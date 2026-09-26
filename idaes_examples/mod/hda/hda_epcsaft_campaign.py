"""ePC-SAFT HDA two-flash simulation, optimization and smoothing campaign (ePC-SAFT issue #73).

For each smoothing factor (1, 0.1, 0.01 x the IDAES SmoothVLE widths): the benzene/toluene
temperature path, the fixed two-flash recycle simulation started from the solved ideal
baseline, and the operating-cost optimization from the three fixed starts. Every solution
gets balances, pressure-equation residuals and a post-solve audit of every dew and _teq
point; every optimization gets independent KKT evidence (hda_epcsaft_kkt).
Run (one process, ~hours on two shared cores):
    OMP_NUM_THREADS=1 PYTHONPATH=. python -m idaes_examples.mod.hda.hda_epcsaft_campaign
"""

import csv
import gzip
import hashlib
import json
import sys
import time
from pathlib import Path

import epcsaft
from pyomo.environ import Block, value
from pyomo.network import Arc
from idaes.core import StateBlockData
from idaes.core.util import from_json, to_json

from idaes_examples.mod.hda import hda_epcsaft_qualification as qual
from idaes_examples.mod.hda import hda_two_flash_epcsaft as flowsheet
from idaes_examples.mod.hda.epcsaft_eos import COMPONENTS, LIBRARY, NASA9, PARAMETERS, audit_row, mixture
from idaes_examples.mod.hda.hda_epcsaft_kkt import kkt
from idaes_examples.mod.hda.hda_epcsaft_VLE import ELEMENTS

RESULTS = Path(__file__).resolve().parents[3] / "hda_epcsaft_results"
FACTORS = (1.0, 0.1, 0.01)
PP = flowsheet.PP
ENVELOPE = (298.0, 800.0)  # K, HDA stream-state envelope
STOICH = {"hydrogen": -1, "methane": 1, "benzene": 1, "toluene": -1}


def envelope_class(T):
    return "below" if T < ENVELOPE[0] else ("above" if T > ENVELOPE[1] else "within")


def balances(m):
    """Independent component, atom and energy residuals from solved stream values (mol/s, W):
    the whole flowsheet (feeds -> products) and every control volume."""
    fs = m.fs
    cv = {u: (getattr(fs, u).control_volume.properties_in[0], getattr(fs, u).control_volume.properties_out[0])
          for u in ("H101", "R101", "F101", "F102", "C101")}
    extent = value(fs.R101.rate_reaction_extent[0, "hydrodealkylation"])
    duty = {u: value(getattr(fs, u).heat_duty[0]) for u in ("H101", "R101", "F101", "F102")}
    duty["C101"] = value(fs.C101.work_mechanical[0])
    units = {u: ([i], [o], extent if u == "R101" else 0.0, duty[u]) for u, (i, o) in cv.items()}
    units["M101"] = ([fs.M101.inlet_1_state[0], fs.M101.inlet_2_state[0], fs.M101.inlet_3_state[0]],
                     [fs.M101.mixed_state[0]], 0.0, 0.0)
    units["S101"] = ([fs.S101.mixed_state[0]], [fs.S101.purge_state[0], fs.S101.recycle_state[0]], 0.0, 0.0)
    units["flowsheet"] = ([fs.I101.properties[0], fs.I102.properties[0]],
                          [fs.P101.properties[0], fs.P102.properties[0], fs.P103.properties[0]],
                          extent, sum(duty.values()))
    out = {}
    for u, (ins, outs, ext, q) in units.items():
        net = {c: value(sum(b.flow_mol_comp[c] for b in ins) - sum(b.flow_mol_comp[c] for b in outs)) + STOICH[c] * ext
               for c in COMPONENTS}
        atoms = {e: sum(ELEMENTS[c].get(e, 0) * (net[c] - STOICH[c] * ext) for c in COMPONENTS) for e in ("C", "H")}
        energy = value(sum(b.flow_mol_phase[p] * b.enth_mol_phase[p] for b in ins for p in b.phase_list)
                       - sum(b.flow_mol_phase[p] * b.enth_mol_phase[p] for b in outs for p in b.phase_list)) + q
        out[u] = {"component_mol_s": max(abs(v) for v in net.values()),
                  "atom_mol_s": max(abs(v) for v in atoms.values()), "energy_W": abs(energy)}
    return out


def pressure_residuals(m):
    """|residual| of every pressure equation (Pa): arcs, unit specifications, mixer minimum."""
    fs = m.fs
    res = {a.local_name: abs(value(a.source.pressure[0] - a.destination.pressure[0]))
           for a in fs.component_data_objects(Arc)}
    for u in ("H101", "R101", "F101", "F102", "C101"):
        c = getattr(fs, u).control_volume
        dp = value(c.deltaP[0]) if hasattr(c, "deltaP") else 0.0
        res[u] = abs(value(c.properties_out[0].pressure - c.properties_in[0].pressure) - dp)
    res["F101_deltaP_spec"] = abs(value(fs.F101.deltaP[0]))
    res["C101_outlet_spec"] = abs(value(fs.C101.outlet.pressure[0]) - 350000.0)
    inlets = [value(getattr(fs.M101, f"inlet_{k}_state")[0].pressure) for k in (1, 2, 3)]
    res["M101_minimum"] = abs(value(fs.M101.mixed_state[0].pressure) - min(inlets))
    res["S101"] = max(abs(value(getattr(fs.S101, f"{o}_state")[0].pressure - fs.S101.mixed_state[0].pressure))
                      for o in ("purge", "recycle"))
    return res


def audit(m, mix):
    """Audit the dew point and the _teq point of every state block (public API only)."""
    rows = []
    for sb in m.component_data_objects(Block, descend_into=True):
        if not (isinstance(sb, StateBlockData) and hasattr(sb, "_teq")):
            continue
        z = [value(sb.mole_frac_comp[c]) for c in COMPONENTS]
        x = [value(sb.mole_frac_phase_comp["Liq", c]) if ("Liq", c) in sb.phase_component_set else 0.0
             for c in COMPONENTS]
        y = [value(sb.mole_frac_phase_comp["Vap", c]) for c in COMPONENTS]
        tdew = [value(sb._mole_frac_tdew[(*PP, c)]) if c in ("benzene", "toluene") else 0.0 for c in COMPONENTS]
        base = {"block": sb.name, "T_K": value(sb.temperature), "P_Pa": value(sb.pressure),
                "T_class": envelope_class(value(sb.temperature)),
                "vapor_fraction": value(sb.phase_frac["Vap"]), "liquid_fraction": value(sb.phase_frac["Liq"])}
        for kind, T, xa, ya in (("dew", value(sb.temperature_dew[PP]), tdew, z), ("teq", value(sb._teq[PP]), x, y)):
            a = audit_row(mix, T, value(sb.pressure), xa, ya)
            rows.append({**base, "point": kind, "T_point_K": T, "T_point_class": envelope_class(T),
                         "x": xa, "y": ya, **a})
    return rows


def streams(m):
    """Stream and phase table of every arc: T, P, flows, phase fractions, compositions, densities."""
    rows = []
    for arc in m.fs.component_data_objects(Arc):
        sb = arc.source.temperature[0].parent_block()
        row = {"stream": arc.local_name, "from": arc.source.parent_block().local_name,
               "to": arc.destination.parent_block().local_name, "T_K": value(sb.temperature),
               "P_Pa": value(sb.pressure), "flow_mol_s": value(sb.flow_mol),
               "vapor_fraction": value(sb.phase_frac["Vap"])}
        row.update({f"F_{c}": value(sb.flow_mol_comp[c]) for c in COMPONENTS})
        for p in ("Liq", "Vap"):
            row.update({f"{p}_x_{c}": value(sb.mole_frac_phase_comp[p, c]) if (p, c) in sb.phase_component_set
                        else 0.0 for c in COMPONENTS})
            row[f"{p}_rho_mol_m3"] = value(sb.dens_mol_phase[p])
            row[f"{p}_h_J_mol"] = value(sb.enth_mol_phase[p])
        rows.append(row)
    return rows


def signature(audit_rows):
    """Topology for branch matching: ordered roots and present phases at every audited point."""
    return [(r["block"], r["point"], r.get("branch_liquid"), r.get("branch_vapor"),
             r.get("rho_liquid", 0) > r.get("rho_vapor", 0), r["vapor_fraction"] > 1e-6, r["liquid_fraction"] > 1e-6)
            for r in audit_rows]


def evidence(m, mix):
    rows = audit(m, mix)
    return {**flowsheet.summary(m), "balances": balances(m), "pressure_residuals_Pa": pressure_residuals(m),
            "audit_accepted": all(r["accepted"] for r in rows),
            "audit_rejected": [f'{r["block"]}:{r["point"]}:{r["reason"]}' for r in rows if not r["accepted"]],
            "audit_max_fugacity_residual": max(r.get("fugacity_residual", float("inf")) for r in rows),
            "_audit": rows, "_streams": streams(m), "_signature": signature(rows)}


def campaign(ref, mix, factor):
    """Fixed simulation from the ideal solution, then the optimization from each fixed start."""
    m = flowsheet.build("epcsaft")
    flowsheet.set_smoothing(m, factor)
    flowsheet.fix_simulation(m)
    flowsheet.initialize_from(m, ref)
    sim = {"factor": factor, **flowsheet.solve(m), **evidence(m, mix)}
    print("simulation", factor, sim["termination"], sim["iterations"], round(sim["wall_s"]), flush=True)
    base = to_json(m, return_dict=True)
    opts = []
    for start in flowsheet.STARTS:
        from_json(m, sd=base)  # initial point: the fixed simulation with the decisions at the start
        flowsheet.fix_simulation(m, start)
        flowsheet.set_optimization(m)
        rec = {"factor": factor, "start": start, **flowsheet.solve(m)}
        rec.update(evidence(m, mix))
        rec["kkt"] = kkt(m, flowsheet.decisions(m))
        print("optimization", factor, start, rec["termination"], rec["iterations"], round(rec["wall_s"]),
              rec["objective_USD_per_yr"], rec["kkt"]["label"], rec["kkt"]["norms"], flush=True)
        opts.append(rec)
    return sim, opts


def compare(ref, other):
    """Smoothing sensitivity of one solution versus the factor-1 solution (issue criteria)."""
    rel = lambda a, b, floor: abs(a - b) <= max(1e-3 * abs(a), floor)  # noqa: E731
    checks = {f"duty_{u}": rel(ref["duties_W"][u], other["duties_W"][u], 1.0) for u in ref["duties_W"]}
    checks["C101_work"] = rel(ref["C101_work_W"], other["C101_work_W"], 1.0)
    checks["operating_cost"] = rel(ref["objective_USD_per_yr"], other["objective_USD_per_yr"], 1.0)
    for k in ("purity", "recovery"):
        checks[k] = abs(ref[k] - other[k]) <= 1e-4
    for u in ref["vapor_fraction"]:
        checks[f"vapor_fraction_{u}"] = abs(ref["vapor_fraction"][u] - other["vapor_fraction"][u]) <= 1e-4
    checks["topology"] = ref["_signature"] == other["_signature"]
    return {"passed": all(checks.values()), "failed": [k for k, ok in checks.items() if not ok],
            "max_abs_change": {"objective_USD_per_yr": abs(ref["objective_USD_per_yr"] - other["objective_USD_per_yr"]),
                               "purity": abs(ref["purity"] - other["purity"]),
                               "recovery": abs(ref["recovery"] - other["recovery"])}}


def compare_binary(ref_rows, rows):
    """Binary path sensitivity: vapor fraction within 1e-4 and the same branches at each point."""
    out = []
    for a, b in zip(ref_rows, rows):
        same = all(a[f"{k}_branch_{p}"] == b[f"{k}_branch_{p}"] for k in ("bubble", "dew", "teq")
                   for p in ("liquid", "vapor"))
        dv = abs(a["vapor_fraction"] - b["vapor_fraction"])
        out.append({"factor": b["factor"], "T_K_ref": a["T_K"], "T_K": b["T_K"], "vapor_fraction_change": dv,
                    "same_branches": same, "passed": same and dv <= 1e-4})
    return out


def write_csv(name, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(RESULTS / name, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def public(rec):
    return {k: v for k, v in rec.items() if not k.startswith("_")}


def main(factors=FACTORS):
    RESULTS.mkdir(exist_ok=True)
    t0 = time.perf_counter()
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()  # noqa: E731
    mix = mixture()
    write_csv("property_transport.csv", qual.transport_rows(mix))
    write_csv("callback_domain.csv", qual.callback_domain_rows())
    ideal_m, ideal = flowsheet.run_ideal_baseline()
    (RESULTS / "ideal_baseline.json").write_text(json.dumps(ideal, indent=1))
    ref = flowsheet.build("ideal")  # the solved fixed ideal simulation initializes every ePC-SAFT run
    flowsheet.initialize_ideal(ref)
    flowsheet.solve(ref)

    binary, sims, opts, sens, audits, stream_rows, primal_dual = {}, {}, [], [], [], [], []
    for f in factors:
        binary[f] = qual.binary_path_rows(mix, f)
        sims[f], o = campaign(ref, mix, f)
        opts += o
        for rec in (sims[f], *o):
            case = "simulation" if "start" not in rec else "start_" + "_".join(str(v) for v in rec["start"])
            audits += [{"factor": f, "case": case, **r} for r in rec["_audit"]]
            stream_rows += [{"factor": f, "case": case, **r} for r in rec["_streams"]]
            if "kkt" in rec:
                pdx = rec["kkt"].pop("_primal_dual")
                primal_dual += [(f, case, "variable", n, v, s) for n, v, s in pdx["variables"]]
                primal_dual += [(f, case, "equality_multiplier_scaled", n, v, s) for n, v, s in pdx["equalities"]]
    for f in factors[1:]:
        sens.append({"factor": f, "case": "simulation", **compare(sims[factors[0]], sims[f])})
        for a, b in zip([o for o in opts if o["factor"] == factors[0]], [o for o in opts if o["factor"] == f]):
            sens.append({"factor": f, "case": f"start {a['start']}", **compare(a, b)})
    binary_sens = [r for f in factors[1:] for r in compare_binary(binary[factors[0]], binary[f])]

    write_csv("binary_path.csv", [r for f in factors for r in binary[f]])
    write_csv("binary_path_sensitivity.csv", binary_sens)
    write_csv("audit.csv", audits)
    write_csv("streams.csv", stream_rows)
    with gzip.open(RESULTS / "optimization_primal_dual.csv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["factor", "case", "kind", "name", "value", "scale"])
        w.writerows(primal_dual)
    wheel = json.loads(next(Path(epcsaft.__file__).parent.parent.glob("epcsaft-*.dist-info"))
                       .joinpath("direct_url.json").read_text())["url"].removeprefix("file://")
    out = {"provenance": {"wheel": wheel, "wheel_sha256": sha(wheel), "asl_library_sha256": sha(LIBRARY),
                          "parameters_sha256": sha(PARAMETERS), "nasa9_sha256": sha(NASA9),
                          "python": sys.version.split()[0]},
           "simulations": [public(sims[f]) for f in factors], "optimizations": [public(o) for o in opts],
           "sensitivity": sens,
           "binary_path_sensitivity_passed": all(r["passed"] for r in binary_sens),
           "ideal_baseline": {"simulation": public(ideal["simulation"]),
                              "optimum": public(ideal["optimizations"][0])},
           "wall_s": time.perf_counter() - t0}
    (RESULTS / "campaign.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main(tuple(float(a) for a in sys.argv[1:]) or FACTORS)
