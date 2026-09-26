"""ePC-SAFT HDA property and phase-root qualification (ePC-SAFT issue #73, study steps 3-4).

Writes hda_epcsaft_results/: IDAES-vs-public property transport, the callback temperature
domain, LogBubbleDew bubble/dew solves on single ePC-SAFT state blocks at the ideal-baseline
stream compositions (audited), the bubble sum sum(K_i z_i) at each stream pressure, and
the corrected ideal-baseline simulation/optimization.
Run: python -m idaes_examples.mod.hda.hda_epcsaft_qualification
"""

import csv
import hashlib
import json
import math
from pathlib import Path

import epcsaft
from pyomo.environ import ConcreteModel, Constraint, value
from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock

from idaes_examples.mod.hda import hda_two_flash_epcsaft as flowsheet
from idaes_examples.mod.hda.epcsaft_eos import COMPONENTS, LIBRARY, NASA9, PARAMETERS, audit_row, mixture
from idaes_examples.mod.hda.hda_epcsaft_VLE import thermo_config

RESULTS = Path(__file__).resolve().parents[3] / "hda_epcsaft_results"
PP = ("Vap", "Liq")
SUBPROBLEM = {"bub": ("temperature_bubble", "tbub"), "dew": ("temperature_dew", "tdew")}
GUESSES = {"bub": ((0.999, 1e-3, 1e-8, 1e-8), (0.5, 0.3, 0.1, 0.1)),  # incipient vapor: light, mixed
           "dew": ((1e-4, 1e-3, 0.5, 0.4989),)}  # incipient liquid: aromatic


def state_block(T, P, z):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = GenericParameterBlock(**thermo_config)
    m.fs.sb = m.fs.props.build_state_block([0], defined_state=True)
    sb = m.fs.sb[0]
    sb.flow_mol.fix(1)
    sb.temperature.fix(T)
    sb.pressure.fix(P)
    for c, v in zip(COMPONENTS, z):
        sb.mole_frac_comp[c].fix(v)
        sb.log_mole_frac_comp[c].set_value(math.log(v))
    return m, sb


def transport_rows(mix):
    """IDAES expressions vs Mixture.state: density, enthalpy, ln phi_i (issue qualification states)."""
    rows = []
    for T, P, w in ((303.2, 350000.0, (0.30, 0.02, 1e-5, 1e-5)), (600.0, 350000.0, (0.10, 0.05, 0.25, 0.60)),
                    (325.0, 350000.0, (0.05, 0.10, 0.45, 0.40)), (375.0, 105000.0, (0.05, 0.10, 0.45, 0.40))):
        x = [v / sum(w) for v in w]
        m, sb = state_block(T, P, x)
        for p, phase in (("Liq", "liquid"), ("Vap", "vapor")):
            for c, v in zip(COMPONENTS, x):
                sb.mole_frac_phase_comp[p, c].set_value(v)
                sb.log_mole_frac_phase_comp[p, c].set_value(math.log(v))
            sb._teq[PP].set_value(T)
            eos = sb.params.get_phase(p).config.equation_of_state
            got = [value(sb.dens_mol_phase[p]), value(sb.enth_mol_phase[p])] + [
                value(eos.log_fug_phase_comp_eq(sb, p, c, PP)) - math.log(v * P / 1e5) for c, v in zip(COMPONENTS, x)]
            st = mix.state(T, P=P, x=x, phase=phase)
            ref = [st.molar_density, st.total_enthalpy, *st.log_fugacity_coefficient]
            rows.append({"T_K": T, "P_Pa": P, "phase": phase, "branch": st.density_branch.name,
                         "stable_roots": st.stable_root_count,
                         "max_rel_diff": max(abs(a - b) / abs(b) for a, b in zip(got, ref))})
    return rows


def callback_domain_rows():
    """Callback availability versus temperature for the gas-feed composition (vapor request)."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = GenericParameterBlock(**thermo_config)
    w = [0.30 / 0.32, 0.02 / 0.32, 1e-5 / 0.32, 1e-5 / 0.32]
    return [{"T_K": T, "observable": obs, "value": getattr(m.fs.props, "_epcsaft_" + kind).evaluate(
        [T, 350000.0, *w, m.fs.props._epcsaft_record, obs, comp, "vapor"])}
        for T in (150.0, 199.9, 200.0, 200.1, 250.0)
        for kind, obs, comp in (("dens", "molar_density", "-"), ("enth", "enthalpy", "-"),
                                ("lnphi", "log_fugacity_coefficient", "0"))]


def bubble_dew_rows(mix, streams):
    """IDAES LogBubbleDew bubble/dew equations alone, from a fixed start grid, then audited."""
    flowsheet.SOLVER_OPTIONS.update(bound_push=1e-8, max_iter=200)  # keep trace starting fractions
    rows = []
    for name, (T, P, z) in streams.items():
        for kind, (tvar, abbrv) in SUBPROBLEM.items():
            for T0, guess in ((T0, g) for T0 in (80, 120, 160, 210, 300, 400) for g in GUESSES[kind]):
                m, sb = state_block(T, P, z)
                keep = {f"eq_{tvar}", f"eq_mole_frac_{abbrv}", f"log_mole_frac_{abbrv}_eqn", "log_mole_frac_comp_eqn"}
                for con in m.component_objects(Constraint, descend_into=True):
                    con.activate() if con.local_name in keep else con.deactivate()
                getattr(sb, tvar)[PP].set_value(T0)
                for c, g in zip(COMPONENTS, guess):
                    getattr(sb, "_mole_frac_" + abbrv)[PP, c].set_value(g)
                    getattr(sb, "log_mole_frac_" + abbrv)[PP, c].set_value(math.log(g))
                res = flowsheet.solve(m)
                Tsol = value(getattr(sb, tvar)[PP])
                inc = [value(getattr(sb, "_mole_frac_" + abbrv)[PP, c]) for c in COMPONENTS]
                x, y = (z, inc) if kind == "bub" else (inc, z)
                audit = audit_row(mix, Tsol, P, x, y) if res["termination"] == "convergenceCriteriaSatisfied" else {}
                print(name, kind, T0, res["termination"], round(Tsol, 3), audit.get("reason"), flush=True)
                rows.append({"stream": name, "kind": kind, "T_stream_K": T, "P_Pa": P, "T_start_K": T0, "incipient_guess": guess,
                             "termination": res["termination"], "T_solution_K": Tsol,
                             "in_stream_envelope": 298 <= Tsol <= 800, **audit})
    return rows


def bubble_sum_rows(mix, streams):
    """Decisive bubble test at the stream pressure: with the liquid at the stream composition z,
    converge the incipient vapor by successive substitution, y_i ~ z_i exp(ln phi_i^L - ln phi_i^V(y)),
    and record S = sum_i K_i z_i. S > 1 means the liquid is unstable to that vapor at T (the bubble
    temperature, if any, is where S = 1); S > 1 at every T means no bubble temperature exists."""
    rows = []
    for name, (_, P, z) in streams.items():
        for T in (12, 15, 20, 25, 30, 40, 60, 80, 100, 120, 150, 180, 210, 240, 270, 300):
            row = {"stream": name, "T_K": T, "P_Pa": P}
            try:
                liquid = mix.state(T, P=P, x=list(z), phase="liquid")
                y = [0.99, 0.0099, 5e-5, 5e-5]  # hydrogen-rich start: y = z is the trivial fixed point
                converged = False
                for iterations in range(1, 201):
                    vapor = mix.state(T, P=P, x=y, phase="vapor")
                    k = [math.exp(a - b) for a, b in zip(liquid.log_fugacity_coefficient, vapor.log_fugacity_coefficient)]
                    total = sum(ki * zi for ki, zi in zip(k, z))
                    new = [max(ki * zi / total, 1e-300) for ki, zi in zip(k, z)]
                    if max(abs(a - b) for a, b in zip(new, y)) < 1e-12:
                        converged = True
                        break
                    y = new
                row.update({"sum_Kz": total, "converged": converged, "iterations": iterations, "y_H2": y[0], "rho_liquid": liquid.molar_density, "rho_vapor": vapor.molar_density,
                            "liquid_roots": liquid.stable_root_count, "trivial": max(abs(a - b) for a, b in zip(y, z)) < 1e-3})
            except Exception as err:  # retained: no density root at this T and P
                row["reason"] = str(err)[:120]
            rows.append(row)
    return rows


def ideal_streams(m):
    """Stream (T, P, z) of the solved ideal baseline, used as HDA composition proxies."""
    fs = m.fs
    blocks = {"gas_feed": fs.I102.properties[0], "liquid_feed": fs.I101.properties[0],
              "M101_out": fs.M101.mixed_state[0], "R101_out": fs.R101.control_volume.properties_out[0],
              "F101_vapor": fs.S101.mixed_state[0], "F102_feed": fs.F102.control_volume.properties_in[0],
              "F102_vapor": fs.P101.properties[0], "F102_liquid": fs.P102.properties[0]}
    out = {}
    for k, b in blocks.items():
        f = [max(value(b.flow_mol_comp[c]), 1e-9) for c in COMPONENTS]
        out[k] = (value(b.temperature), value(b.pressure), [v / sum(f) for v in f])
    return out


def write_csv(name, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(RESULTS / name, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    RESULTS.mkdir(exist_ok=True)
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()  # noqa: E731
    mix = mixture()
    _, ideal = flowsheet.run_ideal_baseline()
    m = flowsheet.build("ideal")  # fresh fixed simulation for the composition proxies
    flowsheet.initialize_ideal(m)
    flowsheet.solve(m)
    streams = ideal_streams(m)
    write_csv("property_transport.csv", transport_rows(mix))
    write_csv("callback_domain.csv", callback_domain_rows())
    write_csv("bubble_dew_idaes.csv", bubble_dew_rows(mix, streams))
    h2_rich = {k: v for k, v in streams.items() if v[2][0] > 0.05}
    write_csv("bubble_sum_at_stream_pressure.csv", bubble_sum_rows(mix, h2_rich))
    wheel = json.loads(next(Path(epcsaft.__file__).parent.parent.glob("epcsaft-*.dist-info"))
                       .joinpath("direct_url.json").read_text())["url"].removeprefix("file://")
    meta = {"wheel": wheel, "wheel_sha256": sha(wheel), "asl_library_sha256": sha(LIBRARY), "parameters_sha256": sha(PARAMETERS),
            "nasa9_sha256": sha(NASA9), "provisional": "benzene/toluene k_ij=-0.003 (Gross pure set)",
            "proxy_streams": {k: {"T_K": T, "P_Pa": P, "z": z} for k, (T, P, z) in streams.items()}}
    (RESULTS / "ideal_baseline.json").write_text(json.dumps(ideal, indent=1))
    (RESULTS / "qualification_meta.json").write_text(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
