"""ePC-SAFT HDA property qualification (ePC-SAFT issue #73, study steps 2-3).

IDAES-vs-public property transport at the issue states (liquid calls carry vapor-only
hydrogen and methane at weight 0.0), the callback temperature domain, and the
benzene/toluene temperature path solved by IDAES SmoothVLE/LogBubbleDew at a given
smoothing factor. hda_epcsaft_campaign.py writes the results.
"""

import math

from pyomo.environ import ConcreteModel, value
from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock

from idaes_examples.mod.hda import hda_two_flash_epcsaft as flowsheet
from idaes_examples.mod.hda.epcsaft_eos import COMPONENTS, PHASES, ROWS, audit_row
from idaes_examples.mod.hda.hda_epcsaft_VLE import thermo_config, thermo_config_aromatics, thermo_config_binary

PP = flowsheet.PP
LIQUID_ABSENT = ("hydrogen", "methane")
TRANSPORT_STATES = (  # (T, P, weights [H2, CH4, benzene, toluene]) from the issue
    ("gas_feed", 303.2, 350000.0, (0.30, 0.02, 1e-5, 1e-5)),
    ("liquid_feed", 303.2, 350000.0, (1e-5, 1e-5, 1e-5, 0.30)),
    ("asymmetric", 600.0, 350000.0, (0.10, 0.05, 0.25, 0.60)),
    ("separator_F101", 325.0, 350000.0, (0.05, 0.10, 0.45, 0.40)),
    ("separator_F102", 375.0, 105000.0, (0.05, 0.10, 0.45, 0.40)),
    # The F102 train's aromatics package (zero hydrogen and methane weights in both phases).
    ("aromatics_F102", 363.3, 105000.0, {"Liq": (0, 0, 0.6447, 0.3553), "Vap": (0, 0, 0.8163, 0.1837)}))
BINARY_P, BINARY_Z = 105000.0, (0.5, 0.5)
BINARY_GRID = [300.0 + 5 * k for k in range(31)]


def state_block(config, defined_state):
    m = ConcreteModel()
    m.package = "epcsaft"  # flowsheet.solve options
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = GenericParameterBlock(**config)
    m.fs.sb = m.fs.props.build_state_block([0], defined_state=defined_state)
    return m, m.fs.sb[0]


def _normalized(w):
    return [v / sum(w) for v in w]


def transport_rows(mix):
    """IDAES density, enthalpy and log-fugacity expressions and the raw ln phi_i callback vs
    Mixture.state, each phase. The liquid weights are (0, 0, w_benzene, w_toluene): the four
    ln phi_i there include the exact infinite-dilution limits for hydrogen and methane."""
    rows = []
    for name, T, P, w in TRANSPORT_STATES:
        m, sb = state_block(thermo_config_aromatics if isinstance(w, dict) else thermo_config, defined_state=False)
        sb.temperature.set_value(T)
        sb.pressure.set_value(P)
        sb._teq[PP].set_value(T)
        for p, phase in (("Liq", "liquid"), ("Vap", "vapor")):
            x = _normalized(w[p] if isinstance(w, dict) else
                            [0.0 if (p == "Liq" and c in LIQUID_ABSENT) else v for c, v in zip(COMPONENTS, w)])
            for c, v in zip(COMPONENTS, x):
                if (p, c) in sb.phase_component_set:
                    sb.mole_frac_phase_comp[p, c].set_value(v)
            st = mix.state(T, P=P, x=x, phase=phase)
            eos = sb.params.get_phase(p).config.equation_of_state
            lnphi = [sb.params._epcsaft_lnphi.evaluate([T, P, *x, sb.params._epcsaft_record,
                                                         str(ROWS["lnphi"][0] + i), str(PHASES[phase])])
                     for i in range(4)]
            pairs = [("molar_density", value(sb.dens_mol_phase[p]), st.molar_density, 1e4),
                     ("enthalpy", value(sb.enth_mol_phase[p]), st.total_enthalpy, 1e4)]
            pairs += [(f"ln_phi_{c}", g, r, 1.0) for c, g, r in zip(COMPONENTS, lnphi, st.log_fugacity_coefficient)]
            pairs += [(f"log_fug_eq_{c}", value(eos.log_fug_phase_comp_eq(sb, p, c, PP)),
                       math.log(x[i] * P / 1e5) + st.log_fugacity_coefficient[i], 1.0)
                      for i, c in enumerate(COMPONENTS) if (p, c) in sb.phase_component_set]
            for quantity, got, ref, scale in pairs:  # A12 fixed scales: rho 1e4 mol/m3, h 1e4 J/mol, ln 1
                diff, limit = abs(got - ref) / scale, 1e-12 + 1e-10 * abs(ref / scale)
                rows.append({"state": name, "T_K": T, "P_Pa": P, "phase": phase, "quantity": quantity,
                             "idaes": got, "public": ref, "scale": scale, "scaled_diff": diff, "limit": limit,
                             "excess": diff - limit, "passed": diff <= limit,
                             "branch": st.density_branch.name, "stable_roots": st.stable_root_count})
    return rows


def callback_domain_rows():
    """Callback availability versus temperature for the gas-feed composition (vapor request)."""
    pb = state_block(thermo_config, defined_state=True)[1].params
    w = _normalized((0.30, 0.02, 1e-5, 1e-5))
    return [{"T_K": T, "row": row, "value": getattr(pb, "_epcsaft_" + kind).evaluate(
        [T, 350000.0, *w, pb._epcsaft_record, str(row), str(PHASES["vapor"])])}
        for T in (150.0, 199.9, 200.0, 200.1, 250.0) for kind, (row, _) in ROWS.items()]


def binary_path_rows(mix, factor):
    """Equimolar benzene/toluene at 105 kPa, 300-450 K every 5 K plus the computed bubble and dew
    temperatures, freshly solved by IDAES SmoothVLE/LogBubbleDew at `factor` x the default
    smoothing widths; the bubble, dew and _teq rows are audited (weights (0, 0, x_b, x_t))."""
    m, sb = state_block(thermo_config_binary, defined_state=True)  # FTPx: always flashes
    sb.flow_mol.fix(1)
    sb.temperature.fix(BINARY_GRID[0])
    sb.pressure.fix(BINARY_P)
    for c, v in zip(sb.component_list, BINARY_Z):
        sb.mole_frac_comp[c].fix(v)
    flowsheet.set_smoothing(m, factor)
    m.fs.sb.initialize()

    def fractions(var, keys):  # four-component vector, hydrogen and methane at 0.0
        return [0.0, 0.0, *(value(var[k]) for k in keys)]

    def solve_at(T):
        sb.temperature.fix(T)
        res = flowsheet.solve(m)
        comps = sb.component_list
        z = fractions(sb.mole_frac_comp, comps)
        x = fractions(sb.mole_frac_phase_comp, [("Liq", c) for c in comps])
        y = fractions(sb.mole_frac_phase_comp, [("Vap", c) for c in comps])
        row = {"factor": factor, "T_K": T, "termination": res["termination"], "iterations": res["iterations"],
               "T_bubble_K": value(sb.temperature_bubble[PP]), "T_dew_K": value(sb.temperature_dew[PP]),
               "T_eq_K": value(sb._teq[PP]), "vapor_fraction": value(sb.phase_frac["Vap"]),
               "x_benzene": x[2], "y_benzene": y[2]}
        points = {"bubble": (row["T_bubble_K"], z, fractions(sb._mole_frac_tbub, [(*PP, c) for c in comps])),
                  "dew": (row["T_dew_K"], fractions(sb._mole_frac_tdew, [(*PP, c) for c in comps]), z),
                  "teq": (row["T_eq_K"], x, y)}
        for kind, (Ta, xa, ya) in points.items():
            a = audit_row(mix, Ta, BINARY_P, xa, ya)
            row.update({f"{kind}_{k}": a.get(k) for k in (
                "accepted", "reason", "fugacity_residual", "branch_liquid", "branch_vapor", "rho_liquid", "rho_vapor")})
        return row

    rows = [solve_at(T) for T in BINARY_GRID]
    return rows + [solve_at(rows[-1][k]) for k in ("T_bubble_K", "T_dew_K")]
