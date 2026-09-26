"""HDA hydrogen/methane/benzene/toluene VLE with ePC-SAFT phases (ePC-SAFT issue #73).

Hydrogen and methane are non-condensable (vapor only), as in the standard IDAES HDA
model; benzene and toluene appear in both phases. Phase equilibrium is IDAES SmoothVLE
with LogBubbleDew and log-fugacity equality for benzene and toluene. With vapor-only
components IDAES builds no bubble temperature, and _teq = smooth_min(T, T_dew, eps_2).
Every density, enthalpy and ln(phi_i) comes from the installed Engine wheel (epcsaft_eos). The Engine enthalpy already
contains the NASA9 formation enthalpy, so IDAES adds no formation term and the
reactor adds no separate heat of reaction.
"""

from pyomo.environ import units as pyunits
from idaes.core import Component, LiquidPhase, PhaseType, VaporPhase
from idaes.models.properties.modular_properties.state_definitions import FpTPxpc, FTPx
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import LogBubbleDew
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity
from idaes.models.properties.modular_properties.pure import RPP5

from idaes_examples.mod.hda.epcsaft_eos import COMPONENTS, EPCSAFT
from idaes_examples.mod.hda.hda_ideal_VLE_modular import thermo_config as ideal_config

ELEMENTS = {"hydrogen": {"H": 2}, "methane": {"C": 1, "H": 4},
            "benzene": {"C": 6, "H": 6}, "toluene": {"C": 7, "H": 8}}


NONCONDENSABLE = ("hydrogen", "methane")


def _component(name):
    ideal = ideal_config["components"][name]["parameter_data"]
    if name in NONCONDENSABLE:
        return {"type": Component, "elemental_composition": ELEMENTS[name],
                "valid_phase_types": [PhaseType.vaporPhase],
                "parameter_data": {"mw": ideal["mw"]}}
    return {
        "type": Component,
        "elemental_composition": ELEMENTS[name],
        # RPP5 vapor pressure and critical constants only seed IDAES's dew initial
        # guesses; no model equation uses them.
        "pressure_sat_comp": RPP5,
        "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
        "parameter_data": {k: ideal[k] for k in (
            "mw", "pressure_sat_comp_coeff", "temperature_crit", "pressure_crit")},
    }


thermo_config = {
    "components": {c: _component(c) for c in COMPONENTS},
    "phases": {
        "Liq": {"type": LiquidPhase, "equation_of_state": EPCSAFT},
        "Vap": {"type": VaporPhase, "equation_of_state": EPCSAFT},
    },
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    # Phase flows and compositions are the state, as in the ideal baseline: flash outlets carry
    # their phase split across arcs, and zero-flow phases are simple bounds.
    "state_definition": FpTPxpc,
    "state_bounds": {
        "flow_mol_phase": (0, 0.5, 100, pyunits.mol / pyunits.s),
        "temperature": (10, 300, 1500, pyunits.K),
        "pressure": (5e4, 1e5, 1e7, pyunits.Pa),
    },
    "pressure_ref": (1e5, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": LogBubbleDew,
    "include_enthalpy_of_formation": False,  # the Engine enthalpy already carries it
}

# The benzene/toluene binary (same Engine mixture, hydrogen and methane at weight 0.0) for the
# temperature-path check: no non-condensables, so SmoothVLE builds both T_bub (eps_1) and T_dew
# (eps_2). FTPx makes one fixed-(F, T, P, z) state block a complete flash.
thermo_config_binary = {
    **thermo_config, "components": {c: _component(c) for c in ("benzene", "toluene")},
    "state_definition": FTPx,
    "state_bounds": {"flow_mol": (0, 1, 100, pyunits.mol / pyunits.s),
                     **{k: thermo_config["state_bounds"][k] for k in ("temperature", "pressure")}}}
