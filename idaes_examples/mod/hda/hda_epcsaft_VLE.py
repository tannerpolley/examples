"""HDA hydrogen/methane/benzene/toluene VLE with ePC-SAFT phases (ePC-SAFT issue #73).

All four species are admitted to both phases. Phase equilibrium is IDAES SmoothVLE
with LogBubbleDew and log-fugacity equality; every density, enthalpy and ln(phi_i)
comes from the installed Engine wheel (epcsaft_eos). The Engine enthalpy already
contains the NASA9 formation enthalpy, so IDAES adds no formation term and the
reactor adds no separate heat of reaction.
"""

from pyomo.environ import units as pyunits
from idaes.core import Component, LiquidPhase, VaporPhase
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import LogBubbleDew
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity
from idaes.models.properties.modular_properties.pure import RPP5

from idaes_examples.mod.hda.epcsaft_eos import COMPONENTS, EPCSAFT
from idaes_examples.mod.hda.hda_ideal_VLE_modular import thermo_config as ideal_config

ELEMENTS = {"hydrogen": {"H": 2}, "methane": {"C": 1, "H": 4},
            "benzene": {"C": 6, "H": 6}, "toluene": {"C": 7, "H": 8}}


def _component(name):
    ideal = ideal_config["components"][name]["parameter_data"]
    return {
        "type": Component,
        "elemental_composition": ELEMENTS[name],
        # RPP5 vapor pressure and critical constants only seed IDAES's bubble/dew initial
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
    "state_definition": FTPx,
    "state_bounds": {
        "flow_mol": (0, 1, 100, pyunits.mol / pyunits.s),
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
