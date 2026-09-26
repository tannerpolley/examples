"""ePC-SAFT equation of state for IDAES modular properties (ePC-SAFT issue #73).

Every thermodynamic value (phase molar density, total molar enthalpy on the NASA9
formation basis at 1 bar, component ln(phi_i)) comes from the installed Engine
wheel's ASL callback ``epcsaft_pressure_property`` with exact first and second
derivatives. This module only maps IDAES variables to the callback arguments and
assembles ln(f_i / P_ref) = ln x_i + ln(P / P_ref) + ln(phi_i), as the IDAES cubic
EoS does. Liquid requests take the highest mechanically stable density root and
vapor requests the lowest; a single stable root is returned for either request.
"""

import json
import math
from pathlib import Path

import epcsaft
from epcsaft import IdealCorrelation, IdealInterval, IdealNASA9, Mixture, Parameters
from epcsaft import ThermochemistryRecord
from pyomo.environ import ExternalFunction, log, units as pyunits
from idaes.models.properties.modular_properties.eos.eos_base import EoSBase

INPUTS = Path(__file__).parent / "epcsaft_inputs"
PARAMETERS = INPUTS / "hda-pc-saft.json"
NASA9 = INPUTS / "nasa9-records.json"
COMPONENTS = ("hydrogen", "methane", "benzene", "toluene")  # fixed callback order
LIBRARY = str(Path(epcsaft.__file__).parent / "libepcsaft_asl.so")
X_MIN = 1e-12  # callback domain: every composition weight >= 1e-12, no clipping
OBSERVABLES = {  # IDAES name -> (callback observable, result units)
    "dens": ("molar_density", pyunits.mol / pyunits.m**3),
    "enth": ("enthalpy", pyunits.J / pyunits.mol),
    "lnphi": ("log_fugacity_coefficient", pyunits.dimensionless),
}


def thermochemistry():
    """NASA9 gas records with the source gas constant converted to the Engine constant."""
    records = json.loads(NASA9.read_text())
    scale = records["gas_constant_J_per_mol_K"] / epcsaft.GAS_CONSTANT_J_PER_MOL_K
    p0 = records["standard_pressure_Pa"]
    species = {s["species_id"]: s["intervals"] for s in records["species"]}
    intervals = []
    for c in COMPONENTS:
        last = len(species[c]) - 1
        intervals.append([
            IdealInterval(iv["lower_K"], iv["upper_K"], True, k == last,
                          IdealCorrelation(IdealNASA9([a * scale for a in iv["coefficients"]]), p0))
            for k, iv in enumerate(species[c])])
    return ThermochemistryRecord(p0, intervals)


def mixture():
    """The single Engine mixture every callback and audit uses."""
    parameters = Parameters.from_json(PARAMETERS, components=list(COMPONENTS))
    return Mixture(parameters, thermochemistry=thermochemistry())


def _phase(b, p):
    return "liquid" if b.params.get_phase(p).is_liquid_phase() else "vapor"


def _call(b, kind, p, T, P, x, j="-"):
    """One callback term; x maps component name -> mole-fraction expression."""
    f = getattr(b.params, "_epcsaft_" + kind)
    for c in COMPONENTS:
        if hasattr(x[c], "setlb"):  # the callback consumes this variable: enforce its domain
            x[c].setlb(X_MIN)
    comp = j if j == "-" else str(COMPONENTS.index(j))
    return f(T, P, *(x[c] for c in COMPONENTS), b.params._epcsaft_record,
             OBSERVABLES[kind][0], comp, _phase(b, p))


def _phase_x(b, p):
    return {c: b.mole_frac_phase_comp[p, c] for c in COMPONENTS}


def _log_fug_bubble_dew(b, p, j, pp, pt_var):
    """ln(f_j / P_ref) at a bubble/dew point, choosing T, P and x as the IDAES cubic EoS does."""
    name = pt_var.local_name
    T, P = (pt_var[pp], b.pressure) if name.startswith("temperature") else (b.temperature, pt_var[pp])
    abbrv = name[0] + ("bub" if name.endswith("bubble") else "dew")
    incipient = b.params.get_phase(p).is_vapor_phase() == abbrv.endswith("bub")
    if incipient:  # vapor at a bubble point, liquid at a dew point
        mf, log_mf = getattr(b, "_mole_frac_" + abbrv), getattr(b, "log_mole_frac_" + abbrv)
        x, log_xj = {c: mf[pp, c] for c in COMPONENTS}, log_mf[pp, j]
    else:
        x, log_xj = {c: b.mole_frac_comp[c] for c in COMPONENTS}, b.log_mole_frac_comp[j]
    return log_xj + log(P / b.params.pressure_ref) + _call(b, "lnphi", p, T, P, x, j)


class EPCSAFT(EoSBase):
    """Vapor and liquid phases from the Engine ASL property callback."""

    @staticmethod
    def build_parameters(b):
        pb = b.parent_block()
        if hasattr(pb, "_epcsaft_record"):
            return
        if tuple(pb.component_list) != COMPONENTS:
            raise ValueError(f"ePC-SAFT HDA callback requires components {COMPONENTS}")
        pb._epcsaft_record = epcsaft.eos.property_asl_record(mixture())
        arg_units = [pyunits.K, pyunits.Pa] + [pyunits.dimensionless] * 4 + [None] * 4
        for kind, (_, units) in OBSERVABLES.items():
            pb.add_component("_epcsaft_" + kind, ExternalFunction(
                library=LIBRARY, function="epcsaft_pressure_property",
                units=units, arg_units=arg_units))

    @staticmethod
    def common(b, pobj):
        pass

    @staticmethod
    def calculate_scaling_factors(b, pobj):
        pass

    @staticmethod
    def dens_mol_phase(b, p):
        return _call(b, "dens", p, b.temperature, b.pressure, _phase_x(b, p))

    @staticmethod
    def enth_mol_phase(b, p):
        return _call(b, "enth", p, b.temperature, b.pressure, _phase_x(b, p))

    @staticmethod
    def log_fug_phase_comp_eq(b, p, j, pp):
        x = _phase_x(b, p)
        return (b.log_mole_frac_phase_comp[p, j] + log(b.pressure / b.params.pressure_ref)
                + _call(b, "lnphi", p, b._teq[pp], b.pressure, x, j))

    @staticmethod
    def log_fug_phase_comp_Tbub(b, p, j, pp):
        return _log_fug_bubble_dew(b, p, j, pp, b.temperature_bubble)

    @staticmethod
    def log_fug_phase_comp_Tdew(b, p, j, pp):
        return _log_fug_bubble_dew(b, p, j, pp, b.temperature_dew)

    @staticmethod
    def log_fug_phase_comp_Pbub(b, p, j, pp):
        return _log_fug_bubble_dew(b, p, j, pp, b.pressure_bubble)

    @staticmethod
    def log_fug_phase_comp_Pdew(b, p, j, pp):
        return _log_fug_bubble_dew(b, p, j, pp, b.pressure_dew)


def molar_masses():
    """Component molar masses (kg/mol) from the parameter document."""
    doc = json.loads(PARAMETERS.read_text())
    mw = {c["component_id"]: c["fixed"]["molar_mass"]["value"]["magnitude"] for c in doc["components"]}
    return [mw[c] for c in COMPONENTS]


def audit_row(mix, T, P, x, y):
    """Independent #73 audit of a declared liquid x / vapor y row at (T, P), public API only.

    Rejects: fractions outside [1e-12, 1] or phase sums off by > 1e-10; unavailable roots; one
    shared root; compositions closer than 1e-3 (trivial y = x branch); a declared liquid that is
    not the denser phase by mass; max |ln x_i + ln phi_i^L - ln y_i - ln phi_i^V| > 1e-8.
    """
    row = {"fractions_ok": all(X_MIN <= v <= 1 for v in (*x, *y))
           and abs(sum(x) - 1) <= 1e-10 and abs(sum(y) - 1) <= 1e-10,
           "sum_error": max(abs(sum(x) - 1), abs(sum(y) - 1))}
    x, y = [v / sum(x) for v in x], [v / sum(y) for v in y]  # sums checked above; roots at normalized x
    try:
        liq, vap = (mix.state(T, P=P, x=z, phase=ph) for z, ph in ((x, "liquid"), (y, "vapor")))
    except Exception as err:  # an unavailable root is a typed hard stop, never a substitute
        return {**row, "accepted": False, "reason": f"root unavailable: {err}"}
    mw = molar_masses()
    row.update(
        rho_liquid=liq.molar_density, rho_vapor=vap.molar_density,
        branch_liquid=liq.density_branch.name, branch_vapor=vap.density_branch.name,
        shared_root=abs(liq.molar_density - vap.molar_density) <= 1e-10 * liq.molar_density,
        composition_distance=max(abs(a - b) for a, b in zip(x, y)),
        mass_ordered=liq.molar_density * sum(a * m for a, m in zip(x, mw))
        > vap.molar_density * sum(a * m for a, m in zip(y, mw)),
        fugacity_residual=max(abs(math.log(x[i]) + liq.log_fugacity_coefficient[i] - math.log(y[i])
                                  - vap.log_fugacity_coefficient[i]) for i in range(len(x))))
    checks = {"fractions": row["fractions_ok"], "distinct_root": not row["shared_root"],
              "composition_distance": row["composition_distance"] >= 1e-3,
              "mass_ordered": row["mass_ordered"], "residual": row["fugacity_residual"] <= 1e-8}
    row["accepted"] = all(checks.values())
    row["reason"] = ",".join(k for k, ok in checks.items() if not ok)
    return row
