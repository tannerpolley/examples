"""HDA two-flash recycle flowsheet: ideal baseline and ePC-SAFT properties (ePC-SAFT issue #73).

Adiabatic stoichiometric reactor at fixed 75 % toluene conversion. Decisions:
H101 outlet T, F101 vapor T, F102 vapor T, F102 outlet P; tutorial operating-cost
objective with overhead-loss, benzene-production and purity constraints.
"""

import time
from pathlib import Path

from pyomo.environ import (Block, ConcreteModel, Constraint, Expression, Objective, Param,
                           TransformationFactory, Var, units as pyunits, value)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from pyomo.contrib.solver.common.results import SolutionStatus
from pyomo.contrib.solver.solvers.ipopt import Ipopt
import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock)
from idaes.core import StateBlockData
from idaes.core.util.initialization import propagate_state
from idaes.models.properties.modular_properties.base.generic_property import _init_Tbub
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.unit_models import (Feed, Flash, Heater, Mixer, PressureChanger, Product,
                                      Separator, StoichiometricReactor, Translator)
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption

from idaes_examples.mod.hda import hda_epcsaft_VLE, hda_ideal_VLE_modular
from idaes_examples.mod.hda.hda_epcsaft_VLE import NONCONDENSABLE
from idaes_examples.mod.hda.hda_flowsheet_extras import fix_inlet_states, manual_propagation
from idaes_examples.mod.hda.hda_reaction_modular import reaction_config

COMPONENTS = ("hydrogen", "methane", "benzene", "toluene")
PP = ("Vap", "Liq")
SMOOTHING_DEFAULTS = {"eps_1": 0.01, "eps_2": 0.0005}  # K, IDAES 2.12 SmoothVLE defaults
UNIT_SCALES = ((pyunits.K, 300.0), (pyunits.Pa, 1e5), (pyunits.mol / pyunits.s, 1.0),
               (pyunits.W, 1e4), (pyunits.dimensionless, 1.0),  # issue #73 fixed KKT scales
               (pyunits.J / pyunits.mol, 1e4))  # the issue's enthalpy comparison scale
OBJECTIVE_SCALE = 1e5  # USD/yr
TRANSLATORS = ("to_gas", "to_aromatics", "to_full")  # full -> gas, full -> aromatics, gas -> full
STARTS = ((550, 325, 375, 107500), (500, 300, 350, 105000), (600, 400, 425, 110000))
SOLVER_OPTIONS = {"nlp_scaling_method": "user-scaling", "ma57_automatic_scaling": "yes",
                  "max_iter": 1000, "tol": 1e-8, "bound_relax_factor": 0.0, "linear_solver": "ma57"}
# ePC-SAFT: primal residuals to 1e-9 in the unit-scaled model (1e-5 W, 1e-4 Pa; the callback's
# ln(phi) precision floor is ~1e-10), and the adaptive barrier update: the fixed update stalls
# with mu = 0.1 while smoothing- and eps-sized phase amounts sit next to their zero bounds.
# Absolute dual and complementarity tolerances: Ipopt's defaults (1, 1e-4) let multipliers of
# 1e9-4e17 on structurally zero rows inflate its relative test (s_d ~ 1.5e8) so that it stops at
# mu ~ 1e-6 before stationarity (ePC-SAFT #73 diagnosis). These tighten the stopping test.
EPCSAFT_OPTIONS = {**SOLVER_OPTIONS, "constr_viol_tol": 1e-9, "mu_strategy": "adaptive",
                   "dual_inf_tol": 1e-6, "compl_inf_tol": 1e-6,
                   "max_iter": 400, "max_cpu_time": 1800}  # caps: 2x the largest converged run


def build(package="ideal"):
    """The tutorial flowsheet with the selected property package ("ideal" | "epcsaft")."""
    m = ConcreteModel()
    m.package = package
    m.fs = FlowsheetBlock(dynamic=False)
    cfg = (hda_ideal_VLE_modular if package == "ideal" else hda_epcsaft_VLE).thermo_config
    fs = m.fs
    fs.thermo_params = props = gas = aromatics = GenericParameterBlock(**cfg)
    if package == "epcsaft":  # split packages joined by translators (owner, 2026-10-13)
        fs.gas_params = gas = GenericParameterBlock(**hda_epcsaft_VLE.thermo_config_gas)
        fs.aromatics_params = aromatics = GenericParameterBlock(**hda_epcsaft_VLE.thermo_config_aromatics)
    fs.reaction_params = GenericReactionParameterBlock(property_package=props, **reaction_config)
    fs.I101, fs.I102 = Feed(property_package=props), Feed(property_package=props)
    fs.M101 = Mixer(property_package=props, num_inlets=3)
    fs.H101 = Heater(property_package=props, has_pressure_change=False, has_phase_equilibrium=True)
    # Formation-based ePC-SAFT enthalpies already carry the reaction enthalpy.
    fs.R101 = StoichiometricReactor(property_package=props, reaction_package=fs.reaction_params,
                                    has_heat_of_reaction=package == "ideal", has_heat_transfer=True,
                                    has_pressure_change=False)
    fs.F101 = Flash(property_package=props, has_heat_transfer=True, has_pressure_change=True)
    fs.S101 = Separator(property_package=gas, ideal_separation=False, outlet_list=["purge", "recycle"])
    fs.C101 = PressureChanger(property_package=gas, compressor=True,
                              thermodynamic_assumption=ThermodynamicAssumption.isothermal)
    fs.F102 = Flash(property_package=aromatics, has_heat_transfer=True, has_pressure_change=True)
    fs.P101, fs.P102 = (Product(property_package=aromatics) for _ in range(2))
    fs.P103 = Product(property_package=gas)
    arcs = {
        "s01": (fs.I101.outlet, fs.M101.inlet_1), "s02": (fs.I102.outlet, fs.M101.inlet_2),
        "s03": (fs.M101.outlet, fs.H101.inlet), "s04": (fs.H101.outlet, fs.R101.inlet),
        "s05": (fs.R101.outlet, fs.F101.inlet), "s06": (fs.F101.vap_outlet, fs.S101.inlet),
        "s07": (fs.F101.liq_outlet, fs.F102.inlet), "s08": (fs.S101.recycle, fs.C101.inlet),
        "s09": (fs.C101.outlet, fs.M101.inlet_3), "s10": (fs.F102.vap_outlet, fs.P101.inlet),
        "s11": (fs.F102.liq_outlet, fs.P102.inlet), "s12": (fs.S101.purge, fs.P103.inlet),
    }
    if package == "epcsaft":
        for arc, name, pin, pout, phase in zip(("s06", "s07", "s09"), TRANSLATORS, (props, props, gas),
                                               (gas, aromatics, props), ("Vap", "Liq", "Vap")):
            setattr(fs, name, tr := Translator(inlet_property_package=pin, outlet_property_package=pout))
            src, dst = arcs.pop(arc)
            arcs[arc + "a"], arcs[arc + "b"] = (src, tr.inlet), (tr.outlet, dst)
            _translate(tr, phase)
    for name, (src, dst) in arcs.items():
        fs.add_component(name, Arc(source=src, destination=dst))
    TransformationFactory("network.expand_arcs").apply_to(m)

    out102, out101 = fs.F102.control_volume.properties_out[0], fs.F101.control_volume.properties_out[0]
    rin, rout = fs.R101.control_volume.properties_in[0], fs.R101.control_volume.properties_out[0]
    fs.purity = Expression(expr=out102.flow_mol_phase_comp["Vap", "benzene"] / (
        out102.flow_mol_phase_comp["Vap", "benzene"] + out102.flow_mol_phase_comp["Vap", "toluene"]))
    fs.recovery = Expression(
        expr=out102.flow_mol_phase_comp["Vap", "benzene"] / rout.flow_mol_comp["benzene"])
    fs.cooling_cost = Expression(expr=0.212e-7 * (-fs.F101.heat_duty[0]) + 0.212e-7 * (-fs.R101.heat_duty[0]))
    fs.heating_cost = Expression(expr=2.2e-7 * fs.H101.heat_duty[0] + 1.9e-7 * fs.F102.heat_duty[0])
    fs.operating_cost = Expression(expr=3600 * 24 * 365 * (fs.heating_cost + fs.cooling_cost))
    fs.R101.control_volume.conversion = conv = Var(initialize=0.75, bounds=(0, 1))
    fs.R101.conv_constraint = Constraint(
        expr=conv * rin.flow_mol_comp["toluene"] == rin.flow_mol_comp["toluene"] - rout.flow_mol_comp["toluene"])
    fs.overhead_loss = Constraint(expr=out101.flow_mol_phase_comp["Vap", "benzene"]
                                  <= 0.20 * rout.flow_mol_comp["benzene"])
    fs.product_flow = Constraint(expr=out102.flow_mol_phase_comp["Vap", "benzene"] >= 0.15)
    fs.product_purity = Constraint(expr=fs.purity >= 0.80)
    fs.objective = Objective(expr=fs.operating_cost)
    for c in (fs.overhead_loss, fs.product_flow, fs.product_purity, fs.objective):
        c.deactivate()
    if package == "epcsaft":
        _epcsaft_numerics(m)
    return m


def _translate(tr, phase):
    """Outlet = inlet for T, P, the present phase's flow and composition; the absent phase is an exact zero flow at
    the equimolar composition (O4). No energy row: both sides evaluate one callback row at identical arguments."""
    i, o = tr.properties_in[0], tr.properties_out[0]
    pairs = [(o.temperature, i.temperature), (o.pressure, i.pressure), (o.flow_mol_phase[phase], i.flow_mol_phase[phase])]
    pairs += [(o.mole_frac_phase_comp[phase, c], i.mole_frac_phase_comp[phase, c]) for c in o.component_list]
    tr.copy = Constraint(range(len(pairs)), rule=lambda b, k: pairs[k][0] == pairs[k][1])
    for p, c in o.phase_component_set:
        if p != phase:
            o.flow_mol_phase[p].fix(0.0)
            o.mole_frac_phase_comp[p, c].fix(0.5)


def _epcsaft_numerics(m):
    """Scaling, absent-phase flows and bounds for the interior-point solve with ePC-SAFT."""
    fs = m.fs
    set_unit_scaling(m)
    # Ideal-separation outlets carry the absent phase as an exact zero, so the F101 liquid's overall
    # hydrogen and methane fractions are exactly 0: they are fixed and their defining rows dropped
    # (z_j F = y_j F_Vap with F_Vap = 0; left free, the row and the active bound make them degenerate).
    fs.F101.split.eps_flow.set_value(0.0)
    fs.F102.split.eps_flow.set_value(0.0)
    for c in NONCONDENSABLE:
        fs.to_aromatics.properties_in[0].mole_frac_comp[c].fix(0.0, skip_validation=True)  # at its 1e-20 floor
        fs.to_aromatics.properties_in[0].mole_frac_comp_eq[c].deactivate()
    # bound_relax_factor = 0 keeps every composition weight inside the callback domain. Phase flows
    # and fractions are never callback inputs; they get Ipopt's default 1e-8 lower-bound relaxation
    # explicitly, because SmoothVLE leaves superheated streams a liquid of O(eps_2^2 / (T - T_dew))
    # (1e-21 to 1e-11 mol/s) whose hard zero bound makes the barrier multipliers diverge.
    for var in m.component_objects(Var, descend_into=True):
        if var.local_name in ("flow_mol_phase", "phase_frac"):
            for v in var.values():
                v.setlb(-1e-8)


def unit_scale(expr):
    """Fixed scale of a variable or constraint body from its physical units; None if unclassified."""
    try:
        u = pyunits.get_units(expr)
    except Exception:
        return None
    u = pyunits.dimensionless if u is None else u
    for ref, s in UNIT_SCALES:
        try:
            if abs(pyunits.convert_value(1.0, from_units=u, to_units=ref) - 1.0) < 1e-12:
                return s
        except Exception:
            continue
    return None


def set_unit_scaling(m):
    """Scaling factors 1/scale for every variable and constraint from its units, and 1e-5 for the
    objective: the fixed scales of the issue's KKT evidence, used by the scaled NL writer."""
    for obj in (*m.component_data_objects(Var, descend_into=True),
                *m.component_data_objects(Constraint, descend_into=True)):
        s = unit_scale(obj if obj.ctype is Var else obj.body)
        if s is not None:
            iscale.set_scaling_factor(obj, 1.0 / s)
    iscale.set_scaling_factor(m.fs.objective, 1.0 / OBJECTIVE_SCALE)


def decisions(m):
    fs = m.fs
    return (fs.H101.outlet.temperature[0], fs.F101.vap_outlet.temperature[0],
            fs.F102.vap_outlet.temperature[0], fs.F102.vap_outlet.pressure[0])


def fix_simulation(m, start=(600, 325, 375, 150000)):
    """Square simulation: tutorial specifications, adiabatic reactor, fixed conversion."""
    fs = m.fs
    fs.R101.control_volume.conversion.fix(0.75)
    fs.R101.heat_duty.fix(0)
    fs.F101.deltaP.fix(0)
    fs.S101.split_fraction[0, "purge"].fix(0.2)
    fs.C101.outlet.pressure.fix(350000)
    fs.F102.deltaP.unfix()
    for var, v in zip(decisions(m), start):
        var.fix(v)
    for c in (fs.overhead_loss, fs.product_flow, fs.product_purity, fs.objective):
        c.deactivate()


def set_optimization(m):
    fs = m.fs
    for var, (lb, ub) in zip(decisions(m), ((500, 600), (298, 450), (298, 450), (105000, 110000))):
        var.unfix()
        var.setlb(lb)
        var.setub(ub)
    fs.F102.deltaP.unfix()
    for c in (fs.overhead_loss, fs.product_flow, fs.product_purity, fs.objective):
        c.activate()
    assert degrees_of_freedom(m) == 4


def set_smoothing(m, factor):
    """Scale every SmoothVLE width (eps_1 for T_bub, eps_2 for T_dew) from the IDAES default."""
    for prm in m.component_objects(Param, descend_into=True):  # scalar Params holding a default
        for eps, default in SMOOTHING_DEFAULTS.items():
            if prm.local_name.startswith(eps + "_"):
                prm.set_value(default * factor)


def initialize_ideal(m):
    tear_guesses = fix_inlet_states(m)
    fix_simulation(m)
    manual_propagation(m, tear_guesses)


def initialize_from(m, ref):
    """Start an ePC-SAFT flowsheet from a solved flowsheet `ref` with the same state definition
    (the ideal baseline): every same-named variable takes the reference value (states, phase
    splits, dew points, _teq, unit variables), and the feeds are fixed as in the tutorial."""
    for v in m.component_data_objects(Var, descend_into=True):
        r = ref.find_component(v.name)
        if r is not None and not v.fixed and r.value is not None:
            v.set_value(r.value)
    for arc in ("s06", "s07", "s09"):  # translators (no ideal counterpart) take both neighbours' states
        propagate_state(arc=getattr(m.fs, arc + "a"))
        propagate_state(arc=getattr(m.fs, arc + "b"), direction="backward")
    for sb in (b for b in m.component_data_objects(Block, descend_into=True) if isinstance(b, StateBlockData)):
        for v, c in [*((sb.mole_frac_comp[j], sb.mole_frac_comp_eq[j]) for j in sb.component_list if not sb.mole_frac_comp[j].fixed),
                     *((sb.phase_frac[p], sb.phase_fraction_constraint[p]) for p in sb.phase_list)]:
            calculate_variable_from_constraint(v, c)
    _init_Tbub(m.fs.F102.control_volume.properties_out[0], pyunits.K)  # F102's bubble operation: IDAES's seed
    SmoothVLE.calculate_teq(m.fs.F102.control_volume.properties_out[0], PP)
    fix_inlet_states(m)


def solve(m, tee=False, log=None):
    """IDAES ipopt_v2 settings (scaled NL writer, linear presolve, MA57); returns solver evidence.
    `log` keeps Ipopt's output file <log>.ipopt and its stdout/stderr (callback failures) <log>.log."""
    options = SOLVER_OPTIONS if m.package == "ideal" else EPCSAFT_OPTIONS
    if log:
        options = {**options, "output_file": f"{log}.ipopt", "file_print_level": 5, "print_user_options": "yes"}
    t0 = time.perf_counter()
    res = Ipopt().solve(m, tee=tee, load_solutions=False, raise_exception_on_nonoptimal_result=False,
                        solver_options=options, writer_config={"scale_model": True, "linear_presolve": True})
    wall = time.perf_counter() - t0
    if res.solution_status != SolutionStatus.noSolution:
        res.solution_loader.load_vars()
    out = {"termination": res.termination_condition.name, "iterations": res.extra_info.iteration_count, "wall_s": wall}
    if log:
        Path(f"{log}.log").write_text(res.solver_log)
        out.update(exit=Path(f"{log}.ipopt").read_text().partition("EXIT: ")[2].split("\n")[0].strip(),
                   command_line=res.extra_info.command_line, return_code=res.extra_info.return_code)
    return out


def summary(m):
    """Key flowsheet numbers of the current solution."""
    fs = m.fs
    vf = {u: value(b.flow_mol_phase["Vap"] / sum(b.flow_mol_phase[p] for p in b.phase_list))
          for u, b in (("F101", fs.F101.control_volume.properties_out[0]),
                       ("F102", fs.F102.control_volume.properties_out[0]))}
    return {"objective_USD_per_yr": value(fs.operating_cost),
            "decisions": [value(v) for v in decisions(m)],
            "duties_W": {u: value(getattr(fs, u).heat_duty[0]) for u in ("H101", "R101", "F101", "F102")},
            "C101_work_W": value(fs.C101.work_mechanical[0]),
            "purity": value(fs.purity), "recovery": value(fs.recovery),
            "reactor_outlet_T_K": value(fs.R101.outlet.temperature[0]), "vapor_fraction": vf}


def run_campaign(m, base_state):
    """Fixed simulation, then the optimization from each fixed start (square solve at the start first)."""
    from idaes.core.util import from_json
    out = {"simulation": {**solve(m), **summary(m)}, "optimizations": []}
    for start in STARTS:
        from_json(m, sd=base_state)
        fix_simulation(m, start)
        square = solve(m)
        set_optimization(m)
        out["optimizations"].append({"start": start, "start_square_solve": square["termination"],
                                     **solve(m), **summary(m)})
    return out


def run_ideal_baseline(factor=1.0):
    from idaes.core.util import to_json
    m = build("ideal")
    set_smoothing(m, factor)
    initialize_ideal(m)
    solve(m)
    return m, run_campaign(m, to_json(m, return_dict=True))
