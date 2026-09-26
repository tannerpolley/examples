"""HDA two-flash recycle flowsheet: ideal baseline and ePC-SAFT properties (ePC-SAFT issue #73).

Adiabatic stoichiometric reactor at fixed 75 % toluene conversion. Decisions:
H101 outlet T, F101 vapor T, F102 vapor T, F102 outlet P; tutorial operating-cost
objective with overhead-loss, benzene-production and purity constraints.
"""

import time

from pyomo.environ import (ConcreteModel, Constraint, Expression, Objective, TransformationFactory,
                           Var, value)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from pyomo.contrib.solver.common.results import SolutionStatus
from pyomo.contrib.solver.solvers.ipopt import Ipopt
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock)
from idaes.models.unit_models import (Feed, Flash, Heater, Mixer, PressureChanger, Product,
                                      Separator, StoichiometricReactor)
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption

from idaes_examples.mod.hda import hda_epcsaft_VLE, hda_ideal_VLE_modular
from idaes_examples.mod.hda.hda_flowsheet_extras import fix_inlet_states, manual_propagation
from idaes_examples.mod.hda.hda_reaction_modular import reaction_config

COMPONENTS = ("hydrogen", "methane", "benzene", "toluene")
STARTS = ((550, 325, 375, 107500), (500, 300, 350, 105000), (600, 400, 425, 110000))
SOLVER_OPTIONS = {"nlp_scaling_method": "user-scaling", "ma57_automatic_scaling": "yes",
                  "max_iter": 1000, "tol": 1e-8, "bound_relax_factor": 0.0,
                  "linear_solver": "ma57"}


def build(package="ideal"):
    """The tutorial flowsheet with the selected property package ("ideal" | "epcsaft")."""
    m = ConcreteModel()
    m.package = package
    m.fs = FlowsheetBlock(dynamic=False)
    cfg = (hda_ideal_VLE_modular if package == "ideal" else hda_epcsaft_VLE).thermo_config
    fs = m.fs
    fs.thermo_params = props = GenericParameterBlock(**cfg)
    fs.reaction_params = GenericReactionParameterBlock(property_package=props, **reaction_config)
    fs.I101, fs.I102 = Feed(property_package=props), Feed(property_package=props)
    fs.M101 = Mixer(property_package=props, num_inlets=3)
    fs.H101 = Heater(property_package=props, has_pressure_change=False, has_phase_equilibrium=True)
    # Formation-based ePC-SAFT enthalpies already carry the reaction enthalpy.
    fs.R101 = StoichiometricReactor(property_package=props, reaction_package=fs.reaction_params,
                                    has_heat_of_reaction=package == "ideal", has_heat_transfer=True,
                                    has_pressure_change=False)
    fs.F101 = Flash(property_package=props, has_heat_transfer=True, has_pressure_change=True)
    fs.S101 = Separator(property_package=props, ideal_separation=False,
                        outlet_list=["purge", "recycle"])
    fs.C101 = PressureChanger(property_package=props, compressor=True,
                              thermodynamic_assumption=ThermodynamicAssumption.isothermal)
    fs.F102 = Flash(property_package=props, has_heat_transfer=True, has_pressure_change=True)
    fs.P101, fs.P102, fs.P103 = (Product(property_package=props) for _ in range(3))
    for name, (src, dst) in {
        "s01": (fs.I101.outlet, fs.M101.inlet_1), "s02": (fs.I102.outlet, fs.M101.inlet_2),
        "s03": (fs.M101.outlet, fs.H101.inlet), "s04": (fs.H101.outlet, fs.R101.inlet),
        "s05": (fs.R101.outlet, fs.F101.inlet), "s06": (fs.F101.vap_outlet, fs.S101.inlet),
        "s07": (fs.F101.liq_outlet, fs.F102.inlet), "s08": (fs.S101.recycle, fs.C101.inlet),
        "s09": (fs.C101.outlet, fs.M101.inlet_3), "s10": (fs.F102.vap_outlet, fs.P101.inlet),
        "s11": (fs.F102.liq_outlet, fs.P102.inlet), "s12": (fs.S101.purge, fs.P103.inlet),
    }.items():
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
    return m


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


def initialize_ideal(m):
    tear_guesses = fix_inlet_states(m)
    fix_simulation(m)
    manual_propagation(m, tear_guesses)


def solve(m, tee=False):
    """IDAES ipopt_v2 settings (scaled NL writer, linear presolve, MA57); returns solver evidence."""
    t0 = time.perf_counter()
    res = Ipopt().solve(m, tee=tee, load_solutions=False, raise_exception_on_nonoptimal_result=False,
                        solver_options=SOLVER_OPTIONS,
                        writer_config={"scale_model": True, "linear_presolve": True})
    wall = time.perf_counter() - t0
    if res.solution_status != SolutionStatus.noSolution:
        res.solution_loader.load_vars()
    return {"termination": res.termination_condition.name, "iterations": res.extra_info.iteration_count,
            "wall_s": wall}


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


def run_ideal_baseline():
    from idaes.core.util import to_json
    m = build("ideal")
    initialize_ideal(m)
    solve(m)
    return m, run_campaign(m, to_json(m, return_dict=True))
