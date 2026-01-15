from pyomo.environ import (
    Constraint,
    Var,
    ConcreteModel,
    Expression,
    TransformationFactory,
    value,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock

from idaes.models.unit_models import (
    PressureChanger,
    Mixer,
    Separator as Splitter,
    Heater,
    StoichiometricReactor,
    Feed,
    Product,
)

from idaes.models.unit_models import Flash

from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.core.solvers import get_solver

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
)
from idaes_examples.mod.hda.hda_ideal_VLE_modular import thermo_config_fxn
from idaes_examples.mod.hda.hda_reaction_modular import reaction_config

from idaes.models.properties.modular_properties.state_definitions.FpTPxpc import FpTPxpc
from idaes.models.properties.modular_properties.state_definitions.FpcTP import FpcTP


from idaes_examples.mod.hda.hda_flowsheet_extras import (manual_propagation, automatic_propagation,
                                                         fix_inlet_states_1, fix_inlet_states_2,
                                                         FpcTP_to_FpTPxpc)

if __name__ == "__main__":

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    state_definition = FpTPxpc
    # state_definition = FpcTP
    thermo_config = thermo_config_fxn(state_definition)

    m.fs.thermo_params = GenericParameterBlock(**thermo_config)

    m.fs.reaction_params = GenericReactionParameterBlock(
        property_package=m.fs.thermo_params, **reaction_config
    )

    m.fs.I101 = Feed(property_package=m.fs.thermo_params)
    m.fs.I102 = Feed(property_package=m.fs.thermo_params)

    m.fs.M101 = Mixer(
        property_package=m.fs.thermo_params,
        num_inlets=3,
    )

    m.fs.H101 = Heater(
        property_package=m.fs.thermo_params,
        has_pressure_change=False,
        has_phase_equilibrium=True,
    )

    m.fs.R101 = StoichiometricReactor(
        property_package=m.fs.thermo_params,
        reaction_package=m.fs.reaction_params,
        has_heat_of_reaction=True,
        has_heat_transfer=True,
        has_pressure_change=False,
    )

    m.fs.F101 = Flash(
        property_package=m.fs.thermo_params,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    m.fs.S101 = Splitter(
        property_package=m.fs.thermo_params,
        ideal_separation=False,
        outlet_list=["purge", "recycle"],
    )

    m.fs.C101 = PressureChanger(
        property_package=m.fs.thermo_params,
        compressor=True,
        thermodynamic_assumption=ThermodynamicAssumption.isothermal,
    )

    m.fs.F102 = Flash(
        property_package=m.fs.thermo_params,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    m.fs.P101 = Product(property_package=m.fs.thermo_params)
    m.fs.P102 = Product(property_package=m.fs.thermo_params)
    m.fs.P103 = Product(property_package=m.fs.thermo_params)

    m.fs.s01 = Arc(source=m.fs.I101.outlet, destination=m.fs.M101.inlet_1)
    m.fs.s02 = Arc(source=m.fs.I102.outlet, destination=m.fs.M101.inlet_2)
    m.fs.s03 = Arc(source=m.fs.M101.outlet, destination=m.fs.H101.inlet)
    m.fs.s04 = Arc(source=m.fs.H101.outlet, destination=m.fs.R101.inlet)
    m.fs.s05 = Arc(source=m.fs.R101.outlet, destination=m.fs.F101.inlet)
    m.fs.s06 = Arc(source=m.fs.F101.vap_outlet, destination=m.fs.S101.inlet)
    m.fs.s07 = Arc(source=m.fs.F101.liq_outlet, destination=m.fs.F102.inlet)
    m.fs.s08 = Arc(source=m.fs.S101.recycle, destination=m.fs.C101.inlet)
    m.fs.s09 = Arc(source=m.fs.C101.outlet, destination=m.fs.M101.inlet_3)
    m.fs.s10 = Arc(source=m.fs.F102.vap_outlet, destination=m.fs.P101.inlet)
    m.fs.s11 = Arc(source=m.fs.F102.liq_outlet, destination=m.fs.P102.inlet)
    m.fs.s12 = Arc(source=m.fs.S101.purge, destination=m.fs.P103.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    print(degrees_of_freedom(m))

    eps = 1e-5
    # eps = 0.0
    flow_mol_phase_comp_101 = {
        ("Vap", "benzene"): eps,
        ("Vap", "toluene"): eps,
        ("Vap", "hydrogen"): eps,
        ("Vap", "methane"): eps,
        ("Liq", "benzene"): eps,
        ("Liq", "toluene"): 0.30,
    }

    flow_mol_phase_101, mole_frac_phase_comp_101 = FpcTP_to_FpTPxpc(flow_mol_phase_comp_101)

    m.fs.I101.flow_mol_phase[0, "Vap"].fix(flow_mol_phase_101["Vap"])
    m.fs.I101.flow_mol_phase[0, "Liq"].fix(flow_mol_phase_101["Liq"])
    m.fs.I101.mole_frac_phase_comp[0, "Vap", "benzene"].fix(mole_frac_phase_comp_101["Vap", "benzene"])
    m.fs.I101.mole_frac_phase_comp[0, "Vap", "toluene"].fix(mole_frac_phase_comp_101["Vap", "toluene"])
    m.fs.I101.mole_frac_phase_comp[0, "Vap", "hydrogen"].fix(mole_frac_phase_comp_101["Vap", "hydrogen"])
    m.fs.I101.mole_frac_phase_comp[0, "Vap", "methane"].fix(mole_frac_phase_comp_101["Vap", "methane"])
    m.fs.I101.mole_frac_phase_comp[0, "Liq", "benzene"].fix(mole_frac_phase_comp_101["Liq", "benzene"])
    m.fs.I101.mole_frac_phase_comp[0, "Liq", "toluene"].fix(mole_frac_phase_comp_101["Liq", "toluene"])
    m.fs.I101.temperature.fix(303.2)
    m.fs.I101.pressure.fix(350000)

    flow_mol_phase_comp_102 = {
        ("Vap", "benzene"): eps,
        ("Vap", "toluene"): eps,
        ("Vap", "hydrogen"): .30,
        ("Vap", "methane"): .02,
        ("Liq", "benzene"): eps,
        ("Liq", "toluene"): eps,
    }

    flow_mol_phase_102, mole_frac_phase_comp_102 = FpcTP_to_FpTPxpc(flow_mol_phase_comp_102)

    m.fs.I102.flow_mol_phase[0, "Vap"].fix(flow_mol_phase_102["Vap"])
    m.fs.I102.flow_mol_phase[0, "Liq"].fix(flow_mol_phase_102["Liq"])
    m.fs.I102.mole_frac_phase_comp[0, "Vap", "benzene"].fix(mole_frac_phase_comp_102["Vap", "benzene"])
    m.fs.I102.mole_frac_phase_comp[0, "Vap", "toluene"].fix(mole_frac_phase_comp_102["Vap", "toluene"])
    m.fs.I102.mole_frac_phase_comp[0, "Vap", "hydrogen"].fix(mole_frac_phase_comp_102["Vap", "hydrogen"])
    m.fs.I102.mole_frac_phase_comp[0, "Vap", "methane"].fix(mole_frac_phase_comp_102["Vap", "methane"])
    m.fs.I102.mole_frac_phase_comp[0, "Liq", "benzene"].fix(mole_frac_phase_comp_102["Liq", "benzene"])
    m.fs.I102.mole_frac_phase_comp[0, "Liq", "toluene"].fix(mole_frac_phase_comp_102["Liq", "toluene"])
    m.fs.I102.temperature.fix(303.2)
    m.fs.I102.pressure.fix(350000)

    tear_guesses = {
        "flow_mol_phase": {
            (0, "Liq"): flow_mol_phase_101["Liq"],
            (0, "Vap"): flow_mol_phase_102["Vap"],
        },
        "mole_frac_phase_comp": {
            (0, "Liq", "benzene"): mole_frac_phase_comp_101["Liq", "benzene"],
            (0, "Liq", "toluene"): mole_frac_phase_comp_101["Liq", "toluene"],
            (0, "Vap", "benzene"): mole_frac_phase_comp_102["Vap", "benzene"],
            (0, "Vap", "toluene"): mole_frac_phase_comp_102["Vap", "toluene"],
            (0, "Vap", "methane"): mole_frac_phase_comp_102["Vap", "hydrogen"],
            (0, "Vap", "hydrogen"): mole_frac_phase_comp_102["Vap", "methane"],
        },
        "temperature": {0: 303},
        "pressure": {0: 350000},
    }

    print('Heater degrees of freedom: ', degrees_of_freedom(m.fs.H101.control_volume) - 10)

    m.fs.H101.outlet.temperature.fix(600)

    print('Heater degrees of freedom: ', degrees_of_freedom(m.fs.H101.control_volume) - 10)

    # m.fs.benzene_purity = Var(initialize=0.75, bounds=(0, 1))
    #
    # m.fs.benzene_purity_constraint = Constraint(
    #     expr=m.fs.benzene_purity ==
    #     (m.fs.F102.control_volume.properties_out[0].flow_mol_phase_comp["Vap", "benzene"] /
    #     (m.fs.F102.control_volume.properties_out[0].flow_mol_phase_comp["Vap", "benzene"] +
    #      m.fs.F102.control_volume.properties_out[0].flow_mol_phase_comp["Vap", "toluene"])))

    # m.fs.benzene_purity.fix(0.80)

    print('Reactor degrees of freedom: ', degrees_of_freedom(m.fs.R101.control_volume) - 10)

    # m.fs.R101.rate_reaction_extent[0, 'hydrodealkylation'].fix(8.457291523913673)

    # m.fs.R101.outlet.temperature.setlb(500)
    # m.fs.R101.outlet.temperature.setub(1000)
    # m.fs.R101.outlet.temperature.fix(725.0)

    # m.fs.benzene_purity = Expression(
    #     expr=m.fs.F102.control_volume.properties_out[0].flow_mol_phase_comp[
    #         "Vap", "benzene"
    #     ]
    #     / (
    #         m.fs.F102.control_volume.properties_out[0].flow_mol_phase_comp["Vap", "benzene"]
    #         + m.fs.F102.control_volume.properties_out[0].flow_mol_phase_comp[
    #             "Vap", "toluene"
    #         ]
    #     )
    # )
    #
    m.fs.R101.control_volume.conversion = Var(initialize=0.75, bounds=(0, 1))

    m.fs.R101.conv_constraint = Constraint(
        expr=m.fs.R101.control_volume.conversion
        * (m.fs.R101.control_volume.properties_in[0].flow_mol_phase_comp["Vap", "toluene"])
        == (
            m.fs.R101.control_volume.properties_in[0].flow_mol_phase_comp["Vap", "toluene"] - m.fs.R101.control_volume.properties_out[0].flow_mol_phase_comp["Vap", "toluene"])
    )
    m.fs.R101.control_volume.conversion.fix(.75)

    m.fs.R101.heat_duty.fix(0.0)
    print('Reactor degrees of freedom: ', degrees_of_freedom(m.fs.R101) - 10)

    print('Flash 1 degrees of freedom: ', degrees_of_freedom(m.fs.F101) - 10)
    m.fs.F101.vap_outlet.temperature.fix(325.0)
    m.fs.F101.deltaP.fix(0.0)
    print('Flash 1 degrees of freedom: ', degrees_of_freedom(m.fs.F101) - 10)

    print('Flash 2 degrees of freedom: ', degrees_of_freedom(m.fs.F102) - 10)
    m.fs.F102.vap_outlet.temperature.fix(375)
    m.fs.F102.deltaP.fix(-200000)
    print('Flash 2 degrees of freedom: ', degrees_of_freedom(m.fs.F102) - 10)

    print('Splitter degrees of freedom: ', degrees_of_freedom(m.fs.S101) - 10)
    m.fs.S101.split_fraction[0, "purge"].fix(0.2)
    print('Splitter degrees of freedom: ', degrees_of_freedom(m.fs.S101) - 10)
    print('Cooler degrees of freedom: ', degrees_of_freedom(m.fs.C101) - 10)
    m.fs.C101.outlet.pressure.fix(350000)
    print('Cooler degrees of freedom: ', degrees_of_freedom(m.fs.C101) - 10)

    m.fs.cooling_cost = Expression(
        expr=0.212e-7 * (-m.fs.F101.heat_duty[0]) + 0.212e-7 * (-m.fs.R101.heat_duty[0])
    )

    m.fs.heating_cost = Expression(
        expr=2.2e-7 * m.fs.H101.heat_duty[0] + 1.9e-7 * m.fs.F102.heat_duty[0]
    )

    m.fs.operating_cost = Expression(
        expr=(3600 * 24 * 365 * (m.fs.heating_cost + m.fs.cooling_cost))
    )

    print(degrees_of_freedom(m))

    # automatic_propagation(m, tear_guesses)
    manual_propagation(m, tear_guesses)

    optarg = {
        "nlp_scaling_method": "user-scaling",
        "OF_ma57_automatic_scaling": "yes",
        "max_iter": 1000,
        "tol": 1e-8,
    }

    solver = get_solver("ipopt_v2", options=optarg)

    # Solve the model
    results = solver.solve(m, tee=True)

    from pyomo.environ import TerminationCondition

    assert results.solver.termination_condition == TerminationCondition.optimal

    print("benzene purity = ", value(m.fs.benzene_purity))
    print("operating cost = $", value(m.fs.operating_cost))
    print("reactor temperature = $", value(m.fs.R101.outlet.temperature[0]))

    # assert value(m.fs.purity) == pytest.approx(0.82429, abs=1e-3)
