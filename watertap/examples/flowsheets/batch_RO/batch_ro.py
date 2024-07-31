import os
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    value,
    Var,
    NonNegativeReals,
    assert_optimal_termination,
    Objective,
    units as pyunits,
)
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor
from watertap.unit_models.reverse_osmosis_1D import (ReverseOsmosis1D, ConcentrationPolarizationType, 
                                                     MassTransferCoefficient, PressureChangeType)
from watertap.unit_models.reverse_osmosis_0D import (ReverseOsmosis0D,
        ConcentrationPolarizationType, MassTransferCoefficient)
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from idaes.core.util.model_statistics import degrees_of_freedom
from watertap.core.solvers import get_solver
from idaes.models.unit_models import Product, Feed
from idaes.models.unit_models.mixer import (
    Mixer,
    MomentumMixingType,
    MaterialBalanceType,
)
from watertap.costing import (
    WaterTAPCosting,
    PumpType,
    MixerType,
    ROType,
)
import idaes.core.util.scaling as iscale
from idaes.core.util.initialization import propagate_state
from pyomo.environ import TransformationFactory
from pyomo.network import Arc
from watertap.unit_models.pressure_changer import Pump

def build_system(time_blk=None):

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()
    m.fs.properties = SeawaterParameterBlock()
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.recirc = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)

    m.fs.P1 = Pump(property_package=m.fs.properties)
    m.fs.P1.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )
    m.fs.P2 = Pump(property_package=m.fs.properties)
    m.fs.P2.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )

    m.fs.M1 = Mixer(
        property_package=m.fs.properties,
        has_holdup=False,
        num_inlets = 2,
        momentum_mixing_type = MomentumMixingType.equality,
    )

    m.fs.RO = ReverseOsmosis0D(
    property_package=m.fs.properties,
    has_pressure_change=True,
    pressure_change_type=PressureChangeType.calculated,
    mass_transfer_coefficient=MassTransferCoefficient.calculated,
    concentration_polarization_type=ConcentrationPolarizationType.calculated,
    has_full_reporting = True,
    )
    m.fs.RO.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    # connect unit models
    m.fs.feed_to_P1 = Arc(source=m.fs.feed.outlet, destination=m.fs.P1.inlet)
    m.fs.recirc_to_P2 = Arc(source=m.fs.recirc.outlet, destination=m.fs.P2.inlet)
    m.fs.P1_to_M1 = Arc(source=m.fs.P1.outlet, destination=m.fs.M1.inlet_1)
    m.fs.P2_to_M1 = Arc(source=m.fs.P2.outlet, destination=m.fs.M1.inlet_2)
    m.fs.M1_to_RO = Arc(source=m.fs.M1.outlet, destination=m.fs.RO.inlet)

    # m.fs.feed_to_P1 = Arc(source=m.fs.feed.outlet, destination=m.fs.P1.inlet)
    # m.fs.P1_to_M1 = Arc(source=m.fs.P1.outlet, destination=m.fs.M1.inlet_1)
    # m.fs.M1_to_RO = Arc(source=m.fs.M1.outlet, destination=m.fs.RO.inlet)

    m.fs.RO_permeate_to_product = Arc(source=m.fs.RO.permeate, destination=m.fs.product.inlet)
    m.fs.RO_retentate_to_disposal = Arc(source=m.fs.RO.retentate, destination=m.fs.disposal.inlet)
    # m.fs.RO_to_P2 = Arc(source=m.fs.RO.retentate, destination=m.fs.P2.inlet)
    # m.fs.P2_to_M1 = Arc(source=m.fs.P2.outlet, destination=m.fs.M1.inlet_2)
    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1e2, index=("Liq", "TDS"))

    # set unit model values
    set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    set_scaling_factor(m.fs.P2.control_volume.work, 1e-3)
    set_scaling_factor(m.fs.RO.area, 1e-2)

    # calculate and propagate scaling factors
    calculate_scaling_factors(m)

    m.fs.feed.properties[0].flow_vol_phase["Liq"]
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"]
    m.fs.recirc.properties[0].flow_vol_phase["Liq"]
    m.fs.recirc.properties[0].mass_frac_phase_comp["Liq", "TDS"]
    m.fs.M1.mixed_state[0].flow_vol_phase["Liq"]
    m.fs.M1.mixed_state[0].conc_mass_phase_comp["Liq", "TDS"]

    return m


def set_operating_conditions(m, Qin = 1, Cin = 35):

    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )

    m.fs.feed_salinity = Var(
        initialize=35,
        bounds=(0, 2000),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )

    m.fs.feed_flow_mass = Var(
        initialize=1,
        bounds=(0.00001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.kg / pyunits.s,
        doc="System Feed Flowrate",
    )
    
    # m.fs.feed_salinity = Var(
    #     initialize=35,
    #     bounds=(0, 2000),
    #     domain=NonNegativeReals,
    #     units=pyunits.dimensionless,
    #     doc="System Water Recovery",
    # )

    # m.fs.feed_flow_mass = Var(
    #     initialize=1,
    #     bounds=(0.00001, 1e6),
    #     domain=NonNegativeReals,
    #     units=pyunits.kg / pyunits.s,
    #     doc="System Feed Flowrate",
    # )

    iscale.set_scaling_factor(m.fs.water_recovery, 10)
    iscale.set_scaling_factor(m.fs.feed_flow_mass, 1)
    iscale.set_scaling_factor(m.fs.feed_salinity, 1)

    if Qin is not None:
        m.fs.feed_flow_mass.fix(Qin)
    if Cin is not None:
        m.fs.feed_salinity.fix(Cin)

    m.fs.feed_flow_constraint = Constraint(
        expr= m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'H2O']
        == m.fs.feed_flow_mass
    )

    m.fs.nacl_mass_constraint = Constraint(
        expr=m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"] * 1000
        == m.fs.feed_flow_mass * m.fs.feed_salinity
    )

    m.fs.feed.properties[0].pressure.fix(101325)                           # pressure (Pa)
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)                   # temperature (K)

    m.fs.recirc.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(0.1)
    m.fs.recirc.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(0.010275)
    m.fs.recirc.properties[0].pressure.fix(101325)                           # pressure (Pa)
    m.fs.recirc.properties[0].temperature.fix(273.15 + 25)                   # temperature (K)

    # high pressure pump, 2 degrees of freedom
    m.fs.P1.efficiency_pump.fix(0.80)                                    # pump efficiency (-)
    m.fs.P1.control_volume.properties_out[0].pressure.fix(75e5)          # pump outlet pressure (Pa)
    m.fs.P2.efficiency_pump.fix(0.80)                                    # pump efficiency (-)
    # m.fs.P2.control_volume.properties_out[0].pressure.fix(75e5)          # pump outlet pressure (Pa)

    # RO unit, 7 degrees of freedom
    m.fs.RO.A_comp.fix(4.2e-12)                                            # membrane water permeability coeff (m/Pa/s)
    m.fs.RO.B_comp.fix(3.5e-8)                                             # membrane salt permeability coeff (m/s)

    # fix 4 module specficiations
    m.fs.RO.area.fix(25)                                                # membrane stage area (m^2)
    m.fs.RO.length.fix(7)                                                # membrane stage width (m)
    m.fs.RO.feed_side.channel_height.fix(1E-3)                          # channel height in membrane stage (m)
    m.fs.RO.feed_side.spacer_porosity.fix(0.97)                         # spacer porosity in membrane stage (-)
    m.fs.RO.permeate.pressure[0].fix(101325)                               # permeate pressure (Pa)

    print("DOF = ", degrees_of_freedom(m))
    assert_no_degrees_of_freedom(m)


def set_operating_constraints(m):

    m.fs.eq_feed_product_balance = Constraint(
        expr=m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
        == m.fs.product.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )

    m.fs.eq_pump_equality = Constraint(
        expr=m.fs.P1.control_volume.properties_out[0].pressure
        == m.fs.P2.control_volume.properties_out[0].pressure
    )


def add_costing(m):

    m.fs.costing.cost_process()
    product_flow_vol_total = m.fs.product.properties[0].flow_vol
    m.fs.costing.add_annual_water_production(product_flow_vol_total)
    m.fs.costing.add_specific_energy_consumption(product_flow_vol_total)
    m.fs.costing.add_LCOW(product_flow_vol_total)
    # m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOW)


def initialize_mixer(m, guess = True):
    if guess:
        recovery_guess = 0.3
        m.fs.M1.inlet_2_state[0.0].flow_mass_phase_comp["Liq", "H2O"].fix((1-recovery_guess)*0.965)
        m.fs.M1.inlet_2_state[0.0].flow_mass_phase_comp["Liq", "TDS"].fix(0.035/(1-recovery_guess))
    m.fs.M1.initialize()
    m.fs.M1.inlet_2_state[0.0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    m.fs.M1.inlet_2_state[0.0].flow_mass_phase_comp["Liq", "TDS"].unfix()


def do_forward_initialization_pass(m, pass_num=1):
    # initialize unit by unit
    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_P1)

    m.fs.recirc.initialize()
    propagate_state(m.fs.recirc_to_P2)

    m.fs.P1.initialize()
    propagate_state(m.fs.P1_to_M1)

    m.fs.P2.initialize()
    propagate_state(m.fs.P2_to_M1)

    # if pass_num > 0:
    #     m.fs.M1.initialize()
    # else:
    #     initialize_mixer(m)
    m.fs.M1.initialize()
    propagate_state(m.fs.M1_to_RO)

    m.fs.RO.initialize()
    propagate_state(m.fs.RO_permeate_to_product)
    propagate_state(m.fs.RO_retentate_to_disposal)

    m.fs.product.initialize()
    m.fs.disposal.initialize()
    # propagate_state(m.fs.RO_to_P2)

    # m.fs.P2.initialize()
    # propagate_state(m.fs.P2_to_M1)


def initialize(m):
    # initialize unit by unit
    for idx, init_pass in enumerate(range(1)):
        print(f'\n\nINITIALIZATION PASS {idx+1}\n\n')
        do_forward_initialization_pass(m, pass_num=idx)
        print_results(m)


def optimize(m, water_recovery=0.1, Q_ro = 0.965):

    m.fs.feed_flow_mass.unfix()
    m.fs.RO.inlet.flow_mass_phase_comp[0,"Liq", "H2O"].fix(Q_ro)
    m.fs.RO.recovery_vol_phase[0, "Liq"].fix(0.1)
    m.fs.P1.control_volume.properties_out[0].pressure.unfix()          # pump outlet pressure (Pa)                                # pump efficiency (-)
    m.fs.P2.control_volume.properties_out[0].pressure.unfix()          # pump outlet pressure (Pa)


def solve(m, solver=None, tee=True, raise_on_failure=True):
    if solver is None:
        solver = get_solver()
        solver.options["max_iter"] = 2000

    print("\n--------- SOLVING ---------\n")

    results = solver.solve(m, tee=tee)
    # store new state so we can see what was chnged before solve
    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        print("\n--------- INFEASIBLE SOLVE!!! ---------\n")

        print("\n--------- CLOSE TO BOUNDS ---------\n")
        print_close_to_bounds(m)

        print("\n--------- INFEASIBLE BOUNDS ---------\n")
        print_infeasible_bounds(m)

        print("\n--------- INFEASIBLE CONSTRAINTS ---------\n")
        print_infeasible_constraints(m)

        print_results(m)


def print_results(m):
    print('\n\n')
    print(f'MIXER INLET 1: {value(m.fs.M1.inlet_1_state[0].flow_mass_phase_comp["Liq", "H2O"]):<5.2f}')
    print(f'MIXER INLET 2: {value(m.fs.M1.inlet_2_state[0].flow_mass_phase_comp["Liq", "H2O"]):<5.2f}')
    print(f'MIXER OUTLET: {value(m.fs.M1.mixed_state[0].flow_mass_phase_comp["Liq", "H2O"]):<5.2f}')
    print(f'MIXER CONC: {value(pyunits.convert(m.fs.M1.mixed_state[0].conc_mass_phase_comp["Liq", "TDS"], to_units=pyunits.g/pyunits.L)):<5.2f} {pyunits.get_units(pyunits.convert(m.fs.M1.mixed_state[0].conc_mass_phase_comp["Liq", "TDS"], to_units=pyunits.g/pyunits.L))}')
    print('\n')
    print(f'PUMP 1 INLET: {value(m.fs.P1.control_volume.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"]):<5.2f}')
    print(f'PUMP 1 OUTLET: {value(m.fs.P1.control_volume.properties_out[0.0].flow_mass_phase_comp["Liq", "H2O"]):<5.2f}')
    print(f'PUMP 1 PRESSURE: {value(pyunits.convert(m.fs.P1.control_volume.properties_out[0.0].pressure, to_units=pyunits.bar)):<5.2f}')
    print('\n')
    print(f'PUMP 2 INLET: {value(m.fs.P2.control_volume.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"]):<5.2f}')
    print(f'PUMP 2 OUTLET: {value(m.fs.P2.control_volume.properties_out[0.0].flow_mass_phase_comp["Liq", "H2O"]):<5.2f}')
    print(f'PUMP 2 PRESSURE: {value(pyunits.convert(m.fs.P2.control_volume.properties_out[0.0].pressure, to_units=pyunits.bar)):<5.2f}')
    print('\n')
    print(f'RO FEED: {value(m.fs.RO.inlet.flow_mass_phase_comp[0,"Liq", "H2O"]):<5.2f}')
    print(f'RO PRODUCT: {value(m.fs.RO.permeate.flow_mass_phase_comp[0,"Liq", "H2O"]):<5.2f}')
    print(f'RO BRINE: {value(m.fs.RO.retentate.flow_mass_phase_comp[0,"Liq", "H2O"]):<5.2f}')
    print('\n')
    print(f'LCOW: {value(m.fs.costing.LCOW):<5.2f} {pyunits.get_units(m.fs.costing.LCOW)}')
    print('\n\n')
    print(m.fs.M1.report())
    print(m.fs.RO.report())


def main():
    m = build_system()
    set_operating_conditions(m)
    add_costing(m)
    initialize(m)
    optimize(m)
    solve(m)
    print_results(m)
    

if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    main()

    #BUG the TDS is going up a ton!
    #BUG Need to move the pump before the mixer