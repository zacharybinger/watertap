import os
import math
from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Constraint,
    Objective,
    NonNegativeReals,
    TransformationFactory,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import FlowsheetBlock
import idaes.logger as idaeslogger
from idaes.core.util.initialization import propagate_state as _prop_state
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import InitializationError
from watertap.core.wt_database import Database

from watertap.property_models.NaCl_prop_pack import NaClParameterBlock

from idaes.models.unit_models import Product, Feed, Separator
from idaes.core.util.model_statistics import *
from watertap.core.util.initialization import *
from watertap.costing import WaterTAPCosting
from watertap.unit_models.pressure_changer import Pump
from idaes.core import FlowsheetBlock, UnitModelCostingBlock

from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
    calculate_operating_pressure,
)

from idaes.models.unit_models.mixer import (
    Mixer,
)

from ro_system import (
    build_ro, 
    display_ro_system_build, 
    init_ro_system, 
    init_ro_stage, 
    calc_scale,
    set_ro_system_operating_conditions)

def propagate_state(arc, debug=False):
    _prop_state(arc)
    if debug:
        print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
        arc.source.display()
        print(arc.destination.name)
        arc.destination.display()
        print('\n')

def main():
    file_dir = os.path.dirname(os.path.abspath(__file__))

    m = build_system()

    # # Connect units and add system-level constraints
    add_connections(m)
    add_constraints(m)
    add_costing(m)
    release_constraints_for_ultrapure_water(m)

    # Set inlet conditions and operating conditions for each unit
    set_operating_conditions(m)
    
    # Initialize system, ititialization routines for each unit in definition for init_system
    initialize(m)

    # # Solve system and display results
    solve(m, tee=True, raise_on_failure=True)
    display_flow_table(m)
    display_system_metrics(m)

def build_system():
    m = ConcreteModel()
    m.db = Database()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)

    m.fs.primary_mixer = Mixer(
        property_package=m.fs.properties,
        has_holdup=False,
        num_inlets = 2,
    )
    
    m.fs.primary_pump = Pump(property_package=m.fs.properties)
    m.fs.booster_pump = Pump(property_package=m.fs.properties)

    m.fs.RO1 = FlowsheetBlock(dynamic=False)
    m.fs.RO2 = FlowsheetBlock(dynamic=False)
    
    build_ro(m,m.fs.RO1, number_of_stages=2, ultra_pute_water=True)
    build_ro(m,m.fs.RO2, number_of_stages=2, ultra_pute_water=True)

    m.fs.splitter = Separator(
        property_package=m.fs.properties,
        outlet_list=["recirc", "downstream"],
    )

    m.fs.brine_mixer = Mixer(
        property_package=m.fs.properties,
        has_holdup=False,
        num_inlets = 2,
    )

    return m

def add_connections(m):

    m.fs.feed_to_primary_mixer = Arc(
        source=m.fs.feed.outlet,
        destination=m.fs.primary_mixer.inlet_1,
    )
    
    m.fs.recirc_mixer_to_primary_pump = Arc(
        source=m.fs.primary_mixer.outlet,
        destination=m.fs.primary_pump.inlet,
    )

    m.fs.primary_pump_to_ro = Arc(
        source=m.fs.primary_pump.outlet,
        destination=m.fs.RO1.feed.inlet,
    )

    m.fs.RO1_product_to_booster = Arc(
        source=m.fs.RO1.product.outlet,
        destination=m.fs.booster_pump.inlet,
    )

    m.fs.RO1_disposal_to_brine_mixer = Arc(
        source=m.fs.RO1.disposal.outlet,
        destination=m.fs.brine_mixer.inlet_1,
    )

    m.fs.booster_pump_to_RO2 = Arc(
        source=m.fs.booster_pump.outlet,
        destination=m.fs.RO2.feed.inlet,
    )

    m.fs.ro_to_product = Arc(
        source=m.fs.RO2.product.outlet,
        destination=m.fs.product.inlet,
    )
    m.fs.ro_to_splitter = Arc(
        source=m.fs.RO2.disposal.outlet,
        destination=m.fs.splitter.inlet,
    )

    m.fs.splitter_to_recirc = Arc(
        source=m.fs.splitter.recirc,
        destination=m.fs.primary_mixer.inlet_2,
    )

    m.fs.splitter_to_brine_mixer = Arc(
        source=m.fs.splitter.downstream,
        destination=m.fs.brine_mixer.inlet_2,
    )
    
    m.fs.brine_mixer_to_disposal = Arc(
        source=m.fs.brine_mixer.outlet,
        destination=m.fs.disposal.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

def add_constraints(m):
    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )
    m.fs.pass_one_recovery = Var(
        initialize=0.5,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="Pass One Water Recovery",
    )
    m.fs.pass_two_recovery = Var(
        initialize=0.5,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="Pass Two Water Recovery",
    )

    m.fs.feed_salinity = Var(
        initialize=35,
        bounds=(0, 5000),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="Salinity",
    )

    m.fs.feed_flow_mass = Var(
        initialize=1,
        bounds=(0.000001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.kg / pyunits.s,
        doc="System Feed Flowrate",
    )
    
    m.fs.perm_flow_mass = Var(
        initialize=1,
        bounds=(0.000001, 1e6),
        domain=NonNegativeReals,
        units=pyunits.kg / pyunits.s,
        doc="System Produce Flowrate",
    )

    m.fs.nacl_mass_constraint = Constraint(
        expr=m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"] * 1000
        == m.fs.feed_flow_mass * m.fs.feed_salinity
    )

    m.fs.h2o_mass_constraint = Constraint(
        expr=m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"]
        == m.fs.feed_flow_mass * (1 - m.fs.feed_salinity / 1000)
    )

    m.fs.eq_water_recovery = Constraint(
        expr=m.fs.feed.properties[0].flow_vol * m.fs.water_recovery
        == m.fs.product.properties[0].flow_vol
    )
    m.fs.eq_pass_one_recovery = Constraint(
        expr=m.fs.RO1.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"] * m.fs.pass_one_recovery
        == m.fs.RO1.product.properties[0.0].flow_mass_phase_comp["Liq", "H2O"]
    )
    m.fs.eq_pass_two_recovery = Constraint(
        expr=m.fs.RO2.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"] * m.fs.pass_two_recovery
        == m.fs.RO2.product.properties[0.0].flow_mass_phase_comp["Liq", "H2O"]
    )


    m.fs.product.properties[0].mass_frac_phase_comp
    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.product.properties[0].conc_mass_phase_comp
    m.fs.disposal.properties[0].conc_mass_phase_comp
    m.fs.primary_pump.control_volume.properties_in[0].conc_mass_phase_comp
    m.fs.primary_pump.control_volume.properties_out[0].conc_mass_phase_comp
    m.fs.booster_pump.control_volume.properties_in[0].conc_mass_phase_comp
    m.fs.booster_pump.control_volume.properties_out[0].conc_mass_phase_comp

def add_costing(m):
    m.fs.primary_pump.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.booster_pump.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.primary_mixer.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.RO1.stage[1].module.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"ro_type": "standard"},
    )
    m.fs.RO1.stage[2].module.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"ro_type": "standard"},
    )
    m.fs.RO2.stage[1].module.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"ro_type": "standard"},
    )
    m.fs.RO2.stage[2].module.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"ro_type": "standard"},
    )

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)

    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

def release_constraints_for_ultrapure_water(m):

    # Release constraints for system
    m.fs.RO1.product.properties[0.0].conc_mass_phase_comp["Liq", "NaCl"].setlb(0)
    m.fs.RO2.product.properties[0.0].conc_mass_phase_comp["Liq", "NaCl"].setlb(0)
    m.fs.product.properties[0.0].conc_mass_phase_comp["Liq", "NaCl"].setlb(0)

def set_inlet_conditions(m, Qin=None, Cin=None, water_recovery=None, primary_pump_pressure=10e5, booster_pump_pressure=7e5):
    """Sets operating condition for the PV-RO system

    Args:
        m (obj): Pyomo model
        flow_in (float, optional): feed volumetric flow rate [m3/s]. Defaults to 1e-2.
        conc_in (int, optional): solute concentration [g/L]. Defaults to 30.
        water_recovery (float, optional): water recovery. Defaults to 0.5.
    """
    print(f'\n{"=======> SETTING OPERATING CONDITIONS <=======":^60}\n')
    solver = get_solver()    
    if Qin is None:
        m.fs.feed_flow_mass.fix(1)
    else:
        m.fs.feed_flow_mass.fix(Qin)

    if Cin is None:
        m.fs.feed_salinity.fix(10)
    else:
        m.fs.feed_salinity.fix(Cin)

    if water_recovery is not None:
        m.fs.water_recovery.fix(water_recovery)
        m.fs.primary_pump.control_volume.properties_out[0].pressure.unfix()
    else:
        m.fs.water_recovery.unfix()
        m.fs.primary_pump.control_volume.properties_out[0].pressure.fix(primary_pump_pressure)
        m.fs.booster_pump.control_volume.properties_out[0].pressure.fix(booster_pump_pressure)

    iscale.set_scaling_factor(m.fs.feed_flow_mass, 1)
    iscale.set_scaling_factor(m.fs.feed_salinity, 1)
    iscale.set_scaling_factor(m.fs.booster_pump.control_volume.work, 1)

    feed_temperature = 273.15 + 20
    pressure_atm = 101325
    supply_pressure = 2.7e5

    # initialize feed
    m.fs.feed.pressure[0].fix(supply_pressure)
    m.fs.feed.temperature[0].fix(feed_temperature)

    m.fs.primary_pump.efficiency_pump.fix(0.85)
    m.fs.booster_pump.efficiency_pump.fix(0.85)
    iscale.set_scaling_factor(m.fs.primary_pump.control_volume.work, 1e-3)

    m.fs.feed.properties[0].flow_vol_phase["Liq"]
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"]

    m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value = (
        m.fs.feed_flow_mass.value * m.fs.feed_salinity.value / 1000
    )
    m.fs.feed.flow_mass_phase_comp[
        0, "Liq", "H2O"
    ].value = m.fs.feed_flow_mass.value * (1 - m.fs.feed_salinity.value / 1000)

    scale_flow = calc_scale(m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].value)
    scale_tds = calc_scale(m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value)

    operating_pressure = calculate_operating_pressure(
    feed_state_block=m.fs.feed.properties[0],
    over_pressure=0.15,
    water_recovery=0.8,
    NaCl_passage=0.01,
    solver=solver,
    )

    operating_pressure_psi = pyunits.convert(
        operating_pressure * pyunits.Pa, to_units=pyunits.psi
    )()
    operating_pressure_bar = pyunits.convert(
        operating_pressure * pyunits.Pa, to_units=pyunits.bar
    )()
    print(
        f"\nOperating Pressure Estimate = {round(operating_pressure_bar, 2)} bar = {round(operating_pressure_psi, 2)} psi\n"
    )

    # REVIEW: Make sure this is applied in the right place
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-scale_flow, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-scale_tds, index=("Liq", "NaCl")
    )

    assert_units_consistent(m)

def set_operating_conditions(m):
    # Set inlet conditions and operating conditions for each unit
    set_inlet_conditions(m, Qin=5, Cin=2, primary_pump_pressure=25e5, booster_pump_pressure=10e5)
    RO1_area = 100
    RO2_area = 50
    set_ro_system_operating_conditions(m, m.fs.RO1, mem_area=[200, 125])
    set_ro_system_operating_conditions(m, m.fs.RO2, mem_area=[200, 125])

    iscale.calculate_scaling_factors(m)

    m.fs.splitter.split_fraction[0, "recirc"].fix(0.15)

def initialize(m, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()
    
    optarg = solver.options
    print(f"Degrees of Freedom @ Initialization = {degrees_of_freedom(m)}")
    assert_no_degrees_of_freedom(m)
    for init in range(3):
        do_forward_initialization_pass(m, verbose=True, solver=solver)

def do_forward_initialization_pass(m, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print(
        "\n\n--------------------START FORWARD INITIALIZATION PASS--------------------\n\n"
    )

    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.feed_to_primary_mixer)

    m.fs.primary_mixer.initialize(optarg=optarg)
    propagate_state(m.fs.recirc_mixer_to_primary_pump)

    m.fs.primary_pump.initialize(optarg=optarg)
    propagate_state(m.fs.primary_pump_to_ro)

    init_ro_system(m, m.fs.RO1)
    propagate_state(m.fs.RO1_product_to_booster)
    propagate_state(m.fs.RO1_disposal_to_brine_mixer)

    m.fs.booster_pump.initialize(optarg=optarg)

    propagate_state(m.fs.booster_pump_to_RO2)
    init_ro_system(m, m.fs.RO2)

    propagate_state(m.fs.ro_to_product)
    propagate_state(m.fs.ro_to_splitter)

    m.fs.splitter.initialize(optarg=optarg)
    
    propagate_state(m.fs.splitter_to_recirc)
    propagate_state(m.fs.splitter_to_brine_mixer)
    
    m.fs.brine_mixer.initialize(optarg=optarg)
    propagate_state(m.fs.brine_mixer_to_disposal)

    m.fs.product.initialize(optarg=optarg)
    m.fs.disposal.initialize(optarg=optarg)

def solve(model, solver=None, tee=False, raise_on_failure=False):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")
    print(f"Degrees of Freedom @ Solve = {degrees_of_freedom(model)}")

    results = solver.solve(model, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        raise RuntimeError(msg)
    else:
        return results

def display_flow_table(m):
    print('\n\n')
    print(f'{"NODE":<34s}{"MASS FLOW RATE H2O (KG/S)":<30s}{"PRESSURE (BAR)":<20s}{"MASS FLOW RATE NACL (KG/S)":<30s}{"CONC. (mg/L)":<20s}')
    print(f'{"Feed":<34s}{m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].value:<30.3f}{value(pyunits.convert(m.fs.feed.pressure[0], to_units=pyunits.bar)):<30.1f}{m.fs.feed.flow_mass_phase_comp[0, "Liq", "NaCl"].value:<20.3e}{value(pyunits.convert(m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg * pyunits.L ** -1)):<20.4f}')
    print(f'{"Product":<34s}{m.fs.product.flow_mass_phase_comp[0, "Liq", "H2O"].value:<30.3f}{value(pyunits.convert(m.fs.product.pressure[0], to_units=pyunits.bar)):<30.1f}{m.fs.product.flow_mass_phase_comp[0, "Liq", "NaCl"].value:<20.3e}{value(pyunits.convert(m.fs.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg * pyunits.L ** -1)):<20.4f}')
    print(f'{"Disposal":<34s}{m.fs.disposal.flow_mass_phase_comp[0, "Liq", "H2O"].value:<30.3f}{value(pyunits.convert(m.fs.disposal.pressure[0], to_units=pyunits.bar)):<30.1f}{m.fs.disposal.flow_mass_phase_comp[0, "Liq", "NaCl"].value:<20.3e}{(value(pyunits.convert(m.fs.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg * pyunits.L ** -1))):<20.4f}')

    print(f'{"Primary Pump Inlet":<34s}{m.fs.primary_pump.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.primary_pump.control_volume.properties_in[0].pressure, to_units=pyunits.bar)():<30.1f}{m.fs.primary_pump.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{value(pyunits.convert(m.fs.primary_pump.control_volume.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg * pyunits.L ** -1)):<20.4f}')
    print(f'{"Primary Pump Outlet":<34s}{m.fs.primary_pump.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.primary_pump.control_volume.properties_out[0].pressure, to_units=pyunits.bar)():<30.1f}{m.fs.primary_pump.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{value(pyunits.convert(m.fs.primary_pump.control_volume.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.mg * pyunits.L ** -1)):<20.4f}')
    
    print(f'{"Booster Pump Inlet":<34s}{m.fs.booster_pump.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.booster_pump.control_volume.properties_in[0].pressure, to_units=pyunits.bar)():<30.1f}{m.fs.booster_pump.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{m.fs.booster_pump.control_volume.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.4f}')
    print(f'{"Booster Pump Outlet":<34s}{m.fs.booster_pump.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(m.fs.booster_pump.control_volume.properties_out[0].pressure, to_units=pyunits.bar)():<30.1f}{m.fs.booster_pump.control_volume.properties_out[0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{m.fs.booster_pump.control_volume.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.4f}')
    

    for idx, unit in enumerate([m.fs.RO1, m.fs.RO2]):
        print(f'\n{"=======> RO Pass " + str(idx + 1) + " <=======":^60}\n')
        print(f'{str(unit.name).split(".")[1] + " Feed":<34s}{unit.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(unit.feed.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{unit.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{unit.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.4f}')
        print(f'{str(unit.name).split(".")[1] + " Product":<34s}{unit.product.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(unit.product.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{unit.product.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{unit.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.4f}')
        print(f'{str(unit.name).split(".")[1] + " Disposal":<34s}{unit.disposal.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(unit.disposal.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{unit.disposal.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{unit.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.4f}')

        for idx, stage in unit.stage.items():
            print(f'{str(unit.name).split(".")[1] + " Stage " + str(idx) + " Feed":<34s}{stage.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.module.feed_side.properties[0, 0].pressure, to_units=pyunits.bar)():<30.1f}{stage.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.feed_side.properties[0,0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.4f}')
        for idx, stage in unit.stage.items():
            print(f'{str(unit.name).split(".")[1] + " Stage " + str(idx) + " Permeate":<34s}{stage.permeate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.permeate.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{stage.permeate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.4f}')
        for idx, stage in unit.stage.items():
            print(f'{str(unit.name).split(".")[1] + " Stage " + str(idx) + " Retentate":<34s}{stage.retentate.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].value:<30.3f}{pyunits.convert(stage.retentate.properties[0.0].pressure, to_units=pyunits.bar)():<30.1f}{stage.retentate.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].value:<20.3e}{stage.module.feed_side.properties[0.0,1.0].conc_mass_phase_comp["Liq", "NaCl"].value:<20.4f}')
    
def display_system_metrics(m):
    print('\n')
    print(f'{"STAGE":<15s}{"RECOVERY %":<15s}{"REJECTION %":<15s}{"AREA (SQ M)":<15s}{"WIDTH (M)":<15s}{"LENGTH (M)":<15s}{"INLET VEL (M/S)":<20s}{"EXIT VEL (M/S)":<20s}{"EXIT DRIVING FORCE (BAR)":<15s}')
    for idx1, unit in enumerate([m.fs.RO1, m.fs.RO2]):
        for idx, stage in unit.stage.items():
            del_pi = value(pyunits.convert(stage.module.feed_side.properties_interface[0.0,1.0].pressure_osm_phase["Liq"] - stage.module.permeate_side[0.0,1.0].pressure_osm_phase["Liq"], to_units=pyunits.bar))
            del_P = value(pyunits.convert(stage.module.feed_side.properties_interface[0.0,1.0].pressure - stage.module.permeate_side[0.0,1.0].pressure, to_units=pyunits.bar))
            print(f'{str(unit.name).split(".")[1] + " Stage " + str(idx):<15s}{100*stage.module.recovery_vol_phase[0.0, "Liq"].value:<15.1f}{stage.module.rejection_phase_comp[0.0, "Liq", "NaCl"].value:<15.3f}{stage.module.area.value:<15.3f}{stage.module.width.value:<15.3f}{stage.module.length.value:<15.3f}{stage.module.feed_side.velocity[0, 0].value:<20.3f}{stage.module.feed_side.velocity[0, 1].value:<20.3f}{del_P-del_pi:<15.3f}')
    print('\n')
    print(f'{"Recirculation %":<20s}{100*(value(m.fs.splitter.split_fraction[0, "recirc"])):<10.1f}{"%":<20s}')
    print(f'{"1st Pass Recovery":<20s}{100*(value(m.fs.RO1.product.properties[0.0].flow_mass_phase_comp["Liq", "H2O"]/m.fs.RO1.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"])):<10.1f}{"%":<20s}')
    print(f'{"2nd Pass Recovery":<20s}{100*(value(m.fs.RO2.product.properties[0.0].flow_mass_phase_comp["Liq", "H2O"]/m.fs.RO2.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"])):<10.1f}{"%":<20s}')
    print(f'{"System Recovery":<20s}{100*m.fs.water_recovery.value:<10.1f}{"%":<20s}')
    print(f'{"System SEC":<20s}{value(m.fs.costing.specific_energy_consumption):<10.3f}{str(pyunits.get_units(m.fs.costing.specific_energy_consumption)):<20s}')
    print(f'{"System LCOW":<20s}{value(m.fs.costing.LCOW):<10.2f}{str(pyunits.get_units(m.fs.costing.LCOW)):<20s}')
    print('\n')

if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    main()
    