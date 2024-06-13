import os
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    value,
    Var,
    NonNegativeReals,
    assert_optimal_termination,
    units as pyunits,
)
from idaes.core import FlowsheetBlock
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

from idaes.core.util.initialization import propagate_state
from pyomo.environ import TransformationFactory
from pyomo.network import Arc
from watertap.unit_models.pressure_changer import Pump

def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    # m.fs.disposal = Product(property_package=m.fs.properties)

    m.fs.P1 = Pump(property_package=m.fs.properties)
    m.fs.P2 = Pump(property_package=m.fs.properties)

    m.fs.M1 = Mixer(
        property_package=m.fs.properties,
        has_holdup=False,
        num_inlets = 2,
    )

    m.fs.RO = ReverseOsmosis0D(
    property_package=m.fs.properties,
    has_pressure_change=True,
    pressure_change_type=PressureChangeType.calculated,
    mass_transfer_coefficient=MassTransferCoefficient.calculated,
    concentration_polarization_type=ConcentrationPolarizationType.calculated,
    has_full_reporting = True,
    )

    # connect unit models
    m.fs.feed_to_M1 = Arc(source=m.fs.feed.outlet, destination=m.fs.M1.inlet_1)
    m.fs.M1_to_P1 = Arc(source=m.fs.M1.outlet, destination=m.fs.P1.inlet)
    m.fs.P1_to_RO = Arc(source=m.fs.P1.outlet, destination=m.fs.RO.inlet)
    m.fs.RO_to_product = Arc(source=m.fs.RO.permeate, destination=m.fs.product.inlet)
    m.fs.RO_to_P2 = Arc(source=m.fs.RO.retentate, destination=m.fs.P2.inlet)
    m.fs.P2_to_M1 = Arc(source=m.fs.P2.outlet, destination=m.fs.M1.inlet_2)
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

    return m


def set_operating_conditions(m):
    # feed, 4 degrees of freedom
    m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'H2O'].fix(0.965)                # volumetric flow rate (m3/s)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(0.035)  # TDS mass fraction (-)
    m.fs.feed.properties[0].pressure.fix(101325)                           # pressure (Pa)
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)                   # temperature (K)

    # high pressure pump, 2 degrees of freedom
    m.fs.P1.efficiency_pump.fix(0.80)                                    # pump efficiency (-)
    m.fs.P1.control_volume.properties_out[0].pressure.fix(75e5)          # pump outlet pressure (Pa)
    m.fs.P2.efficiency_pump.fix(0.80)                                    # pump efficiency (-)
    m.fs.P2.control_volume.properties_out[0].pressure.fix(75e5)          # pump outlet pressure (Pa)

    # RO unit, 7 degrees of freedom
    m.fs.RO.A_comp.fix(4.2e-12)                                            # membrane water permeability coeff (m/Pa/s)
    m.fs.RO.B_comp.fix(3.5e-8)                                             # membrane salt permeability coeff (m/s)

    # fix 4 module specficiations
    m.fs.RO.area.fix(50)                                                # membrane stage area (m^2)
    m.fs.RO.width.fix(5)                                                # membrane stage width (m)
    m.fs.RO.feed_side.channel_height.fix(1E-3)                          # channel height in membrane stage (m)
    m.fs.RO.feed_side.spacer_porosity.fix(0.97)                         # spacer porosity in membrane stage (-)
    m.fs.RO.permeate.pressure[0].fix(101325)                               # permeate pressure (Pa)

    print("DOF = ", degrees_of_freedom(m))
    assert_no_degrees_of_freedom(m)

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
    propagate_state(m.fs.feed_to_M1)

    if pass_num > 1:
        m.fs.M1.initialize()
    else:
        initialize_mixer(m)
    propagate_state(m.fs.M1_to_P1)

    m.fs.P1.initialize()
    propagate_state(m.fs.P1_to_RO)

    m.fs.RO.initialize()
    propagate_state(m.fs.RO_to_product)
    propagate_state(m.fs.RO_to_P2)

    m.fs.P2.initialize()
    propagate_state(m.fs.P2_to_M1)


def initialize(m):
    # initialize unit by unit
    for init_pass in range(5):
        do_forward_initialization_pass(m)


def solve(m):
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)

def print_results(m):
    print(m.fs.RO.report())
    print(m.fs.M1.report())

def main():
    m = build()
    set_operating_conditions(m)
    initialize(m)
    print_results(m)
    solve(m)


if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    main()