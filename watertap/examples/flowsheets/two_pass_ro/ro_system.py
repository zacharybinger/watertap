import os
import math
from pyomo.environ import (
    ConcreteModel,
    value,
    Param,
    Var,
    Constraint,
    NonNegativeReals,
    Block,
    RangeSet,
    units as pyunits,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state as _prop_state
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslogger
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)

from watertap.unit_models.pressure_changer import Pump
from watertap.core.util.initialization import *
from idaes.models.unit_models.mixer import (
    Mixer,
)
from watertap.costing import (
    WaterTAPCosting,
)
from idaes.models.unit_models import Product, Feed, StateJunction
from idaes.core.util.model_statistics import *
from watertap.costing import WaterTAPCosting
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock


def propagate_state(arc):
    _prop_state(arc)
    # print(f"Propogation of {arc.source.name} to {arc.destination.name} successful.")
    # arc.source.display()
    # print(arc.destination.name)
    # arc.destination.display()
    # print('\n')

def _initialize(m, blk, optarg):
    try:
        blk.initialize()
        blk.report()
        # print_close_to_bounds(m)
    except:
        print("----------------------------------\n")
        print(f"Initialization of {blk.name} failed.")
        print("\n----------------------------------\n")
        
        # blk.display()
        blk.report()
        print_infeasible_bounds(blk)
        print_close_to_bounds(blk)
        assert False

_log = idaeslog.getModelLogger("my_model", level=idaeslog.DEBUG, tag="model")

def build_ro(m, blk, number_of_stages=3, ultra_pute_water=False) -> None:
    print(f'\n{"=======> BUILDING RO SYSTEM <=======":^60}\n')

    blk.feed = StateJunction(property_package=m.fs.properties)
    blk.product = StateJunction(property_package=m.fs.properties)
    blk.disposal = StateJunction(property_package=m.fs.properties)
    blk.numberOfStages = Param(initialize=number_of_stages)
    blk.Stages = RangeSet(blk.numberOfStages)
    blk.booster_pumps = False

    blk.FirstStage = blk.Stages.first()
    blk.LastStage = blk.Stages.last()
    blk.NonFinalStages = RangeSet(number_of_stages - 1)

    blk.primary_mixer = Mixer(
        property_package=m.fs.properties,
        has_holdup=False,
        num_inlets = number_of_stages,
    )

    blk.stage = FlowsheetBlock(
        RangeSet(number_of_stages),
        dynamic=False)
    
    for idx, stage in blk.stage.items():
        if stage.index() > 1:
            build_ro_stage(m, stage, booster_pump=blk.booster_pumps)
        else:
            build_ro_stage(m, stage)

    blk.ro_feed_to_first_stage = Arc(
        source=blk.feed.outlet,
        destination=blk.stage[1].feed.inlet,
    )

    blk.stage_retentate_to_next_stage = Arc(
        blk.NonFinalStages,
        rule=lambda blk, n: {
            "source": blk.stage[n].retentate.outlet,
            "destination": blk.stage[n + 1].feed.inlet,
        },
    )

    blk.stage_permeate_to_mixer = Arc(
        blk.Stages,
        rule=lambda blk, n: {
            "source": blk.stage[n].permeate.outlet,
            "destination": getattr(blk.primary_mixer, "inlet_" + str(n)),
        },
    )

    blk.primary_mixer_to_product = Arc(
        source=blk.primary_mixer.outlet,
        destination=blk.product.inlet,
    )

    blk.last_stage_retentate_to_ro_retentate = Arc(
        source=blk.stage[number_of_stages].retentate.outlet, destination=blk.disposal.inlet
    )

    blk.feed.properties[0].conc_mass_phase_comp
    blk.product.properties[0].conc_mass_phase_comp
    blk.disposal.properties[0].conc_mass_phase_comp

    if ultra_pute_water:
        release_constraints_for_ultrapure_water(m, blk)

def build_ro_stage(m, blk, booster_pump=False):
    # Define IO
    blk.feed = StateJunction(property_package=m.fs.properties)
    blk.permeate = StateJunction(property_package=m.fs.properties)
    blk.retentate = StateJunction(property_package=m.fs.properties)
    blk.has_booster_pump = booster_pump

    if booster_pump:
        blk.booster_pump = Pump(property_package=m.fs.properties)

    blk.module = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=10,
        has_full_reporting = True
    )

    blk.module.eq_area.deactivate()
    @blk.module.Constraint(doc="Total Membrane area")
    def eq_area_new(b):
        return b.area == b.length * 2 * b.width
    
    blk.module.eq_mass_flux_equal_mass_transfer.deactivate()
    @blk.module.Constraint(
        blk.flowsheet().config.time,
        blk.module.difference_elements,
        blk.module.config.property_package.phase_list,
        blk.module.config.property_package.component_list,
        doc="Mass transfer term",
        )
    def eq_mass_flux_equal_mass_transfer_new(b, t, x, p, j):
        return (
            b.flux_mass_phase_comp[t, x, p, j] * (b.area / b.length)
            == -b.feed_side.mass_transfer_term[t, x, p, j]
        )

    if booster_pump:
        blk.stage_feed_to_booster_pump = Arc(
            source=blk.feed.outlet,
            destination=blk.booster_pump.inlet,
        )
        blk.stage_booster_pump_to_module = Arc(
            source=blk.booster_pump.outlet,
            destination=blk.module.inlet,
        )
    else:
        blk.stage_feed_to_module = Arc(
            source=blk.feed.outlet,
            destination=blk.module.inlet,
        )
    
    
    blk.stage_module_to_permeate = Arc(
        source=blk.module.permeate,
        destination=blk.permeate.inlet,
    )
    
    blk.stage_module_to_retentate = Arc(
        source=blk.module.retentate,
        destination=blk.retentate.inlet,
    )

def init_ro_system(m, blk, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING RO SYSTEM --------------------\n\n")
    # print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    # print(f"RO Degrees of Freedom: {degrees_of_freedom(blk)}")
    # for stage in blk.stage.values():
    #     print(f"RO Stage {stage} Degrees of Freedom: {degrees_of_freedom(stage)}")
    # print('\n\n')
    # assert_no_degrees_of_freedom(m)

    blk.feed.initialize(optarg=optarg)
    propagate_state(blk.ro_feed_to_first_stage)

    for stage in blk.stage.values():
        init_ro_stage(m, stage, solver=solver)
        if stage.index() < blk.numberOfStages:
            propagate_state(blk.stage_retentate_to_next_stage[stage.index()])
            propagate_state(blk.stage_permeate_to_mixer[stage.index()])
        else:
            propagate_state(blk.last_stage_retentate_to_ro_retentate)
            propagate_state(blk.stage_permeate_to_mixer[stage.index()])

    blk.disposal.initialize(optarg=optarg)
    blk.primary_mixer.initialize(optarg=optarg)
    propagate_state(blk.primary_mixer_to_product)
    blk.product.initialize(optarg=optarg)
    # _initialize(m, blk.product, optarg)

    print("\n\n-------------------- RO INITIALIZATION COMPLETE --------------------\n\n")
    # print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")
    # print(f"RO Degrees of Freedom: {degrees_of_freedom(blk)}")
    # for stage in blk.stage.values():
    #     print(f"RO Stage {stage} Degrees of Freedom: {degrees_of_freedom(stage)}")
    # print('\n\n')

def init_ro_stage(m, stage, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    if stage.has_booster_pump:
        stage.feed.initialize(optarg=optarg)
        propagate_state(stage.stage_feed_to_booster_pump)
        stage.booster_pump.initialize(optarg=optarg)
        propagate_state(stage.stage_booster_pump_to_module)
    else:
        stage.feed.initialize(optarg=optarg)
        propagate_state(stage.stage_feed_to_module)

    _initialize(m, stage.module, optarg)
    propagate_state(stage.stage_module_to_retentate)
    propagate_state(stage.stage_module_to_permeate)

    stage.permeate.initialize(optarg=optarg)
    stage.retentate.initialize(optarg=optarg)

def release_constraints_for_ultrapure_water(m, blk):
    for idx, stage in blk.stage.items():
            # Release constraints related to low velocity and low concentration
            stage.module.feed_side.velocity[0, 1].setlb(0.0)
            stage.module.recovery_mass_phase_comp.setlb(1e-7)
            stage.module.mixed_permeate[0.0].conc_mass_phase_comp["Liq", "NaCl"].setlb(0)

            # Release constraints related to low concentration
            for item in [stage.module.permeate_side, stage.module.feed_side.properties_interface]:
                for idx, param in item.items():
                    if idx[1] > 0:
                        param.molality_phase_comp["Liq", "NaCl"].setlb(0)
                        param.pressure_osm_phase["Liq"].setlb(0)
                        param.conc_mass_phase_comp["Liq", "NaCl"].setlb(0)

            # Release constraints related to low velocity and low flux
            for item in [stage.module.flux_mass_phase_comp, stage.module.feed_side.K, stage.module.feed_side.N_Re]:
                for idx, param in item.items():
                    if idx[1] > 0:
                        param.setlb(0)

            # Release constraints related to low velocity
            for idx, param in stage.module.feed_side.friction_factor_darcy.items():
                if idx[1] > 0:
                    param.setub(100)

def calc_scale(value):
    return math.floor(math.log(value, 10))

def set_ro_system_operating_conditions(m, blk, mem_area=100, booster_pump_pressure=80e5):
    # parameters
    mem_A = 4.0 / 3.6e11  # membrane water permeability coefficient [m/s-Pa]
    mem_B = 0.1 / 1000.0 / 3600.0  # membrane salt permeability coefficient [m/s]
    height = 1e-3  # channel height in membrane stage [m]
    spacer_porosity = 0.8  # spacer porosity in membrane stage [-]
    width = 500  # effective membrane width [m]
    area = mem_area # membrane area [m^2]
    length = 1.016*6  # membrane length [m]
    pressure_atm = 101325  # atmospheric pressure [Pa]
    pump_efi = 0.8  # pump efficiency [-]

    # blk.stage[1].module.feed_side.velocity[0, 0].fix(0.35)

    for idx, stage in blk.stage.items():
        stage.module.A_comp.fix(mem_A)
        stage.module.B_comp.fix(mem_B)
        stage.module.area.fix(area[idx-1])
        stage.module.length.fix(length)
        stage.module.width.unfix()
        stage.module.mixed_permeate[0].pressure.fix(pressure_atm)

        if stage.has_booster_pump:
            stage.booster_pump.control_volume.properties_out[0].pressure.fix(booster_pump_pressure)
            iscale.set_scaling_factor(stage.module.area, 1e-3)

        # stage.module.feed_side.velocity[0, 0].fix(0.25)
        if (
            stage.module.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
        ) or stage.module.config.pressure_change_type == PressureChangeType.calculated:
            stage.module.feed_side.channel_height.fix(height)
            stage.module.feed_side.spacer_porosity.fix(spacer_porosity)

        iscale.set_scaling_factor(stage.module.area, 1e-3)
        iscale.set_scaling_factor(stage.module.feed_side.area, 1e-3)
        iscale.set_scaling_factor(stage.module.width, 1e-3)

    # ---checking model---
    assert_units_consistent(m)

def get_sub_blocks(block, decend = False, report=False):
    blocks = []
    for v in block.component_data_objects(ctype=Block, active=True, descend_into=decend):
        print(v)
        if report:
            try:
                table = v._get_stream_table_contents()
                for item in table:
                    print(table[item])
                # print(v._get_stream_table_contents())
            except:
                pass
        
def display_ro_system_build(m):
    get_sub_blocks(m.fs)
    # get_sub_blocks(m.fs.ro)
    # for stage in m.fs.ro.stage.values():
    #     get_sub_blocks(stage)
    print('\n')

if __name__ == "__main__":
    file_dir = os.path.dirname(os.path.abspath(__file__))
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.disposal = Product(property_package=m.fs.properties)
    
    m.fs.ro = FlowsheetBlock(dynamic=False)
    build_ro(m,m.fs.ro, number_of_stages=2)
    display_ro_system_build(m)