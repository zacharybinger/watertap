import os
import math
import numpy as np
from pyomo.environ import (
    ConcreteModel,
    value,
    Param,
    Var,
    Constraint,
    Set,
    Expression,
    Objective,
    NonNegativeReals,
    Block,
    RangeSet,
    check_optimal_termination,
    units as pyunits,
)

from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent

from idaes.models.unit_models import CSTR, Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *

from watertap.unit_models.uv_aop import Ultraviolet0D, UVDoseType
from watertap.unit_models.gac import (
    GAC,
    FilmTransferCoefficientType,
    SurfaceDiffusionCoefficientType,
)
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state as _prop_state

def propagate_state(arc):
    _prop_state(arc)

def build_UV_GAC(m, blk) -> None:
    print(f'\n{"=======> BUILDING POLISHING SYSTEM <=======":^60}\n')

    blk.feed = StateJunction(property_package=m.fs.properties)
    blk.product = StateJunction(property_package=m.fs.properties)
    
    blk.uv_aop = Ultraviolet0D(
        property_package=m.fs.properties,
        database=m.db,
        process_subtype="hydrogen_peroxide",
    )

    blk.gac = GAC(
            property_package=m.fs.properties,
            film_transfer_coefficient_type="fixed",
            surface_diffusion_coefficient_type="fixed",
        )
    
    blk.feed_to_UV = Arc(
        source=blk.feed.outlet,
        destination=blk.uv_aop.inlet,
    )

    blk.UV_to_GAC = Arc(
        source=blk.uv_aop.outlet,
        destination=blk.gac.inlet,
    )

    blk.GAC_to_product = Arc(
        source=blk.gac.outlet,
        destination=blk.product,
    )
    
def init_UV_GAC(m, blk, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING UV GAC Polishing --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(blk)}")
    print('\n\n')

    blk.feed.initialize(optarg=optarg)
    propagate_state(blk.feed_to_UV)
    blk.uv_aop.initialize(optarg=optarg)
    propagate_state(blk.UV_to_GAC)
    blk.gac.initialize(optarg=optarg)
    propagate_state(blk.GAC_to_product)
    blk.product.initialize(optarg=optarg)

def set_UV_GAC_operating_conditions(m,blk):
    blk.feed.pressure.fix(101325)
    blk.feed.temperature.fix(298.15)

    uv_intensity = 1 * pyunits.mW / pyunits.cm**2
    exporure_time = 500 * pyunits.s
    inactivation_rate = 2.3 * pyunits.cm**2 / pyunits.J
    EEO = 0.25 * pyunits.kWh / pyunits.m**3
    lamp_efficiency = 0.8

    m.fs.unit.uv_intensity.fix(uv_intensity)
    m.fs.unit.exposure_time.fix(exporure_time)
    # m.fs.unit.inactivation_rate["Liq", "NDMA"].fix(inactivation_rate)
    m.fs.unit.electrical_efficiency_phase_comp[0, "Liq", "NDMA"].fix(EEO)
    m.fs.unit.lamp_efficiency.fix(lamp_efficiency) 

    blk.gac.freund_k.fix(37.9e-6 * (1e6**0.8316))
    blk.gac.freund_ninv.fix(0.8316)
    blk.gac.particle_dens_app.fix(722)
    blk.gac.particle_dia.fix(0.00106)
    blk.gac.ebct.fix(300)  # seconds
    blk.gac.bed_voidage.fix(0.449)
    blk.gac.bed_length.fix(6)  # assumed
    blk.gac.conc_ratio_replace.fix(0.50)
    blk.gac.kf.fix(3.29e-5)
    blk.gac.ds.fix(1.77e-13)
    blk.gac.a0.fix(3.68421)
    blk.gac.a1.fix(13.1579)
    blk.gac.b0.fix(0.784576)
    blk.gac.b1.fix(0.239663)
    blk.gac.b2.fix(0.484422)
    blk.gac.b3.fix(0.003206)
    blk.gac.b4.fix(0.134987)