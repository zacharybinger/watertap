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

# Import anaerobic digestor model
from watertap.unit_models.anaerobic_digestor import AD

# import translator models from WaterTAP
from watertap.unit_models.translators.translator_asm1_adm1 import Translator_ASM1_ADM1
from watertap.unit_models.translators.translator_adm1_asm1 import Translator_ADM1_ASM1

def propagate_state(arc):
    _prop_state(arc)

def build_anaerobic_anoxic(m, blk) -> None:
    print(f'\n{"=======> BUILDING AEROBIC ANOXIC SYSTEM <=======":^60}\n')

    blk.feed = StateJunction(property_package=m.fs.properties)
    blk.product = StateJunction(property_package=m.fs.properties)
    blk.disposal = StateJunction(property_package=m.fs.properties)
    

    blk.anaerobic_stage =  AD(
    liquid_property_package=m.fs.props_ADM1,
    vapor_property_package=m.fs.props_vap,
    reaction_package=m.fs.ADM1_rxn_props,
    has_heat_transfer=True,
    has_pressure_change=False,
    )

    # First reactor (anoxic) - standard CSTR
    blk.anoxic_stage = CSTR(
    property_package=m.fs.props_ASM1, reaction_package=m.fs.ASM1_rxn_props
    )

    blk.feed_to_anaerobic = Arc(
        source=blk.feed.outlet,
        destination=blk.anaerobic_stage.inlet,
    )
    
    blk.anaerobic_to_anoxic = Arc(
        source=blk.anaerobic_stage.outlet,
        destination=blk.anoxic_stage.inlet,
    )

    blk.anoxic_to_product = Arc(
        source=blk.anoxic_stage.outlet,
        destination=blk.product,
    )

def init_anaerobic_anoxic(m, blk, optarg):
    
    blk.feed.initialize(optarg=optarg)
    propagate_state(blk.feed_to_anaerobic)
    blk.anaerobic_stage.initialize(optarg=optarg)
    propagate_state(blk.anaerobic_to_anoxic)
    blk.anoxic_stage.initialize(optarg=optarg)
    propagate_state(blk.anoxic_to_product)
    blk.product.initialize(optarg=optarg)
    propagate_state(blk.product)