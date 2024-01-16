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
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state as _prop_state
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    calculate_scaling_factors,
    set_scaling_factor,
)
import idaes.logger as idaeslogger
from idaes.core.util.exceptions import InitializationError
from watertap.core.util.model_diagnostics.infeasible import *

from watertap.unit_models.pressure_changer import Pump
from watertap.core.util.initialization import *


# from analysisWaterTAP.utils import flowsheet_utils as fsTool
# from analysisWaterTAP.flowsheets.lssro_oaro.costing.LSRRO_ORARO_costing import *
from idaes.models.unit_models import CSTR, Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *
from watertap.costing import WaterTAPCosting
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock

# Import Anaerobic Digertor Model properties 
from watertap.property_models.anaerobic_digestion.adm1_properties import (
    ADM1ParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_reactions import (
    ADM1ReactionParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_properties_vapor import (
    ADM1_vaporParameterBlock,
)

# Import Activated Sludge Model properties 
from watertap.property_models.activated_sludge.asm1_properties import (
    ASM1ParameterBlock,
)
from watertap.property_models.activated_sludge.asm1_reactions import (
    ASM1ReactionParameterBlock,
)
from watertap.unit_models.uv_aop import Ultraviolet0D, UVDoseType
from watertap.unit_models.gac import (
    GAC,
    FilmTransferCoefficientType,
    SurfaceDiffusionCoefficientType,
)
# Import anaerobic digestor model
from watertap.unit_models.anaerobic_digestor import AD
# Import CSTR with oxygen injection model
from watertap.unit_models.cstr_injection import CSTR_Injection
from watertap.unit_models.zero_order import UltraFiltrationZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting

# import translator models from WaterTAP
from watertap.unit_models.translators.translator_asm1_adm1 import Translator_ASM1_ADM1
from watertap.unit_models.translators.translator_adm1_asm1 import Translator_ADM1_ASM1

def propagate_state(arc):
    _prop_state(arc)

def _initialize(m, blk, optarg):
    try:
        blk.initialize()
    except:
        print("----------------------------------\n")
        print(f"Initialization of {blk.name} failed.")
        print("\n----------------------------------\n")
        
        # blk.display()
        blk.report()
        print_infeasible_bounds(m)
        print_close_to_bounds(m)
        print_infeasible_constraints(m)
        assert False

_log = idaeslog.getModelLogger("my_model", level=idaeslog.DEBUG, tag="model")

def build_MBR(m, blk) -> None:
    print(f'\n{"=======> BUILDING MBR SYSTEM <=======":^60}\n')

    blk.feed = StateJunction(property_package=m.fs.properties)
    blk.product = StateJunction(property_package=m.fs.properties)
    blk.disposal = StateJunction(property_package=m.fs.properties)
    

    blk.aerobic_stage =  CSTR_Injection(
    property_package=m.fs.props_ASM1, reaction_package=m.fs.ASM1_rxn_props)

    blk.UF = UltraFiltrationZO(property_package=m.fs.params, database=m.db)

 
def add_connections(m):
    print(f'\n{"=======> ADDING CONNECTIONS <=======":^60}\n')
    m.fs.feed_to_translator = Arc(source=m.fs.feed, destination=m.fs.asm_adm.inlet)
    m.fs.translator_to_aerobic = Arc(source=m.fs.asm_adm.outlet, destination=m.fs.aerobic_stage.inlet)
    m.fs.aerobic_to_uf = Arc(source=m.fs.aerobic_stage.outlet, destination=m.fs.UF.inlet)
    m.fs.uf_to_product = Arc(source=m.fs.UF.retentate, destination=m.fs.product)

def init_MBR(m, blk, optarg):
    blk.feed.initialize(optarg=optarg)
    propagate_state(blk.feed_to_translator)
    propagate_state(blk.translator_to_aerobic)
    blk.aerobic_stage.initialize(optarg=optarg)
    propagate_state(blk.aerobic_to_uf)
    blk.UF.initialize(optarg=optarg)
    propagate_state(blk.uf_to_product)
    blk.product.initialize(optarg=optarg)


def main():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props_ASM1 = ASM1ParameterBlock()
    m.fs.props_ADM1 = ADM1ParameterBlock()
    m.fs.props_vap = ADM1_vaporParameterBlock()
    m.fs.ADM1_rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props_ADM1)
    m.fs.ASM1_rxn_props = ASM1ReactionParameterBlock(property_package=m.fs.props_ASM1)
    m.fs.params = WaterParameterBlock(
            solute_list=["eeq", "toc", "tss", "cryptosporidium"]
        )
    
    m.fs.feed = Feed(property_package=m.fs.props_ASM1)
    m.fs.product = Product(property_package=m.fs.props_ASM1)
    # m.fs.disposal = Product(property_package=m.fs.properties)

    m.fs.MBR = FlowsheetBlock(dynamic=False)
    build_MBR(m, m.fs.MBR)
