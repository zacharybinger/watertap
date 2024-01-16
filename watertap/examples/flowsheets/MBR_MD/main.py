import os
import math
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
    TransformationFactory,
    Block,
    RangeSet,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util.initialization import *
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from watertap.core.wt_database import Database

from watertap.property_models.NaCl_prop_pack import NaClParameterBlock

from idaes.models.unit_models import Product, Feed, StateJunction, Separator
from idaes.core.util.model_statistics import *
from watertap.costing import WaterTAPCosting
from watertap.unit_models.pressure_changer import Pump
from idaes.core import FlowsheetBlock, UnitModelCostingBlock

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

# import translator models from WaterTAP
from watertap.unit_models.translators.translator_asm1_adm1 import Translator_ASM1_ADM1
from watertap.unit_models.translators.translator_adm1_asm1 import Translator_ADM1_ASM1

from components.anaerobic_anoxic import (
    build_anaerobic_anoxic,
    init_anaerobic_anoxic,
)
from components.MBR import (
    build_MBR,
    init_MBR)
from components.uv_aop_gac import (
    build_UV_GAC,
    init_UV_GAC,
)

def main():
    m = build_system()
    add_connections(m)
    add_translator(m)
    add_constraints(m)
    set_inlet_conditions(m)
    # init_system(m)
    return m

def build_system():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props_ASM1 = ASM1ParameterBlock()
    m.fs.props_ADM1 = ADM1ParameterBlock()
    m.fs.props_vap = ADM1_vaporParameterBlock()
    m.fs.ADM1_rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props_ADM1)
    m.fs.ASM1_rxn_props = ASM1ReactionParameterBlock(property_package=m.fs.props_ASM1)
    # m.fs.params = WaterParameterBlock(
    #         solute_list=["eeq", "toc", "tss", "cryptosporidium"]
    #     )
    
    m.fs.feed = Feed(property_package=m.fs.props_ASM1)
    m.fs.product = Product(property_package=m.fs.props_ASM1)
    m.fs.disposal = Product(property_package=m.fs.properties)

    m.fs.anaerobic_anoxic = FlowsheetBlock(dynamic=False)
    build_anaerobic_anoxic(m, m.fs.anaerobic_anoxic)

    m.fs.MBR = FlowsheetBlock(dynamic=False)
    build_MBR(m, m.fs.MBR)

    m.fs.UV_GAC = FlowsheetBlock(dynamic=False)
    build_UV_GAC(m, m.fs.UV_GAC)

    return m

def add_connections(m):

    m.fs.feed_to_anaerobic = Arc(
        source=m.fs.feed, 
        destination=m.fs.anaerobic_anoxic.feed
    )

    m.fs.anaerobic_anoxic_to_MBR = Arc(
        source=m.fs.anaerobic_anoxic.product,
        destination=m.fs.MBR.feed
    )

    m.fs.MBR_to_UV_GAC = Arc(
        source=m.fs.MBR.product,
        destination=m.fs.UV_GAC.feed
    )

    m.fs.UV_GAC_to_product = Arc(
        source=m.fs.UV_GAC.product,
        destination=m.fs.product
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

def add_constraints(m):
    pass

def add_translator(m) -> None:
    print(f'\n{"=======> BUILDING TRANSLATOR <=======":^60}\n')

    m.fs.asm_adm = Translator_ASM1_ADM1(
        inlet_property_package=m.fs.props_ASM1,
        outlet_property_package=m.fs.props_ADM1,
        reaction_package=m.fs.ADM1_rxn_props,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    m.fs.adm_asm = Translator_ADM1_ASM1(
        inlet_property_package=m.fs.props_ADM1,
        outlet_property_package=m.fs.props_ASM1,
        reaction_package=m.fs.ADM1_rxn_props,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )
   
def set_inlet_conditions(m):
    pass

def init_system(m, verbose=True, solver=None):
    if solver is None:
        solver = get_solver()

    optarg = solver.options

    print("\n\n-------------------- INITIALIZING SYSTEM --------------------\n\n")
    print(f"System Degrees of Freedom: {degrees_of_freedom(m)}")

    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.feed_to_anaerobic)
    
    init_anaerobic_anoxic(m, m.fs.anaerobic_anoxic, optarg)
    propagate_state(m.fs.anaerobic_anoxic_to_MBR)

    init_MBR(m, m.fs.MBR, optarg)
    propagate_state(m.fs.MBR_to_UV_GAC)

    init_UV_GAC(m, m.fs.UV_GAC, optarg)
    propagate_state(m.fs.UV_GAC_to_product)

    m.fs.product.initialize(optarg=optarg)
    m.fs.disposal.initialize(optarg=optarg)


