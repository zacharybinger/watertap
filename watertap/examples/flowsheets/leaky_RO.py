from pyomo.environ import Param,TransformationFactory
from pyomo.network import Arc
from pyomo.environ import units as pyunits
from idaes.core.util.scaling import constraint_scaling_transform
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import Product, Feed
from idaes.core import UnitModelCostingBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.costing import WaterTAPCosting
from watertap.unit_models.reverse_osmosis_0D import (ReverseOsmosis0D, ConcentrationPolarizationType, MassTransferCoefficient, PressureChangeType)

# Import concrete model from Pyomo
from pyomo.environ import ConcreteModel, Var, Reals, Objective, Constraint, value, units, NonNegativeReals
# Import flowsheet block from IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor
# Import NaCl property model
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock


# flowsheet set up
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.properties = SeawaterParameterBlock() # or NaClParameterBlock

# Control volume flow blocks
m.fs.feed = m.fs.properties.build_state_block([0])
m.fs.perm = m.fs.properties.build_state_block([0])

# fix state variables
m.fs.feed[0].temperature.fix(273 + 25)                      # temperature (K)
m.fs.feed[0].pressure.fix(75e5)                           # pressure (Pa)
m.fs.feed[0].flow_mass_phase_comp['Liq', 'H2O'].fix(0.965)  # mass flowrate of H2O (kg/s)
m.fs.feed[0].flow_mass_phase_comp['Liq', 'TDS'].fix(0.035)  # mass flowrate of TDS (kg/s)
m.fs.feed[0].pressure_osm_phase

# fix state variables
m.fs.perm[0].temperature.fix(273 + 25)                      # temperature (K)
m.fs.perm[0].pressure.fix(101325)                           # pressure (Pa)
m.fs.perm[0].flow_mass_phase_comp['Liq', 'H2O'].fix(1)  # mass flowrate of H2O (kg/s)
m.fs.perm[0].flow_mass_phase_comp['Liq', 'TDS'].fix(0.00035)  # mass flowrate of TDS (kg/s)
m.fs.perm[0].pressure_osm_phase

# set default property values
m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
m.fs.properties.set_default_scaling(
    "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
)

m.fs.flux_mass_phase_comp = Var(
            initialize= 1e-6,
            units=pyunits.kg
            * pyunits.m ** -2
            * pyunits.second ** -1,
            doc="Mass flux across membrane at inlet and outlet",
        )

m.fs.A_prime = Var(
        initialize=1e-8,
        domain=NonNegativeReals,
        units=pyunits.m
                * pyunits.second ** -1
                * pyunits.kPa ** -1,
        doc="Water permeability coefficient of the membrane",
    )

m.fs.A_prime.fix(1e-9)

m.fs.flux_constraint = Constraint(expr=m.fs.flux_mass_phase_comp == m.fs.A_prime * (
        (m.fs.feed[0].pressure - m.fs.perm[0].pressure) - (
        m.fs.feed[0].pressure_osm_phase['Liq'] - m.fs.perm[0].pressure_osm_phase['Liq']))
        )

solver = get_solver()
results = solver.solve(m)

liq = "Liq"
h2o = "H2O"

print('\n')
print('=============== SIMULATION RESULTS ===============\n')
print(
        f'{"Water Permeability Coeff":<30s}{f"{value(m.fs.A_prime):<10,.2e}"}{"m/s*kPa":<10s}'
    )
print('\n')
print(
        f'{"Feed Pressure":<30s}{f"{value(units.convert(m.fs.feed[0].pressure() * pyunits.Pa, to_units=pyunits.bar)):<10,.2f}"}{"bar":<10s}'
    )
print(
         f'{"Feed Flow Rate":<30s}{f"{value(units.convert(m.fs.feed[0].flow_mass_phase_comp[liq, h2o], to_units=pyunits.kg / pyunits.second)):<10,.2f}"}{"kg/s":<10s}'
    )
print(
         f'{"Feed Osm Pressure":<30s}{f"{value(units.convert(m.fs.feed[0].pressure_osm_phase[liq], to_units=pyunits.bar)):<10,.2f}"}{"bar":<10s}'
    )
print('\n')
print(
        f'{"Perm Pressure":<30s}{f"{value(units.convert(m.fs.perm[0].pressure() * pyunits.Pa, to_units=pyunits.bar)):<10,.2f}"}{"bar":<10s}'
    )
print(
         f'{"Perm Flow Rate":<30s}{f"{value(units.convert(m.fs.perm[0].flow_mass_phase_comp[liq, h2o], to_units=pyunits.kg / pyunits.second)):<10,.2f}"}{"kg/s":<10s}'
    )
print(
         f'{"Perm Osm Pressure":<30s}{f"{value(units.convert(m.fs.perm[0].pressure_osm_phase[liq], to_units=pyunits.bar)):<10,.2f}"}{"bar":<10s}'
    )
print('\n')
print(
         f'{"Perm Density":<30s}{f"{value(m.fs.perm[0].dens_mass_solvent()):<10,.2f}"}{str(pyunits.get_units(m.fs.perm[0].dens_mass_solvent)):<10s}'
    )
print(
         f'{"Water Flux":<30s}{f"{value(m.fs.flux_mass_phase_comp):<10,.2f}"}{str(pyunits.get_units(m.fs.flux_mass_phase_comp)):<10s}'
    )
print(
         f'{"Water Flux":<30s}{f"{value(units.convert(m.fs.flux_mass_phase_comp / m.fs.perm[0].dens_mass_solvent, to_units=pyunits.liter * pyunits.m ** -2 * pyunits.hour ** -1)):<10,.2f}"}{"LMH":<10s}'
    )
print('\n')
# print(
#          f'{"Water Flux":<30s}{f"{value(units.convert((5.262E-7 * pyunits.m * pyunits.second ** -1), to_units=pyunits.liter * pyunits.m ** -2 * pyunits.hour ** -1)):<10,.2f}"}{"LMH":<10s}'
#     )
