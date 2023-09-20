import os
# Import concrete model from Pyomo
from pyomo.environ import ConcreteModel, Var, Reals, Objective, Constraint, Expression, value, exp, NonNegativeReals, units as pyunits
# Import flowsheet block from IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from watertap.core.util.initialization import assert_no_degrees_of_freedom
# Import NaCl property model
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from skk_util import *

import numpy as np
import matplotlib.pyplot as plt

par_dir = os.path.dirname(os.path.abspath(__file__))

def add_vars(m):

    m.fs.flux_mass_phase_comp = Var(
                initialize= 0,
                units=pyunits.kg
                * pyunits.m ** -2
                * pyunits.second ** -1,
                doc="Mass flux across membrane at inlet and outlet",
            )

    m.fs.salt_flux_mass_phase_comp = Var(
                initialize= 0,
                units=pyunits.kg
                * pyunits.m ** -2
                * pyunits.second ** -1,
                doc="Mass flux across membrane at inlet and outlet",
            )

    m.fs.A = Var(
            initialize=1e-8,
            domain=NonNegativeReals,
            units=pyunits.m
                    * pyunits.second ** -1
                    * pyunits.kPa ** -1,
            doc="Water permeability coefficient of the membrane",
        )

    m.fs.B = Var(
            initialize=1e-6,
            domain=NonNegativeReals,
            units=pyunits.m
                    * pyunits.second ** -1,
            doc="Salt permeability coefficient of the membrane",
        )

    m.fs.dens_solvent = Var(
            initialize=1000,
            units=pyunits.kg * pyunits.m ** -3,
            doc="Pure water density",
        )

    m.fs.reflect_coeff = Var(
            initialize=0.9,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Reflection coefficient of the membrane",
        )

    m.fs.alpha = Var(
                initialize=1,
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Alpha coefficient of the membrane",
        )
        
    m.fs.rejection_skk = Var(
                initialize=0.9,
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Theoretical solute rejection",
        )

def add_constraints(m):

    m.fs.alpha_constraint = Constraint(expr=m.fs.alpha == (1 - m.fs.reflect_coeff) / m.fs.B)
    
    m.fs.flux_constraint = Constraint(expr=m.fs.flux_mass_phase_comp == m.fs.A *
            m.fs.dens_solvent * (
                (
                m.fs.feed[0].pressure - 
                m.fs.perm[0].pressure
                ) - m.fs.reflect_coeff * (m.fs.feed[0].pressure_osm_phase['Liq'])
                ))

    m.fs.salt_flux_constraint = Constraint(expr=m.fs.salt_flux_mass_phase_comp == m.fs.B * (
            (m.fs.feed[0].conc_mass_phase_comp['Liq', 'TDS'])) +
            (1 - m.fs.reflect_coeff) * (
                (m.fs.flux_mass_phase_comp / m.fs.dens_solvent) * 
                m.fs.feed[0].conc_mass_phase_comp['Liq', 'TDS']))

    m.fs.rejection = Expression(expr= 1 - (
        m.fs.salt_flux_mass_phase_comp / (m.fs.flux_mass_phase_comp/1000)) / (
            m.fs.feed[0].conc_mass_phase_comp['Liq', 'TDS']
        )
        )

    m.fs.rejection_skk_constraint = Constraint(expr=m.fs.rejection_skk == m.fs.reflect_coeff*(
        1-exp(-1*m.fs.alpha*m.fs.salt_flux_mass_phase_comp)) / (
            1-exp(-1*m.fs.alpha*m.fs.salt_flux_mass_phase_comp)*m.fs.reflect_coeff
        ))
    
    m.fs.min_flux_constraint = Constraint(expr=m.fs.flux_mass_phase_comp >= 0)
    m.fs.min_salt_flux_constraint = Constraint(expr=m.fs.salt_flux_mass_phase_comp >= 0)

def build(transport_model='SD'):
    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock() # or NaClParameterBlock

    # Control volume flow blocks
    m.fs.feed = m.fs.properties.build_state_block([0])
    m.fs.perm = m.fs.properties.build_state_block([0])

    add_vars(m)
    add_constraints(m)

    if transport_model == 'SD':
        m.fs.reflect_coeff.fix(1)

    return m

def set_operating_conditions(m, 
                             tds=0.035,
                             pressure = 50e5,
                             A_LMH = 10,
                             B_LMH = 1,
                             reflect_coeff=None):
    # fix state variables
    m.fs.feed[0].temperature.fix(273 + 25)                      # temperature (K)
    m.fs.feed[0].pressure.fix(pressure)                           # pressure (Pa)
    m.fs.feed[0].flow_mass_phase_comp['Liq', 'H2O'].fix(1-tds)  # mass flowrate of H2O (kg/s)
    m.fs.feed[0].flow_mass_phase_comp['Liq', 'TDS'].fix(tds)  # mass flowrate of TDS (kg/s)
    # m.fs.feed[0].pressure_osm_phase
    # m.fs.feed[0].mass_frac_phase_comp

    m.fs.dens_solvent.fix(1000)
    # m.fs.perm[0].flow_mass_phase_comp['Liq', 'H2O'].fix(0)
    # m.fs.perm[0].flow_mass_phase_comp['Liq', 'TDS'].fix(0)
    m.fs.perm[0].temperature.fix(273 + 25)                      # temperature (K)
    m.fs.perm[0].pressure.fix(101325)                           # pressure (Pa)
    # m.fs.perm[0].pressure_osm_phase
    
    m.fs.A.fix(A_LMH/(1000*3600*1e5))                                      # membrane water permeability (m/Pa/s)
    m.fs.B.fix(B_LMH/(1000*3600))
    # m.fs.A.fix(4.2e-12)                                      # membrane water permeability (m/Pa/s)
    # m.fs.B.fix(3.5e-8)

    if reflect_coeff != None:
        m.fs.reflect_coeff.fix(reflect_coeff)

    return m

def initialize(m):
    m.fs.feed.initialize()
    m.fs.perm.initialize()

def set_scaling(m):
    # set default property values
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )

def param_sweep_2D():
    solver = get_solver()
    m = build(transport_model='SKK')
    set_operating_conditions(m, pressure=30e5, tds = 0.035, A_LMH = 5, B_LMH = 13.7)
    fig, ax = plt.subplots(figsize=(8, 6))
    tds_range = np.linspace(0.03, 0.07, 3)
    pressure_range = np.linspace(30e5, 50e5, 3)
    flux_results = []
    salt_flux_results = []
    reports = []
    for idx, conc in enumerate(tds_range):
        for idx, pressure in enumerate(pressure_range):
            set_operating_conditions(m, pressure=pressure, tds = conc, A_LMH = 5, B_LMH = 13.7)
            m.fs.reflect_coeff.fix(0.9)
            results = solver.solve(m)
            flux_results.append(value(pyunits.convert(m.fs.flux_mass_phase_comp, to_units=pyunits.kg /pyunits.m **2 / pyunits.hr)))
            salt_flux_results.append(value(pyunits.convert(m.fs.salt_flux_mass_phase_comp, to_units=pyunits.kg /pyunits.m **2 / pyunits.hr)))
            reports.append(create_report(m))
    report = pd.concat(reports)
    report.to_csv(os.path.join(par_dir, 'reports/', 'SKK_stateblock_conc_vs_pressure.csv'))
    print(report)

def reflection_sensitivity():
    solver = get_solver()
    m = build(transport_model='SKK')
    # set_operating_conditions(m, pressure=30e5, tds = 0.035, A_LMH = 5, B_LMH = 13.7)
    set_operating_conditions(m, pressure=65e5, tds = 0.1, A_LMH = 10.1, B_LMH = 13.7, reflect_coeff=0.95)
    fig, ax = plt.subplots(figsize=(8, 6))
    ref_coeff_range = np.linspace(0.7, 0.95, 10)
    alpha_range = np.linspace(1E9, 1E10, 10)
    flux_results = []
    salt_flux_results = []
    rejection_obs_results = []
    rejection_skk_results = []
    reports = []
    for idx, ref_coeff in enumerate(ref_coeff_range):
        m.fs.reflect_coeff.fix(ref_coeff)
        results = solver.solve(m)
#         # print_results(m)
        flux_results.append(value(pyunits.convert(m.fs.flux_mass_phase_comp, to_units=pyunits.kg /pyunits.m **2 / pyunits.hr)))
        salt_flux_results.append(value(pyunits.convert(m.fs.salt_flux_mass_phase_comp, to_units=pyunits.kg /pyunits.m **2 / pyunits.hr)))
        rejection_obs_results.append(value(m.fs.rejection*100))
        rejection_skk_results.append(value(m.fs.rejection_skk*100))
        reports.append(create_report(m))

    report = pd.concat(reports)
    report.to_csv(os.path.join(par_dir, 'reports/', 'SKK_stateblock_sensitivity.csv'))
    print(report)
    # plot_data(ref_coeff_range, flux_results, salt_flux_results, ax=ax, xlabel=r'Reflect Coeff. $(\sigma)$', ylabel=r'Water Flux $(\frac{L}{m^{2}  hr})$', ylabel2=r'Salt Flux $(\frac{kg}{m^{2}  hr})$')
    plot_data(ref_coeff_range, flux_results, salt_flux_results, rej_skk=rejection_skk_results, ax=ax, xlabel=r'Reflect Coeff. $(\frac{1-\sigma}{\omega^\prime})$', ylabel=r'Water Flux $(\frac{L}{m^{2}  hr})$', ylabel2=r'Salt Flux $(\frac{kg}{m^{2}  hr})$')
    fig.tight_layout()
    plt.show()

def simple_reflection_sensitivity():
    solver = get_solver()
    m = build(transport_model='SKK')
    ref_coeff_range = np.linspace(0.7, 1, 10)
    frames = []
    A = 1
    B = 0.5
    for idx, ref_coeff in enumerate(ref_coeff_range):
        set_operating_conditions(m, pressure=65e5, tds = 0.07, A_LMH = A, B_LMH = B, reflect_coeff=ref_coeff)
        results = solver.solve(m)
        frames.append(create_report(m))
    report = pd.concat(frames)
    report.to_csv(os.path.join(par_dir, 'reports/', 'SKK_stateblock_reflect_sensitivity_simple.csv'))
    print(report)
    print(report[['Feed Pressure',"Feed Osm Pressure",'Water Flux [LMH]', 'Salt Flux [LMH]', 'Reflection Coeff.', 'Rejection_Obs']])
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_data(report['Reflection Coeff.'], report['Water Flux [LMH]'], report['Salt Flux [LMH]'], rej_obs=report['Rejection_Obs'].to_list(), rej_skk=report['Rejection'].to_list(), ax=ax, xlabel=r'Reflection Coeff. $(\sigma)$', ylabel=r'Water Flux $(\frac{L}{m^{2}  hr})$', ylabel2=r'Salt Flux $(\frac{kg}{m^{2}  hr})$', A_LMH=A, B_LMH=B)
    fig.tight_layout()
    # plt.savefig(os.path.join(par_dir, 'reports/', 'SKK_stateblock_reflect_sensitivity_simple.png'), dpi=300)
    plt.show()

def main():
    solver = get_solver()
    m = build(transport_model='SKK')
    set_operating_conditions(m, pressure=70e5, tds = 0.035, A_LMH = 1, B_LMH = 0.5, reflect_coeff=0.85)
    # set_operating_conditions(m, pressure=65e5, tds = 0.1, A_LMH = 1, B_LMH = 1, reflect_coeff=0.99)
    results = solver.solve(m)
    df = create_report(m)
    df.to_csv(os.path.join(par_dir, 'reports/', 'SKK_stateblock_report.csv'))
    print(df)

if __name__ == "__main__":
    # main()
    # param_sweep_2D()
    # reflection_sensitivity()
    simple_reflection_sensitivity()
