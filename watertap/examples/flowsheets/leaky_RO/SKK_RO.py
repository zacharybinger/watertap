import os
# Imports from Pyomo, including "value" for getting the 
# value of Pyomo objects
from pyomo.environ import ConcreteModel, Objective, Expression, check_optimal_termination, assert_optimal_termination, value, units as pyunits
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor
from watertap.core.util.infeasible import *
from watertap.core.util.initialization import assert_no_degrees_of_freedom
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.unit_models.reverse_osmosis_0D import (ReverseOsmosis0D,
        ConcentrationPolarizationType, MassTransferCoefficient)
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from skk_util import plot_data, plot_contour, plot_fixed_alpha, SKK_RO_report

par_dir = os.path.dirname(os.path.abspath(__file__))
liq = "Liq"
h2o = "H2O"
nacl = 'NaCl'
def build():
    # Create a Pyomo concrete model, flowsheet, and NaCl property parameter block.
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()

    m.fs.ROunit_SD_0D = ReverseOsmosis0D(
        property_package=m.fs.properties,
        concentration_polarization_type=ConcentrationPolarizationType.none,
        mass_transfer_coefficient=MassTransferCoefficient.none,
        has_pressure_change=False,
        transport_model = 'SD',
        )
    
    m.fs.ROunit_SD_1D = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=10,
        transport_model = 'SD',
        )
    
    m.fs.ROunit_SKK_0D = ReverseOsmosis0D(
        property_package=m.fs.properties,
        concentration_polarization_type=ConcentrationPolarizationType.none,
        mass_transfer_coefficient=MassTransferCoefficient.none,
        has_pressure_change=False,
        transport_model = 'SKK',
        )
    
    m.fs.ROunit_SKK_1D = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=10,
        transport_model = 'SKK'
        )
    
    # m.fs.units = [m.fs.ROunit_SD_0D, m.fs.ROunit_SKK_0D]
    m.fs.units = [m.fs.ROunit_SD_0D, m.fs.ROunit_SKK_0D, m.fs.ROunit_SD_1D, m.fs.ROunit_SKK_1D]

    set_scaling(m)

    return m

def set_operating_conditions(m, 
                             tds=0.035,
                             pressure = 50e5,
                             A_LMH = 10,
                             B_LMH = 1,
                             reflect_coeff=1):

    for blk in m.fs.units:
        blk.inlet.pressure[0].fix(pressure)                              # feed pressure (Pa)
        blk.inlet.temperature[0].fix(298.15)                         # feed temperature (K)
        blk.A_comp.fix(A_LMH/(1000*3600*1e5))                             # membrane water permeability (m/Pa/s)
        blk.B_comp.fix(B_LMH/(1000*3600))  
        blk.permeate.pressure[0].fix(101325)                         # permeate pressure (Pa)                               # membrane salt permeability (m/s)
        blk.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(tds)  # mass flow rate of NaCl (kg/s)
        blk.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(1-tds)   # mass flow rate of water (kg/s)
        blk.area.fix(13.8)                                             # membrane area (m^2)



        if (blk.name == 'fs.ROunit_SKK_1D') | (blk.name == 'fs.ROunit_SD_1D'): 
            blk.length.fix(1)
            blk.feed_side.channel_height.fix(0.001)
            blk.feed_side.spacer_porosity.fix(0.97)

        if blk.config.transport_model == 'SKK':
            blk.reflect_coeff.fix(reflect_coeff)

def optimize(m, reflect_coeff=0.9):
    for blk in m.fs.units:
        if blk.config.transport_model == 'SKK':
            blk.reflect_coeff.fix(reflect_coeff)

def set_scaling(m):
    # Set scaling factors for component mass flowrates.
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))

    # Set scaling factor for membrane area.
    set_scaling_factor(m.fs.ROunit_SD_0D.area, 1e-2)
    set_scaling_factor(m.fs.ROunit_SKK_0D.area, 1e-2)
    set_scaling_factor(m.fs.ROunit_SKK_1D.feed_side.area, 1e-2)
    set_scaling_factor(m.fs.ROunit_SKK_1D.area, 1e-2)
    set_scaling_factor(m.fs.ROunit_SKK_1D.width, 1e-2)
    set_scaling_factor(m.fs.ROunit_SD_1D.feed_side.area, 1e-2)
    set_scaling_factor(m.fs.ROunit_SD_1D.area, 1e-2)
    set_scaling_factor(m.fs.ROunit_SD_1D.width, 1e-2)

    # Calculate scaling factors for all other variables.
    calculate_scaling_factors(m)

def initialize(m):
    for blk in m.fs.units:
        blk.initialize()

    assert_no_degrees_of_freedom(m)

def solve(m, solver=None, tee=True):
    if solver is None:
        solver = get_solver()

    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)

    if check_optimal_termination(results):
        print('Model solved optimally')
        # print_close_to_bounds(m)
        # print_variables_close_to_bounds(m)
    else:
        print_infeasible_bounds(m)
        print_infeasible_constraints(m)

    return results

def reflection_sensitivity(m):
    solver = get_solver()
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ref_coeff_range = np.linspace(0.7, 0.9, 5)
    salt_perm_range = np.logspace(-1,1,5)
    flux_results = []
    salt_flux_results = []
    rejection_results = []
    rej_diff = np.zeros((len(salt_perm_range), len(ref_coeff_range)))
    salt_flux_matrix = np.zeros((len(salt_perm_range), len(ref_coeff_range)))
    reports = []
    m.fs.ROunit_SKK_0D.alpha.unfix()
    for idx1, ref_coeff in enumerate(ref_coeff_range):
        for idx2, salt_perm in enumerate(salt_perm_range):
            # m.fs.ROunit_SKK_0D.reflect_coeff.fix(ref_coeff)
            # m.fs.ROunit_SKK_0D.B_comp.fix(salt_perm/(1000*3600))
            # m.fs.ROunit_SD_0D.B_comp.fix(salt_perm/(1000*3600))
            # A = (salt_perm/0.01331)**(1/3)
            # m.fs.ROunit_SKK_0D.B_comp.fix(A/(1000*3600*1e5))
            # m.fs.ROunit_SD_0D.B_comp.fix(A/(1000*3600*1e5))
            set_operating_conditions(m, pressure=65e5, tds = 0.01, A_LMH = ((salt_perm/0.01331)**(1/3)), B_LMH = salt_perm, reflect_coeff=ref_coeff)
            results = solver.solve(m)
            assert_optimal_termination(results)
            try:
                assert_optimal_termination(results)
            except:
                print_infeasible_bounds(m)
                print_infeasible_constraints(m)

            # print_results(m)
            skk_rejection  = (1 - value(m.fs.ROunit_SKK_0D.mixed_permeate[0.0].flow_mass_phase_comp['Liq', 'NaCl'] / 
                                m.fs.ROunit_SKK_0D.feed_side.properties_in[0.0].flow_mass_phase_comp['Liq', 'NaCl']))
            sd_rejection  = (1 - value(m.fs.ROunit_SD_0D.mixed_permeate[0.0].flow_mass_phase_comp['Liq', 'NaCl'] / 
                                m.fs.ROunit_SD_0D.feed_side.properties_in[0.0].flow_mass_phase_comp['Liq', 'NaCl']))
            
            flux_results.append(value(pyunits.convert(m.fs.ROunit_SKK_0D.flux_mass_phase_comp_avg[0, 'Liq', 'H2O'] / (1000 * pyunits.kg / pyunits.m ** 3),
                                to_units=pyunits.mm / pyunits.hr)))
            salt_flux_results.append(value(pyunits.convert(m.fs.ROunit_SKK_0D.flux_mass_phase_comp_avg[0, 'Liq', 'NaCl'], 
                                    to_units=pyunits.kg / pyunits.m ** 2 / pyunits.hr)))
            rejection_results.append(100* (1 - value(m.fs.ROunit_SKK_0D.mixed_permeate[0.0].flow_mass_phase_comp['Liq', 'NaCl'] / 
                                m.fs.ROunit_SKK_0D.feed_side.properties_in[0.0].flow_mass_phase_comp['Liq', 'NaCl']
                        )
                        ))
            
            rej_diff[idx2,idx1]=(abs(skk_rejection - sd_rejection)/sd_rejection)*100
            salt_flux_matrix[idx2,idx1] = value(pyunits.convert(m.fs.ROunit_SKK_0D.flux_mass_phase_comp_avg[0, 'Liq', 'H2O'] / (1000 * pyunits.kg / pyunits.m ** 3),
                                to_units=pyunits.mm / pyunits.hr))

    low = min(rej_diff.flatten())
    high = max(rej_diff.flatten())
    mid = (low + high)/2
    levels = np.linspace(low, high, 5).tolist()
    levels = [0.1, 0.5, 1, 5]
    plot_contour(ref_coeff_range, salt_perm_range, rej_diff, levels=levels, x_label=r'Reflection Coeff. $(\sigma)$', y_label=r'Solute Permeability B $(LMH)$', z_label='SD vs SKK \nRejection Difference $(\%)$', low=low, mid=mid, high=high)
    plt.savefig(os.path.join(par_dir, 'reports/', 'SKK_RO_reflect_vs_B_2.png'), dpi=300)
    plt.show()

    return flux_results, salt_flux_results, rejection_results, rej_diff

def skk_contour(m):
    solver = get_solver()

    pressure_range = np.linspace(60e5, 70e5, 2)
    feed_salinity_range = np.linspace(0.05, 0.07, 2)
    flux_results = []
    salt_flux_results = []
    rejection_results = []
    rej_diff = np.zeros((len(feed_salinity_range), len(pressure_range)))
    salt_flux_matrix = np.zeros((len(feed_salinity_range), len(pressure_range)))
    for idx1, pressure_val in enumerate(pressure_range):
        for idx2, conc in enumerate(feed_salinity_range):
            set_operating_conditions(m, pressure=pressure_val, tds = conc, A_LMH = 10.1, B_LMH = 13.7, reflect_coeff=1)
            results = solver.solve(m)
            report = pd.concat([SKK_RO_report(unit) for unit in m.fs.units])
            print(report)

def fixed_alpha_sensitivity(m):
    solver = get_solver()
    alpha = 1E6

    rej_diff = []
    skk_B_calc = []
    skk_rej = []
    sd_rej = []
    m.fs.ROunit_SKK.alpha.fix(alpha)
    m.fs.ROunit_SD.B_comp.fix(1/alpha)
    m.fs.ROunit_SKK.B_comp.unfix()

    fig, ax = plt.subplots(figsize=(8, 6))
    ref_coeff_range = np.linspace(0.5, 0.99, 10)
    for idx1, ref_coeff in enumerate(ref_coeff_range):
        m.fs.ROunit_SKK.reflect_coeff.fix(ref_coeff)

        results = solver.solve(m)

        skk_rejection  = (1 - value(m.fs.ROunit_SKK.mixed_permeate[0.0].flow_mass_phase_comp['Liq', 'NaCl'] / 
                                m.fs.ROunit_SKK.feed_side.properties_in[0.0].flow_mass_phase_comp['Liq', 'NaCl']))
        sd_rejection  = (1 - value(m.fs.ROunit_SD.mixed_permeate[0.0].flow_mass_phase_comp['Liq', 'NaCl'] / 
                                m.fs.ROunit_SD.feed_side.properties_in[0.0].flow_mass_phase_comp['Liq', 'NaCl']))
        print(skk_rejection, sd_rejection)

        rej_diff.append(abs(skk_rejection - sd_rejection)/sd_rejection)
        skk_B_calc.append(value(m.fs.ROunit_SKK.B_comp[0, 'NaCl']))

        skk_rej.append(skk_rejection)
        sd_rej.append(sd_rejection)

    plot_fixed_alpha(ref_coeff_range, skk_B_calc, rej_diff, sd_rej, skk_rej, ax=ax, xlabel=r'Reflection Coeff. $(\sigma)$', ylabel=r'Calc. Solute Perm. Coeff. $(m/s)$', ylabel2='SD vs SKK \nRejection Difference $(\%)$')
    fig.tight_layout()
    # plt.savefig('SKK_RO_fixed_alpha_vs_reflect.png', dpi=300)
    plt.show()
    return rej_diff, skk_B_calc, sd_rej, skk_rej
            
def sensitivity_study():
    m = build()
    # set_operating_conditions(m, pressure=65e5, tds = 0.035, A_LMH = 1, B_LMH = 13.7, reflect_coeff=1)
    set_operating_conditions(m, pressure=65e5, tds = 0.1, A_LMH = 10.1, B_LMH = 13.7, reflect_coeff=1)
    initialize(m)
    # optimize(m, reflect_coeff=0.95)
    # results = solve(m)

    flux, salt_flux, rejection, rej_diff = reflection_sensitivity(m)
    # rej_diff, skk_B_calc, sd_rej, skk_rej = fixed_alpha_sensitivity(m)
    # skk_contour(m)

def main():
    m = build()
    set_operating_conditions(m, pressure=65e5, tds = 0.1, A_LMH = 10.1, B_LMH = 13.7, reflect_coeff=1)
    initialize(m)

    optimize(m, reflect_coeff=0.90)
    results = solve(m)
    report = pd.concat([SKK_RO_report(unit) for unit in m.fs.units])
    # print(m.fs.units[0].display())
    # sd_0D_df = SKK_RO_report(m.fs.ROunit_SD_0D)
    # sd_1D_df = SKK_RO_report(m.fs.ROunit_SD_1D)
    # skk_0D_df =SKK_RO_report(m.fs.ROunit_SKK_0D)
    # skk_1D_df =SKK_RO_report(m.fs.ROunit_SKK_1D)
   
    # report = pd.concat([sd_0D_df, skk_0D_df, sd_1D_df, skk_1D_df])
    report.to_csv(os.path.join(par_dir, 'reports/', 'SKK_unit_model_report.csv'))
    print(report)

    # for unit in m.fs.units:
    #     print(unit.name),":"
    #     print(f'{"Jw":<15s}{f"{value(unit.flux_mass_phase_comp_avg[0, liq, h2o]):<10,.2e}"}{f"{pyunits.get_units(unit.flux_mass_phase_comp_avg[0, liq, h2o])}":<10s}')
    #     print(f'{"A":<15s}{f"{value(pyunits.convert(unit.A_comp[0, h2o], to_units=pyunits.m/pyunits.second/pyunits.Pa)):<10,.2e}"}{f"{pyunits.get_units(pyunits.convert(unit.A_comp[0, h2o], to_units=pyunits.m/pyunits.second/pyunits.Pa))}":<10s}')
    #     print(f'{"Delta P":<15s}{f"{value(unit.inlet.pressure[0] - unit.permeate.pressure[0]):<10,.2e}"}{f"{pyunits.get_units(unit.inlet.pressure[0])}":<10s}')
    #     print(f'{"Reflect Coeff":<15s}{f"{value(unit.reflect_coeff):<10,.2e}"}{f"{pyunits.get_units(unit.reflect_coeff)}":<10s}')
    #     print(f'{"Feed Osm P":<15s}{f"{value(unit.feed_side.properties_interface[0,0].pressure_osm_phase[liq]):<10,.2e}"}{f"{pyunits.get_units(unit.feed_side.properties_interface[0,0].pressure_osm_phase[liq])}":<10s}')
    #     print(f'{"Perm Osm P":<15s}{f"{value(unit.permeate_side[0,0].pressure_osm_phase[liq]):<10,.2e}"}{f"{pyunits.get_units(unit.feed_side.properties_interface[0,0].pressure_osm_phase[liq])}":<10s}')
    #     print(f'{"Delta Osm P":<15s}{f"{value(unit.feed_side.properties_interface[0,0].pressure_osm_phase[liq] - unit.permeate_side[0,0].pressure_osm_phase[liq]):<10,.2e}"}{f"{pyunits.get_units(unit.feed_side.properties_interface[0,0].pressure_osm_phase[liq])}":<10s}')
    #     print('\n')

    # for unit in m.fs.units:
    #     print(unit.name),":"
    #     print(f'{"Js":<15s}{f"{value(unit.flux_mass_phase_comp_avg[0, liq, nacl]):<10,.2e}"}{f"{pyunits.get_units(unit.flux_mass_phase_comp_avg[0, liq, nacl])}":<10s}')
    #     print(f'{"B":<15s}{f"{value(pyunits.convert(unit.B_comp[0, nacl], to_units=pyunits.m/pyunits.second)):<10,.2e}"}{f"{pyunits.get_units(pyunits.convert(unit.B_comp[0, nacl], to_units=pyunits.m/pyunits.second))}":<10s}')
    #     print(f'{"Delta C":<15s}{f"{value(unit.feed_side.properties_in[0.0].conc_mass_phase_comp[liq,nacl]):<10,.2f}"}{f"{pyunits.get_units(unit.feed_side.properties_in[0.0].conc_mass_phase_comp[liq,nacl])}":<10s}')
    #     print(f'{"Reflect Coeff":<15s}{f"{value(unit.reflect_coeff):<10,.2e}"}{f"{pyunits.get_units(unit.reflect_coeff)}":<10s}')
    #     print(f'{"Jw":<15s}{f"{value(unit.flux_mass_phase_comp_avg[0, liq, h2o]):<10,.2e}"}{f"{pyunits.get_units(unit.flux_mass_phase_comp_avg[0, liq, h2o])}":<10s}')
    #     print('\n')

if __name__ == "__main__":
    # main()
    sensitivity_study()
