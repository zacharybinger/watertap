# Imports from Pyomo, including "value" for getting the 
# value of Pyomo objects
from pyomo.environ import ConcreteModel, Objective, Expression, value, units as pyunits
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

salt_perm_mag = -7

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
    
    set_scaling(m)

    return m

def set_operating_conditions(m, 
                             tds=0.035,
                             pressure = 50e5,
                             A_LMH = 10,
                             B_LMH = 1,
                             reflect_coeff=0.9):

    for blk in [m.fs.ROunit_SD_0D, m.fs.ROunit_SD_1D, m.fs.ROunit_SKK_0D, m.fs.ROunit_SKK_1D]:
        blk.inlet.pressure[0].fix(pressure)                              # feed pressure (Pa)
        blk.inlet.temperature[0].fix(298.15)                         # feed temperature (K)
        blk.A_comp.fix(A_LMH/(1000*3600*1e5))                             # membrane water permeability (m/Pa/s)
        blk.B_comp.fix(B_LMH/(1000*3600))  
        blk.permeate.pressure[0].fix(101325)                         # permeate pressure (Pa)                               # membrane salt permeability (m/s)
        blk.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(tds)  # mass flow rate of NaCl (kg/s)
        blk.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(1-tds)   # mass flow rate of water (kg/s)
        blk.area.fix(1)                                             # membrane area (m^2)



        if (blk.name == 'fs.ROunit_SKK_1D') | (blk.name == 'fs.ROunit_SD_1D'): 
            blk.width.fix(1)
            blk.feed_side.channel_height.fix(0.001)
            blk.feed_side.spacer_porosity.fix(0.97)

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
    for blk in [m.fs.ROunit_SD_0D, m.fs.ROunit_SD_1D, m.fs.ROunit_SKK_0D, m.fs.ROunit_SKK_1D]:
        blk.initialize()

    assert_no_degrees_of_freedom(m)

def solve(m, solver=None, tee=True):
    if solver is None:
        solver = get_solver()

    results = solver.solve(m, tee=False)

    return results

def reflection_sensitivity(m):
    solver = get_solver()
    
    # fig, ax = plt.subplots(figsize=(8, 6))
    ref_coeff_range = np.linspace(0.5, 1, 10)
    alpha_range = np.linspace(1.4E6, 2.6E6, 10)
    salt_perm_range = np.linspace(1E-7, 1E-5, 10)
    flux_results = []
    salt_flux_results = []
    rejection_results = []
    rej_diff = np.zeros((len(ref_coeff_range), len(salt_perm_range)))
    reports = []
    m.fs.ROunit_SKK.alpha.unfix()
    for idx1, ref_coeff in enumerate(ref_coeff_range):
        for idx2, salt_perm in enumerate(salt_perm_range):
            m.fs.ROunit_SKK.reflect_coeff.fix(ref_coeff)
            m.fs.ROunit_SKK.B_comp.fix(salt_perm)
            m.fs.ROunit_SD.B_comp.fix(salt_perm)
            results = solver.solve(m)
            # print_results(m)
            skk_rejection  = (1 - value(m.fs.ROunit_SKK.mixed_permeate[0.0].flow_mass_phase_comp['Liq', 'NaCl'] / 
                                m.fs.ROunit_SKK.feed_side.properties_in[0.0].flow_mass_phase_comp['Liq', 'NaCl']))
            sd_rejection  = (1 - value(m.fs.ROunit_SD.mixed_permeate[0.0].flow_mass_phase_comp['Liq', 'NaCl'] / 
                                m.fs.ROunit_SD.feed_side.properties_in[0.0].flow_mass_phase_comp['Liq', 'NaCl']))
            
            flux_results.append(value(pyunits.convert(m.fs.ROunit_SKK.flux_mass_phase_comp_avg[0, 'Liq', 'H2O'] / (1000 * pyunits.kg / pyunits.m ** 3),
                                to_units=pyunits.mm / pyunits.hr)))
            salt_flux_results.append(value(pyunits.convert(m.fs.ROunit_SKK.flux_mass_phase_comp_avg[0, 'Liq', 'NaCl'], 
                                    to_units=pyunits.kg / pyunits.m ** 2 / pyunits.hr)))
            rejection_results.append(100* (1 - value(m.fs.ROunit_SKK.mixed_permeate[0.0].flow_mass_phase_comp['Liq', 'NaCl'] / 
                                m.fs.ROunit_SKK.feed_side.properties_in[0.0].flow_mass_phase_comp['Liq', 'NaCl']
                        )
                        ))
            
            rej_diff[idx2,idx1]=(abs(skk_rejection - sd_rejection)/sd_rejection)*100

    low = min(rej_diff.flatten())
    high = max(rej_diff.flatten())
    mid = (low + high)/2
    levels = np.linspace(low, high, 5).tolist()
    levels = [0.1, 0.5, 1, 5, 10]
    print(levels)
    plot_contour(ref_coeff_range, salt_perm_range, rej_diff, levels=levels, x_label=r'Reflection Coeff. $(\sigma)$', y_label=r'Solute Permeability B $(m/s)$', z_label='SD vs SKK \nRejection Difference $(\%)$', low=low, mid=mid, high=high)
    # plot_data(alpha_range, flux_results, salt_flux_results, rejection_results, ax=ax, xlabel=r'Alpha $(\frac{1-\sigma}{B})$', ylabel=r'Water Flux $(\frac{L}{m^{2}  hr})$', ylabel2=r'Salt Flux $(\frac{kg}{m^{2}  hr})$')
    # plot_data(ref_coeff_range, flux_results, salt_flux_results, rejection_results, ax=ax, xlabel=r'Reflection Coeff. $(\sigma)$', ylabel=r'Water Flux $(\frac{L}{m^{2}  hr})$', ylabel2=r'Salt Flux $(\frac{kg}{m^{2}  hr})$')
    # fig.tight_layout()
    # plt.savefig('SKK_RO_reflect_vs_B.png', dpi=300)
    plt.show()
    return flux_results, salt_flux_results, rejection_results, rej_diff

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
    set_operating_conditions(m)
    # set_scaling(m)
    initialize(m)
    results = solve(m)

    # flux, salt_flux, rejection, rej_diff = reflection_sensitivity(m)
    rej_diff, skk_B_calc, sd_rej, skk_rej = fixed_alpha_sensitivity(m)
    
    print(rej_diff, skk_rej, sd_rej)

def main():
    m = build()
    set_operating_conditions(m, pressure=70e5, tds = 0.035, A_LMH = 1, B_LMH = 1, reflect_coeff=0.9)
    
    initialize(m)
    results = solve(m)

    # print(m.fs.ROunit_SKK_1D.display())

    sd_0D_df = SKK_RO_report(m.fs.ROunit_SD_0D)
    sd_1D_df = SKK_RO_report(m.fs.ROunit_SD_1D)
    skk_0D_df =SKK_RO_report(m.fs.ROunit_SKK_0D)
    skk_1D_df =SKK_RO_report(m.fs.ROunit_SKK_1D)
   
    report = pd.concat([sd_0D_df, sd_1D_df, skk_0D_df, skk_1D_df])
    report.to_csv('SKK_RO_report.csv')
    print(report)

if __name__ == "__main__":
    main()
    # sensitivity_study()
    # concentration_sensitivity()
