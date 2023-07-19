from pyomo.environ import (
    value,
    units as pyunits,
)
import pandas as pd

liq = 'Liq'
h20 = 'H2O'
nacl = 'NaCl'
def report_inlet_condition(m,stage):
    print('\n\n====================', f'Stage = {stage}', '====================\n\n')
    # m.fs.OAROUnits[stage].display()
    print('Feed:')
    print(f'{"Pressure":<20s}  {f"{value(pyunits.convert(m.fs.OAROUnits[stage].feed_inlet.pressure[0], to_units=pyunits.bar)):<10,.1f}"} {"bar":<10s}')
    print(f'{"Flow H20":<20s}  {f"{value(m.fs.OAROUnits[stage].feed_side.properties_in[0.0].flow_mass_phase_comp[liq, h20]):<10,.1f}"} {"kg/s":<10s}')
    print(f'{"Flow NaCl":<20s}  {f"{value(m.fs.OAROUnits[stage].feed_side.properties_in[0.0].flow_mass_phase_comp[liq, nacl]):<10,.2f}"} {"kg/s":<10s}')
    print(f'{"Mass Frac.":<20s}  {f"{(100*value(m.fs.OAROUnits[stage].feed_side.properties_in[0.0].flow_mass_phase_comp[liq, nacl]/m.fs.OAROUnits[stage].feed_side.properties_in[0.0].flow_mass_phase_comp[liq, h20])):<10,.2f}"} {"%":<10s}')
    # print(f'{"Feed Reynolds":<20s}  {f"{value(m.fs.OAROUnits[stage].permeate_side.dh()):<10,.1f}"} {"m**2":<10s}')
    # print(f'{"Feed Reynolds":<20s}  {f"{value(m.fs.OAROUnits[stage].feed_side.deltaP):<10,.1f}"} {"kg/m/s**2":<10s}')
    print('\nPermeate:')
    print(f'{"Pressure":<20s}  {f"{value(pyunits.convert(m.fs.OAROUnits[stage].permeate_inlet.pressure[0], to_units=pyunits.bar)):<10,.1f}"} {"bar":<10s}')
    print(f'{"Flow H20":<20s}  {f"{value(m.fs.OAROUnits[stage].permeate_side.properties_in[0.0].flow_mass_phase_comp[liq, h20]):<10,.1f}"} {"kg/s":<10s}')
    print(f'{"Flow NaCl":<20s}  {f"{value(m.fs.OAROUnits[stage].permeate_side.properties_in[0.0].flow_mass_phase_comp[liq, nacl]):<10,.2f}"} {"kg/s":<10s}')
    print(f'{"Mass Frac":<20s}  {f"{100*value(m.fs.OAROUnits[stage].permeate_side.properties_in[0.0].flow_mass_phase_comp[liq, nacl]/m.fs.OAROUnits[stage].permeate_side.properties_in[0.0].flow_mass_phase_comp[liq, h20]):<10,.2f}"} {"%":<10s}')
    print('\n')
    # Other Possible Outputs
    # m.fs.OAROUnits[stage].length
    # m.fs.OAROUnits[stage].width
    

def report_outlet_condition(m,stage):
    # m.fs.OAROUnits[stage].display()
    print('\nMembrane:')
    print(f'{"Mem. Length":<20s}  {f"{value(m.fs.OAROUnits[stage].length):<10,.1f}"} {"m":<10s}')
    print(f'{"Mem. Width":<20s}  {f"{value(m.fs.OAROUnits[stage].width):<10,.1f}"} {"m":<10s}')
    print(f'{"Mem. Area":<20s}  {f"{value(m.fs.OAROUnits[stage].area):<10,.1f}"} {"m":<10s}')
    print(f'{"Water Flux":<20s}  {f"{value(pyunits.convert(m.fs.OAROUnits[stage].flux_mass_phase_comp[0.0, 0.0, liq, h20], to_units=pyunits.kg / pyunits.m**2 / pyunits.hour)):<10,.1f}"} {"LMH":<10s}')
    print(f'{"Salt Flux":<20s}  {f"{value(pyunits.convert(m.fs.OAROUnits[stage].flux_mass_phase_comp[0.0, 0.0, liq, nacl], to_units=pyunits.gram / pyunits.m**2 / pyunits.hour)):<10,.1f}"} {"kg/m**2/hr":<10s}')
    print(f'{"Water Flux":<20s}  {f"{value(pyunits.convert(m.fs.OAROUnits[stage].flux_mass_phase_comp[0.0, 1.0, liq, h20], to_units=pyunits.kg / pyunits.m**2 / pyunits.hour)):<10,.1f}"} {"LMH":<10s}')
    print(f'{"Salt Flux":<20s}  {f"{value(pyunits.convert(m.fs.OAROUnits[stage].flux_mass_phase_comp[0.0, 1.0, liq, nacl], to_units=pyunits.gram / pyunits.m**2 / pyunits.hour)):<10,.1f}"} {"kg/m**2/hr":<10s}')


    print('\nFeed:')
    print(f'{"Pressure Loss":<20s}  {f"{-1*value(pyunits.convert(m.fs.OAROUnits[stage].feed_side.deltaP[0], to_units=pyunits.bar)):<10,.2f}"} {"bar":<10s}')
    print(f'{"Velocity In":<20s}  {f"{value(m.fs.OAROUnits[stage].feed_side.velocity[0,0]):<10,.2f}"} {"m/s":<10s}')
    print(f'{"Reynolds In":<20s}  {f"{value(m.fs.OAROUnits[stage].feed_side.N_Re[0,0]):<10,.1f}"} {"dimless":<10s}')
    print(f'{"Reynolds Out":<20s}  {f"{value(m.fs.OAROUnits[stage].feed_side.N_Re[0,1]):<10,.1f}"} {"dimless":<10s}')
    print(f'{"Hydr. Dia.":<20s}  {f"{value(pyunits.convert(m.fs.OAROUnits[stage].feed_side.dh, to_units=pyunits.mm)):<10,.1f}"} {"mm":<10s}')
    print(f'{"Conc. In":<20s}  {f"{value(m.fs.OAROUnits[stage].feed_side.properties_in[0.0].conc_mass_phase_comp[liq, nacl]):<10,.1f}"} {"kg/m**3":<10s}')
    print(f'{"Conc. Out":<20s}  {f"{value(m.fs.OAROUnits[stage].feed_side.properties_out[0.0].conc_mass_phase_comp[liq, nacl]):<10,.1f}"} {"kg/m**3":<10s}')

    print('\nPermeate:')
    print(f'{"Pressure Loss":<20s}  {f"{-1*value(pyunits.convert(m.fs.OAROUnits[stage].permeate_side.deltaP[0], to_units=pyunits.bar)):<10,.2f}"} {"bar":<10s}')
    print(f'{"Velocity In":<20s}  {f"{value(m.fs.OAROUnits[stage].permeate_side.velocity[0,0]):<10,.2f}"} {"m/s":<10s}')
    print(f'{"Reynolds In":<20s}  {f"{value(m.fs.OAROUnits[stage].permeate_side.N_Re[0,0]):<10,.1f}"} {"dimless":<10s}')
    print(f'{"Reynolds Out":<20s}  {f"{value(m.fs.OAROUnits[stage].permeate_side.N_Re[0,1]):<10,.1f}"} {"dimless":<10s}')
    print(f'{"Hydr. Dia.":<20s}  {f"{value(pyunits.convert(m.fs.OAROUnits[stage].permeate_side.dh, to_units=pyunits.mm)):<10,.1f}"} {"mm":<10s}')
    print(f'{"Conc. In":<20s}  {f"{value(m.fs.OAROUnits[stage].permeate_side.properties_in[0.0].conc_mass_phase_comp[liq, nacl]):<10,.1f}"} {"kg/m**3":<10s}')
    print(f'{"Conc. Out":<20s}  {f"{value(m.fs.OAROUnits[stage].permeate_side.properties_out[0.0].conc_mass_phase_comp[liq, nacl]):<10,.1f}"} {"kg/m**3":<10s}')
    print('\n\n====================', '=========', '====================\n\n')

def create_system_report(m):
    report = {}
    frames = []
    def get_state(blk):
        stage_report = {}
        if blk == m.fs.RO:
            # m.fs.RO.display()
            stage_report['Mem Area'] = blk.area.value
            stage_report['Mem Length'] = blk.length.value
            stage_report['Mem Width'] = blk.width.value
            stage_report['Water Flux'] = value(pyunits.convert(blk.flux_mass_phase_comp[0.0, 1.0, liq, h20], to_units=pyunits.kg / pyunits.m**2 / pyunits.hour))
            stage_report['Salt Flux'] = value(pyunits.convert(blk.flux_mass_phase_comp[0.0, 1.0, liq, nacl], to_units=pyunits.kg / pyunits.m**2 / pyunits.hour))
            stage_report['Feed Pressure'] = value(pyunits.convert(blk.feed_side.properties_in[0.0].pressure, to_units=pyunits.bar))
            stage_report['Feed Q In'] = value(pyunits.convert(blk.feed_side.properties_in[0.0].flow_mass_phase_comp[liq, h20], to_units=pyunits.kg / pyunits.min))
            stage_report['Feed Q Out'] = value(pyunits.convert(blk.feed_side.properties_out[0.0].flow_mass_phase_comp[liq, h20], to_units=pyunits.kg / pyunits.min))
            stage_report['Perm Q Out'] = value(pyunits.convert(blk.mixed_permeate[0.0].flow_mass_phase_comp[liq, h20], to_units=pyunits.kg / pyunits.min))
            stage_report['Feed Conc'] = value(pyunits.convert(blk.feed_side.properties_in[0.0].conc_mass_phase_comp[liq, nacl], to_units=pyunits.gram / pyunits.liter))
            stage_report['Perm Conc'] = value(pyunits.convert(blk.mixed_permeate[0.0].conc_mass_phase_comp[liq, nacl], to_units=pyunits.gram / pyunits.liter))
            stage_report['Feed P Loss'] = value(pyunits.convert(blk.feed_side.deltaP[0], to_units=pyunits.bar))
            stage_report['Feed MTC'] = value((blk.feed_side.K[0.0, 0.0, nacl] + blk.feed_side.K[0.0, 0.0, nacl])/2  * 1000 * 3600)

        else:
            stage_report['Mem Area'] = blk.area.value
            stage_report['Mem Length'] = blk.length.value
            stage_report['Mem Width'] = blk.width.value
            stage_report['Water Flux'] = value(pyunits.convert(blk.flux_mass_phase_comp[0.0, 1.0, liq, h20], to_units=pyunits.kg / pyunits.m**2 / pyunits.hour))
            stage_report['Salt Flux'] = value(pyunits.convert(blk.flux_mass_phase_comp[0.0, 1.0, liq, nacl], to_units=pyunits.kg / pyunits.m**2 / pyunits.hour))
            stage_report['Feed Pressure'] = value(pyunits.convert(blk.feed_inlet.pressure[0], to_units=pyunits.bar))
            stage_report['Perm Pressure'] = value(pyunits.convert(blk.permeate_inlet.pressure[0], to_units=pyunits.bar))
            stage_report['Feed Q In'] = value(pyunits.convert(blk.feed_side.properties_in[0.0].flow_mass_phase_comp[liq, h20], to_units=pyunits.kg / pyunits.min))
            stage_report['Feed Q Out'] = value(pyunits.convert(blk.feed_side.properties_out[0.0].flow_mass_phase_comp[liq, h20], to_units=pyunits.kg / pyunits.min))
            stage_report['Perm Q In'] = value(pyunits.convert(blk.permeate_side.properties_in[0.0].flow_mass_phase_comp[liq, h20], to_units=pyunits.kg / pyunits.min))
            stage_report['Perm Q Out'] = value(pyunits.convert(blk.permeate_side.properties_out[0.0].flow_mass_phase_comp[liq, h20], to_units=pyunits.kg / pyunits.min))
            stage_report['Feed Conc'] = value(pyunits.convert(blk.feed_side.properties_in[0.0].conc_mass_phase_comp[liq, nacl], to_units=pyunits.gram / pyunits.liter))
            stage_report['Perm Conc'] = value(pyunits.convert(blk.permeate_side.properties_in[0.0].conc_mass_phase_comp[liq, nacl], to_units=pyunits.gram / pyunits.liter))
            stage_report['Feed P Loss'] = value(pyunits.convert(blk.feed_side.deltaP[0], to_units=pyunits.bar))
            stage_report['Perm P Loss'] = value(pyunits.convert(blk.permeate_side.deltaP[0], to_units=pyunits.bar))
            stage_report['Feed MTC'] = value(blk.feed_side.K_avg[0, nacl]  * 1000 * 3600)
            stage_report['Perm MTC'] = value(blk.permeate_side.K_avg[0, nacl]  * 1000 * 3600)

        return stage_report

    for idx, stage in enumerate(m.fs.NonFinalStages):
        print('Generating report for stage', stage)
        report = get_state(m.fs.OAROUnits[stage])
        frames.append(pd.DataFrame(report, index=['Stage '+str(stage)]))
        
    report = get_state(m.fs.RO)
    frames.append(pd.DataFrame(report, index=['RO']))
    df = pd.concat(frames)
    print(df)
    df.to_csv('watertap/examples/flowsheets/oaro/stage_breakdown/breakdown_'+str(value(m.fs.NumberOfStages))+'_stages.csv')
