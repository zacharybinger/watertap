#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import os

from pyomo.environ import (
    units as pyunits,
    check_optimal_termination,
    value,
    Expression,
    Param,
)

from watertap.tools.parameter_sweep import LinearSample, parameter_sweep
from watertap.examples.flowsheets.oaro import oaro_multi as oaro
from idaes.core.util.misc import StrEnum


class ERDtype(StrEnum):
    pump_as_turbine = "pump_as_turbine"


def erd_type_not_found(erd_type):
    raise NotImplementedError(
        "erd_type was {}, but can only " "be pump_as_turbine" "".format(erd_type.value)
    )


def _oaro_presweep(number_of_stages=2):
    m = oaro.build(
        number_of_stages=number_of_stages,
        erd_type=ERDtype.pump_as_turbine,
    )
    oaro.set_operating_conditions(m, number_of_stages)
    oaro.initialize(m, number_of_stages)
    oaro.solve(m)
    m.fs.feed.flow_mass_phase_comp.unfix()
    m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix()
    m.fs.feed.properties[0].flow_vol_phase["Liq"].fix()
    oaro.optimize_set_up(m, number_of_stages)
    oaro.solve(m)

    return m


def run_case(number_of_stages, nx, output_filename=None):
    """
    Run the parameter sweep tool on the LSRRO flowsheet, sweeping over feed
    concentration from 5 to 250 kg/m^3 and water recovery from 30% to 90%.

    Arguments
    ---------
    number_of_stages (int) : The number of LSRRO Stages (including the initial RO stage).
    nx (int) : The number of points for both feed concentration and water recovery. The
               total number of points swept will be nx^2.
    output_filename (str, optional): The place to write the parameter sweeep results
               csv file. By default it is
               ./param_sweep_output/{number_of_stages}_stage/results_LSRRO.csv

    Returns
    -------
    global_results (numpy array) : The raw values from the parameter sweep
    sweep_params (dict) : The dictionary of samples
    model (Pyomo ConcreteModel) : The LSRRO flowsheet used for parameter sweeps

    """

    if output_filename is None:
        output_filename = (
            f"param_sweep_output/{number_of_stages}_stage/results_OARO.csv"
        )

    sweep_params = {}
    outputs = {}

    m = _oaro_presweep(number_of_stages=number_of_stages)

    # Sweep parameters ------------------------------------------------------------------------

    # sweep_params["Feed Concentration"] = LinearSample(
    #     m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"], 5, 120, nx
    # )
    #
    # sweep_params["Volumetric Recovery Rate"] = LinearSample(
    #     m.fs.water_recovery, 0.3, 0.9, nx
    # )
    sweep_params["Max Pressure for Recycle Pump"] = LinearSample(
        m.fs.recycle_pump_max_pressure, 3e5, 100e5, nx
    )

    # Outputs  -------------------------------------------------------------------------------
    outputs["LCOW"] = m.fs.costing.LCOW
    outputs["LCOW wrt Feed Flow"] = m.fs.costing.LCOW_feed
    outputs["SEC"] = m.fs.costing.specific_energy_consumption
    outputs["SEC wrt Feed"] = m.fs.costing.specific_energy_consumption_feed
    outputs["Number of Stages"] = m.fs.NumberOfStages
    outputs["Final Brine Concentration"] = m.fs.disposal.properties[
        0
    ].conc_mass_phase_comp["Liq", "NaCl"]
    outputs["Final Permeate Concentration (ppm)"] = (
        m.fs.product.properties[0].mass_frac_phase_comp["Liq", "NaCl"] * 1e6
    )
    outputs["Annual Feed Flow"] = m.fs.annual_feed
    outputs["Annual Water Production"] = m.fs.costing.annual_water_production

    outputs["Pump Work In (kW)"] = m.fs.total_pump_work
    outputs["Pump Work Recovered (kW)"] = m.fs.recovered_pump_work
    outputs["Net Pump Work In (kW)"] = m.fs.net_pump_work
    outputs["Energy Recovery (%)"] = (
        -m.fs.recovered_pump_work / m.fs.total_pump_work * 100
    )

    outputs["Mass Water Recovery Rate (%)"] = m.fs.mass_water_recovery * 100
    outputs["System Salt Rejection (%)"] = m.fs.system_salt_rejection * 100

    outputs["Total Membrane Area"] = m.fs.total_membrane_area
    outputs["Total Capex LCOW"] = (
        m.fs.costing.total_capital_cost
        * m.fs.costing.factor_capital_annualization
        / m.fs.costing.annual_water_production
    )
    outputs["Total Opex LCOW"] = (
        m.fs.costing.total_operating_cost / m.fs.costing.annual_water_production
    )

    outputs["Primary Pump Capex LCOW"] = m.fs.costing.primary_pump_capex_lcow
    outputs["Recycle Pump Capex LCOW"] = m.fs.costing.recycle_pump_capex_lcow

    outputs["ERD Capex LCOW"] = m.fs.costing.erd_capex_lcow
    outputs["Membrane Capex LCOW"] = m.fs.costing.membrane_capex_lcow
    outputs["Indirect Capex LCOW"] = m.fs.costing.indirect_capex_lcow
    outputs["Electricity LCOW"] = m.fs.costing.electricity_lcow
    outputs["Membrane Replacement LCOW"] = m.fs.costing.membrane_replacement_lcow
    outputs[
        "Chem-labor-maintenance LCOW"
    ] = m.fs.costing.chemical_labor_maintenance_lcow

    outputs["Pumping Energy Agg LCOW"] = m.fs.costing.pumping_energy_aggregate_lcow
    outputs["Membrane Agg LCOW"] = m.fs.costing.membrane_aggregate_lcow

    outputs.update(
        {
            f"Feed Pressure (bar)-Stage {idx}": pyunits.convert(
                pump.control_volume.properties_out[0].pressure, to_units=pyunits.bar
            )
            for idx, pump in m.fs.PrimaryPumps.items()
        }
    )

    outputs.update(
        {
            f"Membrane Area-Stage {idx}": stage.area
            for idx, stage in m.fs.OAROUnits.items()
        }
    )

    outputs.update(
        {
            f"Observed Rejection (%)-Stage {idx}": stage.rejection_phase_comp[
                0, "Liq", "NaCl"
            ]
            * 100
            for idx, stage in m.fs.OAROUnits.items()
        }
    )

    outputs.update(
        {
            f"Observed Salt Passage (%)-Stage {idx}": (
                1 - stage.rejection_phase_comp[0, "Liq", "NaCl"]
            )
            * 100
            for idx, stage in m.fs.OAROUnits.items()
        }
    )

    outputs.update(
        {
            f"Mass Salt Passage (%)-Stage {idx}": (
                stage.permeate_side.properties[0, 0].flow_mass_phase_comp["Liq", "NaCl"]
                / stage.feed_side.properties[0, 0].flow_mass_phase_comp["Liq", "NaCl"]
            )
            * 100
            for idx, stage in m.fs.OAROUnits.items()
        }
    )

    outputs.update(
        {
            f"Volumetric Module Recovery Rate (%)-Stage {idx}": stage.recovery_vol_phase[
                0, "Liq"
            ]
            * 100
            for idx, stage in m.fs.OAROUnits.items()
        }
    )

    outputs.update(
        {
            f"Mass Water Module Recovery Rate (%)-Stage {idx}": stage.recovery_mass_phase_comp[
                0, "Liq", "H2O"
            ]
            * 100
            for idx, stage in m.fs.OAROUnits.items()
        }
    )

    # outputs.update(
    #     {
    #         f"Volumetric Stage Recovery Rate (%)-Stage {idx}": m.fs.stage_recovery_vol[
    #             idx
    #         ]
    #         * 100
    #         for idx in m.fs.Stages
    #     }
    # )

    # outputs.update(
    #     {
    #         f"Mass Water Stage Recovery Rate (%)-Stage {idx}": m.fs.stage_recovery_mass_H2O[
    #             idx
    #         ]
    #         * 100
    #         for idx in m.fs.Stages
    #     }
    # )

    outputs.update(
        {
            f"A-Value (LMH/bar)-Stage {idx}": stage.A_comp[0, "H2O"] * 3.6e11
            for idx, stage in m.fs.OAROUnits.items()
        }
    )

    outputs.update(
        {
            f"B-Value (LMH)-Stage {idx}": stage.B_comp[0, "NaCl"] * 3.6e6
            for idx, stage in m.fs.OAROUnits.items()
        }
    )

    outputs.update(
        {
            f"Average Water Flux (LMH)-Stage {idx}": stage.flux_mass_phase_comp_avg[
                0, "Liq", "H2O"
            ]
            * 3.6e3
            for idx, stage in m.fs.OAROUnits.items()
        }
    )

    outputs.update(
        {
            f"Average NaCl Flux (GMH)-Stage {idx}": stage.flux_mass_phase_comp_avg[
                0, "Liq", "NaCl"
            ]
            * 3.6e6
            for idx, stage in m.fs.OAROUnits.items()
        }
    )

    # outputs.update(
    #     {
    #         f"Pressure Drop (bar)-Stage {idx}": -pyunits.convert(
    #             stage.deltaP[0], to_units=pyunits.bar
    #         )
    #         for idx, stage in m.fs.OAROUnits.items()
    #     }
    # )
    outputs.update(
        {
            f"Inlet Reynolds Number-Stage {idx}": stage.feed_side.N_Re[0, 0]
            for idx, stage in m.fs.OAROUnits.items()
        }
    )
    outputs.update(
        {
            f"Outlet Reynolds Number-Stage {idx}": stage.feed_side.N_Re[0, 1]
            for idx, stage in m.fs.OAROUnits.items()
        }
    )
    outputs.update(
        {
            f"Inlet Crossflow Velocity-Stage {idx}": stage.feed_side.velocity[0, 0]
            for idx, stage in m.fs.OAROUnits.items()
        }
    )
    outputs.update(
        {
            f"Outlet Crossflow Velocity-Stage {idx}": stage.feed_side.velocity[0, 1]
            for idx, stage in m.fs.OAROUnits.items()
        }
    )

    outputs.update(
        {
            f"Primary Pump SEC-Stage {idx}": pyunits.convert(
                pump.work_mechanical[0], to_units=pyunits.kW
            )
            / pyunits.convert(
                m.fs.product.properties[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.hr,
            )
            for idx, pump in m.fs.PrimaryPumps.items()
        }
    )

    m.fs.zero_expression = Expression(expr=0.0)

    outputs["Recycle Pump SEC-Stage 1"] = m.fs.zero_expression
    outputs.update(
        {
            f"Recycle Pump SEC-Stage {idx}": pyunits.convert(
                pump.work_mechanical[0], to_units=pyunits.kW
            )
            / pyunits.convert(
                m.fs.product.properties[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.hr,
            )
            for idx, pump in m.fs.RecyclePumps.items()
        }
    )
    outputs.update(
        {
            f"ERD SEC-Stage {idx}": pyunits.convert(
                erd.work_mechanical[0], to_units=pyunits.kW
            )
            / pyunits.convert(
                m.fs.product.properties[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.hr,
            )
            for idx, erd in m.fs.EnergyRecoveryDevices.items()
        }
    )

    outputs.update(
        {
            f"Net SEC-Stage {idx}": outputs[f"Primary Pump SEC-Stage {idx}"]
            + outputs[f"Recycle Pump SEC-Stage {idx}"]
            + outputs[f"ERD SEC-Stage {idx}"]
            for idx in m.fs.Stages
        }
    )

    global_results = parameter_sweep(
        m,
        sweep_params,
        outputs,
        csv_results_file_name=output_filename,
        optimize_function=oaro.solve,
        interpolate_nan_outputs=True,
    )

    return global_results, sweep_params, m


if __name__ == "__main__":
    for n in range(2, 4):
        global_results, sweep_params, m = run_case(number_of_stages=n, nx=17)
        print(global_results)
