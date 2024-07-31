# General python imports
import os
import numpy as np
import pandas as pd
import logging
from os.path import join, dirname
# Pyomo imports
from pyomo.environ import Set, Expression, value, Objective
# IDAES imports
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel
from idaes.core.solvers.get_solver import get_solver
import idaes.logger as idaeslog
# Flowsheet function imports

from pyomo.environ import Var, value, units as pyunits
from watertap.examples.flowsheets.batch_RO.batch_ro import *
from watertap.examples.flowsheets.batch_RO.batch_ro_config import *

solver = get_solver()


def get_variable_pairs(t1, t2):
    """
    This function returns pairs of variables that need to be connected across two time periods

    Args:
        t1: current time block
        t2: next time block

    Returns:
        None
    """
    # return []
    return [
        (t2.fs.recirc.properties[0].flow_mass_phase_comp["Liq", "H2O"], t1.fs.disposal.properties[0].flow_mass_phase_comp["Liq", "H2O"]),
        (t2.fs.recirc.properties[0].flow_mass_phase_comp["Liq", "TDS"], t1.fs.disposal.properties[0].flow_mass_phase_comp["Liq", "TDS"]),
        (t2.fs.recirc.properties[0].pressure, t1.fs.disposal.properties[0].pressure),
        (t2.fs.recirc.properties[0].temperature, t1.fs.disposal.properties[0].temperature),
    ]

def unfix_dof(
        m,
        water_recovery=None,
        Q_ro=None,
        time_idx = 0):
    """
    This function unfixes a few degrees of freedom for optimization

    Args:
        water_recovery: 
        Q_ro:

    Returns:
        None
    """

    optimize(m, water_recovery=water_recovery, Q_ro = Q_ro)

    if time_idx > 0:
        m.fs.recirc.properties[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
        m.fs.recirc.properties[0].flow_mass_phase_comp["Liq", "TDS"].unfix()
        m.fs.recirc.properties[0].pressure.unfix()
        m.fs.recirc.properties[0].temperature.unfix()


# def add_operating_constraints(mp):
#     n_time_points = len(mp.blocks)
    # for t in range(n_time_points):
    #     # Set operating constraints
    #     # mp.blocks[t].process.fs.P1.control_volume.properties_in[0].pressure.fix(8e5)
    #     # mp.blocks[t].process.fs.P2.control_volume.properties_in[0].pressure.fix(8e5)
    #     # mp.blocks[t].process.fs.feed_flow_mass.fix(77.77)
    #     # mp.blocks[t].process.fs.feed_concentration.fix(0.95)
    #     mp.blocks[t].process.fs.water_recovery.fix(mp.water_recovery)
    #     mp.blocks[t].process.fs.duty_cycle.fix(mp.duty_cycle)

def create_multiperiod_ffrro_model(
        n_time_points=2,
        recovery=0.1,
):
    """
    This function creates a multi-period pv battery flowsheet object. This object contains 
    a pyomo model with a block for each time instance.

    Args:
        n_time_points: Number of time blocks to create

    Returns:
        Object containing multi-period vagmd batch flowsheet model
    """
    watertap_solver = get_solver()

    mp = MultiPeriodModel(
        n_time_points=n_time_points,
        process_model_func=build_system,
        linking_variable_func=get_variable_pairs,
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
        solver=watertap_solver,
        outlvl=logging.WARNING,
    )

    # mp.build_multi_period_model()

    # modes = ['System A', 'System B']
    # flowsheet_options={ t: {"configuration":modes[t],} 
    #                         for t in range(n_time_points)
    # }

    flowsheet_options={ t: {"time_blk":t,} 
                            for t in range(n_time_points)
    }

    mp.build_multi_period_model(
        model_data_kwargs=flowsheet_options,
        unfix_dof_options={'water_recovery':0.1, 'Q_ro':0.965})
    

    active_blks = mp.get_active_process_blocks()
    print(active_blks)
    # Initialize and unfix dof for each period
    # solver = get_solver()
    for idx, blk in enumerate(active_blks):
        fix_dof_and_initialize(
            m=blk
        )
        result = solve(blk, solver=watertap_solver, tee=False, raise_on_failure=True)
        unfix_dof(m=blk, water_recovery=recovery, Q_ro=0.965, time_idx=idx)
    
    # add_costing_constraints(mp)
    # add_operating_constraints(mp)

    return mp
    
def solve(model, solver=None, tee=False, raise_on_failure=False):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")
    print(f'Degrees of Freedom: {degrees_of_freedom(model)}')
    # add_operating_constraints(model)
    results = solver.solve(model, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        print("\n--------- INFEASIBLE SOLVE!!! ---------\n")

        print("\n--------- CLOSE TO BOUNDS ---------\n")
        print_close_to_bounds(model)

        print("\n--------- INFEASIBLE BOUNDS ---------\n")
        print_infeasible_bounds(model)

        print("\n--------- INFEASIBLE CONSTRAINTS ---------\n")
        print_infeasible_constraints(model)

        raise RuntimeError(msg)
    else:
        print(msg)
        # debug(model, solver=solver, automate_rescale=False, resolve=False)
        # check_jac(model)
        assert False

if __name__ == "__main__":
    mp = create_multiperiod_ffrro_model(n_time_points=10)
    results = solve(mp, raise_on_failure=True)
    print("Solve Complete")

    for blk in mp.get_active_process_blocks():
        print(blk.fs.M1.report())

    ro_feed_salinity = [value(blk.fs.M1.mixed_state[0.0].flow_mass_phase_comp["Liq", "TDS"]) for blk in mp.get_active_process_blocks()]
    ro_feed_conc = [value(blk.fs.M1.mixed_state[0.0].conc_mass_phase_comp["Liq", "TDS"]) for blk in mp.get_active_process_blocks()]
    print(ro_feed_salinity)
    print(ro_feed_conc)