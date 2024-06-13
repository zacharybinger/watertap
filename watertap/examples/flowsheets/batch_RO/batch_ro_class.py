# General python imports
import os
import numpy as np
import pandas as pd
import logging
from collections import deque
from os.path import join, dirname
# Pyomo imports
from pyomo.environ import Set, Expression, value, Objective
import datetime
# IDAES imports
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel
from idaes.core.solvers.get_solver import get_solver
import idaes.logger as idaeslog
# Flowsheet function imports
from analysisWaterTAP.flowsheets.FFR_RO.multiperiod.ffrro_config_mp import (
    fix_dof_and_initialize,
)
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.ticker import PercentFormatter
import matplotlib.dates as mdates
# import seaborn as sns
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate
from pyomo.environ import Var, value, units as pyunits
from analysisWaterTAP.flowsheets.FFR_RO.ffrro_system import *
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
    return []

def unfix_dof(
        m,
        water_recovery=None,
        fixed_pressure=None,
        Qin = None,
        Qout = None,
        Qprod=None,
        mem_width = None,
        mem_length = None,
        velocity = None,
        objective=None):
    """
    This function unfixes a few degrees of freedom for optimization

    Args:
        feed_flow_rate: 

    Returns:
        None
    """

    if objective == 'LCOW':
        m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOW)

    print("\n\nDOF before optimization: ", degrees_of_freedom(m))
    if water_recovery is not None:
        print(f"\n------- Fixed Recovery at {100*water_recovery}% -------")
        # product mass flow rate fraction of feed [-]
        m.fs.water_recovery.fix(water_recovery)
        m.fs.feed_flow_mass.fix(Qprod/water_recovery)
        print(f"------- Fixed Feed Flow Rate at {Qprod/water_recovery:<3.1f} kg/s -------")
    else:
        m.fs.water_recovery.unfix()
        m.fs.water_recovery.setlb(0.01)
        m.fs.water_recovery.setub(0.99)

    if Qprod is not None:
        print(f"------- Fixed Product Flow at {Qprod} kg/s -------")
        m.fs.perm_flow_mass.fix(Qprod)
        m.fs.feed_flow_mass.unfix()

    if fixed_pressure is not None:
        print(f"\n------- Fixed Pressure at {fixed_pressure} -------\n")
        m.fs.P1.control_volume.properties_out[0].pressure.fix(fixed_pressure)
    else:
        print(f"------- Unfixed Pressure -------")
        m.fs.P1.control_volume.properties_out[0].pressure.unfix()
        m.fs.P2.control_volume.properties_out[0].pressure.unfix()

def add_costing_constraints(mp):
    n_time_points = len(mp.blocks)

    # Total cost
    @mp.Expression(doc="total cost")
    def adj_LCOW(b):
        # The annualized capital cost is evenly distributed to the multiperiod
        return (b.blocks[0].process.fs.costing.LCOW + b.blocks[1].process.fs.costing.LCOW)/2
    
    @mp.Expression(doc="total cost")
    def LCOW(b):
        return b.adj_LCOW
    
    # Set objective
    mp.obj = Objective(expr=mp.LCOW)

def add_operating_constraints(mp):
    n_time_points = len(mp.blocks)
    for t in range(n_time_points):
        # Set operating constraints
        # mp.blocks[t].process.fs.P1.control_volume.properties_in[0].pressure.fix(8e5)
        # mp.blocks[t].process.fs.P2.control_volume.properties_in[0].pressure.fix(8e5)
        # mp.blocks[t].process.fs.feed_flow_mass.fix(77.77)
        # mp.blocks[t].process.fs.feed_concentration.fix(0.95)
        mp.blocks[t].process.fs.water_recovery.fix(mp.water_recovery)
        mp.blocks[t].process.fs.duty_cycle.fix(mp.duty_cycle)

def create_multiperiod_ffrro_model(
        n_time_points=2,
        recovery=0.8,
        duty_cycle=0.99,
):
    """
    This function creates a multi-period pv battery flowsheet object. This object contains 
    a pyomo model with a block for each time instance.

    Args:
        n_time_points: Number of time blocks to create

    Returns:
        Object containing multi-period vagmd batch flowsheet model
    """
    mp = MultiPeriodModel(
        n_time_points=n_time_points,
        process_model_func=build_system,
        linking_variable_func=get_variable_pairs,
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
        outlvl=logging.WARNING,
    )

    modes = ['System A', 'System B']
    flowsheet_options={ t: {"configuration":modes[t],} 
                            for t in range(n_time_points)
    }

    mp.build_multi_period_model(
        model_data_kwargs=flowsheet_options,
        flowsheet_options=flowsheet_options,
        initialization_options=None,
        unfix_dof_options={
            "water_recovery": recovery,
            "fixed_pressure":None, 
            "Qprod":77.77, 
            "objective":None}
        )
    
    mp.water_recovery = Var(
        initialize=0.8,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )
    
    mp.duty_cycle = Var(
        initialize=duty_cycle,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="3rd Stage Duty Cycle",
    )

    # @mp.Constraint(doc="Water Recovery")
    # def eq_water_recovery(b):
    #     return (b.total_capital_cost / 20)

    active_blks = mp.get_active_process_blocks()
    # Initialize and unfix dof for each period
    solver = get_solver()
    for blk in active_blks:
        fix_dof_and_initialize(
            m=blk
        )
        result = solver.solve(blk)
        unfix_dof(m=blk, water_recovery=value(mp.water_recovery), fixed_pressure=None, Qprod = 77.77, mem_length=7.0)
    
    add_costing_constraints(mp)
    add_operating_constraints(mp)

    return mp
    
def solve(model, solver=None, tee=False, raise_on_failure=False):
    # ---solving---
    if solver is None:
        solver = get_solver()

    print("\n--------- SOLVING ---------\n")
    add_operating_constraints(model)
    results = solver.solve(model, tee=tee)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        return results
    msg = (
        "The current configuration is infeasible. Please adjust the decision variables."
    )
    if raise_on_failure:
        debug(model, solver=solver, automate_rescale=False, resolve=False)
        # debug(model, solver=solver, automate_rescale=False, resolve=False)
        # check_jac(model)
        raise RuntimeError(msg)
    else:
        print(msg)
        # debug(model, solver=solver, automate_rescale=False, resolve=False)
        # check_jac(model)
        return results

if __name__ == "__main__":
    mp = create_multiperiod_ffrro_model()
    mp.water_recovery.fix(0.86)
    mp.duty_cycle.fix(0.95)
    results = solve(mp)
    print("Solve Complete")

    print(f'\n\n{"Normal Operation:":<30}')
    display_flow_table(mp.blocks[0].process)
    # display_system_metrics(mp.blocks[0].process)
    
    print(f'\n\n{"Split Flow Operation:":<30}')
    display_flow_table(mp.blocks[1].process)
    # display_system_metrics(mp.blocks[1].process)