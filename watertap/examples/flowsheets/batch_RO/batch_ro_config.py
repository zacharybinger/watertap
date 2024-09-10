from pyomo.environ import (
    units as pyunits,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import *
from watertap.core.util.initialization import *
from watertap.examples.flowsheets.batch_RO.batch_ro import *

solver = get_solver()


def fix_dof_and_initialize(
    m,
    objective=None,
    outlvl=idaeslog.WARNING,
):
    """Fix degrees of freedom and initialize the flowsheet

    This function fixes the degrees of freedom of each unit and initializes the entire flowsheet.

    Args:
        m: Pyomo `Block` or `ConcreteModel` containing the flowsheet
        outlvl: Logger (default: idaeslog.WARNING)
    """
    set_operating_conditions(m, Qin=1, Cin=1)
    add_costing(m)
    initialize(m)

    return


if __name__ == "__main__":
    m = build_system()
    fix_dof_and_initialize(m)
    solve(m)
    print_results(m)