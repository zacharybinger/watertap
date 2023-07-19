from watertap.core.util.model_diagnostics.infeasible import *
from idaes.core.util.scaling import *

def automate_rescale_variables(self, rescale_factor=1, default=1, verbose=True):
        if rescale_factor is None:
            rescale_factor = 1
        for var, sv in badly_scaled_var_generator(self):
            sf = get_scaling_factor(var)
            if get_scaling_factor(var) is None:
                if verbose:
                    print(f"{var} is missing a scaling factor")
                sf = default
                set_scaling_factor(var, sf, data_objects=False)

            set_scaling_factor(var, sf / sv * rescale_factor)
            calculate_scaling_factors(self)

def debug(m, solver=None, verbose=True, automate_rescale=False, resolve=False):
    badly_scaled_vars = list(badly_scaled_var_generator(m))
    if verbose:
        print(f'\n{"=======> DEBUGGING <=======":^60}\n')
        print(f'\n{"=======> BADLY SCALED VARIABLES <=======":^60}\n')
        print([print(i[0], i[1]) for i in badly_scaled_vars])
        print(f'\n{"=======> INFEASIBLE BOUNDS <=======":^60}\n')
        print_infeasible_bounds(m)
        print(f'\n{"=======> INFEASIBLE CONSTRAINTS <=======":^60}\n')
        print_infeasible_constraints(m)
        print(f'\n{"=======> CONSTRAINTS CLOSE TO BOUNDS <=======":^60}\n')
        print_close_to_bounds(m)
    if automate_rescale:
        if verbose:
            print(
                f"\n{len(badly_scaled_vars)} poorly scaled "
                f"variable(s) will be rescaled so that each scaled variable value = 1\n"
            )
        automate_rescale_variables(m, verbose=False)
        # badly_scaled_vars = list(badly_scaled_var_generator(m))
        if verbose:
            print(
                    f"\nNow {len(badly_scaled_vars)} poorly scaled\n"
                )

    if resolve:
        results = solver.solve(m, tee=True)
        return results
