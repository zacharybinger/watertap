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

import os, pathlib, warnings
import h5py
import itertools
import pprint
import copy
import numpy as np

from scipy.interpolate import griddata

from pyomo.common.config import ConfigDict, ConfigValue


class ParameterSweepWriter:

    CONFIG = ConfigDict()

    CONFIG.declare(
        "debugging_data_dir",
        ConfigValue(
            default=None,
            domain=str,
            description="directory path to output debugging data.",
        ),
    )

    CONFIG.declare(
        "csv_results_file_name",
        ConfigValue(
            default=None, domain=str, description="filepath to the output CSV file."
        ),
    )

    CONFIG.declare(
        "h5_results_file_name",
        ConfigValue(
            default=None, domain=str, description="filepath to the output H5 file."
        ),
    )

    CONFIG.declare(
        "interpolate_nan_outputs",
        ConfigValue(
            default=False,
            domain=bool,
            description="Bool to decide whether to interpolate NaN outputs.",
        ),
    )

    CONFIG.declare(
        "h5_parent_group_name",
        ConfigValue(
            default=None,
            domain=str,
            description="Parent group name (container like objects, similar to a folder/directory in a file system) for parameter sweep outputs to be saved.",
        ),
    )

    def __init__(
        self,
        parallel_manager,
        **options,
    ):

        self.parallel_manager = parallel_manager
        self.config = self.CONFIG(options)

        if self.parallel_manager.is_root_process():
            if (
                self.config.h5_results_file_name is None
                and self.config.csv_results_file_name is None
            ):
                warnings.warn(
                    "No results will be writen to disk as h5_results_file_name and csv_results_file_name are both None"
                )

    @staticmethod
    def _strip_extension(file_name, extension):
        if file_name.lower().endswith(extension):
            return file_name[: -len(extension)], extension
        else:
            return file_name, None

    def _process_results_filename(self, results_file_name):
        # Get the directory path
        dirname = os.path.dirname(results_file_name)
        # Get the file name without the extension
        known_extensions = [".h5", ".csv"]
        for ext in known_extensions:
            fname_no_ext, extension = self._strip_extension(results_file_name, ext)
            if extension is not None:
                break

        return dirname, fname_no_ext, extension

    @staticmethod
    def _interp_nan_values(global_values, global_results):

        global_results_clean = np.copy(global_results)

        n_vals = np.shape(global_values)[1]
        n_outs = np.shape(global_results)[1]

        # Build a mask of all the non-nan saved outputs
        # i.e., where the optimzation succeeded
        mask = np.isfinite(global_results[:, 0])

        # Create a list of points where good data is available
        x0 = global_values[mask, :]

        if np.sum(mask) >= 4:
            # Interpolate to get a value for nan points where possible
            for k in range(n_outs):
                y0 = global_results[mask, k]
                yi = griddata(
                    x0, y0, global_values, method="linear", rescale=True
                ).reshape(-1)
                global_results_clean[~mask, k] = yi[~mask]

        else:
            warnings.warn("Too few points to perform interpolation.")

        return global_results_clean

    def _write_debug_data(
        self,
        sweep_params,
        local_values,
        local_results_dict,
        write_h5,
        write_csv,
        process_number,
    ):

        if write_h5:
            fname_h5 = f"local_results_{process_number:03}.h5"
            self._write_output_to_h5(
                local_results_dict,
                os.path.join(self.config["debugging_data_dir"], fname_h5),
            )
        if write_csv:
            fname_csv = f"local_results_{process_number:03}.csv"

            data_header = ",".join(itertools.chain(sweep_params))
            local_results = np.zeros(
                (np.shape(local_values)[0], len(local_results_dict["outputs"])),
                dtype=float,
            )
            for i, (key, item) in enumerate(local_results_dict["outputs"].items()):
                data_header = ",".join([data_header, key])
                local_results[:, i] = item["value"][:]

            local_save_data = np.hstack((local_values, local_results))

            # Save the local data
            np.savetxt(
                os.path.join(self.config["debugging_data_dir"], fname_csv),
                local_save_data,
                header=data_header,
                delimiter=", ",
                fmt="%.6e",
            )

    def _write_outputs(self, output_dict, txt_options="metadata"):

        self._write_output_to_h5(output_dict, self.config["h5_results_file_name"])

        # We will also create a companion txt file by default which contains
        # the metadata of the h5 file in a user readable format.
        txt_fname = self.config["h5_results_file_name"] + ".txt"

        if (
            self.config.h5_parent_group_name is None
            and os.path.isfile(txt_fname) is False
        ):
            if txt_options == "metadata":
                my_dict = copy.deepcopy(output_dict)
                for key, value in my_dict.items():
                    for subkey, subvalue in value.items():
                        subvalue.pop("value")
            elif txt_options == "keys":
                my_dict = {}
                for key, value in output_dict.items():
                    if key in ["sweep_params", "outputs"]:
                        my_dict[key] = list(value.keys())
            else:
                my_dict = output_dict

            with open(txt_fname, "w") as log_file:
                pprint.pprint(my_dict, log_file)

    def _write_output_to_h5(self, output_dict, h5_results_file_name):

        if self.config.h5_parent_group_name is None:
            # No parent groups exists, a new file will be created regardless
            f = h5py.File(h5_results_file_name, "w")
            parent_grp = f
        else:
            if os.path.isfile(h5_results_file_name):
                # File exists, we only need to add the new parent group
                f = h5py.File(h5_results_file_name, "a")
                parent_grp = f.require_group(self.config.h5_parent_group_name)
            else:
                # Create a new file since none exists and add the parent group
                f = h5py.File(h5_results_file_name, "w")
                parent_grp = f.require_group(self.config.h5_parent_group_name)

        for key, item in output_dict.items():
            grp = parent_grp.create_group(key)
            if key in ["sweep_params", "outputs"]:
                for subkey, subitem in item.items():
                    subgrp = grp.create_group(subkey)
                    for subsubkey, subsubitem in subitem.items():
                        if subsubkey[0] != "_":
                            if subsubkey == "lower bound" and subsubitem is None:
                                subgrp.create_dataset(subsubkey, data=np.finfo("d").min)
                            elif subsubkey == "upper bound" and subsubitem is None:
                                subgrp.create_dataset(subsubkey, data=np.finfo("d").max)
                            else:
                                subgrp.create_dataset(subsubkey, data=subsubitem)
            elif key in ["solve_successful", "nominal_idx", "differential_idx"]:
                grp.create_dataset(key, data=output_dict[key])

        f.close()

    def _write_to_csv(
        self,
        sweep_params,
        global_values,
        global_results_dict,
        global_results_arr,
        process_number,
    ):

        # Create the dataframe that is going to be written to a CSV
        global_save_data = np.hstack((global_values, global_results_arr))

        if process_number == self.parallel_manager.ROOT_PROCESS_RANK:
            data_header = ",".join(itertools.chain(sweep_params))
            for i, (key, item) in enumerate(global_results_dict["outputs"].items()):
                data_header = ",".join([data_header, key])

            if self.config["csv_results_file_name"] is not None:
                # Write the CSV
                np.savetxt(
                    self.config["csv_results_file_name"],
                    global_save_data,
                    header=data_header,
                    delimiter=",",
                    fmt="%.6e",
                )

                # If we want the interpolated output_list in CSV
                if self.config["interpolate_nan_outputs"]:
                    global_results_clean = self._interp_nan_values(
                        global_values, global_results_arr
                    )
                    global_save_data_clean = np.hstack(
                        (global_values, global_results_clean)
                    )

                    head, tail = os.path.split(self.config["csv_results_file_name"])

                    if head == "":
                        interp_file = "interpolated_%s" % (tail)
                    else:
                        interp_file = "%s/interpolated_%s" % (head, tail)

                    np.savetxt(
                        interp_file,
                        global_save_data_clean,
                        header=data_header,
                        delimiter=",",
                        fmt="%.6e",
                    )

        return global_save_data

    def save_results(
        self,
        sweep_params,
        local_values,
        global_values,
        local_results_dict,
        global_results_dict,
        global_results_arr,
        process_number,
    ):

        if process_number == self.parallel_manager.ROOT_PROCESS_RANK:
            if self.config["debugging_data_dir"] is not None:
                os.makedirs(self.config["debugging_data_dir"], exist_ok=True)
            if self.config["h5_results_file_name"] is not None:
                pathlib.Path(self.config["h5_results_file_name"]).parent.mkdir(
                    parents=True, exist_ok=True
                )
            if self.config["csv_results_file_name"] is not None:
                pathlib.Path(self.config["csv_results_file_name"]).parent.mkdir(
                    parents=True, exist_ok=True
                )

        self.parallel_manager.sync_with_peers()

        # Handle values in the debugging data_directory
        if self.config["debugging_data_dir"] is not None:
            self._write_debug_data(
                sweep_params,
                local_values,
                local_results_dict,
                self.config["h5_results_file_name"] is not None,
                self.config["csv_results_file_name"] is not None,
                process_number,
            )

        global_save_data = self._write_to_csv(
            sweep_params,
            global_values,
            global_results_dict,
            global_results_arr,
            process_number,
        )

        if (
            process_number == self.parallel_manager.ROOT_PROCESS_RANK
            and self.config["h5_results_file_name"] is not None
        ):
            # Save the data of output dictionary
            self._write_outputs(global_results_dict, txt_options="keys")

        return global_save_data
