import netCDF4 as nc
import numpy as np
import os


def copy_netcdf_file(
    filename,
    input_folder,
    output_folder,
    scenario_string,
    compress=False,
    remove_string=None,
):
    src = nc.Dataset(input_folder + filename)
    if remove_string:
        trg = nc.Dataset(
            output_folder + filename.replace(remove_string, "") + scenario_string,
            mode="w",
        )
    else:
        trg = nc.Dataset(output_folder + filename + scenario_string, mode="w")

    # Create the dimensions of the file
    for name, dim in src.dimensions.items():
        trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # Copy the global attributes
    trg.setncatts({a: src.getncattr(a) for a in src.ncattrs()})

    # Create the variables in the file
    for name, var in src.variables.items():
        if compress:
            trg.createVariable(
                name, var.dtype, var.dimensions, complevel=9, shuffle=True, zlib=True
            )
        else:
            trg.createVariable(name, var.dtype, var.dimensions)
        # Copy the variable attributes
        trg.variables[name].setncatts({a: var.getncattr(a) for a in var.ncattrs()})

        # Copy the variables values (as 'f4' eventually)
        trg.variables[name][:] = src.variables[name][:]

    # Return the data
    src.close()
    return trg


def insert_interpolated_point(db, time_to_add, ind_before=1, ind_after=1):
    times = db.variables["time"][:]
    if any(times == time_to_add):
        raise ValueError("time already in the database")
    after_time_ind = np.where(times > time_to_add)[0].min()
    assert after_time_ind > 0, "Cannot add element before the start by interpolation"
    step_before = (time_to_add - db.variables["time"][after_time_ind - ind_before]) / (
        db.variables["time"][after_time_ind + ind_after - 1]
        - db.variables["time"][after_time_ind - ind_before]
    )
    assert step_before < 1
    db.variables["time"][after_time_ind : len(times) + 1] = db.variables["time"][
        after_time_ind - 1 : len(times)
    ]
    db.variables["time"][after_time_ind] = time_to_add
    for var in db.variables:
        if db.variables[var].dimensions[0] == "time" and var != "time":
            db.variables[var][-1, ...] = db.variables[var][-2, ...]
            db.variables[var][after_time_ind:-1, ...] = db.variables[var][
                after_time_ind - 1 : -2, ...
            ]
            # We weight the values before and after in proportion to closeness,
            # which is 1 - the fraction of the step between the points.
            db.variables[var][after_time_ind, ...] = (
                db.variables[var][after_time_ind - ind_before, ...] * (1 - step_before)
                + db.variables[var][after_time_ind + ind_after, ...] * step_before
            )
        elif (
            var != "time"
            and "time" in db.variables[var].dimensions
            and db.variables[var].dimensions[1] == "time"
        ):
            db.variables[var][:, -1, ...] = db.variables[var][:, -2, ...]
            db.variables[var][:, after_time_ind:-1, ...] = db.variables[var][
                :, after_time_ind - 1 : -2, ...
            ]
            # We weight the values before and after in proportion to closeness,
            # which is 1 - the fraction of the step between the points.
            db.variables[var][:, after_time_ind, ...] = (
                db.variables[var][:, after_time_ind - ind_before, ...]
                * (1 - step_before)
                + db.variables[var][:, after_time_ind + ind_after, ...]
                * step_before
            )


def cutoff_netcdf_time(
    input_folder,
    output_folder,
    filename,
    tcutoff,
    scenario_string="_cropped.nc",
    compress=True,
    remove_string=None,
    tstart=None,
):
    # This function cuts off data after a particular time and also compresses it if
    # compress == True.
    db = nc.Dataset(input_folder + filename)
    if remove_string:
        trg = nc.Dataset(
            output_folder
            + "cut_"
            + filename.replace(remove_string, "")
            + scenario_string,
            mode="w",
        )
    else:
        trg = nc.Dataset(output_folder + "cut_" + filename + scenario_string, mode="w")
    times = db.variables["time"][:]
    assert tcutoff > min(times)
    if tcutoff > max(times):
        return
    if tstart:
        valid_times = np.where((times <= tcutoff) & (times >= tstart))[0]
    else:
        valid_times = np.where(times <= tcutoff)[0]

    # Create the dimensions of the file
    for name, dim in db.dimensions.items():
        trg.createDimension(
            name, len(dim) if not (dim.isunlimited()) and not (name == "time") else None
        )

    # Copy the global attributes
    trg.setncatts({a: db.getncattr(a) for a in db.ncattrs()})

    # Create the variables in the file
    for name, var in db.variables.items():
        if compress:
            trg.createVariable(
                name, var.dtype, var.dimensions, complevel=9, shuffle=True, zlib=True
            )
        else:
            trg.createVariable(name, var.dtype, var.dimensions)

        # Copy the variable attributes
        trg.variables[name].setncatts({a: var.getncattr(a) for a in var.ncattrs()})

        # Copy the variables values, removing some times if needed
        if "time" in var.dimensions[:]:
            if len(var.dimensions) == 1:  # We assume time is the 2nd dimension if mult
                trg.variables[name][:] = db.variables[name][valid_times]
            elif var.dimensions[0] == "time":
                trg.variables[name][:] = db.variables[name][valid_times, ...]
            else:
                trg.variables[name][:] = db.variables[name][:, valid_times, ...]
        else:
            trg.variables[name][:] = db.variables[name][:]

    # Return the data
    db.close()
    return trg


def cleanup_files(output_folder, working_string, remove_scenario_string=None):
    # This function removes files from the output folder if the filename ends in ".nc_",
    # working_string or ".nc" + remove_scenario_string.
    output_files = os.listdir(output_folder)
    deletable_files = [
        file
        for file in output_files
        if (file[-4:] == ".nc_")
        or (file[-len(working_string) :] == working_string)
        or (
            remove_scenario_string
            and (
                file[-len(remove_scenario_string) - 3 :]
                == ".nc" + remove_scenario_string
            )
        )
    ]
    print("Deleting files {}".format(deletable_files))
    for file in deletable_files:
        os.remove(output_folder + file)
