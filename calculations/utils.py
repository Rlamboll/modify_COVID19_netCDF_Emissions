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
    """
    Copies a netcdf file from one folder to another, optionally appending a string, and
    returns the file open in memory for modification.
    :param filename:
        name of the file to read.
    :param input_folder: string
        Folder where the input file is.
    :param output_folder: string
        Folder where the output should go.
    :param scenario_string:
        String to append on the name of the scenario.
    :param compress:
        If true, the file is compressed (default:false).
    :param remove_string:
        String to remove from the name of the file.
    :return: :obj:`netcdf4.DataFrame`
        The copied file, open in memory.
    """
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
    """
    Performs an interpolation to add an additional timepoint to a netcdf file,
    optionally not using the timepoints immediately either side of the added point.
    :param db: :obj:`netcdf4.DataFrame`
        Structure into which points should be interpolated
    :param time_to_add: numerical
        The time at which a point should be added
    :param ind_before: int
        The number of already existing points before the new point to use for
        interpolating. E.g. to do month-to-month interpolation, use 12 here to reference
        a month of the same name a year ago in monthly data.
        Defaults to 1, i.e. point immediately before the new point.
    :param ind_after: The number of already existing points after the new point to use
        for interpolating. Defaults to 1, i.e. point immediately after the new point.
    :return: :obj:`netcdf4.DataFrame`
        The first input with an additional point interpolated into it
    """

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
    """
    This function cuts off data after a particular time and also compresses it if
    compress == True. The file has "cut_" prepended on the name and saves it to the
    output folder. The file is open in memory in the returned object.

    :param input_folder: string
        Folder where the input file is.
    :param output_folder: string
        Folder where the output should go.
    :param filename:
        name of the file to read.
    :param tcutoff: numerical
        End time - times after this will be cut off.
    :param scenario_string: string
        String to append on the end of the filename.
    :param compress: bool
        If true, (default) the file is compressed.
    :param remove_string: string
        If provided, will remove the given string from the name of the copied file.
    :param tstart: numerical
        Times before this will be cut off.
    :return: :obj:`netcdf4.DataFrame`
        The copied and cropped file, open in memory.
    """

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
    """
    This function removes files from the output folder if the filename ends in ".nc_",
    working_string or ".nc" + remove_scenario_string.
    :param output_folder: string
        Folder where files should be removed from.
    :param working_string: string
        Files ending in this string will be removed.
    :param remove_scenario_string: string
        Files ending in ".nc" followed by this string will be removed.
    :return:
    """
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
