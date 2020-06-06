import netCDF4 as nc
import numpy as np

def copy_netcdf_file(filename, input_folder, output_folder, scenario_string):
    src = nc.Dataset(input_folder + filename)
    trg = nc.Dataset(output_folder + filename + scenario_string, mode='w')

    # Create the dimensions of the file
    for name, dim in src.dimensions.items():
        trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # Copy the global attributes
    trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})

    # Create the variables in the file
    for name, var in src.variables.items():
        trg.createVariable(name, var.dtype, var.dimensions, complevel=9)

        # Copy the variable attributes
        trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})

        # Copy the variables values (as 'f4' eventually)
        trg.variables[name][:] = src.variables[name][:]

    # Return the data
    src.close()
    return trg

def insert_interpolated_point(db, time_to_add):
    times = db.variables["time"][:]
    if any(times == time_to_add):
        raise ValueError("time already in the database")
    after_time_ind = np.where(times > time_to_add)[0].min()
    assert after_time_ind > 0, "Cannot add element before the start by interpolation"
    step_before = (time_to_add - db.variables["time"][after_time_ind - 1]) / (
        db.variables["time"][after_time_ind] - db.variables["time"][after_time_ind - 1]
    )
    assert step_before < 1
    db.variables["time"][after_time_ind:len(times) + 1] = db.variables["time"][after_time_ind - 1:len(times)]
    db.variables["time"][after_time_ind] = time_to_add
    for var in db.variables:
        if "time" in db.variables[var].dimensions and var != "time":
            db.variables[var][
                after_time_ind:len(times), ...
            ] = db.variables[var][after_time_ind - 1:len(times) - 1, ...]
            db.variables[var][after_time_ind, ...] = (
                db.variables[var][after_time_ind - 1, ...] * step_before
                + db.variables[var][after_time_ind + 1, ...] * (1 - step_before)
            )