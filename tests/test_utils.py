import netCDF4 as nc
import numpy as np
import os

from ..calculations.utils import copy_netcdf_file, \
    insert_interpolated_point, cutoff_netcdf_time

nsect = 5
latlen = 12
lonlen = 15
timeslen = 10
test_file = "test_file.nc"
test_db = nc.Dataset(test_file, mode="w", format="NETCDF4")
test_db.createDimension("time", None)
test_db.createDimension("sector", nsect)
test_db.createDimension("lat", latlen)
test_db.createDimension("lon", lonlen)
times = test_db.createVariable("time", float, ("time",))
sectors = test_db.createVariable("sector", int, ("sector",))
lats = test_db.createVariable("lat", float, ("lat",))
lons = test_db.createVariable("lon", float, ("lon",))
temp = test_db.createVariable("temp", float, ("time", "sector", "lat", "lon"))
times[:] = np.arange(timeslen)
sectors[:] = list(range(nsect))
lats[:] = np.arange(latlen)
lons[:] = np.arange(lonlen)
temp[:, :, :, :] = (
        times[:].reshape(timeslen, 1, 1, 1) * sectors[:].reshape(1, nsect, 1, 1) *
        lats[:].reshape(1, 1, latlen, 1) * lons[:].reshape(1, 1, 1, lonlen)
)
test_db.close()

name_append = "_clone"
folder = "./"

def test_copy_netcdf():
    new_scc = copy_netcdf_file(test_file, folder, folder, name_append)
    assert os.path.isfile(folder + test_file + name_append)
    new_scc.close()
    os.remove(folder + test_file + name_append)

def test_insert_interpolated_point():
    new_scc = copy_netcdf_file(test_file, folder, folder, name_append)
    startsize = len(new_scc.variables["temp"][:, :, :, :])
    startshape = new_scc.variables["temp"][:, :, :, :].shape
    newtime = 5.5
    expected_results = (new_scc.variables["temp"][5, :, :, :] +
            new_scc.variables["temp"][6, :, :, :]) / 2
    orig_vals = new_scc.variables["temp"][:, :, :, :]
    insert_interpolated_point(new_scc, newtime)
    assert len(new_scc.variables["temp"][:, :, :, :]) == startsize + 1
    assert all(
        new_scc.variables["temp"][
            :, :, :, :
        ].shape[i] == startshape[i] for i in range(1, 3)
    )
    assert new_scc.variables["time"][6] == newtime
    assert np.allclose(
        new_scc.variables["temp"][6, :, :, :], expected_results
    )
    assert np.allclose(new_scc.variables["temp"][-3:-1, :, :, :], orig_vals[-3:-1, :, :, :])
    secondtime = 6.5
    # This is between
    expected_results_2 = (expected_results * 2.5 + new_scc.variables["temp"][10, :, :, :] * 1) / 3.5
    insert_interpolated_point(new_scc, secondtime, 2, 3)
    assert len(new_scc.variables["temp"][:, :, :, :]) == startsize + 2
    assert np.allclose(new_scc.variables["temp"][8, :, :, :], expected_results_2)
    new_scc.close()
    os.remove(folder + test_file + name_append)

def test_cutoff_time():
    new_scc = copy_netcdf_file(test_file, folder, folder, name_append)
    orig_size = new_scc["temp"].shape
    new_scc.close()
    tcutoff = 6
    new_scc = cutoff_netcdf_time(folder, folder, test_file, tcutoff)
    assert new_scc["temp"].shape[0] == 7
    assert new_scc["temp"].shape[1:] == orig_size[1:]
    new_scc.close()
    os.remove(folder + test_file + "_cropped.nc")