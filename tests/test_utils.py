import netCDF4 as nc
import numpy as np
import os
import pytest

from ..calculations.utils import copy_netcdf_file, \
    insert_interpolated_point, cutoff_netcdf_time

sectlen = 5
latlen = 12
lonlen = 15
timeslen = 10
folder = os.path.join(os.path.abspath(__file__), "..", "testData/")
test_file = "test_file.nc"
test_db = nc.Dataset(folder + test_file, mode="w", format="NETCDF4")
test_db.createDimension("time", None)
test_db.createDimension("sector", sectlen)
test_db.createDimension("lat", latlen)
test_db.createDimension("lon", lonlen)
times = test_db.createVariable("time", float, ("time",))
sectors = test_db.createVariable("sector", int, ("sector",))
lats = test_db.createVariable("lat", float, ("lat",))
lons = test_db.createVariable("lon", float, ("lon",))
temp = test_db.createVariable("temp", float, ("time", "sector", "lat", "lon"))
time_bands = test_db.createVariable("time_bands", float, ("time", "sector"))
times[:] = np.arange(timeslen)
sectors[:] = list(range(sectlen))
lats[:] = np.arange(latlen)
lons[:] = np.arange(lonlen)
temp[:, :, :, :] = (
        times[:].reshape(timeslen, 1, 1, 1) * sectors[:].reshape(1, sectlen, 1, 1) *
        lats[:].reshape(1, 1, latlen, 1) * lons[:].reshape(1, 1, 1, lonlen)
)
time_bands[:, :] = times[:].reshape(timeslen, 1) * sectors[:].reshape(1, sectlen)
test_db.close()

name_append = "_clone"


@pytest.mark.parametrize("compress", [True, False])
def test_copy_netcdf(compress):
    new_scc = copy_netcdf_file(
        test_file, folder, folder, name_append, compress=compress
    )
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


@pytest.mark.parametrize("compress", [True, False])
def test_cutoff_time(compress):
    new_scc = copy_netcdf_file(test_file, folder, folder, name_append)
    orig_size = new_scc["temp"].shape
    new_scc.close()
    tcutoff = 6
    new_scc = cutoff_netcdf_time(folder, folder, test_file, tcutoff, compress=compress)
    assert new_scc["temp"].shape[0] == 7
    assert new_scc["temp"].shape[1:] == orig_size[1:]
    assert new_scc["temp"][-1, 0, 0, 0] == 0
    assert new_scc["time_bands"][-1, 0] == 0
    new_scc.close()
    os.remove(folder + "cut_" + test_file + "_cropped.nc")
    new_scc = cutoff_netcdf_time(
        folder, folder, test_file, tcutoff, compress=compress, remove_string="il"
    )
    new_scc.close()
    assert os.path.isfile(folder + "cut_" + "test_fe.nc_cropped.nc")
    os.remove(folder + "cut_" + "test_fe.nc_cropped.nc")
    os.remove(folder + "test_file.nc_clone")

