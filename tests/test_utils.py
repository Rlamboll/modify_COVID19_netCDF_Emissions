import netCDF4 as nc
import numpy as np
import os

from ..calculations.utils import copy_netcdf_file

nsect = 5
latlen = 20
lonlen = 10
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
times[:] = np.arange(10)
sectors[:] = list(range(nsect))
lats[:] = np.arange(latlen)
lons[:] = np.arange(lonlen)
test_db.close()

def test_copy_netcdf():
    name_append = "_clone"
    folder = "./"
    new_scc = copy_netcdf_file(test_file, folder, folder, name_append)
    assert os.path.isfile(folder + test_file + name_append)
    new_scc.close()
    os.remove(folder + test_file + name_append)

# def test_insert_interpolated_point()
