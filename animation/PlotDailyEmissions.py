import cartopy.crs as ccrs
import datetime
import matplotlib
import netCDF4 as nc
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from matplotlib import pyplot as plt
import matplotlib.animation as anim


# This script plots animations of the NOx and SO2 data

calc_nox = True  # This switches between SO2 if false and nox if true.
# Load data
if not calc_nox:
    baseline_file = "../output/aerosols/daily/cut_SO2-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.ncdaily_v4.7.nc_baseline.nc"
    covid_file = "../output/aerosols/daily/cut_SO2-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.ncdaily_v4.7.nc_1_year.nc"
    varname = "SO2_em_anthro"
    title_str = "Fraction of usual SO$_2$ emissions due to COVID-19 \n Date: {} {}."
else:
    baseline_file = "../output/aerosols/daily/cut_NOx-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.ncdaily_v4.7.nc_baseline.nc"
    covid_file = "../output/aerosols/daily/cut_NOx-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.ncdaily_v4.7.nc_1_year.nc"
    av_varname = "NOx_em_AIR_anthro"
    aviation_base_file = "../output/aviation/cut_NOx-em-AIR-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.nc_baseline_v4.7.nc"
    aviation_covid_file = "../output/aviation/cut_NOx-em-AIR-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.nc_v4.7.nc_flightrd_mp06_daily_v4.7.nc"
    # Which dimension index is the height (to be summed over)?
    height_ax = 1
    # Read the aviation data
    base_av_ds = nc.Dataset(aviation_base_file)
    covid_av_ds = nc.Dataset(aviation_covid_file)
    covid_av = covid_av_ds.variables[av_varname][:]
    base_av = base_av_ds.variables[av_varname][:]
    assert np.allclose(
        covid_av_ds.variables["time"][:], base_av_ds.variables["time"][:]
    )
    varname = "NOx_em_anthro"
    title_str = "Fraction of usual NO$_x$ emissions due to COVID-19 \n Date: {} {}."

savename = "output/animated_COVID_rel_lin_{}_v7_light.gif".format(varname)

# Load
covid_data = nc.Dataset(covid_file)
baseline_data = nc.Dataset(baseline_file)
covid = covid_data.variables[varname][:]
base = baseline_data.variables[varname][:]
lats = covid_data.variables["lat"][:]
lons = covid_data.variables["lon"][:]
img_extent = [min(lons), max(lons), min(lats), max(lats)]
sect_ax = 1


def data_transform(data, _):
    return data.filled(1)


def date_conv(days):
    date = datetime.datetime(year=2015, month=1, day=2) + datetime.timedelta(days=days)
    return date.strftime("%b"), date.day


def interpolate_missing_times(ds, have_times, want_times, lats, lons):
    interpolating_fn = RegularGridInterpolator((have_times, lats, lons), ds)
    return interpolating_fn(
        [
            (x, y, z)
            for x in want_times.filled(np.nan)
            for y in lats.filled(np.nan)
            for z in lons.filled(np.nan)
        ]
    )


base_time = baseline_data.variables["time"][:]
startind = 0
endind = 210
if not calc_nox:
    basesum = covid.sum(axis=sect_ax) / base.sum(axis=sect_ax)
else:
    # Need to include aviation emissions
    base_av_int = interpolate_missing_times(
        base_av.sum(axis=height_ax),
        base_av_ds.variables["time"][:],
        baseline_data.variables["time"][startind:endind],
        lats,
        lons,
    ).reshape(endind - startind, len(lats), len(lons))
    covid_av_int = interpolate_missing_times(
        covid_av.sum(axis=height_ax),
        covid_av_ds.variables["time"][:],
        baseline_data.variables["time"][startind:endind],
        lats,
        lons,
    ).reshape(endind - startind, len(lats), len(lons))
    basesum = (covid.sum(axis=sect_ax)[startind:endind, ...] + covid_av_int) / (
        base.sum(axis=sect_ax)[startind:endind, ...] + base_av_int
    )
# We want to ensure the max and min values always fit on our scale
max_data = basesum.max()
min_data = basesum.min()

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure(figsize=(8, 4))
ax = plt.axes(projection=ccrs.PlateCarree())

cmap = matplotlib.cm.get_cmap("inferno")  # "gist_earth"
cmap.set_bad("dimgrey", 1.0)
plot_options = {
    "vmax": max_data,
    "vmin": min_data,
    "cmap": cmap,
    "origin": "upper",
    "extent": img_extent,
    "transform": ccrs.PlateCarree(),
    "aspect": 'auto'
}
line = ax.imshow(
    data_transform(basesum[startind + 160, ::-1, :], max_data), **plot_options
)
ax.coastlines(color="white")
plt.title(title_str)
month, day = date_conv(base_time[startind])
plt.title(title_str.format(month, day))
cb = fig.colorbar(line, ax=ax, format="%g")
cb.ax.set_ylabel("Emissions proportion")

# initialization function: plot the background of each frame
def init():
    line = ax.imshow(
        data_transform(basesum[startind, ::-1, :], max_data), **plot_options
    )
    month, day = date_conv(base_time[startind])
    plt.title(title_str.format(month, day))
    return (line,)


# animation function.  This is called sequentially
def animate(i):
    line = ax.imshow(
        data_transform(basesum[startind + i, ::-1, :], max_data), **plot_options
    )
    month, day = date_conv(base_time[startind + i])
    plt.title(title_str.format(month, day))
    return (line,)


writer = anim.writers["pillow"]
writer = writer(fps=4, metadata=dict(artist="Robin Lamboll"), bitrate=-1)


# Call the animator.  blit=True means only re-draw the parts that have changed.
animation_inst = anim.FuncAnimation(
    fig, animate, init_func=init, frames=endind - startind, interval=20, blit=True
)
# Save the animation
animation_inst.save(savename, writer=writer)
