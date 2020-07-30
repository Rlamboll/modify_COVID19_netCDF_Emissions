import cartopy.crs as ccrs
import datetime
import matplotlib
import netCDF4 as nc
from matplotlib import pyplot as plt
import matplotlib.animation as anim


# Load data
include_aviation = False
baseline_file = "../output/aerosols/daily/cut_SO2-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.ncdaily_v4.nc_baseline.nc"
covid_file = "../output/aerosols/daily/cut_SO2-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.ncdaily_v4.nc_1_year.nc"
if include_aviation:
    aviation_base_file = "../output/aviation/"
varname = "SO2_em_anthro"
title_str = "Reduction in SO$_2$ emissions resulting from COVID-19 \n Date: {} {}."
var_label = "SO$_2$ emissions ({})"
units = "kg/km$^2$/day"
covid_data = nc.Dataset(covid_file)
baseline_data = nc.Dataset(baseline_file)
covid = covid_data.variables[varname][:]
base = baseline_data.variables[varname][:]
lats = covid_data.variables['lat'][:]
lons = covid_data.variables['lon'][:]
img_extent = [min(lons), max(lons), min(lats), max(lats)]
sect_ax = 1
basesum = (base.sum(axis=sect_ax) - covid.sum(axis=sect_ax)) * 10 ** 6 * 3600 * 24
base_time = baseline_data.variables["time"][:]
startind = 0
endind = 180
lowlim = 5.5
savename = "output/animated_COVID_rel_{}_v4.gif".format(varname)

# We want to plot some transform of the data to prevent nans from log(0)
max_data = basesum.max()

def data_transform(data, max_data):
    return data + max_data * 1e-15


def date_conv(days):
    date = datetime.datetime(year=2015, month=1, day=2) + datetime.timedelta(days=days)
    return date.strftime("%b"), date.day


# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure(figsize=(8, 4))
ax = plt.axes(projection=ccrs.PlateCarree())

cmap = matplotlib.cm.get_cmap("gist_earth")
cmap.set_bad('dimgrey', 1.)
plot_options = {
    "vmax": max_data,
    "vmin": max_data / 10**lowlim,
    "norm": matplotlib.colors.LogNorm(vmin=max_data / 10**lowlim, vmax=max_data),
    "cmap": cmap,
    "origin": "upper",
    "extent": img_extent,
    "transform": ccrs.PlateCarree(),
}
line = ax.imshow(data_transform(basesum[startind+160, ::-1, :], max_data), **plot_options)
ax.coastlines(color="lightgrey")
plt.title(title_str)
month, day = date_conv(base_time[startind])
plt.title(title_str.format(month, day))
cb = fig.colorbar(line, ax=ax, format="%g   ")
cb.ax.set_ylabel("Reduction in " + var_label.format(units))
cb.ax.invert_yaxis()
plt.text(-35, -105, "Dark grey implies a small emissions increase")

# initialization function: plot the background of each frame
def init():
    line = ax.imshow(data_transform(basesum[startind, ::-1, :], max_data),
                     **plot_options)
    month, day = date_conv(base_time[startind])
    plt.title(title_str.format(month, day))
    return line,

# animation function.  This is called sequentially
def animate(i):
    line = ax.imshow(data_transform(basesum[startind + i, ::-1, :], max_data),
                     **plot_options)
    month, day = date_conv(base_time[startind + i])
    plt.title(title_str.format(month, day))
    return line,


writer = anim.writers['pillow']
writer = writer(fps=4, metadata=dict(artist='Robin Lamboll'), bitrate=-1)


# call the animator.  blit=True means only re-draw the parts that have changed.
animation_inst = anim.FuncAnimation(
    fig, animate, init_func=init, frames=endind - startind,
    interval=20, blit=True
)
# Save the animation
animation_inst.save(savename, writer=writer)