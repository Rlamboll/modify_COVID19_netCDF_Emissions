import datetime
import math
import matplotlib.colors
import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as anim
import seaborn


# Load data
baseline_file = "../output/aerosols/daily/cut_SO2-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.ncdaily_v4.nc_baseline.nc"
covid_file = "../output/aerosols/daily/cut_SO2-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.ncdaily_v4.nc_1_year.nc"
varname = "SO2_em_anthro"
var_label = "SO$_2$ emissions ({})"
units = "kg/km$^2$/s"
covid_data = nc.Dataset(covid_file)
baseline_data = nc.Dataset(baseline_file)
covid = covid_data.variables[varname][:]
base = baseline_data.variables[varname][:]

sect_ax = 1
basesum = (base.sum(axis=sect_ax) - covid.sum(axis=sect_ax)) * 10 ** 6
base_time = baseline_data.variables["time"][:]
startind = 0
endind = 180
lowlim = 6
savename = "animated_COVID_{}_v2.gif".format(varname)

# We want to plot some transform of the data to prevent nans from log(0)
max_data = basesum.max()
def data_transform(data, max_data):
    return data + max_data * 1e-15

def date_conv(days):
    date = datetime.datetime(year=2015, month=1, day=2) + datetime.timedelta(days=days)
    return date.strftime("%b"), date.day

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure(figsize=(8, 4))
cbar_ticks = [
    math.pow(10, i) for i in range(
        math.floor(np.log10(max_data) - lowlim), math.ceil(np.log10(max_data))
    )
]

plot_options = {
    "xticklabels": False, "yticklabels": False,
    "vmax": max_data,
    "vmin": max_data / 10**lowlim,
    "norm": matplotlib.colors.LogNorm(vmin=max_data / 10**lowlim, vmax=max_data)
}
line = seaborn.heatmap(
    data_transform(basesum[startind, ::-1, :], max_data),
    cbar_kws={"ticks": cbar_ticks, 'label': var_label.format(units)},
    **plot_options
)
month, day = date_conv(base_time[startind])
plt.title("Date: {} {}".format(month, day))

# initialization function: plot the background of each frame
def init():
    line = seaborn.heatmap(
        data_transform(basesum[startind, ::-1, :], max_data),
        cbar=False,
        **plot_options
    )
    month, day = date_conv(base_time[startind])
    plt.title("Date: {} {}".format(month, day))
    return line,

# animation function.  This is called sequentially
def animate(i):
    line = seaborn.heatmap(
        data_transform(basesum[startind + i, ::-1, :], max_data),
        cbar=False,
        **plot_options
    )
    month, day = date_conv(base_time[startind + i])
    plt.title("Date: {} {}".format(month, day))
    return line,


writer = anim.writers['pillow']
writer = writer(fps=3, metadata=dict(artist='Robin Lamboll'), bitrate=-1)


# call the animator.  blit=True means only re-draw the parts that have changed.
animation_inst = anim.FuncAnimation(
    fig, animate, init_func=init, frames=endind - startind,
    interval=20, blit=False
)
# Save the animation
animation_inst.save(savename, writer=writer)