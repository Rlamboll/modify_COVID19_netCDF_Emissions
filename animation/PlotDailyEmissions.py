import datetime
import math
import matplotlib.colors
import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as anim
import seaborn


# Load data
baseline_file = "../output/aerosols/cut_CO2-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-MESSAGE-GLOBIOM-ssp245-1-1_gn_201501-210012.ncv4.nc_baseline.nc"
baseline_data = nc.Dataset(baseline_file)
base = baseline_data.variables["CO2_em_anthro"][:]
sect_ax = 0
startind = 12
basesum = base.sum(axis=sect_ax)
base_time = baseline_data.variables["time"][:]
lowlim = 6

# We want to plot some transform of the data really
max_data = basesum.max()
def data_transform(data, max_data):
    return data + max_data * 1e-10

def date_conv(days):
    date = datetime.datetime(year=2020, month=1, day=1) + datetime.timedelta(days=days)
    return date.strftime("%b"), date.day

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
cbar_ticks = [
    math.pow(10, i) for i in range(
        math.floor(np.log10(max_data) - lowlim), 1+math.ceil(np.log10(max_data))
    )
]
line = seaborn.heatmap(
    data_transform(basesum[startind, ::-1, :], max_data),
    xticklabels=False, yticklabels=False,
    vmax=max_data,
    vmin=max_data / 10**lowlim,
    norm=matplotlib.colors.LogNorm(vmin=max_data / 10**lowlim, vmax=max_data),
    cbar_kws={"ticks": cbar_ticks}
)
m, d = date_conv(0)
plt.title("Date: {} {}".format(m, d))

# initialization function: plot the background of each frame
def init():
    line = seaborn.heatmap(
        data_transform(basesum[startind, ::-1, :], max_data),
        xticklabels=False, yticklabels=False,
        vmax=max_data,
        vmin=max_data / 10 ** lowlim,
        norm=matplotlib.colors.LogNorm(vmin=max_data / 10 ** lowlim, vmax=max_data),
        cbar=False
    )
    month, day = date_conv(0)
    plt.title("Date: {} {}".format(month, day))
    return line,

# animation function.  This is called sequentially
def animate(i):
    line = seaborn.heatmap(
        data_transform(basesum[startind + i, ::-1, :], max_data),
        xticklabels=False, yticklabels=False,
        vmax=max_data,
        vmin=max_data / 10 ** lowlim,
        norm=matplotlib.colors.LogNorm(vmin=max_data / 10 ** lowlim, vmax=max_data),
        cbar=False
    )
    month, day = date_conv(i)
    plt.title("Date: {} {}".format(month, day))
    return line,


writer = anim.writers['pillow']
writer = writer(fps=1, metadata=dict(artist='Robin Lamboll'), bitrate=1800)


# call the animator.  blit=True means only re-draw the parts that have changed.
animation_inst = anim.FuncAnimation(
    fig, animate, init_func=init, frames=5, # len(basesum) - startind,
    interval=20, blit=False
)
# Save the animation
animation_inst.save('basic_animation.gif', writer=writer)