import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from preprocessing import download, binning
sys.path.append('../../')
from transition_characterization import estimate_transition
from transition_characterization import combined_transition


#################################################################
# data download and import                                      #
#################################################################

data = download()

#################################################################
# omitting nan and binning to equally spaced time axis whereby  #
# dt determines the time step                                   #
#################################################################

d18o_idx = data['d18o'].dropna().index.values
d18o = data['d18o'].iloc[d18o_idx].values
d18o_age = data['age'].iloc[d18o_idx].values

dt = 5

d18o_bins = np.arange(int(d18o_age[0]-dt/2),
                      int(d18o_age[-1] + dt),
                      dt)

d18o_age, d18o = binning(d18o_age, d18o, d18o_bins)


#################################################################
# selecting data window around onset of GI-8c                   # 
#################################################################

d18o_time = -d18o_age[::-1]
d18o = d18o[::-1]

ti = -39000
tf = -37300

mask = (d18o_time > ti) & (d18o_time < tf)

trans = d18o[mask]
time = d18o_time[mask]

#################################################################
# running transition onset characterization                     #
#################################################################

traces = estimate_transition(time,
                             trans,
                             nwalkers=20, nsamples=600, nthin=1)

t0_dist = np.histogram(traces['t0'], time, density=True)
t1_dist = np.histogram(traces['t0'] + traces['dt'], time, density=True)

# make an axis for the distributions of y0 and dy
y_ax = np.linspace(np.min(trans), np.max(trans), 1000)
y0_dist = np.histogram(traces['y0'], y_ax, density=True)
dy_dist = np.histogram(traces['dy'], y_ax, density=True)


p5, p50, p95 = combined_transition(time, traces)

fig = plt.figure(figsize=(8, 5))
height_ratios = [3, 1, 1]
gs = fig.add_gridspec(3, 1,
                      hspace=0.,
                      left=0.1,
                      bottom=0.15,
                      top=0.85,
                      right=0.9,
                      height_ratios=height_ratios)

ax1 = fig.add_subplot(gs[:2])
ax1.spines['bottom'].set_visible(False)
ax1.xaxis.set_visible(False)
ax1.patch.set_visible(False)
ax2 = fig.add_subplot(gs[1:])
ax2.spines['top'].set_visible(False)
ax2.patch.set_visible(False)


ax1.plot(time, trans, color='C0', label='observations', lw=0.5)
ax1.plot(time, p50, color='k', lw=1.2, label='50th percentile', zorder=12)
ax1.plot(time, p5, color='slategray', lw=0.8)
ax1.plot(time, p95, color='slategray',  lw=0.8)
ax1.fill_between(time, p5, p95, color='C1', alpha=.8,
                 zorder=10, label='90th percentile range')
ax1.set_ylabel('$x$')
ax1.legend()

hist_ax = t0_dist[1][:-1] + (t0_dist[1][1] - t0_dist[1][0])/2
mask0 = t0_dist[0] > 0.001
mask1 = t1_dist[0] > 0.001
ax2.fill_between(hist_ax[mask0], t0_dist[0][mask0],
                 color='C4', lw=1, label='$t_0$', alpha = 0.5)
ax2.fill_between(hist_ax[mask1], t1_dist[0][mask1],
                 color='C6', lw=1, label='$t_1$', alpha = 0.5)
ax2.set_xlabel('time')
ax2.set_ylabel('rel. probability')
ax2.yaxis.set_ticks_position('right')
ax2.yaxis.set_label_position('right')
ax2.set_xlim(ax1.get_xlim())
ax2.legend()

