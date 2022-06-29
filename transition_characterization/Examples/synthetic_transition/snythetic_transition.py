import numpy as np
import matplotlib.pyplot as plt
import joblib as jl
import sys
import statsmodels.api as sm
sys.path.append('../../')
from model import linear_ramp, fit_rmse
from distributions import sample_ar1
from transition_characterization import estimate_transition, combined_transition


#################################################################
# test the GLSAR from stats - to be removed from this file      #
#################################################################

# delta = 0.1
# time = np.arange(-100, 100, delta)
# xx = sample_ar1(len(time), 0.3, sigma=0.3, x0=0)
# dx = np.diff(xx)
# xx = sm.add_constant(xx)
# model = sm.GLSAR(dx, xx[:-1], rho=1)
# results = model.iterative_fit(maxiter=10)

# a = results.params[1]


#################################################################
# create synthetic transition                                   #
#################################################################

delta = 0.1
time = np.arange(-100, 100, delta)
trans = linear_ramp(time, t0=-10, dt=20)
noise = sample_ar1(len(time), 0.5, sigma=0.2, x0=0)
synt_trans = trans + noise

#################################################################
# least square fit of a linear ramp model to the noisy data     #
#################################################################

popt = fit_rmse(time, synt_trans)
cntrl = linear_ramp(time,
                    t0=popt[0], dt=np.exp(popt[1]),
                    y0=popt[2], dy=popt[3])

fig, ax = plt.subplots()
ax.plot(time, synt_trans, color='C0', label='synthetic transition')
ax.plot(time, cntrl, color='C1', label='rmse fit')
ax.axvline(popt[0], color='k', lw=0.5)
ax.axvline(popt[0] + np.exp(popt[1]), color='k', lw=0.5)
tax = ax.twiny()
tax.plot([], [])
tax.set_xlim(ax.get_xlim())
tax.set_xticks([popt[0], popt[0] + np.exp(popt[1])])
tax.set_xticklabels(['$t_0 = %0.2f$' % popt[0],
                     '$t_0 = %0.2f$' % (popt[1] + np.exp(popt[1]))],
                    rotation=45,
                    ha='left',
                    rotation_mode='anchor')

ax.set_xlabel('time')
ax.set_ylabel('$x$')
fig.subplots_adjust(top=0.8)
ax.legend()


#################################################################
# characterize the transition with the Bayesian linear ramp fit #
#################################################################

traces = estimate_transition(time,
                             synt_trans,
                             nwalkers=20, nsamples=600, nthin=1)

t0_dist = np.histogram(traces['t0'], time, density=True)
t1_dist = np.histogram(traces['t0'] + traces['dt'], time, density=True)

# make an axis for the distributions of y0 and dy
y_ax = np.linspace(np.min(synt_trans), np.max(synt_trans), 1000)
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


ax1.plot(time, synt_trans, color='C0', label='synthetic transition', lw=0.5)
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
                 color='C4', lw=1, label='$t_0$')
ax2.fill_between(hist_ax[mask1], t1_dist[0][mask1],
                 color='C6', lw=1, label='$t_1$')
ax2.set_xlabel('time')
ax2.set_ylabel('rel. probability')
ax2.yaxis.set_ticks_position('right')
ax2.yaxis.set_label_position('right')
ax2.set_xlim(ax1.get_xlim())
ax2.legend()
