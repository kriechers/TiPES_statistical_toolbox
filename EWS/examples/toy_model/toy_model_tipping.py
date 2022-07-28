import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sys
sys.path.append('../../')
from EWS_functions import EWS, bayesian_AR1

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

#################################################################
# defining a dynamical system which can switch between mono and #
# bistability depending on the value of a control parameter     #
# alpha                                                         #
#################################################################

# noise level
sigma = 0.2

# time step 
dt = 0.01

# default value of the control
alpha = 0


def xdot(x,sigma=sigma,dt=dt,alpha=alpha):
    drift = (-x**3 + x + alpha)*dt
    diffusion = np.random.normal()*sigma * np.sqrt(dt)
    return drift + diffusion


#################################################################
# compute the bifurcation diagram                               #
#################################################################

alpha_ax = np.arange(-2,2, 0.01)
x_ax = np.arange(-3,3,0.01)
xfp_list = []
for a in alpha_ax:

    mask = np.diff(xdot(x_ax, sigma = 0, alpha = a) < 0)
    starters = x_ax[:-1][mask]

    def fun(x): return xdot(x, sigma = 0, alpha = a)
    
    for s in starters:
        solution = fsolve(fun, s, maxfev=1000)
        xfp_list.append((a, solution[0]))

xfp_list = np.array(xfp_list)


### the value of the local minimum and maximum of xdot at     ###
### alpha = 0 namely x_min = -sqrt(1/3) and x_max = sqrt(1/3) ###
### provide a upper and lower bound for the lower and upper   ###
### stable branch of the nullcline                            ###

upper = xfp_list[xfp_list[:,1] > np.sqrt(1/3)]
lower = xfp_list[xfp_list[:,1] < -np.sqrt(1/3)]
unstable = xfp_list[(xfp_list[:,1] < np.sqrt(1/3)) &
                    (xfp_list[:,1] > -np.sqrt(1/3))]


fig, ax = plt.subplots()
ax.plot(*upper.T, color = 'C0', label = 'stable')
ax.plot(*lower.T, color = 'C0')
ax.plot(*unstable.T, color = 'C1', label = 'unstable', ls = '--')
ax.legend()


### compute the critical values for alpha                     ###

val, counts = np.unique(xfp_list[:,0], return_counts = True)
alpha_crit1 = np.min(val[counts == 3])
alpha_crit2 = np.max(val[counts == 3])


#################################################################
# create time series with bifurcation induced tipping           #
#################################################################

t_f = 600
t_1 = t_f /4
t_2 = 3 * t_f/4
time = np.arange(0, t_f, dt)
cond_list = [time <= t_1,
             (t_1 < time) & (time < t_2),
             t_2 <= time]

alpha0 = -1.5
alpha1 = 1.5
linear_ramp = lambda t: alpha0 + (alpha1-alpha0) * (t-t_1)/(t_2-t_1)

alpha_ramp = np.piecewise(time,
                          cond_list,
                          [alpha0,linear_ramp, alpha1])

# compute t_crit at which alpha_crit2 is exceeded

t_crit = time[np.sum(alpha_ramp < alpha_crit2)]


# integrate the system with alpha ramp

xx = np.zeros_like(time)
xx[0] = -1

for i,t in enumerate(time[:-1]):
    xx[i+1] = xx[i] + xdot(xx[i], alpha = alpha_ramp[i])

# chose running window size in time
rw = 10

# apply EWS analysis to data before alpha crosses alpha_crit1
mask = time < t_crit - rw

rv, rac, rrf, ct, pval = EWS(time[mask], xx[mask], int(rw / dt), 1)


#################################################################
# creating the plot                                             #
#################################################################

fig = plt.figure(figsize=(16,9))
gs = fig.add_gridspec(14, 1,
                      hspace=0.5,
                      left=0.1,
                      bottom=0.15,
                      top=0.85,
                      right=0.9)


ax1 = fig.add_subplot(gs[0:5])
ax2 = fig.add_subplot(gs[3:8])
ax3 = fig.add_subplot(gs[6:11])
ax4 = fig.add_subplot(gs[9:14])
rax1 = ax1.twinx()
axbg = fig.add_subplot(gs[:])

ax_list = [ax1,ax2,ax3, ax4, rax1, axbg]

label_list = ['x', 'rv', 'lag1ac', 'rf']
color_list = ['C0', 'C1', 'C2', 'C3']


for i,ax in enumerate(ax_list):
    make_patch_spines_invisible(ax)
    ax.set_xlim(0, t_f)
    
    ax.xaxis.set_visible(False)
    if i<4:
        ax.set_ylabel(label_list[i])
        for spine in ax.spines.values():
            spine.set_color(color_list[i])
            ax.yaxis.label.set_color(color_list[i])
            ax.tick_params(axis='y', colors=color_list[i])

ax1.plot(time, xx)
rax1.plot(time, alpha_ramp, color = 'slategray')
rax1.spines['right'].set_visible(True)
rax1.spines['right'].set_color('slategray')
ax1.spines['left'].set_visible(True)

ax2.yaxis.set_ticks_position('right')
ax2.yaxis.set_label_position('right')
ax2.spines['right'].set_visible(True)
ax2.spines["right"].set_position(("axes", 0.7))
ax2.yaxis.set_label_coords(0.75, 0.5)

# ax3.set_ylim(ax3.get_ylim()[::-1])
ax3.spines['left'].set_visible(True)
#ax3.set_ylim((7, 2))

ax4.yaxis.set_ticks_position('right')
ax4.yaxis.set_label_position('right')
# #ax4.set_ylim((5.5, 2))
ax4.spines['right'].set_visible(True)
ax4.spines['bottom'].set_visible(True)
ax4.spines['bottom'].set_color('k')
ax4.xaxis.set_visible(True)
ax4.spines["right"].set_position(("axes", 0.7))
ax4.yaxis.set_label_coords(0.75, 0.5)


for i, (ax, data) in enumerate(zip([ax2, ax3, ax4], [rv, rac, rrf])):
    a,b = np.polyfit(ct, data, 1)
    ax.plot(ct, data, color = color_list[i+1])
    ls = '--'
    if i<2:
        if pval[i] <0.05:
            ls = '-'
    else:
        if pval[i] >0.95:
            ls = '-'
            
    ax.plot(ct, ct * a + b,
            color = 'k',
            ls = ls,
            lw = 0.5)

axbg.axvline(t_crit, lw = 0.2, color ='k')
axbg.yaxis.set_visible(False)

#################################################################
# BAYESIAN FIT OF AN AR1 PROCESS                                #
#################################################################

# marginal_alpha, marginal_sigma = bayesian_AR1(time[mask], xx[mask], int(rw / dt), res = 10)

# fig, ax = plt.subplots()
# ax.contour(marginal_alpha)
