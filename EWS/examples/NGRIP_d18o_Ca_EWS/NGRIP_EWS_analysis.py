import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import sys
sys.path.append('../../')
from preprocessing import download, binning, NGRIP_stadial_mask
from EWS_functions import runvar
from EWS_functions import EWS
from cheby import cheby_highpass_filter

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

        
def vmarker(x0, x1, ax0, ax1, **kwargs):
    xy0 = (x0, ax0.get_ylim()[0])
    xy1 = (x1, ax1.get_ylim()[1])
    ax0.axvline(x0, **kwargs)
    ax1.axvline(x1, **kwargs)
    con = ConnectionPatch(xy0, xy1, 'data', 'data',
                          ax0, ax1, **kwargs)
    ax0.add_artist(con)




data = download()

d18o_idx = data['d18o'].dropna().index.values
d18o = data['d18o'].iloc[d18o_idx].values
d18o_age = data['age'].iloc[d18o_idx].values

dust_idx = data['dust'].dropna().index.values
dust = data['dust'].iloc[dust_idx].values
dust_age = data['age'].iloc[dust_idx].values

dt = 5
rw = 40 # given in data points => 200y
d18o_bins = np.arange(int(d18o_age[0]-dt/2),
                      int(d18o_age[-1] + dt),
                      dt)

d18o_age, d18o = binning(d18o_age, d18o, d18o_bins)

dust_bins = np.arange(int(dust_age[0]-dt/2),
                      int(dust_age[-1] + dt/2),
                      dt)

dust_age, dust = binning(dust_age, dust, dust_bins)


#################################################################
# analysis of the d18o                                          #
#################################################################

GS, transition_ages, transition_names = NGRIP_stadial_mask(d18o_age)

### transforming age into time ###

d18o_time = -d18o_age[::-1]
d18o = d18o[::-1]
transition_times = -transition_ages[::-1]
GS = GS[::-1]

### apply 100 year highpass filter to the entire time series ###

d18o_100 = cheby_highpass_filter(d18o,1/100, 1/5, 8, 0.05)

### check if time series starts in stadial state ###

if GS[0]:
    shift = 0
else:
    shift = 0

    
#check

# fig, ax = plt.subplots()
# ax.plot(d18o_time[GS], d18o[GS], lw = 0.5)
# ax.plot(d18o_time[~GS], d18o[~GS], color = 'C1',
#         lw = 0.5)
# for t in transition_times:
#     ax.axvline(t)

rv_list = []
rac_list = []
rrf_list = []
significance_list = []
centered_time_list = []

#################################################################
# run EWS analysis and store results and corresponidng lists    #
#################################################################

for t_i, t_f in zip(transition_times[shift:-1:2], transition_times[shift+1::2]):
#for t_i, t_f in zip([transition_times[15]], [transition_times[16]]):
    mask = (d18o_time > t_i) & (d18o_time < t_f)
    xx = d18o_time[mask]
    yy = d18o_100[mask]
    if len(xx) < 2* rw:
        continue
    rv, rac, rrf, centered_time, significance = EWS(xx,yy,rw,1)
    #test_xx, test_yy = runvar(xx, yy, 40, 1)
    rv_list.append(rv)
    rac_list.append(rac)
    rrf_list.append(rrf)
    centered_time_list.append(centered_time)
    significance_list.append(significance)


fig = plt.figure(figsize=(16,9))
height_ratios = [1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,2,2,2,2]
gs = fig.add_gridspec(16, 1,
                      hspace=0.5,
                      left=0.1,
                      bottom=0.15,
                      top=0.85,
                      right=0.9,
                      height_ratios=height_ratios)

ax1 = fig.add_subplot(gs[2:7])
ax2 = fig.add_subplot(gs[5:10])
ax3 = fig.add_subplot(gs[8:13])
ax4 = fig.add_subplot(gs[11:16])
axe = fig.add_subplot(gs[0])
axa = fig.add_subplot(gs[2:])

ax_list = [ax1,ax2,ax3, ax4, axe,axa]

label_list = ['d18o', 'rv', 'lag1ac', 'rf']
color_list = ['C0', 'C1', 'C2', 'C3']

#################################################################
# plot NGRIP data                                               #
#################################################################

for i,ax in enumerate(ax_list):
    make_patch_spines_invisible(ax)
    ax.set_xlim(-60000, -10000)
    if i<4:
        ax.xaxis.set_visible(False)
        ax.set_ylabel(label_list[i])
        for spine in ax.spines.values():
            spine.set_color(color_list[i])
            ax.yaxis.label.set_color(color_list[i])
            ax.tick_params(axis='y', colors=color_list[i])


ax1.plot(d18o_time, d18o)
    

for rv, rac, rrf, centered_time, significance in zip(rv_list, rac_list, rrf_list, centered_time_list, significance_list):

    for i, (ax, data) in enumerate(zip([ax2, ax3, ax4], [rv, rac, rrf])):
        a,b = np.polyfit(centered_time, data, 1)
        ax.plot(centered_time, data, color = color_list[i+1])
        if significance[i] <0.05:
            ls = '-'

        else:
            ls = '--'
        ax.plot(centered_time, centered_time * a + b,
                color = 'k',
                ls = ls,
                lw = 0.5)

    

axe.yaxis.set_visible(False)
axe.xaxis.set_label_position('top')
axe.xaxis.set_ticks_position('top')
axe.spines['top'].set_visible(True)
axe.set_xticks(np.arange(0,len(transition_names[::2])))
axe.set_xlim((0, len(transition_names[::2])))
axe.set_xticklabels(transition_names[::-2],
                    rotation=-45)

make_patch_spines_invisible(axa)
axa.yaxis.set_visible(False)
axa.set_xlim(ax1.get_xlim())
axa.spines['bottom'].set_visible(True)
axa.set_xlabel('GICC05 Age [ky b2k]')
axa.set_ylim((0, 1.05))

####################################################
# all spines invisible, now make individual spines #
# visible and customize outer appearence           #
####################################################

ax1.set_ylabel('$\delta^{18}$O$\;$[\u2030]')
#ax1.yaxis.set_label_coords(-0.04, 0.7)
#ax1.spines['left'].set_position(('axes', 0.035))
ax1.spines['left'].set_visible(True)

# ax2.set_ylim(ax2.get_ylim()[::-1])
ax2.yaxis.set_ticks_position('right')
ax2.yaxis.set_label_position('right')
ax2.spines['right'].set_visible(True)
#ax2.spines["right"].set_position(("axes", 0.95))
#ax2.yaxis.set_label_coords(1.04, 0.5)

# ax3.set_ylim(ax3.get_ylim()[::-1])
ax3.spines['left'].set_visible(True)
#ax3.set_ylim((7, 2))

# ax4.yaxis.set_ticks_position('right')
# ax4.yaxis.set_label_position('right')
# #ax4.set_ylim((5.5, 2))
# ax4.spines['right'].set_visible(True)


for i, a in enumerate(transition_times[1::2]):
    vmarker(i, a, axe, axa,
            lw=0.5,
            ls='solid',
            zorder=-100,
            color='lightsteelblue')

# for i in range(0, len(transition_times[:-1]), 2):
#     axa.axvspan(transition_times[i],
#                 transition_times[i+1],
#                 alpha=0.2,
#                 edgecolor=None,
#                 facecolor='slategray',
#                 zorder=1)
    


# check significance-test_xx

# rv = rv_list[12]
# rac = rac_list[12]
# rac_a, rac_b = np.polyfit(np.arange(0,len(rac)), rac, 1)

# test = fourrier_surrogates(rac, 10000)
# a, b = np.polyfit(np.arange(0,test.shape[1]), test.T, 1)

# fig, ax = plt.subplots()
# ax.plot(rac)
# ax.plot(np.arange(0,len(rac)), rac_a * np.arange(0,len(rac)) + rac_b)
# for i in range(10):
#     ax.plot(test[i], color = 'gray', lw = 0.2)
#     ax.plot(np.arange(0,test.shape[1]),
#             np.arange(0,test.shape[1]) * a[i] +b[i],
#             color = 'k',
#             ls = '--',
#             lw = 0.5)

# sig = np.sum(a > rac_a) / 10000