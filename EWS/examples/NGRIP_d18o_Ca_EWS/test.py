import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('../../')
from preprocessing import download, binning, NGRIP_stadial_mask
from EWS_functions import runvar

data = download()

d18o_idx = data['d18o'].dropna().index.values
d18o = data['d18o'].iloc[d18o_idx].values
d18o_age = data['age'].iloc[d18o_idx].values

dust_idx = data['dust'].dropna().index.values
dust = data['dust'].iloc[dust_idx].values
dust_age = data['age'].iloc[dust_idx].values

dt = 5
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

GS, transition_ages = NGRIP_stadial_mask(d18o_age)

### transforming age into time ### 
d18o_time = -d18o_age[::-1]
d18o = d18o[::-1]
transition_times = -transition_ages[::-1]
GS = GS[::-1]

### check if time series starts in stadial state ###

if GS[0]:
    shift = 0
else:
    shift = 1

rv_list = []

for t_i, t_f in zip(transition_times[shift::2], transition_times[1+shift::2]):
    mask = (d18o_time > t_i) & (d18o_time < t_f)
    xx = d18o_time[mask]
    yy = d18o[mask]
    test_xx, test_yy = runvar(xx, yy, 50, 10)
    rv_list.append((test_xx, test_yy))

fig, ax = plt.subplots()
ax.plot(d18o_time, d18o)
for i in range(0, len(transition_times[:-1]), 2):
    ax.axvspan(transition_times[i],
               transition_times[i+1],
               alpha=0.2,
               edgecolor=None,
               facecolor='slategray',
               zorder=1)

rax = ax.twinx()
for xx, yy in rv_list:
    rax.plot(xx,yy, color = 'C2')
    
    

