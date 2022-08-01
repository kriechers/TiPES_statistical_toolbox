import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
from scipy.signal import filtfilt, cheby1, lombscargle
from scipy.interpolate import interp1d, UnivariateSpline
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

#----------------------------------------------------------------
# data preprocessing according to NB

def cheby_lowpass(cutoff, fs, order, rp):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = cheby1(order, rp, normal_cutoff, btype='low', analog=False)
    return b, a

def cheby_lowpass_filter(x, cutoff, fs, order, rp):
    b, a = cheby_lowpass(cutoff, fs, order, rp)
    y = filtfilt(b, a, x)
    return y

def cheby_highpass(cutoff, fs, order, rp):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = cheby1(order, rp, normal_cutoff, btype='high', analog=False)
    return b, a

def cheby_highpass_filter(x, cutoff, fs, order, rp):
    b, a = cheby_highpass(cutoff, fs, order, rp)
    y = filtfilt(b, a, x)
    return y

# set parameters for chebychef filter:
order = 8
# the maximum ripple allowed below unity gain in the passband (in dB)
rp = .05
# sampling frequency, should be equal to 1 because we have 1-yr steps?
fss = 1.
# the sampling width to which we want to interpolate the data
samp = 5
# MTM time half-bandwidth
NW = 4

sec_factor = .95
sw = 100

print("samp: ", samp)
# the corresponding cutoff frequency
cutoff = .5 * 1. / samp
# offset at beginning and end:
offset = 19

int_method = 'cubic'
print(int_method)

prox = 'd18o'

data_file = 'original_data/icecore/NGRIP_d18O_and_dust_5cm/NGRIP_d18O_and_dust_5cm.xls'
colnames = ['depth', 'd18o', 'dust', 'age', 'age_err']

age_int = np.array(d18o_age, dtype = 'int')
end_age = 10277
start = np.argmin(np.abs(d18o_age - end_age))
end = -1
age = d18o_age[start : end]
age_int = age_int[start : end]
age_ipd = np.arange(age_int.min(), age_int.max())
# age_ipd = np.arange(age_int.min() + 5, age_int.max())
print(age_ipd)
#d18o = ngrip[:, 2][start : end]
#dust = ngrip[:, 3][start : end]
d18o = d18o[start : end]


raw_samp = np.diff(age)
print("sampling steps > 5yr: ", np.where(raw_samp > 5)[0].shape[0] / float(raw_samp.shape[0]))
print("sampling steps > 4yr: ", np.where(raw_samp > 4)[0].shape[0] / float(raw_samp.shape[0]))
bins = np.linspace(raw_samp.min(), raw_samp.max(), 20)
rs_hist = np.histogram(raw_samp, bins = bins, density = True)
print(rs_hist)


time = np.arange(age.min(), age.max(), samp)[::-1]

if prox == 'dust':
    dust = -np.log(dust)
    nan_idx = np.where(np.isnan(dust) == True)[0]
    dust_filled = dust.copy()
    for i in nan_idx:
        swwi = np.min((i, sw))
        swwe = np.min((age.shape[0] - i, sw))
        nonan_idx_temp = np.where(np.isnan(dust[i - swwi : i + swwe]) == False)[0]
        nan_idx_temp = np.where(np.isnan(dust[i - swwi : i + swwe]) == True)[0]
        jdat_temp = np.zeros((2, nonan_idx_temp.shape[0]))
        jdat_temp[0] = d18o[i - swwi : i + swwe][nonan_idx_temp]
        jdat_temp[1] = dust[i - swwi : i + swwe][nonan_idx_temp]
        kernel_temp = st.gaussian_kde(jdat_temp)
        x_temp = d18o[i]
        y_temp = np.linspace(dust[i - swwi : i + swwe][nonan_idx_temp].min(), dust[i - swwi : i + swwe][nonan_idx_temp].max(), 101)
        xx_temp, yy_temp = np.meshgrid(x_temp, y_temp)
        positions_temp = np.vstack([xx_temp.ravel(), yy_temp.ravel()])
        f_temp = np.reshape(kernel_temp(positions_temp).T, xx_temp.shape)
        dust_filled[i] = y_temp[np.argmax(f_temp, axis = 0)]



if int_method == 'cheby':
    dat = np.loadtxt('data/NGRIP_d18O_samp%d_Cheby1_iter_interpolated_order%d_rp%d_at.txt'%(samp, order, rp))[::-1][::samp]
elif int_method == 'cubic':
    # dat = np.loadtxt('data/NGRIP_d18O_samp%d_cubic_spline_interpolated_order%d_filt_at.txt'%(samp, order))[::-1][::samp]
    # dat = np.loadtxt('data/NGRIP_d18o_samp5_cubic_spline_interpolated_end10277_sw100.txt')
    if prox == 'd18o':

        f_d18o = UnivariateSpline(age, d18o, s = 0.)
        d18o_ipd_cub_temp = f_d18o(age_ipd)
        d18o_ipd_cub = cheby_lowpass_filter(d18o_ipd_cub_temp, sec_factor * cutoff, fss, order, rp)
        dat = d18o_ipd_cub[offset : -offset][::samp][::-1]
    elif prox == 'dust':
        f_dust = UnivariateSpline(age, dust_filled, s = 0.)
        dust_ipd_cub_temp = f_dust(age_ipd)
        dust_ipd_cub = cheby_lowpass_filter(dust_ipd_cub_temp, sec_factor * cutoff, fss, order, rp)
        dat = dust_ipd_cub[offset : -offset][::samp][::-1]

        f_d18o = UnivariateSpline(age, d18o, s = 0.)
        d18o_ipd_cub_temp = f_d18o(age_ipd)
        d18o_ipd_cub = cheby_lowpass_filter(d18o_ipd_cub_temp, sec_factor * cutoff, fss, order, rp)
        dat2 = d18o_ipd_cub[offset : -offset][::samp][::-1]
        print("Cor(d18o, dust) = ", np.corrcoef(dat, dat2)[0, 1])

time = age_ipd[offset : -offset][::samp][::-1]
print(time)
print(time.shape)
print(dat.shape)

time = age_ipd[offset : -offset][::samp][::-1]
print(time)
print(time.shape)
print(dat.shape)

# plt.figure()
# plt.plot(time, dat)
# plt.show()

#----------C-O-M-P-U-T-A-T-I-O-N------S-T-A-R-T-S------H-E-R-E------------------------------------------------------

# normalize by standard deviation (not necessary, but makes it easier
# to compare with plot on Interactive Wavelet page, at
# "http://paos.colorado.edu/research/wavelets/plot/"

# dat = (dat - np.mean(dat)) / np.std(dat, ddof=1)
mdat = np.mean(dat)
sdat = np.std(dat, ddof = 1)
# dat = (dat - mdat)
dat = (dat - mdat) / sdat



filt = 100
ft = 'cheby'
stdar = 1

s1 = 10
s2 = 50

scavf = 200 / samp

o = 2
psn = 2000

surr_mode = 'fourier'

lpf = 800

offset = 200 / samp

if filt > 0:
    if ft == 'cheby':
        filt_dat = cheby_highpass_filter(dat, .95 * 1. / filt, 1. / samp, 8, .05)
    elif ft == 'ma':
        filt_dat = dat - runmean(dat, filt / samp)

variance = np.std(dat, ddof=1)**2

d18o_age = time[::-1]
d18o = filt_dat[::-1]


#-----------------------------------------------------------------

#------------- Data Preprocessing according to KR ---------------

# dt = 5
# rw = 40 # given in data points => 200y
# d18o_bins = np.arange(int(d18o_age[0]-dt/2),
#                       int(d18o_age[-1] + dt),
#                       dt)

# d18o_age, d18o = binning(d18o_age, d18o, d18o_bins)

# dust_bins = np.arange(int(dust_age[0]-dt/2),
#                       int(dust_age[-1] + dt/2),
#                       dt)

# dust_age, dust = binning(dust_age, dust, dust_bins)

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