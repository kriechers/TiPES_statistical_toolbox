from scipy.signal import filtfilt, cheby1, lombscargle
from scipy.signal import freqs
import matplotlib.pyplot as plt


# set parameters for chebychef filter:
# order = 8
# the higher n, the steeper the decline in amplification for frequencies
# above the cutt off.

# the maximum ripple allowed below unity gain in the passband (in
# dB) rp = .05
# rp = 0.05 guarantees, that the gain over the
# passed frequency band does not fall below 0.9 approx.

# fss := sampling frequency
# fss = 1/5

# data_40 = cheby_lowpass_filter(data, 1/40, fss, order, rp)
# data_40_100 = cheby_highpass_filter(data_40, 1/100, fss, order, rp)
# data_filt = data_40_100[offset: -offset]


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


# illustrating how the cutoff and the sampling rate work
# in combination

# b,a = cheby1(8, 0.05, 1/100, btype='high', analog=False, fs = 1/5)
# w, h = freqs(b, a)

# plt.semilogx(w, 20 * np.log10(abs(h)), lw = 3)
# plt.title('Chebyshev Type I frequency response (rp=0.05)')
# plt.xlabel('Frequency [radians / second]')
# plt.ylabel('Amplitude [dB]')
# plt.margins(0, 0.1)
# plt.grid(which='both', axis='both')
# plt.axvline(100, color='green') # cutoff frequency
# plt.axhline(-5, color='green') # rp

# b,a = cheby1(8, 0.05, 1/10, btype='high', analog=False)
# w, h = freqs(b, a)
# plt.semilogx(w, 20 * np.log10(abs(h)), color = 'C1')


