import numpy as np
import statsmodels.api as sm
import scipy.stats as st
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def fourrier_surrogates(ts, ns):
    ts_fourier  = np.fft.rfft(ts)
    random_phases = np.exp(
        np.random.uniform(0, 2 * np.pi,
                          (ns, ts.shape[0] // 2 + 1)) * 1.0j)
    ts_fourier_new = ts_fourier * random_phases
    new_ts = np.real(np.fft.irfft(ts_fourier_new))
    return new_ts

def significance_test(ts, ns, mode = 'linear'):
    tlen = ts.shape[0]
    tau = st.linregress(np.arange(tlen), ts)[0]
    tsf = ts - ts.mean()
    nts = fourrier_surrogates(tsf, ns)
    stat = np.zeros(ns)
    tlen = nts.shape[1]
    if mode == 'linear':
        for i in range(ns):
            stat[i] = st.linregress(np.arange(tlen), nts[i])[0]
    elif mode == 'kt':
        for i in range(ns):
            stat[i] = st.kendalltau(np.arange(tlen), nts[i])[0]
    p = 1 - st.percentileofscore(stat, tau) / 100.
    return p

def runvar(xx, yy, rw, res):
    n=len(xx)
    rw_idx=np.arange(0,rw)
    step_idx =np.arange(0,n-rw, res)
    idx_matrix=rw_idx[:,None]+step_idx[None,:]
    rv=np.var(yy[idx_matrix],axis=0)
    rv_xx = xx[(step_idx + rw/2).astype(int)]
    return rv_xx, rv


# for the bias in the estimation of the autocorrelation please
# see: https://stats.stackexchange.com/questions/278703/how-to-
# estimate-the-autocorrelation-function

#def autocorr(x):
#    result = np.correlate(x, x, mode='full')
#    return result[int(result.size/2):]/result[int(result.size/2)]
def exp_dec(x, gamma):
    return np.exp(-gamma * x)

def EWS(xx,yy,rw,res, visual_ctrl = True):
    
    '''
    EWS computes the three early warning indicators 1) variance, 
    2) lag-1 autocorrelation and 3) the restoring force from data 
    yy given on a time axis xx over running windows of length rw
    given in data points. 
    '''
    # subtract the mean from the data
    yy = yy- np.mean(yy)
    dt = xx[1] - xx[0]
    n=len(xx)
    rw_idx=np.arange(0,rw)
    step_idx =np.arange(0,n-rw, res)
    m = len(step_idx)
    idx_matrix=rw_idx[:,None]+step_idx[None,:]

    # each column of rw_matrix represents one running window of
    # of the time series
    rw_matrix = yy[idx_matrix]

    # subtract a linear trend from each running window
    p = np.polyfit(rw_idx, rw_matrix, 1)
    rw_matrix = rw_matrix - (p[0][None,:] * rw_idx[:,None]
                             + p[1][None,:])

    running_variance=np.var(rw_matrix,axis=0)
    autocorr = np.zeros((int(rw/2), m))
    running_restoring_force = np.zeros(m)
    running_lag1ac = np.zeros(m)
    for i in range(m):
        autocorr[:,i] = sm.tsa.acf(rw_matrix[:,i], nlags = int(rw/2)-1)
        
        running_restoring_force[i],_ = curve_fit(exp_dec,
                                                 np.arange(0,int(rw/2)),
                                                 autocorr[:,i]) / dt
        running_lag1ac[i] = autocorr[1,i]

        # if visual_ctrl: 
        #     fig, ax = plt.subplots()
        #     ax.plot(autocorr[:,i])
        #     ax.plot(exp_dec(np.arange(0,int(rw/2)),
        #                     running_restoring_force[i]))
            
    new_xx = xx[(step_idx + rw/2).astype(int)]
    rv_sig = significance_test(running_variance, 10000)
    ac_sig = significance_test(running_lag1ac, 10000)
    rf_sig = significance_test(running_restoring_force, 10000)
    significance = (rv_sig, ac_sig, rf_sig)
    return running_variance, running_lag1ac, running_restoring_force, new_xx, significance
    

def runac1(xx,yy,rw,res):
    n=len(xx)
    rw_idx=np.arange(0,rw)
    step_idx =np.arange(0,n-rw, res)
    idx_matrix=rw_idx[:,None]+step_idx[None,:]
    rv=np.var(yy[idx_matrix],axis=0)
    rv_xx = xx[(step_idx + rw/2).astype(int)]



def fourrier_surrogates(ts, ns):
    ts_fourier  = np.fft.rfft(ts)
    random_phases = np.exp(np.random.uniform(0, 2 * np.pi, (ns, ts.shape[0] // 2 + 1)) * 1.0j)
    ts_fourier_new = ts_fourier * random_phases
    new_ts = np.real(np.fft.irfft(ts_fourier_new))
    return new_ts

def kendall_tau_test(ts, ns, tau, mode1 = 'fourier', mode2 = 'linear'):
    tlen = ts.shape[0]

    if mode1 == 'fourier':
        tsf = ts - ts.mean()
        nts = fourrier_surrogates(tsf, ns)
    elif mode1 == 'shuffle':
        nts = shuffle_surrogates(ts, ns)
    stat = np.zeros(ns)
    tlen = nts.shape[1]
    if mode2 == 'linear':
        for i in range(ns):
            stat[i] = st.linregress(np.arange(tlen), nts[i])[0]
    elif mode2 == 'kt':
        for i in range(ns):
            stat[i] = st.kendalltau(np.arange(tlen), nts[i])[0]
    p = 1 - st.percentileofscore(stat, tau) / 100.
    return p


def runmean(x, w):
   n = x.shape[0]
   xs = np.zeros_like(x)
   for i in range(w // 2):
      xs[i] = np.nanmean(x[: i + w // 2 + 1])
   for i in range(n - w // 2, n):
      xs[i] = np.nanmean(x[i - w // 2 + 1:])

   for i in range(w // 2, n - w // 2):
      xs[i] = np.nanmean(x[i - w // 2 : i + w // 2 + 1])
   return xs


def runstd(x, w):
   n = x.shape[0]
   xs = np.zeros_like(x)
   for i in range(w // 2):
      xw = x[: i + w // 2 + 1]
      xw = xw - xw.mean()
      if np.std(xw) > 0:
          lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
          p0 = lg[0]
          p1 = lg[1]
          xw = xw - p0 * np.arange(xw.shape[0]) - p1

          xs[i] = np.std(xw)
      else:
          xs[i] = np.nan
   for i in range(n - w // 2, n):
      xw = x[i - w // 2 + 1:]
      xw = xw - xw.mean()
      if np.std(xw) > 0:
          lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
          p0 = lg[0]
          p1 = lg[1]

          xw = xw - p0 * np.arange(xw.shape[0]) - p1


          xs[i] = np.std(xw)
      else:
          xs[i] = np.nan

   for i in range(w // 2, n - w // 2):
      xw = x[i - w // 2 : i + w // 2 + 1]
      xw = xw - xw.mean()
      if np.std(xw) > 0:
          lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
          p0 = lg[0]
          p1 = lg[1]
          xw = xw - p0 * np.arange(xw.shape[0]) - p1

          xs[i] = np.std(xw)
      else:
          xs[i] = np.nan

   return xs

def runlf(x, y, w):
   n = x.shape[0]
   xs = np.zeros_like(y)
   for i in range(w // 2):
      xw = x[: i + w // 2 + 1]
      yw = y[: i + w // 2 + 1]
      yw = yw - yw.mean()

      xs[i] = np.polyfit(xw, yw, 1)[0]

   for i in range(n - w // 2, n):
      yw = y[i - w // 2 + 1:]
      xw = x[i - w // 2 + 1:]
      yw = yw - yw.mean()
      xs[i] = st.linregress(xw, yw)[0]

   for i in range(w // 2, n - w // 2):
      yw = y[i - w // 2 : i + w // 2 + 1]
      xw = x[i - w // 2 : i + w // 2 + 1]
      yw = yw - yw.mean()
      xs[i] = st.linregress(xw, yw)[0]
   return xs


def runac(x, w):
   n = x.shape[0]
   xs = np.zeros_like(x)
   for i in range(w // 2):
      xw = x[: i + w // 2 + 1]
      xw = xw - xw.mean()
      if np.std(xw) > 0:
          lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
          p0 = lg[0]
          p1 = lg[1]
          xw = xw - p0 * np.arange(xw.shape[0]) - p1

          xs[i] = np.corrcoef(xw[1:], xw[:-1])[0,1]
      else:
          xs[i] = np.nan

   for i in range(n - w // 2, n):
      xw = x[i - w // 2 + 1:]
      xw = xw - xw.mean()
      if np.std(xw) > 0:
          lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
          p0 = lg[0]
          p1 = lg[1]

          xw = xw - p0 * np.arange(xw.shape[0]) - p1

          xs[i] = np.corrcoef(xw[1:], xw[:-1])[0,1]
      else:
          xs[i] = np.nan

   for i in range(w // 2, n - w // 2):
      xw = x[i - w // 2 : i + w // 2 + 1]
      xw = xw - xw.mean()
      if np.std(xw) > 0:
          lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
          p0 = lg[0]
          p1 = lg[1]

          xw = xw - p0 * np.arange(xw.shape[0]) - p1
          xs[i] = np.corrcoef(xw[1:], xw[:-1])[0,1]
      else:
          xs[i] = np.nan

   return xs



def run_fit_a(x, w):
  n = x.shape[0]
  xs = np.zeros_like(x)

  for i in range(w // 2):
     xw = x[: i + w // 2 + 1]
     xw = xw - xw.mean()
     lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
     p0 = lg[0]
     p1 = lg[1]

     xw = xw - p0 * np.arange(xw.shape[0]) - p1

     dxw = xw[1:] - xw[:-1]
     lg = st.linregress(xw[:-1], dxw)[:]
     a = lg[0]
     b = lg[1]


     xs[i] = a

  for i in range(n - w // 2, n):
     xw = x[i - w // 2 + 1:]
     xw = xw - xw.mean()
     lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
     p0 = lg[0]
     p1 = lg[1]

     xw = xw - p0 * np.arange(xw.shape[0]) - p1

     dxw = xw[1:] - xw[:-1]
     lg = st.linregress(xw[:-1], dxw)[:]
     a = lg[0]
     b = lg[1]
     xs[i] = a

  for i in range(w // 2, n - w // 2):
     xw = x[i - w // 2 : i + w // 2 + 1]
     xw = xw - xw.mean()

     lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
     p0 = lg[0]
     p1 = lg[1]

     xw = xw - p0 * np.arange(xw.shape[0]) - p1


     dxw = xw[1:] - xw[:-1]
     lg = st.linregress(xw[:-1], dxw)[:]
     a = lg[0]
     b = lg[1]

     xs[i] = a
  return xs


def run_fit_a_ar1(x, w):
  n = x.shape[0]
  xs = np.zeros_like(x)

  for i in range(w // 2):
     xs[i] = np.nan

  for i in range(n - w // 2, n):
     xs[i] = np.nan

  for i in range(w // 2, n - w // 2):
     xw = x[i - w // 2 : i + w // 2 + 1]
     xw = xw - xw.mean()

     p0, p1 = np.polyfit(np.arange(xw.shape[0]), xw, 1)
     xw = xw - p0 * np.arange(xw.shape[0]) - p1


     dxw = xw[1:] - xw[:-1]

     xw = sm.add_constant(xw)
     model = sm.GLSAR(dxw, xw[:-1], rho=1)
     results = model.iterative_fit(maxiter=10)

     a = results.params[1]

     xs[i] = a
  return xs
