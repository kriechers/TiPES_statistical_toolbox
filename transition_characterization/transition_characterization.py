from distributions import sample_ar1
from model import fit_mcmc, linear_ramp, fit_rmse
import numpy as np
import matplotlib.pyplot as plt
import joblib as jl
import sys
sys.path.append('../../')


def combined_transition(time, traces):
    m = len(traces[0])
    n = len(time)
    ramps = np.zeros((m, n))
    for i, p in enumerate(zip(*traces)):
        ramps[i] = linear_ramp(time, p[0], p[1]-p[0], p[2], p[3])
    p5, p50, p95 = np.percentile(ramps, [5, 50, 95], axis=0)
    return p5, p50, p95


def find_trans_time(time, obs_norm, rw=None):
    n = len(time)
    if rw is None:
        rw = int(n/10)

    score = np.zeros(int(n-2*rw))
    for i in range(n-2*rw):
        score[i] = np.abs(np.mean(obs_norm[i:i+rw])
                          - np.mean(obs_norm[i+rw:i+2*rw]))

    return time[np.argmax(score)+rw], score


def estimate_transition(time,
                        obs,
                        ptrans=None,
                        pnoise=None,
                        nwalkers=60, nsamples=60000, nthin=600):

    # normalize the time series
    obs_norm = (obs - np.mean(obs)) / np.std(obs)

    if ptrans is None:
        # initial t0 is computed by maximizing the difference
        # between means of two neighbouring sliding windows
        t0, score = find_trans_time(time, obs_norm)

        # no easy way to estimate dt - taking 1/10 of the total
        # time of the time series should be of the right magnitude

        dt = 1/10 * (time[-1] - time[0])

        # initial dy is estimated by subtracting means after
        # and before t0
        dy = (np.mean(obs_norm[time > t0 + 0.5*dt])
              - np.mean(obs_norm[time < t0 - 0.5*dt]))

    else:
        # unpack the initial guess for the parameter
        t0 = ptrans[0]
        dt = ptrans[1]
        dy = ptrans[2]

    if pnoise is None:
        mask1 = time < t0 - 0.5*dt
        mask2 = time > t0 + 0.5*dt
        # sigma is computed as the mean of obs_norm std before and
        # after the estimated transition. Note that this sigma is
        # related to the noise level of the AR1 process used to
        # model fluctuations around the linear ramp by
        # sigma_eff = np.sqrt(sigma ** 2 * (1 - alpha ** 2))

        sigma = (np.std(obs_norm[mask1])/2 +
                 np.std(obs_norm[mask2])/2)

        # the autocorrelation time tau is defined as the inverse
        # restoring force of an OU process: tau = 1/gamma.
        # alpha = exp(-1/tau dt) of an AR1 process can be
        # estimated before and after the transition by
        # liniearly regressing the increments dx = x_i+1 - x_i
        # against the values of the time series x_i.
        # This must be done on either side of the transition

        dobs1 = np.diff(obs_norm[mask1])
        a1, b1 = np.polyfit(obs_norm[mask1][:-1], dobs1, 1)

        dobs2 = np.diff(obs_norm[mask2])
        a2, b2 = np.polyfit(obs_norm[mask2][:-1], dobs2, 1)

        alpha = 1 + (a1+a2)/2
        tau = -(time[1] - time[0])/np.log(alpha)

    # shift time to center the transition
    stime = time - t0

    # take ln of the dt, tau and sigma

    lndt = np.log(dt)
    lntau = np.log(tau)
    lnsigma = np.log(sigma)

    # define thet0
    theta0 = (t0, lndt, -1, dy, lntau, lnsigma)

    # run mcmc machinery

    sampler = fit_mcmc(stime, obs_norm, theta0, nwalkers=20, nsample=6000)

    # load every nthin'th sample from the walkers and reshape to
    # final dimensions
    trace = sampler.chain[:, ::nthin, :].reshape(-1, sampler.ndim).copy()

    # convert from sample space to meaningfull space
    trace[:, [1, 4, 5]] = np.exp(trace[:, [1, 4, 5]])

    t0_tr = trace[:, 0] + t0
    t1_tr = trace[:, 0] + trace[:, 1] + t0
    y0_tr = trace[:, 2] * np.std(obs) + np.mean(obs)
    dy_tr = trace[:, 3] * np.std(obs)

    return t0_tr, t1_tr, y0_tr, dy_tr
