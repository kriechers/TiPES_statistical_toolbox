"""Model and fitting related functions

This file contains all the functions related to the deterministic model
the probabilistic model and the fitting of the model to data
"""
import numpy as np
import emcee
from distributions import ar1ue_like, normal_like, gamma_like
from scipy.optimize import fmin
import warnings


def linear_ramp(t, t0=0.0, dt=1.0, y0=0.0, dy=1.0):
    """Linear Ramp Function

    This function describes the linear transition between two constant values.

    Parameter
    ---------
    t : np.ndarray
        Time variable
    t0 : float
        Start time of the ramp
    dt : float
        Transition length
    y0 : float
        Function value before the transition
    dy : float
        Hight of the transion

    Return
    ------
    y : np.ndarray
        Function values of the linear transiton
    """
    lt_t0 = t < t0
    gt_t1 = t > t0 + dt
    condlist = [lt_t0,
                ~np.logical_or(lt_t0, gt_t1),
                gt_t1]
    funclist = [lambda t: y0,
                lambda t: y0 + dy * (t - t0) / dt,
                lambda t: y0 + dy]
    y = np.piecewise(t, condlist, funclist)
    return y


def lnpost(theta, t, yobs):
    """Posterior log-probability function

    Log-posterior function of the linear ramp function with an AR(1) noise
    model.
    This function also includes the priors for all the parameters

    Parameter
    ---------
    theta : np.ndarray
        Transformed model paramters (t0, ln(dt), y0, dy, ln(sigma), ln(tau))
        where the first for are for the linear ramp and the two last parameters
        are the standard deviation and autocorrelation time of the AR(1) noise
        model
    t : np.ndarray
        Time variable
    yobs : np.ndarray
        Observations of the ramp. Must be the same length as t

    Returns
    -------
    lnp : flaot
        log-posterior probability of the parameters given the data
        lnp(theta|t, yobs)
    """
    # Unpack and transform parameters
    t0, y0, dy = theta[[0, 2, 3]]
    dt, tau, sigma = np.exp(theta[[1, 4, 5]])
    # Priors go here
    # Additional terms are needed because dt and tau are
    # sampled as log-transformed variables
    # this makes sure that the priors are realy the ones specified with
    # the respective distributions
    #
    # Start and stop of transition:
    tmin = np.min(t)   # start of data
    tmax = np.max(t)   # end of data
    if t0 < tmin or t0 > tmax or (t0 + dt) < tmin or (t0 + dt) > tmax:
        return -np.inf
    # if dt is smaller than ~1 days
    if dt <= 0.0027:
        return -np.inf
    if sigma == 0 or sigma > 10:
        return -np.inf
    lnp = normal_like(t0, 0.0, 50.0)                   # start time
    lnp += gamma_like(dt, 2.0, 0.02) + np.log(dt)      # transition length
    lnp += 1.0                                         # start height
    lnp += normal_like(dy, 0.0, 5.0)                   # step height
    lnp += gamma_like(tau, 2.5, 0.15) + np.log(tau)    # autocorrelation time
    # Ramp function
    ypred = linear_ramp(t, t0, dt, y0, dy)
    # likelihood function
    lnp += ar1ue_like(yobs - ypred, t, tau, sigma)
    return lnp


def neglnpost(theta, t, yobs):
    """Negative log-posterior, used for optimization
    Returns the negative of the log-posterior calculated with
    lnpost(tehta, t, yobs).

    For ducumentation see lnpost
    """
    lnp = lnpost(theta, t, yobs)
    if np.isfinite(lnp):
        return -1.0 * lnp
    else:
        return 1e25


def rmse(yobs, ypred):
    """Root mean squared error between two variables

    Parameter
    --------
    yobs : np.ndarray
        Observed values (same length as ypred)
    ypred : np.ndarray
        Predicted values (same length as yobs)

    Return
    ------
    rmse : float
        Root mean squared error between yobs and ypred
    """
    return np.sqrt(np.mean((yobs - ypred) ** 2))


def fit_rmse(t, y, p0=None):
    """Fit linear ramp to observation

    Uses RMSE minimization to fit a linear ramp to observations

    Parameter
    ---------
    t : np.ndarray
        Time variable
    y : np.ndarray
        Observations of the linear ramp
    p0 : None (default) or tuple of four parameters
        Starting parameters for the observation, if p0=None, than a starting
        position for the optimization will be guessed from the data

    Return
    ------
    p : np.ndarray
        Optimal parameter set for the linear ramp (t0, dt, y0, dy)

    See also
    --------
    linear_ramp : function that is fitted to the data
    """
    sort = np.argsort(t)
    if p0 is None:
        p0 = (0.0, np.log(10), np.mean(y[sort][:10]),
              np.mean(y[sort][-10:]) - np.mean(y[sort][:10]))
    p, *_, flag = fmin(lambda p: rmse(linear_ramp(t[sort], p[0], np.exp(p[1]),
                                                  p[2], p[3]), y[sort]),
                       p0, ftol=1e-5, xtol=1e-5, maxfun=1e5, maxiter=1e5,
                       disp=False, full_output=True)
    if flag != 0:
        warnings.warn('RMSE optimisation did not converge, returning guess',
                      RuntimeWarning)
        return p0
    return p


def fit_map(t, y, theta0=None, **kwargs):
    """Fit linear ramp to observations using MAP estimation

    Uses maximization of the log-posterior to fit a linear ramp to obsrvations
    Assumes an AR(1) noise model for the deviations from the ramp

    Parameter
    ---------
    t : np.ndarray
        Time variable
    y : np.ndarray
        Observations of the linear ramp
    theta0 : None (default) or np.ndarray of 6 parameters
        Starting parameters for the observation, if theta0=None,
        than a starting position for the optimization will be guessed from
        the data using some simple heuristics and fitting via RMSE minimization

    Return
    ------
    theta_map : np.ndarray
        Optimal transformed parameters for the linear ramp and the AR(1) noise
        model (t0, ln(dt), y0, dy, ln(sigma), ln(tau)).
        If the optimization fails, an array of np.nan will be returned

    See also
    --------
    lnpost : posterior function used for the fitting
    linear_ramp : model that is fitted
    fit_rmse : start parameter estimation
    """
    if theta0 is None:
        t0, lndt, y0, dy = fit_rmse(t, y)
        sigma = np.std(np.diff(y))
        theta0 = np.array((t0, lndt, y0, dy, np.log(2.0), np.log(sigma)))
    # Find high propability parameters by optimization
    theta_map, _, _, _, flag = fmin(neglnpost, theta0, args=(t, y),
                                    maxfun=1e5, maxiter=1e5,
                                    full_output=True, disp=False, **kwargs)
    if flag != 0:
        warnings.warn('MAP optimisation did not converge, returning nan',
                      RuntimeWarning)
        return theta0 * np.nan
    else:
        return theta_map


def fit_mcmc(t, y, theta0=None,
             nwalkers=60, nsample=20000, nthin=1, nburnin=5000):
    """Run MCMC sampler for linear ramp model with AR(1) noise

    This function sets up and runs an emcee.EnsembleSampler for the
    linear ramp model with AR(1) noise

    Parameter
    ---------
    t : np.ndarray
        Time variable of the observations
    y : np.ndarray
        Observations of the linear ramp
    theta0 : None (default) or np.ndarray of 6 parameters
        Starting parameters for the observation, if theta0=None,
        than a starting position for the optimization will be guessed from
        the data using some simple heuristics and fitting via RMSE minimization
    nwalkers : int
        Number of ensemble walkers used in the MCMC sampler
    nsample : int
        Number of samples drawn during the sample run
    nthin : int
        Thinning of the MCMC chains (number of samples per walker is
        nsample / nthin)
    nburnin : int
        Number of samples run before the nsample samples are drawn
        These samples are not saved

    Return
    ------
    sampler : emcee.EnsembleSampler
        Sampler object after the MCMC run.

    See also
    -------
    linear_ramp : deterministic part of the model
    lnpost : posterior function from which the samples are drawn from
    """
    # Sampler parameters
    ndim = 6
    # Fit MAP as start point
    mask = np.isfinite(y)
    if not np.any(mask):
        return None
    if theta0 is None:
        theta0 = fit_map(t[mask], y[mask])
    # Initialize sampler
    pos0 = theta0 + 0.01 * np.random.randn(nwalkers, ndim)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnpost,
                                    args=(t[mask], y[mask]))
    # Run burnin
    if nburnin > 0:
        pos0, *_ = sampler.run_mcmc(pos0, nburnin, store=False)
        sampler.reset()
    # Run sampler
    sampler.run_mcmc(pos0, nsample, thin=nthin)
    return sampler


def fit_gridsearch(t, y):
    """Fit ramp using a grid search

    Uses a brute force search of all possible starting and
    ending positions of the ramp.

    WARNING: If the data contains mainy observations this
    can take a long time.

    Parameter
    ---------
    t : np.ndarray
        Time variable of the observations
    y : np.ndarray
        Observations of the linear ramp

    Return
    ------
    p : np.ndarray
        Optimal parameter set for the linear ramp (t0, dt, y0, dy)

    See also
    --------
    linear_ramp : function that is fitted to the data
    """
    mask = np.isfinite(y)
    tm = t[mask]
    ym = y[mask]
    pars = (np.nan, np.nan, np.nan, np.nan)
    rmsmin = 1e15
    for i in range(len(tm)):
        t0 = t[i]
        y0 = np.mean(ym[:i + 1])
        for j in range(i + 1, len(tm)):
            dt = t[j] - t0
            dy = np.mean(ym[(j - 1):]) - y0
            ypred = linear_ramp(tm, t0, dt, y0, dy)
            rms = np.sqrt(np.mean((ym - ypred) ** 2))
            if rms < rmsmin:
                pars = (t0, dt, y0, dy)
                rmsmin = rms
    return pars
