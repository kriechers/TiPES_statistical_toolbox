"""Distribution and likelyhood functions """
import numpy as np
from scipy.special import gammaln


def normal_like(x, mu=0.0, sigma=1.0):
    """Normal log-likelyhood function

        lnp({x}|mu, sigma)
    """
    k = 2 * np.pi * sigma ** 2
    d = (x - mu) / sigma
    return -0.5 * np.nansum(np.log(k) + d ** 2)


def exp_like(x, beta):
    """Exponential log-likelihood"""
    return np.log(beta) - beta * np.sum(x)


def gamma_like(x, alpha, beta):
    """Gamma log-likelyhood

        lnp({x}|alpha, beta)

    Note
    ----
    This is the same as gamma.logpdf(x, alpha, scale=1.0/beta) from scipy.stats
    """
    lnp = (alpha * np.log(beta) - gammaln(alpha) +
           (alpha - 1) * np.sum(np.log(x)) - beta * np.sum(x))
    return lnp


def ar1_like(x, alpha, sigma=1.0):
    """AR(1) log-likelihood for evenly sampled series x

    Parameter
    ---------
    x : array
        Series of observations/noise
    alpha : float
        autocorrelation factor of AR(1) process (0 <=k < 1)
    sigma : float
        standard deviation of AR(1) process

    Returns
    -------
    lnp : array
        Log-likelyhood of the observations given the parameters of the
        AR(1) process \ln p({x}|k, sigma)
    """
    xim1 = x[:-1]
    xi   = x[1:]
    lnp = normal_like(xi, alpha * xim1, sigma=sigma)
    sigma_lim = np.sqrt(sigma**2 / (1 - alpha**2))
    lnp += normal_like(x[0], mu=0.0, sigma=sigma_lim)
    return lnp


def ar1ue_like(x, t, tau, sigma=1.0):
    """Log-likelihood of unevenly sampled AR(1) process

    Parameter
    ---------
    x : np.ndarray
        Observations / noise observed at time t
    t : np.ndarray
        Observation times (len(t) == len(x))
    tau : float
        Autocorrelation time
    sigma : float
        Standard deviation of x
    x0 : float
        Value of x[0]

    Returns
    -------
    lnp : np.ndarray
        Log-likelyhood of the observations given the parameters of the
        AR(1) process \ln p(x|t, tau, sigma)
    """
    alpha = np.exp(-np.abs(np.diff(t)) / tau)
    sigma_e = np.sqrt(sigma ** 2 * (1 - alpha ** 2))
    xim1 = x[:-1]
    xi = x[1:]
    lnp = normal_like(xi, alpha * xim1, sigma=sigma_e)
    lnp += normal_like(x[0], mu=0.0, sigma=sigma)
    return lnp


def sample_ar1(n, alpha, sigma=1.0, x0=0):
    """Generate AR(1) noise for evenely sampled series"""
    x = np.zeros(n)
    x[0] = x0 + sigma * np.random.randn()
    for i in range(1, n):
        x[i] = alpha * x[i - 1] + sigma * np.random.randn()
    return x


def sample_ar1ue(t, tau, sigma=1.0):
    """Sample AR(1) process at given times

    Produces a realization of unevenly spaced observation
    of an AR(1) process.

    Parameter
    ---------
    t : np.ndarray
        Times at which the AR(1) process is observed
    tau : float
        Autocorrelation time of the AR(1) process
    sigma : float
        Standard deviation of the observations
    x0 : float
        Value of x at t=0

    Return
    ------
    x : np.ndarray
        Values of the AR(1) process at time t (len(x) == len(t))
    """
    n = len(t)
    x = np.zeros(n)
    x[0] = sigma * np.random.randn()
    alpha = np.exp(-np.abs(np.diff(t)) / tau)
    sigma_e = np.sqrt(sigma ** 2 * (1 - alpha ** 2))
    for i, a, s in zip(range(1, n), alpha, sigma_e):
        x[i] = a * x[i - 1] + s * np.random.randn()
    return x
