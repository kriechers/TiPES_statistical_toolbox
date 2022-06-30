# Transition Characterization

The directory *transition_characterization* is comprised of this
readme.md, three python files and the *Eamples* directory.

The three python files combine to a Bayesian ramp fit, that can
be used to estimate the onset, duration, and amplitude of an
abrupt transition evident in a given time series in an
uncertainty sensitive manner. The algorithm is not designed for
the detection of abrupt transitions.

Please note that the method was published in conjunction with the
article 'Erhardt, T. et al. Decadal-scale progression of the
onset of Dansgaard-Oeschger warming events. Clim. Past 15,
811–825 (2019)' and which is available in its original form from
https://github.com/terhardt/DO-progression (last access:
28.06.22)

## python files

* **model.py**: this file was copied from
    https://github.com/terhardt/DO-progression. It contains the
    specification of the model which is fitted to the data in a
    Bayesian sense. For a description of the model, please see
    below. Furthermore, it contains the specification of the
    posterior distributions and a function which ultimately
    executes the MCMC sampling procedure.

* **distributions.py**: this file was copied from
    https://github.com/terhardt/DO-progression. It contains
    several probability distributions which combine to the
    posterior distribution.

* **transition_characterization.py**: this files contains the
    function *estimate_transition*, which is a top-level
    interface function. It takes only the time series under study
    as an input and returns MCMC-sampled values from the Bayesian
    posterior distribution, which serve as an
    uncertainty-sensitive estimate of the respecitve transition
    onset 't0', duration 'dt' and amplitude 'dy'. This function
    uses the functions defined in the other two .py files.
    Optionally, an initial guess for the model paramters can be
    specified. Also, the parameters of the MCMC sampler can be
    chosen upon a function call. 

    The function *combined_transition* can be used to further
    assess the sampled transition parameters: it computes
    piecewise linear ramp functions from for each parameter
    sample (t0, dt, y0, dy). Subsequently, for each point on the
    given time axis, it returns the 5th, the 50th and th 95th
    percentiles accross all piecewise linear ramps.

    The function *find_trans_time* is only used if not initial
    guess for the transition parameters is passed to the
    *estimate_transition* function. 

## the model 

*Theory: The Bayesian ramp-fit is based on the assumption, that
the time-series under investigation can be modelled as an AR1
process that fluctuates around a time dependent mean value. The
time dependency of the mean is modelled as a piecewise linear
ramp of the form:

          y0                      if t_i<t0
y(t_i) =  dy (t_i-t0) / dt + y0   if t0<t_i<t0+dt
          y0 + dy                 if t_i>t0+dt

The fluctuations around that mean values are mathematically 
described as an AR1 process: 

f(t_i) = alpha * f_i-1 + sqrt(sigma² * (1-alpha²)) * epsilon_i, 

with epsilon_i being a standard Gaussian random variable, alpha
being the AR(1) coefficient and sigma being the variance of the
AR(1) process (not of the gaussian noise). The algorithm works
with the autocorrelation time tau of the AR1 process which is
give by tau = -log(alpha) / delta t, alpha = exp(-1/tau *delta t)
where delta t denotes the sampling time step of the time series.

The overall model thus reads: 
ypred(t_i) = y(t_i) + f(t_i)

and is uniquely defined by the set of 6 parameters: 
(t0, dt, y0, dy, alpha, sigma)

## likelihood function
The likelihood function p(D|M) indicates the probability that a
stochastic model prediction exactly matches the observations:

ypred = yobs

and is thus given by: 

p(D|M) = product_i=1^N 1/sqrt(2 * pi * sigma_eff²) 
     exp(-(1/2) * (delta_i-delta_i-1 * alpha)² / sigma_eff²)
      * p(ypred_0 = yobs_0), 

with

sigma_eff = sqrt(sigma² * (1-alpha²)),

delta_i = yobs_i - y_i,

p(ypred_0 = yobs_0) = 1/sqrt(2 * pi * sigma_eff) 
                               exp(-(1/2) * (yobs_0)² / sigma_eff²).


## priors
the priors for the paramters are defined in the file model.py as
part of the function lnpost

*note that all probabilities are ln transformed for computational
reasons*


### start time t0
lnp = normal_like(t0, 0.0, 50.0)                
Gaussian distribution with std of 50 and mean 0

The *estimate_transition* function shifts the given time series
by the inital guess for t0, such that the time series is
approximately centered around the transition for the application
of the MCMC. The output is transformed back to the original time
axis. 

### transition length dt
lnp += gamma_like(dt, 2.0, 0.02) + np.log(dt)   
This is the same as gamma.logpdf(x, alpha = 2.0, scale=1/0.02) 
from scipy.stats

### start height y0 
lnp += 1.0        
uniform prior 

### step height dy                              
lnp += normal_like(dy, 0.0, 5.0)

### autocorrelation time tau
lnp += gamma_like(tau, 2.5, 0.15) + np.log(tau) 
gamma distribution weighted with a linear function 

### variance sigma
no prior defined, this is equivalent to a uniform prior. 
