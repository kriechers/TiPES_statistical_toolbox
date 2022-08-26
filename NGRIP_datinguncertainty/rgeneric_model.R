# Non-standard models in INLA requires specification using the rgeneric modeling framework. This includes key functions that provide information on precision matrix, mean vector, priors, graph etc. This is intended for internal use only, but the documentation is included here in case someone want to change something. See example below for how rgeneric can be used.

rgeneric_model = function( #specifies necessary functions for INLA to define the linear ramp model
  cmd = c("graph", "Q","mu", "initial", "log.norm.const", "log.prior", "quit"),
  theta = NULL)
{
  linramp = function(t,t0=0,dt=1,y0=0,dy=1){ #linear ramp function
    y = numeric(length(t))
    y = y0 + dy*(t-t0)/dt

    y[t<t0]=y0
    y[t>t0+dt]=y0+dy

    return(y)
  }


  interpret.theta = function() { #helpful function to transform back from internal parametrization
    y0 = theta[3]
    dy = (theta[4])
    y1 = y0+dy
    t0 = theta[1]
    Dt = exp(theta[2])
    t1=t0+Dt


    prec = exp(theta[6])
    tau = exp(theta[5])
    return(list(y0=y0,dy=dy,y1=y1,t0=t0,t1=t1,prec=prec,tau=tau,Dt=Dt))
  }

  mu = function() { #mean vector defined as a linear ramp
    param = interpret.theta()
    y0=param$y0; dy=param$dy; y1=param$y1; t0=param$t0; t1=param$t1; Dt = param$Dt; dt=Dt

    mvek = linramp(timepoints,t0=t0,dt=dt,y0=y0,dy=dy)

    return(mvek)
  }


  graph = function(){ #graphs of conditional dependence structure. 1 where Q[i,j] != 0, 0 where Q[i,j] = 0
    G = Q()
    G[G != 0] = 1
    return (G)
  }

  Q = function(){ #inverse covariance matrix: AR(1) that allow for unequal spacing
    param=interpret.theta()
    ii = 1:n
    jj = 1:n


    rhos=rep(NA,n)
    rhos[2:n] = exp(-diff(timepoints)/param$tau)
    kappa0 = param$prec

    xx = rep(NA,2*n-1)
    xx[1] = 1+rhos[2]^2/(1-rhos[2]^2)
    xx[2] = 1/(1-rhos[n]^2)
    xx[3:n] = 1/(1-rhos[2:(n-1)]^2) + rhos[3:n]^2/(1-rhos[3:n]^2)
    xx[(n+1):(2*n-1)] = -rhos[2:n]/(1-rhos[2:n]^2)

    xx = kappa0*xx

    i = c(1L, n, 2L:(n - 1L), 1L:(n - 1L))
    j = c(1L, n, 2L:(n - 1L), 2L:n)

    Q = Matrix::sparseMatrix(
      i = i,
      j = j,
      x = xx,
      symmetric = TRUE
    )


    return (Q)
  }

  log.norm.const = function(){ #INLA computes this automatically
    return(numeric(0))
  }

  log.prior = function(){
    params = interpret.theta()

    #log-priors are given for internal parametrisation using the change of variables theorem

    lprior = dnorm(theta[1],mean=round(0.5*(tslutt+tstart)),sd=50,log=TRUE) #t0
    lprior = lprior + dgamma(exp(theta[2]),shape=1.0,rate=0.02,log=TRUE) + theta[2] #dt
    lprior = lprior + dnorm(theta[3],mean=ystart,sd=5,log=TRUE) #y0
    lprior = lprior + dnorm(theta[4],mean=0,sd=10.0,log=TRUE) #dy
    lprior = lprior + dgamma(exp(theta[5]),2.5,rate = 0.15,log=TRUE) + theta[6] #tau/rho
    lprior = lprior + dgamma(exp(theta[6]),2,rate = 0.15,log=TRUE) + theta[6] #sigma/kappa
    return (lprior)
  }

  initial = function(){
    ini = c(0,0,0,0,1,2)
    return (ini)
  }

  quit = function(){
    return ()
  }
  if(is.null(theta)){
    theta = initial()
  }

  cmd = match.arg(cmd)
  val = do.call(cmd, args = list())
  return (val)
}
