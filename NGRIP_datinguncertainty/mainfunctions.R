## This is the main wrapper function for the base model of the dating uncertainties.
## It performs all analyses in the correct order
main = function(age,depth,proxy, events=NULL,nsims=10000, eventmeasure = "depth",
                reference.label=NULL,proxy.type="d18O",
                transform="identity",reg.model = list(
                  const=FALSE,depth1=FALSE,depth2=TRUE,proxy=TRUE,psi0=TRUE,psi1=TRUE),
                noise="ar1", method="inla", CI.type="quantiles",
                event.estimation = NULL, bias = NULL,store.everything=FALSE,print.progress=FALSE){

  time.start = Sys.time()
  main.call = sys.call(which=1)
  if(print.progress) cat("Initiating data formatting...",sep="")

  #prepare object by formatting dataset, writing formulastrings etc
  object = prepare(age,depth,proxy, events=events,nsims=nsims, eventmeasure = eventmeasure,
                   reg.model = reg.model,
                   noise=noise, method=method,reference.label=reference.label,transform=transform,proxy.type=proxy.type)
  if(print.progress) cat(" completed!\n",sep="")
  #fit the data, first by least squares, then by INLA (if specified)
  object = modelfitter(object, method=method,print.progress=print.progress)

  #produce samples from the chronologies
  object = chronology_simulation(object, nsims=nsims, method=method,store.means=store.everything,print.progress=print.progress)

  #compute posterior marginal mean, quantiles and other summary statistics
  object = simulationsummarizer(object,CI.type=CI.type,print.progress=print.progress)

  ##compute posterior marginal mean, quantiles and other summary statistics for synchronous time scale
  #object = synchronoussummarizer(object,CI.type=CI.type,print.progress=print.progress)

  #if event.estimation list object (containing specifications) is included, perform dating estimation
  if(!is.null(event.estimation)){
    #find onset depth posterior by fitting linear ramp model with INLA
    object = linrampfitter(object,interval=event.estimation$interval,
                           h=event.estimation$h,t1.sims=event.estimation$t1.sims,
                           rampsims=event.estimation$rampsims,label=event.estimation$label,
                           depth.reference = event.estimation$depth.reference,
                           print.progress=print.progress,log.ramp=event.estimation$log.ramp)

    #perform Monte Carlo simulations to produce samples for onset age of warming transition
    object = event_depth_to_age(object, nsims = nsims, print.progress=print.progress,label=event.estimation$label,age.reference = event.estimation$age.reference)
  }
  #if bias list object is included, perform this analysis
  if(!is.null(bias)){
    object = biased_chronologies(object,bias.model=bias$bias.model,biasparams = bias$biasparams,nsims=nsims,store.samples=store.everything)
  }
  time.total = difftime(Sys.time(), time.start,units="secs")[[1]]
  object$.args$call = main.call
  object$time$total = time.total

  if(print.progress){
    cat("Completed everything in ",time.total, "seconds\n",sep="")
  }
  return(object)
}



## This function generates chronologies that have been synchronized with the tie-points
## sampled above
syncsampler = function(object,nsims=10000, tie_depths=NULL,tiepointsims=NULL,
                       free_indexes=NULL,
                       tie_indexes=NULL){

  tie_indexes = which.index(tie_depths,object$data$z)
  if(is.null(tiepointsims)&& is.null(template)) stop("No tie-points specified")

  n = length(object$data$y)
  m = length(tie_indexes)
  free_n = n-m

  latentselection = list()
  if(reg.model$const) latentselection$`(Intercept)`=1
  if(reg.model$depth1) latentselection$z=1
  if(reg.model$depth2) latentselection$z2=1
  if(reg.model$proxy) latentselection$x = 1

  for(i in 2:object$.args$nevents){
    if(reg.model$psi0) latentselection[[paste0("a",i-1)]] = 1
    if(reg.model$psi1)latentselection[[paste0("c",i-1)]] = 1
  }
  latentsamples = inla.posterior.sample(nsims,object$fitting$fit,selection=latentselection,verbose=FALSE,add.names=FALSE)
  samples = matrix(NA,nrow=n,ncol=nsims)


  time.start = Sys.time()
  for(r in 1:nsims){
    sigma_sample = object$simulation$sigma[r]
    phi_sample = object$simulation$phi[r]

    if(noise=="ar1"){
      Qfull = Qmaker_ar1cum(n,sigma_sample,phi_sample)
    }else{
      #'noise' is the precision matrix of the layer differences Q_x
      Qfull = Qymaker(noise)
    }

    Qa = Qfull[-tie_indexes,-tie_indexes]
    Qab = Qfull[-tie_indexes,tie_indexes]
    if(m==1){
      Qab = as.matrix(Qab,ncol=1)
    }


    coefs = latentsamples[[r]]$latent
    dmeansim = meanmaker( coefs, reg.model, nevents=object$.args$nevents,data = object$data )

    y_mu = cumsum(c(object$preceeding$y0,dmeansim))[2:(n+1)]
    free_mu = y_mu[free_indexes]
    tie_mu = y_mu[tie_indexes]


    La = t(chol(Qa))
    b_temp = (-Qab%*%(tiepointsims[,r]-tie_mu))[,1]
    w = solve(La,b_temp)
    mu_temp = solve(t(La),w)[,1]
    z0 = rnorm(free_n)
    v = solve(t(La),z0)[,1]
    samples[free_indexes,r] = mu_temp + free_mu + v

    #mu_amidb = (free_mu - solve(Qa)%*%Qab%*%(skewsamples[r]-tie_mu))[,1]

    #samples[free_indexes,r] = Qsimmer(1,Qa,mu_amidb)
    if((r %% 1000) == 0){
      cat("Simulation ",r,"/",nsims,". Elapsed time: ",Sys.time()-time.start ,"\n",sep="")
    }
  }

  time.end = Sys.time()-time.start

  object$simulation$age = samples
  object$time$tiepoints = time.end

  return(object)
}
