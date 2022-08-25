
# Fits the linear regression model to observations.
modelfitter = function(object, method="inla",noise=NULL,
                       print.progress=FALSE,verbose=FALSE){

# 
   time.start = Sys.time()

  formula = object$formula
  if(print.progress) cat("Performing least squares fit...",sep="")
  
  fit = lm(object$ls.formula,object$data) #fits model using least squares
  
  dmean = fit$fitted.value
  resi = fit$residuals

  ## extract parameter estimates
  object$LS.fitting = list(fit=fit, params=list(meanvector=as.numeric(dmean)))
  
  if(tolower(object$.args$noise) %in% c("iid","independent","ar(0)")){
    sigma = sd(resi)
    object$LS.fitting$params[["sigma"]] = sigma
  }else if(tolower(object$.args$noise) %in% c("ar1",1,"ar(1)")){
    noisefit = arima(resi,order = c(1,0,0))
    phi = noisefit$coef[1]
    sigma = sd(resi)
    object$LS.fitting$params[["sigma"]] = sigma
    object$LS.fitting$params[["phi"]] = phi
    object$LS.fitting$noisefit=noisefit
  }else if(tolower(object$.args$noise) %in% c("ar2",2,"ar(2)")){
    noisefit = arima(resi,order = c(2,0,0))
    phi1 = noisefit$coef[1]
    phi2 = noisefit$coef[2]
    sigma = sd(resi)
    object$LS.fitting$params[["sigma"]] = sigma
    object$LS.fitting$params[["phi1"]] = phi1
    object$LS.fitting$params[["phi2"]] = phi2
    object$LS.fitting$noisefit=noisefit
  }
  time.ls=Sys.time()
  if(print.progress) cat(" completed!\n",sep="")
  #}

  ### move to separate function

  if(tolower(method) == "inla"){

    ## will use results from least squares fit as starting point in INLA optimization. Requires proper parametrization
    if(print.progress) cat("Performing INLA fit...\n",sep="")

    #set initial values for fixed parameters based on least squares 'fit'
    my.control.fixed = control.fixed.priors(reg.model, fit, object$.args$nevents)

    resi = object$LS.fitting$fit$residuals

    initialmodes = log(1/object$LS.fitting$params$sigma^2)
    
    if(tolower(object$.args$noise) %in% c(1,"ar1","ar(1)")){
      phi.ls = object$LS.fitting$params$phi
      initialmodes = c(initialmodes, log( (1+phi)/(1-phi) ))
      
    }else if(tolower(object$.args$noise) %in% c(2,"ar2","ar(2)")){
      phi1.ls = object$LS.fitting$params$phi1
      phi2.ls = object$LS.fitting$params$phi2
      
      initialmodes=c(initialmodes, log( (1+phi1/(1-phi2))/(1-phi1/(1-phi2)) ), log( (1+phi2)/(1-phi2) ) )
    }

    num.threads=1 #rgeneric can sometimes be more stable if more than one core is used
    
    object$data$idy=1:nrow(object$data) #create covariate for random effect in INLA
    
    ## fit using INLA
    inlafit = inla(object$formula, family="gaussian",data=object$data, 
                   control.family=list(hyper=list(prec=list(initial=12, fixed=TRUE))) ,
                   control.fixed=my.control.fixed,num.threads = num.threads,
                   control.compute=list(config=TRUE),verbose=FALSE,
                   control.inla=list(restart=TRUE,h=0.1), 
                   control.mode=list(theta=initialmodes,restart=TRUE)  )

    
    object$fitting = list(fit=inlafit)
    
    
    if(print.progress){
      cat("Computing remaining posteriors using Monte Carlo simulation...\n",sep="")
    }
    
    ## extract posteriors for hyperparameters
  posterior_sigma = inla.tmarginal(function(x)1/sqrt(x),inlafit$marginals.hyperpar$`Precision for idy`); zmarg_sigma=inla.zmarginal(posterior_sigma,silent=TRUE)
  object$fitting$hyperparameters = list(posteriors=list(sigma_epsilon=posterior_sigma))
  object$fitting$hyperparameters$results$sigma_epsilon = zmarg_sigma
  if(tolower(object$.args$noise)%in% c(1,"ar1","ar(1)") ){
    posterior_phi = inlafit$marginals.hyperpar$`Rho for idy`; zmarg_phi = inla.zmarginal(posterior_phi,silent=TRUE)
    object$fitting$hyperparameters$posteriors$phi = posterior_phi
    object$fitting$hyperparameters$results$phi = zmarg_phi
    
  }else if(object$.args$noise %in% c(2,"ar2","ar(2)") ){
    hypersamples = inla.hyperpar.sample(50000,inlafit)
    
    p=2
    phii = hypersamples[, 2L:(2L+(p-1L))]
    phis = apply(phii, 1L, inla.ar.pacf2phi)
    posterior_phi1 = cbind(density(phis[1,])$x,density(phis[1,])$y);colnames(posterior_phi1)=c("x","y"); zmarg_phi1 = inla.zmarginal(posterior_phi1,silent=TRUE)
    posterior_phi2 = cbind(density(phis[2,])$x,density(phis[2,])$y);colnames(posterior_phi2)=c("x","y"); zmarg_phi2=inla.zmarginal(posterior_phi2,silent=TRUE)
    object$fitting$hyperparameters$posteriors$phi1 = posterior_phi1
    object$fitting$hyperparameters$results$phi1 = zmarg_phi1
    object$fitting$hyperparameters$posteriors$phi2 = posterior_phi2
    object$fitting$hyperparameters$results$phi2 = zmarg_phi2
  }

    
    time.inla=Sys.time()
    elapsed.inla = difftime(time.inla,time.ls,units="secs")[[1]]
    if(print.progress) cat("INLA fit completed in ",elapsed.inla, " seconds!\n",sep="")

    object$time$fit$inla = elapsed.inla
    object$time$fit$ls = difftime(time.ls,time.start,units="secs")[[1]]
    object$time$fit$total = time.inla-time.start
  }

  return(object)
}
