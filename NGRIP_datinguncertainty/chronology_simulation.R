# Simulates chronologies based on the fitted regression model.

chronology_simulation = function(object, nsims=10000, method="inla",
                                 store.means=FALSE,print.progress=FALSE){
  
  ## sample hyperparameters

  time.start = Sys.time()

  if(tolower(method) == "inla"){ #if INLA is used
    noise= object$.args$noise
    reg.model = object$.args$reg.model

    if(print.progress) cat("Simulating ",nsims, " hyperparameters from INLA posterior...",sep="")
    
    
  if(tolower(noise) %in% c(0,"ar(0)","ar0","iid","independent")){
    hypersamples = inla.hyperpar.sample(nsims,object$fitting$fit)
    object$simulation = list(sigma = 1/sqrt(hypersamples[,1]))
  }else if (tolower(noise) %in% c(1,"ar1","ar(1)")){
    hypersamples = inla.hyperpar.sample(nsims,object$fitting$fit)
    object$simulation = list(sigma = 1/sqrt(hypersamples[,1]), phi=hypersamples[,2])
    
  }else if (tolower(noise) %in% c(2,"ar2","ar(2)")){
    hypersamples = inla.hyperpar.sample(nsims,object$fitting$fit)
    p=2
    hypersamplesar2 = inla.hyperpar.sample(nsims,object$fitting$fit)
    phii = hypersamplesar2[, 2L:(2L+(p-1L))]
    phis = apply(phii, 1L, inla.ar.pacf2phi)
    object$simulation = list(sigma = 1/sqrt(hypersamples[,1]),phi1=phis[1,],phi2=phis[2,])
  }

    
    if(print.progress) cat(" completed!\n",sep="")
  }


  ## sample mean ("fixed") vector and ("stochastic") noise component

  if(tolower(method) == "inla"){
    if(print.progress) cat("Simulating mean vector from fitted coefficients...",sep="")
    time.startmean = Sys.time()
    
    ##sample fixed parameters first
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

    n=dim(object$data)[1]
    if(store.means) object$simulation$dmean = matrix(NA,nrow=n,ncol=nsims)
    time.endmean = Sys.time()
    if(print.progress) cat(" completed in ",difftime(time.endmean,time.startmean,units="secs")[[1]]," seconds!\n",sep="")

    object$simulation$age = matrix(NA,nrow=n,ncol=nsims)

    if(print.progress) cat("Simulating chronologies...\n",sep="")
    time.startage=Sys.time()
    
    
    is.biascorrected = TRUE
    for(i in 1:nsims){
      if(i %% 1000 == 0 && print.progress){
        cat("Age simulation ",i,"/",nsims,". Elapsed time: ",difftime(Sys.time(),time.startage,units="secs")[[1]]," seconds...","\n",sep="")
      }
      ## from fixed parameters compute fixed model component using 'meanmaker' function
      coefs = latentsamples[[i]]$latent
      dmeansim = meanmaker( coefs, reg.model, nevents=object$.args$nevents,data = object$data )
      
      
      
      if(store.means) object$simulation$dmean[,i] = dmeansim #store mean if we want
      
      ##sample noise component
        if(tolower(noise) %in% c(0,"iid","independent","ar0","ar(0)")){
          noisesim = rnorm(n,mean=0,sd=object$simulation$sigma[i])
          
        }else if(tolower(noise) %in% c(1,"ar1","ar(1)")){
          noisesim = arima.sim(n=n,list(ar=c(object$simulation$phi[i])),
                               sd = object$simulation$sigma[i]*sqrt(1-object$simulation$phi[i]^2))
          
        }else if(tolower(noise) %in% c(2,"ar2","ar(2)")){
          gamma0 = (1-phis[2,i])/((1+phis[2,i])*(1-phis[1,i]-phis[2,i])*(1+phis[1,i]-phis[2,i]))
          
          noisesim = arima.sim(n = n, list(ar = c( phis[1,i],phis[2,i])),
                               sd = object$simulation$sigma[i]*sqrt(1/gamma0))
        }
      
      
    ## Take cumulatives. If log transformation is used, transform back first
    if(object$.args$transform == "log"){
      object$simulation$age[,i] = object$preceeding$y0+ cumsum(exp(dmeansim+noisesim))
    }else{
      object$simulation$age[,i] = object$preceeding$y0+ cumsum(dmeansim+noisesim)
    }
    
      
    }
    }
    
     time.endage = Sys.time()
    if(print.progress) cat("Completed in ",difftime(time.endage,time.startage,units="secs")[[1]]," seconds!\n",sep="")
    object$time$simulation$mean = difftime(time.endmean,time.startmean,"secs")[[1]]
    object$time$simulation$age = difftime(time.endage,time.startage,"secs")[[1]]
    object$time$simulation$total = difftime(time.endage,time.start,units="secs")[[1]]
  

  return(object)
}

