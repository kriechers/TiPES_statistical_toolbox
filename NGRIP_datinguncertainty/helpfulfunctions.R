# Contains functions used by other functions

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


# This function samples from a skew normal distribution, which is essentially two Gaussians of
# equal mean, but different variance merged together at the mean.
skewsampler = function(n,mean=0,sdL=1,sdU=1,log=FALSE,plothist=list(compute=FALSE,breaks=50,col="orange")){
  uppers = runif(n) < sdU/(sdL+sdU)
  samples = numeric(n)
  for(i in 1:n){
    if(uppers[i]){
      sample = abs(rnorm(1))
      samples[i] = sdU*sample+mean
    }else{
      sample = -abs(rnorm(1))
      samples[i] = sdL*sample+mean
    }
  }
  if(log) samples = log(samples)

  return(samples)
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

  object$simulation$age_sync = samples
  object$time$tiepoints = time.end

  return(object)
}

meanmaker = function(coefs,reg.model,nevents=69,data){
  ## Computes mean vector from given fixed effects 'coefs'.
  ## Requires specification of which effects to include ('reg.model'), the number of climate transitions ('nevents') and a data.frame with covariates ('data')
  coefcounter=1
  fitted=numeric(dim(data)[1])
  if(reg.model$const){
    fitted=coefs[1]
    coefcounter=coefcounter+1
  }
  if(reg.model$depth1){
    fitted = fitted + coefs[coefcounter]*data$z
    coefcounter=coefcounter+1
  }
  if(reg.model$depth2){
    fitted = fitted + coefs[coefcounter]*data$z2
    coefcounter=coefcounter+1
  }
  if(reg.model$proxy){
    fitted=fitted + coefs[coefcounter]*data$x
    coefcounter=coefcounter+1
  }
  if(nevents>0){
    for(i in 2:nevents){
      if(reg.model$psi0){

        fitted = fitted + coefs[coefcounter]*data[[paste0("a",i-1)]]
        coefcounter=coefcounter+1
      }
      if(reg.model$psi1){

        fitted = fitted + coefs[coefcounter]*data[[paste0("c",i-1)]]
        coefcounter=coefcounter+1
      }
    }
  }
  return(fitted)
}

# Computes the (noiseless) linear ramp function.
linramp = function(t,t0=0,dt=1,y0=0,dy=1){
  y = numeric(length(t))
  y = y0 + dy*(t-t0)/dt
  y[t<t0]=y0
  y[t>t0+dt]=y0+dy
  return(y)
}


# Computes the (noiseless) linear ramp function, but in reverse.
linramprev = function(t,t0=0,dt=1,y0=0,dy=1){
  y = numeric(length(t))
  y = y0 + dy*(t-t0)/dt

  y[t>t0]=y0
  y[t<t0+dt]=y0+dy
  return(y)
}


# Finds the indices where 'events' best match values in a given 'record'.
which.index = function(events, record){
  eventindexes = numeric(length(events))
  for(i in 1:length(events)){
    if(events[i] < min(record) || events[i] > max(record)){ #Gives NA if located outside range of 'record'
      warning(paste0("Event ",i,", located at ",events[i]," is outside the interval covered by 'record' (",min(record),", ",max(record),"). The event will be omitted!"))
      eventindexes[i] = NA
    }else{
      eventindexes[i] = which(abs(events[i]-record) == min(abs(events[i]-record)))
    }

  }
  #eventindexes = unique(c(1,eventindexes[!is.na(eventindexes)])) #Placing transition at the start of record. Removing NA and duplicates
  return(eventindexes)
}

# Computes posterior marginal mean and uncertainty intervals from simulations.
simulationsummarizer = function(object,CI.type="hpd",print.progress=FALSE){
  if(print.progress) cat("Computing posterior marginal mean and 95% credible intervals from chronology samples...\n",sep="")
  time.start = Sys.time()
  n = dim(object$simulation$age)[1]
  nsims = dim(object$simulation$age)[2]

  meanvek = rowMeans2(object$simulation$age)
  sdvek = sqrt(rowVars(object$simulation$age))

  lower = numeric(n); upper = numeric(n)
  if(CI.type=="hpd"){
    modevek = numeric(n)
    for(i in 1:n){
      dens = density(object$simulation$age[i,])
      modevek[i]=dens$x[which(dens$y == max(dens$y))]

      lower[i] = inla.hpdmarginal(0.95,dens)[1]
      upper[i] = inla.hpdmarginal(0.95,dens)[2]

    }
  }else{
    lower = meanvek-1.96*sdvek
    upper = meanvek+1.96*sdvek
  }

  time.summary = Sys.time()
  object$simulation$summary = list(mean=meanvek,sd=sdvek,lower=lower,upper=upper,
                                   .args=list(interval=cbind(lower,upper),print.progress=print.progress,CI.type=CI.type))
  if(CI.type=="hpd") object$simulation$summary$mode = modevek
  if(print.progress) cat(" completed in ",difftime(time.summary,time.start,units="secs")[[1]],"\n",sep="")

  object$time$samplesummary = list(total=difftime(time.summary,time.start,units="secs")[[1]])
  object$simulation$summary$sim.sum.time = difftime(time.summary,time.start,units="secs")[[1]]

  return(object)
}



## sets initial values for fixed parameters equal to least squares 'fit'
control.fixed.priors = function(reg.model, fit, nevents){

  my.control.fixed = list(mean=list(  ))

  if(reg.model$depth1) my.control.fixed$mean[["z1"]] = fit$coefficients[["z1"]]
  if(reg.model$depth2) my.control.fixed$mean[["z2"]] = fit$coefficients[["z2"]]
  if(reg.model$proxy) my.control.fixed$mean[["x"]] = fit$coefficients[["x"]]

  if(reg.model$psi0 || reg.model$psi1){
    for(i in 1:(nevents-1)){
      if(reg.model$psi0){
        my.control.fixed$mean[[paste0("a",i)]] = fit$coefficients[[paste0("a",i)]]
      }
      if(reg.model$psi0){
        my.control.fixed$mean[[paste0("c",i)]] = fit$coefficients[[paste0("c",i)]]
      }
    }
  }

  return(my.control.fixed)
}



# Simulates from a Gaussian process given the precision matrix (and mean vector)
Qsimmer = function(nsims, Q, muvek){
  nn=dim(Q)[1]
  samples = matrix(NA,nrow=nn,ncol=nsims)

  La = chol(Q)

  for(i in 1:nsims){
    samples_z = rnorm(nn)
    v = solve(La,samples_z)[,1]
    samples[,i] = muvek[(length(muvek)-nn+1):length(muvek)] + v
  }
  return(samples)
}

# Construct the precision matrix for the partial sum AR(1) process
Qmaker_ar1cum = function(n,sigma,rho){

  noise = sigma^2*(1-rho^2)
  ii = c(1L, n, 2L:(n - 1L), 1L:(n - 1L),1L:(n-2L)); jj = c(1L, n, 2L:(n - 1L), 2L:n,3L:n)
  if(n==3){
    ii = c(1L, n, n-1L,  1L:(n - 2L), n-1L, 1L:(n-2L))
    jj = c(1L, n, n-1L,  2L:(n - 1L), n    ,3L:n)
    values = 1/(noise)*c(2+2*rho+1*rho^2, 1, rep(2+2*rho+1*rho^2,n-2L), rep(-(1+2*rho+rho^2),n-2),-(1+rho), rep(rho,n-2))
  }else{
    ii = c(1L, n, n-1L, 2L:(n-2L), 1L:(n-2L), n-1L, 1L:(n-2L))
    jj = c(1L, n, n-1L, 2L:(n-2L), 2L:(n-1),  n,    3L:n)
    values = 1/(noise)*c( 2+2*rho+rho^2, 1, rho^2+2*rho+2,rep(2+2*rho+2*rho^2,n-3L), rep(-(1+2*rho+rho^2),n-2),-(1+rho), rep(rho,n-2)  )
  }
  return(sparseMatrix(i=ii, j=jj, x=values, giveCsparse = TRUE, symmetric = TRUE))
}

# Construct the precision matrix Qy of the partial sum of a Gaussian process with given precision matrix Qx
Qymaker = function(Qx){
  nn = dim(Qx)[1]
  ii = c(1:nn,2:nn); jj = c(1:nn,1:(nn-1)); xx = c(rep(1,nn),rep(-1,nn-1))
  S = sparseMatrix(i=ii,j=jj,x=xx)
  return(t(S)%*%Qx%*%S)
}


## This function imports and formats the Adolphi (2018) and Muscheler (2020) tie-points
adolphi_tiepoint_pdfs = function(plotdens=FALSE,tieshifts = numeric(5),x.ref=NULL){

  tiepointdata1 = scan("datasets_used/Adolphi18_pdf_tie_1-3.txt",skip=6,nlines=1,quiet=TRUE)
  tie1_range = tiepointdata1[1:2]
  tie1_y = tiepointdata1[3:303]
  tie1_x = c((-150):150) + tieshifts[1]
  tie1 = cbind(tie1_x,tie1_y)

  tiepointdata2 = scan("datasets_used/Adolphi18_pdf_tie_1-3.txt",skip=7,nlines=1,quiet=TRUE)
  tie2_range = tiepointdata2[1:2]
  tie2_y = tiepointdata2[3:303]
  tie2_x = c((-150):150)+ tieshifts[2]
  tie2 = cbind(tie2_x,tie2_y)

  tiepointdata3 = scan("datasets_used/Adolphi18_pdf_tie_1-3.txt",skip=8,nlines=1,quiet=TRUE)
  tie3_range = tiepointdata3[1:2]
  tie3_y = tiepointdata3[3:303]
  tie3_x = c((-150):150)+ tieshifts[3]
  tie3 = cbind(tie3_x,tie3_y)


  tiepointdata4 =  scan("datasets_used/Adolphi18_pdf_tie_4_update_Muscheler2020.txt",skip=6,quiet=TRUE)
  tie4_range = tiepointdata4[1:2]
  tie4_y = tiepointdata4[3:2003]
  tie4_x = (-1000):1000+ tieshifts[4]
  tie4 = cbind(tie4_x,tie4_y)

  tiepointdata5 = scan("datasets_used/Adolphi18_pdf_tie_5_update_Muscheler2020.txt",skip=6,quiet=TRUE)
  tie5_range = tiepointdata5[1:2]
  tie5_y = tiepointdata5[3:3003]
  tie5_x = (-1500):1500+ tieshifts[5]
  tie5 = cbind(tie5_x,tie5_y)

  if(plotdens){
    la=layout(mat=matrix(c(1,4,1,4,2,4,2,5,3,5,3,5) ,nrow=2))
    if(sum(tieshifts)==0){
      xlab = "x - GICC05 (years)"
    }else{
      xlab = "Age (y b2k)"
    }

    plot(tie1,type="l",main="Tie-point 1",xlab=xlab,ylab="Density")

    zmarg = inla.zmarginal(tie1,silent=TRUE)
    abline(v=zmarg$mean)
    abline(v=zmarg$quant0.025,col="gray")
    abline(v=zmarg$quant0.975,col="gray")
    if(!is.null(x.ref)){
      abline(v=x.ref[1],col="blue")
    }
    plot(tie2,type="l",main="Tie-point 2",xlab=xlab,ylab="Density")

    zmarg = inla.zmarginal(tie2,silent=TRUE)
    abline(v=zmarg$mean,lwd=0.8)
    abline(v=zmarg$quant0.025,col="gray",lwd=0.8)
    abline(v=zmarg$quant0.975,col="gray",lwd=0.8)
    if(!is.null(x.ref)){
      abline(v=x.ref[2],col="blue")
    }
    plot(tie3,type="l",main="Tie-point 3",xlab=xlab,ylab="Density")

    zmarg = inla.zmarginal(tie3,silent=TRUE)
    abline(v=zmarg$mean,lwd=0.8)
    abline(v=zmarg$quant0.025,col="gray",lwd=0.8)
    abline(v=zmarg$quant0.975,col="gray",lwd=0.8)

    if(!is.null(x.ref)){
      abline(v=x.ref[3],col="blue")
    }

    plot(tie4,type="l",main="Tie-point 4",xlab=xlab,ylab="Density")

    zmarg = inla.zmarginal(tie4,silent=TRUE)
    abline(v=zmarg$mean,lwd=0.8)
    abline(v=zmarg$quant0.025,col="gray",lwd=0.8)
    abline(v=zmarg$quant0.975,col="gray",lwd=0.8)
    if(!is.null(x.ref)){
      abline(v=x.ref[4],col="blue")
    }
    plot(tie5,type="l",main="Tie-point 5",xlab=xlab,ylab="Density")

    zmarg = inla.zmarginal(tie5,silent=TRUE)
    abline(v=zmarg$mean,lwd=0.8)
    abline(v=zmarg$quant0.025,col="gray",lwd=0.8)
    abline(v=zmarg$quant0.975,col="gray",lwd=0.8)

    if(!is.null(x.ref)){
      abline(v=x.ref[5],col="blue")
    }
    par(mfrow=c(1,1))
  }


  tie_pdfs = list(tie1,tie2,tie3,tie4,tie5)
  return(tie_pdfs)
}

# This function produces samples from the Adolphi (2018) and Muscheler (2020) samples
adolphi_tiepoint_simmer = function(nsims=10000,tieshifts = numeric(5), plotdens=FALSE,x.ref=NULL, plothist=list(compute=FALSE,breaks=50,col="orange")){

  tie_pdfs = adolphi_tiepoint_pdfs(plotdens=plotdens,tieshifts = tieshifts,x.ref=x.ref)


  tie_samples = matrix(NA,ncol=nsims,nrow=5)

  for(i in 1:5){
    sims = inla.rmarginal(nsims,tie_pdfs[[i]])
    tie_samples[i,] = sims
  }

  if(sum(tieshifts)==0){
    xlab = "x - GICC05 (years)"
  }else{
    xlab = "Age (y b2k)"
  }

  if(!is.null(plothist) && plothist$compute==TRUE){
    la=layout(mat=matrix(c(1,4,1,4,2,4,2,5,3,5,3,5) ,nrow=2))
    hist(tie_samples[1,],freq=0,col=plothist$col,breaks=plothist$breaks,main="Tie-point 1",xlab=xlab)
    hist(tie_samples[2,],freq=0,col=plothist$col,breaks=plothist$breaks,main="Tie-point 2",xlab=xlab)
    hist(tie_samples[3,],freq=0,col=plothist$col,breaks=plothist$breaks,main="Tie-point 3",xlab=xlab)
    hist(tie_samples[4,],freq=0,col=plothist$col,breaks=plothist$breaks,main="Tie-point 4",xlab=xlab)
    hist(tie_samples[5,],freq=0,col=plothist$col,breaks=plothist$breaks,main="Tie-point 5",xlab=xlab)
    par(mfrow=c(1,1))
  }


  #par(mfrow=c(1,1))
  return(tie_samples)
}






# This function samples from a skew normal distribution, which is essentially two Gaussians of
# equal mean, but different variance merged together at the mean.
skewsampler = function(n,mean=0,sdL=1,sdU=1,log=FALSE,plothist=list(compute=FALSE,breaks=50,col="orange")){
  uppers = runif(n) < sdU/(sdL+sdU)
  samples = numeric(n)
  for(i in 1:n){
    if(uppers[i]){
      sample = abs(rnorm(1))
      samples[i] = sdU*sample+mean
    }else{
      sample = -abs(rnorm(1))
      samples[i] = sdL*sample+mean
    }
  }
  if(log) samples = log(samples)

  return(samples)
}





