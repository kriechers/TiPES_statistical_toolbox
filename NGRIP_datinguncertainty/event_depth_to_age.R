# Combines linear ramp model fit and Bayesian regression modeling to estimate complete dating uncertainty of the onset of Dansgaard-Oeschger events using Monte Carlo simulation.

event_depth_to_age = function(object, nsims = 10000, print.progress=FALSE, label=NULL,age.reference=NULL){
  time.start=Sys.time()
  n=length(object$data$y)
  if(is.null(label)) label=object$linramp$.args$label
  if(print.progress){
    if(!is.null(label) && label !=""){
      cat("Initiating simulation from complete dating uncertainty of abrupt warming event: ",label,"\n",sep="")
    }else{
      cat("Initiating simulation from complete dating uncertainty of unlabeled abrupt warming event...\n",sep="")
    }
  }


  if(nsims > dim(object$simulation$age)[2]) stop(paste0("Not enough samples from age-depth model (",dim(object$simulation$age)[2],", ",nsims, "needed)"))
  z.sample = inla.rmarginal(nsims,object$linramp$param$t0$marg.t0)

  if(is.null(object$simulation)){
    warning("Cannot find simulations from age-depth model in 'agedepthmodel'!")
  }else{

    time.startsim = Sys.time()
    ysamps=numeric(nsims)
    z = object$data$z
    for(i in 1:nsims){
      if(print.progress==TRUE && i%%1000==0){
        cat("Onset age simulation ",i,"/",nsims,"\n",sep="")
      }
      zgrid = sort(c(z.sample[i],z))
      sampleindex = which(zgrid==z.sample[i] )
      
      #interpolation to find dating uncertainty for in-between depths
      distance = (z.sample[i]-zgrid[sampleindex-1])/(zgrid[sampleindex+1]-zgrid[sampleindex-1])
      ysamps[i] = (1-distance)*object$simulation$age[sampleindex-1,i] +distance*object$simulation$age[sampleindex,i]
    }
  }
  time.endsim = Sys.time()
  if(print.progress) cat("Completed in ", difftime(Sys.time(),time.startsim,units="secs")[[1]]," seconds!\n",sep="")



  dens = density(ysamps); Y0marg= cbind(dens$x,dens$y); colnames(Y0marg) = c("x","y"); Y0marg = inla.smarginal(Y0marg)
  z.Y0 = inla.zmarginal(Y0marg,silent=TRUE)

  object$event_dating = list(samples = ysamps, mean = z.Y0$mean,sd=z.Y0$sd,q0.025=z.Y0$quant0.025,q0.5=z.Y0$quant0.5,q0.975=z.Y0$quant0.975)

  object$event_dating$.args = list(nsims = nsims,label=label,age.reference=age.reference)
  object$time$event_age = list(simulation = difftime(time.endsim,time.startsim,units="secs")[[1]],total = difftime(Sys.time(),time.start,units="secs")[[1]])

  return(object)
}


