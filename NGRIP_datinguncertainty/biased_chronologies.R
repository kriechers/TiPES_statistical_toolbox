## computes chronologies where an additional stochastic bias following a uniform distribution with parameters given by biasparams
biased_chronologies = function(object,bias.model="uniform",biasparams = c(0.99,1.01),nsims=10000,store.samples=FALSE){
  if(nsims > dim(object$simulation$age)[2]) stop("Number of simulated biases exceeds number of simulated chronologies! Shutting down...")
  time.start = Sys.time()
  n = dim(object$simulation$age)[1]
  m = dim(biasparams)[2]
  for(iter in 1:m){
    biasparam = biasparams[,iter]
    biases = runif(nsims,min=biasparam[1],max=biasparam[2])

    if(store.samples){
      biasedages = matrix(NA,nrow=n,ncol=nsims)
    }
    bias.x1=numeric(n)
    bias.x2=numeric(n)
    for(i in 1:nsims){
      sample = biases[i]*object$simulation$age[,i]
      if(store.samples){
        biasedages[,i] = sample
      }
      bias.x1 = bias.x1 + sample
      bias.x2 = bias.x2 + sample^2
    }
    biasmean = bias.x1/nsims
    biassd = sqrt(  1/(nsims-1) * (bias.x2 - 2*bias.x1*biasmean + nsims*biasmean**2)   )

    listr = paste0("bias",iter)
    object$biases[[listr]] = list(mean = biasmean, sd = biassd, quant0.025 = biasmean-1.96*biassd,quant0.975=biasmean+1.96*biassd
    )
    if(store.samples){
      object$biases[[listr]]$simulations = biasedages
    }

  }

  object$biases$.args = list(bias.model=bias.model,biasparam=biasparams,nsims=nsims,store.samples=store.samples,nbiases = m)
  time.end = Sys.time()
  time.full = difftime(time.end,time.start,units="secs")[[1]]
  object$biases$time = time.full
  object$time$biases = time.full
  return(object)
}
