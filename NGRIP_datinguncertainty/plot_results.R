# Plots results.

plot_results = function(x,
                       plot.proxydata=list(age=TRUE,depth=FALSE,xrev=FALSE,label=NULL),
                       plot.ls = list(fitted=TRUE,legend=NULL,residuals=TRUE,histogram=TRUE,qqplot=TRUE,acf=TRUE,xrev=FALSE,
                                      label.fit=NULL,label.res=NULL,label.hist=NULL,label.qq=NULL,label.acf=NULL),
                       plot.inla.posterior = list(posteriors=TRUE,label=NULL),
                       plot.inlasims = list(nsims=30,legend=NULL,xrev=FALSE,label=NULL),
                       plot.bias = list(MCE=NULL,legend=NULL,xrev=FALSE,label=NULL),
                       plot.linramp = list(depth.reference=NULL,show.t0=TRUE,show.t1=TRUE,xrev=TRUE,label=NULL),
                       plot.event_depth = list(depth.reference=NULL,xrev=TRUE,label=NULL),
                       plot.event_age = list(age.reference=NULL,xrev=TRUE,label=NULL),
                       postscript=FALSE,
                       pdf=FALSE,
                       prefix = "results.plots/figure-",
                       ...){
  if(!postscript && !pdf){
    dev=getOption("device")
    if(!is.character(dev) || dev != "RStudioGD"){
      dev.new(...)
    }
  }else{
    dir = dirname(prefix)
    if (!file.exists(dir) && nchar(dir) > 0L) {
      dir.create(dir, recursive=TRUE)
    } else {
      stopifnot(file.info(dir)$isdir)
    }
  }

  figure.count = 1L
  ageref = x$data$y; z = x$data$z; xx = x$data$x; n=length(ageref)
  reference.label = x$.args$reference.label
  if(is.null(reference.label)){
    reference.label = "reference"
  }
  eventindexes = x$.args$eventindexes; nevents = length(eventindexes)
  oldpar = par()

  if(!is.null(plot.proxydata)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    par(mfrow=c(1,1),mar=(c(5,4.5,4,2)+0.1))
    if(plot.proxydata$age){
      xdata = ageref; xlab=paste0(reference.label," (yr b2k)")
    }else if(plot.proxydata$depth){
      xdata = z; xlab = "Depth (m)"
    }
    xlim = c(xdata[1],xdata[n])
    if(plot.proxydata$xrev) xlim=rev(xlim)
    if(tolower(x$.args$proxy.type) %in% c("d18o","d18","o","oxygen","d","del")){
      plot(xdata,xx,type="l",xlab=xlab,ylab=expression(paste(delta^18,"O (permil)")),main=plot.proxydata$label,xlim=xlim)
    }else if(tolower(x$.args$proxy.type) %in% c("calcium","ca","ca2","ca2+","dust","logcalcium")){
      plot(xdata,xx,type="l",xlab=xlab,ylab=expression(paste("log(",Ca^"2+",") (mL"^"-1",")")),main=plot.proxydata$label,xlim=xlim)
    }

    if(nevents>0) abline(v=xdata[eventindexes],lwd=0.6,col="gray")#rgb(red=0.5,green=0.5,blue=0.5,alpha=1),lwd=0.8)

    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }

  if(tolower(x$.args$transform) %in% "log"){
    dy = x$data$logdy
    ylab = "Log-layer increments per 5 cm"
  }else{
    dy=x$data$dy
    ylab = "Layer increments per 5cm"
  }

  if(!is.null(plot.ls) && !is.null(x$LS.fitting$fit)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
    if(plot.ls$fitted){
      xlim = range(x$data$z)
      if(plot.ls$xrev) xlim=rev(xlim)
      if(is.null(plot.ls$label.fit)){
        plot.label="Least squares fit"
      }else{
        plot.label=plot.ls$label.fit
      }
      plot(x$data$z,dy,type="l",xlab=paste0("Depth (m)"),xlim=xlim,ylab=ylab,main=plot.label)
      lines(x$data$z,x$LS.fitting$fit$fitted.values,col="red")
      abline(v=x$data$z[eventindexes],col="gray",lwd=0.8)
    }
    if(!is.null(plot.ls$legend)) legend(x=leg$x,y=leg$y,legend=leg$legend,col=leg$col,lty=leg$lty,cex=leg$cex,pch=leg$pch,lwd=leg$lwd,pt.cex=leg$pt.cex,bty=leg$bty)

    if(plot.ls$residuals){
      if(is.null(plot.ls$label.res)){
        plot.label="Least squares residuals"
      }else{
        plot.label=plot.ls$label.res
      }
      plot(z,x$LS.fitting$fit$residuals,type="l",xlab="Depth (m)",ylab="Residual errors per 5 cm",main=plot.label,xlim=xlim); abline(h=0,lty=3,col="gray")
    }
    if(plot.ls$histogram){
      if(is.null(plot.ls$label.hist)){
        plot.label="Histogram"
      }else{
        plot.label=plot.ls$label.hist
      }
      hist(x$LS.fitting$fit$residuals,freq=0,col="orange",breaks=20,xlab="Residual errors per 5 cm", main=plot.label)
    }
    if(plot.ls$qqplot){
      if(is.null(plot.ls$label.qq)){
        plot.label="Q-Q Plot"
      }else{
        plot.label=plot.ls$label.qq
      }
      qqnorm(x$LS.fitting$fit$residuals,main=plot.label); qqline(x$LS.fitting$fit$residuals)
    }
    if(plot.ls$acf){
      if(is.null(plot.ls$label.acf)){
        plot.label="Autocorrelation function"
      }else{
        plot.label=plot.ls$label.acf
      }
      acf(x$LS.fitting$fit$residuals,lag.max = 30,main=plot.label,lwd=2)
    }
    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }

  if(!is.null(plot.inla.posterior) && !is.null(x$fitting$fit)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    if(plot.inla.posterior$posteriors){
      if(length(plot.inla.posterior$label)<=1){
        plot.label1 = plot.inla.posterior$label
        plot.label2 = plot.inla.posterior$label
        plot.label3 = plot.inla.posterior$label
      }else{
        plot.label1 = plot.inla.posterior$label[1]
      }
      plot(x$fitting$hyperparameters$posteriors$sigma_epsilon,type="l",xlab=expression(paste(sigma[epsilon])),ylab="Density",lwd=2,main=plot.label1)
      abline(v=x$fitting$hyperparameters$results$sigma_epsilon$mean)
      abline(v=c(x$fitting$hyperparameters$results$sigma_epsilon$quant0.025,x$fitting$hyperparameters$results$sigma_epsilon$quant0.975),col="gray")

      if(tolower(x$.args$noise) %in% c(1,"ar1","ar(1)")){
        if(length(plot.inla.posterior$label)>=2){
          plot.label2 = plot.inla.posterior$label[2]
        }
        plot(x$fitting$hyperparameters$posteriors$phi,xlab=expression(paste(phi)),ylab="Density",lwd=2,type="l",main=plot.label2)
        abline(v=x$fitting$hyperparameters$results$phi$mean)
        abline(v=c(x$fitting$hyperparameters$results$phi$quant0.025,x$fitting$hyperparameters$results$phi$quant0.975),col="gray")
      }else if(tolower(x$.args$noise) %in% c(2,"ar2","ar(2)")){
        if(length(plot.inla.posterior$label)>=3){
          plot.label2 = plot.inla.posterior$label[2]
          plot.label3 = plot.inla.posterior$label[3]
        }
        plot(x$fitting$hyperparameters$posteriors$phi1,xlab=expression(paste(phi[1])),ylab="Density",lwd=2,type="l",main=plot.label2)
        abline(v=x$fitting$hyperparameters$results$phi1$mean)
        abline(v=c(x$fitting$hyperparameters$results$phi1$quant0.025,x$fitting$hyperparameters$results$phi1$quant0.975),col="gray")

        plot(x$fitting$hyperparameters$posteriors$phi2,xlab=expression(paste(phi[2])),ylab="Density",lwd=2,type="l",main=plot.label3)
        abline(v=x$fitting$hyperparameters$results$phi2$mean)
        abline(v=c(x$fitting$hyperparameters$results$phi2$quant0.025,x$fitting$hyperparameters$results$phi2$quant0.975),col="gray")
      }
    }
    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }

  if(!is.null(plot.inlasims) && !is.null(x$simulation)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    nsims = plot.inlasims$nsims
    xlim=range(z)
    if(plot.inlasims$xrev) xlim=rev(xlim)

    ylim=range(x$simulation$summary$lower-ageref,x$simulation$summary$upper-ageref)*1.1
    plot(NA,xlim=xlim,ylim=ylim,xlab="Depth (m)",ylab=paste0("Simulated timescale - ",reference.label," (years)"),main=plot.inlasims$label)
    for(iter in 1:nsims){
      lines(z,x$simulation$age[,iter]-ageref,col="gray",lwd=0.8)
    }
    lines(z,x$simulation$summary$lower-ageref,col="red",lwd=2)
    lines(z,x$simulation$summary$upper-ageref,col="red",lwd=2)
    abline(h=0,lty=3)
    if(x$simulation$summary$.args$CI.type == "hpd"){
      lines(z,x$simulation$summary$median-ageref,col="blue",lwd=2)
    }else{
      lines(z,x$simulation$summary$mean-ageref,col="blue",lwd=2)
    }

    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }

  if(!is.null(plot.bias) && !is.null(x$biases)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    nbiases = x$biases$.args$nbiases

    yrange = c(0,0)
    for( iter in 1:nbiases){
      yrange = range(yrange,x$biases[[paste0("bias",iter)]]$quant0.975-ageref,x$biases[[paste0("bias",iter)]]$quant0.025-ageref)
    }
    if(!is.null(plot.bias$MCE)){
      yrange = range(yrange,-plot.bias$MCE,plot.bias$MCE)
    }
    xlim=range(z)
    if(!is.null(plot.bias$legend)){
      yrange = range(yrange,plot.bias$legend$y)
      xlim=range(xlim,plot.bias$legend$y)
    }

    if(plot.bias$xrev) xlim=rev(xlim)
    plot(NA, xlim=xlim,xlab="Depth (m)", ylab=paste0("Simulated timescale - ",reference.label," (years)"),type="l",col="blue",ylim=yrange,main=plot.bias$label)
    abline(h=0,lty=3,lwd=1,col="gray")
    for(iter in 1:nbiases){
      lines(z,x$biases[[paste0("bias",iter)]]$quant0.975-ageref,col="blue",lty=iter)
      lines(z,x$biases[[paste0("bias",iter)]]$quant0.025-ageref,col="blue",lty=iter)
    }
    if(!is.null(plot.bias$MCE)){
      lines(z,plot.bias$MCE,lwd=2)
      lines(z,-plot.bias$MCE,lwd=2)
    }
    if(!is.null(plot.bias$legend)){
      leg=plot.bias$legend
      legend(leg$x,leg$y,leg$legend,col=leg$col,lty=leg$lty,cex=leg$cex,pch=leg$pch,lwd=leg$lwd,pt.cex=leg$pt.cex)
    }
    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }

  }

  if(!is.null(plot.linramp) && !is.null(x$linramp)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    par(mfrow=c(1,1),mar=(c(5,4.5,4,2)+0.1))
    yrange = range(x$linramp$data$y,x$linramp$linrampfit$q0.025,x$linramp$linrampfit$q0.975)
    xlim = range(x$linramp$data$x)
    if(plot.linramp$xrev) xlim=rev(xlim)
    xval = x$linramp$data$x
    if(!is.null(plot.linramp$label)){
      plot.label=plot.linramp$label
    }else{
      plot.label=x$linramp$.args$label
    }
    if(tolower(x$.args$proxy.type) %in% c("d18o","d18","o","oxygen","d","del")){
    plot(xval,x$linramp$data$y,type="l",lwd=1.25,col="gray",xlim=xlim,ylim=yrange,xlab="Depth (m)",ylab=expression(paste(delta^18, "O (permil)")),main=plot.label)
    }else if(tolower(x$.args$proxy.type) %in% c("calcium","ca","ca2","ca2+","dust","logcalcium")){
      plot(xval,x$linramp$data$y,type="l",lwd=1.25,col="gray",xlim=xlim,ylim=yrange,xlab="Depth (m)",ylab=expression(paste("log(",Ca^"2+",") (mL"^"-1",")")),main=plot.label)
    }
    lines(xval,x$linramp$linrampfit$q0.025,col="red",lwd=2)
    lines(xval,x$linramp$linrampfit$q0.975,col="red",lwd=2)
    lines(xval,x$linramp$linrampfit$mean,col="black",lwd=2)

    if(plot.linramp$show.t0){
      ybottom = min(x$linramp$data$y)
      ytop = min(x$linramp$linrampfit$q0.025)
      margt0=x$linramp$param$t0$marg.t0
      normt0.y = margt0[,2]/diff(range(margt0[,2]))*(ytop-ybottom)-min(margt0[,2])+ybottom
      lines(x=margt0[,1],y=normt0.y,col="blue",lwd=2)
    }
    if(plot.linramp$show.t0 && !is.null(x$linramp$param$t1)){
      ybottom = min(x$linramp$data$y)
      ytop = min(x$linramp$linrampfit$q0.025)
      margt1 = x$linramp$param$t1$marginal
      normt1.y = margt1[,2]/diff(range(margt1[,2]))*(ytop-ybottom)-min(margt1[,2])+ybottom
      lines(x=margt1[,1],y=normt1.y,col="blue",lwd=2,lty=3)
    }

    if(!is.null(plot.event_depth$depth.reference)){
      abline(v=plot.event_age$age.reference,lwd=2,lty=3)
    }else{
      if(!is.null(x$linramp$.args$depth.reference)){
        abline(v= x$linramp$.args$depth.reference,lwd=2,lty=3)
      }
    }
  }

  if(!is.null(plot.event_depth) && !is.null(plot.linramp)){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    par(mfrow=c(1,1),mar=(c(5,4,4,2)+0.1))
    if(is.null(plot.event_depth$label)){
      plot.label = x$linramp$.args$label
    }else{
      plot.label = plot.event_depth$label
    }
    xlim = range(x$linramp$param$t0$marg.t0[,1])
    if(!is.null(plot.event_depth$depth.reference)){
      xlim=range(xlim,plot.event_depth$depth.reference)
    }else{
      if(!is.null(x$linramp$.args$depth.reference)){
        xlim=range(xlim,x$linramp$.args$depth.reference)
      }
    }
    if(plot.event_age$xrev) xlim=rev(xlim)
    plot(x$linramp$param$t0$marg.t0,type="l",lwd=2,xlab="Onset depth (m)",ylab="Density",main=plot.label,xlim=xlim)
    abline(v=x$linramp$param$t0$mean,lwd=2)
    abline(v=c(x$linramp$param$t0$q0.025,x$linramp$param$t0$q0.975),lwd=2,col="gray")

    if(!is.null(plot.event_depth$depth.reference)){
      abline(v=plot.event_age$depth.reference,lwd=2,lty=3)
    }else{
      if(!is.null(x$linramp$.args$depth.reference)){
        abline(v= x$linramp$.args$depth.reference,lwd=2,lty=3)
      }
    }
  }

  if(!is.null(plot.event_age) && !is.null(x$event_dating)){
    xlim = range(x$event_dating$samples)
    if(plot.event_age$xrev) xlim=rev(xlim)
    hist(x$event_dating$samples,xlab="Onset age (y b2k)",freq=FALSE,col="orange",breaks=30,main=x$event_dating$.args$label,xlim=xlim)
    abline(v=x$event_dating$mean,lwd=2)
    abline(v=c(x$event_dating$q0.025,x$event_dating$q0.975),col="gray",lwd=2)
    if(!is.null(plot.event_age$age.reference)){
      abline(v=plot.event_age$age.reference,lwd=2,lty=3)
    }else{
      if(!is.null(x$event_dating$.args$age.reference)){
        abline(v=x$event_dating$.args$age.reference,lwd=2,lty=3)
      }
    }
  }

  return(invisible(x))
}


new.plot = function(postscript,pdf,prefix,figure.count,...)
{

  #dev = getOption("device")
  if(postscript && pdf){
    stop("Multiple file types have been seleced.")
  }
  if(postscript) {
    ending=".eps"
  }else if(pdf){
    ending=".pdf"
  }
  if(!postscript && !pdf){
    dev=getOption("device")
    if(!is.character(dev) || dev != "RStudioGD"){
      dev.new(...)
    }
  }else{

    file.found=FALSE
    while(!file.found){
      filename=paste(prefix,figure.count,ending,sep="")

      if(file.exists(filename)){
        figure.count <- figure.count +1L
      }else{
        file.found=TRUE
      }
    }
    if(postscript){
      postscript(file=filename,...)
    }else if(pdf){
      pdf(file=filename,...)
    }
  }
  return (invisible(figure.count))
}
