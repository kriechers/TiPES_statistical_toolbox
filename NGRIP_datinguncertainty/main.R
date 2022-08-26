
###################################################################################
################################## Control panel ##################################
###################################################################################
#Set proper directory
setwd("/Users/emy016/Dropbox/Postdoc/Submission files/Uncertaintypaper/NGRIP_datinguncertainty/")

nsims=200 #number of simulated chronologies, make sure you have enough memory

do.sync=TRUE #should chronologies be synchronized? Currently only supports noise="ar1",
#do.dust=FALSE.

tiepoint_pdf = "adolphi" #if synchronized, which distribution is assumed on tie-points.
# Possible choices ("adolphi", "skewgaussian", "gaussian", "fixed")
# The Adolphi distribution is taken from Adolphi et al (2018) and Muscheler et al (2020)

do.dust=FALSE #should dust proxy be used?

do.bias = FALSE #should unknown bias also be included (not implemented for do.sync=TRUE)

compare.noise = FALSE #should AR(1) model be compared to iid and AR(2) models?

do.transitions = TRUE #should dating of abrupt warming transitions be included?

eventnumber=13 #If you want a single transition to be computed (between 1 and 29)

print.progress=TRUE #Should progression reports of each function be printed to screen?

### To compute everything, run: source("main.R") in console

###################################################################################
###################################################################################
###################################################################################

#load functions and packages
source("biased_chronologies.R"); source("chronology_simulation.R"); source("event_depth_to_age.R")
source("helpfulfunctions.R"); source("linrampfitter.R"); source("modelfitter.R")
source("plot_results.R"); source("prepare.R"); source("rgeneric_model.R")
source("summary_results.R"); source("mainfunctions.R")
library("INLA"); library("matrixStats"); library("numDeriv")

library(readODS)
library(readxl)

#read data
maindata = read_excel("datasets_used/NGRIP_d18O_and_dust_5cm.xls",
                      sheet="NGRIP-2 d18O and Dust",col_types="numeric" )

####################################################
##### Fit model and sample chronology ensemble #####
####################################################


startindex = 2921 #11703yb2k (holocene onset)

depth = maindata$`NGRIP-2 depth (m)`[startindex:nrow(maindata)]
age = maindata$`GICC05 age (yr b2k)`[startindex:nrow(maindata)]
MCE = maindata$`GICC05 MCE (yr)`[startindex:nrow(maindata)]
water = maindata$`Delta O18 (permil)`[startindex:nrow(maindata)]#dust, NA until 1167
dust0 = maindata$`Dust count (ml^-1)`[startindex:nrow(maindata)]#dust, NA until 1167

library(zoo)
dust = na.approx(dust0) #Impute missing data using linear interpolation

dust = log(dust)

do.dust = FALSE #set FALSE for d18O proxy, and true for Ca2+ proxy

if(do.dust){
  proxy.type = "ca"
  proxy=dust

}else{
  proxy.type = "d18O"
  proxy=water
}


#load abrupt warming transitions
eventdata = read_ods("datasets_used/GISevents.ods")
GISevents = eventdata[(eventdata$`NGRIP depth (m)`>min(depth))&(eventdata$`NGRIP depth (m)`<max(depth)),]
event_intervals = read.table("datasets_used/event_intervals.txt")


eventdepths = GISevents$`NGRIP depth (m)`
eventindexes = c(1,which.index(eventdepths, depth[2:length(depth)]) )


#plot d18O proxies as a function of depth and age (GICC05), respectively
par(mfrow=c(1,1),mar=c(5,4.5,2,2)+0.1)
plot(depth,water,type="l",xlab="Depth (m)",ylab=expression(paste(delta^18,"O (permil)")),xlim=rev(range(depth))); abline(v=eventdepths,lwd=0.7,col="gray")
plot(age,water,type="l",xlab="Age (yb2k)",ylab=expression(paste(delta^18,"O (permil)")),xlim=rev(range(age))); abline(v=age[eventindexes],lwd=0.7,col="gray")




#events used in regression model
events=eventdepths #locations of transitions used in regression model. Pairs with 'eventmeasure' for finding the corresponding indices
eventmeasure="depth"

#reg.model specifies the structure for which the formula string should be created
# const for intercept, depth1 for linear response to depth, depth2 for quadratic response
# to depth, proxy for linear response to proxy, psi0 for intercept for each climate period
# psi1 for linear trend for each climate period
reg.model = list(
  const=FALSE,depth1=FALSE,depth2=TRUE,proxy=TRUE,psi0=TRUE,psi1=TRUE); method="inla";





#load data window and specifics to transition
lowerints = which.index(event_intervals$depth_int_lower.m, depth[2:length(depth)])
upperints = which.index(event_intervals$depth_int_upper.m, depth[2:length(depth)])
depth.reference = event_intervals$NGRIP_depth_m[eventnumber]
age.reference = event_intervals$GICC_age.yb2k[eventnumber]
interval = lowerints[eventnumber]:upperints[eventnumber]



#list object specifying parameters for abrupt warming transition dating
event.estimation = list(interval=interval,t1.sims=50000,rampsims=50000,label="GI-11",
                     depth.reference=event_intervals$NGRIP_depth_m[eventnumber],
                     age.reference=event_intervals$NGRIP_agee_m[eventnumber])

reference.label="GICC05" #this string is used in plot_results




transform = "identity" #set equal to "log" for logarithmic transformation
noise = "ar1" #iid, ar1 and ar2 supported

if(do.bias){
  bias = list(bias.model="uniform",biasparams=cbind( c(1,1),c(0.98,1.02),c(0.96,1.04) ),
              store.samples=FALSE)
}else{
  bias = NULL
}

store.everything=FALSE #should fixed model components and stochastic both be stored (TRUE), or just the sum (FALSE)

#run main function wrapper
object = main(age,depth,proxy, events=eventdepths,nsims=nsims, eventmeasure = eventmeasure,proxy.type=proxy.type,
              reference.label=reference.label,transform=transform,reg.model = reg.model,
              noise=noise,
  event.estimation = event.estimation,store.everything=store.everything,
  print.progress=print.progress,bias = bias
  )

summary_results(object) #print out summary statistics

plot_results(object) #plot results



#############################
###### Synchronization ######
#############################



plotcenter = "median" #can be mean, mode or median. For symmetric tie-points these are all equal

# set tie-point distribution. Possible choices ("adolphi", "skewgaussian", "gaussian", "fixed")
# The Adolphi distribution is taken from Adolphi et al (2018) and Muscheler et al (2020)
tiepoint_pdf = tiepoint_pdf


z=object$data$z; y=object$data$y; n=length(z)

#Set location of tie-points in GICC05 age
tie_gicc05 = c(11050,12050,13050,22050,42050)

# We omit first Adolphi point since it takes place after the Holocene (hence it is not considered)
tie_used = 2:5

tie_n = length(tie_gicc05)

# Find indexes corresponding to the tie-point datings
tie_indexes = numeric(tie_n)
for(i in 1:tie_n){
  #tie_indexes[i] = which( abs( tie_depths[i]-z ) == min(abs(  tie_depths[i]-z  )) )
  tie_indexes[i] = which( abs( tie_gicc05[i]-y ) == min(abs(  tie_gicc05[i]-y  )) )
}
#tie_depths = c(1455.80,1503.25,1533.15,1767.70,2132.85)[1:5]

tie_depths = (z[tie_indexes]) #Tie-point location in NGRIP depth
free_indexes = (1:n)[-tie_indexes] #which indexes are not associated with tie-points
free_n = length(free_indexes) #number of free variables

plotdens=TRUE #should tie-point distributions be plotted?

## The code block below computes the tie-point distributions and generates samples from them
if(tolower(tiepoint_pdf) == "adolphi"){
  tieshifts = c(11050,12050,13050,22050,42050)
  tiepointsims = adolphi_tiepoint_simmer(nsims=nsims,tieshifts = tieshifts,x.ref=tieshifts,plotdens=plotdens,
                                         plothist=list(compute=TRUE,breaks=50,col="orange"))
  #tiepointsims = tiepointsims[tie_used,]
  tie_marginals = adolphi_tiepoint_pdfs()
  tie_summary = numeric(0)
  for(i in 1:5){
    zmarg = inla.zmarginal(tie_marginals[[i]],silent=TRUE)
    q1 = inla.qmarginal(0.1586553,tie_marginals[[i]])
    q2 = inla.qmarginal(0.8413447,tie_marginals[[i]])
    tiemode=((tie_marginals[[i]])[,1])[which( (tie_marginals[[i]])[,2] == max((tie_marginals[[i]])[,2]))]
    tie_summary = rbind(tie_summary,c(tie_depths[i],tie_gicc05[i],zmarg$mean,zmarg$sd,tiemode,zmarg$quant0.025,q1,zmarg$quant0.5,q2,zmarg$quant0.975))
  }
  colnames(tie_summary) = c("NGRIP depth","GICC05 age","mean","sd","mode","q0.025","q0.158","median","q0.841","q0.975")
  rownames(tie_summary) = 1:tie_n
  tie_summary

}else if(tolower(tiepoint_pdf) == "fixed"){
  #tie_offsets = c(-60.1,-6.42,2.64,532.36,243.57)#[tie_used]
  tie_offsets = c(-60.0,-7.0,3.3,521.4,275.2)#[tie_used]
  tie_age = y[tie_indexes]+tie_offsets
  tiepointsims = t(matrix(rep(tie_age,each=nsims),ncol=tie_n))
  tie_summary = numeric(0)
  for(i in 1:tie_n){
    tie_summary = rbind(tie_summary,c(tie_depths[i],tie_age[i],tie_age[i]))
  }
  colnames(tie_summary) = c("NGRIP depth","GICC05 age","offset")
  rownames(tie_summary) = 1:tie_n
}else if(tolower(tiepoint_pdf) == "skewgaussian"){
  tie_offsets = c(-60,-6,1,538,244)[tie_used]
  tie_age = y[tie_indexes]+tie_offsets
  tie_sigmaL = c(7.79,10.33,6.25,145.25,321.99)[tie_used]
  tie_sigmaU = c(9,   9,    10,  63.26, 213.5)[tie_used]
  tiepointsims = matrix(NA,nrow=tie_n,ncol=nsims)
  for(i in 1:tie_n){
    tiepointsims[i,] = skewsampler(nsims,mean=tie_age[i],sdL=tie_sigmaL[i],sdU=tie_sigmaU[i])
  }
  tie_summary = numeric(0)
  for(i in 1:tie_n){
    tie_marginal = density(tiepointsims[i,])
    zmarg = inla.zmarginal(tie_marginal,silent=TRUE)
    tiemode=((tie_marginal)[,1])[which( (tie_marginal)[,2] == max((tie_marginal)[,2]))]

    ###fill in theoretical values
    tie_summary = rbind(tie_summary,c(tie_depths[i],tie_age[i],zmarg$mean,zmarg$sd,tiemode,zmarg$quant0.025,zmarg$quant0.5,zmarg$quant0.975))
  }
  colnames(tie_summary) = c("NGRIP depth","GICC05 age","mean","sd","mode","q0.025","median","q0.975")
  rownames(tie_summary) = 1:tie_n
  if(plotdens){
    if(!is.null(plothist) && plothist$compute==TRUE){
      la=layout(mat=matrix(c(1,4,1,4,2,4,2,5,3,5,3,5) ,nrow=2))
      layout.show(la)

      plot(density(tiepointsims[1,]-y[tie_indexes[1]]),type="l",lwd=2,xlab="Age (yb2k)",ylab="Density",main="(a) Tie-point 1")
      abline(v=tie_offsets[1]); abline(v=c(tie_offsets[1]+2*tie_sigmaU[1],tie_offsets[1]-2*tie_sigmaL[1]),col="gray")
      plot(density(tiepointsims[2,]-y[tie_indexes[2]]),type="l",lwd=2,xlab="Age (yb2k)",ylab="Density",main="(b) Tie-point 2")
      abline(v=tie_offsets[2]); abline(v=c(tie_offsets[2]+2*tie_sigmaU[2],tie_offsets[2]-2*tie_sigmaL[2]),col="gray")
      plot(density(tiepointsims[3,]-y[tie_indexes[3]]),type="l",lwd=2,xlab="Age (yb2k)",ylab="Density",main="(c) Tie-point 3")
      abline(v=tie_offsets[3]); abline(v=c(tie_offsets[3]+2*tie_sigmaU[3],tie_offsets[3]-2*tie_sigmaL[3]),col="gray")
      plot(density(tiepointsims[4,]-y[tie_indexes[4]]),type="l",lwd=2,xlab="Age (yb2k)",ylab="Density",main="(d) Tie-point 4")
      abline(v=tie_offsets[4]); abline(v=c(tie_offsets[4]+2*tie_sigmaU[4],tie_offsets[4]-2*tie_sigmaL[4]),col="gray")

      par(mfrow=c(1,1))

    }
  }
}else if(tolower(tiepoint_pdf) == "gauss"){
  tie_sigma = c(8.42,9.69,8.34,112.03,273.68)[2:5]

}else{
  stop("Could not recognize tie-point uncertainty distribution")
}



# This code makes sure that the tie-points and free variables are correctly represented
# in case some tie-points have been chosen to be omitted earlier in the code
tie_depths_used = tie_depths[tie_used]
tie_indexes_used = tie_indexes[tie_used]
free_indexes_used = (1:n)[-tie_indexes_used]
free_n_used = length(free_indexes_used)
tie_n_used = length(tie_indexes_used)
tie_marginals_used = tie_marginals[tie_used]
free_depths_used = object$data$z[free_indexes_used]

tie_gicc05_used = tie_gicc05[tie_used]
tiepointsims_used = tiepointsims[tie_used,]



## This code executes the syncsampler function and generates an ensemble of plausible
## synchronized chronologies.
object_sync = syncsampler(object,nsims=nsims,tie_depths=tie_depths_used,
                     tiepointsims=tiepointsims_used,
                     free_indexes=free_indexes_used,
                     tie_indexes=tie_indexes_used)



# ##### fixed tie-points #####
#
#
# syncsampler(object,
#             tiepoints = list(depths=tie_depths,
#                              samples=c(rep(tie_age[1],10000),
#                                        rep(tie_age[2],10000),
#                                        rep(tie_age[3],10000),
#                                        rep(tie_age[4],10000),
#                                        rep(tie_age[5],10000)),
#                              template=NULL))
#
#



## Summarize the synchronized chronology ensemble ##


tie_n_used = length(tie_used)
free_n=free_n_used;tie_n=tie_n_used
free_indexes=free_indexes_used;tie_indexes=tie_indexes_used

emp.modes = numeric(free_n)
emp.hpdupper = numeric(free_n)
emp.hpdlower = numeric(free_n)
emp.medians = numeric(free_n)

if(tiepoint_pdf != "fixed"){
  tie_means = numeric(tie_n_used)
  tie_sds = numeric(tie_n_used)
  tie_medians = numeric(tie_n_used)
  tie_hpdlower= numeric(tie_n_used)
  tie_hpdupper= numeric(tie_n_used)
  for(i in 1:tie_n_used){
    dens=tie_marginals_used[[i]]
    dens = data.frame(x=dens[,1],y=dens[,2])
    zm=inla.zmarginal(dens,silent=TRUE)
    tie_means[i] = zm$mean
    tie_sds[i] = zm$sd
    tie_medians[i]=zm$quant0.5
    hpd = inla.hpdmarginal(0.95,dens)
    tie_hpdlower[i] = hpd[1]
    tie_hpdupper[i] = hpd[2]
  }
}


my.var = function(vek,mode){
  return(  1/(length(vek)-1)*sum( (vek-mode)^2 )  )
}

for( i in 1:free_n){
  dens = density(object_sync$simulation$age[free_indexes[i],])
  if(plotcenter == "mean"){
    emp.modes[i] = mean(object_sync$simulation$age[free_indexes[i],]) #mean
  }else if(plotcenter == "median"){
    emp.modes[i] = median(object_sync$simulation$age[free_indexes[i],]) #median
  }else{
    emp.modes[i]=dens$x[which(dens$y == max(dens$y))] #mode
  }

  dens0 = data.frame(x=dens$x,y=dens$y)
  zm=inla.zmarginal(dens0,silent = TRUE)
  emp.medians[i] =zm$quant0.5

  emp.hpdlower[i] = inla.hpdmarginal(0.95,dens)[1]
  emp.hpdupper[i] = inla.hpdmarginal(0.95,dens)[2]
}
object_sync$simulation$summary = list(mean=rep(NA,n),sd=rep(NA,n),lower=rep(NA,n),upper=rep(NA,n))
object_sync$simulation$summary$mean[free_indexes] = rowMeans(object_sync$simulation$age[free_indexes,])
object_sync$simulation$summary$sd[free_indexes] = rowSds(object_sync$simulation$age[free_indexes,])
object_sync$simulation$summary$median[free_indexes] = rowMedians(object_sync$simulation$age[free_indexes,])
object_sync$simulation$summary$lower[free_indexes] = emp.hpdlower
object_sync$simulation$summary$upper[free_indexes] = emp.hpdupper

object_sync$simulation$summary$mean[tie_indexes] = tie_means+tie_gicc05_used
object_sync$simulation$summary$sd[tie_indexes] = tie_sds+tie_gicc05_used
object_sync$simulation$summary$median[tie_indexes] = tie_medians+tie_gicc05_used
object_sync$simulation$summary$lower[tie_indexes] = tie_hpdlower+tie_gicc05_used
object_sync$simulation$summary$upper[tie_indexes] = tie_hpdupper+tie_gicc05_used

object_sync$simulation$summary$.args$CI.type = "hpd"
object_sync$simulation$summary$.args$interval = cbind(emp.hpdlower,emp.hpdupper)
colnames(object_sync$simulation$summary$.args$interval)=c("lower","upper")

object_sync$reference.age="GICC05"

## Visualize synchronized chronology

## Plot synchronized time scale

library(ggplot2)
gicc = object_sync$data$y[free_indexes_used]
fullpd = data.frame(depth = object_sync$data$z[free_indexes_used], medians=emp.medians-gicc,
                    lower=emp.hpdlower-gicc,upper=emp.hpdupper-gicc)
gg2 = ggplot(data=fullpd,aes(x=depth)) +
  geom_line(aes(y=0),color="blue",linetype="dotted",size=0.2)+
  geom_line(aes(y=medians))+
  geom_ribbon(aes(ymin=lower,ymax=upper),color="red",fill="red",alpha=0.3)+
  theme_bw()+
  xlab("NGRIP depth (m)")+ ylab("Estimated age - GICC05 (years)")

if(tiepoint_pdf != "fixed"){
  gg2 = gg2 + geom_segment(data=data.frame(depth=tie_depths_used,median=tie_medians,
                                           lower=tie_hpdlower,upper=tie_hpdupper),
                           aes(x=tie_depths_used,y=lower,xend=tie_depths_used,yend=upper),
                           col="magenta")
}else{
  gg2 = gg2 + geom_point(data=data.frame(tiedepths = tie_depths_used,
                                         tiemid = tie_offsets[tie_used]),
                         aes(x=tiedepths,y=tiemid),
                         color="blue")
}

print(gg2)




##############################################
### Comparing models: iid, AR(1) and AR(2) ###
##############################################

transform = "identity" #set equal to "log" for logarithmic transformation, "identity" for no transformation
do.dust = TRUE #TRUE for Ca2+ proxy, FALSE for d18O proxy
if(do.dust){
  proxy.type = "ca"
  proxy=dust

}else{
  proxy.type = "d18O"
  proxy=water
}


#perform analysis on iid, ar1 and ar2 models
objectiid= main(age,depth,proxy,eventdepths,transform=transform,noise="iid")
objectar1= main(age,depth,proxy,eventdepths,transform=transform,noise="ar1")
objectar2= main(age,depth,proxy,eventdepths,transform=transform,noise="ar2")


#summary results
summary_results(objectiid); summary_results(objectar1); summary_results(objectar2)


## plot uncertainties and compare
xlim = rev(range(objectiid$data$z))
gicc05=objectiid$data$y
ylim=range(objectar2$simulation$summary$lower-gicc05,objectar2$simulation$summary$upper-gicc05)

plot(x=objectiid$data$z,y=objectiid$simulation$summary$mean-gicc05,xlim=xlim,ylim=ylim,xlab="Depth (m)",ylab="Simulated time scale - GICC05 (years)",
     type="l",col="gray",lwd=1.5)
lines(x=objectiid$data$z,y=objectiid$simulation$summary$lower-gicc05,col="black",lwd=1.5)
lines(x=objectiid$data$z,y=objectiid$simulation$summary$upper-gicc05,col="black",lwd=1.5)
lines(x=objectiid$data$z,y=objectar1$simulation$summary$lower-gicc05,col="blue",lwd=1.5)
lines(x=objectiid$data$z,y=objectar1$simulation$summary$upper-gicc05,col="blue",lwd=1.5)
lines(x=objectiid$data$z,y=objectar2$simulation$summary$lower-gicc05,col="red",lwd=1.5)
lines(x=objectiid$data$z,y=objectar2$simulation$summary$upper-gicc05,col="red",lwd=1.5)
abline(h=0,lty=3,col="blue",lwd=0.8)


if(transform == "log"){
  legend(1600,310,legend=c("Mean","iid CI", "AR(1) CI", "AR(2) CI"),
         col=c("gray","black","blue","red"),
         lty=c(1,1), cex=0.5)
}else{
  legend(1570,250,legend=c("Mean","iid CI", "AR(1) CI", "AR(2) CI"),
         col=c("gray","black","blue","red"),
         lty=c(1,1), cex=0.4)
}



#plot posterior marginal distributions of model parameters for the different models
layout(mat=matrix(c(1,2,4,7,3,5,8,9,6),nrow=3))
par(mar=c(4,4,2,1))
xrange_sigma = c(0.420,0.442)
xrange_sigma = range(objectiid$fitting$hyperparameters$posteriors$sigma_epsilon[,1],
                     objectar1$fitting$hyperparameters$posteriors$sigma_epsilon[,1],
                     objectar2$fitting$hyperparameters$posteriors$sigma_epsilon[,1])
xrange_phi = range(objectar1$fitting$hyperparameters$posteriors$phi[,1],
                   objectar2$fitting$hyperparameters$posteriors$phi1[,1],
                   objectar2$fitting$hyperparameters$posteriors$phi2[,1])

plot(objectiid$fitting$hyperparameters$posteriors$sigma_epsilon,xlab=expression(paste(sigma[epsilon])),ylab="Density",lwd=2,type="l",main="(a) Independent residuals",xlim=xrange_sigma)
abline(v=objectiid$fitting$hyperparameters$results$sigma_epsilon$mean);abline(v=c(objectiid$fitting$hyperparameters$results$sigma_epsilon$quant0.025,objectiid$fitting$hyperparameters$results$sigma_epsilon$quant0.975),col="gray")

plot(objectar1$fitting$hyperparameters$posteriors$sigma_epsilon,xlab=expression(paste(sigma[epsilon])),ylab="Density",lwd=2,type="l",main="(b) AR(1) residuals",xlim=xrange_sigma)
abline(v=objectar1$fitting$hyperparameters$results$sigma_epsilon$mean);abline(v=c(objectar1$fitting$hyperparameters$results$sigma_epsilon$quant0.025,objectar1$fitting$hyperparameters$results$sigma_epsilon$quant0.975),col="gray")
plot(objectar1$fitting$hyperparameters$posteriors$phi,xlab=expression(phi),ylab="Density",lwd=2,type="l",main="(c) AR(1) residuals",xlim=xrange_phi)
abline(v=objectar1$fitting$hyperparameters$results$phi$mean);abline(v=c(objectar1$fitting$hyperparameters$results$phi$quant0.025,objectar1$fitting$hyperparameters$results$phi$quant0.975),col="gray")

plot(objectar2$fitting$hyperparameters$posteriors$sigma_epsilon,xlab=expression(paste(sigma[epsilon])),ylab="Density",lwd=2,type="l",main="(d) AR(2) residuals",xlim=xrange_sigma)
abline(v=objectar2$fitting$hyperparameters$results$sigma_epsilon$mean);abline(v=c(objectar2$fitting$hyperparameters$results$sigma_epsilon$quant0.025,objectar2$fitting$hyperparameters$results$sigma_epsilon$quant0.975),col="gray")
plot(objectar2$fitting$hyperparameters$posteriors$phi1,xlab=expression(phi[1]),ylab="Density",lwd=2,type="l",main="(e) AR(2) residuals",xlim=xrange_phi)
abline(v=objectar2$fitting$hyperparameters$results$phi1$mean);abline(v=c(objectar2$fitting$hyperparameters$results$phi1$quant0.025,objectar2$fitting$hyperparameters$results$phi1$quant0.975),col="gray")
plot(objectar2$fitting$hyperparameters$posteriors$phi2,xlab=expression(phi[2]),ylab="Density",lwd=2,type="l",main="(f) AR(2) residuals",xlim=xrange_phi)
abline(v=objectar2$fitting$hyperparameters$results$phi2$mean);abline(v=c(objectar2$fitting$hyperparameters$results$phi2$quant0.025,objectar2$fitting$hyperparameters$results$phi2$quant0.975),col="gray")

par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)





#############################
### Abrupt warming events ###
#############################



transform = "identity" ##"log" for logaritmic transformation, "identity" for ordinary
#nsims = 10000

eventdata = read_ods("datasets_used/GISevents.ods")
event_intervals = read.table("datasets_used/event_intervals.txt")


if(!do.sync){
  #d18O
  depth = maindata$`NGRIP-2 depth (m)`[2979:nrow(maindata)]
  proxy = maindata$`Delta O18 (permil)`[2979:nrow(maindata)]#dust, NA until 1167
  age = maindata$`GICC05 age (yr b2k)`[2979:nrow(maindata)]
  MCE = maindata$`GICC05 MCE (yr)`[2979:nrow(maindata)]

  GISevents = eventdata[(eventdata$`NGRIP depth (m)`>min(depth))&(eventdata$`NGRIP depth (m)`<max(depth)),]

  eventdepths = GISevents$`NGRIP depth (m)`
  eventindexes = c(1,which.index(eventdepths, depth[2:length(depth)]) )

  plot(depth,proxy,type="l",xlab="Depth",ylab="d18O"); abline(v=eventdepths,lwd=0.7,col="gray")

  reg.model = list( const=FALSE,depth1=FALSE,depth2=TRUE,proxy=TRUE,psi0=TRUE,psi1=TRUE)
  object_d18O= main(age,depth,proxy,eventdepths,nsims=nsims,transform=transform,
                    noise="ar1",reg.model=reg.model,proxy.type="d18O")
}else{
  object_d18O = object_sync
  do.dust=FALSE #dust is currently disabled for synchronized chronologies
}
#dust
if(do.dust){

  depth = maindata$`NGRIP-2 depth (m)`[1167:nrow(maindata)]
  proxy0 = maindata$`Dust count (ml^-1)`[1167:nrow(maindata)]#dust, NA until 1167
  age = maindata$`GICC05 age (yr b2k)`[1167:nrow(maindata)]
  MCE = maindata$`GICC05 MCE (yr)`[1167:nrow(maindata)]
  library(zoo)
  proxy = log(na.approx(proxy0))

  GISevents = eventdata[(eventdata$`NGRIP depth (m)`>min(depth))&(eventdata$`NGRIP depth (m)`<max(depth)),]

  eventdepths = GISevents$`NGRIP depth (m)`
  eventindexes = c(1,which.index(eventdepths, depth[2:length(depth)]) )

  plot(depth,proxy,type="l",xlab="Depth",ylab="log(Ca2+)"); abline(v=eventdepths,lwd=0.7,col="gray")

  reg.model = list(const=FALSE,depth1=FALSE,depth2=TRUE,proxy=TRUE,psi0=TRUE,psi1=TRUE)


  object_dust= main(age,depth,proxy,eventdepths,nsims=nsims,transform=transform,
                    noise="ar1",reg.model=reg.model,proxy.type="calcium")


}



  #allocate data

  if(do.dust){

    agesimmatrix_dust = matrix(NA,29,nsims)
    depthlist_dust = c()
    depthstats = as.data.frame(matrix(NA,7,29))
    colnames(depthstats) = c("true", "mean", "CIL","CIU","dustmean","dustCIL","dustCIU")
    agestats = as.data.frame(matrix(NA,7,29))
    colnames(agestats) = c("true", "mean", "CIL","CIU","dustmean","dustCIL","dustCIU")

  }else{
    depthstats = as.data.frame(matrix(NA,4,29))
    colnames(depthstats) = c("true", "mean", "CIL","CIU")
    agestats = as.data.frame(matrix(NA,4,29))
    colnames(agestats) = c("true", "mean", "CIL","CIU")

  }





do.plot.rampfit = TRUE #if do.transitions, should each linear ramp model fit be plotted?
do.plot.agehist = TRUE #if do.transitions, should each histogram of onset age be plotted?


agesimmatrix_d18O = matrix(NA,29,nsims)

depthlist_d18O = c()


steplengths = rep(0.01,29)
steplengths[21]=0.001 #sometimes changing the steplength slightly might improve convergence (default is 0.01)

#iterate over all 29 transitions
for(eventnumber in 1:29){
  #find data windows
  lowerints = which.index(event_intervals$depth_int_lower.m, object_d18O$data$z)
  upperints = which.index(event_intervals$depth_int_upper.m, object_d18O$data$z)
  interval_d18O = lowerints[eventnumber]:upperints[eventnumber]




  depth.reference = event_intervals$NGRIP_depth_m[eventnumber]
  age.reference = event_intervals$GICC_age.yb2k[eventnumber]

  #specify parameters for transition onset age estimation
  event.estimation = list(interval_dust=interval_dust,interval_d18O=interval_d18O,
                          t1.sims=50000,rampsims=50000,h=steplengths[eventnumber],
                          label=paste0(eventnumber," (",event_intervals[eventnumber,1],"): ",event_intervals[eventnumber,2]),
                          depth.reference=event_intervals$NGRIP_depth_m[eventnumber],
                          age.reference=event_intervals$NGRIP_agee_m[eventnumber])

  #fit linear ramp


  if(do.plot.rampfit && do.dust){ #plot linear ramp fit
    plot_results(eventobject_dust,plot.proxydata=NULL,plot.ls = NULL,
                 plot.inla.posterior = NULL,plot.inlasims = NULL,plot.bias = NULL,
                 plot.linramp = list(depth.reference=NULL,show.t0=TRUE,show.t1=TRUE,xrev=TRUE,
                                     paste0("Dust: ",event.estimation$label)),
                 plot.event_depth = NULL,plot.event_age = NULL)
  }
  eventobject_d18O = linrampfitter(object_d18O,interval=event.estimation$interval_d18O,
                                h=event.estimation$h,t1.sims=event.estimation$t1.sims,
                                rampsims=event.estimation$rampsims,
                                label=paste0("d18O: ",event.estimation$label),
                                log.ramp=FALSE,
                                depth.reference = event.estimation$depth.reference,
                                print.progress=print.progress)
  if(do.plot.rampfit){
    plot_results(eventobject_d18O,plot.proxydata=NULL,plot.ls = NULL,
                 plot.inla.posterior = NULL,plot.inlasims = NULL,plot.bias = NULL,
                 plot.linramp = list(depth.reference=NULL,show.t0=TRUE,show.t1=TRUE,xrev=TRUE,
                                     paste0("d18O: ",event.estimation$label)),
                 plot.event_depth = NULL,plot.event_age = NULL)
  }

  #sample onset age

  agesimmatrix_dust[eventnumber,] = eventobject_dust$event_dating$samples
  if(do.plot.agehist && do.dust){ #plot histogram of simulated onset ages
    hist(agesimmatrix_dust[eventnumber,],freq=0,breaks=50,col="orange",main=paste0("Dust: ",event.estimation$label))
  }
  eventobject_d18O = event_depth_to_age(eventobject_d18O, nsims = nsims, print.progress=print.progress,
                                  label=event.estimation$label,
                                  age.reference = event.estimation$age.reference)
  agesimmatrix_d18O[eventnumber,] = eventobject_d18O$event_dating$samples
  if(do.plot.agehist){
    hist(agesimmatrix_d18O[eventnumber,],freq=0,breaks=50,col="orange",
         main=paste0("d18O: ",event.estimation$label),
         xlim=range(agesimmatrix_d18O[eventnumber,],age.reference))
    abline(v=age.reference,lwd=3,lty=2)
  }

  #compute summary statistics
  depthstats[1,eventnumber] = event_intervals$NGRIP_depth_m[eventnumber]
  depthstats[2,eventnumber] = eventobject_d18O$linramp$param$t0$mean
  depthstats[3,eventnumber] = eventobject_d18O$linramp$param$t0$q0.025
  depthstats[4,eventnumber] = eventobject_d18O$linramp$param$t0$q0.975


  zage1 = inla.zmarginal(densage1,silent=TRUE)


  agestats[1,eventnumber] = event_intervals$GICC_age.yb2k[eventnumber]
  agestats[2,eventnumber] = mean(agesimmatrix_d18O[eventnumber,])
  agestats[3,eventnumber] = zage1$quant0.025
  agestats[4,eventnumber] = zage1$quant0.975

  if(do.dust){
    lowerints = which.index(event_intervals$depth_int_lower.m, object_dust$data$z)
    upperints = which.index(event_intervals$depth_int_upper.m, object_dust$data$z)
    interval_dust = lowerints[eventnumber]:upperints[eventnumber]

    eventobject_dust = linrampfitter(object_dust,interval=event.estimation$interval_dust,
                                     h=event.estimation$h,t1.sims=event.estimation$t1.sims,
                                     rampsims=event.estimation$rampsims,
                                     label=paste0("Dust: ",event.estimation$label),
                                     depth.reference = event.estimation$depth.reference,
                                     print.progress=print.progress)

    eventobject_dust = event_depth_to_age(eventobject_dust, nsims = nsims, print.progress=print.progress,
                                          label=event.estimation$label,
                                          age.reference = event.estimation$age.reference)


    depthstats[5,eventnumber] = eventobject_d18O$linramp$param$t0$mean
    depthstats[6,eventnumber] = eventobject_d18O$linramp$param$t0$q0.025
    depthstats[7,eventnumber] = eventobject_d18O$linramp$param$t0$q0.975
    densage1 = density(agesimmatrix_d18O[eventnumber,]); densage1 = data.frame(x=densage1$x,y=densage1$y)
    densage2 = density(agesimmatrix_dust[eventnumber,]); densage2 = data.frame(x=densage2$x,y=densage2$y)
    zage2 = inla.zmarginal(densage2,silent=TRUE)
    agestats[5,eventnumber] = mean(agesimmatrix_dust[eventnumber])
    agestats[6,eventnumber] = zage1$quant0.025
    agestats[7,eventnumber] = zage1$quant0.975


  }
}

#import onsets from other sources to compare with: Buizert et al. (2015) and Capron et al. (2021)
rasmussen_ages_full = read_excel("datasets_used/Rasmussen_et_al_2014_QSR_Table_2.xlsx",
                                 skip=21)

rasmussen_ages = event_intervals$GICC_age.yb2k #GISevents$`Age (a b2k)`

buizert_onsets = read_excel("datasets_used/Buizert_onsets.xlsx",
                            col_types=c("text","numeric","numeric","numeric","numeric"))
buizert_depths = buizert_onsets$Depth; buizert_ages = buizert_onsets$Age
#buizert_depths[c(2:6,8:9,12:13,15:17,19:21)]
buizert_ages = buizert_ages[c(2:6,8:9,12:13,15:17,19:21)]

capron_onsets_NGRIP_d18O = read_excel("datasets_used/Capron_onsets.xls",
                           sheet="d18O",n_max=25)
capron_ages_NGRIP_d18O = capron_onsets_NGRIP_d18O$`t1 (50%)`

capron_onsets_NGRIP_dust = read_excel("datasets_used/Capron_onsets.xls",
                                sheet="Ca2+",n_max=24)
capron_ages_NGRIP_dust = capron_onsets_NGRIP_dust$`t1 (50%)`


capron_onsets_NEEM_d18O = read_excel("datasets_used/Capron_onsets.xls",
                                      sheet="d18O",skip=26)
capron_ages_NEEM_d18O = capron_onsets_NEEM_d18O$`t1 (50%)`

capron_onsets_NEEM_dust = read_excel("datasets_used/Capron_onsets.xls",
                                      sheet="Ca2+",skip=25)
capron_ages_NEEM_dust = capron_onsets_NEEM_dust$`t1 (50%)`


#visualize dating uncertainty with comparisons to Rasmussen, Buizert and Capron
library(stringr)
par(mfrow=c(5,6),mar=c(2,2,1.5,1))
capind = numeric(29); buiind = numeric(29)

capind = c(NA,2,3, 4,5,6,NA,NA,7,8,NA,9,10,11,NA,NA,NA,NA,NA,12,13,14,NA,NA,15,NA,NA,16,17)
buiind = c(NA,2,NA,3,4,6,NA,NA,8,9,NA,11,12,13,NA,NA,NA,NA,NA,15,16,17,NA,NA,19,NA,NA,20,21)
for(i in 1:29){

  dens2 = density(agesimmatrix_d18O[i,]); dens2 = data.frame(x=dens2$x,y=dens2$y)#/max(dens2$y))
  if(do.dust){
    dens1 = density(agesimmatrix_dust[i,]); dens1 = data.frame(x=dens1$x,y=dens1$y)#/max(dens1$y))
    xlim = range(dens1$x,dens2$x);ylim = range(dens1$y,dens2$y)
    plot(dens2,type="l",col=1,xlab="Onset year (b2k)",ylab="Density",
         main=str_sub(event_intervals[i,2],str_locate(event_intervals[i,2],"GI-")[1]),
         xlim=xlim,ylim=ylim)
    lines(dens1,type="l",col="gray")
  }else{
    xlim = range(dens2$x);ylim = range(dens2$y)
    plot(dens2,type="l",col=1,xlab="Onset year (b2k)",ylab="Density",
         main=str_sub(event_intervals[i,2],str_locate(event_intervals[i,2],"GI-")[1]),
         xlim=xlim,ylim=ylim)
  }


  abline(v=event_intervals$GICC_age.yb2k[i],col="blue")
  if(!is.na(buiind[i])){
    abline(v=buizert_onsets$Age[buiind[i]],col="green")
    #abline(v=buizert_ages[buiind[i]],col="green")
  }


  # abline(v=capron_ages_NGRIP_d18O,col="orange")
  # abline(v=capron_ages_NGRIP_dust,col="orange",lty=3)
  if(!is.na(capind[i])){
    abline(v=capron_ages_NGRIP_d18O[capind[i]],col="red")
    abline(v=capron_ages_NGRIP_dust[capind[i]],col="pink",lty=1)
  }


}
par(mar=c(0,0,0,0));plot(-1,axes=FALSE,xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
if(do.dust){
  legend(0,1,legend=c("d18O onset posterior","Ca2+ onset posterior",
                      "Rasmussen onset","Buizert onset",
                      "Capron d18O onset","Capron Ca2+ onset"),
         col=c("black","gray", "blue", "green", "red", "pink"),
         lty=c(1,1,1,1,1,1), cex=0.65,bty="n")

}else{
  legend(0,1,legend=c("d18O onset posterior",
                      "Rasmussen onset","Buizert onset",
                      "Capron d18O onset","Capron Ca2+ onset"),
         col=c("black", "blue", "green", "red", "pink"),
         lty=c(1,1,1,1,1), cex=0.65,bty="n")

}

par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)

## Format all results to a data.frame
rownames=c()
for(i in 1:29){
  rownames=c(rownames,str_sub(event_intervals[i,2],str_locate(event_intervals[i,2],"GI-")[1]))
}

if(do.dust){
  colnames = colnames(depthstats)[1:7]
  colnames = c("d18O_mean","d18O_q0.025","d18O_q0.975",
               "Ca_mean","Ca_q0.025","Ca_q0.975")
  depthresults = t(depthstats)[,2:7]
  ageresults = t(agestats[2:7,])
}else{
  colnames = colnames(depthstats)[1:4]
  colnames = c("d18O_mean","d18O_q0.025","d18O_q0.975")
  depthresults = t(depthstats)[,2:4]
  ageresults = t(agestats[2:4,])
}

colnames(depthresults)=colnames
rownames(depthresults)=rownames

colnames(ageresults)=colnames
rownames(ageresults)=rownames

