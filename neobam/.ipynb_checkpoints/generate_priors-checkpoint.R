generate_priors=function(parameters){

   
    
library(settings, lib.loc = "/nas/cee-water/cjgleason/r-lib/",quietly = TRUE,warn.conflicts=FALSE)
library(whisker, lib.loc = "/nas/cee-water/cjgleason/r-lib/",quietly = TRUE,warn.conflicts=FALSE)

Wobs=parameters$swot_data$width
Sobs=parameters$swot_data$slope2
Q_priors=parameters$sos_data$Qpriors
reachid=parameters$reachid
date=parameters$swot_data$time

    
    #need dAdW and dHdW
    #dH, dW, and dA are relative to minimum. simple finite difference should work



    #create a bamdata object for the prior generation
      Qmean=exp(Q_priors$logQ_hat)
      bamdata = bam_data(h=Hobs,w=Wobs, dA=dAobs,Qhat=as.vector(Qmean),s=Sobs) 
  #make full BAM priors (increase knowledge content as much as possible)
  makeBamPriors = function(bamdata,Q_priors){
    bampriors = bam_priors(bamdata = bamdata, logQ_sd = Q_priors$logQ_sd, variant = "amhg",classification='expert')
    return(bampriors)
  }
  
    #create other BAM priors (can be adjusted)
    bampriors = makeBamPriors(bamdata,Q_priors)

    #print(bampriors$other_priors)
    norm_to_lognorm=function(mu,sigma){
    lognorm_mu=2*log(mu)-0.5*log(mu^2+sigma^2)
    lognorm_sigma= log(mu^2 + sigma^2) - 2*(log(mu))
    return(list('mu'=lognorm_mu,'sigma'=lognorm_sigma))}


# #nudge priors
# #the max width is less than bankfull given the precentile filter.
# bampriors$river_type_priors$lowerbound_logn= log(0.02)
# bampriors$river_type_priors$upperbound_logn= log(0.04)# nudging
# bampriors$river_type_priors$logn_hat=rep(log(0.05),times=nrow(Wobs))
# bampriors$river_type_priors$logn_sd = rep(2,times=nrow(Hobs)) #bampriors$river_type_priors$logr_hat /2
    
# Q_priors$logQ_hat=rep(Q_priors$logQ_hat,times=nrow(Hobs))
# Q_priors$logQ_sd=rep(Q_priors$logQ_sd,times=nrow(Hobs))

if(any(is.na(c(Q_priors$logQ_sd,Q_priors$logQ_hat)))){return(list('neo_priors'=NA))}
                
hasdat=matrix(as.integer(Hobs>0),nrow=nrow(Hobs),ncol=ncol(Hobs))
ntot=sum(hasdat)

#for stan, NA to zero
Sobs[is.na(Sobs)]=0
Wobs[is.na(Wobs)]=0
    
                     
   
   
return(list('neo_priors'=bampriors,
            'date'=date,
            'Hobs'=Hobs,
            'Sobs'=Sobs,
            'Wobs'=Wobs,
            'dAdWobs'=dAdWobs,
            'dHdWobs'=dHdWobs,
            'Wdiffobs'=Wdiffobs,
            'dHobs'=dHobs,
            'dAobs'=dAobs,
            'reachid'=reachid,
            'ntot'=ntot,
            'hasdat'=hasdat,
            'Qobs'=Qobs,
            'Q_priors'=Q_priors))
    
}
