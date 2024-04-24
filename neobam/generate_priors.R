generate_priors=function(parameters){

   
    
library(settings, lib.loc = "/nas/cee-water/cjgleason/r-lib/",quietly = TRUE,warn.conflicts=FALSE)
library(whisker, lib.loc = "/nas/cee-water/cjgleason/r-lib/",quietly = TRUE,warn.conflicts=FALSE)

Wobs=parameters$Wobs
Sobs=parameters$Sobs
Hobs=parameters$Hobs
dAdWobs=parameters$dAdWobs
dHdWobs=parameters$dHdWobs   
Wdiffobs=parameters$Wdiffobs
dAobs=parameters$dAobs
dHobs=parameters$dHobs
Qobs=parameters$Qobs
Q_priors=parameters$Qpriors
reachid=parameters$reachid
date=parameters$date
sample_size=parameters$sample_size
    
    #need dAdW and dHdW
    #dH, dW, and dA are relative to minimum. simple finite difference should work


    
     #if enough data to be subsampled, do that
      if(nrow(Hobs)>parameters$sample_size){
          order_frame=data.frame(row=1:nrow(Hobs),nt=rowSums(!is.na(Hobs)))
          keep_index=order(order_frame$nt, decreasing=TRUE)[1:sample_size]
          Hobs=Hobs[keep_index,] 
          Wobs=Wobs[keep_index,] 
          Qobs=Qobs[keep_index,] 
          dAobs=dAobs[keep_index,] 
          dHobs=dHobs[keep_index,] 
          dWobs=dWobs[keep_index,] 
      }    
    
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
Hobs[is.na(Hobs)]=0
Wobs[is.na(Wobs)]=0
    
dHobs[is.na(dHobs)]=0
dAobs[is.na(dAobs)]=0
dAdWobs[is.na(dAdWobs)]=0
    dHdWobs[is.na(dHdWobs)]=0
    Wdiffobs[is.na(Wdiffobs)]=0
                     
   
   
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
