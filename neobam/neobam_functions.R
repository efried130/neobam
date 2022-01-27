window_mode=function(Wobs,Sobs,date,perc_lower,perc_upper,nt_window,nx_sample,priors){

    # ?Handled valid observation detection in input
    # Sobsmat= matrix(Sobs,nrow=nrow(Wobs),ncol=ncol(Wobs))
    # Sobs=Sobsmat

    # #cleanup
    #   #more cleanup needed- want at least 5 nx for every nt and 5nt for every nx
    #   index_of_nodata_row=rowSums(Wobs==0) >= (ncol(Wobs)-5)
    #   Wobs=Wobs[!index_of_nodata_row,]
    #   Sobs=Sobs[!index_of_nodata_row,]
    #
    #   index_of_nodata_col=colSums(Wobs==0) >= (nrow(Wobs)-5)
    #   index_of_nodata_col[is.na(index_of_nodata_col)] = TRUE
    #
    #   Wobs=Wobs[,!index_of_nodata_col]
    #   Sobs=Sobs[,!index_of_nodata_col]
    #   date=date[!index_of_nodata_col]
    #
    #   #stop if there aren't enough data
    #   if (nrow(Wobs)<5 || ncol(Wobs) < 5 ) {return(flag=1)}

    #if there are less times than the window size
    if (ncol(Wobs)<nt_window){
      hasdat=matrix(as.integer(Wobs>0),nrow=nrow(Wobs),ncol=ncol(Wobs))
      ntot=sum(hasdat)
      W_index=1:ncol(Wobs)
    } else {
      W_mean=colMeans(Wobs, na.rm=TRUE)
      W_percentiles=quantile(W_mean, c(perc_lower,perc_upper), na.rm=TRUE)

      W_index= which(W_mean>=W_percentiles[1] & W_mean<=W_percentiles[2] )
    }
    if (length(W_index)>nt_window){
      #get a random sample of those in that percentile
      W_index=sample(W_index,nt_window,replace=FALSE)
    }

    if (nrow(Wobs)>nx_sample){
      which_rows=sample(1:nrow(Wobs),nx_sample,replace=FALSE)
    } else {
      which_rows=1:nrow(Wobs)
    }

    # Apply window to observations and priors
    Wobs=Wobs[which_rows,W_index]
    Sobs=Sobs[which_rows,W_index]
    hasdat=matrix(as.integer(Wobs>0),nrow=nrow(Wobs),ncol=ncol(Wobs))
    ntot=sum(hasdat)
    priors=lapply(priors, window_prior, which_rows)

 return(list('hasdat'=hasdat,'ntot'=ntot,'Wobs'=Wobs,'Sobs'=Sobs,'row_filter'=which_rows,'column_filter'=W_index, 'priors'=priors))
}

window_prior = function(prior, which_rows) {
  if (length(prior) > 1) {
    return(prior[which_rows])
  } else {
    return(prior)
  }
}

generate_neobam_priors_and_data = function(windowed, Q_priors,
                               neobam_prep_data,global_params,window_params) {
  #Q priors assume the user has a list of pre-calculated priors
  # logQ_hat= Q_priors$logQ_hat[windowed$column_filter] #because we sampled in time
  logQ_hat = rep(Q_priors$logQ_hat, times=length(windowed$column_filter))
  lowerbound_logQ= Q_priors$lowerbound_logQ
  upperbound_logQ=Q_priors$upperbound_logQ
  # logQ_sd= Q_priors$logQ_sd[windowed$column_filter]
  logQ_sd = rep(Q_priors$logQ_sd, times=length(windowed$column_filter))

  #performative variables
  nx=nrow(windowed$Wobs)
  nt=ncol(windowed$Wobs)
  ntot=windowed$ntot
  hasdat=windowed$hasdat
  Wobs=windowed$Wobs
  Sobs=windowed$Sobs

neobam_data_and_priors = list(

    Wobs=Wobs,
    Sobs=Sobs,
logQ_hat=logQ_hat,
lowerbound_logQ=lowerbound_logQ,
upperbound_logQ=upperbound_logQ,
logQ_sd=logQ_sd,

nx=nx,
nt=nt,
ntot=ntot,
hasdat=hasdat,
iter=global_params$iter,

Werr_sd=global_params$Werr_sd, #expected uncertainty in width observations
Serr_sd=global_params$Serr_sd, #expected uncerainty in slope observations
iter=global_params$iter, #how many iterations in the markov chain?
sigma_man=matrix(global_params$sigma_man,nrow=nrow(Wobs),ncol=ncol(Wobs)),

logWb_hat = window_params$logWb_hat,
logWb_sd = window_params$logWb_sd,
lowerbound_logWb = window_params$lowerbound_logWb,
upperbound_logWb = window_params$upperbound_logWb,

logDb_hat = window_params$logDb_hat,
logDb_sd = window_params$logDb_sd,
lowerbound_logDb = window_params$lowerbound_logDb,
upperbound_logDb = window_params$upperbound_logDb,

r_hat = window_params$r_hat,
r_sd = window_params$r_sd,
lowerbound_r = window_params$lowerbound_r,
upperbound_r = window_params$upperbound_r,

logn_hat = window_params$logn_hat,
logn_sd = window_params$logn_sd,
lowerbound_logn = window_params$lowerbound_logn,
upperbound_logn = window_params$upperbound_logn


  )


}

run_neobam = function(neobam_data_and_priors,sourcefile){
  #run it here
  library(rstan,quietly=TRUE,warn.conflicts = FALSE)
  library(dplyr,quietly=TRUE,warn.conflicts = FALSE)
  library(tidyr, quietly=TRUE,warn.conflicts = FALSE)
  library(stringr,quietly=TRUE,warn.conflicts = FALSE)

  return_posterior_mean_sd= function(parameter,stan_output){


    var_names=gsub("\\[.*\\]" ,"",   row.names(stan_output) )
    parameter_index=!is.na(str_match(var_names,parameter))
    string_length_index=unlist(lapply(var_names,nchar)) == nchar(parameter)

    combined_index=which(parameter_index*string_length_index==1)

    posterior_mean=mean(stan_output[combined_index,]$mean)
    posterior_sd=sqrt(mean(stan_output[combined_index,]$sd^2))

    return(list('mean'=posterior_mean,'sd'=posterior_sd))

  }

    iter=neobam_data_and_priors$iter

  fit1 <- stan(

    file = sourcefile,  # Stan program
    data = neobam_data_and_priors,    # named list of data
    chains = 1,             # number of Markov chains
    warmup = floor(0.1*iter),          # number of warmup iterations per chain
    iter = iter,   # total number of iterations per chain
    cores = 1,              # number of cores (could use one per chain)
    refresh = 0             # no progress shown

  )



  output= data.frame(rstan::summary(fit1)$summary)
  #OK! we now need a list of all the posteriors so we can return them to stuff back into the flow law
  posterior_list=c('r','logn','logWb','logDb')

  posteriors=lapply(posterior_list,return_posterior_mean_sd,output)
  names(posteriors)=posterior_list

  return(posteriors)

}

remake_discharge =function (Wobs,Sobs,posteriors){
  r=posteriors$r$mean
  logWb=posteriors$logWb$mean
  logDb=posteriors$logDb$mean
  logn=posteriors$logn$mean



delta=(1.67*r) +1
logQ=  (delta*log(Wobs))-(1.67*r*logWb)+(1.67*logDb)-(1.67*log((r+1)/r))-(logn)+(0.5*log(Sobs))

logQ[is.infinite(logQ)]=NA

Q=exp(colMeans(logQ,na.rm=T))



}


norm_to_lognorm=function(mu,sigma){
    lognorm_mu=2*log(mu)-0.5*log(mu^2+sigma^2)
    lognorm_sigma= log(mu^2 + sigma^2) - 2*(log(mu))
    return(list('mu'=lognorm_mu,'sigma'=lognorm_sigma))

}











