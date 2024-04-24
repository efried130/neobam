
run_neobam_stan = function(neobam_data_and_priors,sourcefile){
    
  library(BH,quietly=TRUE,warn.conflicts = FALSE)
  library(rstan,quietly=TRUE,warn.conflicts = FALSE)
  library(dplyr,quietly=TRUE,warn.conflicts = FALSE)
  library(tidyr, lib.loc = "/nas/cee-water/cjgleason/r-lib/",quietly=TRUE,warn.conflicts = FALSE)
  library(stringr,quietly=TRUE,warn.conflicts = FALSE)
    
    rstan_options(auto_write = TRUE)


    
  return_posterior_mean_sd= function(parameter,stan_output){
      
      
    
    var_names=gsub("\\[.*\\]" ,"",   row.names(stan_output) )
    parameter_index=!is.na(str_match(var_names,parameter))
    string_length_index=unlist(lapply(var_names,nchar)) == nchar(parameter)
    
    combined_index=which(parameter_index*string_length_index==1)
      
      #print(stan_output[combined_index,])
    
    posterior_mean=stan_output[combined_index,]$mean
    posterior_sd=mean(stan_output[combined_index,]$sd^2)
    
    return(list('mean'=posterior_mean,'sd'=posterior_sd))
    
  }
    
   
    
  
iter=neobam_data_and_priors$iter
    

        
  fit1 <- stan(
    
    file = sourcefile,  # Stan program
    data = neobam_data_and_priors,    # named list of data
    chains = 3,             # number of Markov chains
    warmup = floor(0.1*iter),          # number of warmup iterations per chain
    iter = iter,   # total number of iterations per chain
    cores = 1,              # number of cores (could use one per chain)
    refresh = 0
    #        boost_lib='/home/cjgleason_umass_edu/R/x86_64-pc-linux-gnu-library/4.0',
   #   eigen_lib='/home/cjgleason_umass_edu/R/x86_64-pc-linux-gnu-library/4.0'
  )
  
  output= data.frame(rstan::summary(fit1)$summary)
  #OK! we now need a list of all the posteriors so we can return them to stuff back into the flow law
  # posterior_list=c('r','logn','betaslope','betaint','logbeta','logQ')
     posterior_list=c('r','logn','logWb','logDb','logQ')
    
  posteriors=lapply(posterior_list,return_posterior_mean_sd,output)
  names(posteriors)=posterior_list
    
    nchains <- fit1@sim$chains
     chain <- 1:nchains
      alpha <- 1 - 0.95
    
  hydrograph_posterior=  rstan::extract(fit1, "logQ", permuted = FALSE) %>%
    reshape2::melt()%>%
    dplyr::mutate(chains = gsub("^chain:", "", .data$chains)) %>%
    dplyr::filter(.data$chains %in% chain) %>%
    dplyr::mutate(value = exp(.data$value)) %>%
    dplyr::group_by(.data$parameters) %>%
    dplyr::summarize(mean = mean(.data$value),
              conf.low = quantile(.data$value, alpha / 2),
              conf.high = quantile(.data$value, 1 - (alpha / 2)),
              sd = sd(.data$value)) %>%
    dplyr::rename(time = .data$parameters) %>%
    dplyr::mutate(time = gsub("^logQ\\[", "", .data$time),
           time = gsub("\\]$", "", .data$time),
           time = as.numeric(.data$time)) %>%
    dplyr::arrange(.data$time)
    
    # posteriors$logbeta$mean = (posteriors$r$mean * posteriors$betaslope$mean) + posteriors$betaint$mean
    # posteriors$logbeta$sd = (posteriors$r$sd * posteriors$betaslope$sd) + posteriors$betaint$sd


  return(list('posteriors'=posteriors,'hydrograph_posterior'=hydrograph_posterior))

}

remake_discharge =function (Wobs,Sobs,posteriors){
    
    

  r=posteriors$r$mean
  logn=posteriors$logn$mean
  logbeta=posteriors$logDb$mean -(r*posteriors$logWb$mean)
  
    logW=log(Wobs)
    logS=log(Sobs)
    
finalQ=matrix(NA,nrow=nrow(logS),ncol=ncol(logS))


    
for (i in 1:nrow(logS)){

logQ= 0.5*logS[i,] -
     logn[i] +
     1.66*logbeta[i] +
     1.66*(log(r[i]/(r[i]+1))) +
     ((1.66*r[i]) +1)*log(Wobs[i,])  
    
    
finalQ[i,]=exp(logQ)
}

    
finalQ[is.infinite(finalQ)]=NA
    
    
finalQ=colMeans(finalQ,na.rm=TRUE) #take the mean for a given time t
    #in linspace

  
return(finalQ)

  
}

norm_to_lognorm=function(mu,sigma){
    lognorm_mu=2*log(mu)-0.5*log(mu^2+sigma^2)
    lognorm_sigma= log(mu^2 + sigma^2) - 2*(log(mu))
    return(list('mu'=lognorm_mu,'sigma'=lognorm_sigma))
    
}








