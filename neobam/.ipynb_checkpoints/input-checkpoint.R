#' neoBAM input operations
#'
#' SOS: max_q, mean_q, min_q
#' SWOT: node/time, node/width, node/slope2

# Libraries
library(RNetCDF,quietly=TRUE,warn.conflicts = FALSE)
#' Title
#'
#' @param swot_file string path to SWOT NetCDF file
#' @param sos_file string path to SOS NetCDF file
#' @param reach_id integer unique reach identifier
#'
#' @return named list of data needed to run neoBAM
#' @export
get_input = function(swot_file, sos_file, reach_id) {


  # Get SWOT
  swot_data = get_swot(swot_file)
    
    # print(swot_data$node$width)

  # Get SOS
  sos_data = get_sos(sos_file, reach_id)

  # Check for valid number of observations
  data = check_observations(swot_data, sos_data)
  # print("here first bs")
  # print(swot_data$obs_times)

  # Return list of valid observations
  if (length(data) == 0) {
    print('in get_input(), neobam has decided the data are invalid')
         # print(swot_data)
    return(list(valid=FALSE, obs_times=swot_data$obs_times,reach_id=reach_id, nx=swot_data$nx, nt=swot_data$nt, thisisdumb=1, node_ids=sos_data$Q_priors$nids))
  } else {
          # print(swot_data)
    # Create a list of data with reach identifier
    return(list(valid=TRUE, reach_id = reach_id,obs_times=swot_data$obs_times ,swot_data=data$swot_data,
                sos_data=data$sos_data, invalid_nodes=data$invalid_nodes,
                invalid_times=data$invalid_times, node_ids=sos_data$Q_priors$nids))
  }
}

#' Retrieve SWOT data.
#'
#' @param swot_file string path to SWOT NetCDF file
#'
#' @return list of width, slope2, and time matrices
get_swot = function(swot_file) {
  swot = open.nc(swot_file)
  nx = var.get.nc(swot, "nx")
  nt = var.get.nc(swot, "nt")

  node_grp = grp.inq.nc(swot, "node")$self
  width = t(var.get.nc(node_grp, "width"))
  slope2 = t(var.get.nc(node_grp, "slope2"))
  time = t(var.get.nc(node_grp, "time"))

  reach_grp = grp.inq.nc(swot, "reach")$self
  obs_times = t(var.get.nc(reach_grp, "time_str"))   


  close.nc(swot)

  return(list(nx=nx, nt=nt, width=width, slope2=slope2, time=time, obs_times=obs_times))

}

reload_sos = function(sos_file, retry_number){
  tries = retry_number
  while (tries > 0){
    tryCatch (
      {
        sos = open.nc(sos_file)
        # print("sos loaded...")
        tries = 0

      },
      error=function(e) {
              message('An Error Occurred')
              # print(e)
              tries = tries - 1
              Sys.sleep(runif(1, min=100, max=500))

          },
      warning=function(w) {
            message('A Warning Occurred')
            # print(w)
            return(NA)
        
  }
    )

  }
  return(sos)

}


#' Retrieve SOS data.
#'
#' @param sos_file string path to SOS NetCDF file
#' @param reach_id integer unique reach identifier
#'
#' @return list of Q priors
#' @export
get_sos = function(sos_file, reach_id) {

  Q_priors = list()
  sos = reload_sos(sos_file, 5)



  tries = 5
  while (tries >0){
    tryCatch ( {
    # print('trying...')
    reach_grp = grp.inq.nc(sos, "reaches")$self
    rids = var.get.nc(reach_grp, "reach_id")
    index = which(rids == reach_id, arr.ind=TRUE)

    node_grp = grp.inq.nc(sos, "nodes")$self
    nrids = var.get.nc(node_grp, "reach_id")
    indexes = which(nrids == reach_id, arr.ind=TRUE)
    nids = var.get.nc(node_grp, "node_id")
    Q_priors$nids = nids[indexes]

    model_grp = grp.inq.nc(sos, "model")$self
    # print("made it to the model group")
    # print(model_grp)
    Q_priors$logQ_hat = log(var.get.nc(model_grp, "monthly_q")[,index])
  
    Q_priors$upperbound_logQ = log(var.get.nc(model_grp, "max_q")[index])
    min_q = var.get.nc(model_grp, "min_q")[index]   # Check action taken
    # if(is.null(min_q)){Q_priors$lowerbound_logQ = NA}
    if(is.na(min_q)){Q_priors$lowerbound_logQ = NA} else {
        Q_priors$lowerbound_logQ = log(min_q)
        if(min_q ==0){
            min_q=0.01
         Q_priors$lowerbound_logQ = log(min_q)}
        
    }
        
      
    

    r_grp = grp.inq.nc(sos, "gbpriors/reach")$self
    Q_priors$logQ_sd = var.get.nc(r_grp, "logQ_sd")[index]

    window_params = list()
    n_grp = grp.inq.nc(sos, "gbpriors/node")$self
    window_params$logWb_hat = var.get.nc(n_grp, "logWb_hat")[indexes]
    window_params$logWb_sd = var.get.nc(n_grp, "logWb_sd")[indexes]
    window_params$lowerbound_logWb = min(var.get.nc(n_grp, "lowerbound_logWb")[index])
    window_params$upperbound_logWb = max(var.get.nc(n_grp, "upperbound_logWb")[index])

    window_params$logDb_hat = var.get.nc(n_grp, "logDb_hat")[indexes]
    window_params$logDb_sd = var.get.nc(n_grp, "logDb_sd")[indexes]
    window_params$lowerbound_logDb = min(var.get.nc(n_grp, "lowerbound_logDb")[index])
    window_params$upperbound_logDb = max(var.get.nc(n_grp, "upperbound_logDb")[index])

    window_params$r_hat = exp(var.get.nc(n_grp, "logr_hat")[indexes])
    window_params$r_sd = exp(var.get.nc(n_grp, "logr_sd")[indexes])
    window_params$lowerbound_r = min(exp(var.get.nc(n_grp, "lowerbound_logr")[index]))
    window_params$upperbound_r = max(exp(var.get.nc(n_grp, "upperbound_logr")[index]))

    window_params$logn_hat = var.get.nc(n_grp, "logn_hat")[indexes]
    window_params$logn_sd = var.get.nc(n_grp, "logn_sd")[indexes]
    window_params$lowerbound_logn = min(var.get.nc(n_grp, "lowerbound_logn")[index])
    window_params$upperbound_logn = max(var.get.nc(n_grp, "upperbound_logn")[index])
    close.nc(sos)
    tries = 0


  },
      error=function(e) {
              message('An Error Occurred, reloading sos')
              close.nc(sos)


              # print(e)
              tries = tries - 1
              Sys.sleep(runif(1, min=10, max=600))
              sos = reload_sos(sos_file, 5)
          },
      warning=function(w) {
            message('A Warning Occurred')
            # print(w)
            return(NA)
      }
  )
}

  # print("Read was successful...")

  return(list(Q_priors=Q_priors, window_params=window_params))

}
#' Checks if observation data is valid.
#'
#' @param swot_data named list of SWOT observations
#' @param sos_data named list of priors (Q and window_params)
#'
#' @return list of valid observations or an empty list if there are none
check_observations = function(swot_data, sos_data) {

  # Q priors
  qhat = sos_data$Q_priors$logQ_hat
  qmax = sos_data$Q_priors$upperbound_logQ
  qmin = sos_data$Q_priors$lowerbound_logQ
  qsd = sos_data$Q_priors$logQ_sd

  if (all(is.na(qhat[[1]])) || is.na(qmax[[1]]) || is.na(qmin[[1]]) || is.na(qsd[[1]])) { return(vector(mode = "list")) }
    
  
  # SWOT data
  swot_data$width[swot_data$width < 0] = NA
  swot_data$slope2[swot_data$slope2 < 0] = NA
  # print(swot_data$width)
  # print(swot_data$slope2)
  # print(swot_data$time)
    
 
  invalid = get_invalid_nodes_times(swot_data$width, swot_data$slope2, swot_data$time)
 

  # Return valid data (or empty list if invalid)
       #here, we pass the inpute data do the file, so this is where we rewrite missing qhats
  return(remove_invalid(swot_data, sos_data, invalid$invalid_nodes, invalid$invalid_times))
}

#' Retrieve invalid nodes and time steps
#'
#' @param width matrix of SWOT width obs
#' @param slope2 matrix of SWOT slope2 obs
#' @param time matrix of SWOT time obs
#'
#' @return named list of invalid nodes and times indexes
get_invalid_nodes_times = function(width, slope2, time) {
    
    # print(width)
    # print(slope2)
    # print(time)


  invalid_width = get_invalid(width)
  # print(invalid_width)
  invalid_slope2 = get_invalid(slope2)
  invalid_time = get_invalid(time)
  invalid_nodes = unique(c(which(invalid_width$invalid_nodes == TRUE),
                           which(invalid_slope2$invalid_nodes == TRUE),
                           which(invalid_time$invalid_nodes == TRUE)))
  invalid_times = unique(c(which(invalid_width$invalid_times == TRUE),
                           which(invalid_slope2$invalid_times == TRUE),
                           which(invalid_time$invalid_times == TRUE)))
  return(list(invalid_nodes=invalid_nodes, invalid_times=invalid_times))
}

#' Remove invalid nodes and times from SWOT observations and SOS node priors
#'
#' @param swot_data named list of SWOT observations
#' @param sos_data named list of priors (Q and window_params)
#' @param invalid_nodes vector of invalid node indexes
#' @param invalid_times vector of invalid time indexes
#'
#' @return named list of SWOT observations, Q priors, and other priors
remove_invalid = function(swot_data, sos_data, invalid_nodes, invalid_times){
  # print(invalid_nodes)
    #make a dataframe with a month keyf ield
    Q_hat_df= data.frame(logQ_hat=sos_data$Q_priors$logQ_hat, month=as.numeric(1:12))
    #change the date to a datetime object  using the !#$!@#$ origin of the data
    #time is nx by nt, we want an nt vector only
    obs_date= as.Date(as.Date(swot_data$time[1,]/86400,origin = '2000-01-01'),format='%Y%m%d')
  
    #make it a month
    obs_month=data.frame(month=as.numeric(format(obs_date,'%m')))
      # print(obs_month)
    #now join the Qhat_df by month
    Q_hat_month=left_join(obs_month,Q_hat_df,by='month')
    #some of these will be blank, so we need to drop the NAs in order to make the data work out
    # print(Q_hat_month)
     Q_hat_month$logQ_hat[is.na(Q_hat_month$logQ_hat)]=mean(Q_hat_month$logQ_hat,na.rm=TRUE)
#     print(Q_hat_month)
#                                 bonk
    
  # All valid
  if (identical(invalid_nodes, integer(0)) && identical(invalid_times, integer(0))) {
    print('all valid')
      sos_data$Q_priors$logQ_hat=Q_hat_month$logQ_hat
    return(list(swot_data=swot_data, sos_data=sos_data,
                invalid_nodes=invalid_nodes, invalid_times=invalid_times))
    # Valid nodes
  } else if (identical(invalid_nodes, integer(0))) {
    print('valid nodes')
    swot_data$width = swot_data$width[, -invalid_times]
    swot_data$slope2 = swot_data$slope2[, -invalid_times]
    swot_data$time = swot_data$time[, -invalid_times]
    sos_data$Q_priors$logQ_hat=Q_hat_month$logQ_hat[-invalid_times]

    # Valid time steps
  } else if (identical(invalid_times, integer(0))) {
    print('valid time')
    swot_data$width = swot_data$width[-invalid_nodes,]
    swot_data$slope2 = swot_data$slope2[-invalid_nodes,]
    swot_data$time = swot_data$time[-invalid_nodes,]
    sos_data$window_params$logWb_hat = sos_data$window_params$logWb_hat[-invalid_nodes]
    sos_data$window_params$logWb_sd = sos_data$window_params$logWb_sd[-invalid_nodes]
    sos_data$window_params$logDb_hat = sos_data$window_params$logDb_hat[-invalid_nodes]
    sos_data$window_params$logDb_sd = sos_data$window_params$logDb_sd[-invalid_nodes]
    sos_data$window_params$r_hat = sos_data$window_params$r_hat[-invalid_nodes]
    sos_data$window_params$r_sd = sos_data$window_params$r_sd[-invalid_nodes]
    sos_data$window_params$logn_hat = sos_data$window_params$logn_hat[-invalid_nodes]
    sos_data$window_params$logn_sd = sos_data$window_params$logn_sd[-invalid_nodes]

    # Both invalid
  } else {
    print('both invalid')
    swot_data$width = swot_data$width[-invalid_nodes, -invalid_times]
    swot_data$slope2 = swot_data$slope2[-invalid_nodes, -invalid_times]
    swot_data$time = swot_data$time[-invalid_nodes, -invalid_times]
    sos_data$Q_priors$logQ_hat=Q_hat_month$logQ_hat[-invalid_times]  
      
      
    sos_data$window_params$logWb_hat = sos_data$window_params$logWb_hat[-invalid_nodes]
    sos_data$window_params$logWb_sd = sos_data$window_params$logWb_sd[-invalid_nodes]
    sos_data$window_params$logDb_hat = sos_data$window_params$logDb_hat[-invalid_nodes]
    sos_data$window_params$logDb_sd = sos_data$window_params$logDb_sd[-invalid_nodes]
    sos_data$window_params$r_hat = sos_data$window_params$r_hat[-invalid_nodes]
    sos_data$window_params$r_sd = sos_data$window_params$r_sd[-invalid_nodes]
    sos_data$window_params$logn_hat = sos_data$window_params$logn_hat[-invalid_nodes]
    sos_data$window_params$logn_sd = sos_data$window_params$logn_sd[-invalid_nodes]

  }
    


  # Return list to indicate invalid data
  if (is.null(dim(swot_data$width)) || nrow(swot_data$width) < 3 || ncol(swot_data$width) < 3 ) { return(vector(mode = "list")) }
  if (is.null(dim(swot_data$slope2)) || nrow(swot_data$slope2) < 3 || ncol(swot_data$slope2) < 3 ) { return(vector(mode = "list")) }
  if (is.null(dim(swot_data$time)) || nrow(swot_data$time) < 3 || ncol(swot_data$time) < 3 ) { return(vector(mode = "list")) }

    
    # print(sos_data$Q_priors$logQ_hat)
    # bonk

    
  # Return list of remaining valid observation data
  return(list(swot_data=swot_data, sos_data=sos_data,
              invalid_nodes=invalid_nodes, invalid_times=invalid_times))

}

#' Checks if observation parameter has valid nx (nodes) and valid nt (time steps)
#'
#' @param obs matrix
#'
#' @return list of invalid nodes and invalid time steps
get_invalid = function(obs) {

  # Determine invalid nx and nt for obs

  invalid_nodes = rowSums(is.na(obs)) >= (ncol(obs) - 3)
  # print('rowsums')
  # print(rowSums(is.na(obs)))
  # print('colsums')
  # print(colSums(is.na(obs)))

  invalid_times = colSums(is.na(obs)) >= (nrow(obs) - 3)
  # invalid_nodes = 0
  # invalid_times = 0


  return(list(invalid_nodes=invalid_nodes, invalid_times=invalid_times))

}