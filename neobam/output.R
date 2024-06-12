#' Write neoBAM output to NetCDF file.
#'
#' discharge time series
#' posteriors (mean and sd): r, logn, logWb, logDb

# Libraries
# library(RNetCDF,lib.loc='/home/cjgleason_umass_edu/.conda/pkgs/r-rnetcdf-2.6_2-r42h498a2f1_0/lib/R/library/',quietly=TRUE,warn.conflicts = FALSE)
library(dplyr)

# Constants
FILL = -999999999999

#' Insert NA values back into discharge list vectors to account for invalid
#' time steps.
#'
#' @param discharge list of discharge time series to insert NA at time steps
#' @param invalid_times list of invalid time indexes
#'
#' @return list of data with NA in place of invalid nodes
concatenate_invalid = function(discharge, invalid_times) {

  # Time-level data
  for (index in invalid_times) {
    discharge = append(discharge, NA, after=index-1)
  }
  return(discharge)
}

#' Write posteriors data to NetCDF file.  lapply(list(1,2,3), write_posteriors, nc_out=nc_out, posteriors=posteriors)
#'
#' @param chain integer number of neoBAM run
#' @param nc_out NetCDF file pointer to write to
#' @param posteriors list of posteriors
write_posteriors = function(nc_out, posteriors, is_valid, out_data) {

  # Chain
  # Posteriors
  r = tryCatch(
    error = function(cond) grp.def.nc(nc_out, "r"),
    grp.inq.nc(nc_out, "r")$self
  )
  var.def.nc(r, "mean", "NC_DOUBLE", "nx")
  att.put.nc(r, "mean", "_FillValue", "NC_DOUBLE", FILL)
  var.def.nc(r, "sd", "NC_DOUBLE", NA)
  att.put.nc(r, "sd", "_FillValue", "NC_DOUBLE", FILL)

  if(is_valid){
    var.put.nc(r, "mean", posteriors$r$mean)
    var.put.nc(r, "sd", posteriors$r$sd)
  }else{
    var.put.nc(r, "mean", rep(FILL, length(out_data$nt)))
    var.put.nc(r, "sd", FILL)

  }

  logn = tryCatch(
    error = function(cond) grp.def.nc(nc_out, "logn"),
    grp.inq.nc(nc_out, "logn")$self
  )
  var.def.nc(logn, "mean", "NC_DOUBLE", "nx")
  att.put.nc(logn, "mean", "_FillValue", "NC_DOUBLE", FILL)
  var.def.nc(logn, "sd", "NC_DOUBLE", NA)
  att.put.nc(logn, "sd", "_FillValue", "NC_DOUBLE", FILL)

  if(is_valid){
    var.put.nc(logn, "mean", posteriors$logn$mean)
    var.put.nc(logn, "sd", posteriors$logn$sd)
  }else{
    var.put.nc(logn, "mean", rep(FILL, length(out_data$nt)))
    var.put.nc(logn, "sd", FILL)

  }

  logWb = tryCatch(
    error = function(cond) grp.def.nc(nc_out, "logWb"),
    grp.inq.nc(nc_out, "logWb")$self
  )
  print('logwb')
  var.def.nc(logWb, "mean", "NC_DOUBLE", "nx")
  att.put.nc(logWb, "mean", "_FillValue", "NC_DOUBLE", FILL)
  var.def.nc(logWb, "sd", "NC_DOUBLE", NA)
  att.put.nc(logWb, "sd", "_FillValue", "NC_DOUBLE", FILL)

  if(is_valid){
    var.put.nc(logWb, "mean", posteriors$logn$mean)
    var.put.nc(logWb, "sd", posteriors$logn$sd)
  }else{
    var.put.nc(logWb, "mean", rep(FILL, length(out_data$nt)))
    var.put.nc(logWb, "sd", FILL)

  }

  logDb = tryCatch(
    error = function(cond) grp.def.nc(nc_out, "logDb"),
    grp.inq.nc(nc_out, "logDb")$self
  )
  var.def.nc(logDb, "mean", "NC_DOUBLE", "nx")
  att.put.nc(logDb, "mean", "_FillValue", "NC_DOUBLE", FILL)
  # var.put.nc(logDb, "mean", posteriors$logDb$mean)
  var.def.nc(logDb, "sd", "NC_DOUBLE", NA)
  att.put.nc(logDb, "sd", "_FillValue", "NC_DOUBLE", FILL)
  # var.put.nc(logDb, "sd", posteriors$logDb$sd)

  if(is_valid){
    var.put.nc(logDb, "mean", posteriors$logDb$mean)
    var.put.nc(logDb, "sd", posteriors$logDb$sd)
  }else{
    var.put.nc(logDb, "mean", rep(FILL, length(out_data$nt)))
    var.put.nc(logDb, "sd", FILL)

  }

}

#' Write discharge data to NetCDF file.
#'
#' @param chain integer number that indicates neoBAM run
#' @param nc_out NetCDF file pointer to write to
#' @param discharge list of discharge values
  # write_output(out_data, neobam_output$posteriors, neobam_output$posterior_Q, neobam_output$posterior_Qsd, in_data$valid)
write_discharge = function(nc_out, discharge,discharge_sd,is_valid) {

  # Chain


  # Discharge
  q = tryCatch(
    error = function(cond) grp.def.nc(nc_out, "q"),
    grp.inq.nc(nc_out, "q")$self
  )
  # var.def.nc(q, "name", "NC_DOUBLE", "nt")
  # att.put.nc(q, "name", "_FillValue", "NC_DOUBLE", FILL)
  # discharge[is.nan(discharge)] = NA
  # discharge <- discharge %>% mutate_all(~ifelse(is.nan(.), NA, .))

  # print(do.call(rbind,discharge))
  # print(dim(discharge))

  var.def.nc(q, "q", "NC_DOUBLE", "nt")
  att.put.nc(q, "q", "_FillValue", "NC_DOUBLE", FILL)
  if(is_valid){
      var.put.nc(q, "q", discharge)

  

  }else{
    var.put.nc(q,"q",FILL)
  }

  var.def.nc(q, "q_sd", "NC_DOUBLE", NA)
  att.put.nc(q, "q_sd", "_FillValue", "NC_DOUBLE", FILL)
  if(is_valid){
      var.put.nc(q, "q_sd", discharge_sd)
  }else{
    var.put.nc(q,"q_sd",FILL)
  }

}

#' Write discharge and posteriors to NetCDF file.
#'
#' @param data named list of metadata
#' @param posteriors list of posterior "chains"
#' @param discharge list of discharge "chains"
#' @param out_dir string to output directory
#'
#' @export
write_output = function(data, posteriors, discharge,discharge_sd, out_dir, is_valid) {

  print('writing output...')

print('here are invalid')
# print(data$invalid_times)
  # Concatenate invalid times back into discharge
  # discharge = lapply(discharge, concatenate_invalid, invalid_times=data$invalid_times)
  discharge = concatenate_invalid(discharge, invalid_times = data$invalid_times)

  # File creation
  nc_file = paste(out_dir, paste0(data$reach_id, "_geobam.nc"), sep=.Platform$file.sep)
  nc_out = create.nc(nc_file, format="netcdf4")

  # Global attributes
  att.put.nc(nc_out, "NC_GLOBAL", "reach_id", "NC_INT64", data$reach_id)
  att.put.nc(nc_out, "NC_GLOBAL", "node_ids", "NC_DOUBLE", unlist(data$node_ids))

  # Dimensions
  dim.def.nc(nc_out, "nt", length(data$nt))
  var.def.nc(nc_out, "nt", "NC_INT", "nt")
  att.put.nc(nc_out, "nt", "units", "NC_STRING", "time")
  var.put.nc(nc_out, "nt", data$nt)

  dim.def.nc(nc_out, "nx", length(posteriors$r$mean))
  var.def.nc(nc_out, "nx", "NC_INT", "nx")
  att.put.nc(nc_out, "nx", "units", "NC_STRING", "num_nodes")
  var.put.nc(nc_out, "nx", seq(from = 0, by = 1, length.out = length(posteriors$r$mean))
)

  # dim.def.nc(nc_out, "_dis?", length(discharge))
  # var.def.nc(nc_out, "_dis?", "NC_DOUBLE", "_dis?")
  # att.put.nc(nc_out, "_dis?", "units", "NC_DOUBLE", "wtf_is_the_discharge_dim")
  # var.put.nc(nc_out, "_dis?", discharge)


  # Write data
  print('writing post...')
  # lapply(list(1,2,3), write_posteriors, nc_out=nc_out, posteriors=posteriors)
  write_posteriors(nc_out = nc_out, posteriors=posteriors, is_valid=is_valid, out_data=data)
  print('writing discharge...')

  # Discharge
  # lapply(list(1,2,3), write_discharge, nc_out=nc_out, discharge=discharge)
  write_discharge(nc_out = nc_out, discharge=discharge,discharge_sd=discharge_sd, is_valid=is_valid)

  # Close file
  close.nc(nc_out)
}
