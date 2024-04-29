#' Write neoBAM output to NetCDF file.
#'
#' discharge time series
#' posteriors (mean and sd): r, logn, logWb, logDb

# Libraries
# library(RNetCDF)


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

#' Write posteriors data to NetCDF file.
#'
#' @param chain integer number of neoBAM run
#' @param nc_out NetCDF file pointer to write to
#' @param posteriors list of posteriors
write_posteriors = function(nc_out, posteriors) {

  # Chain

  # Posteriors
  r = tryCatch(
    error = function(cond) grp.def.nc(nc_out, "r"),
    grp.inq.nc(nc_out, "r")$self
  )
  var.def.nc(r, mean, "NC_DOUBLE", NA)
  att.put.nc(r, mean, "_FillValue", "NC_DOUBLE", FILL)
  var.put.nc(r, mean, posteriors$r$mean)
  var.def.nc(r, sd, "NC_DOUBLE", NA)
  att.put.nc(r, sd, "_FillValue", "NC_DOUBLE", FILL)
  var.put.nc(r, sd, posteriors$r$sd)

  logn = tryCatch(
    error = function(cond) grp.def.nc(nc_out, "logn"),
    grp.inq.nc(nc_out, "logn")$self
  )
  var.def.nc(logn, mean, "NC_DOUBLE", NA)
  att.put.nc(logn, mean, "_FillValue", "NC_DOUBLE", FILL)
  var.put.nc(logn, mean, posteriors$logn$mean)
  var.def.nc(logn, sd, "NC_DOUBLE", NA)
  att.put.nc(logn, sd, "_FillValue", "NC_DOUBLE", FILL)
  var.put.nc(logn, sd, posteriors$logn$sd)

  logWb = tryCatch(
    error = function(cond) grp.def.nc(nc_out, "logWb"),
    grp.inq.nc(nc_out, "logWb")$self
  )
  var.def.nc(logWb, mean, "NC_DOUBLE", NA)
  att.put.nc(logWb, mean, "_FillValue", "NC_DOUBLE", FILL)
  var.put.nc(logWb, mean, posteriors$logWb$mean)
  var.def.nc(logWb, sd, "NC_DOUBLE", NA)
  att.put.nc(logWb, sd, "_FillValue", "NC_DOUBLE", FILL)
  var.put.nc(logWb, sd, posteriors$logWb$sd)

  logDb = tryCatch(
    error = function(cond) grp.def.nc(nc_out, "logDb"),
    grp.inq.nc(nc_out, "logDb")$self
  )
  var.def.nc(logDb, mean, "NC_DOUBLE", NA)
  att.put.nc(logDb, mean, "_FillValue", "NC_DOUBLE", FILL)
  var.put.nc(logDb, mean, posteriors$logDb$mean)
  var.def.nc(logDb, sd, "NC_DOUBLE", NA)
  att.put.nc(logDb, sd, "_FillValue", "NC_DOUBLE", FILL)
  var.put.nc(logDb, sd, posteriors$logDb$sd)

}

#' Write discharge data to NetCDF file.
#'
#' @param chain integer number that indicates neoBAM run
#' @param nc_out NetCDF file pointer to write to
#' @param discharge list of discharge values
write_discharge = function(chain, nc_out, discharge) {

  # Chain


  # Discharge
  q = tryCatch(
    error = function(cond) grp.def.nc(nc_out, "q"),
    grp.inq.nc(nc_out, "q")$self
  )
  var.def.nc(q, name, "NC_DOUBLE", "nt")
  att.put.nc(q, name, "_FillValue", "NC_DOUBLE", FILL)
  discharge[is.nan(discharge)] = NA
  var.put.nc(q, name, discharge)
}

#' Write discharge and posteriors to NetCDF file.
#'
#' @param data named list of metadata
#' @param posteriors list of posterior "chains"
#' @param discharge list of discharge "chains"
#' @param out_dir string to output directory
#'
#' @export
write_output = function(data, posteriors, discharge, out_dir) {

  # Concatenate invalid times back into discharge
  discharge = lapply(discharge, concatenate_invalid, invalid_times=data$invalid_times)

  # File creation
  nc_file = paste(out_dir, paste0(data$reach_id, "_geobam.nc"), sep=.Platform$file.sep)
  nc_out = create.nc(nc_file, format="netcdf4")

  # Global attributes
  att.put.nc(nc_out, "NC_GLOBAL", "reach_id", "NC_INT64", data$reach_id)

  # Dimensions
  dim.def.nc(nc_out, "nt", length(data$nt))
  var.def.nc(nc_out, "nt", "NC_INT", "nt")
  att.put.nc(nc_out, "nt", "units", "NC_STRING", "time")
  var.put.nc(nc_out, "nt", data$nt)

  # Write data
  lapply(list(1,2,3), write_posteriors, nc_out=nc_out, posteriors=posteriors)

  # Discharge
  lapply(list(1,2,3), write_discharge, nc_out=nc_out, discharge=discharge)

  # Close file
  close.nc(nc_out)
}
