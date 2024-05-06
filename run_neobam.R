#' Run neoBAM
#'
#' Execute neoBAM on input files and write to container output directory.
#'
#' Commandline arguments (optional):
#' reach_files
#'    - name of JSON file that contains associated reach file data


# Functions
source("/app/neobam/input.R")
source("/app/neobam/neobam_functions.R")
source("/app/neobam/output.R")
source("/app/neobam/process.R")

# source("neobam/input.R")
# source("neobam/neobam_functions.R")
# source("neobam/output.R")
# source("neobam/process.R")

# Constants
# IN_DIR = file.path("/nas/cee-water/cjgleason/SWOT_Q_UMASS/mnt",  "input")
# OUT_DIR = file.path("/nas/cee-water/cjgleason/SWOT_Q_UMASS/mnt", "output")
IN_DIR = file.path("/mnt/data/input")
OUT_DIR = file.path("/mnt/data/output")
STAN_FILE = file.path("/app", "neobam", "neobam_stan_engine.stan")
# STAN_FILE = file.path( "neobam", "neobam_stan_engine.stan")


#' Identify reach and locate SWOT and SoS files.
#'
#' @param reaches_json name of json reach file
#'
#' @return list of swot file and sos file
get_reach_files = function(reaches_json){
  # Get reach identifier from array environment variable
  index = strtoi(Sys.getenv("AWS_BATCH_JOB_ARRAY_INDEX")) + 1
  json_data = rjson::fromJSON(file=file.path(reaches_json))[[index]]
  return(list(reach_id = json_data$reach_id,
              swot_file = file.path(IN_DIR, "swot", json_data$swot),
              sos_file = file.path(IN_DIR, "sos", json_data$sos)))
}

#' Create output data structure for invalid observations
#'
#' @param nt number of time steps
#'
#' @return named list of discharge and posteriors
create_invalid_out = function(nt) {
  nt_vector = rep(NA_real_, nt)
  base_discharge = list(nt_vector, nt_vector, nt_vector)
  base_posteriors = list(
    r = list(mean=NA_real_, sd=NA_real_),
    logn = list(mean=NA_real_, sd=NA_real_),
    logWb = list(mean=NA_real_, sd=NA_real_),
    logDb = list(mean=NA_real_, sd=NA_real_)
  )
  return(list(discharge=base_discharge, posteriors=list(base_posteriors, base_posteriors, base_posteriors)))
}

#' Execute neoBAM
main = function() {

  # Identify reach files to process
  start = Sys.time()
  args = commandArgs(trailingOnly=TRUE)
  reaches_json = ifelse(identical(args, character(0)), "/mnt/data/input/reaches.json", args[1])
      # reaches_json= '/nas/cee-water/cjgleason/SWOT_Q_UMASS/mnt/input/reaches.json'
  io_data = get_reach_files(reaches_json)
    
  # Get Input
  in_data = get_input(io_data$swot_file, io_data$sos_file, io_data$reach_id)

   
  # Process
  if (in_data$valid != FALSE) {
    neobam_output = process_data(in_data, STAN_FILE)
    out_data = list(reach_id = io_data$reach_id,
                    nt = in_data$swot_data$nt,
                    invalid_nodes = in_data$invalid_nodes,
                    invalid_times = in_data$invalid_times)
  } else {
   
    neobam_output = create_invalid_out(length(in_data$nt))
    out_data = list(reach_id = io_data$reach_id,
                    nt = in_data$nt,
                    invalid_nodes = vector(mode = "list"),
                    invalid_times = vector(mode = "list"))
  }
    
   

  # Write output
  # write_output(out_data, neobam_output$posteriors, neobam_output$posterior_Q, OUT_DIR)
  end = Sys.time()
  print(paste("Total execution time for reach", io_data$reach_id, ":", (end - start), "seconds."))
    
    return(list(c(neobam_output,out_data)))
}

neobam_output=main()
