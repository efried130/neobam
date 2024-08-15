#' Run neoBAM
#'
#' Execute neoBAM on input files and write to container output directory.
#'
#' Commandline arguments (optional):
#' reach_files
#'    - name of JSON file that contains associated reach file data

# Imports
library(reticulate)
library(optparse)

# Functions
# source("/app/neobam/input.R")
# source("/app/neobam/neobam_functions.R")
# source("/app/neobam/output.R")
# source("/app/neobam/process.R")

# source("neobam/input.R")
# source("neobam/neobam_functions.R")
# source("neobam/output.R")
# source("neobam/process.R")

source("/Users/tebaldi/Documents/workspace/confluence/workspace/neobam/neobam/input.R")
source("/Users/tebaldi/Documents/workspace/confluence/workspace/neobam/neobam/neobam_functions.R")
source("/Users/tebaldi/Documents/workspace/confluence/workspace/neobam/neobam/output.R")
source("/Users/tebaldi/Documents/workspace/confluence/workspace/neobam/neobam/process.R")

# Constants
# # IN_DIR = file.path("/nas/cee-water/cjgleason/SWOT_Q_UMASS/mnt",  "input")
# # OUT_DIR = file.path("/nas/cee-water/cjgleason/SWOT_Q_UMASS/mnt", "output")
# IN_DIR = file.path("/mnt/data/input")
# OUT_DIR = file.path("/mnt/data/output")
# STAN_FILE = file.path("/app", "neobam", "neobam_stan_engine.stan")
# # STAN_FILE = file.path( "neobam", "neobam_stan_engine.stan")
# PYTHON_EXE = "/usr/bin/python3"
# PYTHON_FILE = "/app/sos_read/sos_read.py"
# TMP_PATH = "/tmp"
IN_DIR = file.path("/Users/tebaldi/Documents/workspace/confluence/data/modules/neobam/input")
OUT_DIR = file.path("/Users/tebaldi/Documents/workspace/confluence/data/modules/neobam/output")
STAN_FILE = file.path("/Users/tebaldi/Documents/workspace/confluence/workspace/neobam", "neobam", "neobam_stan_engine.stan")
PYTHON_EXE = "/Users/tebaldi/Documents/workspace/environments/sos_read/bin/python3"
PYTHON_FILE = "/Users/tebaldi/Documents/workspace/confluence/workspace/neobam/sos_read/sos_read.py"
TMP_PATH = "/Users/tebaldi/Documents/workspace/confluence/data/modules/neobam/tmp"


#' Identify reach and locate SWOT and SoS files.
#'
#' @param reaches_json name of json reach file
#'
#' @return list of swot file and sos file
get_reach_files = function(reaches_json, index, bucket_key){
  # Get reach data from index
  json_data = rjson::fromJSON(file=file.path(IN_DIR, reaches_json))[[index]]

  if (bucket_key != "") {
    # Download the SoS file and reference the file path
    use_python(PYTHON_EXE)
    source_python(PYTHON_FILE)

    sos_filepath = file.path(TMP_PATH, json_data$sos)
    download_sos(bucket_key, sos_filepath)
    reach_list = list(reach_id = json_data$reach_id,
                      swot_file = file.path(IN_DIR, "swot", json_data$swot),
                      sos_file = sos_filepath)
  } else {
    reach_list = list(reach_id = json_data$reach_id,
                      swot_file = file.path(IN_DIR, "swot", json_data$swot),
                      sos_file = file.path(IN_DIR, "sos", json_data$sos))
  }

  return(reach_list)
}

#' Create output data structure for invalid observations
#'
#' @param nt number of time steps
#'
#' @return named list of discharge and posteriors
create_invalid_out = function(nt) {
  nt_vector = rep(NA_real_, nt)
  # base_discharge = list(nt_vector, nt_vector, nt_vector)
  base_discharge = list(nt_vector)
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

  option_list <- list(
    make_option(c("-i", "--index"), type = "integer", default = NULL, help = "Index to run on"),
    make_option(c("-b", "--bucket_key"), type = "character", default = "", help = "Bucket key to find the sos"),
    make_option(c("-r", "--reaches_json"), type = "character", default = "reaches.json", help = "Name of reaches.json")
  )

  opt_parser <- OptionParser(option_list = option_list)
  opts <- parse_args(opt_parser)
  bucket_key <- opts$bucket_key
  index <- opts$index + 1    # Add 1 to AWS 0-based index
  reaches_json <- opts$reaches_json
  print(paste("bucket_key: ", bucket_key))
  print(paste("index: ", index))
  print(paste("reaches_json: ", reaches_json))

  io_data = get_reach_files(reaches_json, index, bucket_key)
  print(paste("reach_id: ", io_data$reach_id))
  print(paste("swot_file: ", io_data$swot_file))
  print(paste("sos_file: ", io_data$sos_file))

  # Get Input
  in_data = get_input(io_data$swot_file, io_data$sos_file, io_data$reach_id)

  print('processing data...')
  # Process
  if (in_data$valid != FALSE) {
    print('data was valid...')
    neobam_output = process_data(in_data, STAN_FILE)
    out_data = list(reach_id = io_data$reach_id,
                    nt = in_data$swot_data$nt,
                    invalid_nodes = in_data$invalid_nodes,
                    invalid_times = in_data$invalid_times,
                    node_ids = in_data$node_ids)
  } else {

    neobam_output = create_invalid_out(length(in_data$nt))
    out_data = list(reach_id = io_data$reach_id,
                    nt = in_data$nt,
                    invalid_nodes = vector(mode = "list"),
                    invalid_times = vector(mode = "list"),
                    node_ids = in_data$node_ids)
  }



  # Write output
  # print(neobam_output)
  write_output(out_data, neobam_output$posteriors, neobam_output$posterior_Q, neobam_output$posterior_Q_sd, OUT_DIR, in_data$valid, in_data$obs_times)
  end = Sys.time()
  print(paste("Total execution time for reach", io_data$reach_id, ":", (end - start), "seconds."))

    return(list(c(neobam_output,out_data)))
}

neobam_output=main()
# print(neobam_output)
