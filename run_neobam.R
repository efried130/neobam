#' Run neoBAM
#'
#' Execute neoBAM on input files and write to container output directory.
#'
#' Commandline arguments (optional):
#' reach_files
#'    - name of JSON file that contains associated reach file data

# Functions
source("/app/neobamdata/neobam/config.R")
source("/app/neobamdata/neobam/input.R")
source("/app/neobamdata/neobam/neobam_functions.R")
source("/app/neobamdata/neobam/output.R")
source("/app/neobamdata/neobam/process.R")

# Constants
IN_DIR = file.path("/mnt", "data", "input")
OUT_DIR = file.path("/mnt", "data", "output")
STAN_FILE = file.path("/app", "neobamdata", "neobam", "neobam_stan_engine.stan")

#' Identify reach and locate SWOT and SoS files.
#'
#' @param reaches_json name of json reach file
#'
#' @return list of swot file and sos file
get_reach_files = function(reaches_json){
  # Get reach identifier from array environment variable
  index = strtoi(Sys.getenv("AWS_BATCH_JOB_ARRAY_INDEX")) + 1
  json_data = rjson::fromJSON(file=file.path(IN_DIR, reaches_json))[[index]]
  return(list(reach_id = json_data$reach_id,
              swot_file = file.path(IN_DIR, "swot", json_data$swot),
              sos_file = file.path(IN_DIR, "sos", json_data$sos)))
}

#' Execute neoBAM
main = function() {

  # Identify reach files to process
  start = Sys.time()
  args = commandArgs(trailingOnly=TRUE)
  reaches_json = ifelse(identical(args, character(0)), "reaches.json", args[1])
  io_data = get_reach_files(reaches_json)

  # Get Input
  in_data = get_input(io_data$swot_file, io_data$sos_file, io_data$reach_id)

  # Process
  process_data = process_data(in_data, STAN_FILE)

  # Write Output
  out_data = list(reach_id = io_data$reach_id,
                  nt = in_data$swot_data$nt,
                  invalid_nodes = in_data$invalid_nodes,
                  invalid_times = in_data$invalid_times)
  write_output(out_data, process_data$posteriors, process_data$discharge, OUT_DIR)
  end = Sys.time()
  print(paste("Total execution time:", (end - start)))

}

main()
