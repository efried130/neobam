


#' Performs the following checks:
#' - types:
#'     - logQ_hat is numeric vector
#'     - everything else matrix
#' - dimensions:
#'     - all matrices have same dims
#'     - logQ_hat has length equal to ncol of matrices
#'
#' @param datalist A list of BAM data inputs
bam_check_args <- function(datalist) {

  dA <- datalist$dA
  logQ_hat <- datalist$logQ_hat
  matlist <- datalist[names(datalist) != "logQ_hat"]

  # Check types
  if (!(is(logQ_hat, "numeric") && is(logQ_hat, "vector")))
    stop("Qhat must be a numeric vector.\n")
  if (!all(vapply(matlist, is, logical(1), "matrix")))
    stop("All data must be a supplied as a matrix.\n")

  # Check dims
  nr <- nrow(matlist[[1]])
  nc <- ncol(matlist[[1]])
  if (!(all(vapply(matlist, nrow, 0L) == nr) &&
        all(vapply(matlist, ncol, 0L) == nc)))
    stop("All data must have same dimensions.\n")
  if (!length(logQ_hat) == nc)
    logQ_hat <- rep(logQ_hat, length.out = nc)

  out <- c(matlist, list(logQ_hat = logQ_hat))

  out
}

#' Add missing-data inputs to data list
#'
#' Binary matrices indicating where data are/aren't missing are
#' added to the data list. This is required in order to run
#' ragged-array data structures in the stanfile.
#'
#' Previously this function omitted any times with missing data,
#' but now that ragged arrays are accommodated in the stanfile the
#' operations are entirely different.
#'
#' @param datalist a list of BAM inputs
#' @param variant which geoBAM is being run
#' @importFrom stats median
bam_check_nas <- function(datalist, variant) {

  mats <- vapply(datalist, is.matrix, logical(1))
  nonas <- lapply(datalist[mats], function(x) !is.na(x))

  # AMHG has-data matrix (needs slopes and widths)
  hasdat_s <- (!is.na(datalist[["Sobs"]])) * 1
  hasdat_w <- (!is.na(datalist[["Wobs"]])) * 1
  hasdat_amhg <- hasdat_w * hasdat_s

  # Replace NA's with zeros so Stan will accept the data
  datalist[["Wobs"]][!hasdat_amhg] <- 0
  datalist[["Sobs"]][!hasdat_amhg] <- 0

  # Manning has-data matrix (only nonzero if all Manning obs present)
  if (identical(setdiff(c("Wobs", "Sobs", "dAobs"),
                        names(datalist[mats])),
                character(0))) {
    #hasdat_s <- (!is.na(datalist[["Sobs"]])) * 1
    hasdat_a <- (!is.na(datalist[["dAobs"]])) * 1

    hasdat_man <- hasdat_amhg * hasdat_a

    # Replace NA's with zeros so Stan will accept the data
    datalist[["dAobs"]][!hasdat_man] <- 0

  } else {
    hasdat_man <- matrix(0, nrow = nrow(hasdat_amhg), ncol = ncol(hasdat_amhg))
  }

  if (!is.null(datalist[["dAobs"]])) {
    dA_shift <- apply(datalist[["dAobs"]], 1, function(x) median(x) - min(x))
  } else {
    dA_shift <- rep(0, nrow(datalist[["Wobs"]]))
  }

  newbits <- list(
    hasdat_man = hasdat_man,
    hasdat_amhg = hasdat_amhg,
    ntot_man = sum(hasdat_man),
    ntot_amhg = sum(hasdat_amhg),
    dA_shift = dA_shift
  )

  out <- c(datalist, newbits)
  out
}

#' Establish prior hyperparameters for BAM estimation
#'
#' Produces a bampriors object that can be passed to bam_estimate function
#'
#' @useDynLib geoBAMr, .registration = TRUE
#' @param bamdata An object of class bamdata, as returned by \code{bam_data}
#' @param variant Which geoBAM variant to use. Options are "manning_amhg" (default),
#'   "manning", or "amhg".
#' @param classification Which classification framework to use. Options are 'expert' (default),
#'   or 'unsupervised'.
#' @param ... Optional manually set parameters. Unquoted expressions are allowed,
#'   e.g. \code{logQ_sd = cv2sigma(0.8)}. Additionally, any variables present in
#'   \code{bamdata} may be referenced, e.g. \code{lowerbound_logQ = log(mean(Wobs)) + log(5)}
#' @export

bam_priors <- function(bamdata,
                       variant = c("manning_amhg", "manning", "amhg"),
                       classification = c('expert', 'unsupervised'),
                       ...) {
  variant <- match.arg(variant)
  classification <- match.arg(classification)
  if (variant != "amhg" && !attr(bamdata, "manning_ready"))
    stop("bamdata must have dA data for non-amhg variants.")

  force(bamdata)
  paramset <- bam_settings("paramnames")

  myparams0 <- rlang::quos(..., .named = TRUE)
  myparams <- do.call(settings::clone_and_merge,
                      args = c(list(options = bam_settings), myparams0))

  quoparams <- myparams()[-1] # first one is parameter set
  params <- lapply(quoparams, rlang::eval_tidy, data = bamdata)

  if (classification == 'unsupervised'){ #recalculate parameter set if classification switched to unsupervised
    paramset <- bam_settings_unsupervised("paramnames")

    myparams0 <- rlang::quos(..., .named = TRUE)
    myparams <- do.call(settings::clone_and_merge,
                        args = c(list(options = bam_settings_unsupervised), myparams0))

    quoparams <- myparams()[-1] # first one is parameter set
    params <- lapply(quoparams, rlang::eval_tidy, data = bamdata)
  }

  if (!length(params[["logQ_sd"]]) == bamdata$nt)
    params$logQ_sd <- rep(params$logQ_sd, length.out = bamdata$nt)

  if (!identical(dim(params[["sigma_man"]]),
                 as.integer(c(bamdata$nx, bamdata$nt)))) {
    params$sigma_man <- matrix(rep(params$sigma_man,
                                   length.out = bamdata$nt * bamdata$nx),
                               nrow = bamdata$nx)
  }

  if (!identical(dim(params[["sigma_amhg"]]),
                 as.integer(c(bamdata$nx, bamdata$nt)))) {
    params$sigma_amhg <- matrix(rep(params$sigma_amhg,
                                    length.out = bamdata$nt * bamdata$nx),
                                nrow = bamdata$nx)
  }

  #just priors the user wants to see
  user_paramset <- c('lowerbound_logQ', 'upperbound_logQ', 'lowerbound_logWc', 'upperbound_logWc', 'lowerbound_logQc', 'upperbound_logQc',
                     'logWc_hat',"logQc_hat",
                     'logQ_sd','logWc_sd','logQc_sd',
                     "Werr_sd", "Serr_sd", "dAerr_sd",
                     "sigma_man", "sigma_amhg")
  user_paramset <- params[user_paramset]

  #total priors needed to run geoBAM
  geoBAM_paramset <- c( "lowerbound_A0", "upperbound_A0", "lowerbound_logn", "upperbound_logn","lowerbound_b","upperbound_b", "lowerbound_logWb", 'upperbound_logWb', 'lowerbound_logDb','upperbound_logDb', 'lowerbound_logr', 'upperbound_logr',
                        "logA0_hat", "logn_hat","b_hat",  'logWb_hat', 'logDb_hat', 'logr_hat',
                       "logA0_sd", "logn_sd", "b_sd",'logWb_sd', 'logDb_sd', 'logr_sd')
  geoBAMparams <- params[geoBAM_paramset]

  riverType <- params[["River_Type"]]

  out <- list( 'River_Type'=riverType, 'river_type_priors'=geoBAMparams, 'other_priors'=user_paramset)
  out <- structure(out,
                   class = c("bampriors"))
  out
}

compose_bam_inputs <- function(bamdata, priors = bam_priors(bamdata)) {

  inps <- c(bamdata, priors)

  out <- inps
  out

}


#' Take a random sample of a bamdata object's cross-sections.
#'
#' @param bamdata a bamdata object, as returned by \code{bam_data()}
#' @param n Number of cross-sections to
#' @param seed option RNG seed, for reproducibility.
#' @importFrom methods is
#' @export
sample_xs <- function(bamdata, n, seed = NULL) {

  stopifnot(is(bamdata, "bamdata"))

  if (n >= bamdata$nx)
    return(bamdata)

  if (!is.null(seed))
    set.seed(seed)
  keepxs <- sort(sample(1:bamdata$nx, size = n, replace = FALSE))

  bamdata$nx <- n
  bamdata$Wobs <- bamdata$Wobs[keepxs, ]

  if (!is.null(bamdata$Sobs)) {
    bamdata$Sobs <- bamdata$Sobs[keepxs, ]
    bamdata$dAobs <- bamdata$dAobs[keepxs, ]
  }

  bamdata
}



#' Calculate lognormal moments based on truncated normal parameters
#'
#' Used to put measurement errors into original log-normal parameterization.
#'
#' @param obs A numeric vector of observations
#' @param err_sigma Standard deviation of measurement error
#' @param a zero-reference for method of moments.
#' @importFrom stats dnorm pnorm
ln_moms <- function(obs, err_sigma, a = 0) {
  alpha <- (a - obs) / err_sigma
  Z <- 1 - pnorm(alpha)

  mean <- obs + (dnorm(alpha)) / Z * err_sigma
  sdquan <- 1 + (alpha * dnorm(alpha)) / Z -
    (dnorm(alpha) / Z)^2
  sd <- err_sigma * sqrt(sdquan)

  out <- list(mean = mean, sd = sd)
  out
}

#' Calculate lognormal sigma parameter based on truncated normal parameters
#'
#' Used to put measurement errors into original log-normal parameterization.
#'
#' @param obs A numeric vector of observations
#' @param err_sigma Standard deviation of measurement error
#' @param a zero-reference for method of moments.
ln_sigsq <- function(obs, err_sigma, a = 0) {
  moms <- ln_moms(obs = obs, err_sigma = err_sigma, a = a)
  mn <- unname(moms[["mean"]])
  sd <- unname(moms[["sd"]])
  mu <- 2 * log(mn) - 0.5 * log(sd^2 + mn^2)
  sigsq <- 2 * log(mn) - 2 * mu

  sigsq
}
                      
                      
                     

#Bounds on Manning's n are manually set to log(0.01) and log(0.05).
  #Functions are written to revert back to standard approach to defining priors, currently commented out.
#Mark's defualt BAM priors are also commented out.



#' Options manager for geoBAMr defaults
#'
#' @param ... (Optional) named settings to query or set.
#' @param .__defaults See \code{?settings::option_manager}
#' @param .__reset See \code{?settings::option_manager}
#' @export
bam_settings <- settings::options_manager(
  paramnames = c("lowerbound_logQ", "upperbound_logQ", "lowerbound_A0",
                 "upperbound_A0", "lowerbound_logn", "upperbound_logn", "lowerbound_logQc",
                 "upperbound_logQc", "lowerbound_logWc", "upperbound_logWc", "lowerbound_b",
                 "upperbound_b", "lowerbound_logWb", 'upperbound_logWb', 'lowerbound_logDb',
                 'upperbound_logDb', 'lowerbound_logr', 'upperbound_logr',
                 "sigma_man", "sigma_amhg",
                 "logQc_hat", "logWc_hat", "b_hat", "logA0_hat", "logn_hat", 'logWb_hat', 'logDb_hat', 'logr_hat',
                 "logQ_sd", "logQc_sd", "logWc_sd", "b_sd", "logA0_sd", "logn_sd", 'logWb_sd', 'logDb_sd', 'logr_sd',
                 "Werr_sd", "Serr_sd", "dAerr_sd",
                 'River_Type'),

  # Bounds on parameters
  lowerbound_logQ = rlang::quo(estimate_lowerbound_logQ(Wobs)),
  upperbound_logQ = rlang::quo(estimate_upperbound_logQ(Wobs)),

  lowerbound_A0 = rlang::quo(estimate_lowerboundA0(Wobs)), #0.72,
  upperbound_A0 = rlang::quo(estimate_upperboundA0(Wobs)), #114500,
  lowerbound_logn = log(0.01), #rlang::quo(estimate_lowerboundlogn(Wobs)), #-4.60517,
  upperbound_logn = log(0.05), #rlang::quo(estimate_upperboundlogn(Wobs)), #-2.995732,

  lowerbound_logQc = 0.01,
  upperbound_logQc = 10,
  lowerbound_logWc = 1,
  upperbound_logWc = 8, # 3 km
  lowerbound_b = rlang::quo(estimate_lowerboundb(Wobs)), #0.000182,
  upperbound_b = rlang::quo(estimate_upperboundb(Wobs)), #0.773758,

  lowerbound_logDb = rlang::quo(estimate_lowerboundlogDb(Wobs)), #-3.02002,
  upperbound_logDb = rlang::quo(estimate_upperboundlogDb(Wobs)), #3.309359,
  lowerbound_logWb = rlang::quo(estimate_lowerboundlogWb(Wobs)),#-0.12273,
  upperbound_logWb = rlang::quo(estimate_upperboundlogWb(Wobs)), #7.006786,
  lowerbound_logr = rlang::quo(estimate_lowerboundlogr(Wobs)), #-2.58047,
  upperbound_logr = rlang::quo(estimate_upperboundlogr(Wobs)), #8.03772,


  # *Known* likelihood parameters
  sigma_man = 0.25,
  sigma_amhg = 0.22, # UPDATE THIS FROM CAITLINE'S WORK


  # Hyperparameters
  logQc_hat = rlang::quo(mean(logQ_hat, na.rm=TRUE)),
  logWc_hat = rlang::quo(estimate_Wc(Wobs)),
  b_hat = rlang::quo(estimate_b(Wobs)),
  logA0_hat = rlang::quo(estimate_logA0(Wobs)),
  logn_hat = rlang::quo(estimate_logn(Sobs, Wobs)),
  logWb_hat = rlang::quo(estimate_logWb(Wobs)),
  logDb_hat = rlang::quo(estimate_logDb(Wobs)),
  logr_hat = rlang::quo(estimate_logr(Wobs)),

  logQ_sd = sqrt(log(1^2 + 1)), # CV of Q equals 1
  logQc_sd = sqrt(log(1^2 + 1)), # CV of Q equals 1; UPDATE THIS
  logWc_sd = sqrt(log(0.01)^2 + 1),

  #set from my model outputs
  b_sd = rlang::quo(estimate_bSD(Wobs)), #0.068077044,
  logA0_sd = rlang::quo(estimate_A0SD(Wobs)), #0.58987527,
  logn_sd = rlang::quo(estimate_lognSD(Wobs)), #0.761673112,
  logWb_sd = rlang::quo(estimate_logWbSD(Wobs)), #0.137381044,
  logDb_sd = rlang::quo(estimate_logDbSD(Wobs)), #0.576212733,
  logr_sd = rlang::quo(estimate_logrSD(Wobs)), #0.67332688,

  # Observation errors.
  Werr_sd = 10,
  Serr_sd = 1e-5,
  dAerr_sd = 10,

  #River type
  River_Type=rlang::quo(apply(Wobs, 1, classify_func))
)


#' Options manager for geoBAMr defaults using unsupervised classification
#'
#' @param ... (Optional) named settings to query or set.
#' @param .__defaults See \code{?settings::option_manager}
#' @param .__reset See \code{?settings::option_manager}
#' @export
bam_settings_unsupervised <- settings::options_manager(
  paramnames = c("lowerbound_logQ", "upperbound_logQ", "lowerbound_A0",
                 "upperbound_A0", "lowerbound_logn", "upperbound_logn", "lowerbound_logQc",
                 "upperbound_logQc", "lowerbound_logWc", "upperbound_logWc", "lowerbound_b",
                 "upperbound_b", "lowerbound_logWb", 'upperbound_logWb', 'lowerbound_logDb',
                 'upperbound_logDb', 'lowerbound_logr', 'upperbound_logr',
                 "sigma_man", "sigma_amhg",
                 "logQc_hat", "logWc_hat", "b_hat", "logA0_hat", "logn_hat", 'logWb_hat', 'logDb_hat', 'logr_hat',
                 "logQ_sd", "logQc_sd", "logWc_sd", "b_sd", "logA0_sd", "logn_sd", 'logWb_sd', 'logDb_sd', 'logr_sd',
                 "Werr_sd", "Serr_sd", "dAerr_sd",
                 'River_Type'),

  # Bounds on parameters
  lowerbound_logQ = rlang::quo(estimate_lowerbound_logQ(Wobs)),
  upperbound_logQ = rlang::quo(estimate_upperbound_logQ(Wobs)),

  lowerbound_A0 = rlang::quo(estimate_lowerboundA0_unsupervised(Wobs)), #0.72,
  upperbound_A0 = rlang::quo(estimate_upperboundA0_unsupervised(Wobs)), #114500,
  lowerbound_logn = log(0.01), #rlang::quo(estimate_lowerboundlogn_unsupervised(Wobs)), #-4.60517,
  upperbound_logn = log(0.05), #rlang::quo(estimate_upperboundlogn_unsupervised(Wobs)), #-2.995732,

  lowerbound_logQc = 0.01,
  upperbound_logQc = 10,
  lowerbound_logWc = 1,
  upperbound_logWc = 8, # 3 km
  lowerbound_b = rlang::quo(estimate_lowerboundb_unsupervised(Wobs)), #0.000182,
  upperbound_b = rlang::quo(estimate_upperboundb_unsupervised(Wobs)), #0.773758,

  lowerbound_logDb = rlang::quo(estimate_lowerboundlogDb_unsupervised(Wobs)), #-3.02002,
  upperbound_logDb = rlang::quo(estimate_upperboundlogDb_unsupervised(Wobs)), #3.309359,
  lowerbound_logWb = rlang::quo(estimate_lowerboundlogWb_unsupervised(Wobs)),#-0.12273,
  upperbound_logWb = rlang::quo(estimate_upperboundlogWb_unsupervised(Wobs)), #7.006786,
  lowerbound_logr = rlang::quo(estimate_lowerboundlogr_unsupervised(Wobs)), #-2.58047,
  upperbound_logr = rlang::quo(estimate_upperboundlogr_unsupervised(Wobs)), #8.03772,


  # *Known* likelihood parameters
  sigma_man = 0.25,
  sigma_amhg = 0.22, # UPDATE THIS FROM CAITLINE'S WORK


  # Hyperparameters
  logQc_hat = rlang::quo(mean(logQ_hat, na.rm=TRUE)),
  logWc_hat = rlang::quo(estimate_Wc(Wobs)),
  b_hat = rlang::quo(estimate_b_unsupervised(Wobs)),
  logA0_hat = rlang::quo(estimate_logA0_unsupervised(Wobs)),
  logn_hat = rlang::quo(estimate_logn_unsupervised(Sobs, Wobs)),
  logWb_hat = rlang::quo(estimate_logWb_unsupervised(Wobs)),
  logDb_hat = rlang::quo(estimate_logDb_unsupervised(Wobs)),
  logr_hat = rlang::quo(estimate_logr_unsupervised(Wobs)),

  logQ_sd = sqrt(log(1^2 + 1)), # CV of Q equals 1
  logQc_sd = sqrt(log(1^2 + 1)), # CV of Q equals 1; UPDATE THIS
  logWc_sd = sqrt(log(0.01^2 + 1)), #CV of Wc equals 0.01

  #set from my model outputs
  b_sd = rlang::quo(estimate_bSD_unsupervised(Wobs)), #0.068077044,
  logA0_sd = rlang::quo(estimate_A0SD_unsupervised(Wobs)), #0.58987527,
  logn_sd = rlang::quo(estimate_lognSD_unsupervised(Wobs)), #0.761673112,
  logWb_sd = rlang::quo(estimate_logWbSD_unsupervised(Wobs)), #0.137381044,
  logDb_sd = rlang::quo(estimate_logDbSD_unsupervised(Wobs)), #0.576212733,
  logr_sd = rlang::quo(estimate_logrSD_unsupervised(Wobs)), #0.67332688,

  # Observation errors.
  Werr_sd = 10,
  Serr_sd = 1e-5,
  dAerr_sd = 10,

  #river type
  River_Type=rlang::quo(apply(Wobs, 1, classify_func_unsupervised))
)

                      # Using advice from read-and-delete-me:
# "Be sure to add useDynLib(mypackage, .registration = TRUE) to the NAMESPACE file
# which you can do by placing the line   #' @useDynLib rstanarm, .registration = TRUE
# in one of your .R files
# see rstanarm's 'rstanarm-package.R' file"

#' Preprocess data for BAM estimation
#'
#' Produces a bamdata object that can be passed to bam_estimate function
#'
#' @useDynLib geoBAMr, .registration = TRUE
#' @param w Matrix (or data frame) of widths: time as columns, space as rows
#' @param s Matrix (or data frame) of slopes: time as columns, space as rows
#' @param dA Matrix of area above base area: time as columns, space as rows
#' @param Qhat Vector of Q estimates. Needed to create prior on Q.
#' @param variant Which geoBAM variant to use. Options are "manning_amhg" (default),
#'   "manning", or "amhg".
#' @param max_xs Maximum number of cross-sections to allow in data. Used to reduce
#'   sampling time. Defaults to 30.
#' @param seed RNG seed to use for sampling cross-sections, if nx > max_xs.
#' @export


bam_data <- function(h,w,dA,Qhat,s){
 
  datalist <- list(Hobs = h,
                Sobs = s,
                dAobs = dA,
                   Wobs=w,
                logQ_hat = log(Qhat))

  nx <- nrow(datalist$Hobs)
  nt <- ncol(datalist$WHobs)

  out <- structure(c(datalist,
                        nx = nx,
                        nt = nt),
                   class = c("bamdata"))
  out
}

#' Performs the following checks:
#' - types:
#'     - logQ_hat is numeric vector
#'     - everything else matrix
#' - dimensions:
#'     - all matrices have same dims
#'     - logQ_hat has length equal to ncol of matrices
#'
#' @param datalist A list of BAM data inputs
bam_check_args <- function(datalist) {

  dA <- datalist$dA
  logQ_hat <- datalist$logQ_hat
  matlist <- datalist[names(datalist) != "logQ_hat"]

  # Check types
  if (!(is(logQ_hat, "numeric") && is(logQ_hat, "vector")))
    stop("Qhat must be a numeric vector.\n")
  if (!all(vapply(matlist, is, logical(1), "matrix")))
    stop("All data must be a supplied as a matrix.\n")

  # Check dims
  nr <- nrow(matlist[[1]])
  nc <- ncol(matlist[[1]])
  if (!(all(vapply(matlist, nrow, 0L) == nr) &&
        all(vapply(matlist, ncol, 0L) == nc)))
    stop("All data must have same dimensions.\n")
  if (!length(logQ_hat) == nc)
    logQ_hat <- rep(logQ_hat, length.out = nc)

  out <- c(matlist, list(logQ_hat = logQ_hat))

  out
}

#' Add missing-data inputs to data list
#'
#' Binary matrices indicating where data are/aren't missing are
#' added to the data list. This is required in order to run
#' ragged-array data structures in the stanfile.
#'
#' Previously this function omitted any times with missing data,
#' but now that ragged arrays are accommodated in the stanfile the
#' operations are entirely different.
#'
#' @param datalist a list of BAM inputs
#' @param variant which geoBAM is being run
#' @importFrom stats median
bam_check_nas <- function(datalist, variant) {

  mats <- vapply(datalist, is.matrix, logical(1))
  nonas <- lapply(datalist[mats], function(x) !is.na(x))

  # AMHG has-data matrix (needs slopes and widths)
  hasdat_s <- (!is.na(datalist[["Sobs"]])) * 1
  hasdat_w <- (!is.na(datalist[["Wobs"]])) * 1
  hasdat_amhg <- hasdat_w * hasdat_s

  # Replace NA's with zeros so Stan will accept the data
  datalist[["Wobs"]][!hasdat_amhg] <- 0
  datalist[["Sobs"]][!hasdat_amhg] <- 0

  # Manning has-data matrix (only nonzero if all Manning obs present)
  if (identical(setdiff(c("Wobs", "Sobs", "dAobs"),
                        names(datalist[mats])),
                character(0))) {
    #hasdat_s <- (!is.na(datalist[["Sobs"]])) * 1
    hasdat_a <- (!is.na(datalist[["dAobs"]])) * 1

    hasdat_man <- hasdat_amhg * hasdat_a

    # Replace NA's with zeros so Stan will accept the data
    datalist[["dAobs"]][!hasdat_man] <- 0

  } else {
    hasdat_man <- matrix(0, nrow = nrow(hasdat_amhg), ncol = ncol(hasdat_amhg))
  }

  if (!is.null(datalist[["dAobs"]])) {
    dA_shift <- apply(datalist[["dAobs"]], 1, function(x) median(x) - min(x))
  } else {
    dA_shift <- rep(0, nrow(datalist[["Wobs"]]))
  }

  newbits <- list(
    hasdat_man = hasdat_man,
    hasdat_amhg = hasdat_amhg,
    ntot_man = sum(hasdat_man),
    ntot_amhg = sum(hasdat_amhg),
    dA_shift = dA_shift
  )

  out <- c(datalist, newbits)
  out
}

#' Establish prior hyperparameters for BAM estimation
#'
#' Produces a bampriors object that can be passed to bam_estimate function
#'
#' @useDynLib geoBAMr, .registration = TRUE
#' @param bamdata An object of class bamdata, as returned by \code{bam_data}
#' @param variant Which geoBAM variant to use. Options are "manning_amhg" (default),
#'   "manning", or "amhg".
#' @param classification Which classification framework to use. Options are 'expert' (default),
#'   or 'unsupervised'.
#' @param ... Optional manually set parameters. Unquoted expressions are allowed,
#'   e.g. \code{logQ_sd = cv2sigma(0.8)}. Additionally, any variables present in
#'   \code{bamdata} may be referenced, e.g. \code{lowerbound_logQ = log(mean(Wobs)) + log(5)}
#' @export

bam_priors <- function(bamdata,
                       variant = c("manning_amhg", "manning", "amhg"),
                       classification = c('expert', 'unsupervised'),
                       ...) {
    
   # library(settings,lib.loc='/nas/cee-water/cjgleason/r-lib/')
   # library(rlang,lib.loc='/nas/cee-water/cjgleason/r-lib/')
    
  variant <- match.arg(variant)
  classification <- match.arg(classification)
  if (variant != "amhg" && !attr(bamdata, "manning_ready"))
    stop("bamdata must have dA data for non-amhg variants.")

  force(bamdata)
  paramset <- bam_settings("paramnames")

  myparams0 <- rlang::quos(..., .named = TRUE)
  myparams <- do.call(settings::clone_and_merge,
                      args = c(list(options = bam_settings), myparams0))

  quoparams <- myparams()[-1] # first one is parameter set
  params <- lapply(quoparams, rlang::eval_tidy, data = bamdata)

  if (classification == 'unsupervised'){ #recalculate parameter set if classification switched to unsupervised
    paramset <- bam_settings_unsupervised("paramnames")

    myparams0 <- rlang::quos(..., .named = TRUE)
    myparams <- do.call(settings::clone_and_merge,
                        args = c(list(options = bam_settings_unsupervised), myparams0))

    quoparams <- myparams()[-1] # first one is parameter set
    params <- lapply(quoparams, rlang::eval_tidy, data = bamdata)
  }



  if (!identical(dim(params[["sigma_man"]]),
                 as.integer(c(bamdata$nx, bamdata$nt)))) {
    params$sigma_man <- matrix(rep(params$sigma_man,
                                   length.out = bamdata$nt * bamdata$nx),
                               nrow = bamdata$nx)
  }

  if (!identical(dim(params[["sigma_amhg"]]),
                 as.integer(c(bamdata$nx, bamdata$nt)))) {
    params$sigma_amhg <- matrix(rep(params$sigma_amhg,
                                    length.out = bamdata$nt * bamdata$nx),
                                nrow = bamdata$nx)
  }

  #just priors the user wants to see
  user_paramset <- c('lowerbound_logQ', 'upperbound_logQ', 'lowerbound_logWc', 'upperbound_logWc', 'lowerbound_logQc', 'upperbound_logQc',
                     'logWc_hat',"logQc_hat",
                     'logQ_sd','logWc_sd','logQc_sd',
                     "Werr_sd", "Serr_sd", "dAerr_sd",
                     "sigma_man", "sigma_amhg")
  user_paramset <- params[user_paramset]

  #total priors needed to run geoBAM
  geoBAM_paramset <- c( "lowerbound_A0", "upperbound_A0", "lowerbound_logn", "upperbound_logn","lowerbound_b","upperbound_b", "lowerbound_logWb", 'upperbound_logWb', 'lowerbound_logDb','upperbound_logDb', 'lowerbound_logr', 'upperbound_logr',
                        "logA0_hat", "logn_hat","b_hat",  'logWb_hat', 'logDb_hat', 'logr_hat',
                       "logA0_sd", "logn_sd", "b_sd",'logWb_sd', 'logDb_sd', 'logr_sd')
  geoBAMparams <- params[geoBAM_paramset]

  riverType <- params[["River_Type"]]

  out <- list( 'River_Type'=riverType, 'river_type_priors'=geoBAMparams, 'other_priors'=user_paramset)
  out <- structure(out,
                   class = c("bampriors"))
  out
}

compose_bam_inputs <- function(bamdata, priors = bam_priors(bamdata)) {

  inps <- c(bamdata, priors)

  out <- inps
  out

}


#' Take a random sample of a bamdata object's cross-sections.
#'
#' @param bamdata a bamdata object, as returned by \code{bam_data()}
#' @param n Number of cross-sections to
#' @param seed option RNG seed, for reproducibility.
#' @importFrom methods is
#' @export
sample_xs <- function(bamdata, n, seed = NULL) {

  stopifnot(is(bamdata, "bamdata"))

  if (n >= bamdata$nx)
    return(bamdata)

  if (!is.null(seed))
    set.seed(seed)
  keepxs <- sort(sample(1:bamdata$nx, size = n, replace = FALSE))

  bamdata$nx <- n
  bamdata$Wobs <- bamdata$Wobs[keepxs, ]

  if (!is.null(bamdata$Sobs)) {
    bamdata$Sobs <- bamdata$Sobs[keepxs, ]
    bamdata$dAobs <- bamdata$dAobs[keepxs, ]
  }

  bamdata
}



#' Calculate lognormal moments based on truncated normal parameters
#'
#' Used to put measurement errors into original log-normal parameterization.
#'
#' @param obs A numeric vector of observations
#' @param err_sigma Standard deviation of measurement error
#' @param a zero-reference for method of moments.
#' @importFrom stats dnorm pnorm
ln_moms <- function(obs, err_sigma, a = 0) {
  alpha <- (a - obs) / err_sigma
  Z <- 1 - pnorm(alpha)

  mean <- obs + (dnorm(alpha)) / Z * err_sigma
  sdquan <- 1 + (alpha * dnorm(alpha)) / Z -
    (dnorm(alpha) / Z)^2
  sd <- err_sigma * sqrt(sdquan)

  out <- list(mean = mean, sd = sd)
  out
}

#' Calculate lognormal sigma parameter based on truncated normal parameters
#'
#' Used to put measurement errors into original log-normal parameterization.
#'
#' @param obs A numeric vector of observations
#' @param err_sigma Standard deviation of measurement error
#' @param a zero-reference for method of moments.
ln_sigsq <- function(obs, err_sigma, a = 0) {
  moms <- ln_moms(obs = obs, err_sigma = err_sigma, a = a)
  mn <- unname(moms[["mean"]])
  sd <- unname(moms[["sd"]])
  mu <- 2 * log(mn) - 0.5 * log(sd^2 + mn^2)
  sigsq <- 2 * log(mn) - 2 * mu

  sigsq
}
                      
                      #Prior calculation for non geoBAM priors-------------------------------------------------------------------------------
#' Estimate AHG logWc_hat using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_Wc <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  logWc_hat <- mean(log(Wobs), na.rm=TRUE)
}

#' Estimate AHG lowerbound_logQ using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerbound_logQ <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lowerbound_logQ <- maxmin(log(Wobs)) + log(0.5) + log(0.5)
}

#' Estimate AHG upperbound_logQ using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperbound_logQ <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  upperbound_logQ <- minmax(log(Wobs)) + log(40) + log(5)
}

# Prior calculation using expert classification framework------------------------------------------------------------------
#class 17 are 'big' rivers

#' Estimate AHG b exponent using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_b <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwsd <- apply(log(Wobs), 1, function(x) sd(x, na.rm = TRUE))

  #expert classification
  temp <- c(0.227633914,
            0.24121913,
            0.204897294,
            0.204451567,
            0.192887656,
            0.188263846,
            0.211646902,
            0.176051924,
            0.164989336,
            0.14444077,
            0.179925708,
            0.151049254,
            0.107294865,
            0.125720128,
            0.117940326,
            0.404257196)

  class <- apply(Wobs, 1, classify_func)
  b_hat <- ifelse(class != 17, temp[class], 0.0569 + 0.3822 * lwsd) #repeat by sptial unit
  #global r2: 0.726
}

#' Estimate AHG b lowerbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  temp <- c(0.009634045,
            0.009206385,
            0.010378626,
            0.00984031,
            0.000182357,
            0.008909195,
            0.025751932,
            0.005703884,
            0.017018127,
            0.013954247,
            0.019210795,
            0.010974071,
            0.016252382,
            0.017584834,
            0.004872167,
            0.029379248)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_b <- ifelse(class != 17, temp[class], 0.000182357)
  lowerbound_b <- min(lowerbound_b, na.rm = TRUE)
}

#' Estimate AHG b upperbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  temp <- c(0.401485695,
            0.450546257,
            0.442897059,
            0.427490668,
            0.415703708,
            0.388788209,
            0.422853699,
            0.41859221,
            0.397420318,
            0.432497404,
            0.404125905,
            0.413263659,
            0.361701335,
            0.328343393,
            0.348032765,
            0.773757585)

  class <- apply(Wobs, 1, classify_func)
  upperbound_b <- ifelse(class != 17, temp[class], 0.773757585)
  upperbound_b <- max(upperbound_b, na.rm = TRUE)
}

#' Estimate AHG b SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_bSD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  temp <- c(0.094625319,
            0.103986506,
            0.110231244,
            0.103910665,
            0.101259677,
            0.098511991,
            0.097920891,
            0.0962566,
            0.100083519,
            0.10412018,
            0.099927194,
            0.098731293,
            0.099282596,
            0.075781266,
            0.091252051,
            0.112180583)

  class <- apply(Wobs, 1, classify_func)
  b_sd <- ifelse(class != 17, temp[class], 0.068077044)
}

#' Estimate base cross-sectional area using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logA0 <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #Expert classification
  temp <- c(3.723280881,
            4.497073804,
            4.753590191,
            4.990432587,
            4.911328517,
            5.350615603,
            5.422742508,
            5.523458921,
            5.774532256,
            6.446836611,
            6.527953643,
            6.873163834,
            7.102499356,
            8.007965013,
            8.937204637,
            4.432006567)

  class <- apply(Wobs, 1, classify_func)
  logA0_hat <- ifelse(class != 17, temp[class], -0.2918 + 1.6930 * lwbar - 0.1887 * lwsd) #repeat for each sptial unit
  #global r2: 0.907
}

#' Estimate base cross-sectional area lowerbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundA0 <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #expert classification
  temp <- c(0.732367894,
            1.508511994,
            0.91027266,
            2.517696473,
            1.199964783,
            2.681021529,
            2.148850993,
            1.545432582,
            2.415913778,
            3.106826321,
            3.874321138,
            2.694627181,
            3.696351469,
            3.593194204,
            4.043928076,
            0.262364264)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_A0 <- ifelse(class != 17, exp(temp[class]), exp(4.540631665))
  lowerbound_A0 <- min(lowerbound_A0, na.rm = TRUE)
}

#' Estimate base cross-sectional area upperbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundA0 <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #expert classification
  temp <- c(7.640123173,
            7.355641103,
            8.997147152,
            9.164296433,
            8.554488976,
            9.417354541,
            7.677863501,
            8.144679183,
            7.863266724,
            8.793308627,
            8.776475789,
            9.014325488,
            8.78186249,
            9.61580548,
            11.6483301,
            11.55214618)

  class <- apply(Wobs, 1, classify_func)
  upperbound_A0 <- ifelse(class != 17, exp(temp[class]), 1000000)
  upperbound <- max(upperbound_A0, na.rm = TRUE)
}

#' Estimate base cross-sectional area SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_A0SD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #expert classification
  temp <- c(1.446820524,
            1.350055707,
            1.363875688,
            1.260787706,
            1.306619211,
            1.315826017,
            1.190542029,
            1.337220271,
            1.21047899,
            1.359096608,
            1.245863462,
            1.287297345,
            1.08535437,
            1.154319081,
            1.5575699,
            2.272342301)

  class <- apply(Wobs, 1, classify_func)
  logA0_sd <- ifelse(class != 17, temp[class], 0.58987527)
}

#'Estimate bankful width using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logWb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  #expert classification
  temp <- c(2.476118144,
            2.864001065,
            3.103015939,
            3.249308032,
            3.284178964,
            3.371669039,
            3.574693905,
            3.664586762,
            3.683922384,
            4.002696788,
            4.032912323,
            4.357733942,
            4.436574004,
            4.921742348,
            5.287893051,
            3.032064203)

  class <- apply(Wobs, 1, classify_func)
  logWb_hat <- ifelse(class != 17, temp[class], 0.0037 + 1.0028 * lwbar) #repeat for each sptial unit
  #global r2: 0.984
}

#'Estimate bankful width lower bound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundlogWb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  #expert classification
  temp <- c(0.316633894,
            0.814036498,
            0.461215123,
            1.377505855,
            0.744077909,
            2.070653036,
            1.73032723,
            0.42199441,
            1.756995477,
            2.396075436,
            2.293796587,
            2.268873179,
            2.99460664,
            1.852305657,
            2.585317436,
            -0.122732765)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_logWb <- ifelse(class != 17, temp[class], 3.089222617)
  lowerbound_logWb <- min(lowerbound_logWb, na.rm = TRUE)
}

#'Estimate bankful width upper bound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundlogWb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  #expert classification
  temp <- c(4.016563185,
            4.410978457,
            4.833579847,
            5.094823245,
            5.029162599,
            5.216130696,
            4.5361416,
            5.091077922,
            5.194372515,
            5.100415058,
            5.372171735,
            5.860073719,
            5.521860838,
            5.749870579,
            7.006785802,
            6.917259966)

  class <- apply(Wobs, 1, classify_func)
  upperbound_logWb <- ifelse(class != 17, temp[class], 10)
  upperbound_logWb <- max(upperbound_logWb, na.rm = TRUE)
}

#'Estimate bankful width SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logWbSD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  #expert classification
  temp <- c(0.750069561,
            0.738069364,
            0.710662053,
            0.669257108,
            0.685680677,
            0.630380021,
            0.563796126,
            0.812622517,
            0.691428554,
            0.654155675,
            0.587954845,
            0.728133987,
            0.540567805,
            0.680678526,
            0.889425931,
            1.268865384)

  class <- apply(Wobs, 1, classify_func)
  logWb_sd <- ifelse(class != 17, temp[class], 0.137381044)
}

#'Estimate bankful depth using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logDb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #expert classification
  temp <- c(-1.035199401,
            -0.738778811,
            -0.759250143,
            -0.714358207,
            -0.834152133,
            -0.629032654,
            -0.471723514,
            -0.547540836,
            -0.480357291,
            0.130087152,
            -0.07180952,
            0.15467964,
            0.299392088,
            0.579044554,
            1.48611948,
            -0.974789558)

  class <- apply(Wobs, 1, classify_func)
  logDb_hat <- ifelse(class != 17, temp[class], -2.6189 - 0.2436 * lwsd + 0.6854 * lwbar) #repeat for each spatial unit
  #global r2: 0.640
}

#'Estimate bankful depth lower bound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundlogDb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #expert classification
  temp <- c(-2.326877786,
            -1.788516004,
            -2.05494407,
            -1.892425141,
            -2.057110918,
            -1.947205432,
            -1.915682003,
            -1.514919292,
            -1.730529988,
            -1.638497383,
            -1.809169822,
            -1.622549192,
            -1.649478962,
            -0.83244782,
            -0.924234286,
            -2.379581849)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_logDb <- ifelse(class != 17, temp[class], -0.951962577)
  lowerbound_logDb <- min(lowerbound_logDb, na.rm = TRUE)
}

#'Estimate bankful depth upper bound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundlogDb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #expert classification
  temp <- c(1.49936159,
            1.045713084,
            1.7933029,
            1.725376856,
            1.127811443,
            1.827833367,
            1.174682984,
            1.449780646,
            1.260032945,
            1.708874134,
            1.765116599,
            1.327316111,
            1.639749354,
            2.102212505,
            3.309358647,
            2.571816749)

  class <- apply(Wobs, 1, classify_func)
  upperbound_logDb <- ifelse(class != 17, temp[class], 5)
  upperbound_logDb <- max(upperbound_logDb, na.rm = TRUE)
}

#'Estimate bankful depth SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logDbSD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #expert classification
  temp <- c(0.754179034,
            0.719470145,
            0.751836561,
            0.699278157,
            0.699865238,
            0.771195464,
            0.69921735,
            0.665457667,
            0.662592139,
            0.810898004,
            0.764674068,
            0.719504763,
            0.690395904,
            0.664391944,
            0.8459272,
            1.132642792)

  class <- apply(Wobs, 1, classify_func)
  upperbound_logDb <- ifelse(class != 17, temp[class], 0.576212733)
}

#'Estimate channel shape using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logr <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #expert classification
  temp <- c(0.479582219,
            0.475424989,
            0.742030366,
            0.747458374,
            0.810066771,
            0.764035146,
            0.728537155,
            0.914703855,
            1.012946569,
            1.012082534,
            0.964521619,
            1.007259927,
            1.368493142,
            1.217111712,
            1.067996942,
            -0.247245593)

  class <- apply(Wobs, 1, classify_func)
  logr_hat <- ifelse(class != 17, temp[class], 1.4241 - 1.9097 * lwsd + 0.0420 * lwbar) #repeat for each spatial unit
  #Global r2: 0.421
}

#'Estimate channel shape lowerbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundlogr <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #expert classification
  temp <- c(0.010687262,
            0.003572302,
            0.030800262,
            0.055237814,
            0.075268048,
            0.005163827,
            0.003031095,
            0.140133893,
            0.014043795,
            0.007552884,
            0.033288292,
            0.017194676,
            0.081171573,
            0.128060852,
            0.020365892,
            -2.580471126)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_logr <- ifelse(class != 17, temp[class], -2.580471126)
  lowerbound_logr <- min(lowerbound_logr, na.rm = TRUE)
}

#'Estimate channel shape upperbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundlogr <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #expert classification
  temp <- c(2.659590969,
            3.660494543,
            3.75765632,
            3.809886281,
            8.037716276,
            3.885278632,
            3.141534762,
            4.058323672,
            2.894177697,
            3.145908752,
            3.291615219,
            3.759676814,
            3.688940701,
            3.293279235,
            4.028976176,
            -0.002840899)

  class <- apply(Wobs, 1, classify_func)
  upperbound_logr <- ifelse(class != 17, temp[class], 8.037716276)
  upperbound_logr <- max(upperbound_logr, na.rm = TRUE)
}

#'Estimate channel shape SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logrSD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #expert classification
  temp <- c(0.574788166,
            0.686443448,
            0.735995097,
            0.660600833,
            1.152463807,
            0.751672132,
            0.646178286,
            0.721545901,
            0.678734433,
            0.73732522,
            0.714799269,
            0.783012838,
            0.824128369,
            0.666709963,
            0.741805733,
            0.418335751)

  class <- apply(Wobs, 1, classify_func)
  logr_sd <- ifelse(class != 17, temp[class], 0.67332688)
}

#'Estimate manning's n using bam data
#'
#' @param Sobs Observed S, as a space-down, time-across matrix
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_logn <- function(Sobs, Wobs) {
  Sobs[Sobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lsbar <- apply(log(Sobs), 1, mean, na.rm = TRUE)
  #expert classification
  temp <- c(-2.956882373,
            -3.218074003,
            -3.073438472,
            -3.316850463,
            -3.318215717,
            -3.15335243,
            -3.456243594,
            -3.545592423,
            -3.240085716,
            -3.401877538,
            -3.269861085,
            -3.372948626,
            -3.404202242,
            -3.274729636,
            -3.405497386,
            -3.318138928)

  class <- apply(Wobs, 1, classify_func)
  logn_hat <- ifelse(class != 17, temp[class], -0.1636 + 0.4077 * lsbar) #repeat for each sptial unit
  #Global r2: 0.631
}

#'Estimate manning's n lowerbound using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_lowerboundlogn <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #expert classification
  temp <- c(-5.678333858,
            -5.467593416,
            -5.67856698,
            -5.768715685,
            -5.682940688,
            -5.841321589,
            -5.56434293,
            -5.957755533,
            -5.646550723,
            -6.064509006,
            -6.029056118,
            -6.405084634,
            -6.678294272,
            -6.211213087,
            -6.06854241,
            -6.246537799)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_logn <- ifelse(class != 17, temp[class], log(0.01))
  lowerbound_logn <- min(lowerbound_logn, na.rm = TRUE)
}

#'Estimate manning's n upperbound using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_upperboundlogn <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #expert classification
  temp <- c(-0.377797559,
            -1.024617813,
            -0.469297566,
            -1.330691863,
            -1.17053266,
            -1.181653683,
            -0.819808414,
            -1.085062606,
            -0.540120016,
            -1.349783005,
            -0.325687448,
            -0.025491812,
            -0.554602291,
            0.501441352,
            1.930473909,
            0.094113848)

  class <- apply(Wobs, 1, classify_func)
  upperbound_logn <- ifelse(class != 17, temp[class], log(0.05))
  upperbound_logn <- max(upperbound_logn, na.rm = TRUE)
}

#'Estimate manning's n SD using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_lognSD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #expert classification
  temp <- c(1.069101032,
            1.004209047,
            1.051391594,
            1.050420395,
            1.103812115,
            1.050696519,
            1.205687421,
            1.229666517,
            1.101795647,
            1.194241604,
            1.268715794,
            1.260820049,
            1.199816945,
            1.345333569,
            1.724262114,
            1.122652969)

  class <- apply(Wobs, 1, classify_func)
  logn_sd <- ifelse(class != 17, temp[class], 0.761673112)
}


# Prior calculation using unsupervised framework------------------------------------------------------------------------------
#class 101 are noisy rivers

#' Estimate AHG b exponent using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_b_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwsd <- apply(log(Wobs), 1, function(x) sd(x, na.rm = TRUE))

  #unsupervised classification
  temp <- c(0.193965612,
            0.186677399,
            0.053815029,
            0.193643212,
            0.171665843,
            0.297094389,
            0.409781961)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  b_hat <- ifelse(class != 101, temp[class], 0.0569 + 0.3822 * lwsd) #repeat by sptial unit
  #global r2: 0.726
}

#' Estimate AHG b lowerbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundb_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  temp <- c(0.004872167,
            0.038484492,
            0.022462362,
            0.029379248,
            0.077325498,
            0.024351422,
            0.103600981)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  lowerbound_b <- ifelse(class != 101, temp[class], 0.000182357)
  lowerbound_b <- min(lowerbound_b, na.rm = TRUE)
}

#' Estimate AHG b upperbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundb_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  temp <- c(0.773757585,
            0.437298948,
            0.478046051,
            0.514297361,
            0.488190912,
            0.455487205,
            0.692551528)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  upperbound_b <- ifelse(class != 101, temp[class], 0.773757585)
  upperbound_b <- max(upperbound_b, na.rm = TRUE)
}

#' Estimate AHG b SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_bSD_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  temp <- c(0.123620812,
            0.133008102,
            0.159164985,
            0.160152755,
            0.151764685,
            0.139129593,
            0.123429852)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  b_sd <- ifelse(class != 101, temp[class], 0.068077044)
}

#' Estimate base cross-sectional area using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logA0_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #unsupervised classification
  temp <- c(5.538317552,
            7.31986493,
            8.157657015,
            7.835974582,
            5.329969058,
            6.223916904,
            3.019662854)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logA0_hat <- ifelse(class != 101, temp[class], -0.2918 + 1.6930 * lwbar - 0.1887 * lwsd) #repeat for each sptial unit
  #global r2: 0.907
}

#' Estimate base cross-sectional area lowerbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundA0_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #unsupervised classification
  temp <- c(0.385262401,
            2.05284086,
            4.952299717,
            5.059425458,
            3.446807893,
            1.596352673,
            0.262364264)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  lowerbound_A0 <- ifelse(class != 101, exp(temp[class]), exp(-0.328504067))
  lowerbound_A0 <- min(lowerbound_A0, na.rm = TRUE)
}

#' Estimate base cross-sectional area upperbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundA0_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #unsupervised classification
  temp <- c(11.6483301,
            9.952277717,
            8.942460927,
            9.472704636,
            7.316548177,
            9.164296433,
            6.200509174)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  upperbound_A0 <- ifelse(class != 101, exp(temp[class]), exp(11.6483301))
  upperbound_A0 <- max(upperbound_A0, na.rm = TRUE)
}

#' Estimate base cross-sectional area SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_A0SD_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #unsupervised classification
  temp <- c(1.83150945,
            2.054790651,
            1.248676431,
            0.954728776,
            1.426273172,
            2.302489877,
            1.252208404)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logA0_sd <- ifelse(class != 101, temp[class], 0.58987527)
}

#'Estimate bankful width using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logWb_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  #unsupervised classification
  temp <- c(3.574730432,
            4.392286324,
            4.361632583,
            4.461530736,
            3.769275439,
            3.890165778,
            2.060762964)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logWb_hat <- ifelse(class != 101, temp[class], 0.0037 + 1.0028 * lwbar) #repeat for each sptial unit
  #global r2: 0.984
}

#'Estimate bankful width lower bound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundlogWb_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  #unsupervised classification
  temp <- c(-0.122732765,
            1.173410499,
            2.971439581,
            3.417726684,
            1.6804554,
            1.042570898,
            0.381172416)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  lowerbound_logWb <- ifelse(class != 101, temp[class], -0.122732765)
  lowerbound_logWb <- min(lowerbound_logWb, na.rm = TRUE)
}

#'Estimate bankful width upper bound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundlogWb_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  #unsupervised classification
  temp <- c(6.917259966,
            5.427282098,
            4.91287545,
            5.860073719,
            4.670489652,
            5.56832542,
            3.948354935)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  upperbound_logWb <- ifelse(class != 101, temp[class], 7.006785802)
  upperbound_logWb <- max(upperbound_logWb, na.rm = TRUE)
}

#'Estimate bankful width SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logWbSD_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  #unsupervised classification
  temp <- c(1.026265531,
            1.051896455,
            0.597547965,
            0.533663213,
            1.093140262,
            1.318621719,
            0.830429076)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logWb_sd <- ifelse(class != 101, temp[class], 0.137381044)
}

#'Estimate bankful depth using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logDb_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #unsupervised classification
  temp <- c(-0.477454537,
            0.570414415,
            1.402823663,
            0.829185144,
            -0.618893436,
            -0.325969172,
            -1.458486344)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logDb_hat <- ifelse(class != 101, temp[class], -2.6189 - 0.2436 * lwsd + 0.6854 * lwbar) #repeat for each spatial unit
  #global r2: 0.640
}

#'Estimate bankful depth lower bound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundlogDb_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #unsupervised classification
  temp <- c(-2.241241775,
            -1.502582421,
            -0.688452336,
            -0.718666115,
            -0.856326031,
            -1.894752997,
            -2.379581849)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  lowerbound_logDb <- ifelse(class != 101, temp[class], -3.020024966)
  lowerbound_logDb <- min(lowerbound_logDb, na.rm = TRUE)
}

#'Estimate bankful depth upper bound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundlogDb_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #unsupervised classification
  temp <- c(3.309358647,
            2.145955751,
            2.009770826,
            1.729110741,
            0.249099039,
            1.7933029,
            -0.221917668)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  upperbound_logDb <- ifelse(class != 101, temp[class], 3.309358647)
  upperbound_logDb <- max(upperbound_logDb, na.rm = TRUE)
}

#'Estimate bankful depth SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logDbSD_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #unsupervised classification
  temp <- c(0.909913045,
            1.11158809,
            0.709944526,
            0.609960942,
            0.465557504,
            1.157661125,
            0.49869352)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logDb_sd <- ifelse(class != 101, temp[class], 0.576212733)
}

#'Estimate channel shape using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logr_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #supervised classification
  temp <- c(0.722592009,
            0.977709482,
            1.230106834,
            0.303821877,
            0.917785312,
            0.141431684,
            -0.248390176) #DBSCAN successfully identified rivers with r < 1!!

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logr_hat <- ifelse(class != 101, temp[class], 1.4241 - 1.9097 * lwsd + 0.0420 * lwbar) #repeat for each spatial unit
  #Global r2: 0.421
}

#'Estimate channel shape lowerbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundlogr_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #unsupervised classification
  temp <- c(-1.694462984,
            -2.580471126,
            -0.188149031,
            -2.297051297,
            -0.374243132,
            -0.636077687,
            -0.827781183)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  lowerbound_logr <- ifelse(class != 101, temp[class], -2.580471126)
  lowerbound_logr <- min(lowerbound_logr, na.rm = TRUE)
}

#'Estimate channel shape upperbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundlogr_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #unsupervised classification
  temp <- c(3.885278632,
            2.701880636,
            2.468097757,
            2.311888633,
            1.777431218,
            2.674955812,
            0.9913607)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  upperbound_logr <- ifelse(class != 101, temp[class], 8.037716276)
  upperbound_logr <- max(upperbound_logr, na.rm = TRUE)
}

#'Estimate channel shape SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logrSD_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #unsupervised classification
  temp <- c(0.796841513,
            1.174459584,
            0.85622092,
            1.205119078,
            0.744722627,
            0.955376128,
            0.429377882)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logr_sd <- ifelse(class != 101, temp[class], 0.67332688)
}

#'Estimate manning's n using bam data
#'
#' @param Sobs Observed S, as a space-down, time-across matrix
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_logn_unsupervised <- function(Sobs, Wobs) {
  Sobs[Sobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lsbar <- apply(log(Sobs), 1, mean, na.rm = TRUE)
  #unsupervised classification
  temp <- c(-3.293181606,
            -3.705388959,
            -0.735356858,
            -3.156147655,
            -3.95179067,
            -3.406376251,
            -3.121679163)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logn_hat <- ifelse(class != 101, temp[class], -0.1636 + 0.4077 * lsbar) #repeat for each sptial unit
  #Global r2: 0.631
}

#'Estimate manning's n lowerbound using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_lowerboundlogn_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #unsupervised classification
  temp <- c(-6.678294272,
            -5.933075927,
            -1.881832288,
            -4.865235717,
            -5.661908145,
            -5.42815632,
            -5.577638217)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  lowerbound_logn <- ifelse(class != 101, temp[class], log(0.01))
  lowerbound_logn <- min(lowerbound_logn, na.rm = TRUE)
}

#'Estimate manning's n upperbound using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_upperboundlogn_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #unsupervised classification
  temp <- c(-0.025491812,
            -1.681204649,
            0.501441352,
            -1.270558784,
            -2.631875724,
            -1.540988288,
            -0.377797559)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  upperbound_logn <- ifelse(class != 101, temp[class], log(0.05))
  upperbound_logn <- max(upperbound_logn, na.rm = TRUE)
}

#'Estimate manning's n SD using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_lognSD_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #unsupervised classification
  temp <- c(1.137136185,
            1.107362891,
            0.679465988,
            0.885903042,
            1.311510864,
            1.096048948,
            0.925151557)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logn_sd <- ifelse(class != 101, temp[class], 0.761673112)
}
                
                # utility functions

#' Convert coefficient of variation to sigma parameter of lognormal diistribution
#'
#' @param cv Coefficient of variation
#' @export

cv2sigma <- function (cv) {
  sqrt(log(cv^2 + 1))
}



# functions for getting Q bounds.

#' Minimum across xs of max across time of width
#'
#' @param x a numeric matrix
minmax <- function(x)
  min(apply(x, 1, max, na.rm=TRUE), na.rm=TRUE)

#' Maximum across xs of min across time of width
#'
#' @param x a numeric matrix
maxmin <- function(x) {
  max(apply(x, 1, min, na.rm=TRUE), na.rm=TRUE)
}

#Functions for classifying rivers----------------------------------------------------------------------

#'Classify river for expert framework
#'
#'@param Wobs observed widths matrix
classify_func <- function(Wobs) {
  lwbar <- mean(log(Wobs), na.rm=TRUE)
  lwsd <- sd(log(Wobs), na.rm= TRUE)

  maxWidth = 6.5
  classes <- c(2.476118144,
               2.864001065,
               3.103015939,
               3.249308032,
               3.284178964,
               3.371669039,
               3.56827873,
               3.664586762,
               3.683922384,
               4.002696788,
               4.031559142,
               4.357733942,
               4.436574004,
               4.921166637,
               5.287893051) #median width of each river type
  index <- ifelse(lwbar > maxWidth, 17, which.min(abs(classes-lwbar))) #17 for big rivers
  index <- ifelse(lwsd >= 0.45, 16, index)  #16 for width-variable rivers, which overrides 'big' rivers
  return(index)
}

#'Classify river for unsupervised framework
#'
#' @param Wobs observed widths matrix
classify_func_unsupervised <- function(Wobs) {
  width <- median((Wobs), na.rm = TRUE) #note these are not log widths!

  #One-vs-Rest Logistic regression model
  #test set accuracy rate of 87%
  p_noise <- 1 / (1.0 + exp(-(-3.11056777 + 0.00189261 * width)))
  p_1 <- 1 / (1.0 + exp(-(1.95763357 + -0.00114445 * width)))
  p_2 <- 1 / (1.0 + exp(-(-4.05274822 + 0.00156913 * width)))
  p_3<- 1 / (1.0 + exp(-(-4.43158034 + -0.00116211 * width)))
  p_4 <- 1 / (1.0 + exp(-(-3.99272449 + 0.00168791* width)))
  p_5 <- 1 / (1.0 + exp(-(-4.28204232 + -0.01493184 * width)))
  p_6 <- 1 / (1.0 + exp(-(-3.87408245 + -0.00175089 * width)))
  p_7 <- 1 / (1.0 + exp(-(-1.16209678 + -0.12664823 * width)))

  probs <- data.frame(p_noise, p_1, p_2, p_3, p_4, p_5, p_6, p_7)

  index <- which.max(probs)
  index <- ifelse(index == 1, 101, index-1) #101 for 'noisey' rivers
  return(index)
}
