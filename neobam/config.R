GLOBAL_PARAMS=list(
  perc_lower=0.3, #lower percentile from which to sample for the window
  perc_upper=0.7, #upper percentie from which to smaple for the window
  nt_window=20, #how many times are in the window? direct tradeoff with runtime
  Werr_sd=45, #expected uncertainty in width observations
  Serr_sd=0.001, #expected uncerainty in slope observations
  iter=2000, #how many iterations in the markov chain? (was 5000)
  sigma_man=0.26, #structural error in flow law
  nx_sample =20
)
