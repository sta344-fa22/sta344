#' @exportS3Method
plot.gplm_fit = function(x, combo = c("dens", "trace"), ...) {
  print(bayesplot::mcmc_combo(x$mcmc, combo = combo, ...))
}

#' @exportS3Method
print.gplm_fit = function(x, ...) {
  nchains = posterior::nchains(x$mcmc)
  ndraws = posterior::ndraws(x$mcmc)
  nvar = posterior::nvariables(x$mcmc)

  cat( stringr::str_glue(
    "# A gplm model (spBayes spLM) with {nchains} chains, {nvar} variables, and {ndraws} iterations."
  ), "\n" )
  summary(x)
}

#' @exportS3Method
summary.gplm_fit = function(object, ...) {
  print(summary(object$mcmc))
}

#' @exportS3Method
print.spLM = function(x, ...) {
  print("spLM model object")
  invisible(x)
}
