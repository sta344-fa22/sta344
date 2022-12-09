get_mcmc = function(rec) {
  posterior::bind_draws(
    posterior::as_draws_matrix(rec$p.theta.recover.samples),
    posterior::as_draws_matrix(rec$p.beta.recover.samples),
    along = "variable"
  )
}

pad_coords_mat = function(coords, ncol=2) {
  stopifnot(ncol(coords) <= ncol)

  cbind(
    coords,
    matrix(0, nrow = nrow(coords), ncol = ncol - ncol(coords))
  )
}

build_coord_mat = function(data, coords, ncol=2) {
  if (!is.matrix(coords)) {
    if (is.numeric(coords) & length(coords) == nrow(data)) {
      coords = matrix(coords, ncol=1)
    } else {
      # Assume select syntax for coords
      coords = data |>
        dplyr::select(dplyr::all_of(coords)) |>
        dplyr::mutate(dplyr::across(everything(), as.numeric)) |>
        as.matrix()
    }
  }

  storage.mode(coords) = "double"

  stopifnot(all(!is.na(coords)))

  pad_coords_mat(coords, ncol)
}


#' @title Fit a gplm model (spBayes spLM)
#'
#' @description A more user friendly (and limited) wrapper around the spBayes' spLM
#' model.
#'
#' @param formula A formula object. A symbolic description of the model to be fitted.
#' @param data A data.frame object containing data of all variables used in the model.
#' @param coords Either a numeric matrix of coordinates or quoted column names from `data` to use for distance calculations.
#' @param chains Numeric. Number of MCMC chains to fit.
#' @param n_batch Numeric. Number of adaptive batches to use when fitting each chain.
#' @param batch_len Numeric. Number of iterations per batch.
#' @param cov_model Character. Name of the covariance model to use, supported values are:
#' "exponential", "matern", "spherical", and "gaussian"
#' @param starting Named list of parameter starting values, allowed names: `beta`, `sigma.sq`, `tau.sq`, `phi`, and `nu`.
#' @param prior Named list of priors, each value is a vector of prior hyperparameter values. See [spBayes::spLM()] for details.
#' @param tuning Names list of variance values for the MH sampler.
#' @param burnin_frac Numeric. Proportion of iterations to discard as burnin.
#' @param accept_rate Numeric. Desired acceptance rate used by the adaptive MCMC algorithm
#' @param thin Numeric. Amount of thinning to apply to posterior samples.
#' @param n_report Numeric. The interval to report Metropolis sampler acceptance and MCMC progress.
#' @param verbose Logical. Should verbose output (sampling progress) be printed.
#'
#' @examples
#' m = gplm(
#'   avg_temp~1, data = avg_temp, coords = "week",
#'   starting=list(
#'     "phi"=sqrt(3)/4, "sigma.sq"=1, "tau.sq"=1
#'   ),
#'   tuning=list(
#'     "phi"=1, "sigma.sq"=1, "tau.sq"=1
#'   ),
#'   priors=list(
#'     "phi.unif"=c(sqrt(3)/52, sqrt(3)/1),
#'     "sigma.sq.ig"=c(2, 1),
#'     "tau.sq.ig"=c(2, 1)
#'   ),
#'   thin=10,
#'   n_batch = 100,
#'   batch_len = 50
#' )
#'
#' newdata = data.frame(
#'   week = seq(0,3.5*52) |> jitter()
#' )
#'
#' (pred = predict(m, newdata=newdata, coords = "week"))
#'
#' @export
#'
gplm = function(
    formula,
    data,
    coords,
    chains = 4,
    n_batch = 200,
    batch_len = 100,
    cov_model = c("gaussian", "exponential", "matern", "spherical"),
    starting = list(
      "phi"=sqrt(3)/1, "sigma.sq"=1, "tau.sq"=1
    ),
    priors = list(
      "phi.unif"=c(sqrt(3)/52, sqrt(3)/1),
      "sigma.sq.ig"=c(2, 1),
      "tau.sq.ig"=c(2, 1)
    ),
    tuning = list(
      "phi"=1, "sigma.sq"=1, "tau.sq"=1
    ),
    burnin_frac = 0.5,
    accept_rate = 0.43,
    thin = 1,
    n_report = 50,
    verbose = FALSE
) {
  args = as.list(environment())
  n_samples = n_batch*batch_len
  burnin = floor(n_samples*burnin_frac + 1)

  coords = build_coord_mat(data, coords)

  cov_model = match.arg(cov_model)

  l = lapply(
    seq_len(chains),
    function(i) {
      message("Fitting chain ", i)

      model = spBayes::spLM(
        formula,
        data = data,
        coords = coords,
        starting = starting,
        priors = priors,
        tuning = tuning,
        amcmc = list(
          "n.batch"=n_batch, "batch.length"=batch_len, "accept.rate"=accept_rate
        ),
        cov.model = cov_model,
        n.report = n_report,
        verbose = verbose
      )

      rec = spBayes::spRecover(
        model,
        start = burnin, thin = thin,
        get.w = FALSE,
        verbose = verbose
      )

      list(
        model = model,
        rec = rec,
        mcmc = get_mcmc(rec)
      )
    }
  )

  structure(
    list(
      models = l,
      args = args,
      mcmc = purrr::map(l, "mcmc") |>
        purrr::reduce(posterior::bind_draws, along="chain")
    ),
    class = "gplm_fit"
  )
}


#' @exportS3Method
predict.gplm_fit = function(
    object,
    newdata,
    coords,
    verbose = FALSE,
    n_report = 100,
    ...
) {
  args = object$args

  coords = build_coord_mat(newdata, coords, ncol=2)

  lapply(
    object$models,
    function(m) {
      n_y = nrow(newdata)

      covars = stats::model.matrix(
        stats::update(object$args$formula, NULL ~ .),
        newdata
      )

      sp_predict_gplm(
        m$rec,
        coords = coords,
        covars = covars,
        verbose = verbose,
        n.report = n_report
      )
    }
  ) |>
    purrr::reduce(posterior::bind_draws, along="chain")
}

# Copied from
# https://github.com/cran/spBayes/blob/e890f7c1e76af4020e9a3e16020f651512e28380/R/spPredict.R#L498-L624
# * Removed all pp stuff
# * Simplified checks and related stuff

sp_predict_gplm = function(
    gplm,
    coords,
    covars,
    verbose = FALSE,
    n.report = 100
) {
  stopifnot(is.matrix(coords))
  stopifnot(ncol(coords) == 2)

  X = gplm$X
  Y = gplm$Y
  p = ncol(X)
  n = nrow(X)
  obs.coords = gplm$coords
  cov.model = gplm$cov.model
  p.theta.samples = gplm$p.theta.recover.samples
  n.samples = nrow(p.theta.samples)
  nugget = gplm$nugget
  beta.prior = gplm$beta.prior
  beta.Norm = gplm$beta.Norm

  beta = gplm$p.beta.recover.samples

  stopifnot(!is.null(gplm$p.beta.recover.samples))

  sigma.sq.indx = 0; tau.sq.indx = 0; phi.indx = 0; nu.indx = 0

  if(!nugget && cov.model != "matern"){
    sigma.sq.indx = 0; phi.indx = 1
  }else if(nugget && cov.model != "matern"){
    sigma.sq.indx = 0; tau.sq.indx = 1; phi.indx = 2
  }else if(!nugget && cov.model == "matern"){
    sigma.sq.indx = 0; phi.indx = 1; nu.indx = 2
  }else{
    sigma.sq.indx = 0; tau.sq.indx = 1; phi.indx = 2; nu.indx = 3
  }

  obs.pred.D = spBayes::iDist(obs.coords, coords)
  obs.D = spBayes::iDist(obs.coords)

  n.pred = nrow(coords)

  storage.mode(X) = "double"
  storage.mode(Y) = "double"
  storage.mode(n) = "integer"
  storage.mode(p) = "integer"
  storage.mode(covars) = "double"
  storage.mode(n.pred) = "integer"
  storage.mode(p.theta.samples) = "double"
  storage.mode(n.samples) = "integer"
  storage.mode(beta) = "double"
  storage.mode(sigma.sq.indx) = "integer"
  storage.mode(tau.sq.indx) = "integer"
  storage.mode(phi.indx) = "integer"
  storage.mode(nu.indx) = "integer"
  storage.mode(verbose) = "integer"
  storage.mode(n.report) = "integer"
  storage.mode(obs.pred.D) = "double"
  storage.mode(obs.D) = "double"

  p = .Call(
    "spLMPredict", X, Y, n, p, covars, n.pred,
    p.theta.samples, n.samples,
    beta.prior, beta.Norm, beta, sigma.sq.indx, tau.sq.indx, phi.indx, nu.indx,
    obs.D, obs.pred.D, cov.model, nugget, verbose, n.report,
    PACKAGE = "spBayes"
  )$p.y.predictive.samples

  rownames(p) = stringr::str_glue("y[{seq_len(nrow(p))}]")
  posterior::as_draws(t(p))
}


#' @exportS3Method
plot.gplm_fit = function(x, combo = c("dens", "trace"), ..., vars=NULL) {
  mcmc = x$mcmc
  if (!is.null(vars))
    mcmc = mcmc[,vars]

  print(bayesplot::mcmc_combo(mcmc, combo = combo, ...))
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
