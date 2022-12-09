#' @title Fit a gpglm model (spBayes spGLM)
#'
#' @description A more user friendly (and limited) wrapper around the spBayes' spLM
#' model.
#'
#' @param formula A formula object. A symbolic description of the model to be fitted.
#' @param family Currently only supports binomial and poisson data using the logit
#' and log link functions, respectively.
#' @param weights An optional vector of weights to be used in the fitting process.
#' Weights correspond to number of trials and offset for each location for the binomial and poisson family, respectively.
#' @param data A data.frame object containing data of all variables used in the model.
#' @param coords Either a numeric matrix of coordinates or quoted column names from `data` to use for distance calculations.
#' @param chains Numeric. Number of MCMC chains to fit.
#' @param n_batch Numeric. Number of adaptive batches to use when fitting each chain.
#' @param batch_len Numeric. Number of iterations per batch.
#' @param cov_model Character. Name of the covariance model to use, supported values are:
#' "exponential", "matern", "spherical", and "gaussian"
#' @param starting Named list of parameter starting values, allowed names: `beta`, `sigma.sq`, `tau.sq`, `phi`, and `nu`.
#' @param prior Named list of priors, each value is a vector of prior hyperparameter values. See [spBayes::spGLM()] for details.
#' @param tuning Names list of variance values for the MH sampler.
#' @param burnin_frac Numeric. Proportion of iterations to discard as burnin.
#' @param accept_rate Numeric. Desired acceptance rate used by the adaptive MCMC algorithm
#' @param n_report Numeric. The interval to report Metropolis sampler acceptance and MCMC progress.
#' @param verbose Logical. Should verbose output (sampling progress) be printed.
#'
#' @examples
#'
#' eta = function(x) 2.5*sin(2.1*pi*(x-0.25))
#' d = data.frame(x=runif(n)) |>
#'   dplyr::mutate(
#'     eta = eta(x),
#'     p = 1/(1+exp(-eta)),
#'     y = rbinom(dplyr::n(), size=1, prob = p)
#'   )
#'
#' m = gpglm(
#'   y~1, family="binomial",
#'   data = d, coords = c("x"),
#'   cov_model = "exponential",
#'   starting=list(
#'     "beta"=0, "phi"=3/0.1, "sigma.sq"=1, "w"=0
#'   ),
#'   tuning=list(
#'     "beta"=0.5, "phi"=0.5, "sigma.sq"=0.5, "w"=0.5
#'   ),
#'   priors=list(
#'     "beta.Normal"=list(0,1),
#'     "phi.unif"=c(sqrt(3)/0.5, sqrt(3)/0.01),
#'     "sigma.sq.ig"=c(2, 1)
#'   ),
#'   n_batch = 100,
#'   batch_len = 100,
#'   verbose = TRUE
#' )
#'
#' newdata = data.frame(
#'   x=seq(0,1,length.out=101)
#' )
#'
#' p = predict(m, newdata=newdata, coords="x", thin=50)
#'
#' # Predicted y
#' p |>
#'   tidybayes::gather_draws(y[i]) |>
#'   ggplot2::ggplot(ggplot2::aes(x=i,y=.value)) +
#'     tidybayes::stat_lineribbon()
#'
#' # Predicted w
#' p |>
#'   tidybayes::gather_draws(w[i]) |>
#'   ggplot2::ggplot(ggplot2::aes(x=i,y=.value)) +
#'     tidybayes::stat_lineribbon()
#'
#' @export
#'
gpglm = function(
    formula,
    family = "binomial",
    weights = rep(1, nrow(data)),
    data,
    coords,
    chains = 4,
    n_batch = 200,
    batch_len = 100,
    cov_model = c("gaussian", "exponential", "matern", "spherical"),
    starting = list(
      beta=0, "phi"=sqrt(3)/1, "sigma.sq"=1
    ),
    priors = list(
      "beta.Normal"=list(0,1),
      "phi.unif"=c(sqrt(3)/52, sqrt(3)/1),
      "sigma.sq.ig"=c(2, 1)
    ),
    tuning = list(
      "beta"=1, "phi"=1, "sigma.sq"=1
    ),
    burnin_frac = 0.5,
    accept_rate = 0.43,
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

      model = spBayes::spGLM(
        formula = formula,
        family = family,
        weights = weights,
        data = data,
        coords = coords,
        starting = starting,
        priors = priors,
        tuning = tuning,
        amcmc = list(
          "n.batch"=n_batch, "batch.length"=batch_len, "accept.rate"=accept_rate
        ),
        cov.model = cov_model,
        verbose = verbose,
        n.report = n_report
      )

      list(
        model = model,
        mcmc = posterior::as_draws_matrix(
          model$p.beta.theta.samples
        )
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
    class = "gpglm_fit"
  )
}

#' @exportS3Method
predict.gpglm_fit = function(
    object,
    newdata,
    coords,
    verbose = FALSE,
    thin=1,
    burnin_frac = 0.5,
    ...
) {
  args = object$args

  end = nrow(object$models[[1]]$mcmc)
  start = floor(end * burnin_frac) + 1


  coords = build_coord_mat(newdata, coords, ncol=2)

  lapply(
    object$models,
    function(m) {
      n_y = nrow(newdata)

      covars = stats::model.matrix(
        stats::update(object$args$formula, NULL ~ .),
        newdata
      )

      p = spBayes::spPredict(
        m$model,
        pred.coords = coords,
        pred.covars = covars,
        verbose = verbose,
        thin = thin,
        start = start,
        end = end,
        ...
      )

      rownames(p$p.w.predictive.samples) = stringr::str_glue("w[{seq_len(n_y)}]")
      rownames(p$p.y.predictive.samples) = stringr::str_glue("y[{seq_len(n_y)}]")

      posterior::bind_draws(
        posterior::as_draws(t(p$p.w.predictive.samples)),
        posterior::as_draws(t(p$p.y.predictive.samples)),
        along = "variable"
      )
    }
  ) |>
    purrr::reduce(posterior::bind_draws, along="chain")
}



#' @exportS3Method
plot.gpglm_fit = function(x, combo = c("dens", "trace"), ..., vars=NULL) {
  mcmc = x$mcmc
  if (!is.null(vars))
    mcmc = mcmc[,vars]

  print(bayesplot::mcmc_combo(mcmc, combo = combo, ...))
}

#' @exportS3Method
print.gpglm_fit = function(x, ...) {
  nchains = posterior::nchains(x$mcmc)
  ndraws = posterior::ndraws(x$mcmc)
  nvar = posterior::nvariables(x$mcmc)

  cat( stringr::str_glue(
    "# A gpglm model (spBayes spGLM) with {nchains} chains, {nvar} variables, and {ndraws} iterations."
  ), "\n" )
  summary(x)
}

#' @exportS3Method
summary.gpglm_fit = function(object, ...) {
  print(summary(object$mcmc))
}

#' @exportS3Method
print.spGLM = function(x, ...) {
  print("spGLM model object")
  invisible(x)
}
