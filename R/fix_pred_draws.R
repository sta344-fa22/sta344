#' @name draws_fix
#'
#' @title Stopgap fixes for tidybayes draws functions
#'
#' @description Wrappers around tidybayes draw functions that maintain the
#' `.iteration`, `.draw`, and `.chain` columns.
#'
#' @param object A supported Bayesian model fit that can provide fits and predictions.
#' @param newdata Data frame to generate predictions from.
#' @param ... Any additional arguments to be passed to the base function.
#'
#' @seealso [tidybayes::predicted_draws()]
#'
NULL


#' @rdname draws_fix
#' @export
predicted_draws_fix = function(object, newdata, ...) {
  fix_draws(object, newdata, ..., func = tidybayes::predicted_draws)
}

#' @rdname draws_fix
#' @export
epred_draws_fix = function(object, newdata, ...) {
  fix_draws(object, newdata, ..., func = tidybayes::epred_draws)
}

#' @rdname draws_fix
#' @export
residual_draws_fix = function(object, newdata, ...) {
  fix_draws(object, newdata, ..., func = tidybayes::residual_draws)
}

fix_draws = function(object, newdata, ..., func = tidybayes::predicted_draws) {
  draws = func(object, newdata, ...)

  n = names(draws)

  dplyr::full_join(
    draws |> dplyr::select(-.chain, -.iteration),
    tidybayes::tidy_draws(object) |>
      dplyr::select(.chain, .iteration, .draw),
    by = ".draw"
  ) |>
    dplyr::select(dplyr::all_of(n)) |>
    dplyr::ungroup()
}
