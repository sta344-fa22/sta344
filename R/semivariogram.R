#' @name sv_functions
#'
#' @title Semivariogram functions
#'
#' @description Functions for calculating various semivariograms
#'
#' @inheritParams covariance_functions
#'
NULL

#' @rdname sv_functions
#' @export
exp_sv  = function(d, sigma2=1, l=1, sigma2_w=0) {
  sigma2 + sigma2_w - exp_cov(d,sigma2,l) - nugget_cov(d,sigma2_w)
}

#' @rdname sv_functions
#' @export
sq_exp_sv  = function(d, sigma2=1, l=1, sigma2_w=0) {
  sigma2 + sigma2_w - sq_exp_cov(d,sigma2,l) - nugget_cov(d,sigma2_w)
}

#' @rdname sv_functions
#' @export
pow_exp_sv = function(d, sigma2=1, l=1, p=2, sigma2_w=0) {
  sigma2 + sigma2_w - pow_exp_cov(d,sigma2,l,p) - nugget_cov(d,sigma2_w)
}


dist_mat = function(d) {
  as.matrix(stats::dist(d))
}

dist_long = function(d) {
  d = as.matrix(d)
  d[upper.tri(d, diag = TRUE)] = NA

  tidyr::expand_grid(
    i=1:nrow(d),
    j=1:nrow(d)
  ) |>
    dplyr::mutate(
      dist = c(d)
    ) |>
    dplyr::filter(!is.na(dist))
}

bin = function(df, var, binwidth, start = NULL, end = NULL) {
  n = nrow(df)

  var = as.character(substitute(var))
  x = df[[var]]

  if (is.null(start)) {
    start = min(x) - (min(x) %% binwidth)
  }
  if (is.null(end)) {
    end = max(x) + binwidth - (max(x) %% binwidth)
  }

  bins = seq(start, end, by = binwidth)

  df |>
    dplyr::mutate(
      bins = cut(.data[[var]], breaks = bins),
      bin_mid = .data[[var]] - .data[[var]] %% binwidth + binwidth/2
    ) |>
    dplyr::group_by(bins)
}

#' @title Estimate the empirical semivariogram
#'
#' @description Calculates the empirical semivariogram for a given
#' binwidth.
#'
#' @param d Data frame containing data
#' @param y Unquoted variable name for y values
#' @param x Unquoted variable name for x values
#' @param bin Logical. Should variance estimatess be binned.
#' @param binwidth Numeric. The width of the bins to use.
#' @param range_max  Numeric. Filter bins with a midpoint above this value.
#' Defaults to 1/3 of the max distance between `x` values.
#'
#' @export
#'
emp_semivariogram = function(d, y, x, bin=FALSE, binwidth, range_max=NULL) {
  y_col = as.character(substitute(y))
  x_col = as.character(substitute(x))

  d = d[[x_col]] |>
    dist() |>
    dist_long() |>
    mutate(
      y_i = d[[y_col]][i],
      y_j = d[[y_col]][j]
    )

  if (is.null(range_max))
    range_max = max(d$dist, na.rm=TRUE)/3

  if (!bin) {
    return( d |> dplyr::mutate(bin_mid = dist) |> dplyr::rowwise() )
  }

  purrr::map_dfr(
    binwidth,
    function(bw) {

      res = d |>
        bin(dist, binwidth = bw, start = 0) |>
        dplyr::summarize(
          gamma = sum( (y_i - y_j)^2 / (2*n()) ),
          d = mean(bin_mid),
          n = n()
        ) |>
        dplyr::filter(d < range_max)

      res |>
        mutate(bw = bw)
    }
  )
}
