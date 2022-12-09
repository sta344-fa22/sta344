#' @name covariance_functions
#'
#' @title Covariance functions
#'
#' @description Functions for calculating various covariances
#'
#' @param d Distance matrix
#' @param sigma2 Scale variance parameter
#' @param l Inverse length scale parameter
#' @param sigma2_w Nugget variance parameter
#' @param p Power parameter (powered exponential covariance)
#' @param a Alpha parameter (rational quartic covariance)
#' @param per Period parameter (periodic covariance)
#' @param nu Nu parameter (Matern covariance)
#'
NULL

#' @rdname covariance_functions
#' @export
nugget_cov = function(d, sigma2=1) {
  ifelse(d==0, sigma2, 0)
}

#' @rdname covariance_functions
#' @export
exp_cov = function(d, sigma2=1, l=1, sigma2_w=0) {
  sigma2 * exp(-abs(d)*l) + nugget_cov(d,sigma2_w)
}

#' @rdname covariance_functions
#' @export
sq_exp_cov = function(d, sigma2=1, l=1, sigma2_w=0) {
  sigma2 * exp(-(abs(d)*l)^2) + nugget_cov(d,sigma2_w)
}

#' @rdname covariance_functions
#' @export
pow_exp_cov = function(d, sigma2=1, l=1, p=2, sigma2_w=0) {
  sigma2 * exp(-(abs(d)*l)^p) + nugget_cov(d,sigma2_w)
}

#' @rdname covariance_functions
#' @export
rquad_cov = function(d, sigma2=1, l=1, a=1) {
  sigma2 * (1+d^2*l^2/a)^(-a)
}

#' @rdname covariance_functions
#' @export
periodic_cov = function(d, sigma2=1, l=1, per=1) {
  sigma2 * exp(-2*l^2*sin(pi*d/per)^2)
}

#' @rdname covariance_functions
#' @export
matern_cov = function(d, sigma2=1, l=1, nu=1/2) {
  fields::Matern(d, alpha=l, nu=nu, phi=sigma2)
}

#' @rdname covariance_functions
#' @export
sphere_cov = function(d, sigma2=1, l=1) {
  ifelse(d > 1/l, 0, sigma2*(1 - 1.5*d*l + 0.5*(d*l)^3))
}

#' @export
linear_cov   = function(x1, x2, sigma2_b=1, sigma2_v=1, c=0) {
  if (!missing(x2)) {
    sigma2_b + sigma2_v * (x1-c) * (x2-c)
  } else {
    d = expand.grid(t_i=x1, t_j=x1)

    linear_cov(d[[1]], d[[2]], sigma2_b=sigma2_b, sigma2_v=sigma2_v, c=c) |>
      matrix(ncol=length(x1))
  }
}











