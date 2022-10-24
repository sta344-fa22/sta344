#' @export
strip_attrs = function(obj) {
  attributes(obj) = NULL
  obj
}

#' @export
normalize_weights = function(w) {
  diag(w) = 0
  rs = rowSums(w)
  rs[rs == 0] = 1
  w/rs
}
