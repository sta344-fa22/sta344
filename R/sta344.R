#' @importFrom rlang ".data"
NULL


#' Average weekly temperatures in Durham, NC
#'
#' Data collected from the Darksky weather API and aggregated to weekly averages.
#'
#' @format A data frame with 156 rows and 3 columns:
#'
#' * `date` - date representing each week (Wednesday)
#' * `avg_temp` - average weekly temperature in degrees Farhenheight
#' * `week` - numerical index of week from the start of the data (0-based)
#'
"avg_temp"
