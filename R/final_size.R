#' Title
#'
#' @param .data
#' @param precision
#'
#' @return
#' @export
#'
#' @examples
final_size <- function(.data, precision = 1E-6) {
  # check that epidemic has run its course
  if (abs(.data$S[length(.data$S) - 10] - min(.data$S)) > precision) {
    warning(
      paste0(
        "susceptible depletion has not yet converged within ",
        precision,
        ". You may want to increase the end time of the epidemic simulation"
      )
    )
  }

  # get model parameters
  params <- attr(.data, "parameters")
  # fraction of the population that got infected in the epidemic
  final <- (max(.data$S) - min(.data$S)) / params$population_size
  return(final)
}
