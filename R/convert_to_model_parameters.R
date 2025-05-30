#' Title
#'
#' @param serial_interval
#' @param r0
#'
#' @returns
#' @export
#'
#' @examples
convert_to_model_parameters <- function(serial_interval, r0) {
  ## use params to calculate latent period, infectious period
  infectiousness_rate <- 2 / serial
  recovery_rate <- 2 / serial
  transmission_rate <- r0 * recovery_rate
  return(list(infectiousness_rate = infectiousness_rate,
              recovery_rate = recovery_rate,
              transmission_rate = transmission_rate))
}


#' Title
#'
#' @param serial_interval
#' @param r0
#' @param time_end
#' @param population_size
#' @param seed_infected
#' @param increment
#'
#' @returns
#' @export
#'
#' @examples
model_reference_with_epi_parms <- function(serial_interval, r0, time_end, population_size, seed_infected, increment) {
  parms <- convert_to_model_parameters(serial_interval, r0)
  data <- model_reference(
    transmission_rate = parms$transmission_rate,
    infectiousness_rate = parms$infectiousness_rate,
    recovery_rate = parms$recovery_rate,
    time_end = time_end,
    population_size = population_size,
    seed_infected = seed_infected,
    increment = increment
  )
  return(data)
}
