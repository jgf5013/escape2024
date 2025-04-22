#' Title
#'
#' @param r_0
#' @param transmission_rate
#' @param infectiousness_rate
#' @param recovery_rate
#' @param time_end
#' @param increment
#' @param population_size
#' @param seed_infected
#'
#' @return
#' @export
#'
#' @examples
model_reference <- function(
    transmission_rate = 0.25,
    infectiousness_rate = 0.25,
    recovery_rate = 0.2,
    time_end = 400,
    increment = 1,
    population_size = 20E6,
    seed_infected = 1E-3) {
  # setup parameters
  params <- list(
    transmission_rate = transmission_rate,
    infectiousness_rate = infectiousness_rate,
    recovery_rate = recovery_rate,
    population_size = population_size,
    seed_infected = seed_infected
  )

  # setup initial values
  initial_values <- c(
    S = population_size * (1 - seed_infected),
    E = population_size * seed_infected,
    I = 0,
    R = 0
  )

  # simulate outbreak
  seir_out <- PBSddesolve::dde(
    y = initial_values,
    times = seq(0, time_end, by = increment),
    func = .ode_model_reference,
    parms = params
  )
  # convert class of output and type of
  seir_out <- seir_out |>
    as.data.frame() |>
    # change deSolve classes to numeric
    type.convert(as.is = TRUE)

  # add parameters as attribute
  attr(seir_out, "parameters") <- params
  attr(seir_out, "infection") <- "SEIR"

  return(seir_out)
}


#' Title
#'
#' @param t
#' @param current_state
#' @param params
#'
#' @return
#' @export
#'
#' @examples
.ode_model_reference <- function(t, current_state, params) {
  with(as.list(c(current_state, params)), {
    # ODEs
    dS <- -transmission_rate * S * I / population_size
    dE <- transmission_rate * S * I / population_size - infectiousness_rate * E
    dI <- infectiousness_rate * E - recovery_rate * I
    dR <- recovery_rate * I
    # output
    return(list(c(dS, dE, dI, dR),
      incidence = transmission_rate * S * I / population_size
    ))
  })
}
