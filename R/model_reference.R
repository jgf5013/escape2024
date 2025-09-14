#' Title
#'
#' @param transmission_rate
#' @param infectiousness_rate
#' @param recovery_rate
#' @param time_end
#' @param increment
#' @param population_size
#' @param seed_infected
#'
#' @returns
#' @export
#'
#' @examples
model_reference <- function(
    transmission_rate = 0.25,
    infectiousness_rate = 0.25,
    recovery_rate = 0.2,
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

  jacobian <- matrix(
    c(
      -infectiousness_rate, transmission_rate, 0,
      infectiousness_rate, -recovery_rate, 0,
      0, recovery_rate, 0
    ),
    byrow = TRUE,
    nrow = 3
  )
  dfe <- rep(0, 3)
  initial_values <- calculate_initial_values(jacobian, dfe, 3, 1, population_size * seed_infected)
  names(initial_values) <- c("E", "I", "R")
  return(initial_values)
}


simulate_outbreak_seir_reference <- function(t, increment, current_state, params) {
    step_result <- PBSddesolve::dde(
      y = current_state,
      times = c(t, t + increment),
      func = .ode_model_reference,
      parms = params
    )
    current_state <- as.numeric(step_result[nrow(step_result), -1])  # Exclude time column
    names(current_state) <- c("E", "I", "R", "S", "incidence")

    print(jsonlite::toJSON(
      list(
        state = as.list(current_state),
        time = unname(t),
        model_type = "model_reference"
      ),
      pretty = FALSE,
      auto_unbox = TRUE
    ))
    current_state <- current_state[c("E", "I", "R")]
}


#' Title
#'
#' @param t
#' @param current_state
#' @param params
#'
#' @returns
#' @export
#'
#' @examples
.ode_model_reference <- function(t, current_state, params) {
  with(as.list(c(current_state, params)), {
    S <- population_size - E - I - R
    # ODEs
    dE <- transmission_rate * S * I / population_size - infectiousness_rate * E
    dI <- infectiousness_rate * E - recovery_rate * I
    dR <- recovery_rate * I
    # output
    return(list(
      c(dE, dI, dR),
      c(
        S = S,
        incidence = transmission_rate * S * I / population_size
      )
    ))
  })
}
