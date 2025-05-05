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

  current_state <- initial_values
  for (t in seq(0, time_end, by = increment)) {
    step_result <- PBSddesolve::dde(
      y = current_state,
      times = c(t, t + increment),
      func = .ode_model_reference,
      parms = params
    )

    current_state <- as.numeric(step_result[nrow(step_result), -1])  # Exclude time column

    names(current_state) <- names(initial_values)

    print(jsonlite::toJSON(
      c(lapply(as.list(current_state), function(x) unname(x)), time = unname(t)),
      pretty = FALSE,
      auto_unbox = TRUE
    ))


    Sys.sleep(0.1)  # Add a delay of 1 second after each iteration
  }
}

.ode_model_reference <- function(t, current_state, params) {
  with(as.list(c(current_state, params)), {
    # ODEs
    dS <- -transmission_rate * S * I / population_size
    dE <- transmission_rate * S * I / population_size - infectiousness_rate * E
    dI <- infectiousness_rate * E - recovery_rate * I
    dR <- recovery_rate * I
    # output
    return(list(c(dS, dE, dI, dR)))
  })
}
