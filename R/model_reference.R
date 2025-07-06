model_reference <- function(
    simulation_id,
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
  # eigensystem of linearization in disease free state
  es <- eigen(jacobian)

  if (sum(Re(es$values) > 0) == 1) {
    # if there is exactly one unstable eigenvector
    es1 <- Re(es$vectors[, Re(es$values) > 0])
    # normalize eigenvector such that first element is -1
    es1 <- es1 / es1[1]
  } else if (sum(Re(es$values) > 0) > 1) {
    es1 <- Re(es$vectors[, Re(es$values) > 0][, 1])
    # normalize eigenvector such that first element is -1
    es1 <- es1 / es1[1]
  } else {
    # unstable eigenvector is set to zero
    es1 <- rep(0, 3)
  }

  # perturb DFE with unstable eigenvector
  initial_values <- dfe + population_size * seed_infected * es1
  names(initial_values) <- c("E", "I", "R")

  current_state <- initial_values

  for (t in seq(0, time_end, by = increment)) {
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
}

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
