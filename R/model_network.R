#' Title
#'
#' @param degree_distribution
#' @param infection
#' @param transmission_rate
#' @param infectiousness_rate
#' @param recovery_rate
#' @param time_end
#' @param increment
#' @param population_size
#' @param seed_infected
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
model_network <- function(
    degree_distribution = c("poisson", "negative_binomial", "constant", "geometric"),
    infection = c("SIR", "SEIR"),
    transmission_rate = 1.25,
    infectiousness_rate = 0.25,
    recovery_rate = 0.2,
    time_end = 400,
    increment = 1,
    population_size = 20E6,
    seed_infected = 1E-3, ...) {
  args <- list(...)
  if (!degree_distribution %in% c("poisson", "negative_binomial", "constant", "geometric")) {
    stop("to be implemented")
  }
  if (infection == "SIR") {
    data <- model_network_sir(
      degree_distribution = degree_distribution,
      transmission_rate = transmission_rate,
      recovery_rate = recovery_rate,
      time_end = time_end,
      increment = increment,
      population_size = population_size,
      seed_infected = seed_infected,
      args
    )
  } else if (infection == "SEIR") {
    data <- model_network_seir(
      degree_distribution = degree_distribution,
      transmission_rate = transmission_rate,
      infectiousness_rate = infectiousness_rate,
      recovery_rate = recovery_rate,
      time_end = time_end,
      increment = increment,
      population_size = population_size,
      seed_infected = seed_infected,
      args
    )
  } else {
    stop("to be implemented")
  }
  return(data)
}


#' Title
#'
#' @param degree_distribution
#' @param transmission_rate
#' @param recovery_rate
#' @param time_end
#' @param increment
#' @param population_size
#' @param seed_infected
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
model_network_sir <- function(degree_distribution = c("poisson", "negative_binomial", "constant", "geometric"),
                              transmission_rate = 0.25,
                              recovery_rate = 0.25,
                              time_end = 400,
                              increment = 1,
                              population_size = 20E6,
                              seed_infected = 1E-3, ...) {
  degree_params <- list(...) |>
    unlist(recursive = FALSE)
  if (degree_distribution == "poisson") {
    stopifnot("lambda" %in% names(degree_params))
  } else if (degree_distribution == "negative_binomial") {
    stopifnot(
      "size" %in% names(degree_params),
      "prob" %in% names(degree_params)
    )
  } else if (degree_distribution == "constant") {
    stopifnot(
      "n" %in% names(degree_params)
    )
  } else if (degree_distribution == "geometric") {
    stopifnot(
      "p" %in% names(degree_params)
    )
  } else {
    stop("to be implemented")
  }

  degree_params <- c(degree_distribution = degree_distribution, degree_params)

  transmission_params <- list(
    transmission_rate = transmission_rate,
    recovery_rate = recovery_rate,
    population_size = population_size,
    seed_infected = seed_infected
  )

  params <- c(
    transmission_params,
    degree_params
  )
  initial_values <- c(
    xbar = initial_x(seed_infected, degree_params),
    R = 0
  )

  # simulate outbreak
  current_state <- initial_values
  for (t in seq(0, time_end, by = increment)) {
    step_result <- PBSddesolve::dde(
      y = initial_values,
      times = 0:time_end,
      func = .ode_model_network_sir,
      parms = params
    )
    current_state <- as.numeric(step_result[nrow(step_result), -1])  # Exclude time column
    names(current_state) <- c("xbar", "R", "incidence", "S", "I")
    current_state[c("S", "I", "R", "incidence")] <- transmission_params$population_size * current_state[c("S", "E", "I", "R", "incidence")]

    print(jsonlite::toJSON(
      list(state = as.list(current_state[c("S", "I", "R", "incidence")]), time = unname(t)),
      pretty = TRUE,
      auto_unbox = TRUE
    ))
    current_state <- current_state[c("xbar", "R")]
  }
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
.ode_model_network_sir <- function(t, current_state, params) {
  with(as.list(c(current_state, params)), {
    dxbar <- transmission_rate * helper_g_function(x = xbar, params) - (transmission_rate + recovery_rate) * xbar + recovery_rate
    dR <- recovery_rate * (1 - probability_generating_function(x = xbar, params) - R)
    return(list(
      c(dxbar, dR),
      c(incidence = -helper_g_function(x = xbar, params) * mean_degree(params) * dxbar,
        S = probability_generating_function(x = xbar, params),
        I = 1 - probability_generating_function(x = xbar, params) - R)
    ))
  })
}


#' Title
#'
#' @param degree_distribution
#' @param transmission_rate
#' @param infectiousness_rate
#' @param recovery_rate
#' @param time_end
#' @param increment
#' @param population_size
#' @param seed_infected
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
model_network_seir <- function(degree_distribution = c("poisson", "negative_binomial", "constant", "geometric"),
                               transmission_rate = 0.25,
                               infectiousness_rate = 0.1,
                               recovery_rate = 0.25,
                               time_end = 400,
                               increment = 1,
                               population_size = 20E6,
                               seed_infected = 1E-3, ...) {
  degree_params <- list(...) |>
    unlist(recursive = FALSE)
  if (degree_distribution == "poisson") {
    stopifnot("lambda" %in% names(degree_params))
  } else if (degree_distribution == "negative_binomial") {
    stopifnot(
      "size" %in% names(degree_params),
      "prob" %in% names(degree_params)
    )
  } else if (degree_distribution == "constant") {
    stopifnot(
      "n" %in% names(degree_params)
    )
  } else if (degree_distribution == "geometric") {
    stopifnot(
      "p" %in% names(degree_params)
    )
  } else {
    stop("to be implemented")
  }

  degree_params <- c(degree_distribution = degree_distribution, degree_params)

  transmission_params <- list(
    transmission_rate = transmission_rate,
    infectiousness_rate = infectiousness_rate,
    recovery_rate = recovery_rate,
    population_size = population_size,
    seed_infected = seed_infected
  )

  f_inf <- .kernel_captial_curly_f_inf_limit(
    infection = "SEIR",
    transmission_rate,
    infectiousness_rate,
    recovery_rate
  )

  params <- c(
    transmission_params,
    degree_params,
    f_inf = f_inf
  )

  # initial value specified as fraction of susceptibles that are infected (randomly in network)
  # initial value specified as fraction of susceptibles that are infected (randomly in network)
  x0 <- initial_x(seed_infected, degree_params)

  # linearization in disease free steady state
  jacobian <- matrix(
    c(
      0, 0, 0, -transmission_rate, 0, 0, 0,
      0, 0, 0, -transmission_rate * helper_g_function_derivative(1, degree_params), 0, 0, 0,
      0, 0, -infectiousness_rate, transmission_rate * helper_g_function_derivative(1, degree_params), 0, 0, 0,
      0, 0, infectiousness_rate, -(transmission_rate + recovery_rate), 0, 0, 0,
      0, 0, 0, transmission_rate * helper_g_function(1, degree_params) * mean_degree(degree_params), -infectiousness_rate, 0, 0,
      0, 0, 0, 0, infectiousness_rate, -recovery_rate, 0,
      0, 0, 0, 0, 0, recovery_rate, 0
    ),
    byrow = TRUE, nrow = 7
  )

  # disease free steady state
  dfe <- c(1, 1, 0, 0, 0, 0, 0)
  # eigensystem of linearization in disease free state
  es <- eigen(jacobian)

  if (sum(es$values > 0) == 1) {
    # if there is exactly one unstable eigenvector
    es1 <- es$vectors[, es$values > 0]
  } else if (sum(es$values > 0) > 1) {
    es1 <- es$vectors[, es$values > 0][, 1]
  } else {
    # unstable eigenvector is set to zero
    es1 <- c(0, 0, 0, 0)
  }
  # normalize eigenvector such that first element is -1
  es1 <- -es1 / es1[1]
  # perturb DFE with unstable eigenvector
  perturb <- dfe + (1 - x0) * es1
  # ensure that xbar and xS remain between 0 and 1

  initial_values <- perturb
  names(initial_values) <- c("xbar", "xS", "xE", "xI", "E", "I", "R")

  # simulate outbreak
  # simulate outbreak
  current_state <- initial_values
  for (t in seq(0, time_end, by = increment)) {
    step_result <- PBSddesolve::dde(
      y = initial_values,
      times = 0:time_end,
      func = .ode_model_network_seir,
      parms = params
    )

    current_state <- as.numeric(step_result[nrow(step_result), -1])  # Exclude time column
    names(current_state) <- c("z1", "z2", "E", "I", "R", "x", "S", "incidence")
    current_state[c("S", "E", "I", "R", "incidence")] <- transmission_params$population_size * current_state[c("S", "E", "I", "R", "incidence")]

    print(jsonlite::toJSON(
      list(state = as.list(current_state[c("S", "E", "I", "R", "incidence")]), time = unname(t)),
      pretty = TRUE,
      auto_unbox = TRUE
    ))
    current_state <- current_state[c("z1", "z2", "E", "I", "R")]
  }
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
.ode_model_network_seir <- function(t, current_state, params) {
  with(as.list(c(current_state, params)), {
    if (transmission_rate > 0) {
      dxbar <- -transmission_rate * xI
      dxS <- -transmission_rate * helper_g_function_derivative(x = xbar, params) * xI
      dxE <- transmission_rate * helper_g_function_derivative(x = xbar, params) * xI - infectiousness_rate * xE
      dxI <- infectiousness_rate * xE - (transmission_rate + recovery_rate) * xI

      S <- probability_generating_function(x = xbar, params)

      dE <- -helper_g_function(x = xbar, params) * mean_degree(params) * dxbar - infectiousness_rate * E
      dI <- infectiousness_rate * E - recovery_rate * I
      dR <- recovery_rate * I
    } else if (transmission_rate == 0) {
      dx <- 0
      dz1 <- 0
      dS <- 0
      dE <- 0
      dI <- 0
      dR <- 0
    }
    return(list(
      c(dz1, dz2, dE, dI, dR),
      c(
        x = x,
        S = S,
        incidence = -helper_g_function(x = x, params) * mean_degree(params) * dx
      )
    ))
  })
}
