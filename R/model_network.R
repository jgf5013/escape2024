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
model_network_sir <- function(
    degree_distribution = c("poisson", "negative_binomial", "constant", "geometric"),
    transmission_rate = 0.25,
    recovery_rate = 0.25,
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

  degree_params <- c(list(degree_distribution = degree_distribution), degree_params)

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
  return(current_state)
}


#' Title
#'
#' @param time_end
#' @param increment
#' @param current_state
#' @param params
#'
#' @returns
#' @export
#'
#' @examples
simulate_outbreak_sir_network <- function(t, increment, current_state, params) {
  step_result <- PBSddesolve::dde(
    y = current_state,
    times = c(t, t + increment),
    func = .ode_model_network_sir,
    parms = params
  )
  current_state <- as.numeric(step_result[nrow(step_result), -1]) # Exclude time column
  names(current_state) <- c("xbar", "R", "incidence", "S", "I")

  print(jsonlite::toJSON(
    list(
      state = as.list(params$population_size * current_state[c("S", "I", "R", "incidence")]),
      time = unname(t),
      # TODO model_reference, model_network + degree
      model_type = paste0("model_network", params$degree_distribution)
    ),
    pretty = FALSE,
    auto_unbox = TRUE
  ))
  current_state <- current_state[c("xbar", "R")]
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
      c(
        incidence = -helper_g_function(x = xbar, params) * mean_degree(params) * dxbar,
        S = probability_generating_function(x = xbar, params),
        I = 1 - probability_generating_function(x = xbar, params) - R
      )
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

  degree_params <- c(list(degree_distribution = degree_distribution), degree_params)
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

  initial_values <- calculate_initial_values(jacobian, dfe, 7, 5, seed_infected)
  names(initial_values) <- c("xbar", "xS", "xE", "xI", "E", "I", "R")
  return(initial_values)
}


#' Title
#'
#' @param time_end
#' @param increment
#' @param current_state
#' @param params
#'
#' @returns
#' @export
#'
#' @examples
simulate_outbreak_seir_network <- function(t, increment, current_state, params) {
  step_result <- PBSddesolve::dde(
    y = current_state,
    times = c(t, t + increment),
    func = .ode_model_network_seir,
    parms = params
  )
  current_state <- as.numeric(step_result[nrow(step_result), -1]) # Exclude time column
  names(current_state) <- c("xbar", "xS", "xE", "xI", "E", "I", "R", "S", "incidence")
  print(jsonlite::toJSON(
    list(
      state = as.list(params$population_size * current_state[c("S", "E", "I", "R", "incidence")]),
      time = unname(t),
      model_type = paste0("model_network_", params$degree_distribution)
    ),
    pretty = FALSE,
    auto_unbox = TRUE
  ))
  current_state <- current_state[c("xbar", "xS", "xE", "xI", "E", "I", "R")]
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
      dxbar <- 0
      dxS <- 0
      dxE <- 0
      dxI <- 0
      dE <- 0
      dI <- 0
      S <- probability_generating_function(x = xbar, params)
    }
    return(list(
      c(dxbar, dxS, dxE, dxI, dE, dI, dR),
      c(
        S = S,
        incidence = -helper_g_function(x = xbar, params) * mean_degree(params) * dxbar
      )
    ))
  })
}
