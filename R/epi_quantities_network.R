#' Title
#'
#' @param .infection
#' @param .params
#' @param .degree
#'
#' @return
#' @export
#'
#' @examples
theory_network_reproduction_number <- function(.infection, .params, .degree) {
  if (!.degree$degree_distribution %in% c("poisson", "negative_binomial", "constant", "geometric")) {
    stop("not implemented")
  }
  var <- var_degree(.degree)
  avg <- mean_degree(.degree)
  c_degree <- (var + avg^2 - avg) / avg

  r_0 <- c_degree * (1 - .kernel_captial_curly_f_inf_limit(
    .infection,
    transmission_rate = .params$transmission_rate,
    infectiousness_rate = .params$infectiousness_rate,
    recovery_rate = .params$recovery_rate
  ))
  return(r_0)
}


#' Title
#'
#' @param .infection
#' @param .params
#' @param .degree
#'
#' @return
#' @export
#'
#' @examples
theory_network_growth_rate <- function(.infection, .params, .degree) {
  # convenience abbreviations for readability of code
  beta <- .params$transmission_rate
  gamma <- .params$recovery_rate

  var <- var_degree(.degree)
  avg <- mean_degree(.degree)
  c_degree <- (var + avg^2 - avg) / avg

  if (.degree$degree_distribution %in% c("poisson", "negative_binomial", "constant", "geometric") & .infection == "SIR") {
    growth <- beta * c_degree - (beta + gamma)
  } else if (.infection == "SEIR") {
    sigma <- .params$infectiousness_rate
    # growth rate is positive root of quadratic equation
    bb <- beta + sigma + gamma
    cc <- -beta * sigma * (1 + c_degree) + gamma * sigma
    growth <- (-bb + sqrt(bb^2 - 4 * cc)) / 2
  } else {
    stop("to be implemented")
  }

  # check that R_0 > 1 iff r > 0, R_0 < 1 iff r < 0, R_0 = 1 iff r=0
  r_0 <- theory_network_reproduction_number(.params, .degree)
  if (r_0 > 1) {
    stopifnot(growth > 0)
  } else if (r_0 < 1) {
    stopifnot(growth < 0)
  } else { # r_0 = 1
    growth <- 0
  }
  return(growth)
}


#' Title
#'
#' @param .params
#' @param .degree
#' @param .infection
#' @param .fraction
#'
#' @return
#' @export
#'
#' @examples
theory_network_final_size <- function(.params, .degree, .infection, .fraction = TRUE) {
  r_0 <- theory_network_reproduction_number(.infection, .params, .degree)

  if (r_0 < 1) { # trivial solution
    final <- 0
    x_final <- 1
  } else if (.degree$degree_distribution %in% c("poisson", "negative_binomial", "constant", "geometric")) {
    # lim_{tau-> infinity} F(tau) as function of infection parameters
    f_inf <- .kernel_captial_curly_f_inf_limit(.infection,
                                               transmission_rate = .params$transmission_rate,
                                               infectiousness_rate = .params$infectiousness_rate,
                                               recovery_rate = .params$recovery_rate
    )
    # initial condition for x given seed infected
    x0 <- initial_x(.params$seed_infected, .degree)
    # helper function
    .root_fun <- function(x) {
      f_inf + (1 - f_inf) * helper_g_function(x = x, .degree) - x
    }
    # hack: use base uniroot and avoid trivial solution 1 from root finding
    x_final <- uniroot(.root_fun, interval = c(0, 1-1E-6), tol = .Machine$double.eps)$root
    # final number of susceptibles
    s_final <- probability_generating_function(x = x_final, .degree)
    # take into account initially infected fraction
    x_final <- x0 - x_final
    final <- (1 - .params$seed_infected) - s_final
  } else {
    stop("not (yet) implemented")
  }

  if (!.fraction) { # convert into population number
    final <- .params$population_size * final
  }

  return(list(
    xbar = x_final,
    final = final
  ))
}
