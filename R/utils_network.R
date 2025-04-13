#' F(tau) function
#'
#' @param transmission_rate
#' @param infectiousness_rate
#' @param recovery_rate
#' @param tau
#' @param infection
#'
#' @return
#' @export
#'
#' @examples
.kernel_capital_curly_f <- function(infection = c("SIR", "SEIR"),
                                    transmission_rate,
                                    infectiousness_rate = NULL,
                                    recovery_rate, tau) {
  sum1 <- transmission_rate + recovery_rate
  if (infection == "SIR") {
    fun <- recovery_rate / sum1 + transmission_rate / sum1 * exp(-sum1 * tau)
  } else if (infection == "SEIR") {
    # readability
    sum2 <- sum1 - infectiousness_rate
    fun <- (recovery_rate - infectiousness_rate) / sum2 + (transmission_rate * infectiousness_rate) / (sum1 * sum2) * exp(-(sum1 * tau)) + transmission_rate / sum2 * exp(-infectiousness_rate * tau)
  } else {
    stop("not implemented")
  }
  return(fun)
}


#' Derivative F'(tau)
#'
#' @param transmission_rate
#' @param infectiousness_rate
#' @param recovery_rate
#' @param tau
#' @param infection
#'
#' @return
#' @export
#'
#' @examples
.kernel_capital_curly_f_derivative <- function(infection = c("SIR", "SEIR"),
                                               transmission_rate, infectiousness_rate = NULL, recovery_rate, tau) {
  sum1 <- transmission_rate + recovery_rate
  if (infection == "SIR") {
    fun <- transmission_rate * exp(-sum1 * tau)
  } else if (infection == "SEIR") {
    sum2 <- sum1 - infectiousness_rate
    fun <- transmission_rate * infectiousness_rate / sum2 * (exp(-sum1 * tau) - exp(-infectiousness_rate * tau))
  } else {
    stop("not implemented")
  }
  return(fun)
}

#' Title
#'
#' @param infection
#' @param transmission_rate
#' @param infectiousness_rate
#' @param recovery_rate
#'
#' @return
#' @export
#'
#' @examples
.kernel_captial_curly_f_inf_limit <- function(infection = c("SIR", "SEIR"), transmission_rate, infectiousness_rate = NULL, recovery_rate) {
  beta <- transmission_rate
  gamma <- recovery_rate
  if (infection == "SIR" | infection == "SEIR") {
    lim <- gamma / (beta + gamma)
  } else {
    stop("not implemented")
  }
  return(lim)
}


#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
probability_generating_function <- function(x, ...) {
  degree_params <- list(...) |>
    unlist(recursive = FALSE)
  stopifnot("degree_distribution" %in% names(degree_params))
  if (degree_params$degree_distribution == "poisson") {
    stopifnot("lambda" %in% names(degree_params))
    lambda <- degree_params$lambda
    fun <- exp(lambda * (x - 1))
  } else if (degree_params$degree_distribution == "negative_binomial") {
    r <- degree_params$size
    p <- degree_params$prob
    fun <- (p / (1 - (1 - p) * x))^r
  } else if (degree_params$degree_distribution == "constant") {
    n <- degree_params$n
    fun <- x^n
  } else if (degree_params$degree_distribution == "geometric") {
    p <- degree_params$p
    fun <- p / (1 - (1 - p) * x)
  } else {
    stop("to be implemented")
  }
  return(fun)
}


#' Title
#'
#' @param seed_infected
#' @param degree_params
#'
#' @return
#' @export
#'
#' @examples
initial_x <- function(seed_infected, degree_params) {
  if (degree_params$degree_distribution == "poisson") {
    x <- 1 + log(1 - seed_infected) / degree_params$lambda
  } else if (degree_params$degree_distribution == "negative_binomial") {
    r <- degree_params$size
    p <- degree_params$prob
    x <- 1 - p * (1 - seed_infected)^(-1 / r)
    x <- x / (1 - p)
  } else if (degree_params$degree_distribution == "constant") {
    n <- degree_params$n
    x <- (1 - seed_infected)^(1 / n)
  } else if (degree_params$degree_distribution == "geometric") {
    p <- degree_params$p
    x <- (1 - seed_infected) - p
    x <- x / ((1 - seed_infected) * (1 - p))
  } else {
    stop("to be implemented")
  }
  return(x)
}


#' G'(x) / G'(1) = G'(x) / <k>, <k> = mean degree, G(x) probability generating function
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
helper_g_function <- function(x, ...) {
  degree_params <- list(...) |>
    unlist(recursive = FALSE)
  stopifnot("degree_distribution" %in% names(degree_params))
  if (degree_params$degree_distribution == "poisson") {
    stopifnot("lambda" %in% names(degree_params))
    lambda <- degree_params$lambda
    fun <- probability_generating_function(x = x, degree_params)
  } else if (degree_params$degree_distribution == "negative_binomial") {
    r <- degree_params$size
    p <- degree_params$prob
    fun <- (p / (1 - (1 - p) * x))^(r + 1)
  } else if (degree_params$degree_distribution == "constant") {
    n <- degree_params$n
    fun <- x^(n - 1)
  } else if (degree_params$degree_distribution == "geometric") {
    p <- degree_params$p
    fun <- p^2 / (1 - (1 - p) * x)^2
  } else {
    stop("to be implemented")
  }
  return(fun)
}


#' Derivative of [helper_g_function]
#' G''(x) / G'(1) = G''(x) / <k>, <k> = mean degree, G(x) probability generating function
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
helper_g_function_derivative <- function(x, ...) {
  degree_params <- list(...) |>
    unlist(recursive = FALSE)
  stopifnot("degree_distribution" %in% names(degree_params))
  if (degree_params$degree_distribution == "poisson") {
    stopifnot("lambda" %in% names(degree_params))
    lambda <- degree_params$lambda
    fun <- lambda * exp(lambda * (x - 1))
  } else if (degree_params$degree_distribution == "negative_binomial") {
    r <- degree_params$size
    p <- degree_params$prob
    fun <- (1 - p) * p * (1 + r) * (p / (1 - (1 - p) * x))^r
    fun <- fun / (1 - (1 - p) * x)^2
  } else if (degree_params$degree_distribution == "constant") {
    n <- degree_params$n
    fun <- (n - 1) * x^(n - 2)
  } else if (degree_params$degree_distribution == "geometric") {
    p <- degree_params$p
    fun <- 2 * (1 - p) * p^2 / (1 + (1 - p) * x)^3
  } else {
    stop("to be implemented")
  }
  return(fun)
}


#' Title
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
mean_degree <- function(...) {
  degree_params <- list(...) |>
    unlist(recursive = FALSE)
  stopifnot("degree_distribution" %in% names(degree_params))
  if (degree_params$degree_distribution == "poisson") {
    stopifnot("lambda" %in% names(degree_params))
    lambda <- degree_params$lambda
    avg <- degree_params$lambda
  } else if (degree_params$degree_distribution == "negative_binomial") {
    r <- degree_params$size
    p <- degree_params$prob
    avg <- r * (1 - p) / p
  } else if (degree_params$degree_distribution == "constant") {
    avg <- degree_params$n
  } else if (degree_params$degree_distribution == "geometric") {
    p <- degree_params$p
    avg <- (1 - p) / p
  } else {
    stop("to be implemented")
  }
  return(avg)
}


#' Title
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
var_degree <- function(...) {
  degree_params <- list(...) |>
    unlist(recursive = FALSE)
  stopifnot("degree_distribution" %in% names(degree_params))
  if (degree_params$degree_distribution == "poisson") {
    stopifnot("lambda" %in% names(degree_params))
    var <- degree_params$lambda
  } else if (degree_params$degree_distribution == "negative_binomial") {
    r <- degree_params$size
    p <- degree_params$prob
    var <- r * (1 - p) / p^2
  } else if (degree_params$degree_distribution == "constant") {
    var <- 0
  } else if (degree_params$degree_distribution == "geometric") {
    p <- degree_params$p
    var <- (1 - p) / p^2
  } else {
    stop("to be implemented")
  }
  return(var)
}
