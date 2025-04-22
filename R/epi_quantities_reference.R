#' Title
#'
#' @param .params
#'
#' @return
#' @export
#'
#' @examples
theory_reference_reproduction_number <- function(.params) {
  return(.params$transmission_rate / .params$recovery_rate)
}


#' Title
#'
#' @param .params
#'
#' @return
#' @export
#'
#' @examples
theory_reference_growth_rate <- function(.params) {
  # convenience abbreviations for readability of code
  beta <- .params$transmission_rate
  sigma <- .params$infectiousness_rate
  gamma <- .params$recovery_rate
  # growth rate is largest root of quadratic equation
  growth <- (-(sigma + gamma) + sqrt((sigma + gamma)^2 - 4 * ((gamma - beta) * sigma))) / 2
  # check that R_0 > 1 iff r > 0, R_0 < 1 iff r < 0, R_0 = 1 iff r=0
  r_0 <- theory_reference_reproduction_number(.params)
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
#' @param .fraction
#' @param r_0
#'
#' @return
#' @export
#'
#' @examples
theory_reference_final_size <- function(.params = NULL, r_0 = NULL, .fraction = TRUE) {
  if (!is.null(.params) & is.null(r_0)) {
    r_0 <- theory_reference_reproduction_number(.params = .params)
    # initial fraction susceptible
    s_0 <- 1 - .params$seed_infected
  } else if ((!is.null(.params) & !is.null(r_0)) |
             (is.null(.params) & is.null(r_0))) {
    stop("provide either parameters or R_0, not both")
  }

  if (r_0 < 1) { # trivial solution
    return(0)
  } else {
    # initial fraction susceptible when only r_0 is provided
    s_0 <- 1
    # helper function
    .root_fun <- function(s) {
      # s(t) = 1-S(t)/S_0, S0 = initial number susceptible population
      log(s) - log(s_0) - r_0 * (s - 1)
    }
    # hack: use base uniroot and avoid trivial solution 1 from root finding
    final <- uniroot(.root_fun, interval = c(0, 1-1E-6), tol = .Machine$double.eps)$root
    # final infected fraction s(0) - s(oo)
    final <- s_0 - final
    if (!.fraction) { # convert into population number
      final <- .params$population_size * final
    }
    return(final)
  }
}
