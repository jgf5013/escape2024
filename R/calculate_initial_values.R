calculate_initial_values <- function(jacobian, disease_free_state, dimension, dim_exposed, seed_infected) {
  # eigensystem of linearization in disease free state
  es <- eigen(jacobian)
  if (sum(Re(es$values) > 0) == 1) {
    # if there is exactly one unstable eigenvector
    es1 <- Re(es$vectors[, Re(es$values) > 0])
  } else if (sum(Re(es$values) > 0) > 1) {
    es1 <- Re(es$vectors[, Re(es$values) > 0][, 1])
  } else {
    # unstable eigenvector is set to zero
    es1 <- rep(0, dimension)
  }
  # normalize eigenvector such that first element is -1
  es1 <- es1 / es1[dim_exposed]
  # perturb DFE with unstable eigenvector
  perturb <- disease_free_state + seed_infected * es1
  return(perturb)
}
