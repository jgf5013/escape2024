% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_network.R
\name{model_network}
\alias{model_network}
\title{Title}
\usage{
model_network(
  simulation_id,
  degree_distribution = c("poisson", "negative_binomial", "constant", "geometric"),
  infection = c("SIR", "SEIR"),
  transmission_rate = 1.25,
  infectiousness_rate = 0.25,
  recovery_rate = 0.2,
  time_end = 400,
  increment = 1,
  population_size = 2e+07,
  seed_infected = 0.001,
  ...
)
}
\arguments{
\item{...}{}
}
\description{
Title
}
