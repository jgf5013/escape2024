################################################################################
#
# Table 2 and Figure 4
#
################################################################################

library(dplyr)
library(ggplot2)
library(escape2024)

fs::dir_create("results")

################################################################################
# Illustrative figures mass action vs configuration network
# (standard degree distributions)
# 1. homogeneous mass action
# 2. heterogeneous degree distribution

### parameter values
# mean degree
mu <- 2.4
infectiousness_rate <- 1 / 5
recovery_rate <- 1 / 5
reproduction_number <- 1.6

# degree parameter values
n <- 5
lambda <- mu
p <- 1 / mu
size <- n
prob <- size / (mu + size)

# variance of degree
var_p <- var_degree(list(degree_distribution = "poisson", lambda = lambda))
var_nb <- var_degree(list(degree_distribution = "negative_binomial", size = size, prob = prob))

# determine transmission rate as function of reproduction number and parameters
beta_ma <- reproduction_number * recovery_rate
beta_p <- reproduction_number / ((var_p + mu^2 - mu) / mu - reproduction_number) * recovery_rate
beta_nb <- reproduction_number / ((var_nb + mu^2 - mu) / mu - reproduction_number) * recovery_rate

# gather models in tibble
df <- tribble(
  ~name,                   ~func,             ~transmission_rate, ~degree_distribution, ~degree_distribution_display, ~degree_params,
  "Mass action",           "model_reference", beta_ma,            NA,                   "-",                          NA,
  "Configuration network", "model_network",   beta_p,             "poisson",            "Poisson",                    list(lambda = lambda),
  "Configuration network", "model_network",   beta_nb,            "negative_binomial",  "Negative Binomial",          list(size = size, prob = prob)
)

# helper function to apply function to each row in tibble
epidemic_model <- function(func, transmission_rate, degree_distribution, degree_params) {
  if (func == "model_reference") {
    data <- model_reference(
      transmission_rate = transmission_rate,
      infectiousness_rate = infectiousness_rate,
      recovery_rate = recovery_rate,
      population_size = 1,
      time_end = 1000
    )
  } else {
    data <- do.call(
      func,
      c(
        list(
          degree_distribution = degree_distribution,
          infection = "SEIR",
          transmission_rate = transmission_rate,
          infectiousness_rate = infectiousness_rate,
          recovery_rate = recovery_rate,
          population_size = 1,
          time_end = 1000
        ),
        degree_params
      )
    )
  }
  return(list(data))
}


# simulate each model and store simulated data in tibble
df <- df |>
  rowwise() |>
  mutate(
    data = epidemic_model(func, transmission_rate, degree_distribution, degree_params),
    final_size_sim = final_size(data)
  )

df_manuscript <- df |>
  bind_cols(
    degree_mean = c("-", mu, mu),
    degree_variance = c("-", var_p, var_nb),
    reproduction_number = reproduction_number
  ) |>
  select(
    "Contact structure" = name,
    "Degree distribution" = degree_distribution_display,
    "Degree mean" = degree_mean,
    "Degree variance" = degree_variance,
    "Transmission rate" = transmission_rate,
    "Reproduction number" = reproduction_number,
    "Final size" = final_size_sim
  )

df_manuscript |> readr::write_csv("results/final_size_ma_vs_network.csv")

# manipulate tibble for plotting
df_plt <- purrr::map(1:3, \(x) {
  p <- df$data[[x]] |> plot_epidemic()
  name <- df$name[x]
  degree_distribution <- df$degree_distribution_display[x]
  df_plt <- p$data |> mutate(
    name = name,
    degree_distribution = degree_distribution
  )
}) |> bind_rows()

# plot for report
# course of epidemic through recovered individuals
plt <- df_plt |>
  filter(
    compartment == "R",
    name == "Configuration network",
    time <= 300
  ) |>
  ggplot(aes(
    x = time, y = value,
    col = degree_distribution
  )) +
  geom_line(linewidth = 1) +
  # homogeneous mass action model
  geom_line(
    data = df_plt |>
      filter(
        compartment == "R", name == "Mass action",
        time <= 300
      ),
    linetype = "dashed",
    linewidth = 1
  ) +
  ylim(c(0, 1)) +
  labs(
    x = "time",
    y = "fraction infected",
    col = "Reference model",
  ) +
  scale_color_discrete(
    limits = c("-", "Poisson", "Negative Binomial"),
    labels = c(
      "-" = "mass action",
      "Poisson" = "configuration network\n(Poisson)",
      "Negative Binomial" = "configuration network\n(negative binomial)"
    )
  ) +
  theme_light()
plt

ggsave("results/configuration_network.png",
       plot = plt, width = 14, height = 8, units = "cm"
)
