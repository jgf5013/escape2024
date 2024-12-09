################################################################################
#
# Figure 5
#
################################################################################


library(dplyr)
library(ggplot2)
library(patchwork)

fs::dir_create("results")

### generate samples from negative binomial distributions -----
# Define the mean as 5 and variance ranging from 5 to 10
mu <- 2.4
variances <- 3.6

# Create a data frame and calculate parameters using dplyr's mutate
params_df <- data.frame(mean = mu, var = variances) %>%
  mutate(
    prob = mean / var,
    size = mean^2 / (var - mean)
  )

# Check the resulting parameters
print(params_df)

# Generate samples from the negative binomial distribution
n_samples <- 10^5
samples_list <- params_df %>%
  rowwise() %>%
  mutate(samples = list(rnbinom(n_samples, size = size, prob = prob))) %>%
  ungroup()

### compute p th moement for each sample set -----
compute_pth_moment <- function(p, samples) {
  # define N
  N <- length(samples)
  # remove 0 from samples (= this means that we define 0^a = 0 for any a, which is necessary to calculate the Lehmer mean for p =< 0)
  samples_nonzero <- samples[samples > 0]
  # compute p-th moment
  moment <- sum(samples_nonzero^p) / N
  # return
  return(moment)
}

### Compute Lp_value for each sample set -----
# Define p_range
p_range <- seq(-10, 50, by = 0.01)

# Lp_value
# Extract the samples
samples <- samples_list$samples[[1]]
# Compute Lp_value for the current sample set
Lp_value <- unlist(lapply(p_range, function(x) {
  compute_pth_moment(p = x, samples = samples) / compute_pth_moment(p = (x - 1), samples = samples)
}))

results_df <- data.frame(p = p_range, Lp_value = Lp_value, variance = samples_list$var[1])

### Plot using ggplot2 -----
Lp_line_plot <- ggplot(
  results_df,
  aes(x = p, y = Lp_value)
) +
  geom_line() +
  xlim(c(min(p_range), max(p_range))) +
  labs(
    title = "L_p Value (NB distribution with mean of 2.4 and variance of 3.6)",
    x = "p",
    y = "L_p Value",
    color = "Variance"
  ) +
  theme_bw() +
  geom_line(size = 1)

Lp_line_zoomed_plot <- ggplot(
  results_df %>% filter(p >= 1) %>% filter(p <= 2),
  aes(x = p, y = Lp_value)
) +
  geom_line() +
  xlim(c(1, 2)) +
  labs(
    title = "",
    x = "p",
    y = "L_p Value",
    color = "Variance"
  ) +
  theme_bw()

Lp_line_plot_merged <- Lp_line_plot +
  inset_element(Lp_line_zoomed_plot, left = 0.55, bottom = 0.01, right = 0.99, top = 0.60)

# Create a data frame for ggplot
sampled_degrees_infectee <- samples # or, rnbinom(n_samples, size = 5, prob = 0.5) # mean 5  var 10 prob 0.5000000  size 5.000000
df_infectee <- data.frame(degree = sampled_degrees_infectee)

# Define f(x) as (x/<x>)*p(x)
f_x <- function(x, size, prob) {
  # Mean <x> (expected value)
  mean_x <- (size * (1 - prob)) / prob

  # Negative binomial probability p(x)
  p_x <- dnbinom(x, size = size, prob = prob)

  # Compute f(x) = (x/<x>) * p(x)
  f_value <- (x / mean_x) * p_x

  # Return the value of f(x)
  return(f_value)
}
x_values <- 0:10^4
f_x_dens <- f_x(x_values, size = params_df$size, prob = params_df$prob)
sampled_degrees_infector <- sample(x_values, size = n_samples, prob = f_x_dens, replace = T)
df_infector <- data.frame(degree = sampled_degrees_infector)

# these values should be close
mean(sampled_degrees_infectee^2) / mean(sampled_degrees_infectee)
mean(sampled_degrees_infector)

# Plot histogram using ggplot2
degree_hist_infectee_plot <- ggplot(df_infectee, aes(x = degree)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", boundary = 0.5) +
  scale_x_continuous(breaks = seq(0, 20, by = 1)) +
  theme_bw() +
  xlim(c(-1, 21)) +
  labs(
    title = "Random individuals' contact degree Distribution",
    x = "Degree (x)",
    y = "Frequency"
  )

degree_hist_infector_plot <- ggplot(df_infector, aes(x = degree)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", boundary = 0.5) +
  scale_x_continuous(breaks = seq(0, 20, by = 1)) +
  theme_bw() +
  xlim(c(-1, 21)) +
  labs(
    title = "Identified cases' contact degree Distribution",
    x = "Degree (x)",
    y = "Frequency"
  )

# combine plots into one figure
p <- ((degree_hist_infectee_plot / degree_hist_infector_plot) | (Lp_line_plot_merged)) + plot_annotation(tag_levels = "A") # 1200 x 577


ggsave("results/Lp_massaction_vs_network.png",
       plot = p, width = 1200, height = 577, units = "px", dpi = 100
)
