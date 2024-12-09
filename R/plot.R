#' Title
#'
#' @param .data
#'
#' @return
#' @import ggplot2
#' @export
#'
#' @examples
plot_epidemic <- function(.data) {
  infection <- attr(.data, "infection")
  if (infection == "SIR") {
    cols <- c("S", "I", "R")
  } else if (infection == "SEIR") {
    cols <- c("S", "E", "I", "R")
  }
  p <- .data |>
    tidyr::pivot_longer(cols = cols, names_to = "compartment") |>
    ggplot(aes(
      x = time,
      y = value,
      col = factor(compartment,
        levels = c("S", "E", "I", "R"),
        labels = c("Susceptible", "Exposed", "Infectious", "Recovered")
      )
    )) +
    geom_line() +
    scale_color_discrete(name = "Compartment") +
    ylab("population size") +
    theme_light()
  return(p)
}
