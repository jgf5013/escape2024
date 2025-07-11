#' Title
#'
#' @param time_end
#' @param increment
#' @param params_p
#' @param params_nb
#' @param params_reference
#'
#' @returns
#' @export
#'
#' @examples
models_combined <- function(simulation_id, time_end, increment, params_p, params_nb, params_reference) {
  current_state_network_nb <- do.call(
    model_network_seir,
    c(list(simulation_id = simulation_id), params_nb)
  )
  current_state_network_p <- do.call(
    model_network_seir,
    c(list(simulation_id = simulation_id), params_p)
  )
  current_state_reference <- do.call(
    model_reference,
    c(list(simulation_id = simulation_id), params_reference)
  )
  for (t in seq(0, time_end, by = increment)) {
    simulate_outbreak_seir_network(t, increment, current_state_network_p, params_p)
    simulate_outbreak_seir_network(t, increment, current_state_network_nb, params_nb)
    simulate_outbreak_seir_reference(t, increment, current_state_reference, params_reference)
  }
}
