
# Convert week numbers to number of weeks since outbreak start (starting at 1)
# Also convert some variables to integers
convert_weeks_etc <- function(df) {
  df %>% mutate(true_duration = as.integer(Week - min(Week) + 1),
                detected_duration = as.integer(duration),
                duration = NULL,
                Week = NULL,
                true_zone_nr = as.integer(true_zone_nr),
                detected_zone_nr = as.integer(detected_zone_nr))
}


add_simulation_nr <- function(df) {
  df %>%
    group_by(mu, p, true_zone_nr, q, true_duration) %>%
    mutate(simulation = 1:N_outbreak_sims) %>%
    ungroup
}
