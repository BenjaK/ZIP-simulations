

# Week of (first) detection ----------------------------------------------------
calc_wod <- function(df, alpha) {
  df %>%
    group_by(Method, mu, p, true_zone_nr, q, simulation) %>%
    summarise(wod = min(which(pvalue <= alpha)),
              detected_zone_nr = detected_zone_nr[wod],
              detected_duration = detected_duration[wod],
              precision = precision[wod],
              recall = recall[wod]) %>%
    ungroup %>%
    mutate(alpha = alpha)
}

calc_wod_multiple <- function(df, alphas) {
  out <- NULL
  for (i in seq_along(alphas)) {
    a <- alphas[i]
    print(paste0("i = ", i, " / ", length(alphas)))
    out <- rbind(out, calc_wod(df, a))
  }
  out %>% arrange(mu, p, true_zone_nr, q)
}


# Proportion detected within time window ---------------------------------------

calculate_FWER <- function(wod_df) {
  wod_df %>%
    filter(q == 1) %>%
    group_by(Method, mu, p, alpha) %>%
    summarize(n_negatives = n(),
              false_positives = sum(wod < Inf)) %>%
    ungroup %>%
    mutate(
      FWER = false_positives / n_negatives,
      specificity = 1 - FWER)
}

count_true_positives <- function(wod_df) {
  wod_df %>%
    filter(q > 1) %>%
    group_by(Method, mu, p, q, true_zone_nr, alpha) %>%
    summarize(n_positives = n(),
              true_positives = sum(wod < Inf)) %>%
    ungroup %>%
    mutate(
      TPR = true_positives / n_positives)
}


roc <- function(wod_df) {
  FWER <- wod_df %>%
    filter(q == 1) %>%
    group_by(mu, p, alpha) %>%
    summarize(negatives = n(),
              false_positives = sum(wod < Inf)) %>%
    ungroup %>%
    mutate(FWER = false_positives / negatives)
  tpr <- wod_df %>%
    filter(q > 1) %>%
    group_by(mu, p, q, true_zone_nr, alpha) %>%
    summarize(positives = n(),
              true_positives = sum(wod < Inf)) %>%
    ungroup %>%
    mutate(TPR = true_positives / positives)
  right_join(FWER, tpr, by = c("mu", "p", "alpha")) %>%
    dplyr::select(mu, p, true_zone_nr, q, alpha,
                  negatives, positives, false_positives, true_positives,
                  FWER, TPR)
}

detected_in_time <- function(wod_df) {
  wod_df %>%
    filter(q > 1) %>%
    group_by(mu, p, q, true_zone_nr, alpha) %>%
    summarize(n_detected = sum(wod < Inf),
              prop_detected = n_detected / n(),
              prop_undetected = 1 - prop_detected) %>%
    ungroup
}

detected_overall <- function(df, alphas) {
  out <- NULL
  for (a in alphas) {
    out <- rbind(out, df %>%
      filter(q > 1) %>%
      group_by(mu, p, q, true_zone_nr, true_duration) %>%
      summarize(n_detected = sum(pvalue < a),
                prop_detected = mean(pvalue < a)) %>%
      mutate(prop_undetected = 1 - prop_detected,
             alpha = a) %>%
      ungroup)
  }
  out
}

bivariate_power <- function(df, zones, alpha, max_zone_length = 25) {
  n <- nrow(df)
  df2 <- df %>% filter(pvalue <= alpha)
  df2 %<>% mutate(detected_zone_length = purrr::map_int(detected_zone_nr, ~ length(zones[[.]])))

  true_zone_nr <- df$true_zone_nr[1]
  true_zone <- zones[[true_zone_nr]]
  max_width <- true_zone_nr # Equal to outbreak zone length
  out <- matrix(0L, max_zone_length, max_width)
  for (i in seq_len(max_zone_length)) {
    df3 <- df2 %>% filter(detected_zone_length == i)
    detected_zones <- df3 %>% pull(detected_zone_nr)

    isects <- purrr:::map_int(detected_zones, ~ length(base::intersect(zones[[.]], true_zone)))
    isects <- isects[isects > 0]

    for (j in isects) {
      out[i, j] <- out[i, j] + 1
    }
  }

  margin_table <- out %>%
    rbind(., colSums(.)) %>%
    cbind(., rowSums(.)) %>%
    apply(., 2, as.integer) %>%
    as.data.frame
  rownames(margin_table) <- c(seq_len(max_zone_length), "Total")
  names(margin_table) <- c(seq_len(max_width), "Total")

  return(list(table = out,
              margin_table = margin_table,
              n = n))
}
