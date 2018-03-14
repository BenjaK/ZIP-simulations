library(magrittr)
library(tidyverse)
library(scanstatistics)

source("./code/performance_metric_functions.R")

# Define simulation parameters -------------------------------------------------

# Generate locations and zones
set.seed(20180101)
n_locations <- 100
geo <- data.frame(location = 1:n_locations,
                  x = runif(n_locations, -1, 1),
                  y = runif(n_locations, -1, 1))

# Time parameters
max_duration <- 10
weeks_scanned <- 20
max_zone_length <- 25

# Create zones
zones <- coords_to_knn(geo[, c("x", "y")], k = max_zone_length) %>% knn_zones

# Baseline parameters
scenarios <- expand.grid(mu = c(1, 5, 10),
                         p = c(0, 0.01, 0.05, 0.15, 0.25, 0.5),
                         q = c(1, 1.1, 1.25, 1.5),
                         obz = c(1, 5, 10, 15))

# Significance levels
alphas <- c(0.001, 0.005, seq(0.01, 0.1, by = 0.005))

# Change values here to explore data
display_scenario <- list(
  mu = 5,
  p = 0.15,
  true_zone_nr = 5,
  true_duration = 4,
  q = 1.5,
  alpha = 0.05)


# Load simulation data ---------------------------------------------------------
sim_df <- readRDS("./data/sim_df.rds")

# Calculate week of first detection for each significance level ----------------

# False positives
fp_wod_df <- calc_wod_multiple(sim_df %>% filter(q == 1), alphas)

# True positives
tp_wod_df <- calc_wod_multiple(sim_df %>% filter(q > 1), alphas)


# Family-wise error rate (FWER) ------------------------------------------------
fwer_df <- fp_wod_df %>% calculate_FWER

(multi_fpr_plot <- fwer_df %>%
    filter(p > 0, alpha <= 0.1) %>%
    mutate(Method = as.character(Method)) %>%
    bind_rows(fwer_df %>% 
                filter(Method == "UC-ZIP") %>%
                filter(p > 0, alpha <= 0.1) %>%
                mutate(Method = "Nominal") %>%
                mutate(FWER = 1 - (1 - alpha)^11)) %>%
    mutate(Method = factor(Method, 
                           levels = c("UC-ZIP", "UC-POI", 
                                      "KC-POI", "Nominal"))) %>%
    ggplot +
    geom_line(aes(x = alpha, y = FWER, 
                  linetype = Method, 
                  size = Method,
                  color = Method)) +
    # geom_line(aes(x = alpha, y = FWER),
              # linetype = "longdash", size = 1.3, alpha = 0.5) +
    # scale_linetype_manual(values = rep("solid", 3)) +
    scale_size_manual(values = c("UC-ZIP" = 0.5, "UC-POI" = 0.5, 
                                 "KC-POI" = 0.5, "Nominal" = 1)) +
    scale_color_manual(values = c("UC-ZIP" = "black", "UC-POI" = "black", 
                                 "KC-POI" = "black", "Nominal" = "grey50")) +
    # xlab(expression(alpha)) +
    xlab(bquote("Significance level" ~ alpha)) +
    ylab("False positive rate") +
    ylim(0, 1) +
    theme_bw() +
    facet_grid(p ~ mu,
               labeller = label_bquote(cols = mu ~ "=" ~.(mu),
                                       rows = p ~ "=" ~ .(p)))
)


# Week of first detection and F-score ------------------------------------------

(joint_wod_Fscore_plot <- tp_wod_df %>%
    filter(mu == display_scenario$mu,
           p == display_scenario$p,
           q == display_scenario$q,
           true_zone_nr == display_scenario$true_zone_nr,
           alpha == display_scenario$alpha) %>%
    mutate(Fscore = 2 / (1 / precision + 1 / recall)) %>%
    # mutate(F = ifelse(is.na(F), 100, F)) %>%
    mutate(segment = cut(Fscore, c(-Inf, 0, 0.25, 0.5, 0.75, Inf))) %>%
    mutate(segment = factor(segment,
                            exclude = NULL,
                            levels = rev(levels(segment)),
                            labels = rev(c("0", "(0, 0.25]", "(0.25, 0.5]",
                                           "(0.5, 0.75]", "(0.75, 1]")))) %>%
    mutate(Week = factor(wod, levels = c(1:11, Inf))) %>% #exclude = NULL)) %>%
    group_by(Week, Method, segment) %>%
    tally %>%
    ungroup %>%
    ggplot +
    geom_col(aes(x = Week, y = n, fill = segment),
             color = "black", size = 0.1) +
    scale_fill_brewer("F-score\nsegment", palette = "Greys") +
    xlab("Week of detection") + ylab("Outbreaks detected") +
    scale_x_discrete(breaks = c(1:11, Inf), labels = c(1:11, "Not\ndetected")) +
    theme_bw() +
    facet_grid(Method ~ ., scales = "free_y")
)

# Limit as p -> 0 --------------------------------------------------------------

(limp0_plot <- sim_df %>%
   filter(Method %in% c("UC-ZIP", "UC-POI")) %>%
   filter(mu == display_scenario$mu,
          q == display_scenario$q,
          true_zone_nr == 10, #display_scenario$true_zone_nr,
          true_duration == display_scenario$true_duration) %>%
   mutate(Fscore = 2 / (1 / precision + 1 / recall)) %>%
   group_by(Method, p) %>%
   summarize(low = quantile(Fscore, 0.05),
             med = median(Fscore),
             high = quantile(Fscore, 0.95)) %>%
   ungroup %>%
   mutate(Method = factor(Method, levels = c("UC-ZIP", "UC-POI"))) %>%
   ggplot +
   geom_line(aes(x = p, y = med)) +
   geom_ribbon(aes(x = p, ymin = low, ymax = high), alpha = 0.2) +
   facet_grid(. ~ Method) +
   xlab(expression(Structural~zero~probability~italic(p))) +
   ylab(expression(italic(F)-score)) +
   theme_bw()
)

# Bivariate power distribution table -------------------------------------------

# Compute the bivariate power distribution for the three scan statistics
bivp_uczip <- sim_df %>%
  filter(Method == "UC-ZIP",
         mu == display_scenario$mu,
         p == display_scenario$p,
         q == display_scenario$q,
         true_zone_nr == display_scenario$true_zone_nr,
         true_duration == display_scenario$true_duration) %>%
  bivariate_power(zones, display_scenario$alpha, max_zone_length)

bivp_ucpoi <- sim_df %>%
  filter(Method == "UC-POI",
         mu == display_scenario$mu,
         p == display_scenario$p,
         q == display_scenario$q,
         true_zone_nr == display_scenario$true_zone_nr,
         true_duration == display_scenario$true_duration) %>%
  bivariate_power(zones, display_scenario$alpha, max_zone_length)

bivp_kcpoi <- sim_df %>%
  filter(Method == "KC-POI",
         mu == display_scenario$mu,
         p == display_scenario$p,
         q == display_scenario$q,
         true_zone_nr == display_scenario$true_zone_nr,
         true_duration == display_scenario$true_duration) %>%
  bivariate_power(zones, display_scenario$alpha, max_zone_length)

# Combine into a single table
bivp_table <- cbind(seq_len(max_zone_length),
                    bivp_uczip$table, rowSums(bivp_uczip$table),
                    bivp_ucpoi$table, rowSums(bivp_kcpoi$table),
                    bivp_ucpoi$table, rowSums(bivp_kcpoi$table))

# Print table in LaTeX format (still needs some manual editing to reproduce the
# table in the paper)
library(kableExtra)
bivp_table %>%
  knitr::kable(format = "latex",
               booktabs = TRUE,
               row.names = TRUE,
               caption = "Bivariate power distribution") %>%
  # kableExtra::kable_styling() %>%
  # column_spec(7, border_left = T, bold = T) %>%
  add_header_above(c("Length" = 1,
                     "Shared locations" = 5,
                     "Sum" = 1,
                     "Shared locations" = 5,
                     "Sum" = 1,
                     "Shared locations" = 5,
                     "Sum" = 1))
