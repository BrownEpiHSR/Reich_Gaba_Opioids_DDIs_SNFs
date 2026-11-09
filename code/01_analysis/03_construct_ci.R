################################################################################################################################
# Purpose: Construct 95% credible intervals using output from 02_bootstrap_estimates.R and create table of results
#          at 30 and 99 days post-SNF admission.
#
# Input datasets:
#
#   The primary input dataset is the final dataset created from the 02_bootstrap_estimates.R file (outcome_name_bootstrap_no_covs_FINAL_rundate_nsuccess.rds OR
#                                                                                                  outcome_name_bootstrap_with_covs_FINAL_rundate_nsuccess.rds).
#     - This dataset has a row for each t_intrv from a given bootstrap replicate.
#       - So the total rows should be # of successful replicates x (the total discrete-time intervals).
#
#   The dataset from 01_outcome_estimation.R is also used to overlay the credible intervals with the point estimates.  
#
################################################################################################################################

## Import libraries
library(data.table)
library(haven)
library(tidyverse)
library(dplyr)
library(modelsummary)
library(fs)
library(here)

## Function to construct 95% credible intervals
construct_ci <- function(bootstrap_in_path, wt_plr_in_path, 
                         outcome_name, include_covs_out_mod, stabilized_wts, model_counter, save_path) {
  
  ## Create directory
  cov_folder <- if (include_covs_out_mod == 1) {
    "with_covs/"
  } else {
    "no_covs/"
  }
  
  wt_folder <- if(stabilized_wts == 1) {
    "stab_wts/"
  } else {
    "nonstab_wts/"
  }
  
  dir <- file.path(save_path, paste0(outcome_name, "/03_construct_ci/model_", model_counter, "/", cov_folder, wt_folder))
  fs::dir_create(dir)
  
  ## Read in data
  # Bootstrap data
  bootstrap_d <- readRDS(bootstrap_in_path)
  
  # Point estimate data
  point_d <- readRDS(wt_plr_in_path)
  
  ## Calculate 95% confidence (credible) intervals
  # Note that bootstrap_d is in a long format (i.e., unique row for t_intrv for each bootstrap iteration), so we need to group by t_intrv
  bootstrap_ci <- bootstrap_d %>%
    group_by(t_intrv) %>%
    summarise(
      pr_ev_lc_1 = quantile(pr_ev_1, 0.025, na.rm=T), # Assign = 1 cumulative incidence lower limit
      pr_ev_lc_0 = quantile(pr_ev_0, 0.025, na.rm=T), # Assign = 0 cumulative incidence lower limit
      cid_lc     = quantile(cid, 0.025, na.rm=T),     # Risk difference lower limit
      cir_lc     = quantile(cir, 0.025, na.rm=T),     # Risk ratio lower limit
      pr_ev_uc_1 = quantile(pr_ev_1, 0.975, na.rm=T), # Assign = 1 cumulative incidence upper limit
      pr_ev_uc_0 = quantile(pr_ev_0, 0.975, na.rm=T), # Assign = 0 cumulative incidence upper limit
      cid_uc     = quantile(cid, 0.975, na.rm=T),     # Risk difference upper limit
      cir_uc     = quantile(cir, 0.975, na.rm=T)      # Risk ratio upper limit
    )
  
  ## Merge 95% confidence (credible) intervals onto point estimates
  final_estimates_day <- point_d %>%
    left_join(bootstrap_ci, by = "t_intrv")
  
  ### Plot weighted point estimates w/ CLs
  
  ## Plot risk difference with 95% CLs

  gg_plr_cid <- ggplot(final_estimates_day, aes(x = time, y = cid)) +
    # Add a reference line at 0 (no difference)
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    # Plot the risk difference curve
    geom_line(aes(y = cid), color = "#2EC4B6") + 
    # Plot the lower and upper confidence limits
    geom_ribbon(aes(ymin = cid_lc, ymax = cid_uc), alpha = 0.15, fill = "#2EC4B6") +
    # Axis scale
    scale_x_continuous(limits = c(0, max(final_estimates_day$time, na.rm = TRUE)), 
                       breaks = seq(0, max(final_estimates_day$time, na.rm = TRUE), 7)) + # Controls interval length between x-axis ticks
    theme_bw() +
    # Labels
    labs(
      title = "Cumulative Risk Difference",
      subtitle = paste0("Outcome: ", outcome_name,
                        "\nCovariates in outcome model? ",
                        ifelse(include_covs_out_mod == 1, "Yes", "No"),
                        "\nWeights used: ",
                        ifelse(stabilized_wts == 1, "Stabilized", "Nonstabilized")),
      x = "Days Since SNF Admission",
      y = "Risk Difference (Assign 1 - Assign 0)"
    ) 

  ## Plot risk ratios with 95% CLs

  gg_plr_cir <- ggplot(final_estimates_day, aes(x = time, y = cir)) +
    # Add a reference line at 1 (no difference)
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    # Plot the risk ratio curve
    geom_line(aes(y = cir), color = "#2EC4B6") +
    # Plot the lower and upper confidence limits
    geom_ribbon(aes(ymin = cir_lc, ymax = cir_uc), alpha= 0.15, fill = "#2EC4B6") +
    # Axis scale
    scale_x_continuous(limits = c(0, max(final_estimates_day$time, na.rm = TRUE)), 
                       breaks = seq(0, max(final_estimates_day$time, na.rm = TRUE), 7)) + # Controls interval length between x-axis ticks
    theme_bw() +
    # Labels
    labs(
      title = "Cumulative Risk Ratio",
      subtitle = paste0("Outcome: ", outcome_name,
                        "\nCovariates in outcome model? ",
                        ifelse(include_covs_out_mod == 1, "Yes", "No"),
                        "\nWeights used: ",
                        ifelse(stabilized_wts == 1, "Stabilized", "Nonstabilized")),
      x = "Days Since SNF Admission",
      y = "Risk Ratio (Assign 1 / Assign 0)"
    )

  ## Plot cumulative incidence curves with 95% CLs
  
  gg_cum_inc <- ggplot(final_estimates_day,
                       aes(x = time)) +
    geom_line(aes(y = pr_ev_1, color = "Assign = 1")) +
    geom_line(aes(y = pr_ev_0, color = "Assign = 0")) +
    geom_ribbon(aes(ymin = pr_ev_lc_1, ymax = pr_ev_uc_1, fill = "Assign = 1"), color = NA, alpha = 0.3) +
    geom_ribbon(aes(ymin = pr_ev_lc_0, ymax = pr_ev_uc_0, fill = "Assign = 0"), color = NA, alpha = 0.3) +
    # Axis scale
    scale_x_continuous(limits = c(0, max(final_estimates_day$time, na.rm = TRUE)), 
                       breaks = seq(0, max(final_estimates_day$time, na.rm = TRUE), 7)) + # Controls interval length between x-axis ticks
    theme_bw() +
    scale_color_manual(values = c("Assign = 1" = "#FF6B8B", "Assign = 0" = "#2EC4B6")) +
    scale_fill_manual(values = c("Assign = 1" = "#FF6B8B", "Assign = 0" = "#2EC4B6")) +
    labs(
      title = "Cumulative incidence",
      subtitle = paste0("Outcome: ", outcome_name,
                        "\nCovariates in outcome model? ",
                        ifelse(include_covs_out_mod == 1, "Yes", "No"),
                        "\nWeights used: ",
                        ifelse(stabilized_wts == 1, "Stabilized", "Nonstabilized")),
      x = "Days since baseline",
      y = "Cumulative incidence",
    )  +
    guides(
      fill = "none"
    )
  
  ### Save dataset and figures
  
  # Save dataset with 95% CLs
  cov_spec <- if (include_covs_out_mod == 1) {
    "_w_covs_"
  } else {
    "_no_covs_"
  }
  
  saveRDS(final_estimates_day, file = file.path(dir, paste0("/final_estimates", cov_spec, outcome_name, ".rds")))
  
  ## Save figures
  
  # RD
  ggsave(filename = file.path(dir, paste0("/rd", cov_spec, outcome_name, ".png")),
         plot = gg_plr_cid, 
         dpi = 300, 
         width = 6, 
         height = 6, 
         units = "in")
  # RR
  ggsave(filename = file.path(dir, paste0("/rr", cov_spec, outcome_name, ".png")),
         plot = gg_plr_cir, 
         dpi = 300, 
         width = 6, 
         height = 6, 
         units = "in")
  
  # Cumulative incidence
  ggsave(filename = file.path(dir, paste0("/cum_inc", cov_spec, outcome_name, ".png")),
         plot = gg_cum_inc, 
         dpi = 300, 
         width = 6, 
         height = 6, 
         units = "in")
  
  ## Return final dataset with 95% CLs for in-session inspection
  return(final_estimates_day)
  
}

## Call function

## Example function call 
# All-cause mortality (nonstabilized weights)
death_nonstab <- construct_ci(
  bootstrap_in_path = "//your/path/death_bootstrap_with_covs_FINAL_20260910_n500.rds", # Bootstrapped data from 02_bootstrap_estimates.R
  wt_plr_in_path = "//your/path/wt_plr_covs_wide_death.rds", # Output dataset from 01_outcome_estimation.R with point estimates; If your bootstrapping included covariates in the outcome model, make sure this dataset did too
  outcome_name = "death",   # used to name folder, datasets, and figures
  include_covs_out_mod = 1, # 1 if you included covariates in the outcome model for the bootstrapping; 0 otherwise (this is used to name the final dataset)
  stabilized_wts = 0,       # 1 if you used stabilized weights in the bootstrapping; 0 if you used nonstabilized weights in the bootstrapping
  model_counter = 1,        # used to name folder
  save_path = "//your/path"
)


## Function to create table with estimates from specific time points of interest
make_results_table <- function(model_dict, outcome_name, t_targets, save_path) {
  
  ## Create directory
  dir <- file.path(save_path, paste0(outcome_name, "/03_construct_ci/"))
  fs::dir_create(dir)
  
  ## Code to process each file and restrict to time points of interest
  process_file <- function(path, outcome_name, model_type) {
    readRDS(path) %>%
      filter(time %in% t_targets) %>%
      mutate(
        Outcome = outcome_name,
        Model = model_type
      ) %>%
      select(Outcome, Model, everything())
  }
  
  # Map through dictionary
  df <- pmap_dfr(model_dict, function(path, model_type) {
    process_file(
      path = path,
      outcome_name = outcome_name,
      model_type = model_type
    )
  })

  # Helper function to format "Estimate (LC, UC)" or just "Estimate" if CIs are missing
  fmt_ci <- function(est, lc, uc, digits = 2) {
    ifelse(
      is.na(lc) | is.na(uc),
      as.character(round(est, digits)),
      paste0(round(est, digits), " (", round(lc, digits), ", ", round(uc, digits), ")")
    )
  }
  
  # Combine and arrange
  final_table <- df %>%
    arrange(Outcome, Model, time) %>%
    mutate(
      pr_ev_1_perc = pr_ev_1 * 100,
      pr_ev_lc_1_perc=pr_ev_lc_1 * 100,
      pr_ev_uc_1_perc=pr_ev_uc_1 * 100,
      pr_ev_0_perc = pr_ev_0 * 100,
      pr_ev_lc_0_perc=pr_ev_lc_0 * 100,
      pr_ev_uc_0_perc=pr_ev_uc_0 * 100,
      cid_perc = cid * 100,
      cid_lc_perc = cid_lc * 100,
      cid_uc_perc = cid_uc * 100,
      
      # Formatted columns
      CumInc_Perc_1 = fmt_ci(pr_ev_1_perc, pr_ev_lc_1_perc, pr_ev_uc_1_perc),
      CumInc_Perc_0 = fmt_ci(pr_ev_0_perc, pr_ev_lc_0_perc, pr_ev_uc_0_perc),
      RD_Perc = fmt_ci(cid_perc, cid_lc_perc, cid_uc_perc),
      RR = fmt_ci(cir, cir_lc, cir_uc)
    ) %>%
    select(Outcome, Model, time, CumInc_Perc_1, CumInc_Perc_0, RD_Perc, RR)

  ## Save dataset as CSV file
  write.csv(final_table, file.path(dir, "final_table_est.csv"))
  
  ## Return table in-session
  return(final_table)
    
}

## Call function
# Create dictionary with all the models you want to include in the table
  # You can just list one model here or multiple
model_dict <- tibble(
  model_type = c("PLR-Weights-covs"), # Give a name for each model. This will appear in the table output.
  path = c("//your/path/final_estimates_w_covs_death.rds")
)

# Example call
# Note: SNF admission is day 0. time = 29 corresponds to t_intrv = 11 and the first 30 follow-up days; time = 98 corresponds to t_intrv = 34 and the first 99 follow-up days.
example_make_results_table <- make_results_table(
  model_dict = model_dict,
  outcome_name = "death", # Name for the outcome. This will appear in the table output
  t_targets = c(29, 98),  # Vector of time points you want to include in the results table. You can include one or more values in this vector, which should be based on the time variable (i.e., day-level)
  save_path = "//your/path/Output/"
)

