################################################################################################################################
# Purpose: Calculate artificial censoring weights for a given outcome-specific person-period dataset.
#
# Input dataset:
#
#    Cloned, person-period dataset (i.e., one row per discretized time interval for each clone).
#       This program assumes the dataset has the following variables:
#           id_var: the unique ID variable for each eligible patient.
#           assign: (0/1) indicates treatment assignment for the clone.
#           t_intrv: Discrete-time interval index (e.g., 0, 1, 2, ...).
#           event: 1 if the outcome occurred in the time interval, 0 otherwise.
#           artificial_cens: 1 if clone was artificially censored in the time interval, 0 otherwise.
#           natural_cens: 1 if clone was naturally censored in the time interval, 0 otherwise.
#           trt_init_prior: (0/1) For assign = 1, flags strategy initiation before the final grace-period day, from the initiation interval onward;
#                           intended for restriction at the final grace-period interval.
#           covariates: Baseline and/or time-varying covariates. Covariate variable names are specified in the covs argument. 
#		- Baseline covariates were added to the person-period datasets in a separate preparation step not shown here. 
#
#       Censoring logic:
#       - Terminal interval: Follow-up for a clone ends with the interval where the outcome or censoring occurs (i.e., event, natural_cens,
#         or artificial_cens = 1). No subsequent intervals should exist for that clone.
#       - Timing assumption: Consequently, clones are censored at the end of the discrete-time interval where the outcome or censoring occurs.
#
# Program notes:
#    - This program fits separate models to calculate artificial censoring probabilities for each treatment assignment.
#    - This program assumes model covariates are complete with no missing values. In cases where covariates are missing for some clones, like max pain rating, create a separate
#      category for "not measured" or "missing." 
#    - Because clones are censored at the end of a given discrete-time interval, the weights calculated in one interval are applied to the subsequent interval.
#    - This program calculates both nonstabilized and stabilized weights. art_cens_ipw is the variable for the nonstabilized weights; art_cens_stab_ipw represents the stabilized weights.
#    - This program aligns with the artificial weight estimation procedure in the primary analysis. That is, in the gabapentinoids + opioids strategy,
#      we include all clones who are alive and uncensored in the last interval of the grace period when modeling the probability of treatment deviation.
#		- For our last-day initiator IPCW specification, we simply restricted the assign = 1 model to those who have yet to initiate the strategy in the last
#		  interval of the grace period (using trt_init_prior)
#			- Then, we manually gave clones who initiated the strategy prior to the last interval of the grace period a weight of 1.
#   
################################################################################################################################

## Import libraries
library(haven)
library(data.table)
library(openxlsx)
library(tidyverse)
library(broom)
library(dplyr)
library(modelsummary)
library(data.table)
library(rlang)
library(fs)
library(here)

## Function to calculate artificial censoring weights
fit_artificial_cens_weights <- function(data_in_path, outcome_name, covs, id_var,
                                        time_func_assign0=NULL, time_func_assign1=NULL, grace_intrv_end, model_counter, save_path) {
  
  ## Create directory
  dir <- file.path(save_path, paste0(outcome_name, "/00_art_cens_weights/model_", model_counter, "/"))
  fs::dir_create(dir)
  
  ## Read in the person-period dataset with covariates (sas or R dataset)
  # Extract the file extension and convert to lowercase
  ext <- tolower(tools::file_ext(data_in_path))
  
  # Read file based on extension
  if (ext == "rds") {
    data <- readRDS(data_in_path)
    
  } else if (ext == "sas7bdat") {
    data <- haven::read_sas(data_in_path)
    
  } else {
    stop("Unsupported file type. Only .rds and .sas7bdat files are accepted.")
  }
  
  ## Covariate prep
  
  # Collapse covariate names into one string
  cov_string <- paste(covs, collapse = " + ")
  
  ### Calculate artificial censoring probabilities for assign = 0 (opioids-only strategy)
  
  # Determine time component
  time_component0 <- if(!is.null(time_func_assign0)) {
    time_func_assign0
  } else {
    1
  }
  
  ## Prep numerator model formula object
  f_string_num0 <- paste("artificial_cens ~", time_component0)
  formula_num_obj0 <- as.formula(f_string_num0)
  
  ## Prep denominator model formula object
  # Combine time + covariates 
  time_part0 <- if (time_component0 == 1 && length(covs) > 0) NULL else time_component0
  cov_part0  <- if (length(covs) > 0) cov_string else NULL
  rhs_terms0 <- c(time_part0, cov_part0)
  
  # Construct formula object
  if (length(rhs_terms0) > 0) { # If there's a time function and/or covariates, add a + sign between each one
    rhs_string0 <- paste(rhs_terms0, collapse = " + ")
  } else {
    rhs_string0 <- "1" # Fallback to intercept-only if no time and no covariates
  }
  
  formula_string0 <- paste("artificial_cens ~", rhs_string0) # Add artificial_cens ~ at the beginning
  formula_den_obj0 <- as.formula(formula_string0)            # Turn the string into a formula object
  
  # Subset person-period dataset to assign = 0 and time period where artificial censoring occurs
  art_cens_0 <- data %>%
    filter(assign == 0,                   # Restrict to assign = 0
           if (time_component0 == 1) {    # If no time function was specified, restrict to the last interval in the grace period
             t_intrv == grace_intrv_end
           } else {
             t_intrv <= grace_intrv_end   # If a time function was specified, restrict to all intervals in the grace period
           }
          )
  
  ## Fit numerator model for assign = 0
  n_glm_wt0 = glm(formula_num_obj0, data = art_cens_0, family = binomial())
  
  ## Fit denominator model for assign = 0
  d_glm_wt0 = glm(formula_den_obj0, data = art_cens_0, family = binomial())
  
  # Retrieve predicted probabilities
  art_cens_0$pr_art_cens_n = n_glm_wt0$fitted.values # numerator; contains the interval-specific marginal probabilities of being artificially censored
  art_cens_0$pr_art_cens_d = d_glm_wt0$fitted.values # denominator
  
  ### Calculate artificial censoring probabilities for assign = 1 (gabapentinoids + opioids strategy)
  
  # Determine time component
  time_component1 <- if(!is.null(time_func_assign1)) {
    time_func_assign1
  } else {
    1
  }
  
  ## Prep numerator model formula object
  f_string_num1 <- paste("artificial_cens ~", time_component1)
  formula_num_obj1 <- as.formula(f_string_num1)
  
  # Combine time + covariates
  time_part1 <- if (time_component1 == 1 && length(covs) > 0) NULL else time_component1
  cov_part1  <- if (length(covs) > 0) cov_string else NULL
  rhs_terms1  <- c(time_part1, cov_part1) 
  
  # Construct formula object
  if (length(rhs_terms1) > 0) { # If there's a time function and/or covariates, add a + sign between each one
    rhs_string1 <- paste(rhs_terms1, collapse = " + ")
  } else {
    rhs_string1 <- "1" # Fallback to intercept-only if no time and no covariates
  }
  
  formula_string1 <- paste("artificial_cens ~", rhs_string1) # Add artificial_cens ~ at the beginning
  formula_den_obj1 <- as.formula(formula_string1)            # Turn the string into a formula object
  
  # Subset person-period dataset to assign = 1 and time period where artificial censoring occurs
  art_cens_1 <- data %>%
    filter(assign == 1,                   # Restrict to assign = 1
           if (time_component1 == 1) {    # If no time function was specified, restrict to the last interval in the grace period
             t_intrv == grace_intrv_end
           } else {
             t_intrv <= grace_intrv_end   # If a time function was specified, restrict to all intervals in the grace period
           }
    )
  
  # Fit numerator model
  n_glm_wt1 = glm(formula_num_obj1, data = art_cens_1, family = binomial())
  
  # Fit denominator model
  d_glm_wt1 = glm(formula_den_obj1, data = art_cens_1, family = binomial())
  
  # Retrieve predicted probabilities
  art_cens_1$pr_art_cens_n = n_glm_wt1$fitted.values # numerator; contains the interval-specific marginal probabilities of being artificially censored
  art_cens_1$pr_art_cens_d = d_glm_wt1$fitted.values # denominator
    
  ## Calculate artificial censoring weights
  
  # Append assign = 0 and assign = 1 artificial censoring datasets
  art_cens_all <- rbind(art_cens_0,
                        art_cens_1)
  
  # Merge art_cens_all with person-period dataset
  outcome_data_join = left_join(data,
                                 select(art_cens_all, {{ id_var }}, t_intrv, pr_art_cens_n, pr_art_cens_d, assign), 
                                 by=join_by({{ id_var }}, t_intrv, assign))
    # time intervals not included in the models above for a given treatment assignment will have pr_art_cens_n = NA and pr_art_cens_d = NA
  
  # Ensure order of variables before calculating artificial censoring weights in each interval
  setDT(outcome_data_join)
  id_str <- rlang::as_name(rlang::enquo(id_var))
  setorderv(outcome_data_join, cols = c(id_str, "assign", "t_intrv"))
  
  # Define interval condition logic for assign = 0
  cond_0 <- if (!is.null(time_func_assign0)) {
    outcome_data_join$t_intrv <= grace_intrv_end # If a time function was specified, calculate weights for every interval in the grace period
  } else {
    outcome_data_join$t_intrv == grace_intrv_end # If a time function was not specified, only calculate weights for the last interval in the grace period
  }
  
  # Define the interval logic for assign = 1
  cond_1 <- if (!is.null(time_func_assign1)) {
    outcome_data_join$t_intrv <= grace_intrv_end # If a time function was specified, calculate weights for every interval in the grace period
  } else {
    outcome_data_join$t_intrv == grace_intrv_end # If a time function was not specified, only calculate weights for the last interval in the grace period
  }
    
  # Initialize weights to 1
  outcome_data_join[, ":=" (art_cens_wt = 1, art_cens_stab_wt = 1)]
  
  # Calculate inverse probability of remaining uncensored for assign = 0
  outcome_data_join[assign == 0 & cond_0,
                    ":="(art_cens_wt = (1 / (1 - pr_art_cens_d)),
                         art_cens_stab_wt = ((1 - pr_art_cens_n) / (1 - pr_art_cens_d)))]
  
  # Calculate inverse probability of remaining uncensored for assign = 1
  outcome_data_join[assign == 1 & cond_1,
                    ":=" (art_cens_wt = (1 / (1 - pr_art_cens_d)),
                         art_cens_stab_wt = ((1 - pr_art_cens_n) / (1 - pr_art_cens_d)))]
    
  # Lag the weights to the next interval (i.e., apply the weights to the next interval)
  outcome_data_join[,  ":="(
    lagged_art_wt = shift(art_cens_wt, n = 1, fill = 1, type = "lag"),             # First row of each clone gets a weight of 1
    lagged_art_stab_wt = shift(art_cens_stab_wt, n = 1, fill = 1, type = "lag")),  # First row of each clone gets a weight of 1
    by = c(id_str, "assign")]
  
  # Calculate cumulative artificial censoring weights
  outcome_data_join[
    , `:=`(art_cens_ipw = cumprod(lagged_art_wt),
           art_cens_stab_ipw = cumprod(lagged_art_stab_wt))
    , by = c(id_str, "assign")
  ]
  
  ## Save person-period dataset with weights (art_cens_ipw and art_cens_stab_ipw)
  saveRDS(outcome_data_join,
          file = file.path(dir, paste0("pp_w_art_weights_", outcome_name, ".rds")))
  
  ###################################################################################
  # Output figures to review
  ###################################################################################
  
  ## Plot the cumulative percentage of clones removed over follow-up (relative to baseline)
  baseline <- outcome_data_join %>%
    filter(t_intrv == 0) %>%
    group_by(assign) %>%
    summarize(n_baseline = n_distinct({{ id_var }}))
  
  remaining <- outcome_data_join %>%
    group_by(assign, t_intrv) %>%
    summarize(n_remaining = n_distinct({{ id_var }}), .groups = "drop")
  
  cumulative_removal <- remaining %>%
    left_join(baseline, by = "assign") %>%
    mutate(
      pct_removed = 100 * (1 - n_remaining / n_baseline)
    )
  
   d_cum_perc_rem <- ggplot(cumulative_removal,
                            aes(x = t_intrv,
                                y = pct_removed,
                                color = factor(assign))) +
                            geom_line(size = 1.0) +
                            geom_point() +
                            scale_color_manual(
                              values = c("0" = "#2EC4B6", "1" = "#FF6B8B")
                            ) +
                            labs(
                              x = "Time Interval (t_intrv)",
                              y = "Cumulative % Removed",
                              color = "Assign"
                            ) +
                            theme_minimal()
  
   ggsave(
     filename = file.path(dir, paste0("cumulative_perc_removed.png")),
     plot = d_cum_perc_rem,
     dpi = 300, 
     width = 6,
     height = 6, 
     units = "in"
   )
  
  ## Save histogram with distribution of pr_art_cens_d by assign status
  dist_prob <- ggplot(art_cens_all, aes(x = pr_art_cens_d, fill = factor(assign))) +
    geom_histogram(
      aes(y = (after_stat(count / ave(count, group, FUN=sum)*100))),
      position = "identity",
      alpha = 0.4,
      bins = 40
    ) +
    scale_fill_manual(
      values = c(
        "0" = "#2EC4B6",
        "1" = "#FF6B8B" 
      ),
      name = "Assign"
    ) +
    labs(
      x = "Artificial Censoring Probability (denominator)",
      y = "Percent",
      title = "Artificial Censoring Probability (denom) Distribution by Assign Status"
    ) +
    theme_minimal() 
  
  ggsave(
    filename = file.path(dir, "prob_art_cens_distribution.png"),
    width = 8,
    height = 6,
    dpi = 300
  )
  
  ## Save boxplot with distribution of art_cens_ipw over follow-up for assign = 0
  png(
    filename = file.path(dir, "dist_wt_followup_assign0.png"),
    width  = 1400,   # in pixels
    height = 600,   # in pixels
    res    = 100    # resolution in ppi
  )
  
  # Configure 1 row, 2 columns
  par(mfrow = c(1, 2))
  
  # Panel 1: Nonstabilized
  boxplot_assign0 <- boxplot(
        art_cens_ipw ~ factor(t_intrv),
        data = outcome_data_join[outcome_data_join$assign == 0, ],
        xlab = "Time Interval (t_intrv)",
        ylab = "Cumulative ART Censoring IPW (Nonstabilized)",
        main = "ART Censoring IPW by Time Interval (Nonstabilized) \n(assign = 0)",
        notch = FALSE,
        varwidth = TRUE,
        outline = FALSE,
        col = "lightblue",
        border = "darkblue"
  )
  
  # Panel 2: Stabilized
  boxplot_stab_assign0 <- boxplot(
    art_cens_stab_ipw ~ factor(t_intrv),
    data = outcome_data_join[outcome_data_join$assign == 0, ],
    xlab = "Time Interval (t_intrv)",
    ylab = "Cumulative ART Censoring IPW (Stabilized)",
    main = "ART Censoring IPW by Time Interval (Stabilized) \n(assign = 0)",
    notch = FALSE,
    varwidth = TRUE,
    outline = FALSE,
    col = "lightblue",
    border = "darkblue"
  )
  
  dev.off()
  
  ## Save boxplot with distribution of art_cens_ipw over follow-up for assign = 1
  png(
    filename = file.path(dir, "dist_wt_followup_assign1.png"),
    width  = 1400,   # in pixels
    height = 600,   # in pixels
    res    = 100    # resolution in ppi
  )
  
  # Configure 1 row, 2 columns
  par(mfrow = c(1, 2))
  
  boxplot_assign1 <- boxplot(
    art_cens_ipw ~ factor(t_intrv),
    data = outcome_data_join[outcome_data_join$assign == 1, ],
    xlab = "Time Interval (t_intrv)",
    ylab = "Cumulative ART Censoring IPW (Nonstabilized)",
    main = "ART Censoring IPW by Time Interval (Nonstabilized)\n(assign = 1)",
    notch = FALSE,
    varwidth = TRUE,
    outline = FALSE,
    col = "lightblue",
    border = "darkblue"
  )
  
  boxplot_stab_assign1 <- boxplot(
    art_cens_stab_ipw ~ factor(t_intrv),
    data = outcome_data_join[outcome_data_join$assign == 1, ],
    xlab = "Time Interval (t_intrv)",
    ylab = "Cumulative ART Censoring IPW (Stabilized)",
    main = "ART Censoring IPW by Time Interval (Stabilized)\n(assign = 1)",
    notch = FALSE,
    varwidth = TRUE,
    outline = FALSE,
    col = "lightblue",
    border = "darkblue"
  )
  
  dev.off()
    
  ###################################################################################
  # Compile output for Excel file
  ###################################################################################
  
  ## Construct metadata and convergence summary table
  metadata_df <- data.frame(
    Parameter = c(
      "Model Run At",
      "Outcome",
      "ID Variable",
      "Grace Interval End",
      "Time Function (Assign = 0)",
      "Time Function (Assign = 1)",
      "Covariates Included",
      "Assign = 0 Numerator Model Converged?",
      "Assign = 0 Denominator Model Converged?",
      "Assign = 1 Numerator Model Converged?",
      "Assign = 1 Denominator Model Converged?"
    ),
    Value = c(
      format(Sys.time(), "%B %d, %Y at %I:%M %p"),
      outcome_name,
      id_str,
      as.character(grace_intrv_end),
      ifelse(is.null(time_func_assign0), "None", time_func_assign0),
      ifelse(is.null(time_func_assign1), "None", time_func_assign1),
      paste(covs, collapse = ", "),
      as.character(isTRUE(n_glm_wt0$converged)),
      as.character(isTRUE(d_glm_wt0$converged)),
      as.character(isTRUE(n_glm_wt1$converged)),
      as.character(isTRUE(d_glm_wt1$converged))
    )
  )
  
  ## Construct table with % of clones artificially censored at each interval during the grace period
  # Calculate instantaneous censoring percentage at each time interval
  censoring_by_time <- outcome_data_join %>%
    group_by(assign) %>%
    mutate(n_clones = n_distinct({{ id_var }}[t_intrv == 0])) %>%
    group_by(assign, t_intrv) %>%
    summarize(
      n_clones = first(n_clones),
      n_censored = sum(artificial_cens == 1),
      pct_censored = n_censored / n_clones * 100,
      .groups = "drop"
    ) %>%
    arrange(assign, t_intrv)
  
  # make wide format
  censoring_table_wide <- censoring_by_time %>%
    filter(t_intrv <= grace_intrv_end) %>%
    select(assign, t_intrv, pct_censored) %>%
    pivot_wider(
      names_from = t_intrv,
      values_from = pct_censored,
      names_prefix = "t"
    ) %>%
    mutate(
      total_perc_cens_grace = rowSums(pick(starts_with("t")), na.rm = TRUE)
    )
  
  ## Construct table with value of pr_art_cens_n for clones at each interval in the grace period
  unique_pr_art_cens_n_table <- outcome_data_join %>%
    filter(t_intrv <= grace_intrv_end) %>%
    distinct(assign, t_intrv, pr_art_cens_n) %>%
    group_by(assign, t_intrv) %>%
    mutate(pr_art_cens_n_id = paste0("pr_art_cens_n_", row_number())) %>%
    pivot_wider(
        names_from = pr_art_cens_n_id,
        values_from = pr_art_cens_n
    )
    
  ## Construct weight summary distribution tables
  calc_weight_stats <- function(df, weight, group_label, period_label) {
    weights <- df[[weight]]
    
    # Map variable name to descriptive label
    weight_label <- if (weight == "art_cens_ipw") {
      "Nonstabilized"
    } else if (weight == "art_cens_stab_ipw") {
      "Stabilized"
    } else {
      weight # Fallback to raw name if it's a different weight
    }
    
    # Compute percentiles
    p <- quantile(weights, probs = c(0.90, 0.95, 0.99, 0.995, 0.998), na.rm = TRUE)
    
    data.frame(
      Weight = weight_label,
      Period = period_label,
      Group = group_label,
      Min = min(weights, na.rm = TRUE),
      Q1 = quantile(weights, 0.25, na.rm = TRUE),
      Median = median(weights, na.rm = TRUE),
      Mean = mean(weights, na.rm = TRUE),
      Q3 = quantile(weights, 0.75, na.rm = TRUE),
      Max = max(weights, na.rm = TRUE),
      SD = sd(weights, na.rm = TRUE),
      P90 = p["90%"],
      P95 = p["95%"],
      P99 = p["99%"],
      P995 = p["99.5%"],
      P998 = p["99.8%"],
      check.names = FALSE
    )
  }
  
  interval_after_grace <- outcome_data_join %>%
    filter(t_intrv == (grace_intrv_end + 1))
  
  df_weight_summary <- bind_rows(
    # Weight distribution across follow-up (nonstabilized)
    calc_weight_stats(outcome_data_join, "art_cens_ipw", "Overall", "All Follow-up"),
    calc_weight_stats(filter(outcome_data_join, assign == 0), "art_cens_ipw", "Assign = 0", "All Follow-up"),
    calc_weight_stats(filter(outcome_data_join, assign == 1), "art_cens_ipw", "Assign = 1", "All Follow-up"),
    # Weight distribution in interval after grace period (nonstabilized)
    calc_weight_stats(interval_after_grace, "art_cens_ipw", "Overall", "Interval after end of grace period"),
    calc_weight_stats(filter(interval_after_grace, assign == 0), "art_cens_ipw", "Assign = 0", "Interval after end of grace period"),
    calc_weight_stats(filter(interval_after_grace, assign == 1), "art_cens_ipw", "Assign = 1", "Interval after end of grace period"),
    # Weight distribution across follow-up (stabilized)
    calc_weight_stats(outcome_data_join, "art_cens_stab_ipw", "Overall", "All Follow-up"),
    calc_weight_stats(filter(outcome_data_join, assign == 0), "art_cens_stab_ipw", "Assign = 0", "All Follow-up"),
    calc_weight_stats(filter(outcome_data_join, assign == 1), "art_cens_stab_ipw", "Assign = 1", "All Follow-up"),
    # Weight distribution in interval after grace period (stabilized)
    calc_weight_stats(interval_after_grace, "art_cens_stab_ipw", "Overall", "Interval after end of grace period"),
    calc_weight_stats(filter(interval_after_grace, assign == 0), "art_cens_stab_ipw", "Assign = 0", "Interval after end of grace period"),
    calc_weight_stats(filter(interval_after_grace, assign == 1), "art_cens_stab_ipw", "Assign = 1", "Interval after end of grace period"),
  )
    
  ## Weight sums relative to baseline N and observed risk-set N (during grace period + interval after grace period)
  baseline_n <- outcome_data_join %>%
    filter(t_intrv == 0) %>%
    group_by(assign) %>%
    summarize(
      n_baseline = n_distinct({{ id_var }}),
      .groups = "drop"
    )
  
  weight_sums <- outcome_data_join %>%
    filter(t_intrv <= (grace_intrv_end + 1)) %>%
    group_by(assign, t_intrv) %>%
    summarize(
      sum_unstab_weights = sum(art_cens_ipw, na.rm = TRUE),
      sum_stab_weights = sum(art_cens_stab_ipw, na.rm = TRUE),
      n_remaining=n_distinct({{ id_var }}),
      .groups = "drop"
    )
  
  weight_check_table <- weight_sums %>%
    left_join(baseline_n, by = "assign") %>%
    mutate(
      ratio_to_baseline_unstab = sum_unstab_weights / n_baseline,
      ratio_to_remaining_stab  = sum_stab_weights / n_remaining
    ) %>%
    arrange(assign, t_intrv)

  ## Bring model summaries together into one table
  # Tidy each model
  t0_n <- tidy(n_glm_wt0) %>%
    mutate(Model = "Assign = 0 Numerator Model \n(estimate [se])")
  t0_d <- tidy(d_glm_wt0) %>%
    mutate(Model = "Assign = 0 Denominator Model \n(estimate [se])")
  
  t1_n <- tidy(n_glm_wt1) %>%
    mutate(Model = "Assign = 1 Numerator Model \n(estimate [se])")
  t1_d <- tidy(d_glm_wt1) %>%
    mutate(Model = "Assign = 1 Denominator Model \n(estimate [se])")

  # Combine
  all_models <- bind_rows(t0_n, t0_d, t1_n, t1_d)
    
  # Create "Estimate (SE)" column
  all_models <- all_models %>%
    mutate(est_se = sprintf("%.3f (%.3f)", estimate, std.error)) %>%
    select(term, Model, est_se)
  
  # Pivot to wide format
  model_output <- all_models %>%
    pivot_wider(names_from = Model,
                values_from = est_se)

  ## Build and style Excel workbook
  wb <- createWorkbook()
  
  header_style <- createStyle(
    fontName= "Calibri", fontSize = 11, fontColour = "#FFFFFF",
    fgFill = "#2EC4B6", halign = "center", textDecoration = "bold"
  )
  title_style <- createStyle(
    fontName= "Calibri", fontSize = 14, fontColour = "#1A2530", textDecoration = "bold"
  )
  
  ## Construct df of select IDs from outcome_data_join to make sure weights were lagged and the cumulative product was taken correctly
  sample_ids <- outcome_data_join %>% 
    pull({{ id_var }}) %>%
    unique() %>%
    head(15) # select 15 unique IDs
  
  outcome_data_join_check <- outcome_data_join %>%
    filter({{ id_var }} %in% sample_ids) %>%
    select(
      {{ id_var }},
      assign,
      t_intrv,
      artificial_cens,
      art_cens_wt,
      art_cens_stab_wt,
      lagged_art_wt,
      lagged_art_stab_wt,
      art_cens_ipw,
      art_cens_stab_ipw
    )
  
  ## Construct df with clones that have nonstabilized weights in the 99.8th percentile or higher
  model_covs_clean <- covs %>%
    str_remove_all("factor\\(") %>%
    str_remove_all("\\)")
  
  df_top998 <- outcome_data_join %>%
    group_by(assign) %>%
    filter(art_cens_ipw >= quantile(art_cens_ipw, 0.998, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by({{ id_var }}, assign) %>%
    slice_max(art_cens_ipw, n = 1, with_ties = FALSE) %>% # Ensures you only take one row per unique ID
    ungroup() %>%
    arrange(assign) %>%
    select(
      {{ id_var }},
      assign,
      t_intrv,
      art_cens_ipw,
      all_of(model_covs_clean)
    )
  
  ## Construct df with clones that have stabilized weights in the 99.8th percentile or higher
  df_top998_stab <- outcome_data_join %>%
    group_by(assign) %>%
    filter(art_cens_stab_ipw >= quantile(art_cens_stab_ipw, 0.998, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by({{ id_var }}, assign) %>%
    slice_max(art_cens_stab_ipw, n = 1, with_ties = FALSE) %>% # Ensures you only take one row per unique ID
    ungroup() %>%
    arrange(assign) %>%
    select(
      {{ id_var }},
      assign,
      t_intrv,
      art_cens_stab_ipw,
      all_of(model_covs_clean)
    )
  
  ## Put all the sheets together
  
  # Run summary sheet
  addWorksheet(wb, "Run Summary")
  
  writeData(wb, "Run Summary", "Model Configuration & Convergence", startRow = 1, startCol = 1)
  addStyle(wb, "Run Summary", title_style, rows = 1, cols = 1)
  
  writeData(wb, "Run Summary", metadata_df, startRow = 2, startCol = 1, headerStyle = header_style)
  
  setColWidths(wb, "Run Summary", cols = 1:2, widths = c(30, 45))
  
  # Percent of clones artificially censored during grace period sheet
  addWorksheet(wb, "% Artificially Censored")
  
  writeData(wb, "% Artificially Censored", "% of Clones Artificially Censored (relative to no. of clones at baseline)", startRow = 1, startCol = 1)
  addStyle(wb, "% Artificially Censored", title_style, rows = 1, cols = 1)
  
  writeData(wb, "% Artificially Censored", censoring_table_wide, startRow = 2, startCol = 1, headerStyle = header_style)
  
  # Unique values of pr_art_cens_n sheet
  addWorksheet(wb, "pr_art_cens_n Values")
  
  writeData(wb, "pr_art_cens_n Values", "Unique values of pr_art_cens_n, by assign and t_intrv", startRow = 1, startCol = 1)
  addStyle(wb, "pr_art_cens_n Values", title_style, rows = 1, cols = 1)
  
  writeData(wb, "pr_art_cens_n Values", unique_pr_art_cens_n_table, startRow = 2, startCol = 1, headerStyle = header_style)
  
  # Weight distribution sheet
  addWorksheet(wb, "Weight Distributions")
  
  writeData(wb, "Weight Distributions", "Weight Distribution (Overall and By Treatment Strategy)", startRow = 1, startCol = 1)
  addStyle(wb, "Weight Distributions", title_style, rows = 1, cols = 1)
  
  writeData(wb, "Weight Distributions", df_weight_summary, startRow = 2, startCol = 1, headerStyle = header_style)
  
  setColWidths(wb, "Weight Distributions", cols = 1:2, widths = c(45, 30))
  
  # Sum of weights sheet
  addWorksheet(wb, "Sum of Weights")
  
  writeData(wb, "Sum of Weights", "Weight sums relative to baseline N and observed risk-set N", startRow = 1, startCol = 1)
  addStyle(wb, "Sum of Weights", title_style, rows = 1, cols = 1)
  writeData(wb, "Sum of Weights", "Nonstabilized weights: A ratio close to 1 indicates that weighting approximately recovers the original baseline population size; outcome events and natural censoring will also reduce the risk set", startRow=2, startCol=1, headerStyle = header_style)
  writeData(wb, "Sum of Weights", "Stabilized weights: Because stabilization is intended to preserve the size of the observed risk set, the ratio of the sum of stabilized weights to the number remaining should generally be close to 1", startRow=3, startCol=1, headerStyle=header_style)
  
  writeData(wb, "Sum of Weights", weight_check_table, startRow = 4, startCol = 1, headerStyle = header_style)
  
  # Model output sheet
  addWorksheet(wb, "Model Output")
  
  writeData(wb, "Model Output", "Assign = 0 and Assign = 1 Model Output", startRow = 1, startCol = 1)
  addStyle(wb, "Model Output", title_style, rows = 1, cols = 1)
  
  writeData(wb, "Model Output", model_output, startRow = 2, startCol = 1, headerStyle = header_style)
  
  setColWidths(wb, "Model Output", cols = 1:3, widths = c(45, 30, 30))
    
  # Weight lag check sheet
  addWorksheet(wb, "Check Weight Lag")
  
  writeData(wb, "Check Weight Lag", "Check Lag and Cumulative Product of Weights for Select IDs", startRow = 1, startCol = 1)
  addStyle(wb, "Check Weight Lag", title_style, rows = 1, cols = 1)
  
  writeData(wb, "Check Weight Lag", outcome_data_join_check, startRow = 2, startCol = 1, headerStyle = header_style)
  
  # Extreme weight sheet (nonstab)
  addWorksheet(wb, "Extreme Weights (Nonstab)")
  
  writeData(wb, "Extreme Weights (Nonstab)", "Clones with nonstabilized weights in the 99.8th percentile", startRow = 1, startCol = 1)
  addStyle(wb, "Extreme Weights (Nonstab)", title_style, rows = 1, cols = 1)
  
  writeData(wb, "Extreme Weights (Nonstab)", df_top998, startRow = 2, startCol = 1, headerStyle = header_style)

  # Extreme weight sheet (stab)
  addWorksheet(wb, "Extreme Weights (Stab)")
  
  writeData(wb, "Extreme Weights (Stab)", "Clones with stabilized weights in the 99.8th percentile", startRow = 1, startCol = 1)
  addStyle(wb, "Extreme Weights (Stab)", title_style, rows = 1, cols = 1)
  
  writeData(wb, "Extreme Weights (Stab)", df_top998_stab, startRow = 2, startCol = 1, headerStyle = header_style)
  
  # Save Excel
  saveWorkbook(wb, file = file.path(dir, "summary_art_wt_estimation.xlsx"), overwrite = TRUE)

  ## Return person-period dataset with art_cens_ipw for in-session inspection
  return(outcome_data_join)

}

## Call function

# Example function call (all-cause mortality)
death_art_weight_est <- fit_artificial_cens_weights(
  data_in_path="//your/path/covs_pp_outcome_death.sas7bdat", # cloned, person-period dataset with covariates included
  outcome_name = "death", # used to name folder, datasets, and figures
  covs = c("factor(cov_pain_max)", "factor(cov_prior_hosp_count2)", "factor(cov_ADRD)", "factor(cov_sex)", "cov_Gagne", "factor(cov_Heart_Failure_Dis)",
           "factor(cov_COPD)", "factor(cov_adv_kid)", "factor(cov_Other_musculo_pain_s)", "factor(cov_loop)", "factor(cov_o2sat_min)", "age_at_start", "factor(cov_Chronic_pain_other)",
           "factor(cov_calendar_year)", "factor(cov_Spine_degen_dis)", "factor(cov_Neuropathic_pain)", "factor(cov_urbrrl_snf)", "factor(cov_antipsych)", "factor(cov_Diabetes)",
           "factor(cov_cardioresp)", "factor(cov_nsaid)", "factor(cov_hyperpolypharm)"), # vector of covariates. Write factor() around any categorical covariates 
  id_var = id_var, # this is the name of the unique ID variable for each eligible patient
  time_func_assign0 = "factor(t_intrv)", # Time function for assign = 0 model
  time_func_assign1 = , # If artificial censoring only occurs in the last interval of the grace period, leave time_func_ empty
  grace_intrv_end = 2,  # last interval in the grace period (e.g., t_intrv = 2)
  model_counter = 1,    # used to name folder
  save_path = "//your/path"
)
