################################################################################################################################
# Purpose: Fit weighted pooled logistic regression models to estimate outcome risks, risk differences, and risk ratios.
# 
# Input dataset:
#
#   Cloned, person-period dataset with covariates and artificial censoring weights added.
#     - Generated from 00_artificial_censoring_weight_estimation.R.
#
# Program notes:
#   - This program fits separate outcome models by treatment strategy
#   - When you run the function (estimate_weighted_plr), it fits the PLR models 2 ways:
#     - Weighted without covariates in the model formula.
#     - Weighted WITH baseline covariates in the model formula.
#       - When we add baseline covariates to the model formula, we standardize counterfactual predictions over the 
#         empirical distribution of these covariates.
#		- In the manuscript, we reported results WITH baseline covariates in the model formula.
#
################################################################################################################################

## Import libraries
library(haven)
library(openxlsx)
library(tidyverse)
library(dplyr)
library(broom)
library(modelsummary)
library(rlang)
library(fs)
library(here)
library(purrr)

## Function for outcome estimation
estimate_weighted_plr <- function(data_in_path, covs, id_var, weight_var, rhs0, rhs1, 
                                  grace_intrv_end, pre_grace_length, post_grace_length, outcome_name, model_counter, save_path) {
  
  ## Create directory
  dir <- file.path(save_path, paste0(outcome_name, "/01_outcome_estimation/model_", model_counter, "/"))
  fs::dir_create(dir)
  
  ## Read in person-period dataset with covariates and artificial censoring weights added
  outcome_dataset <- readRDS(data_in_path)
  
  ## Build and style Excel workbook
  wb <- createWorkbook()
  
  header_style <- createStyle(
    fontName= "Calibri", fontSize = 11, fontColour = "#FFFFFF",
    fgFill = "#2EC4B6", halign = "center", textDecoration = "bold"
  )
  title_style <- createStyle(
    fontName= "Calibri", fontSize = 14, fontColour = "#1A2530", textDecoration = "bold"
  )
  
  ### Construct formula for models
  ## Make assign a categorical variable
  outcome_dataset$assign_cat <- factor(
    outcome_dataset$assign,
    levels = c(0, 1) 
  )
    
  ## Prep weight variable
  weight_str <- deparse(substitute(weight_var))
  
  weight_type <- if (weight_str == "art_cens_ipw") {
    "Nonstabilized"
  } else if (weight_str == "art_cens_stab_ipw") {
    "Stabilized"
  } else {
    "Other"
  }
  
  ## Assign = 0 model (opioids-only strategy)
  # Restrict dataset to assign = 0
  outcome_dataset_assign0 <- outcome_dataset %>%
    filter(assign_cat == "0")
  
  # No covariates
  form0 <- as.formula(rhs0)
  
  # With covariates
  # Collapse covariate names into one string
  cov_string <- paste(covs, collapse = " + ")
  
  # Combine
  form_covs0 <- update(rhs0, paste(". ~ . +", cov_string))

  ## Assign = 1 model (gabapentinoids + opioids strategy)
  # Restrict dataset to assign = 1
  outcome_dataset_assign1 <- outcome_dataset %>%
    filter(assign_cat == "1")
  
  # No covariates
  form1 <- as.formula(rhs1)
  
  # With covariates
  form_covs1 <- update(rhs1, paste(". ~ . +", cov_string))
    
  ## Outcome estimation
  fit_outcome_specs <- function(outcome_model_form0, outcome_model_form1, include_covs) {
    
    # Initialize empty vectors to store warnings
    model_warnings0 <- character() # Assign = 0 model
    model_warnings1 <- character() # Assign = 1 model
    
    ## Assign = 0 Outcome model (opioids-only strategy)
    outcome_dataset_assign0$.weights <- outcome_dataset_assign0[[weight_str]] 
    
    d_glm_pe0 <- withCallingHandlers(
      glm(
        outcome_model_form0,
        data = outcome_dataset_assign0,
        family = binomial(),
        weights = .weights
      ),
      warning = function(w) {
        # Store the warning message string
        model_warnings0 <<- c(model_warnings0, w$message)
      }
    )
    
    ## Assign = 1 Outcome model (gabapentinoids + opioids strategy)
    outcome_dataset_assign1$.weights <- outcome_dataset_assign1[[weight_str]] 
       
    d_glm_pe1 <- withCallingHandlers(
      glm(
       outcome_model_form1,
        data = outcome_dataset_assign1,
        family = binomial(),
        weights = .weights
      ),
      warning = function(w) {
        # Store the warning message string
        model_warnings1 <<- c(model_warnings1, w$message)
      }
    )
        
    ## Construct baseline pseudo-population
    # Extract baseline covariates for unique residents
      # To avoid duplicate rows, we only keep one clone per resident
    baseline_cohort <- outcome_dataset %>%
      filter(t_intrv == 0 & assign_cat == "0") %>%
      select(-t_intrv, -assign_cat)
    
    # Expand baseline cohort across all potential follow-up time intervals to create a fully balanced grid for counterfactual prediction
    rectangular_base <- crossing(
      baseline_cohort,
      t_intrv = min(outcome_dataset$t_intrv):max(outcome_dataset$t_intrv)
    )
    
    ## Counterfactual assignment and outcome prediction
    # Force all residents to have assign = 1
    data_all_1 <- rectangular_base %>% mutate(assign_cat = "1")
    # Force all residents to have assign = 0
    data_all_0 <- rectangular_base %>% mutate(assign_cat = "0")
    
    # Predict discrete interval-specific survival probabilities under each treatment strategy
    rectangular_base$pr_surv_1 = predict(d_glm_pe1, newdata = data_all_1, type = "response")
    rectangular_base$pr_surv_0 = predict(d_glm_pe0, newdata = data_all_0, type = "response")
        
    indiv_counterfactual_surv <- rectangular_base %>%
      # Sort by each resident's follow-up
      arrange({{ id_var }}, t_intrv) %>%
      
      # Calculate cumulative survival under each treatment strategy per resident
      group_by({{ id_var }}) %>%
      mutate(
        pr_cum_1 = cumprod(pr_surv_1), # Cumulative survival under assign = 1
        pr_cum_0 = cumprod(pr_surv_0)  # Cumulative survival under assign = 0
      ) %>%
      ungroup()
    
    # Average survival probability by treatment and time interval
    marginal_survival_by_time <- indiv_counterfactual_surv %>%
      group_by(t_intrv) %>%
      summarize(
        surv_1 = mean(pr_cum_1, na.rm = TRUE),
        surv_0 = mean(pr_cum_0, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Derive cumulative incidence, risk difference, and risk ratio
    final_effect_estimates <- marginal_survival_by_time %>%
      mutate(
        # Cumulative incidence under each treatment strategy
        pr_ev_1 = 1 - surv_1,
        pr_ev_0 = 1 - surv_0,
        # Risk difference
        cid = pr_ev_1 - pr_ev_0,
        # Risk ratio
        cir = pr_ev_1 / pr_ev_0,
        Est = paste("PLR-Weights", include_covs, sep = "-")
      )
  
    # Scale to day-level
    final_effect_estimates2 <- final_effect_estimates %>%
      mutate(
        # Calculate the exact last day of the grace period (0-indexed)
        grace_end_day = (grace_intrv_end + 1) * pre_grace_length - 1,
        time = ifelse(
          t_intrv <= grace_intrv_end,
          # Pre-grace: (t_intrv + 1) * length - 1
          (t_intrv + 1) * pre_grace_length - 1,
          # Post-grace: add post_grace_length for each interval after grace period
          grace_end_day + (t_intrv - grace_intrv_end) * post_grace_length
        ),
        # Cap maximum follow-up duration
        time = pmin(time, 98) # if time is greater than the end of follow-up, cap at the last day of follow-up
      ) %>%
      select(-grace_end_day) # drop temporary column
    
    ### Output info to Excel tabs
    
    ## Metadata
    id_str <- rlang::as_name(rlang::enquo(id_var))
    
    metadata_df <- data.frame(
      Parameter = c(
        "Models Run At",
        "Outcome",
        "ID Variable",
        "Weight Variable",
        "Grace Interval End",
        "Covariates in the Model?",
        "Assign = 0 Model Converged?",
        "Assign = 0 Model Warnings?",
        "Assign = 1 Model Converged?",
        "Assign = 1 Model Warnings?"
      ),
      Value = c(
        format(Sys.time(), "%B %d, %Y at %I:%M %p"),
        outcome_name,
        id_str,
        weight_str,
        as.character(grace_intrv_end),
        include_covs,
        as.character(isTRUE(d_glm_pe0$converged)),
        if (length(model_warnings0) > 0) paste(model_warnings0, collapse = " | ") else "No",
        as.character(isTRUE(d_glm_pe1$converged)),
        if (length(model_warnings1) > 0) paste(model_warnings1, collapse = " | ") else "No"
      )
    )
        
    # Output metadata to Excel
    addWorksheet(wb, paste0("Run Summary (", include_covs, ")"))
    
    writeData(wb, paste0("Run Summary (", include_covs, ")"), paste0("Model Configuration & Convergence (", include_covs, ")"), startRow = 1, startCol = 1)
    addStyle(wb, paste0("Run Summary (", include_covs, ")"), title_style, rows = 1, cols = 1)
    
    writeData(wb, paste0("Run Summary (", include_covs, ")"), metadata_df, startRow = 2, startCol = 1, headerStyle = header_style)
    
    setColWidths(wb, paste0("Run Summary (", include_covs, ")"), cols = 1:2, widths = c(30, 45))
    
    ## Model summary
    # Tidy each model
    coef_df0 <- tidy(d_glm_pe0) %>%
      mutate(Model = "Assign = 0 Outcome Model")
    coef_df1 <- tidy(d_glm_pe1) %>%
      mutate(Model = "Assign = 1 Outcome Model")
    
    # Combine
    all_models <- bind_rows(coef_df0, coef_df1)
        
    # Create "Estimate (SE)" column
    all_models <- all_models %>%
      mutate(est_se = sprintf("%.3f (%.3f)", estimate, std.error)) %>%
      select(term, Model, est_se)
    
    # Pivot to wide format
    model_output <- all_models %>%
      pivot_wider(names_from = Model,
                  values_from = est_se)
    
    # Output model summary to Excel
    addWorksheet(wb, paste0("Model Summary (", include_covs, ")"))
    
    writeData(wb, paste0("Model Summary (", include_covs, ")"), paste0("Model Coefficient Estimates and SE (", include_covs, ")"), startRow = 1, startCol = 1)
    addStyle(wb, paste0("Model Summary (", include_covs, ")"), title_style, rows = 1, cols = 1)
    
    writeData(wb, paste0("Model Summary (", include_covs, ")"), model_output, startRow = 2, startCol = 1, headerStyle = header_style)
    
    ## Distribution of predicted survival probabilities for each treatment strategy (should be between 0 and 1)
    calc_pr_stats <- function(assign) {
      pr_surv <- rectangular_base[[paste0("pr_surv_", assign)]]
      
      data.frame(
        Treatment_Strategy = assign,
        Min = min(pr_surv, na.rm = TRUE),
        Q1 = quantile(pr_surv, 0.25, na.rm = TRUE),
        Median = median(pr_surv, na.rm = TRUE),
        Mean = mean(pr_surv, na.rm = TRUE),
        Q3 = quantile(pr_surv, 0.75, na.rm = TRUE),
        Max = max(pr_surv, na.rm = TRUE),
        SD = sd(pr_surv, na.rm = TRUE)
      )
    }
    
    df_pr_surv_summary <- bind_rows(
      calc_pr_stats("1"),
      calc_pr_stats("0")
    )
        
    # Output distribution of predicted survival probabilities to Excel
    addWorksheet(wb, paste0("Predicted Probs (", include_covs, ")"))
    
    writeData(wb, paste0("Predicted Probs (", include_covs, ")"), 
              paste0("Distribution of predicted survival probabilities under each treatment strategy (", include_covs, ")"), startRow = 1, startCol = 1)
    addStyle(wb, paste0("Predicted Probs (", include_covs, ")"), title_style, rows = 1, cols = 1)
    
    writeData(wb, paste0("Predicted Probs (", include_covs, ")"), df_pr_surv_summary, startRow = 2, startCol = 1, headerStyle = header_style)
    
    ### Return dataset with final effect estimates
    return(final_effect_estimates2)
    
  }
  
  # Call function
  model_nocovs <- fit_outcome_specs(
    outcome_model_form0=form0,
    outcome_model_form1=form1,
    include_covs="nocovs"
  )
  model_covs <- fit_outcome_specs(
    outcome_model_form0=form_covs0,
    outcome_model_form1=form_covs1,
    include_covs="covs"
  )
    
  ## Save Excel and datasets
  
  # Save weighted PLR results without covariates
  saveRDS(model_nocovs,
          file = file.path(dir, paste0("wt_plr_nocovs_wide_", outcome_name, ".rds")))
  
  # Save weighted PLR results with covariates
  saveRDS(model_covs,
          file = file.path(dir, paste0("wt_plr_covs_wide_", outcome_name, ".rds")))
  
  # Save Excel
  saveWorkbook(wb, file = file.path(dir, "summary_outcome_estimation.xlsx"), overwrite = TRUE)

  ## Return weighted PLR results for in-session inspection
  return(list(model_nocovs=model_nocovs, model_covs=model_covs))
  
}

## Call function

## Example function calls
# Nonstabilized weights 
death_nonstab_model1 <- estimate_weighted_plr(
  data_in_path = "//your/path/pp_w_art_weights_death.rds", # Dataset from 00_artificial_censoring_weight_estimation.R
  covs = c("factor(cov_pain_max)", "factor(cov_prior_hosp_count2)", "factor(cov_ADRD)", "factor(cov_sex)", "cov_Gagne", "factor(cov_Heart_Failure_Dis)",
           "factor(cov_COPD)", "factor(cov_adv_kid)", "factor(cov_Other_musculo_pain_s)", "factor(cov_loop)", "factor(cov_o2sat_min)", "age_at_start", "factor(cov_Chronic_pain_other)",
           "factor(cov_calendar_year)", "factor(cov_Spine_degen_dis)", "factor(cov_Neuropathic_pain)", "factor(cov_urbrrl_snf)", "factor(cov_antipsych)", "factor(cov_Diabetes)",
           "factor(cov_cardioresp)", 
           "factor(cov_nsaid)", "factor(cov_hyperpolypharm)"), # baseline covariates you want to standardize over
  id_var = id_var, # this is the name of the unique ID variable for each eligible patient
  weight_var = art_cens_ipw,   # art_cens_ipw for unstabilized weights; art_cens_stab_ipw for stabilized weights
  rhs0 = event == 0 ~ splines::ns(t_intrv, knots = c(2, 6, 17)), # time function for assign = 0 model
  rhs1 = event == 0 ~ splines::ns(t_intrv, knots = c(2, 6, 17)), # time function for assign = 1 model
  grace_intrv_end = 2,   # last interval in the grace period (e.g., t_intrv = 2)
  pre_grace_length = 1,  # length of discrete-time intervals prior to the end of the grace period (e.g., 1-day intervals)
  post_grace_length = 3, # length of discrete-time intervals after the end of the grace period (e.g., 3-day intervals)
  outcome_name = "death",# used to name folder, datasets, and figures
  model_counter = 1,     # used to name folder
  save_path = "//your path"
)

# Stabilized weights (all-cause mortality)
death_stab_model2 <- estimate_weighted_plr(
  data_in_path = "//your/path/pp_w_art_weights_death.rds", # Dataset from 00_artificial_censoring_weight_estimation.R
  covs = c("factor(cov_pain_max)", "factor(cov_prior_hosp_count2)", "factor(cov_ADRD)", "factor(cov_sex)", "cov_Gagne", "factor(cov_Heart_Failure_Dis)",
           "factor(cov_COPD)", "factor(cov_adv_kid)", "factor(cov_Other_musculo_pain_s)", "factor(cov_loop)", "factor(cov_o2sat_min)", "age_at_start", "factor(cov_Chronic_pain_other)",
           "factor(cov_calendar_year)", "factor(cov_Spine_degen_dis)", "factor(cov_Neuropathic_pain)", "factor(cov_urbrrl_snf)", "factor(cov_antipsych)", "factor(cov_Diabetes)",
           "factor(cov_cardioresp)", 
           "factor(cov_nsaid)", "factor(cov_hyperpolypharm)"), # baseline covariates you want to standardize over
  id_var = id_var, # this is the name of the unique ID variable for each eligible patient
  weight_var = art_cens_stab_ipw,   # art_cens_ipw for unstabilized weights; art_cens_stab_ipw for stabilized weights
  rhs0 = event == 0 ~ splines::ns(t_intrv, knots = c(2, 6, 17)), # time function for assign = 0 model
  rhs1 = event == 0 ~ splines::ns(t_intrv, knots = c(2, 6, 17)), # time function for assign = 1 model
  grace_intrv_end = 2,   # last interval in the grace period (e.g., t_intrv = 2)
  pre_grace_length = 1,  # length of discrete-time intervals prior to the end of the grace period (e.g., 1-day intervals)
  post_grace_length = 3, # length of discrete-time intervals after the end of the grace period (e.g., 3-day intervals)
  outcome_name = "death",# used to name folder, datasets, and figures
  model_counter = 2,     # used to name folder
  save_path = "//your/path"
)
