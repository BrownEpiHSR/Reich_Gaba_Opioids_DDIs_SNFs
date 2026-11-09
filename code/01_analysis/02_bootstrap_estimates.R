################################################################################################################################
# Purpose: Generate and save estimates from bootstrap replicates for subsequent calculation of 
#          confidence (credible) intervals.
#
# Input dataset:
#
# Cloned, person-period dataset with covariates (but no artificial censoring weights).
#     - This should be the same dataset you used in the data_in_path argument in 00_artificial_censoring_weight_estimation.R
#
# Program notes:
#
#   - This program uses the same artificial censoring weight estimation and outcome estimation procedures as 00_artificial_censoring_weight_estimation.R
#     and 01_outcome_estimation.R.
#   - This program lets you choose between two bootstrapping methods: Bayesian and percentile-based Poisson bootstrapping.
#     - If you choose Bayesian, you'll technically be calculating credible intervals
#	  - For the Gabapentinoids + Opioids study, we used Bayesian bootstrapping
#   - Warnings are muffled because the job-submission environment used for this study stops on warnings. Convergence failures, fitted-probability
#     warnings, and errors are recorded in the run log; affected replicates are excluded.
#
################################################################################################################################
  
## Suppress warnings so the program runs on the job board
options(warn = -1)

## Import libraries
library(data.table)
library(haven)
library(tidyverse)
library(dplyr)
library(data.table)
library(rlang)
library(modelsummary)
library(fs)
library(here)

## Set seed
set.seed(as.integer(ymd('2026-09-03'))) 

## Function to bootstrap
boot_est <- function(data_in_path, outcome_name, art_cens_covs, out_mod_covs=NULL, id_var, weight_var,
                     art_cens_time_func_assign0=NULL, art_cens_time_func_assign1=NULL, grace_intrv_end, 
                     outc_rhs0, outc_rhs1, include_covs_out_mod, nboot, max_attempts, bayesian_boot=1, model_counter, save_path) {
     
  ## Initialize counters
  success_count <- 0
  attempt_count <- 0
  boot_results <- list() # store results here
  
  ## Create directory
  dir <- file.path(save_path, paste0(outcome_name, "/02_bootstrap_estimates/model_", model_counter, "/"))
  fs::dir_create(dir)
  # Create data directory in dir
  dir_data <- file.path(dir, "data/")
  fs::dir_create(dir_data)
  
  ## Covariate prep for artificial censoring weights
  
  # Collapse covariate names into one string
  art_cov_string <- paste(art_cens_covs, collapse = " + ")
  
  ## Assign = 0 model prep (opioids-only strategy)
  # Determine time component
  time_component0 <- if(!is.null(art_cens_time_func_assign0)) {
    art_cens_time_func_assign0
  } else {
    1
  }
  
  ## Prep numerator model formula object
  f_string_num0 <- paste("artificial_cens ~", time_component0)
  formula_num_obj0 <- as.formula(f_string_num0)
  
  ## Prep denominator model formula object
  # Combine time + covariates 
  time_part0 <- if (time_component0 == 1 && length(art_cens_covs) > 0) NULL else time_component0
  cov_part0  <- if (length(art_cens_covs) > 0) art_cov_string else NULL
  rhs_terms0 <- c(time_part0, cov_part0)
  
  # Construct formula object
  if (length(rhs_terms0) > 0) { # If there's a time function and/or covariates, add a + sign between each one
    rhs_string0 <- paste(rhs_terms0, collapse = " + ")
  } else {
    rhs_string0 <- "1" # Fallback to intercept-only if no time and no covariates
  }
  
  formula_string0 <- paste("artificial_cens ~", rhs_string0) # Add artificial_cens ~ at the beginning
  formula_den_obj0 <- as.formula(formula_string0)            # Turn the string into a formula object
  
  ## Assign 1 model prep (gabapentinoids + opioids strategy)
  # Determine time component
  time_component1 <- if(!is.null(art_cens_time_func_assign1)) {
    art_cens_time_func_assign1
  } else {
    1
  }
  
  ## Prep numerator model formula object
  f_string_num1 <- paste("artificial_cens ~", time_component1)
  formula_num_obj1 <- as.formula(f_string_num1)
  
  # Combine time + covariates
  time_part1 <- if (time_component1 == 1 && length(art_cens_covs) > 0) NULL else time_component1
  cov_part1  <- if (length(art_cens_covs) > 0) art_cov_string else NULL
  rhs_terms1  <- c(time_part1, cov_part1) 
  
  # Construct formula object
  if (length(rhs_terms1) > 0) { # If there's a time function and/or covariates, add a + sign between each one
    rhs_string1 <- paste(rhs_terms1, collapse = " + ")
  } else {
    rhs_string1 <- "1" # Fallback to intercept-only if no time and no covariates
  }
  
  formula_string1 <- paste("artificial_cens ~", rhs_string1) # Add artificial_cens ~ at the beginning
  formula_den_obj1 <- as.formula(formula_string1)            # Turn the string into a formula object
    
  ## Weighted outcome model prep
  
  # Weight variable prep
  weight_str <- deparse(substitute(weight_var))
  
  # Check whether including covariates in the weighted outcome model or not
  # Collapse covariate names into one string
  out_cov_string <- paste(out_mod_covs, collapse = " + ")
  
  if (include_covs_out_mod == 1 & length(out_mod_covs) > 0) {
    # Assign = 0
    rhs_string_out0 <- format(outc_rhs0)
    rhs_final_out0  <- paste(rhs_string_out0, "+", out_cov_string)
    
    # Assign = 1
    rhs_string_out1 <- format(outc_rhs1)
    rhs_final_out1  <- paste(rhs_string_out1, "+", out_cov_string)
  } else {
    # Assign = 0
    rhs_final_out0 <- outc_rhs0
    
    # Assign = 1
    rhs_final_out1 <- outc_rhs1
  }
  
  ## Log file prep
  # Create a clean log name
  log_name <- paste0(dir, "/", outcome_name,
                     if(include_covs_out_mod == 1) "_with_covs" else "_no_covs",
                     "_runlog_", format(Sys.Date(), "%Y%m%d"), ".txt")
  
  # Start sending console output to log file
  log_con <- file(log_name, open = "wt")
  sink(log_con, type = "message") # Captures messages/warnings
  sink(log_con, type = "output")  # Captures printed output
  
  # Prep out_rhs arguments to be printed to log
  outc_rhs1_string <- deparse(substitute(outc_rhs1))
  outc_rhs0_string <- deparse(substitute(outc_rhs0))
  
  # Write the Header
  start_time <- Sys.time()
  cat("=======================================================================\n")
  cat("BOOTSTRAP LOG\n")
  cat("Run Date: ", as.character(start_time), "\n")
  cat("Model Counter: ", model_counter, "\n")
  cat("Outcome: ", outcome_name, "\n")
  cat("Weight variable: ", weight_str, "\n")
  cat("Grace Period Day End: ", grace_intrv_end, "\n")
  cat("Time function for artificial censoring assign = 1 models: ", if(is.null(art_cens_time_func_assign1)) "None" else art_cens_time_func_assign1, "\n")
  cat("Time function for artificial censoring assign = 0 models: ", if(is.null(art_cens_time_func_assign0)) "None" else art_cens_time_func_assign0, "\n")
  cat("Covariates in Weight Model: ", art_cens_covs, "\n")
  cat("Specification of assign = 1 outcome model: ", outc_rhs1_string, "\n")
  cat("Specification of assign = 0 outcome model: ", outc_rhs0_string, "\n")
  cat("Covariates in Outcome Model?", if(include_covs_out_mod == 1) "Yes" else "No", "\n")
  cat("Covariates in Outcome Model: ", out_mod_covs, "\n")
  cat("Bayesian bootstrap?", if(bayesian_boot == 1) "Yes" else "No", "\n")
  cat("Nboot: ", nboot, "\n")
  cat("Max Attempts: ", max_attempts, "\n")
  cat("=======================================================================\n")
  
  ### Read in outcome dataset
  # Extract the file extension and convert to lowercase
  ext <- tolower(tools::file_ext(data_in_path))
  
  # Read file based on extension
  if (ext == "rds") {
    data_in <- readRDS(data_in_path)
    
  } else if (ext == "sas7bdat") {
    data_in <- haven::read_sas(data_in_path)
    
  } else {
    stop("Unsupported file type. Only .rds and .sas7bdat files are accepted.")
  }
  setDT(data_in)
  
  ## Make assign a categorical variable
  data_in$assign_cat <- factor(
    data_in$assign,
    levels = c(0, 1)
  )
  data_in$assign_cat <- relevel(data_in$assign_cat, ref = "0")
  
  ### Bootstrap loop
  while (success_count < nboot && attempt_count < max_attempts) {
    # Initialize separation detection variable (i.e., "glm.fit: fitted probabilities numerically 0 or 1 occurred")
    separation_detected <- FALSE
    
    # Increase attempt count by 1
    attempt_count <- attempt_count + 1
    
    ### Conduct sampling procedure
    if (bayesian_boot == 1) { 
      ## Draw independent exponential values with rate = 1 and normalize them to obtain patient-level Bayesian bootstrap weights that sum to 1
      # To respect cloning, we assign the frequency weight at the person-level
      id_list = distinct(data_in, {{ id_var }}) %>%
        mutate(prior = rexp(n=n(), rate=1),
               freqwt = prior / sum(prior))
      
      bootstrap_data = inner_join(data_in, id_list, by = join_by({{ id_var }}))
    
    } else {
      ## Sample IDs from a Poisson distribution
      # Note that we need to sample by the unique person to account for cloning
      data_in[, freqwt:=rpois(n=1, lambda=1), by=factor({{ id_var }})]
      
      # Drop stays with frequency weight of 0 for efficiency
      bootstrap_data <- data_in[freqwt > 0]
      
    }
  
    ### Artificial censoring weights
    
    iteration_result <- tryCatch({
      
      ### Assign = 0 model (opioids-only strategy)
      # Subset bootstrap dataset to assign = 0 and time period where artificial censoring occurs
      art_cens_0 <- bootstrap_data %>%
        filter(assign == 0,                     # Restrict to assign = 0
               if (time_component0 == 1) {      # If no time function was specified, restrict to the last interval in the grace period
                 t_intrv == grace_intrv_end
               } else {
                 t_intrv <= grace_intrv_end     # If a time function was specified, restrict to all intervals in the grace period
               }
        )
      
      ## Fit models
      # Numerator model
      options(warn = 0)
      n_glm_wt0 <- withCallingHandlers(
        expr = {
          glm(formula_num_obj0, data = art_cens_0, family = binomial(), weights = freqwt)
        },
        warning = function(w) {
          # Check if the warning message contains the specific 'separation' string
          if (grepl("glm.fit: fitted probabilities numerically 0 or 1 occurred", w$message)) {
            separation_detected <<- TRUE # Use <<- to update the variable outside this function
          }
          invokeRestart("muffleWarning")
        }
      )
      options(warn = -1)
    
      # Denominator model
      options(warn = 0)
      d_glm_wt0 <- withCallingHandlers(
        expr = {
          glm(formula_den_obj0, data = art_cens_0, family = binomial(), weights = freqwt) # Note that we've applied freqwt here
        },
        warning = function(w) {
          # Check if the warning message contains the specific 'separation' string
          if (grepl("glm.fit: fitted probabilities numerically 0 or 1 occurred", w$message)) {
            separation_detected <<- TRUE # Use <<- to update the variable outside this function
          }
          invokeRestart("muffleWarning")
        }
      )
      options(warn = -1)
      
      # Retrieve predicted probabilities
      art_cens_0$pr_art_cens_n = n_glm_wt0$fitted.values
      art_cens_0$pr_art_cens_d = d_glm_wt0$fitted.values
      
      ### Assign = 1 model (gabapentinoids + opioids strategy)
      # Subset bootstrap dataset to assign = 1 and time period where artificial censoring occurs
      art_cens_1 <- bootstrap_data %>%
        filter(assign == 1,                     # Restrict to assign = 1
               if (time_component1 == 1) {      # If no time function was specified, restrict to the last interval in the grace period
                 t_intrv == grace_intrv_end
               } else {
                 t_intrv <= grace_intrv_end     # If a time function was specified, restrict to all intervals in the grace period
               }
        )
      
      ## Fit models
      # Numerator model
      options(warn = 0)
      n_glm_wt1 = withCallingHandlers(
        expr = {
          glm(formula_num_obj1, data = art_cens_1, family = binomial(), weights = freqwt) # Note that we've applied freqwt here
        },
        warning = function(w) {
          # Check if the warning message contains the specific 'separation' string
          if (grepl("glm.fit: fitted probabilities numerically 0 or 1 occurred", w$message)) {
            separation_detected <<- TRUE # Use <<- to update the variable outside this function
          }
          invokeRestart("muffleWarning")
        }  
      )
      options(warn = -1)
      
      # Denominator model
      options(warn = 0)
      d_glm_wt1 = withCallingHandlers(
        expr = {
          glm(formula_den_obj1, data = art_cens_1, family = binomial(), weights = freqwt) # Note that we've applied freqwt here
        },
        warning = function(w) {
          # Check if the warning message contains the specific 'separation' string
          if (grepl("glm.fit: fitted probabilities numerically 0 or 1 occurred", w$message)) {
            separation_detected <<- TRUE # Use <<- to update the variable outside this function
          }
          invokeRestart("muffleWarning")
        }  
      )
      options(warn = -1)
      
      # Retrieve predicted probabilities
      art_cens_1$pr_art_cens_n = n_glm_wt1$fitted.values
      art_cens_1$pr_art_cens_d = d_glm_wt1$fitted.values
            
      ## Calculate artificial censoring weights
      
      # Append assign = 0 and assign = 1 artificial censoring datasets
      art_cens_all <- rbind(art_cens_0,
                            art_cens_1)
      
      # Merge art_cens_all with the bootstrap dataset
      outcome_data_join = left_join(bootstrap_data,
                                    select(art_cens_all, {{ id_var }}, t_intrv, pr_art_cens_n, pr_art_cens_d, assign), 
                                    by=join_by({{ id_var }}, t_intrv, assign))
        # time intervals not included in the models above for a given treatment assignment will have pr_art_cens = NA

      # Ensure order of variables before calculating artificial censoring weights in each interval
      setDT(outcome_data_join)
      id_str <- rlang::as_name(rlang::enquo(id_var))
      setorderv(outcome_data_join, cols = c(id_str, "assign", "t_intrv"))
      
      # Define interval condition logic for assign = 0
      cond_0 <- if (!is.null(art_cens_time_func_assign0)) {
        outcome_data_join$t_intrv <= grace_intrv_end # If a time function was specified, calculate weights for every interval in the grace period
      } else {
        outcome_data_join$t_intrv == grace_intrv_end # If a time function was not specified, only calculate weights for the last interval in the grace period
      }
      
      # Define the interval logic for assign = 1
      cond_1 <- if (!is.null(art_cens_time_func_assign1)) {
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
      
      ### Outcome estimation
      
      ## Fit models
      
      ## Assign = 0 model (opioids-only strategy)
      # Restrict to assign = 0
      outcome_data_join_assign0 <- outcome_data_join %>%
        filter(assign_cat == "0")
      
      # Fit model
      outcome_data_join_assign0$.weights <- outcome_data_join_assign0[[weight_str]]
      options(warn = 0)
      d_glm_pe0 <- withCallingHandlers(
        expr = {
          glm(rhs_final_out0, data = outcome_data_join_assign0, family=binomial(), weights = .weights*freqwt) # Note that weights are the product of the artificial censoring weight and freqwt
        },
        warning = function(w) {
          # Check if the warning message contains the specific 'separation' string
          if (grepl("glm.fit: fitted probabilities numerically 0 or 1 occurred", w$message)) {
            separation_detected <<- TRUE # Use <<- to update the variable outside this function
          }
          invokeRestart("muffleWarning")
        }
      ) 
      options(warn = -1)
      
      ## Assign = 1 model (gabapentinoids + opioids strategy)
      # Restrict to assign = 1
      outcome_data_join_assign1 <- outcome_data_join %>%
        filter(assign_cat == "1")
      
      # Fit model
      outcome_data_join_assign1$.weights <- outcome_data_join_assign1[[weight_str]]
      options(warn = 0)
      d_glm_pe1 <- withCallingHandlers(
        expr = {
          glm(rhs_final_out1, data = outcome_data_join_assign1, family=binomial(), weights = .weights*freqwt) # Note that weights are the product of the artificial censoring weight and freqwt
        },
        warning = function(w) {
          # Check if the warning message contains the specific 'separation' string
          if (grepl("glm.fit: fitted probabilities numerically 0 or 1 occurred", w$message)) {
            separation_detected <<- TRUE # Use <<- to update the variable outside this function
          }
          invokeRestart("muffleWarning")
        }
      ) 
      options(warn = -1)
      
      ## Construct baseline pseudo-population
      # Extract baseline covariates for unique residents
      # To avoid duplicate rows, we only keep one clone per resident (and each clone has the same freqwt for a given individual)
      baseline_cohort <- outcome_data_join %>%
        filter(t_intrv == 0 & assign_cat == "0") %>%
        select(-t_intrv, -assign_cat)
      
      # Expand baseline cohort across all potential follow-up time intervals to create a fully balanced grid for counterfactual prediction
      rectangular_base <- crossing(
        baseline_cohort,
        t_intrv = min(outcome_data_join$t_intrv):max(outcome_data_join$t_intrv)
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
          # Take the weighted mean of cumulative survival across the bootstrap sample
          surv_1 = weighted.mean(pr_cum_1, w = freqwt, na.rm = TRUE),
          surv_0 = weighted.mean(pr_cum_0, w = freqwt, na.rm = TRUE),
          .groups = 'drop'
        )
      
      # Derive cumulative incidence, risk difference, and risk ratio in bootstrap sample
      final_effect_estimates <- marginal_survival_by_time %>%
        mutate(
          # Cumulative incidence under each treatment strategy
          pr_ev_1 = 1 - surv_1,
          pr_ev_0 = 1 - surv_0,
          # Risk difference
          cid = pr_ev_1 - pr_ev_0,
          # Risk ratio
          cir = pr_ev_1 / pr_ev_0
        )
      
      ### Check if all models converged
      if (n_glm_wt0$converged && d_glm_wt0$converged && n_glm_wt1$converged && d_glm_wt1$converged && d_glm_pe0$converged && d_glm_pe1$converged) {
        
        if (separation_detected) {
          if (separation_detected) message("Attempt ", attempt_count, "glm.fit: fitted probabilities numerically 0 or 1 occurred: Skipping.")
          
          next
        }
        
        # Update the success count
        success_count <- success_count + 1
        
        # Add the success ID and convert to a standard format
        iteration_data <- final_effect_estimates %>%
          mutate(boot_id = success_count) %>%
          select(boot_id, t_intrv, everything()) # Puts identifiers at the front
        
        # Store results
        boot_results[[success_count]] <- iteration_data
        
        # Only log every 25th success
        if (success_count %% 25 == 0) {
          current_time <- Sys.time()
          message(paste0(">>> Milestone reached: ", success_count, " boots completed at ",
                         format(current_time, "%H:%M:%S")))
        }
        
        # Checkpoint logic: save every 50 successful runs
        if (success_count %% 50 == 0) {
          message(paste("--- CHECKPOINT: Saving success #", success_count, "---"))
          
          # Save as an RDS file
          # Include the success_count in the filename so we don't overwrite
          cov_suffix <- if(include_covs_out_mod == 1) "with_covs" else "no_covs"
          date_stamp <- format(Sys.Date(), "%Y%m%d")
          
          current_filename <- paste0(
            outcome_name,
            "_bootstrap_",
            cov_suffix,
            "_checkpoint_",
            date_stamp,
            "_n",
            success_count,
            ".rds"
          )
          
          full_save_path <- file.path(dir_data, current_filename)
          
          iter_boot_df <- dplyr::bind_rows(boot_results)
          
          saveRDS(iter_boot_df, full_save_path)
        }
        
        } else {
          reasons <- c()
          if (!n_glm_wt0$converged) reasons <- c(reasons, "Numerator Weight Model Assign = 0")
          if (!d_glm_wt0$converged) reasons <- c(reasons, "Denominator Weight Model Assign = 0")
          if (!n_glm_wt1$converged) reasons <- c(reasons, "Numerator Weight Model Assign = 1")
          if (!d_glm_wt1$converged) reasons <- c(reasons, "Denominator Weight Model Assign = 1")     
          if (!d_glm_pe0$converged) reasons <- c(reasons, "Assign = 0 Outcome Model")
          if (!d_glm_pe1$converged) reasons <- c(reasons, "Assign = 1 Outcome Model")
          
          message(paste0("Attempt ", attempt_count, " failed convergence in: ", paste(reasons, collapse = ", ")))
        }
      
      }, error = function(e) {
        message("--- ERROR DETECTED ---")
        print(e) # prints the full error object
        
        return(NULL)
      })
  }
  
  if (success_count < nboot) {
    message("Reached max_attempts before hitting nboot target!")
    
  }
  
  ### Create the final filename
  
  cov_suffix <- if(include_covs_out_mod == 1) "with_covs" else "no_covs"
  date_stamp <- format(Sys.Date(), "%Y%m%d")
  
  final_filename <- paste0(
    outcome_name,
    "_bootstrap_",
    cov_suffix,
    "_FINAL_",
    date_stamp,
    "_n",
    success_count,
    ".rds"
  )
  
  final_save_path <- file.path(dir_data, final_filename)
  
  ### Stack bootstrap results vertically
  final_boot_df <- dplyr::bind_rows(boot_results)
  
  saveRDS(final_boot_df, final_save_path)
  
  # Ending log message
  end_time <- Sys.time()
  cat("=======================================================================\n")
  cat("BOOTSTRAP COMPLETE\n")
  cat("End Time:        ", as.character(end_time), "\n")
  cat("Total Duration:  ", format(end_time - start_time), "\n")
  cat("Total Successes: ", success_count, "\n")
  cat("Total Attempts:  ", attempt_count, "\n")
  cat("=======================================================================\n")
  
  # Close out log file
  sink(type = "message")
  sink(type = "output")
  close(log_con)
  
  message("Analysis complete. Log saved to: ", log_name)
  
  return(final_boot_df)
  
}
  
## Call Function

## Example function call
# All-cause mortality (nonstabilized weights)
death_boot <- boot_est(
  data_in_path = "//your/path/covs_pp_g2_outcome_death.sas7bdat", # cloned, person-period dataset with covariates (can be a SAS or R dataset)
  outcome_name = "death", # used to name folder, datasets, and figures
  art_cens_covs = c("factor(cov_pain_max)", "factor(cov_prior_hosp_count2)", "factor(cov_ADRD)", "factor(cov_sex)", "cov_Gagne", "factor(cov_Heart_Failure_Dis)",
                    "factor(cov_COPD)", "factor(cov_adv_kid)", "factor(cov_Other_musculo_pain_s)", "factor(cov_loop)", "factor(cov_o2sat_min)", "age_at_start", "factor(cov_Chronic_pain_other)",
                    "factor(cov_calendar_year)", "factor(cov_Spine_degen_dis)", "factor(cov_Neuropathic_pain)", "factor(cov_urbrrl_snf)", "factor(cov_antipsych)", "factor(cov_Diabetes)",
                    "factor(cov_cardioresp)", 
                    "factor(cov_nsaid)", "factor(cov_hyperpolypharm)"), # covariates included in artificial censoring weight estimation
  out_mod_covs = c("factor(cov_pain_max)", "factor(cov_prior_hosp_count2)", "factor(cov_ADRD)", "factor(cov_sex)", "cov_Gagne", "factor(cov_Heart_Failure_Dis)",
                   "factor(cov_COPD)", "factor(cov_adv_kid)", "factor(cov_Other_musculo_pain_s)", "factor(cov_loop)", "factor(cov_o2sat_min)", "age_at_start", "factor(cov_Chronic_pain_other)",
                   "factor(cov_calendar_year)", "factor(cov_Spine_degen_dis)", "factor(cov_Neuropathic_pain)", "factor(cov_urbrrl_snf)", "factor(cov_antipsych)", "factor(cov_Diabetes)",
                   "factor(cov_cardioresp)", 
                   "factor(cov_nsaid)", "factor(cov_hyperpolypharm)"),  # baseline covariates included in outcome model formula (you can leave this blank if include_covs_out_mod = 0)
  id_var = id_var, # this is the name of the unique ID variable for each eligible patient
  weight_var = art_cens_ipw,   # art_cens_ipw for nonstabilized weights; art_cens_stab_ipw for stabilized weights
  art_cens_time_func_assign0="factor(t_intrv)", 
  art_cens_time_func_assign1=, # If artificial censoring only occurs in the last interval of the grace period, leave time_func_ empty
  grace_intrv_end=2,  # last interval in the grace period (e.g., t_intrv = 2)
  # Outcome model specifications
  outc_rhs0 = event == 0 ~ splines::ns(t_intrv, knots = c(2, 6, 17)), # Assign = 0 model
  outc_rhs1 = event == 0 ~ splines::ns(t_intrv, knots = c(2, 6, 17)), # Assign = 1 model
  include_covs_out_mod = 1, # 1 if you want to include baseline covariates in the outcome models.
  nboot = 500,              # Total number of successful bootstrap replicates you want.
  max_attempts = 800,       # Max attempts you want to make to get your total number of successful bootstrap replicates.
  bayesian_boot = 1,        # 1 if you want to use Bayesian bootstrap; 0 if you want to assign frequency weights to IDs using a Poisson distribution.
  model_counter = 1,        # used to name folder
  save_path = "//your/path" 
)

