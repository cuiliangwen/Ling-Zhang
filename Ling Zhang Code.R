# ============================================================================== #
# Supplementary R code: s-EASIX and mortality in sepsis-associated AHRF          #
# ============================================================================== #
# This script was reconstructed from the archived dataset, result tables, figures, #
# and analysis reports because the original analysis script was not retained.     #
# Validate every regenerated result against the locked manuscript before upload.  #

# 1. Reproducibility settings --------------------------------------------------- #
set.seed(123) # Use the random seed documented in the archived analysis reports.
options(stringsAsFactors = FALSE) # Prevent automatic character-to-factor conversion.

required_packages <- c( # Declare all packages needed by the analysis workflow.
  "readxl", "dplyr", "tidyr", "purrr", "tibble", "ggplot2", # Data handling and plots.
  "survival", "survminer", "rms", "mice", "car", "pROC", # Epidemiologic analyses.
  "Boruta", "caret", "glmnet", "randomForest", "xgboost", # Feature selection and ML.
  "gbm", "e1071", "rpart" # Gradient boosting, SVM, and decision tree.
) # End the package vector.

missing_packages <- required_packages[ # Identify packages that are not installed.
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE) # Check packages.
] # End the missing-package subset.

if (length(missing_packages) > 0L) { # Stop early if the environment is incomplete.
  stop("Install the following packages before running: ", paste(missing_packages, collapse = ", ")) # Report missing packages.
} # End the package check.

invisible(lapply(required_packages, library, character.only = TRUE)) # Attach required packages.

data_file <- "s-EASIX-汇合（2026-7-6）——原始资料.xlsx" # Set the portable input filename.
output_dir <- "reproduced_outputs" # Set the directory for regenerated files.
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE) # Create the output directory.

# 2. Import the source spreadsheet and create an analysis dataset --------------- #
raw_names <- names(readxl::read_excel(data_file, n_max = 0, .name_repair = "minimal")) # Retain original duplicate headers.
raw_data <- readxl::read_excel(data_file, .name_repair = "minimal") # Read all patient-level records.

pull_occurrence <- function(data, original_names, variable, occurrence = 1L) { # Extract a specified duplicate column.
  matching_columns <- which(original_names == variable) # Locate all columns with the requested header.
  if (length(matching_columns) < occurrence) stop("Column not found: ", variable) # Guard against schema changes.
  data[[matching_columns[[occurrence]]]] # Return the requested occurrence.
} # End the column extractor.

analysis_data <- tibble::tibble( # Build a uniquely named analysis-ready dataset.
  subject_id = pull_occurrence(raw_data, raw_names, "subject_id"), # Patient identifier.
  age = pull_occurrence(raw_data, raw_names, "age"), # Age in years.
  gender = factor(pull_occurrence(raw_data, raw_names, "gender"), levels = c("F", "M")), # Sex.
  race = factor(ifelse(toupper(pull_occurrence(raw_data, raw_names, "race")) == "WHITE", "White", "Other"), levels = c("Other", "White")), # Race group.
  platelet = pull_occurrence(raw_data, raw_names, "lab_platelet_value_seq1", 1L), # Platelet count.
  ldh = pull_occurrence(raw_data, raw_names, "lab_ld_ldh_value_seq1", 1L), # Lactate dehydrogenase.
  s_EASIX = pull_occurrence(raw_data, raw_names, "s-EASIX"), # Archived s-EASIX value.
  charlson_cci = pull_occurrence(raw_data, raw_names, "charlson_cci"), # Charlson comorbidity index.
  myocardial_infarct = pull_occurrence(raw_data, raw_names, "charlson_myocardial_infarct"), # Myocardial infarction.
  heart_failure = pull_occurrence(raw_data, raw_names, "charlson_congestive_heart_failure"), # Heart failure.
  cerebrovascular = pull_occurrence(raw_data, raw_names, "charlson_cerebrovascular"), # Cerebrovascular disease.
  chronic_pulmonary = pull_occurrence(raw_data, raw_names, "charlson_chronic_pulmonary"), # Chronic pulmonary disease.
  diabetes = pull_occurrence(raw_data, raw_names, "charlson_diabetes_all"), # Diabetes mellitus.
  renal_disease = pull_occurrence(raw_data, raw_names, "charlson_renal"), # Chronic renal disease.
  malignant_cancer = pull_occurrence(raw_data, raw_names, "charlson_malignant_cancer"), # Malignant cancer.
  hemoglobin = pull_occurrence(raw_data, raw_names, "lab_hemoglobin_value_seq1"), # Hemoglobin.
  wbc = pull_occurrence(raw_data, raw_names, "lab_wbc_value_seq1"), # White blood cell count.
  creatinine = pull_occurrence(raw_data, raw_names, "lab_creatinine_value_seq1"), # Serum creatinine.
  bun = pull_occurrence(raw_data, raw_names, "lab_bun_value_seq1"), # Blood urea nitrogen.
  albumin = pull_occurrence(raw_data, raw_names, "lab_albumin_value_seq1"), # Serum albumin.
  pt = pull_occurrence(raw_data, raw_names, "lab_pt_value_seq1"), # Prothrombin time.
  ptt = pull_occurrence(raw_data, raw_names, "lab_ptt_value_seq1"), # Activated partial thromboplastin time.
  inr = pull_occurrence(raw_data, raw_names, "lab_inr_value_seq1"), # International normalized ratio.
  bilirubin = pull_occurrence(raw_data, raw_names, "lab_bilirubin_total_value_seq1"), # Total bilirubin.
  alt = pull_occurrence(raw_data, raw_names, "lab_alt_value_seq1"), # Alanine aminotransferase.
  ast = pull_occurrence(raw_data, raw_names, "lab_ast_value_seq1"), # Aspartate aminotransferase.
  lactate = pull_occurrence(raw_data, raw_names, "lab_lactate_value_seq1"), # Lactate.
  pco2 = pull_occurrence(raw_data, raw_names, "lab_pco2_value_seq1", 2L), # Arterial PaCO2.
  ph = pull_occurrence(raw_data, raw_names, "lab_ph_value_seq1", 2L), # Arterial pH.
  pao2 = pull_occurrence(raw_data, raw_names, "lab_po2_value_seq1", 2L), # Arterial PaO2.
  pao2_fio2 = pull_occurrence(raw_data, raw_names, "lab_pao2fio2ratio_value_seq1"), # PaO2/FiO2 ratio.
  ventilation = pull_occurrence(raw_data, raw_names, "vent_flag"), # Invasive mechanical ventilation.
  crrt = pull_occurrence(raw_data, raw_names, "crrt_flag"), # Continuous renal replacement therapy.
  sofa = pull_occurrence(raw_data, raw_names, "sofa_score_max"), # SOFA score.
  apsiii = pull_occurrence(raw_data, raw_names, "apsiii_score_max"), # APS III score.
  sapsii = pull_occurrence(raw_data, raw_names, "sapsii_score_max"), # SAPS II score.
  aki = pull_occurrence(raw_data, raw_names, "aki_flag"), # Acute kidney injury.
  weight = pull_occurrence(raw_data, raw_names, "vital_weight_value_min"), # Body weight.
  spo2_min = pull_occurrence(raw_data, raw_names, "vital_spo2_value_min"), # Minimum oxygen saturation.
  temperature_max = pull_occurrence(raw_data, raw_names, "vital_temperature_value_max"), # Maximum temperature.
  heart_rate_max = pull_occurrence(raw_data, raw_names, "vital_heart_rate_value_max"), # Maximum heart rate.
  resp_rate_max = pull_occurrence(raw_data, raw_names, "vital_resp_rate_value_max"), # Maximum respiratory rate.
  vasopressor = pull_occurrence(raw_data, raw_names, "vaso_any_used"), # Any vasopressor use.
  time_30d = pull_occurrence(raw_data, raw_names, "hosp_time_30d"), # Follow-up time for 30-day outcome.
  death_30d = pull_occurrence(raw_data, raw_names, "hosp_outcome_30d"), # 30-day mortality indicator.
  time_60d = pull_occurrence(raw_data, raw_names, "hosp_time_60d"), # Follow-up time for 60-day outcome.
  death_60d = pull_occurrence(raw_data, raw_names, "hosp_outcome_60d"), # 60-day mortality indicator.
  time_90d = pull_occurrence(raw_data, raw_names, "hosp_time_90d"), # Follow-up time for 90-day outcome.
  death_90d = pull_occurrence(raw_data, raw_names, "hosp_outcome_90d", 2L) # 90-day mortality indicator.
) # Finish creation of the analysis dataset.

analysis_data <- analysis_data %>% # Validate and derive exposure categories.
  mutate( # Create derived analysis variables.
    s_EASIX_check = ldh / platelet, # Recalculate s-EASIX for quality control.
    s_EASIX_quartile = factor(ntile(s_EASIX, 4L), levels = 1:4, labels = paste0("Q", 1:4)), # Define quartiles.
    quartile_num = as.integer(s_EASIX_quartile), # Create ordinal quartile score for trend tests.
    across(c(myocardial_infarct, heart_failure, cerebrovascular, chronic_pulmonary, diabetes, renal_disease, malignant_cancer, ventilation, crrt, aki, vasopressor), as.integer) # Standardize binary indicators.
  ) # End derivation.

stopifnot(nrow(analysis_data) == 3220L) # Verify the documented cohort size.
stopifnot(sum(analysis_data$death_30d, na.rm = TRUE) == 1029L) # Verify 30-day deaths.
stopifnot(sum(analysis_data$death_60d, na.rm = TRUE) == 1162L) # Verify 60-day deaths.
stopifnot(sum(analysis_data$death_90d, na.rm = TRUE) == 1239L) # Verify 90-day deaths.

quartile_cutoffs <- quantile(analysis_data$s_EASIX, probs = c(0.25, 0.50, 0.75), na.rm = TRUE) # Calculate quartile cutoffs.
write.csv(data.frame(percentile = names(quartile_cutoffs), cutoff = as.numeric(quartile_cutoffs)), file.path(output_dir, "quartile_cutoffs.csv"), row.names = FALSE) # Export cutoffs.

# 3. Missing-data handling ------------------------------------------------------ #
analysis_variables <- setdiff(names(analysis_data), c("subject_id", "s_EASIX_check")) # Select variables eligible for analysis.
missing_rate <- vapply(analysis_data[analysis_variables], function(x) mean(is.na(x)), numeric(1)) # Calculate missing proportions.
excluded_for_missingness <- names(missing_rate[missing_rate >= 0.30]) # Identify variables with at least 30% missingness.
imputation_variables <- setdiff(analysis_variables, excluded_for_missingness) # Retain variables below the threshold.
imputation_frame <- analysis_data[imputation_variables] # Create the imputation dataset.

imputation_methods <- mice::make.method(imputation_frame) # Initialize variable-specific imputation methods.
imputation_methods[] <- "rf" # Use random-forest imputation for incomplete eligible variables.
imputation_methods[vapply(imputation_frame, function(x) !anyNA(x), logical(1))] <- "" # Skip complete variables.
imputation_methods[grepl("^(time_|death_)", names(imputation_frame))] <- "" # Never impute outcome or follow-up variables.

imputation_predictors <- mice::make.predictorMatrix(imputation_frame) # Initialize the imputation predictor matrix.
imputation_predictors[, grepl("^(time_|death_)", names(imputation_frame))] <- 0L # Exclude outcomes as predictors.
diag(imputation_predictors) <- 0L # Prevent self-prediction.

imputed_object <- mice::mice( # Perform the archived single-completed-dataset workflow.
  imputation_frame, # Supply the incomplete data.
  m = 1L, # Generate one completed dataset to match the archived non-pooled outputs.
  maxit = 5L, # Run five chained-equation iterations.
  method = imputation_methods, # Apply random-forest imputation where needed.
  predictorMatrix = imputation_predictors, # Apply the outcome-excluding predictor matrix.
  seed = 123, # Make imputation reproducible.
  printFlag = FALSE # Suppress iterative console output.
) # End imputation.

completed_data <- mice::complete(imputed_object, action = 1L) # Extract the completed dataset.
completed_data$subject_id <- analysis_data$subject_id # Restore patient identifiers.
completed_data$s_EASIX_check <- analysis_data$s_EASIX_check # Restore the quality-control variable.

# 4. Baseline table by s-EASIX quartile ---------------------------------------- #
continuous_baseline <- c("s_EASIX", "age", "charlson_cci", "sofa", "apsiii", "sapsii", "wbc", "hemoglobin", "platelet", "albumin", "ldh", "lactate", "bilirubin", "alt", "ast", "bun", "creatinine", "inr", "pt", "ptt", "pao2", "pao2_fio2") # Continuous baseline variables.
categorical_baseline <- c("gender", "race", "myocardial_infarct", "heart_failure", "cerebrovascular", "chronic_pulmonary", "diabetes", "renal_disease", "aki", "malignant_cancer", "ventilation", "vasopressor", "crrt", "death_30d", "death_60d", "death_90d") # Categorical baseline variables.

summarize_continuous <- function(data, variable) { # Summarize a continuous variable by quartile.
  data %>% # Start from the analysis data.
    group_by(s_EASIX_quartile) %>% # Stratify by exposure quartile.
    summarise(variable = variable, median = median(.data[[variable]], na.rm = TRUE), q1 = quantile(.data[[variable]], 0.25, na.rm = TRUE), q3 = quantile(.data[[variable]], 0.75, na.rm = TRUE), .groups = "drop") %>% # Calculate median and IQR.
    mutate(p_value = kruskal.test(data[[variable]] ~ data$s_EASIX_quartile)$p.value) # Calculate the Kruskal-Wallis P value.
} # End continuous summarizer.

summarize_categorical <- function(data, variable) { # Summarize a binary/categorical variable by quartile.
  tabulation <- table(data$s_EASIX_quartile, data[[variable]], useNA = "no") # Build the contingency table.
  p_value <- if (any(tabulation < 5L)) fisher.test(tabulation)$p.value else chisq.test(tabulation, correct = FALSE)$p.value # Select Fisher or chi-square test.
  tibble::tibble(variable = variable, p_value = p_value) # Return the categorical test result.
} # End categorical summarizer.

baseline_continuous_results <- purrr::map_dfr(continuous_baseline, ~ summarize_continuous(completed_data, .x)) # Analyze continuous variables.
baseline_categorical_results <- purrr::map_dfr(categorical_baseline, ~ summarize_categorical(completed_data, .x)) # Analyze categorical variables.
write.csv(baseline_continuous_results, file.path(output_dir, "baseline_continuous.csv"), row.names = FALSE) # Export continuous summaries.
write.csv(baseline_categorical_results, file.path(output_dir, "baseline_categorical_tests.csv"), row.names = FALSE) # Export categorical tests.

# 5. Kaplan-Meier analyses ------------------------------------------------------ #
endpoints <- tibble::tribble( # Define the three mortality endpoints.
  ~label, ~time, ~event, # Name endpoint mapping columns.
  "30-day mortality", "time_30d", "death_30d", # Define the 30-day endpoint.
  "60-day mortality", "time_60d", "death_60d", # Define the 60-day endpoint.
  "90-day mortality", "time_90d", "death_90d" # Define the 90-day endpoint.
) # End endpoint definitions.

run_km_analysis <- function(endpoint_label, time_variable, event_variable) { # Fit and plot one Kaplan-Meier analysis.
  survival_formula <- as.formula(paste0("Surv(", time_variable, ", ", event_variable, ") ~ s_EASIX_quartile")) # Build the survival formula.
  km_fit <- survival::survfit(survival_formula, data = completed_data) # Fit quartile-specific survival curves.
  log_rank <- survival::survdiff(survival_formula, data = completed_data) # Perform the log-rank test.
  log_rank_p <- pchisq(log_rank$chisq, df = length(log_rank$n) - 1L, lower.tail = FALSE) # Convert chi-square statistic to P value.
  km_plot <- survminer::ggsurvplot(km_fit, data = completed_data, risk.table = TRUE, pval = TRUE, conf.int = FALSE, xlab = "Follow-up, days", ylab = "Survival probability", legend.title = "s-EASIX quartile", ggtheme = theme_bw()) # Draw the survival plot.
  ggsave(file.path(output_dir, paste0(gsub("[^0-9]", "", endpoint_label), "_KM.png")), km_plot$plot, width = 7, height = 6, dpi = 300) # Export the plot.
  tibble::tibble(outcome = endpoint_label, log_rank_chisq = unname(log_rank$chisq), log_rank_p = log_rank_p) # Return test results.
} # End Kaplan-Meier function.

km_results <- purrr::pmap_dfr(endpoints, run_km_analysis) # Run all Kaplan-Meier analyses.
write.csv(km_results, file.path(output_dir, "kaplan_meier_log_rank.csv"), row.names = FALSE) # Export log-rank results.

# 6. Cox proportional hazards models ------------------------------------------- #
model_2_covariates <- c("age", "gender", "race", "weight", "charlson_cci") # Define the basic adjustment set.
model_3_covariates <- c(model_2_covariates, "myocardial_infarct", "cerebrovascular", "renal_disease", "malignant_cancer", "apsiii") # Define the manuscript multivariable model.

tidy_cox <- function(model, outcome, exposure, adjustment) { # Convert a Cox model to a publication-ready table.
  coefficients <- summary(model)$coefficients # Extract coefficient statistics.
  confidence_intervals <- summary(model)$conf.int # Extract exponentiated confidence intervals.
  tibble::tibble(outcome = outcome, exposure = exposure, model = adjustment, term = rownames(coefficients), HR = confidence_intervals[, "exp(coef)"], CI_low = confidence_intervals[, "lower .95"], CI_high = confidence_intervals[, "upper .95"], p_value = coefficients[, "Pr(>|z|)"], n = model$n, events = model$nevent) # Return tidy results.
} # End Cox tidier.

fit_cox_set <- function(endpoint_label, time_variable, event_variable) { # Fit continuous, quartile, and trend Cox models.
  survival_response <- paste0("Surv(", time_variable, ", ", event_variable, ")") # Construct the survival response.
  adjustment_sets <- list("Model 1" = character(0), "Model 2" = model_2_covariates, "Model 3" = model_3_covariates) # Define all adjustment sets.
  purrr::imap_dfr(adjustment_sets, function(covariates, model_name) { # Iterate over adjustment models.
    suffix <- if (length(covariates) == 0L) "" else paste0(" + ", paste(covariates, collapse = " + ")) # Build the covariate suffix.
    continuous_fit <- survival::coxph(as.formula(paste0(survival_response, " ~ s_EASIX", suffix)), data = completed_data, x = TRUE) # Fit continuous exposure model.
    quartile_fit <- survival::coxph(as.formula(paste0(survival_response, " ~ s_EASIX_quartile", suffix)), data = completed_data, x = TRUE) # Fit quartile exposure model.
    trend_fit <- survival::coxph(as.formula(paste0(survival_response, " ~ quartile_num", suffix)), data = completed_data, x = TRUE) # Fit ordinal trend model.
    bind_rows(tidy_cox(continuous_fit, endpoint_label, "Continuous", model_name), tidy_cox(quartile_fit, endpoint_label, "Quartile", model_name), tidy_cox(trend_fit, endpoint_label, "Trend", model_name)) # Combine exposure parameterizations.
  }) # End adjustment-model iteration.
} # End Cox model-set function.

cox_results <- purrr::pmap_dfr(endpoints, fit_cox_set) # Fit Cox models for all endpoints.
write.csv(cox_results, file.path(output_dir, "cox_models_30_60_90.csv"), row.names = FALSE) # Export Cox estimates.

# 7. Proportional-hazards assumption ------------------------------------------- #
run_ph_test <- function(endpoint_label, time_variable, event_variable) { # Test the fully adjusted Cox model.
  ph_formula <- as.formula(paste0("Surv(", time_variable, ", ", event_variable, ") ~ s_EASIX_quartile + ", paste(model_3_covariates, collapse = " + "))) # Build Model 3 formula.
  ph_model <- survival::coxph(ph_formula, data = completed_data, x = TRUE) # Fit the model used for diagnostics.
  ph_result <- survival::cox.zph(ph_model, transform = "km") # Test scaled Schoenfeld residuals.
  tibble::as_tibble(as.data.frame(ph_result$table), rownames = "term") %>% mutate(outcome = endpoint_label, .before = 1L) # Return diagnostic results.
} # End proportional-hazards test function.

ph_results <- purrr::pmap_dfr(endpoints, run_ph_test) # Test all endpoints.
write.csv(ph_results, file.path(output_dir, "proportional_hazards_tests.csv"), row.names = FALSE) # Export diagnostic tests.

# 8. Restricted cubic splines -------------------------------------------------- #
rms_distribution <- rms::datadist(completed_data) # Store predictor distributions for rms functions.
options(datadist = "rms_distribution") # Register the distribution object.
rcs_knots <- as.numeric(quantile(completed_data$s_EASIX, probs = c(0.05, 0.35, 0.65, 0.95), na.rm = TRUE)) # Place four knots at prespecified percentiles.
rcs_reference <- median(completed_data$s_EASIX, na.rm = TRUE) # Use the median s-EASIX as reference.

run_rcs <- function(endpoint_label, time_variable, event_variable) { # Fit one restricted cubic spline model.
  rcs_formula <- as.formula(paste0("Surv(", time_variable, ", ", event_variable, ") ~ rcs(s_EASIX, rcs_knots) + ", paste(model_3_covariates, collapse = " + "))) # Build the RCS formula.
  rcs_fit <- rms::cph(rcs_formula, data = completed_data, x = TRUE, y = TRUE, surv = TRUE) # Fit the spline Cox model.
  rcs_test <- rms::anova(rcs_fit) # Obtain Wald tests for overall and nonlinear effects.
  prediction_grid <- rms::Predict(rcs_fit, s_EASIX, ref.zero = TRUE, fun = exp) # Estimate adjusted hazard ratios.
  rcs_plot <- ggplot(as.data.frame(prediction_grid), aes(x = s_EASIX, y = yhat)) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#9ECAE1", alpha = 0.40) + geom_line(color = "#8C7A2B", linewidth = 1.1) + geom_hline(yintercept = 1, linetype = 2, color = "grey40") + coord_cartesian(xlim = quantile(completed_data$s_EASIX, c(0.01, 0.99), na.rm = TRUE)) + labs(title = paste("Restricted cubic spline:", endpoint_label), x = "s-EASIX", y = "Hazard ratio (95% CI)") + theme_bw() # Draw the spline curve.
  ggsave(file.path(output_dir, paste0(gsub("[^0-9]", "", endpoint_label), "_RCS.png")), rcs_plot, width = 8, height = 6, dpi = 300) # Export the spline plot.
  tibble::tibble(outcome = endpoint_label, knots = paste(round(rcs_knots, 4), collapse = ", "), reference = rcs_reference, overall_p = rcs_test["s_EASIX", "P"], nonlinear_p = rcs_test[" Nonlinear", "P"]) # Return spline tests.
} # End RCS function.

rcs_results <- purrr::pmap_dfr(endpoints, run_rcs) # Fit all RCS models.
write.csv(rcs_results, file.path(output_dir, "rcs_tests.csv"), row.names = FALSE) # Export RCS test results.

# 9. ROC curves for s-EASIX alone ---------------------------------------------- #
run_single_marker_roc <- function(endpoint_label, event_variable) { # Evaluate s-EASIX alone for one outcome.
  roc_object <- pROC::roc(completed_data[[event_variable]], completed_data$s_EASIX, direction = "<", quiet = TRUE) # Construct the ROC curve.
  auc_ci <- pROC::ci.auc(roc_object, method = "delong") # Calculate the DeLong 95% confidence interval.
  tibble::tibble(outcome = endpoint_label, AUC = as.numeric(pROC::auc(roc_object)), CI_low = as.numeric(auc_ci[1]), CI_high = as.numeric(auc_ci[3])) # Return AUC results.
} # End ROC function.

single_marker_auc <- purrr::map2_dfr(endpoints$label, endpoints$event, run_single_marker_roc) # Analyze all endpoints.
write.csv(single_marker_auc, file.path(output_dir, "sEASIX_auc.csv"), row.names = FALSE) # Export AUC estimates.

# 10. Candidate predictors and Boruta selection -------------------------------- #
candidate_predictors <- c("s_EASIX", "age", "gender", "race", "charlson_cci", "aki", "vasopressor", "hemoglobin", "wbc", "creatinine", "bun", "albumin", "bilirubin", "alt", "ast", "pt", "ptt", "lactate", "ph", "pco2", "pao2_fio2", "weight", "temperature_max", "heart_rate_max", "resp_rate_max") # Define the 25 candidate predictors.

set.seed(123) # Reset the seed before cohort splitting.
training_index <- sample(seq_len(nrow(completed_data)), size = floor(0.70 * nrow(completed_data)), replace = FALSE) # Reproduce the documented non-stratified 7:3 split.
training_data <- completed_data[training_index, , drop = FALSE] # Create the training cohort.
test_data <- completed_data[-training_index, , drop = FALSE] # Create the independent test cohort.

split_summary <- endpoints %>% # Summarize event counts in both cohorts.
  transmute(outcome = label, training_n = nrow(training_data), training_events = map_int(event, ~ sum(training_data[[.x]])), test_n = nrow(test_data), test_events = map_int(event, ~ sum(test_data[[.x]]))) # Calculate counts.
write.csv(split_summary, file.path(output_dir, "train_test_split.csv"), row.names = FALSE) # Export split summary.

run_boruta <- function(event_variable) { # Run Boruta for one binary mortality outcome.
  boruta_frame <- training_data[, c(event_variable, candidate_predictors)] # Assemble outcome and candidate predictors.
  boruta_frame[[event_variable]] <- factor(boruta_frame[[event_variable]], levels = c(0, 1), labels = c("Survived", "Died")) # Convert outcome to a classification factor.
  set.seed(123) # Reset the Boruta seed.
  boruta_fit <- Boruta::Boruta(as.formula(paste0(event_variable, " ~ .")), data = boruta_frame, doTrace = 0, maxRuns = 500) # Run all-relevant feature selection.
  boruta_fixed <- Boruta::TentativeRoughFix(boruta_fit) # Resolve tentative variables reproducibly.
  boruta_stats <- Boruta::attStats(boruta_fixed) # Extract importance and decisions.
  tibble::as_tibble(boruta_stats, rownames = "variable") # Return Boruta results.
} # End Boruta function.

boruta_results <- purrr::map_dfr(endpoints$event, ~ run_boruta(.x) %>% mutate(outcome = .x, .before = 1L)) # Run Boruta for all outcomes.
write.csv(boruta_results, file.path(output_dir, "boruta_feature_selection.csv"), row.names = FALSE) # Export Boruta results.

# 11. Consensus core features -------------------------------------------------- #
# The original code defining model-specific selection thresholds was not retained. #
# Therefore, the locked consensus sets below are taken from the archived selection #
# ledger, where every retained variable was selected by Boruta plus all seven      #
# machine-learning selection procedures.                                           #
core_features <- list( # Define the locked endpoint-specific consensus features.
  death_30d = c("age", "bilirubin", "charlson_cci", "hemoglobin", "s_EASIX", "temperature_max"), # Six 30-day features.
  death_60d = c("age", "aki", "bilirubin", "pt", "s_EASIX", "vasopressor", "weight"), # Seven 60-day features.
  death_90d = c("bun", "charlson_cci", "hemoglobin", "lactate", "pt", "ptt", "s_EASIX", "weight") # Eight 90-day features.
) # End consensus feature definitions.

core_feature_table <- purrr::imap_dfr(core_features, ~ tibble::tibble(outcome = .y, variable = .x)) # Convert feature sets to a table.
write.csv(core_feature_table, file.path(output_dir, "consensus_core_features.csv"), row.names = FALSE) # Export the feature table.

# 12. Seven prediction models and internal validation -------------------------- #
training_control <- caret::trainControl(method = "cv", number = 5L, classProbs = TRUE, summaryFunction = caret::twoClassSummary, savePredictions = "final", allowParallel = TRUE) # Define five-fold cross-validation.

fit_prediction_models <- function(event_variable) { # Train seven models for one endpoint.
  predictors <- core_features[[event_variable]] # Retrieve endpoint-specific core predictors.
  training_frame <- training_data[, c(event_variable, predictors)] # Create the modeling training frame.
  test_frame <- test_data[, c(event_variable, predictors)] # Create the modeling test frame.
  training_frame[[event_variable]] <- factor(training_frame[[event_variable]], levels = c(1, 0), labels = c("Died", "Survived")) # Make death the event class.
  test_frame[[event_variable]] <- factor(test_frame[[event_variable]], levels = c(1, 0), labels = c("Died", "Survived")) # Apply identical outcome coding.
  modeling_formula <- as.formula(paste0(event_variable, " ~ .")) # Build the model formula.
  set.seed(123) # Reset the tuning seed.
  fitted_models <- list( # Fit the seven prespecified prediction algorithms.
    "SVM" = caret::train(modeling_formula, data = training_frame, method = "svmRadial", metric = "ROC", trControl = training_control, preProcess = c("center", "scale"), tuneLength = 5L), # Support vector machine.
    "XGBoost" = caret::train(modeling_formula, data = training_frame, method = "xgbTree", metric = "ROC", trControl = training_control, tuneLength = 5L, verbose = FALSE), # Extreme gradient boosting.
    "Gradient Boosting" = caret::train(modeling_formula, data = training_frame, method = "gbm", metric = "ROC", trControl = training_control, tuneLength = 5L, verbose = FALSE), # Gradient boosting machine.
    "Random forest" = caret::train(modeling_formula, data = training_frame, method = "rf", metric = "ROC", trControl = training_control, tuneLength = 5L, importance = TRUE), # Random forest.
    "Logistic regression" = caret::train(modeling_formula, data = training_frame, method = "glm", family = binomial(), metric = "ROC", trControl = training_control), # Logistic regression.
    "LASSO" = caret::train(modeling_formula, data = training_frame, method = "glmnet", metric = "ROC", trControl = training_control, preProcess = c("center", "scale"), tuneGrid = expand.grid(alpha = 1, lambda = 10^seq(-4, 0, length.out = 50))), # LASSO regression.
    "Decision tree" = caret::train(modeling_formula, data = training_frame, method = "rpart", metric = "ROC", trControl = training_control, tuneLength = 10L) # Classification tree.
  ) # End model fitting.
  list(models = fitted_models, train = training_frame, test = test_frame) # Return fitted models and analysis frames.
} # End prediction-model function.

evaluate_prediction_model <- function(model, model_name, training_frame, test_frame, event_variable) { # Evaluate one fitted model.
  training_probability <- predict(model, newdata = training_frame, type = "prob")[, "Died"] # Obtain training-set death probabilities.
  test_probability <- predict(model, newdata = test_frame, type = "prob")[, "Died"] # Obtain test-set death probabilities.
  training_roc <- pROC::roc(training_frame[[event_variable]], training_probability, levels = c("Survived", "Died"), direction = "<", quiet = TRUE) # Build training ROC.
  optimal_threshold <- as.numeric(pROC::coords(training_roc, x = "best", best.method = "youden", ret = "threshold", transpose = FALSE)[1]) # Select training Youden threshold.
  test_roc <- pROC::roc(test_frame[[event_variable]], test_probability, levels = c("Survived", "Died"), direction = "<", quiet = TRUE) # Build test ROC.
  auc_ci <- pROC::ci.auc(test_roc, method = "delong") # Calculate DeLong AUC confidence interval.
  predicted_class <- factor(ifelse(test_probability >= optimal_threshold, "Died", "Survived"), levels = c("Died", "Survived")) # Apply the locked threshold.
  observed_class <- factor(test_frame[[event_variable]], levels = c("Died", "Survived")) # Standardize observed levels.
  confusion <- caret::confusionMatrix(predicted_class, observed_class, positive = "Died") # Calculate confusion-matrix metrics.
  sensitivity <- unname(confusion$byClass["Sensitivity"]) # Extract sensitivity.
  specificity <- unname(confusion$byClass["Specificity"]) # Extract specificity.
  ppv <- unname(confusion$byClass["Pos Pred Value"]) # Extract positive predictive value.
  npv <- unname(confusion$byClass["Neg Pred Value"]) # Extract negative predictive value.
  accuracy <- unname(confusion$overall["Accuracy"]) # Extract accuracy.
  f1_score <- if ((ppv + sensitivity) == 0) NA_real_ else 2 * ppv * sensitivity / (ppv + sensitivity) # Calculate F1 score.
  tibble::tibble(model = model_name, AUC = as.numeric(pROC::auc(test_roc)), CI_low = as.numeric(auc_ci[1]), CI_high = as.numeric(auc_ci[3]), Sensitivity = sensitivity, Specificity = specificity, Accuracy = accuracy, PPV = ppv, NPV = npv, F1_score = f1_score, Youden_J = sensitivity + specificity - 1, Cutoff = optimal_threshold) # Return performance metrics.
} # End model evaluator.

all_model_performance <- purrr::map_dfr(endpoints$event, function(event_variable) { # Iterate over mortality endpoints.
  fitted_set <- fit_prediction_models(event_variable) # Train all seven models.
  purrr::imap_dfr(fitted_set$models, ~ evaluate_prediction_model(.x, .y, fitted_set$train, fitted_set$test, event_variable)) %>% mutate(outcome = event_variable, .before = 1L) # Evaluate and label models.
}) # End endpoint iteration.

write.csv(all_model_performance, file.path(output_dir, "Supplementary_Table_S2_metrics.csv"), row.names = FALSE) # Export comprehensive test-set performance.

# 13. Prespecified subgroup analyses ------------------------------------------- #
subgroup_definitions <- list( # Define clinically prespecified subgroup variables.
  age_group = factor(ifelse(completed_data$age < 65, "<65", ">=65")), # Age cutoff at 65 years.
  gender = completed_data$gender, # Female versus male.
  race = completed_data$race, # White versus other race.
  pf_group = factor(ifelse(completed_data$pao2_fio2 < 100, "<100", ">=100")), # PaO2/FiO2 cutoff at 100 mmHg.
  pco2_group = factor(ifelse(completed_data$pco2 < 45, "<45", ">=45")), # PaCO2 cutoff at 45 mmHg.
  ventilation = factor(completed_data$ventilation), # Mechanical ventilation status.
  vasopressor = factor(completed_data$vasopressor), # Vasopressor status.
  aki = factor(completed_data$aki), # AKI status.
  crrt = factor(completed_data$crrt), # CRRT status.
  cci_group = factor(ifelse(completed_data$charlson_cci < 5, "<5", ">=5")) # CCI cutoff at 5.
) # End subgroup definitions.

completed_data <- bind_cols(completed_data, as_tibble(subgroup_definitions)[setdiff(names(subgroup_definitions), names(completed_data))]) # Append newly derived subgroup fields.

subgroup_map <- c(age_group = "age", gender = "gender", race = "race", pf_group = "pao2_fio2", pco2_group = "pco2", ventilation = "ventilation", vasopressor = "vasopressor", aki = "aki", crrt = "crrt", cci_group = "charlson_cci") # Map strata to adjustment variables.

run_subgroup <- function(endpoint_label, time_variable, event_variable, subgroup_variable) { # Estimate stratum-specific effects and interaction.
  adjustment_variables <- setdiff(model_3_covariates, subgroup_map[[subgroup_variable]]) # Exclude the stratifying covariate from adjustment.
  interaction_formula <- as.formula(paste0("Surv(", time_variable, ", ", event_variable, ") ~ quartile_num * ", subgroup_variable, " + ", paste(adjustment_variables, collapse = " + "))) # Build interaction formula.
  interaction_fit <- survival::coxph(interaction_formula, data = completed_data) # Fit the interaction model.
  interaction_terms <- grep(paste0("quartile_num:", subgroup_variable), names(coef(interaction_fit)), value = TRUE) # Locate interaction coefficients.
  interaction_p <- if (length(interaction_terms) == 1L) summary(interaction_fit)$coefficients[interaction_terms, "Pr(>|z|)"] else anova(update(interaction_fit, paste0(". ~ . - quartile_num:", subgroup_variable)), interaction_fit, test = "LRT")[2, "P(>|Chi|)"] # Obtain interaction P value.
  purrr::map_dfr(levels(factor(completed_data[[subgroup_variable]])), function(level_value) { # Iterate over strata.
    stratum_data <- completed_data[completed_data[[subgroup_variable]] == level_value, , drop = FALSE] # Subset one stratum.
    stratum_formula <- as.formula(paste0("Surv(", time_variable, ", ", event_variable, ") ~ quartile_num + ", paste(adjustment_variables, collapse = " + "))) # Build stratum-specific formula.
    stratum_fit <- survival::coxph(stratum_formula, data = stratum_data) # Fit the stratum-specific model.
    coefficient <- summary(stratum_fit)$coefficients["quartile_num", ] # Extract the ordinal s-EASIX coefficient.
    confidence <- summary(stratum_fit)$conf.int["quartile_num", ] # Extract its confidence interval.
    tibble::tibble(outcome = endpoint_label, subgroup = subgroup_variable, level = level_value, n = nrow(stratum_data), events = sum(stratum_data[[event_variable]]), HR = confidence["exp(coef)"], CI_low = confidence["lower .95"], CI_high = confidence["upper .95"], p_value = coefficient["Pr(>|z|)"], p_interaction = interaction_p) # Return stratum results.
  }) # End stratum iteration.
} # End subgroup function.

subgroup_results <- purrr::pmap_dfr(endpoints, function(label, time, event) purrr::map_dfr(names(subgroup_map), ~ run_subgroup(label, time, event, .x))) # Run every endpoint and subgroup.
write.csv(subgroup_results, file.path(output_dir, "subgroup_analysis_30_60_90.csv"), row.names = FALSE) # Export subgroup results.

# 14. Reproducibility record ---------------------------------------------------- #
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt")) # Save exact R and package versions.
saveRDS(list(imputation = imputed_object, cox_results = cox_results, rcs_results = rcs_results, model_performance = all_model_performance, subgroup_results = subgroup_results), file.path(output_dir, "analysis_objects.rds")) # Preserve key analysis objects.

# End of supplementary analysis code ------------------------------------------- #
