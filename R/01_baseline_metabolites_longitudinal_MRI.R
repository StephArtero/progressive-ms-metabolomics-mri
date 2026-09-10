############################################################
# Baseline Metabolites and Longitudinal MRI Outcomes
# Discovery analysis: SPRINT-MS
#
# This script evaluates associations between baseline plasma
# metabolites and longitudinal MRI outcomes using linear mixed-
# effects models.
#
# Primary model:
#   MRI outcome ~ time * baseline metabolite + covariates + (1 | PID)
#
# The metabolite-by-time interaction estimates whether baseline
# metabolite abundance is associated with longitudinal MRI change.
############################################################

#===========================================================
# 1. Load packages
#===========================================================

library(dplyr)
library(purrr)
library(lmerTest)
library(broom.mixed)
library(stringr)
library(tibble)
library(readr)

#===========================================================
# 2. User-defined input files
#===========================================================

# These are example relative paths for GitHub.
# Replace with the correct file names used in your project.
# Raw individual-level data are not provided in this repository.

metabolite_file <- "data/baseline_metabolites.csv"
mri_file        <- "data/longitudinal_mri.csv"
covariate_file  <- "data/baseline_covariates.csv"

# Output directory
output_dir <- "results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#===========================================================
# 3. Read input data
#===========================================================

# Expected structure:
# baseline_metabolites:
#   PID + one column per baseline metabolite
#
# longitudinal_mri:
#   PID, SCAN_WEEK, BPF, WMF, GMF, CTH
#
# baseline_covariates:
#   PID, SEX, RACE, BMI, TOBACCO, DMT, MSTYPE,
#   RANDASSIGNMENT, AGE_YRS_AT_COLLECTION,
#   DISEASE_DURATION_DIAGNOSIS

baseline_metabolites <- read_csv(metabolite_file, show_col_types = FALSE)
longitudinal_mri     <- read_csv(mri_file, show_col_types = FALSE)
baseline_covariates  <- read_csv(covariate_file, show_col_types = FALSE)

#===========================================================
# 4. Define analysis variables
#===========================================================

id_var <- "PID"

mri_outcomes <- c("BPF", "WMF", "GMF", "CTH")

continuous_covariates <- c(
  "AGE_YRS_AT_COLLECTION",
  "BMI",
  "DISEASE_DURATION_DIAGNOSIS"
)

categorical_covariates <- c(
  "SEX",
  "RACE",
  "TOBACCO",
  "DMT",
  "MSTYPE",
  "RANDASSIGNMENT"
)

metabolite_names <- setdiff(colnames(baseline_metabolites), id_var)

#===========================================================
# 5. Prepare longitudinal analysis dataset
#===========================================================

analysis_data <- longitudinal_mri %>%
  left_join(baseline_covariates, by = id_var) %>%
  left_join(baseline_metabolites, by = id_var) %>%
  mutate(
    PID = factor(PID),

    # Convert scheduled MRI week to continuous time in days.
    # If your dataset already contains MRI days, replace this line
    # with the corresponding variable.
    scan_week_numeric = suppressWarnings(as.numeric(as.character(SCAN_WEEK))),
    mri_days_raw = scan_week_numeric * 7,

    # Standardized time improves numerical stability.
    mri_days = as.numeric(scale(mri_days_raw)),

    # Categorical covariates
    SEX = factor(SEX),
    RACE = factor(RACE),
    TOBACCO = factor(TOBACCO),
    DMT = factor(DMT),
    MSTYPE = factor(MSTYPE),
    RANDASSIGNMENT = relevel(factor(RANDASSIGNMENT), ref = "Placebo"),

    # Standardized continuous covariates
    AGE_YRS_AT_COLLECTION = as.numeric(scale(AGE_YRS_AT_COLLECTION)),
    BMI = as.numeric(scale(BMI)),
    DISEASE_DURATION_DIAGNOSIS = as.numeric(scale(DISEASE_DURATION_DIAGNOSIS)),

    # Standardized MRI outcomes
    BPF = as.numeric(scale(BPF)),
    WMF = as.numeric(scale(WMF)),
    GMF = as.numeric(scale(GMF)),
    CTH = as.numeric(scale(CTH))
  ) %>%
  filter(!is.na(mri_days))

# Keep participants with at least two MRI observations
analysis_data <- analysis_data %>%
  group_by(PID) %>%
  filter(n_distinct(mri_days_raw[!is.na(mri_days_raw)]) >= 2) %>%
  ungroup()

#===========================================================
# 6. Helper function: extract interaction term
#===========================================================

extract_metabolite_time_interaction <- function(fit) {

  tidy_fit <- broom.mixed::tidy(
    fit,
    effects = "fixed",
    conf.int = TRUE,
    conf.level = 0.95,
    conf.method = "Wald"
  )

  interaction_row <- tidy_fit %>%
    filter(term %in% c(
      "mri_days:baseline_metabolite_z",
      "baseline_metabolite_z:mri_days"
    )) %>%
    slice(1)

  if (nrow(interaction_row) == 0) {
    interaction_row <- tibble(
      term = "mri_days:baseline_metabolite_z",
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_
    )
  }

  return(interaction_row)
}

#===========================================================
# 7. Primary model: random intercept
#===========================================================

run_primary_metabolite_model <- function(outcome_name, metabolite_name) {

  model_data <- analysis_data %>%
    mutate(
      baseline_metabolite_z = as.numeric(scale(.data[[metabolite_name]]))
    ) %>%
    select(
      PID,
      all_of(outcome_name),
      mri_days,
      baseline_metabolite_z,
      all_of(continuous_covariates),
      all_of(categorical_covariates)
    ) %>%
    na.omit()

  if (
    nrow(model_data) < 20 ||
    n_distinct(model_data$PID) < 10 ||
    isTRUE(sd(model_data$baseline_metabolite_z, na.rm = TRUE) == 0)
  ) {
    return(tibble(
      outcome = outcome_name,
      metabolite = metabolite_name,
      model_status = "skipped_insufficient_data"
    ))
  }

  model_formula <- as.formula(
    paste0(
      outcome_name,
      " ~ mri_days * baseline_metabolite_z + ",
      paste(c(continuous_covariates, categorical_covariates), collapse = " + "),
      " + (1 | PID)"
    )
  )

  fit <- tryCatch(
    lmerTest::lmer(
      model_formula,
      data = model_data,
      REML = TRUE,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    ),
    error = function(e) {
      message(
        "Model failed: outcome = ", outcome_name,
        "; metabolite = ", metabolite_name,
        "; error = ", e$message
      )
      return(NULL)
    }
  )

  if (is.null(fit)) {
    return(tibble(
      outcome = outcome_name,
      metabolite = metabolite_name,
      model_status = "model_failed"
    ))
  }

  interaction_result <- extract_metabolite_time_interaction(fit)

  tibble(
    outcome = outcome_name,
    metabolite = metabolite_name,
    model_status = "ok",
    n_obs = nrow(model_data),
    n_participants = n_distinct(model_data$PID),
    beta = interaction_result$estimate,
    se = interaction_result$std.error,
    lower_CI = interaction_result$conf.low,
    upper_CI = interaction_result$conf.high,
    p_value = interaction_result$p.value,
    AIC = AIC(fit),
    BIC = BIC(fit),
    singular_fit = isSingular(fit, tol = 1e-4)
  )
}

#===========================================================
# 8. Run primary analyses for all metabolites and MRI outcomes
#===========================================================

primary_results <- map_dfr(mri_outcomes, function(outcome_name) {

  message("Running outcome: ", outcome_name)

  map_dfr(metabolite_names, function(metabolite_name) {

    run_primary_metabolite_model(
      outcome_name = outcome_name,
      metabolite_name = metabolite_name
    )

  })
})

# FDR correction within each MRI outcome
primary_results <- primary_results %>%
  group_by(outcome) %>%
  mutate(
    q_value = p.adjust(p_value, method = "fdr")
  ) %>%
  ungroup() %>%
  arrange(outcome, p_value)

write_csv(
  primary_results,
  file.path(output_dir, "baseline_metabolites_longitudinal_MRI_primary_results.csv")
)

#===========================================================
# 9. Random-slope sensitivity analysis
#===========================================================

run_random_slope_sensitivity <- function(outcome_name, metabolite_name) {

  model_data <- analysis_data %>%
    mutate(
      baseline_metabolite_z = as.numeric(scale(.data[[metabolite_name]]))
    ) %>%
    select(
      PID,
      all_of(outcome_name),
      mri_days,
      baseline_metabolite_z,
      all_of(continuous_covariates),
      all_of(categorical_covariates)
    ) %>%
    na.omit()

  if (
    nrow(model_data) < 20 ||
    n_distinct(model_data$PID) < 10 ||
    isTRUE(sd(model_data$baseline_metabolite_z, na.rm = TRUE) == 0)
  ) {
    return(tibble(
      outcome = outcome_name,
      metabolite = metabolite_name,
      model_status = "skipped_insufficient_data"
    ))
  }

  fixed_part <- paste(
    paste0(
      outcome_name,
      " ~ mri_days * baseline_metabolite_z"
    ),
    paste(c(continuous_covariates, categorical_covariates), collapse = " + "),
    sep = " + "
  )

  formula_random_intercept <- as.formula(
    paste0(fixed_part, " + (1 | PID)")
  )

  # Uncorrelated random-slope model.
  # This evaluates participant-specific slopes for time without
  # estimating the intercept-slope covariance.
  formula_random_slope <- as.formula(
    paste0(fixed_part, " + (1 | PID) + (0 + mri_days | PID)")
  )

  fit_ri <- tryCatch(
    lmerTest::lmer(
      formula_random_intercept,
      data = model_data,
      REML = TRUE,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    ),
    error = function(e) NULL
  )

  fit_rs <- tryCatch(
    lmerTest::lmer(
      formula_random_slope,
      data = model_data,
      REML = TRUE,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit_ri) || is.null(fit_rs)) {
    return(tibble(
      outcome = outcome_name,
      metabolite = metabolite_name,
      model_status = "sensitivity_model_failed"
    ))
  }

  lrt <- anova(fit_ri, fit_rs, refit = FALSE)

  ri_interaction <- extract_metabolite_time_interaction(fit_ri)
  rs_interaction <- extract_metabolite_time_interaction(fit_rs)

  vc <- as.data.frame(VarCorr(fit_rs))

  random_slope_variance <- vc %>%
    filter(grp == "PID", var1 == "mri_days", is.na(var2)) %>%
    pull(vcov)

  if (length(random_slope_variance) == 0) {
    random_slope_variance <- NA_real_
  }

  tibble(
    outcome = outcome_name,
    metabolite = metabolite_name,
    model_status = "ok",
    n_obs = nrow(model_data),
    n_participants = n_distinct(model_data$PID),

    beta_random_intercept = ri_interaction$estimate,
    p_random_intercept = ri_interaction$p.value,

    beta_random_slope = rs_interaction$estimate,
    p_random_slope = rs_interaction$p.value,

    random_slope_variance = random_slope_variance,
    LRT_chisq_random_slope = lrt$Chisq[2],
    LRT_df_random_slope = lrt$`Chi Df`[2],
    LRT_p_random_slope = lrt$`Pr(>Chisq)`[2],

    AIC_random_intercept = AIC(fit_ri),
    AIC_random_slope = AIC(fit_rs),
    delta_AIC_slope_minus_intercept = AIC(fit_rs) - AIC(fit_ri),

    BIC_random_intercept = BIC(fit_ri),
    BIC_random_slope = BIC(fit_rs),
    delta_BIC_slope_minus_intercept = BIC(fit_rs) - BIC(fit_ri),

    singular_random_slope = isSingular(fit_rs, tol = 1e-4),

    direction_consistent = case_when(
      beta_random_intercept > 0 & beta_random_slope > 0 ~ "Yes",
      beta_random_intercept < 0 & beta_random_slope < 0 ~ "Yes",
      TRUE ~ "No"
    )
  )
}

# Run random-slope sensitivity analysis.
# This may take longer than the primary analysis.

random_slope_sensitivity_results <- map_dfr(mri_outcomes, function(outcome_name) {

  message("Running random-slope sensitivity for outcome: ", outcome_name)

  map_dfr(metabolite_names, function(metabolite_name) {

    run_random_slope_sensitivity(
      outcome_name = outcome_name,
      metabolite_name = metabolite_name
    )

  })
})

random_slope_sensitivity_results <- random_slope_sensitivity_results %>%
  group_by(outcome) %>%
  mutate(
    q_random_slope_LRT = p.adjust(LRT_p_random_slope, method = "fdr")
  ) %>%
  ungroup()

write_csv(
  random_slope_sensitivity_results,
  file.path(output_dir, "baseline_metabolites_longitudinal_MRI_random_slope_sensitivity.csv")
)

#===========================================================
# 10. Summary tables
#===========================================================

primary_summary <- primary_results %>%
  filter(model_status == "ok") %>%
  group_by(outcome) %>%
  summarise(
    n_metabolites_tested = n(),
    n_nominal_p_lt_0.05 = sum(p_value < 0.05, na.rm = TRUE),
    n_fdr_q_lt_0.05 = sum(q_value < 0.05, na.rm = TRUE),
    n_singular_primary = sum(singular_fit, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  primary_summary,
  file.path(output_dir, "baseline_metabolites_longitudinal_MRI_primary_summary.csv")
)

random_slope_summary <- random_slope_sensitivity_results %>%
  filter(model_status == "ok") %>%
  group_by(outcome) %>%
  summarise(
    n_metabolites_tested = n(),
    n_random_slope_LRT_p_lt_0.05 = sum(LRT_p_random_slope < 0.05, na.rm = TRUE),
    n_random_slope_LRT_q_lt_0.05 = sum(q_random_slope_LRT < 0.05, na.rm = TRUE),
    n_singular_random_slope = sum(singular_random_slope, na.rm = TRUE),
    n_direction_consistent = sum(direction_consistent == "Yes", na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  random_slope_summary,
  file.path(output_dir, "baseline_metabolites_longitudinal_MRI_random_slope_summary.csv")
)

#===========================================================
# 11. Optional diagnostic plots
#===========================================================

plot_lmm_diagnostics <- function(outcome_name, metabolite_name,
                                 model_type = c("random_intercept", "random_slope")) {

  model_type <- match.arg(model_type)

  model_data <- analysis_data %>%
    mutate(
      baseline_metabolite_z = as.numeric(scale(.data[[metabolite_name]]))
    ) %>%
    select(
      PID,
      all_of(outcome_name),
      mri_days,
      baseline_metabolite_z,
      all_of(continuous_covariates),
      all_of(categorical_covariates)
    ) %>%
    na.omit()

  fixed_part <- paste(
    paste0(
      outcome_name,
      " ~ mri_days * baseline_metabolite_z"
    ),
    paste(c(continuous_covariates, categorical_covariates), collapse = " + "),
    sep = " + "
  )

  random_part <- ifelse(
    model_type == "random_intercept",
    "(1 | PID)",
    "(1 | PID) + (0 + mri_days | PID)"
  )

  model_formula <- as.formula(
    paste0(fixed_part, " + ", random_part)
  )

  fit <- lmerTest::lmer(
    model_formula,
    data = model_data,
    REML = TRUE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )

  cat("\nOutcome:", outcome_name, "\n")
  cat("Metabolite:", metabolite_name, "\n")
  cat("Model type:", model_type, "\n")
  cat("AIC:", AIC(fit), "\n")
  cat("BIC:", BIC(fit), "\n")
  cat("Singular fit:", isSingular(fit, tol = 1e-4), "\n\n")

  if (model_type == "random_intercept") {
    par(mfrow = c(2, 2))
  } else {
    par(mfrow = c(2, 3))
  }

  hist(
    resid(fit),
    breaks = 30,
    main = "Residual distribution",
    xlab = "Residuals"
  )

  qqnorm(
    resid(fit),
    main = "Residual QQ plot"
  )
  qqline(resid(fit))

  plot(
    fitted(fit),
    resid(fit),
    xlab = "Fitted values",
    ylab = "Residuals",
    main = "Residuals vs fitted"
  )
  abline(h = 0, lty = 2)

  random_effects <- ranef(fit)$PID

  qqnorm(
    random_effects[, 1],
    main = "Random intercept QQ plot"
  )
  qqline(random_effects[, 1])

  if (model_type == "random_slope") {
    if ("mri_days" %in% colnames(random_effects)) {
      qqnorm(
        random_effects[, "mri_days"],
        main = "Random slope QQ plot"
      )
      qqline(random_effects[, "mri_days"])
    }
  }

  par(mfrow = c(1, 1))

  invisible(fit)
}

# Example diagnostic plot:
# Replace with a metabolite name from your results table.

# plot_lmm_diagnostics(
#   outcome_name = "BPF",
#   metabolite_name = metabolite_names[1],
#   model_type = "random_intercept"
# )

# plot_lmm_diagnostics(
#   outcome_name = "BPF",
#   metabolite_name = metabolite_names[1],
#   model_type = "random_slope"
# )

#===========================================================
# 12. Session information
#===========================================================

sessionInfo()
