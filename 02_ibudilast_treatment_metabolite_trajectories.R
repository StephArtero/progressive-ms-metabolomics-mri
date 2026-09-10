############################################################
# Ibudilast Treatment and Longitudinal Metabolomic Trajectories
# Discovery analysis: SPRINT-MS
#
# This script evaluates whether randomized ibudilast assignment is
# associated with longitudinal changes in metabolomic module eigengenes
# and individual metabolites.
#
# Primary model:
#   metabolomic outcome ~ time * treatment + covariates + (1 | PID)
#
# Sensitivity model:
#   metabolomic outcome ~ time * treatment + covariates + (1 + time | PID)
#
# The treatment-by-time interaction estimates whether metabolomic
# trajectories differ between ibudilast and placebo groups.

############################################################

#===========================================================
# 1. Load packages
#===========================================================

library(dplyr)
library(purrr)
library(readr)
library(lmerTest)
library(lme4)
library(broom.mixed)
library(stringr)
library(tibble)
library(WGCNA)

# Recommended for WGCNA in scripted analyses
options(stringsAsFactors = FALSE)
WGCNA::allowWGCNAThreads()

#===========================================================
# 2. User-defined input files
#===========================================================

# Expected input files are examples/placeholders. Replace with the
# appropriate relative file names used in your project.
#
# longitudinal_metabolomics_file should contain:
#   PID, sample_date, timepoint, and one column per metabolite.
#
# baseline_covariates_file should contain one row per participant:
#   PID, RANDASSIGNMENT, SEX, RACE, BMI, DMT,
#   DISEASE_DURATION_DIAGNOSIS, TOBACCO,
#   AGE_YRS_AT_COLLECTION, MSTYPE.
#
# module_membership_file should contain:
#   Metabolite, Module.
#
# metabolite_metadata_file is optional and may contain:
#   metabolite plus biochemical annotation columns.

longitudinal_metabolomics_file <- "data/sprintms_longitudinal_plasma_metabolites.csv"
baseline_covariates_file      <- "data/sprintms_baseline_covariates.csv"
module_membership_file        <- "data/metabolite_module_membership.csv"
metabolite_metadata_file      <- "data/plasma_metabolite_metadata.csv"

output_dir <- "results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#===========================================================
# 3. Read input data
#===========================================================

long_metabolomics <- read_csv(longitudinal_metabolomics_file, show_col_types = FALSE)
baseline_covariates <- read_csv(baseline_covariates_file, show_col_types = FALSE)
modules <- read_csv(module_membership_file, show_col_types = FALSE)

if (file.exists(metabolite_metadata_file)) {
  metabolite_metadata <- read_csv(metabolite_metadata_file, show_col_types = FALSE)
} else {
  metabolite_metadata <- NULL
}

#===========================================================
# 4. Define analysis variables
#===========================================================

id_var <- "PID"
timepoint_var <- "timepoint"
treatment_var <- "RANDASSIGNMENT"
time_var <- "days_since_baseline"

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
  "MSTYPE"
)

model_covariates <- c(continuous_covariates, categorical_covariates)

# Set this to the modules used for individual metabolite follow-up.
# Use NULL to analyze all metabolites in all non-grey modules.
selected_modules_for_metabolite_followup <- c("greenyellow", "salmon")

exclude_grey_module <- TRUE

#===========================================================
# 5. Helper functions
#===========================================================

scale_z <- function(x) {
  as.numeric(scale(as.numeric(x), center = TRUE, scale = TRUE))
}

median_impute_matrix <- function(x) {
  x <- as.data.frame(x)
  x <- x %>% mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
  for (j in seq_along(x)) {
    col_j <- x[[j]]
    if (all(is.na(col_j))) {
      x[[j]] <- 0
    } else {
      x[[j]][is.na(col_j)] <- median(col_j, na.rm = TRUE)
    }
  }
  x
}

extract_treatment_time_interaction <- function(fit) {
  tidy_fit <- broom.mixed::tidy(
    fit,
    effects = "fixed",
    conf.int = TRUE,
    conf.level = 0.95,
    conf.method = "Wald"
  )

  interaction_row <- tidy_fit %>%
    filter(str_detect(term, fixed(time_var)) &
             str_detect(term, fixed(treatment_var))) %>%
    slice(1)

  if (nrow(interaction_row) == 0) {
    return(tibble(
      term = paste0(time_var, ":", treatment_var),
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_
    ))
  }

  interaction_row
}

make_model_formula <- function(outcome_var,
                               random_effect = c("random_intercept",
                                                 "random_slope_correlated",
                                                 "random_slope_uncorrelated")) {
  random_effect <- match.arg(random_effect)

  fixed_part <- paste0(
    outcome_var,
    " ~ ", time_var, " * ", treatment_var,
    " + ", paste(model_covariates, collapse = " + ")
  )

  random_part <- switch(
    random_effect,
    random_intercept = paste0("(1 | ", id_var, ")"),
    random_slope_correlated = paste0("(1 + ", time_var, " | ", id_var, ")"),
    random_slope_uncorrelated = paste0("(1 | ", id_var, ") + (0 + ", time_var, " | ", id_var, ")")
  )

  as.formula(paste(fixed_part, random_part, sep = " + "))
}

#===========================================================
# 6. Prepare longitudinal metabolomics data
#===========================================================

# Calculate days since baseline from sample_date if not already available.
if (!(time_var %in% colnames(long_metabolomics))) {
  long_metabolomics <- long_metabolomics %>%
    mutate(sample_date = as.Date(sample_date)) %>%
    group_by(.data[[id_var]]) %>%
    mutate(
      baseline_date = min(sample_date, na.rm = TRUE),
      days_since_baseline = as.numeric(sample_date - baseline_date)
    ) %>%
    ungroup()
}

# Identify metabolite columns.
non_metabolite_cols <- c(
  id_var,
  "sample_date",
  "baseline_date",
  timepoint_var,
  time_var
)

metabolite_cols <- setdiff(colnames(long_metabolomics), non_metabolite_cols)

# Prepare module membership.
modules <- modules %>%
  rename(
    Metabolite = any_of(c("Metabolite", "metabolite", "CHEMICAL_NAME")),
    Module = any_of(c("Module", "module", "ModuleColor", "module_color"))
  ) %>%
  filter(!is.na(Metabolite), !is.na(Module))

if (exclude_grey_module) {
  modules <- modules %>% filter(Module != "grey")
}

common_metabolites <- intersect(metabolite_cols, modules$Metabolite)

if (length(common_metabolites) == 0) {
  stop("No overlapping metabolites were found between the longitudinal metabolomics data and module membership table.")
}

# Keep only needed variables and participants with treatment/covariates.
analysis_metabolomics <- long_metabolomics %>%
  select(all_of(c(id_var, timepoint_var, time_var, common_metabolites))) %>%
  left_join(baseline_covariates, by = id_var)

# Clean factors and standardize time/continuous covariates.
analysis_metabolomics <- analysis_metabolomics %>%
  mutate(
    PID = factor(.data[[id_var]]),
    RANDASSIGNMENT = relevel(factor(.data[[treatment_var]]), ref = "Placebo"),
    SEX = factor(SEX),
    RACE = factor(RACE),
    TOBACCO = factor(TOBACCO),
    DMT = factor(DMT),
    MSTYPE = factor(MSTYPE),
    days_since_baseline = scale_z(.data[[time_var]]),
    AGE_YRS_AT_COLLECTION = scale_z(AGE_YRS_AT_COLLECTION),
    BMI = scale_z(BMI),
    DISEASE_DURATION_DIAGNOSIS = scale_z(DISEASE_DURATION_DIAGNOSIS)
  )

#===========================================================
# 7. Compute module eigengenes at each metabolomics timepoint
#===========================================================

compute_module_eigengenes_by_timepoint <- function(dat, metabolite_names, module_table) {
  module_table <- module_table %>%
    filter(Metabolite %in% metabolite_names)

  metabolite_names <- metabolite_names[metabolite_names %in% module_table$Metabolite]
  module_colors <- module_table$Module[match(metabolite_names, module_table$Metabolite)]

  dat %>%
    group_split(.data[[timepoint_var]]) %>%
    map_dfr(function(time_df) {
      current_timepoint <- unique(time_df[[timepoint_var]])[1]

      x <- time_df[, metabolite_names, drop = FALSE]
      x <- median_impute_matrix(x)

      module_eigengenes <- WGCNA::moduleEigengenes(
        expr = as.data.frame(x),
        colors = module_colors
      )$eigengenes

      colnames(module_eigengenes) <- sub("^ME", "", colnames(module_eigengenes))

      module_eigengenes <- as_tibble(module_eigengenes) %>%
        mutate(
          PID = time_df[[id_var]],
          timepoint = current_timepoint,
          days_since_baseline = time_df[[time_var]]
        ) %>%
        relocate(PID, timepoint, days_since_baseline)

      module_eigengenes
    })
}

module_long <- compute_module_eigengenes_by_timepoint(
  dat = analysis_metabolomics,
  metabolite_names = common_metabolites,
  module_table = modules
)

# Add treatment and baseline covariates to module eigengene dataset.
module_long <- module_long %>%
  left_join(baseline_covariates, by = id_var) %>%
  mutate(
    PID = factor(PID),
    RANDASSIGNMENT = relevel(factor(RANDASSIGNMENT), ref = "Placebo"),
    SEX = factor(SEX),
    RACE = factor(RACE),
    TOBACCO = factor(TOBACCO),
    DMT = factor(DMT),
    MSTYPE = factor(MSTYPE),
    days_since_baseline = scale_z(days_since_baseline),
    AGE_YRS_AT_COLLECTION = scale_z(AGE_YRS_AT_COLLECTION),
    BMI = scale_z(BMI),
    DISEASE_DURATION_DIAGNOSIS = scale_z(DISEASE_DURATION_DIAGNOSIS)
  )

module_names <- setdiff(
  colnames(module_long),
  c(id_var, timepoint_var, time_var, treatment_var, model_covariates)
)

#===========================================================
# 8. Module-level treatment trajectory analysis
#===========================================================

run_module_treatment_model <- function(module_name,
                                       random_effect = c("random_intercept",
                                                         "random_slope_correlated",
                                                         "random_slope_uncorrelated"),
                                       REML_value = TRUE) {
  random_effect <- match.arg(random_effect)

  dat <- module_long %>%
    mutate(module_z = scale_z(.data[[module_name]])) %>%
    select(
      PID,
      module_z,
      all_of(c(time_var, treatment_var, model_covariates))
    ) %>%
    na.omit()

  if (nrow(dat) < 20 || n_distinct(dat$PID) < 10 ||
      isTRUE(sd(dat$module_z, na.rm = TRUE) == 0)) {
    return(tibble(
      module = module_name,
      model_status = "skipped_insufficient_data"
    ))
  }

  fml <- make_model_formula("module_z", random_effect = random_effect)

  fit <- tryCatch(
    lmerTest::lmer(
      fml,
      data = dat,
      REML = REML_value,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    ),
    error = function(e) {
      message("Module model failed for ", module_name, ": ", e$message)
      return(NULL)
    }
  )

  if (is.null(fit)) {
    return(tibble(module = module_name, model_status = "model_failed"))
  }

  interaction <- extract_treatment_time_interaction(fit)

  tibble(
    module = module_name,
    model_status = "ok",
    random_effect = random_effect,
    n_obs = nrow(dat),
    n_participants = n_distinct(dat$PID),
    beta = interaction$estimate,
    se = interaction$std.error,
    lower_CI = interaction$conf.low,
    upper_CI = interaction$conf.high,
    p_value = interaction$p.value,
    AIC = AIC(fit),
    BIC = BIC(fit),
    singular_fit = lme4::isSingular(fit, tol = 1e-4)
  )
}

module_treatment_results <- map_dfr(module_names, function(module_name) {
  run_module_treatment_model(
    module_name = module_name,
    random_effect = "random_intercept",
    REML_value = TRUE
  )
}) %>%
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

write_csv(
  module_treatment_results,
  file.path(output_dir, "ibudilast_module_trajectory_results.csv")
)

# Formatted version for manuscript tables.
module_treatment_results_formatted <- module_treatment_results %>%
  mutate(
    beta_CI = paste0(
      sprintf("%.3f", beta),
      " (", sprintf("%.3f", lower_CI),
      ", ", sprintf("%.3f", upper_CI), ")"
    ),
    p_value_formatted = formatC(p_value, format = "e", digits = 2),
    q_value_formatted = formatC(q_value, format = "e", digits = 2)
  )

write_csv(
  module_treatment_results_formatted,
  file.path(output_dir, "ibudilast_module_trajectory_results_formatted.csv")
)

#===========================================================
# 9. Random-effects comparison for module-level treatment models
#===========================================================

compare_module_random_effects <- function(module_name) {
  dat <- module_long %>%
    mutate(module_z = scale_z(.data[[module_name]])) %>%
    select(
      PID,
      module_z,
      all_of(c(time_var, treatment_var, model_covariates))
    ) %>%
    na.omit()

  if (nrow(dat) < 20 || n_distinct(dat$PID) < 10 ||
      isTRUE(sd(dat$module_z, na.rm = TRUE) == 0)) {
    return(tibble(
      module = module_name,
      model_status = "skipped_insufficient_data"
    ))
  }

  f_intercept <- make_model_formula("module_z", random_effect = "random_intercept")
  f_slope <- make_model_formula("module_z", random_effect = "random_slope_correlated")

  fit_intercept <- tryCatch(
    lmerTest::lmer(
      f_intercept,
      data = dat,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    ),
    error = function(e) NULL
  )

  fit_slope <- tryCatch(
    lmerTest::lmer(
      f_slope,
      data = dat,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    ),
    error = function(e) NULL
  )

  if (is.null(fit_intercept) || is.null(fit_slope)) {
    return(tibble(module = module_name, model_status = "model_failed"))
  }

  lrt <- anova(fit_intercept, fit_slope)
  interaction_intercept <- extract_treatment_time_interaction(fit_intercept)
  interaction_slope <- extract_treatment_time_interaction(fit_slope)

  vc <- as.data.frame(VarCorr(fit_slope))
  slope_variance <- vc %>%
    filter(grp == id_var, var1 == time_var, is.na(var2)) %>%
    pull(vcov)

  if (length(slope_variance) == 0) {
    slope_variance <- NA_real_
  }

  tibble(
    module = module_name,
    model_status = "ok",
    n_obs = nrow(dat),
    n_participants = n_distinct(dat$PID),

    beta_random_intercept = interaction_intercept$estimate,
    se_random_intercept = interaction_intercept$std.error,
    p_random_intercept = interaction_intercept$p.value,

    beta_random_slope = interaction_slope$estimate,
    se_random_slope = interaction_slope$std.error,
    p_random_slope = interaction_slope$p.value,

    AIC_random_intercept = AIC(fit_intercept),
    AIC_random_slope = AIC(fit_slope),
    delta_AIC_slope_minus_intercept = AIC(fit_slope) - AIC(fit_intercept),

    BIC_random_intercept = BIC(fit_intercept),
    BIC_random_slope = BIC(fit_slope),
    delta_BIC_slope_minus_intercept = BIC(fit_slope) - BIC(fit_intercept),

    LRT_chisq_random_slope = lrt$Chisq[2],
    LRT_df_random_slope = lrt$`Chi Df`[2],
    LRT_p_random_slope = lrt$`Pr(>Chisq)`[2],

    random_slope_variance = slope_variance,
    singular_random_slope = lme4::isSingular(fit_slope, tol = 1e-4),

    direction_consistent = case_when(
      beta_random_intercept > 0 & beta_random_slope > 0 ~ "Yes",
      beta_random_intercept < 0 & beta_random_slope < 0 ~ "Yes",
      TRUE ~ "No"
    ),

    AIC_favors = case_when(
      delta_AIC_slope_minus_intercept < 0 ~ "Random slope",
      delta_AIC_slope_minus_intercept > 0 ~ "Random intercept",
      TRUE ~ "Equivalent"
    ),

    BIC_favors = case_when(
      delta_BIC_slope_minus_intercept < 0 ~ "Random slope",
      delta_BIC_slope_minus_intercept > 0 ~ "Random intercept",
      TRUE ~ "Equivalent"
    )
  )
}

module_random_effects_comparison <- map_dfr(
  module_names,
  compare_module_random_effects
) %>%
  mutate(
    q_random_intercept = p.adjust(p_random_intercept, method = "BH"),
    q_random_slope = p.adjust(p_random_slope, method = "BH"),
    q_random_slope_LRT = p.adjust(LRT_p_random_slope, method = "BH")
  ) %>%
  arrange(q_random_slope)

write_csv(
  module_random_effects_comparison,
  file.path(output_dir, "ibudilast_module_random_effects_comparison.csv")
)

#===========================================================
# 10. Individual metabolite treatment trajectory analysis
#===========================================================

if (is.null(selected_modules_for_metabolite_followup)) {
  selected_metabolites <- modules %>%
    distinct(Metabolite, Module)
} else {
  selected_metabolites <- modules %>%
    filter(Module %in% selected_modules_for_metabolite_followup) %>%
    distinct(Metabolite, Module)
}

selected_metabolites <- selected_metabolites %>%
  filter(Metabolite %in% common_metabolites)

run_metabolite_treatment_model <- function(metabolite_name,
                                           module_name,
                                           random_effect = c("random_intercept",
                                                             "random_slope_correlated",
                                                             "random_slope_uncorrelated"),
                                           REML_value = TRUE) {
  random_effect <- match.arg(random_effect)

  dat <- analysis_metabolomics %>%
    mutate(metabolite_z = scale_z(.data[[metabolite_name]])) %>%
    select(
      PID,
      metabolite_z,
      all_of(c(time_var, treatment_var, model_covariates))
    ) %>%
    na.omit()

  if (nrow(dat) < 20 || n_distinct(dat$PID) < 10 ||
      isTRUE(sd(dat$metabolite_z, na.rm = TRUE) == 0)) {
    return(tibble(
      metabolite = metabolite_name,
      module = module_name,
      model_status = "skipped_insufficient_data"
    ))
  }

  fml <- make_model_formula("metabolite_z", random_effect = random_effect)

  fit <- tryCatch(
    lmerTest::lmer(
      fml,
      data = dat,
      REML = REML_value,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    ),
    error = function(e) {
      message("Metabolite model failed for ", metabolite_name, ": ", e$message)
      return(NULL)
    }
  )

  if (is.null(fit)) {
    return(tibble(
      metabolite = metabolite_name,
      module = module_name,
      model_status = "model_failed"
    ))
  }

  interaction <- extract_treatment_time_interaction(fit)

  tibble(
    metabolite = metabolite_name,
    module = module_name,
    model_status = "ok",
    random_effect = random_effect,
    n_obs = nrow(dat),
    n_participants = n_distinct(dat$PID),
    beta = interaction$estimate,
    se = interaction$std.error,
    lower_CI = interaction$conf.low,
    upper_CI = interaction$conf.high,
    p_value = interaction$p.value,
    AIC = AIC(fit),
    BIC = BIC(fit),
    singular_fit = lme4::isSingular(fit, tol = 1e-4)
  )
}

metabolite_treatment_results <- pmap_dfr(
  list(selected_metabolites$Metabolite, selected_metabolites$Module),
  function(Metabolite, Module) {
    run_metabolite_treatment_model(
      metabolite_name = Metabolite,
      module_name = Module,
      random_effect = "random_intercept",
      REML_value = TRUE
    )
  }
) %>%
  group_by(module) %>%
  mutate(q_value_within_module = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(q_value_all_selected = p.adjust(p_value, method = "BH")) %>%
  arrange(module, p_value)

# Add metabolite metadata if available.
if (!is.null(metabolite_metadata) && "metabolite" %in% colnames(metabolite_metadata)) {
  metabolite_treatment_results <- metabolite_treatment_results %>%
    left_join(metabolite_metadata, by = "metabolite")
}

write_csv(
  metabolite_treatment_results,
  file.path(output_dir, "ibudilast_individual_metabolite_trajectory_results.csv")
)

metabolite_treatment_results_formatted <- metabolite_treatment_results %>%
  mutate(
    beta_CI = paste0(
      sprintf("%.3f", beta),
      " (", sprintf("%.3f", lower_CI),
      ", ", sprintf("%.3f", upper_CI), ")"
    ),
    p_value_formatted = formatC(p_value, format = "e", digits = 2),
    q_value_within_module_formatted = formatC(q_value_within_module, format = "e", digits = 2),
    q_value_all_selected_formatted = formatC(q_value_all_selected, format = "e", digits = 2)
  )

write_csv(
  metabolite_treatment_results_formatted,
  file.path(output_dir, "ibudilast_individual_metabolite_trajectory_results_formatted.csv")
)

#===========================================================
# 11. Summary tables
#===========================================================

module_summary <- module_treatment_results %>%
  filter(model_status == "ok") %>%
  summarise(
    n_modules_tested = n(),
    n_nominal_p_lt_0.05 = sum(p_value < 0.05, na.rm = TRUE),
    n_fdr_q_lt_0.05 = sum(q_value < 0.05, na.rm = TRUE),
    n_singular = sum(singular_fit, na.rm = TRUE)
  )

metabolite_summary <- metabolite_treatment_results %>%
  filter(model_status == "ok") %>%
  group_by(module) %>%
  summarise(
    n_metabolites_tested = n(),
    n_nominal_p_lt_0.05 = sum(p_value < 0.05, na.rm = TRUE),
    n_within_module_q_lt_0.05 = sum(q_value_within_module < 0.05, na.rm = TRUE),
    n_singular = sum(singular_fit, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(module_summary, file.path(output_dir, "ibudilast_module_trajectory_summary.csv"))
write_csv(metabolite_summary, file.path(output_dir, "ibudilast_metabolite_trajectory_summary.csv"))

#===========================================================
# 12. Optional diagnostic plots shown in R, not saved
#===========================================================

plot_treatment_model_diagnostics <- function(outcome_type = c("module", "metabolite"),
                                             outcome_name,
                                             random_effect = c("random_intercept",
                                                               "random_slope_correlated",
                                                               "random_slope_uncorrelated")) {
  outcome_type <- match.arg(outcome_type)
  random_effect <- match.arg(random_effect)

  if (outcome_type == "module") {
    dat <- module_long %>%
      mutate(outcome_z = scale_z(.data[[outcome_name]])) %>%
      select(PID, outcome_z, all_of(c(time_var, treatment_var, model_covariates))) %>%
      na.omit()
  } else {
    dat <- analysis_metabolomics %>%
      mutate(outcome_z = scale_z(.data[[outcome_name]])) %>%
      select(PID, outcome_z, all_of(c(time_var, treatment_var, model_covariates))) %>%
      na.omit()
  }

  fml <- make_model_formula("outcome_z", random_effect = random_effect)

  fit <- lmerTest::lmer(
    fml,
    data = dat,
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )

  cat("\nOutcome type:", outcome_type, "\n")
  cat("Outcome name:", outcome_name, "\n")
  cat("Random-effect structure:", random_effect, "\n")
  cat("AIC:", AIC(fit), "\n")
  cat("BIC:", BIC(fit), "\n")
  cat("Singular fit:", isSingular(fit, tol = 1e-4), "\n\n")

  if (random_effect == "random_intercept") {
    par(mfrow = c(2, 2))
  } else {
    par(mfrow = c(2, 3))
  }

  hist(resid(fit), breaks = 30, main = "Residual distribution", xlab = "Residuals")

  qqnorm(resid(fit), main = "Residual QQ plot")
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

  qqnorm(random_effects[, 1], main = "Random intercept QQ plot")
  qqline(random_effects[, 1])

  if (random_effect != "random_intercept") {
    slope_col <- colnames(random_effects)[str_detect(colnames(random_effects), fixed(time_var))]
    if (length(slope_col) > 0) {
      qqnorm(random_effects[, slope_col[1]], main = "Random slope QQ plot")
      qqline(random_effects[, slope_col[1]])
    }
  }

  par(mfrow = c(1, 1))
  invisible(fit)
}

# Example diagnostic plots:
# plot_treatment_model_diagnostics(
#   outcome_type = "module",
#   outcome_name = module_names[1],
#   random_effect = "random_intercept"
# )
#
# plot_treatment_model_diagnostics(
#   outcome_type = "metabolite",
#   outcome_name = selected_metabolites$Metabolite[1],
#   random_effect = "random_intercept"
# )


