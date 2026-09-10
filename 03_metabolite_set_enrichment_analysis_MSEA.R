############################################################
# Metabolite Set Enrichment Analysis (MSEA)
# Discovery analysis: SPRINT-MS
#
#
# Primary mixed-effects model:
#   MRI outcome ~ time * baseline metabolite + covariates + (1 | PID)
#
# MSEA input:
#   ranked metabolite-level statistics from the primary model
############################################################

#===========================================================
# 1. Load packages
#===========================================================

library(dplyr)
library(purrr)
library(readr)
library(tibble)
library(stringr)
library(forcats)
library(lmerTest)
library(broom.mixed)
library(fgsea)
library(ggplot2)

#===========================================================
# 2. User-defined input and output files
#===========================================================

# Expected input files:
# 1) baseline_metabolites.csv:
#      PID + one column per baseline metabolite
# 2) longitudinal_mri.csv:
#      PID, SCAN_WEEK, BPF, WMF, GMF, CTH
# 3) baseline_covariates.csv:
#      PID, SEX, RACE, BMI, TOBACCO, DMT, MSTYPE,
#      RANDASSIGNMENT, AGE_YRS_AT_COLLECTION,
#      DISEASE_DURATION_DIAGNOSIS
# 4) metabolite_pathway_annotation.csv:
#      Metabolite, PATHWAY

metabolite_file <- "data/baseline_metabolites.csv"
mri_file        <- "data/longitudinal_mri.csv"
covariate_file  <- "data/baseline_covariates.csv"
pathway_file    <- "data/metabolite_pathway_annotation.csv"

output_dir <- "results"
plot_dir   <- "results/figures"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#===========================================================
# 3. Read input data
#===========================================================

baseline_metabolites <- read_csv(metabolite_file, show_col_types = FALSE)
longitudinal_mri     <- read_csv(mri_file, show_col_types = FALSE)
baseline_covariates  <- read_csv(covariate_file, show_col_types = FALSE)
pathway_annotation   <- read_csv(pathway_file, show_col_types = FALSE)

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

# Optional: keep only metabolites that have pathway annotation.
metabolite_names <- intersect(
  metabolite_names,
  unique(pathway_annotation$Metabolite)
)

#===========================================================
# 5. Prepare longitudinal MRI analysis dataset
#===========================================================

analysis_data <- longitudinal_mri %>%
  left_join(baseline_covariates, by = id_var) %>%
  left_join(baseline_metabolites, by = id_var) %>%
  mutate(
    PID = factor(PID),

    # Convert scheduled MRI week to continuous time in days.
    # If the dataset already contains MRI days, replace this line
    # with the corresponding variable.
    scan_week_numeric = suppressWarnings(as.numeric(as.character(SCAN_WEEK))),
    mri_days_raw = scan_week_numeric * 7,

    # Standardized time improves numerical stability and makes
    # coefficients comparable across models.
    mri_days = as.numeric(scale(mri_days_raw)),

    # Categorical covariates.
    SEX = factor(SEX),
    RACE = factor(RACE),
    TOBACCO = factor(TOBACCO),
    DMT = factor(DMT),
    MSTYPE = factor(MSTYPE),
    RANDASSIGNMENT = relevel(factor(RANDASSIGNMENT), ref = "Placebo"),

    # Continuous covariates.
    AGE_YRS_AT_COLLECTION = as.numeric(scale(AGE_YRS_AT_COLLECTION)),
    BMI = as.numeric(scale(BMI)),
    DISEASE_DURATION_DIAGNOSIS = as.numeric(scale(DISEASE_DURATION_DIAGNOSIS)),

    # MRI outcomes.
    BPF = as.numeric(scale(BPF)),
    WMF = as.numeric(scale(WMF)),
    GMF = as.numeric(scale(GMF)),
    CTH = as.numeric(scale(CTH))
  ) %>%
  filter(!is.na(mri_days))

# Keep participants with at least two non-missing MRI observations.
analysis_data <- analysis_data %>%
  group_by(PID) %>%
  filter(n_distinct(mri_days_raw[!is.na(mri_days_raw)]) >= 2) %>%
  ungroup()

#===========================================================
# 6. Helper function: extract metabolite-by-time interaction
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
# 7. Run metabolite-level mixed-effects models
#===========================================================

run_one_metabolite_mri_model <- function(outcome_name, metabolite_name) {

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
    all(!is.finite(model_data$baseline_metabolite_z)) ||
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
    statistic = interaction_result$statistic,
    lower_CI = interaction_result$conf.low,
    upper_CI = interaction_result$conf.high,
    p_value = interaction_result$p.value,
    AIC = AIC(fit),
    BIC = BIC(fit),
    singular_fit = isSingular(fit, tol = 1e-4)
  )
}

metabolite_mri_results <- map_dfr(mri_outcomes, function(outcome_name) {

  message("Running metabolite-level MRI models for outcome: ", outcome_name)

  map_dfr(metabolite_names, function(metabolite_name) {
    run_one_metabolite_mri_model(
      outcome_name = outcome_name,
      metabolite_name = metabolite_name
    )
  })
})

metabolite_mri_results <- metabolite_mri_results %>%
  group_by(outcome) %>%
  mutate(q_value = p.adjust(p_value, method = "fdr")) %>%
  ungroup() %>%
  arrange(outcome, p_value)

write_csv(
  metabolite_mri_results,
  file.path(output_dir, "metabolite_level_longitudinal_MRI_results_for_MSEA.csv")
)

#===========================================================
# 8. Build pathway sets for MSEA
#===========================================================

pathway_df <- pathway_annotation %>%
  transmute(
    Metabolite = as.character(Metabolite),
    PATHWAY = as.character(PATHWAY)
  ) %>%
  filter(
    !is.na(Metabolite), Metabolite != "",
    !is.na(PATHWAY), PATHWAY != ""
  ) %>%
  distinct()

pathway_list <- split(pathway_df$Metabolite, pathway_df$PATHWAY)

stopifnot(
  is.list(pathway_list),
  length(pathway_list) > 0,
  all(map_lgl(pathway_list, ~ is.character(.x)))
)

#===========================================================
# 9. Run fgsea for one MRI outcome
#===========================================================

run_msea_one_outcome <- function(outcome_name,
                                 rank_metric = c("statistic", "beta", "signed_logp"),
                                 min_size = 5,
                                 max_size = 300) {

  rank_metric <- match.arg(rank_metric)

  outcome_results <- metabolite_mri_results %>%
    filter(outcome == outcome_name, model_status == "ok") %>%
    mutate(
      beta = suppressWarnings(as.numeric(beta)),
      se = suppressWarnings(as.numeric(se)),
      statistic = suppressWarnings(as.numeric(statistic)),
      p_value = suppressWarnings(as.numeric(p_value)),
      rank_score = case_when(
        rank_metric == "statistic" & is.finite(statistic) ~ statistic,
        rank_metric == "beta" & is.finite(beta) ~ beta,
        rank_metric == "signed_logp" & is.finite(p_value) & is.finite(beta) ~
          -log10(pmax(p_value, .Machine$double.xmin)) * sign(beta),
        TRUE ~ NA_real_
      )
    ) %>%
    filter(
      !is.na(metabolite), metabolite != "",
      is.finite(rank_score)
    ) %>%
    group_by(metabolite) %>%
    slice_max(order_by = abs(rank_score), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(rank_score))

  ranks <- outcome_results$rank_score
  names(ranks) <- outcome_results$metabolite

  stopifnot(
    is.numeric(ranks),
    !is.null(names(ranks)),
    all(names(ranks) != ""),
    all(!is.na(names(ranks)))
  )

  universe <- intersect(names(ranks), unique(unlist(pathway_list)))

  if (length(universe) < 10) {
    stop(
      "Very small overlap between ranked metabolites and pathway annotation for ",
      outcome_name,
      ". Check metabolite names in the model results and pathway annotation."
    )
  }

  ranks_use <- ranks[universe]
  ranks_use <- sort(ranks_use, decreasing = TRUE)

  pathway_list_use <- lapply(pathway_list, function(x) intersect(x, universe))
  pathway_list_use <- pathway_list_use[
    lengths(pathway_list_use) >= min_size & lengths(pathway_list_use) <= max_size
  ]

  if (length(pathway_list_use) == 0) {
    stop("No pathways meet min_size/max_size after metabolite matching for ", outcome_name)
  }

  fgsea_result <- fgsea::fgsea(
    pathways = pathway_list_use,
    stats = ranks_use,
    minSize = min_size,
    maxSize = max_size
  ) %>%
    arrange(padj) %>%
    mutate(
      outcome = outcome_name,
      rank_metric = rank_metric,
      ES = as.numeric(ES),
      NES = as.numeric(NES),
      pval = as.numeric(pval),
      padj = as.numeric(padj),
      size = as.integer(size)
    )

  return(fgsea_result)
}

#===========================================================
# 10. Run MSEA for all MRI outcomes
#===========================================================

msea_results <- map_dfr(mri_outcomes, function(outcome_name) {
  message("Running MSEA for outcome: ", outcome_name)
  run_msea_one_outcome(
    outcome_name = outcome_name,
    rank_metric = "statistic",
    min_size = 5,
    max_size = 300
  )
})

# Collapse leading-edge list column before writing to CSV.
msea_results_flat <- msea_results %>%
  mutate(
    leadingEdge = map_chr(leadingEdge, function(x) {
      if (length(x) == 0 || all(is.na(x))) return(NA_character_)
      paste(as.character(x), collapse = ";")
    })
  )

write_csv(
  msea_results_flat,
  file.path(output_dir, "MSEA_fgsea_results_all_MRI_outcomes.csv")
)

# Also write separate files by MRI outcome.
walk(mri_outcomes, function(outcome_name) {

  msea_results_flat %>%
    filter(outcome == outcome_name) %>%
    write_csv(
      file.path(output_dir, paste0("MSEA_fgsea_results_", outcome_name, ".csv"))
    )
})

#===========================================================
# 11. Summary table
#===========================================================

msea_summary <- msea_results_flat %>%
  group_by(outcome) %>%
  summarise(
    n_pathways_tested = n(),
    n_pathways_p_lt_0.05 = sum(pval < 0.05, na.rm = TRUE),
    n_pathways_fdr_lt_0.05 = sum(padj < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  msea_summary,
  file.path(output_dir, "MSEA_summary_by_MRI_outcome.csv")
)

#===========================================================
# 12. Dot plot of significant MSEA pathways
#===========================================================

plot_msea_dotplot <- function(msea_df,
                              padj_threshold = 0.05,
                              top_n_if_none = 10,
                              output_file = file.path(plot_dir, "MSEA_dotplot_all_MRI_outcomes.png")) {

  plot_df <- msea_df %>%
    mutate(
      outcome = factor(outcome, levels = mri_outcomes),
      pathway_label = str_wrap(pathway, width = 35),
      padj = as.numeric(padj),
      pval = as.numeric(pval),
      NES = as.numeric(NES),
      ES = as.numeric(ES)
    )

  significant_df <- plot_df %>%
    filter(!is.na(padj), padj < padj_threshold)

  # If no pathways pass FDR, plot top pathways by p-value for visualization only.
  if (nrow(significant_df) == 0) {
    message(
      "No pathways passed FDR < ", padj_threshold,
      ". Plotting top ", top_n_if_none,
      " pathways per outcome by nominal p-value for visualization."
    )

    significant_df <- plot_df %>%
      group_by(outcome) %>%
      slice_min(order_by = pval, n = top_n_if_none, with_ties = FALSE) %>%
      ungroup()
  }

  # Cap very small adjusted p-values for bubble-size stability.
  cap_value <- quantile(significant_df$padj, probs = 0.05, na.rm = TRUE)
  cap_value <- max(cap_value, .Machine$double.xmin, na.rm = TRUE)

  significant_df <- significant_df %>%
    mutate(
      padj_capped = pmax(padj, cap_value),
      neg_log10_padj = -log10(padj_capped)
    )

  pathway_order <- significant_df %>%
    group_by(pathway_label) %>%
    summarise(median_NES = median(NES, na.rm = TRUE), .groups = "drop") %>%
    arrange(median_NES) %>%
    pull(pathway_label)

  significant_df <- significant_df %>%
    mutate(pathway_label = factor(pathway_label, levels = rev(pathway_order)))

  xmax <- max(abs(significant_df$NES), na.rm = TRUE)
  xlim_use <- c(-xmax - 0.2, xmax + 0.2)

  p <- ggplot(significant_df, aes(x = NES, y = pathway_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_point(aes(size = neg_log10_padj, color = ES), alpha = 0.95) +
    facet_grid(. ~ outcome, scales = "free_x", space = "free_x") +
    scale_x_continuous(limits = xlim_use, expand = expansion(mult = 0.07)) +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    scale_size_continuous(
      name = expression(-log[10](FDR)),
      range = c(2.0, 8.0)
    ) +
    scale_color_gradient2(
      name = "ES",
      low = "#3B4CC0",
      mid = "white",
      high = "#FCA636",
      midpoint = 0
    ) +
    labs(
      title = "Metabolite Set Enrichment Analysis (MSEA)",
      x = "Normalized enrichment score (NES)",
      y = "Pathway"
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      axis.title.y = element_text(face = "bold"),
      axis.text.y = element_text(size = 8),
      panel.spacing.x = unit(8, "pt"),
      plot.margin = margin(10, 20, 10, 10)
    )

  print(p)

  ggsave(
    filename = output_file,
    plot = p,
    width = 12,
    height = 8,
    dpi = 300
  )

  return(p)
}

msea_dotplot <- plot_msea_dotplot(
  msea_df = msea_results_flat,
  padj_threshold = 0.05,
  top_n_if_none = 10
)
