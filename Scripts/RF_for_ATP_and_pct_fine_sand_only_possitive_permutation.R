# =========================
# 0) Packages
# =========================
rm(list = ls())

library(tidyverse)
library(tidymodels)
library(ranger)
library(vip)
library(janitor)
library(caret)

library(patchwork)
library(writexl)
library(stringr)

tidymodels::tidymodels_prefer()
set.seed(123)

theme_set(theme_bw())

# =========================
# 1) Read data + clean names
# =========================
cube <- read_csv("Data/all_cube_variables_used_in_Lasso.csv") %>% clean_names()
meta <- read_csv("Data/EC_Field_metadata_and_Geospatial.csv") %>% clean_names()

# Match your naming expectations
names(meta)[2] <- "ph"
names(cube)[1] <- "parent_id"
cube$parent_id <- gsub("_all", "", cube$parent_id)

dat0 <- cube %>% left_join(meta, by = "parent_id")

# =========================
# 2) Define outcomes + predictors
# =========================
outcomes <- c(
  "cube_median_atp_picomoles_per_g",
  "cube_percent_fine_sand"
)

num_preds <- c(
  "water_depth_cm",
  "stream_order",
  "total_drainage_area_sq_km",
  "slope",
  "elevation",
  "aridity_ws",
  "pct_fst",
  "pct_ag"
)

# Nominal (k-1 dummies) — available in data, but we will MANUALLY select per outcome below
nominal_all <- c(
  "hydrogeomorphology",
  "sediment"
)

# Ordinal coverage (0/1/2) — optional columns
ordinal_cov_all <- c(
  "canopy_coverage",
  "macrophyte_coverage",
  "algal_mat_coverage"
)

# =========================
# MANUAL predictor selection based on INITIAL VIP plots
# (raw-column level; factors kept as whole factors)
# =========================
keep_map <- list(
  "cube_median_atp_picomoles_per_g" = list(
    num = c(
      "aridity_ws",
      "slope",
      "total_drainage_area_sq_km",
      "pct_ag"
      # dropping pct_fst based on negative perm importance
    ),
    nominal = c(
      "hydrogeomorphology",
      "sediment"
    ),
    ordinal = c(
      "canopy_coverage"
    )
  ),
  
  "cube_percent_fine_sand" = list(
    num = c(
      "pct_ag",
      "slope",
      "elevation",
      "aridity_ws",
      "stream_order",
      "water_depth_cm"
      # dropping total_drainage_area_sq_km and pct_fst based on ~0/negative
    ),
    nominal = c(
      "sediment"
      # dropping hydrogeomorphology based on negative/near 0 i
    ),
    ordinal = c(
      "canopy_coverage"
      # dropping algal_mat_coverage based on negative 
    )
  )
)

# =========================
# Required columns check 
# =========================
# Ensure all outcomes exist
missing_outcomes <- setdiff(outcomes, names(dat0))
if (length(missing_outcomes) > 0) {
  stop("These OUTCOME columns were not found:\n", paste(missing_outcomes, collapse = ", "))
}

# Ensure manual selections exist (numeric + nominal); ordinal optional
for (o in outcomes) {
  miss_num <- setdiff(keep_map[[o]]$num, names(dat0))
  miss_nom <- setdiff(keep_map[[o]]$nominal, names(dat0))
  if (length(miss_num) > 0 || length(miss_nom) > 0) {
    stop(
      "\nFor outcome: ", o,
      "\nMissing numeric predictors: ", paste(miss_num, collapse = ", "),
      "\nMissing nominal predictors: ", paste(miss_nom, collapse = ", "),
      "\n(Ordinal predictors are optional.)\n"
    )
  }
}

# Keep only ordinal columns that exist
available_ordinal <- intersect(ordinal_cov_all, names(dat0))
missing_optional <- setdiff(ordinal_cov_all, names(dat0))
if (length(missing_optional) > 0) {
  cat("\nNOTE: These OPTIONAL ordinal coverage columns were not found and will be skipped:\n")
  print(missing_optional)
  cat("\n")
}

# Update keep_map ordinal lists to available columns only
for (o in outcomes) {
  keep_map[[o]]$ordinal <- intersect(keep_map[[o]]$ordinal, available_ordinal)
}

# =========================
# Label helper functions
# =========================
pretty_outcome <- function(x) {
  dplyr::case_when(
    x == "cube_effect_size_atp_picomoles_per_g" ~ "ATP effect size",
    x == "cube_median_atp_picomoles_per_g"      ~ "ATP median",
    x == "cube_percent_fine_sand"               ~ "% fine sand",
    TRUE ~ x
  )
}

pretty_var <- function(x) {
  x %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_replace_all("percent", "%") %>%
    stringr::str_replace_all("\\bpct\\b", "%") %>%
    stringr::str_replace_all("sq km", "km^2") %>%
    stringr::str_trim()
}

wrap_lab <- function(x, width = 22) stringr::str_wrap(x, width = width)

# =========================
# 3) Create a dataset containing only columns we might use
# =========================
all_num <- unique(unlist(lapply(keep_map, `[[`, "num")))
all_nom <- unique(unlist(lapply(keep_map, `[[`, "nominal")))
all_ord <- unique(unlist(lapply(keep_map, `[[`, "ordinal")))

keep_cols <- unique(c("parent_id", outcomes, all_num, all_nom, all_ord))
dat1 <- dat0 %>% select(all_of(keep_cols))

# =========================
# 4) Remove rows with missing data & print removed Parent_IDs
# =========================
rows_with_na <- dat1 %>%
  mutate(any_na = if_any(everything(), is.na)) %>%
  filter(any_na) %>%
  pull(parent_id) %>%
  unique()

cat("\n============================\n")
cat("Rows removed due to missing data (Parent_ID):\n")
if (length(rows_with_na) == 0) cat("None\n") else print(rows_with_na)
cat("============================\n\n")

dat_complete <- dat1 %>%
  filter(if_all(everything(), ~ !is.na(.x)))

cat("Complete-case n =", nrow(dat_complete), "\n\n")

# =========================
# 5) Correlation screen on numeric predictors (ONLY those used anywhere)
# =========================
if (length(all_num) >= 2) {
  
  num_df <- dat_complete %>% select(all_of(all_num))
  cor_mat <- cor(num_df, use = "pairwise.complete.obs")
  
  thr <- 0.75
  
  high_pairs <- which(abs(cor_mat) >= thr & abs(cor_mat) < 1, arr.ind = TRUE)
  
  cat("\n============================\n")
  cat("Highly correlated numeric pairs (|r| >=", thr, "):\n")
  if (nrow(high_pairs) == 0) {
    cat("None\n")
  } else {
    high_tbl <- tibble(
      var1 = rownames(cor_mat)[high_pairs[, 1]],
      var2 = colnames(cor_mat)[high_pairs[, 2]],
      r = cor_mat[high_pairs]
    ) %>%
      rowwise() %>%
      mutate(pair = paste(sort(c(var1, var2)), collapse = " ~~ ")) %>%
      ungroup() %>%
      distinct(pair, .keep_all = TRUE) %>%
      arrange(desc(abs(r)))
    print(high_tbl)
  }
  cat("============================\n\n")
  
  to_drop_idx <- findCorrelation(cor_mat, cutoff = thr, exact = TRUE)
  to_drop <- if (length(to_drop_idx) > 0) colnames(cor_mat)[to_drop_idx] else character(0)
  num_kept_global <- setdiff(all_num, to_drop)
  
  cat("\n============================\n")
  cat("Numeric variables dropped due to correlation:\n")
  if (length(to_drop) == 0) cat("None\n") else print(to_drop)
  
  cat("\nNumeric variables kept (global):\n")
  print(num_kept_global)
  cat("============================\n\n")
  
  # Apply the drops to keep_map per outcome
  for (o in outcomes) {
    keep_map[[o]]$num <- intersect(keep_map[[o]]$num, num_kept_global)
  }
  
} else {
  cat("\nNOTE: <2 numeric predictors total; skipping correlation screen.\n")
}

# =========================
# Helper: Bootstrap stability of variable importance
# =========================
bootstrap_importance <- function(final_wf, d, n_boot = 200, seed = 123) {
  set.seed(seed)
  
  imps <- purrr::map_dfr(seq_len(n_boot), function(b) {
    d_boot <- d %>% slice_sample(n = nrow(d), replace = TRUE)
    fit_b <- fit(final_wf, data = d_boot)
    ranger_b <- extract_fit_parsnip(fit_b)$fit
    
    tibble(
      boot = b,
      variable = names(ranger_b$variable.importance),
      importance = as.numeric(ranger_b$variable.importance)
    )
  })
  
  imp_summary <- imps %>%
    group_by(variable) %>%
    summarise(
      median_imp = median(importance, na.rm = TRUE),
      q25 = quantile(importance, 0.25, na.rm = TRUE),
      q75 = quantile(importance, 0.75, na.rm = TRUE),
      prop_gt0 = mean(importance > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(median_imp))
  
  rank_summary <- imps %>%
    group_by(boot) %>%
    mutate(rank = rank(-importance, ties.method = "average")) %>%
    ungroup() %>%
    group_by(variable) %>%
    summarise(
      median_rank = median(rank, na.rm = TRUE),
      q25_rank = quantile(rank, 0.25, na.rm = TRUE),
      q75_rank = quantile(rank, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(median_rank)
  
  list(imps = imps, imp_summary = imp_summary, rank_summary = rank_summary)
}

# =========================
# 6) Fit RF for ONE outcome using repeated CV + OOB reporting
# =========================
fit_rf_explain <- function(
    data,
    outcome,
    keep_map,
    v = 5,
    repeats = 50,
    imp_boots = 200,
    preview_n = 6
) {
  
  cov_levels <- c("No coverage", "Partial coverage", "Full coverage")
  
  sel_num     <- keep_map[[outcome]]$num
  sel_nominal <- keep_map[[outcome]]$nominal
  sel_ordinal <- keep_map[[outcome]]$ordinal
  
  d <- data %>%
    select(
      parent_id,
      all_of(outcome),
      all_of(sel_num),
      all_of(sel_nominal),
      any_of(sel_ordinal)
    )
  
  # factor setup only for factors you kept
  if ("hydrogeomorphology" %in% sel_nominal) d <- d %>% mutate(hydrogeomorphology = factor(hydrogeomorphology))
  if ("sediment" %in% sel_nominal) d <- d %>% mutate(sediment = factor(sediment))
  
  # ordinal vars if present
  if (length(sel_ordinal) > 0) {
    d <- d %>%
      mutate(across(any_of(sel_ordinal), ~ factor(.x, levels = cov_levels, ordered = TRUE)))
  }
  
  # reference levels (only if present AND levels exist)
  if ("hydrogeomorphology" %in% sel_nominal && "Single-channel straight" %in% levels(d$hydrogeomorphology)) {
    d <- d %>% mutate(hydrogeomorphology = relevel(hydrogeomorphology, ref = "Single-channel straight"))
  }
  if ("sediment" %in% sel_nominal && "Sand" %in% levels(d$sediment)) {
    d <- d %>% mutate(sediment = relevel(sediment, ref = "Sand"))
  }
  
  cat("\n============================\n")
  cat("Outcome:", outcome, "\n\n")
  cat("Numeric predictors USED:\n")
  print(sel_num)
  
  cat("\nNominal predictors USED (k-1 dummy encoding; one_hot = FALSE):\n")
  print(sel_nominal)
  
  if ("hydrogeomorphology" %in% sel_nominal) {
    cat("\nReference (dropped) hydrogeomorphology:", levels(d$hydrogeomorphology)[1], "\n")
  }
  if ("sediment" %in% sel_nominal) {
    cat("Reference (dropped) sediment:", levels(d$sediment)[1], "\n")
  }
  
  cat("\nOrdinal coverage vars USED (step_ordinalscore -> 0/1/2):\n")
  if (length(sel_ordinal) == 0) cat("None present/selected\n") else print(sel_ordinal)
  cat("============================\n")
  
  rec <- recipe(as.formula(paste(outcome, "~ .")), data = d) %>%
    update_role(parent_id, new_role = "id") %>%
    step_rm(has_role("id"))
  
  # dummies only for nominal vars kept
  if (length(sel_nominal) > 0) {
    rec <- rec %>% step_dummy(all_of(sel_nominal), one_hot = FALSE)
  }
  
  # ordinalscore only for ordinal vars kept
  if (length(sel_ordinal) > 0) {
    rec <- rec %>% step_ordinalscore(all_of(sel_ordinal))
  }
  
  # preview baked columns
  rec_prep <- prep(rec, training = d, verbose = FALSE)
  baked <- bake(rec_prep, new_data = NULL)
  
  dummy_cols <- if (length(sel_nominal) > 0) {
    unlist(lapply(sel_nominal, function(vv) {
      grep(paste0("^", vv, "_"), names(baked), value = TRUE)
    }))
  } else character(0)
  
  score_cols <- if (length(sel_ordinal) > 0) paste0(sel_ordinal, "_score") else character(0)
  
  cat("\n============================\n")
  cat("Dummy columns created:\n")
  if (length(dummy_cols) == 0) cat("None\n") else print(dummy_cols)
  cat("\nOrdinal score columns created:\n")
  if (length(score_cols) == 0) cat("None\n") else print(score_cols)
  cat("\nPreview of encoded columns (first ", preview_n, " rows):\n", sep = "")
  print(baked %>% select(any_of(dummy_cols), any_of(score_cols)) %>% head(preview_n))
  cat("============================\n\n")
  
  # RF spec
  rf_spec <- rand_forest(
    trees = 1000,
    mtry = tune(),
    min_n = tune()
  ) %>%
    set_engine("ranger", importance = "permutation") %>%
    set_mode("regression")
  
  wf <- workflow() %>% add_recipe(rec) %>% add_model(rf_spec)
  
  rs <- vfold_cv(d, v = v, repeats = repeats)
  
  # predictor count post-encoding for mtry
  p <- ncol(baked) - 1
  
  grid <- grid_regular(
    mtry(range = c(2, max(3, floor(sqrt(p))))),
    min_n(range = c(12, 25)),
    levels = 5
  )
  
  tuned <- tune_grid(
    wf,
    resamples = rs,
    grid = grid,
    metrics = metric_set(rsq, rmse)
  )
  
  best <- select_best(tuned, metric = "rsq")
  cat("\nBest params for", outcome, ":\n")
  print(best)
  
  final_wf <- finalize_workflow(wf, best)
  
  final_rs <- fit_resamples(
    final_wf,
    resamples = rs,
    metrics = metric_set(rsq, rmse),
    control = control_resamples(save_pred = FALSE)
  )
  
  perf <- collect_metrics(final_rs)
  
  rsq_dist <- collect_metrics(final_rs, summarize = FALSE) %>%
    filter(.metric == "rsq") %>%
    summarise(
      q10 = quantile(.estimate, 0.10, na.rm = TRUE),
      q50 = quantile(.estimate, 0.50, na.rm = TRUE),
      q90 = quantile(.estimate, 0.90, na.rm = TRUE)
    )
  
  cat("\nRepeated CV R^2 quantiles (10/50/90%):\n")
  print(rsq_dist)
  
  # Fit on ALL data
  final_fit <- fit(final_wf, data = d)
  ranger_fit <- extract_fit_parsnip(final_fit)$fit
  
  oob_r2 <- ranger_fit$r.squared
  oob_rmse <- sqrt(ranger_fit$prediction.error)
  
  cat("\nOOB metrics from final ranger fit:\n")
  cat("OOB R^2:", round(oob_r2, 4), "\n")
  cat("OOB RMSE:", round(oob_rmse, 4), "\n\n")
  
  vip_tbl <- tibble(
    variable = names(ranger_fit$variable.importance),
    importance = as.numeric(ranger_fit$variable.importance)
  ) %>% arrange(desc(importance))
  
  # Importance stability
  imp_stab <- bootstrap_importance(final_wf, d, n_boot = imp_boots, seed = 123)
  
  cat("\n============================\n")
  cat("Importance stability (top 15 by median importance):\n")
  print(imp_stab$imp_summary %>% slice_head(n = 15))
  cat("\nRank stability (top 15 by median rank; lower is better):\n")
  print(imp_stab$rank_summary %>% slice_head(n = 15))
  cat("============================\n\n")
  
  list(
    outcome = outcome,
    best = best,
    cv_metrics = perf,
    rsq_quantiles = rsq_dist,
    oob_r2 = oob_r2,
    oob_rmse = oob_rmse,
    vip = vip_tbl,
    importance_boot = imp_stab,
    baked_preview = baked %>% select(any_of(dummy_cols), any_of(score_cols)) %>% head(preview_n),
    dummy_columns = dummy_cols,
    score_columns = score_cols,
    final_fit = final_fit
  )
}

# =========================
# 7) Run models
# =========================
results <- purrr::map(outcomes, ~fit_rf_explain(
  data = dat_complete,
  outcome = .x,
  keep_map = keep_map,
  v = 5,
  repeats = 50,
  imp_boots = 200
))

# =========================
# 8) Print summary performance
# =========================
summary_tbl <- purrr::map_dfr(results, function(res) {
  res$cv_metrics %>%
    mutate(outcome = res$outcome) %>%
    select(outcome, .metric, mean, std_err)
})

cat("\n============================\n")
cat("Repeated CV performance summary (mean ± std_err):\n")
print(summary_tbl)
cat("============================\n\n")

# =========================
# 9) Export GOF + tuning params + importance stability to Excel
# =========================
perf_mean_se <- purrr::map_dfr(results, function(res) {
  res$cv_metrics %>%
    mutate(outcome = pretty_outcome(res$outcome)) %>%
    select(outcome, .metric, mean, std_err)
})

rsq_quantiles_tbl <- purrr::map_dfr(results, function(res) {
  res$rsq_quantiles %>%
    mutate(outcome = pretty_outcome(res$outcome))
}) %>%
  select(outcome, everything())

best_params_tbl <- purrr::map_dfr(results, function(res) {
  res$best %>%
    mutate(outcome = pretty_outcome(res$outcome))
}) %>%
  select(outcome, everything())

oob_tbl <- purrr::map_dfr(results, function(res) {
  tibble(
    outcome = pretty_outcome(res$outcome),
    oob_r2 = res$oob_r2,
    oob_rmse = res$oob_rmse
  )
})

gof_wide <- perf_mean_se %>%
  tidyr::pivot_wider(names_from = .metric, values_from = c(mean, std_err)) %>%
  left_join(rsq_quantiles_tbl, by = "outcome") %>%
  left_join(oob_tbl, by = "outcome") %>%
  left_join(best_params_tbl, by = "outcome")

imp_summary_tbl <- purrr::map_dfr(results, function(res) {
  res$importance_boot$imp_summary %>%
    mutate(outcome = pretty_outcome(res$outcome)) %>%
    select(outcome, everything())
})

rank_summary_tbl <- purrr::map_dfr(results, function(res) {
  res$importance_boot$rank_summary %>%
    mutate(outcome = pretty_outcome(res$outcome)) %>%
    select(outcome, everything())
})

write_xlsx(
  list(
    "gof_wide" = gof_wide,
    "perf_mean_se" = perf_mean_se,
    "rsq_quantiles" = rsq_quantiles_tbl,
    "oob_metrics" = oob_tbl,
    "best_params" = best_params_tbl,
    "imp_summary_prop_gt0" = imp_summary_tbl,
    "rank_summary" = rank_summary_tbl
  ),
  path = "rf_model_gof_manual_positive_only.xlsx"
)

cat("\nWrote model fit summary to: rf_model_gof_manual_positive_only.xlsx\n")

# =========================
# 10) Variable importance plots (2-panel tagged A/B)
# =========================
vip_plots <- list()
stab_plots <- list()

for (res in results) {
  
  vip_plots[[res$outcome]] <-
    res$vip %>%
    slice_head(n = 25) %>%
    mutate(variable_pretty = wrap_lab(pretty_var(variable))) %>%
    ggplot(aes(x = reorder(variable_pretty, importance), y = importance)) +
    geom_col() +
    coord_flip() +
    labs(
      title = pretty_outcome(res$outcome),
      x = NULL,
      y = "Permutation importance"
    )
  
  top_vars <- res$importance_boot$imp_summary %>%
    slice_head(n = 15) %>%
    pull(variable)
  
  stab_plots[[res$outcome]] <-
    res$importance_boot$imps %>%
    filter(variable %in% top_vars) %>%
    mutate(variable_pretty = wrap_lab(pretty_var(variable))) %>%
    ggplot(aes(x = reorder(variable_pretty, importance, FUN = median), y = importance)) +
    geom_boxplot(outlier.size = 0.7) +
    coord_flip() +
    labs(
      title = pretty_outcome(res$outcome),
      x = NULL,
      y = "Bootstrap importance"
    )
}

vip_plots_2 <- vip_plots[outcomes]
stab_plots_2 <- stab_plots[outcomes]

vip_panel_2 <-
  wrap_plots(vip_plots_2, ncol = 2) +
  plot_annotation(
    title = "Permutation importance (single fit on all data)",
    tag_levels = "A"
  )

stab_panel_2 <-
  wrap_plots(stab_plots_2, ncol = 2) +
  plot_annotation(
    title = "Importance stability (bootstrap distributions)",
    tag_levels = "A"
  )

print(vip_panel_2)
print(stab_panel_2)

ggsave("rf_vip_2panel_AB_manual_positive_only.png", vip_panel_2, width = 12, height = 5, dpi = 300)
ggsave("rf_importance_stability_2panel_AB_manual_positive_only.png", stab_panel_2, width = 12, height = 5, dpi = 300)

cat("\nSaved figures:\n")
cat(" - rf_vip_2panel_AB_manual_positive_only.png\n")
cat(" - rf_importance_stability_2panel_AB_manual_positive_only.png\n")

# =========================
# SHAP + Partial Dependence
# =========================

library(fastshap)
library(pdp)

get_baked_x_from_fitted_workflow <- function(final_fit_wf, raw_d, outcome) {
  
  # extract the recipe that is already trained/prepped inside the fitted workflow
  rec_trained <- workflows::extract_recipe(final_fit_wf)
  
  baked <- bake(rec_trained, new_data = raw_d)
  
  X <- baked %>% dplyr::select(-all_of(outcome))
  X
}

# ---- helper: SHAP for ranger regression ----
compute_shap <- function(ranger_fit, X, nsim = 300, seed = 123) {
  set.seed(seed)
  
  pred_wrapper <- function(object, newdata) {
    predict(object, data = newdata)$predictions
  }
  
  shap <- fastshap::explain(
    object = ranger_fit,
    X = X,
    pred_wrapper = pred_wrapper,
    nsim = nsim,
    adjust = TRUE
  )
  
  as_tibble(shap)
}

# ---- helper: summarize SHAP directionality per feature ----
summarize_shap_direction <- function(shap_tbl, X) {
  feats <- colnames(X)
  
  purrr::map_dfr(feats, function(v) {
    sv <- shap_tbl[[v]]
    xv <- X[[v]]
    
    tibble(
      variable = v,
      mean_abs_shap = mean(abs(sv), na.rm = TRUE),
      mean_shap = mean(sv, na.rm = TRUE),
      prop_positive_shap = mean(sv > 0, na.rm = TRUE),
      cor_shap_x = suppressWarnings(cor(sv, xv, use = "pairwise.complete.obs"))
    )
  }) %>%
    arrange(desc(mean_abs_shap)) %>%
    mutate(
      direction = case_when(
        is.na(cor_shap_x) ~ NA_character_,
        cor_shap_x >  0.10 ~ "increasing",
        cor_shap_x < -0.10 ~ "decreasing",
        TRUE ~ "weak/none"
      )
    )
}

# ---- helper: PD for a set of variables + slope sign ----
compute_pd_slopes <- function(ranger_fit, X, vars, grid.resolution = 30) {
  
  pred.fun <- function(object, newdata) {
    predict(object, data = newdata)$predictions
  }
  
  purrr::map_dfr(vars, function(v) {
    
    pd <- pdp::partial(
      object = ranger_fit,
      pred.var = v,
      train = X,
      pred.fun = pred.fun,
      grid.resolution = grid.resolution,
      progress = "none"
    )
    
    x <- pd[[v]]
    yhat <- pd$yhat
    
    slope <- suppressWarnings(coef(lm(yhat ~ x))[2])
    
    tibble(
      variable = v,
      pd_slope = slope,
      pd_direction = case_when(
        is.na(slope) ~ NA_character_,
        slope > 0 ~ "increasing",
        slope < 0 ~ "decreasing",
        TRUE ~ "flat"
      )
    )
  })
}

# =========================
# Run SHAP + PD per outcome
# =========================
shap_summaries <- list()
pd_summaries <- list()

for (res in results) {
  
  outcome <- res$outcome
  
  # Recreate the same raw modeling dataframe (same logic as fit_rf_explain)
  sel_num     <- keep_map[[outcome]]$num
  sel_nominal <- keep_map[[outcome]]$nominal
  sel_ordinal <- keep_map[[outcome]]$ordinal
  
  cov_levels <- c("No coverage", "Partial coverage", "Full coverage")
  
  d <- dat_complete %>%
    dplyr::select(
      parent_id,
      all_of(outcome),
      all_of(sel_num),
      all_of(sel_nominal),
      any_of(sel_ordinal)
    )
  
  # factor setup (same as in fit_rf_explain)
  if ("hydrogeomorphology" %in% sel_nominal) d <- d %>% mutate(hydrogeomorphology = factor(hydrogeomorphology))
  if ("sediment" %in% sel_nominal) d <- d %>% mutate(sediment = factor(sediment))
  
  if (length(sel_ordinal) > 0) {
    d <- d %>%
      mutate(across(any_of(sel_ordinal), ~ factor(.x, levels = cov_levels, ordered = TRUE)))
  }
  
  if ("hydrogeomorphology" %in% sel_nominal && "Single-channel straight" %in% levels(d$hydrogeomorphology)) {
    d <- d %>% mutate(hydrogeomorphology = relevel(hydrogeomorphology, ref = "Single-channel straight"))
  }
  if ("sediment" %in% sel_nominal && "Sand" %in% levels(d$sediment)) {
    d <- d %>% mutate(sediment = relevel(sediment, ref = "Sand"))
  }
  
  # Ranger model object
  ranger_fit <- extract_fit_parsnip(res$final_fit)$fit
  
  # Bake X using the already-trained recipe inside the workflow (NO prep())
  X <- get_baked_x_from_fitted_workflow(res$final_fit, d, outcome)
  
  # ---- SHAP ----
  shap_tbl <- compute_shap(ranger_fit, X, nsim = 300, seed = 123)
  
  shap_sum <- summarize_shap_direction(shap_tbl, X) %>%
    mutate(outcome = pretty_outcome(outcome))
  
  shap_summaries[[outcome]] <- shap_sum
  
  # ---- Partial dependence on top SHAP vars ----
  top_vars <- shap_sum %>% slice_head(n = 8) %>% pull(variable)
  
  pd_sum <- compute_pd_slopes(ranger_fit, X, vars = top_vars, grid.resolution = 30) %>%
    mutate(outcome = pretty_outcome(outcome))
  
  pd_summaries[[outcome]] <- pd_sum
}

shap_summary_tbl <- bind_rows(shap_summaries)
pd_summary_tbl   <- bind_rows(pd_summaries)

# =========================
# Export to Excel
# =========================
write_xlsx(
  list(
    "shap_direction_summary" = shap_summary_tbl,
    "pd_direction_summary"   = pd_summary_tbl
  ),
  path = "rf_directionality_shap_pd.xlsx"
)

cat("\nWrote directionality summaries to: rf_directionality_shap_pd.xlsx\n")