# =========================
# Random Forest (Explanation-focused) with n ~ 46
# - Complete-case only (prints removed Parent_IDs)
# - Drops highly correlated numeric predictors (prints pairs + drops)
# - Encoding:
#     * hydrogeomorphology + sediment: k-1 dummy (reference level dropped)
#     * canopy/macrophyte/algal coverage: ordinal score 0/1/2 (if columns exist)
# - Performance:
#     * repeated v-fold CV for tuning + GOF (less optimistic than bootstrap)
#     * ALSO prints OOB R^2 and OOB RMSE from final ranger fit
# - Importance:
#     * permutation importance (ranger)
#     * bootstrap stability for importance (distributions + rank stability + prop_gt0)
# - PLOTTING:
#     * plots combined into multi-panel figures tagged a/b/c automatically
#     * theme_bw() globally
#     * shorter labels
# - OUTPUT:
#     * GOF + quantiles + best params to Excel
#     * importance + stability summaries to Excel
# =========================

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
  # add more outcomes here if desired
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

# Nominal (k-1 dummies)
nominal_kminus1 <- c(
  "hydrogeomorphology",
  "sediment"
)

# Ordinal coverage (0/1/2) — use if present in data
ordinal_cov_all <- c(
  "canopy_coverage",
  "macrophyte_coverage",
  "algal_mat_coverage"
)

# Keep columns:
# - outcomes and numeric predictors are required
# - nominal predictors are required
# - ordinal predictors are optional (handled gracefully)
required_cols <- c("parent_id", outcomes, num_preds, nominal_kminus1)
optional_cols <- ordinal_cov_all

missing_required <- setdiff(required_cols, names(dat0))
if (length(missing_required) > 0) {
  stop(
    "These REQUIRED columns were not found after clean_names():\n",
    paste(missing_required, collapse = ", ")
  )
}

available_ordinal <- intersect(optional_cols, names(dat0))
missing_optional <- setdiff(optional_cols, names(dat0))

if (length(missing_optional) > 0) {
  cat("\nNOTE: These OPTIONAL ordinal coverage columns were not found and will be skipped:\n")
  print(missing_optional)
  cat("\n")
}

keep_cols <- c(required_cols, available_ordinal)

dat1 <- dat0 %>% select(all_of(keep_cols))

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
# 3) Remove rows with missing data & print removed Parent_IDs
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
# 4) Correlation screen on numeric predictors
# =========================
num_df <- dat_complete %>% select(all_of(num_preds))
cor_mat <- cor(num_df, use = "pairwise.complete.obs")

thr <- 0.75  # less aggressive default; change if you want

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
num_kept <- setdiff(num_preds, to_drop)

cat("\n============================\n")
cat("Numeric variables dropped due to correlation:\n")
if (length(to_drop) == 0) cat("None\n") else print(to_drop)

cat("\nNumeric variables kept:\n")
print(num_kept)
cat("============================\n\n")

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
# 5) Fit RF for ONE outcome using repeated CV + OOB reporting
# =========================
fit_rf_explain <- function(
    data,
    outcome,
    num_kept,
    ordinal_cov_present,
    v = 5,
    repeats = 50,
    imp_boots = 200,
    preview_n = 6
) {
  
  cov_levels <- c("No coverage", "Partial coverage", "Full coverage")  # 0 < 1 < 2
  
  d <- data %>%
    select(
      parent_id,
      all_of(outcome),
      all_of(num_kept),
      all_of(nominal_kminus1),
      any_of(ordinal_cov_present)
    ) %>%
    mutate(
      hydrogeomorphology = factor(hydrogeomorphology),
      sediment = factor(sediment)
    ) %>%
    mutate(
      # if present, coerce ordinal vars to ordered factor with defined levels
      across(
        any_of(ordinal_cov_present),
        ~ factor(.x, levels = cov_levels, ordered = TRUE)
      )
    ) %>%
    mutate(
      hydrogeomorphology = relevel(hydrogeomorphology, ref = "Single-channel straight"),
      sediment = relevel(sediment, ref = "Sand")
    )
  
  cat("\n============================\n")
  cat("Outcome:", outcome, "\n\n")
  cat("Nominal vars (k-1 dummy encoding; one_hot = FALSE):\n")
  print(nominal_kminus1)
  cat("\nReference (dropped) hydrogeomorphology:", levels(d$hydrogeomorphology)[1], "\n")
  cat("Reference (dropped) sediment:", levels(d$sediment)[1], "\n\n")
  
  cat("Ordinal coverage vars used (step_ordinalscore -> 0/1/2):\n")
  if (length(ordinal_cov_present) == 0) cat("None present\n") else print(ordinal_cov_present)
  cat("============================\n")
  
  rec <- recipe(as.formula(paste(outcome, "~ .")), data = d) %>%
    update_role(parent_id, new_role = "id") %>%
    step_rm(has_role("id")) %>%
    step_dummy(all_of(nominal_kminus1), one_hot = FALSE)
  
  # only add ordinalscore step if those columns exist
  if (length(ordinal_cov_present) > 0) {
    rec <- rec %>% step_ordinalscore(all_of(ordinal_cov_present))
  }
  
  # preview baked columns
  rec_prep <- prep(rec, training = d, verbose = FALSE)
  baked <- bake(rec_prep, new_data = NULL)
  
  dummy_cols <- unlist(lapply(nominal_kminus1, function(vv) {
    grep(paste0("^", vv, "_"), names(baked), value = TRUE)
  }))
  
  score_cols <- if (length(ordinal_cov_present) > 0) paste0(ordinal_cov_present, "_score") else character(0)
  
  cat("\n============================\n")
  cat("Dummy columns created:\n")
  print(dummy_cols)
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
  
  # repeated v-fold CV (less optimistic at small n than bootstrap resamples)
  rs <- vfold_cv(d, v = v, repeats = repeats)
  
  # predictor count post-encoding for mtry
  p <- ncol(baked) - 1
  
  # conservative grid for small n
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
  
  # CV performance distribution (rsq)
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
  
  # OOB reporting (important for small n)
  # ranger: r.squared is OOB R^2, prediction.error is OOB MSE for regression
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
# 6) Run models
# =========================
results <- purrr::map(outcomes, ~fit_rf_explain(
  data = dat_complete,
  outcome = .x,
  num_kept = num_kept,
  ordinal_cov_present = available_ordinal,
  v = 5,
  repeats = 50,
  imp_boots = 200
))

# =========================
# 7) Print summary performance
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
# 8) Export GOF + tuning params + importance stability to Excel
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

# importance stability summaries (including prop_gt0)
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
  path = "rf_model_gof.xlsx"
)

cat("\nWrote model fit summary to: rf_model_gof.xlsx\n")

# =========================
# 9) Variable importance plots (2-panel tagged A/B)
# =========================
vip_plots <- list()
stab_plots <- list()

for (res in results) {
  
  # VIP plot (top 25)
  vip_plots[[res$outcome]] <-
    res$vip %>%
    slice_head(n = 25) %>%
    mutate(variable_pretty = wrap_lab(pretty_var(variable))) %>%
    ggplot(aes(x = reorder(variable_pretty, importance), y = importance)) +
    geom_col() +
    coord_flip() +
    labs(
      title = pretty_outcome(res$outcome),   # outcome label shown in each panel
      x = NULL,
      y = "Permutation importance"
    )
  
  # Stability plot (boxplots) for top 15 variables
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
      title = pretty_outcome(res$outcome),   # outcome label shown in each panel
      x = NULL,
      y = "Bootstrap importance"
    )
}

# --- Force exactly 2 panels in a fixed order matching outcomes ---
vip_plots_2 <- vip_plots[outcomes]
stab_plots_2 <- stab_plots[outcomes]

vip_panel_2 <-
  wrap_plots(vip_plots_2, ncol = 2) +
  plot_annotation(
    title = "Permutation importance (single fit on all data)",
    tag_levels = "A"   # <-- makes tags A, B
  )

stab_panel_2 <-
  wrap_plots(stab_plots_2, ncol = 2) +
  plot_annotation(
    title = "Importance stability (bootstrap distributions)",
    tag_levels = "A"   # <-- makes tags A, B
  )

print(vip_panel_2)
print(stab_panel_2)

ggsave("rf_vip_2panel_AB.png", vip_panel_2, width = 12, height = 5, dpi = 300)
ggsave("rf_importance_stability_2panel_AB.png", stab_panel_2, width = 12, height = 5, dpi = 300)

cat("\nSaved figures:\n")
cat(" - rf_vip_2panel_AB.png\n")
cat(" - rf_importance_stability_2panel_AB.png\n")