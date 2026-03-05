# =========================
# 0) Packages
# =========================
# install.packages(c("tidyverse", "tidymodels", "ranger", "vip", "janitor", "caret"))
rm(list = ls())
library(tidyverse)
library(tidymodels)
library(ranger)
library(vip)
library(janitor)
library(caret)

tidymodels::tidymodels_prefer()
set.seed(123)

# =========================
# 1) Read data + clean names
# =========================
cube <- read_csv("Data/all_cube_variables_used_in_Lasso.csv") %>% clean_names()
meta <- read_csv("Data/EC_Field_metadata_and_Geospatial.csv") %>% clean_names() %>%
  dplyr::select(-pct_shrub) # drop shrub for correlation (caret will pick one to drop anyway)
names(meta)[2] = 'ph'
names(cube)[1] = 'parent_id'
cube$parent_id = gsub('_all','',cube$parent_id)

dat0 <- cube %>%
  left_join(meta, by = "parent_id")

# =========================
# 2) Define outcomes + predictors (EXACTLY per agreement)
# =========================
outcomes <- c(
  "cube_effect_size_atp_picomoles_per_g",
  "cube_median_atp_picomoles_per_g",
  "cube_percent_fine_sand"
)

num_preds <- c(
  "ph",
  "water_depth_cm",
  "stream_order",
  "total_drainage_area_sq_km",
  "slope",
  "elevation",
  "aridity_ws",
  "precipitation",
  "pct_fst",
  "pct_ag"
)

cat_preds_to_encode <- c(
  "hydrogeomorphology",
  "sediment",
  "canopy_coverage",
  "macrophyte_coverage",
  "algal_mat_coverage"
)

keep_cols <- c("parent_id", outcomes, num_preds, cat_preds_to_encode)

# Keep only those columns (some might not exist if names differ; catch early)
missing_cols <- setdiff(keep_cols, names(dat0))
if (length(missing_cols) > 0) {
  stop("These expected columns were not found after clean_names():\n",
       paste(missing_cols, collapse = ", "))
}

dat1 <- dat0 %>% select(all_of(keep_cols)) 
# =========================
# 3) Remove rows with missing data & print removed Parent_IDs
# =========================
rows_with_na <- dat1 %>%
  mutate(any_na = if_any(everything(), is.na)) %>%
  filter(any_na) %>%
  pull(parent_id)

cat("\n============================\n")
cat("Rows removed due to missing data (Parent_ID):\n")
if (length(rows_with_na) == 0) {
  cat("None\n")
} else {
  print(unique(rows_with_na))
}
cat("============================\n\n")

dat_complete <- dat1 %>%
  filter(if_all(everything(), ~ !is.na(.x)))

cat("Complete-case n =", nrow(dat_complete), "\n\n")

# =========================
# 4) Correlation screen on numeric predictors
#    - Print high-corr pairs
#    - Drop variables (caret::findCorrelation)
# =========================
num_df <- dat_complete %>% select(all_of(num_preds)) 

cor_mat <- cor(num_df, use = "pairwise.complete.obs")

# print high correlation pairs
thr <- 0.60
high_pairs <- which(abs(cor_mat) >= thr & abs(cor_mat) < 1, arr.ind = TRUE)

cat("\n============================\n")
cat("Highly correlated numeric pairs (|r| >=", thr, "):\n")
if (nrow(high_pairs) == 0) {
  cat("None\n")
} else {
  high_tbl <- tibble(
    var1 = rownames(cor_mat)[high_pairs[,1]],
    var2 = colnames(cor_mat)[high_pairs[,2]],
    r = cor_mat[high_pairs]
  ) %>%
    # keep each pair once
    rowwise() %>%
    mutate(pair = paste(sort(c(var1, var2)), collapse = " ~~ ")) %>%
    ungroup() %>%
    distinct(pair, .keep_all = TRUE) %>%
    arrange(desc(abs(r)))
  print(high_tbl)
}
cat("============================\n\n")

# drop correlated variables (caret chooses ones that reduce overall redundancy)
to_drop_idx <- findCorrelation(cor_mat, cutoff = thr, names = FALSE, exact = TRUE)
to_drop <- if (length(to_drop_idx) > 0) colnames(cor_mat)[to_drop_idx] else character(0)

num_kept <- setdiff(num_preds, to_drop)

cat("\n============================\n")
cat("Numeric variables dropped due to correlation:\n")
if (length(to_drop) == 0) {
  cat("None\n")
} else {
  print(to_drop)
}

cat("\nNumeric variables kept:\n")
print(num_kept)
cat("============================\n\n")

# =========================
# 5) Build a function to fit RF for ONE outcome using bootstrap resampling
#    - No train/test split
#    - Encode ONLY your 9 categorical variables
#    - Use permutation importance
# =========================
fit_rf_bootstrap <- function(data, outcome, num_kept, boots = 200, preview_n = 6) {
  
  # ---- define encoding groups (ONLY these 5 variables are used from categoricals) ----
  nominal_kminus1 <- c("hydrogeomorphology", "sediment")
  
  ordinal_cov <- c("canopy_coverage", "macrophyte_coverage", "algal_mat_coverage")
  cov_levels  <- c("No coverage", "Partial coverage", "Full coverage")  # 0 < 1 < 2
  
  # ---- build modeling dataset ----
  d <- data %>%
    select(parent_id, all_of(outcome), all_of(num_kept),
           all_of(nominal_kminus1), all_of(ordinal_cov)) %>%
    mutate(
      # Nominal: factor (and set baseline/reference level explicitly)
      hydrogeomorphology = factor(hydrogeomorphology),
      sediment = factor(sediment),
      
      # Ordinal: ordered factor with explicit ordering
      canopy_coverage = factor(canopy_coverage, levels = cov_levels, ordered = TRUE),
      macrophyte_coverage = factor(macrophyte_coverage, levels = cov_levels, ordered = TRUE),
      algal_mat_coverage = factor(algal_mat_coverage, levels = cov_levels, ordered = TRUE)
    )
  
  # ---- choose reference levels for nominal vars (the level that will be DROPPED) ----
  # Change these baselines if you want a different dropped category.
  # Good practice: pick a common or conceptually "neutral" baseline.
  d <- d %>%
    mutate(
      hydrogeomorphology = relevel(hydrogeomorphology, ref = "Single-channel straight"),
      sediment = relevel(sediment, ref = "Sand")
    )
  
  # ---- print encoding details BEFORE modeling ----
  cat("\n============================\n")
  cat("Outcome:", outcome, "\n\n")
  
  cat("Nominal vars (k-1 dummy encoding; one_hot = FALSE):\n")
  print(nominal_kminus1)
  cat("\nReference (dropped) level for hydrogeomorphology:", levels(d$hydrogeomorphology)[1], "\n")
  cat("Reference (dropped) level for sediment:", levels(d$sediment)[1], "\n\n")
  
  cat("Levels for nominal vars:\n")
  cat("\n--- hydrogeomorphology ---\n"); print(levels(d$hydrogeomorphology))
  cat("\n--- sediment ---\n");         print(levels(d$sediment))
  
  cat("\nOrdinal coverage vars (ordinal score 0/1/2 via step_ordinalscore):\n")
  print(ordinal_cov)
  cat("\nOrdered levels (0 < 1 < 2):\n")
  print(cov_levels)
  cat("============================\n")
  
  # ---- recipe: remove id, encode nominal as k-1 dummies, ordinal as score ----
  rec <- recipe(as.formula(paste(outcome, "~ .")), data = d) %>%
    update_role(parent_id, new_role = "id") %>%
    step_rm(has_role("id")) %>%
    
    # k-1 dummy encoding (drops the reference level)
    step_dummy(all_of(nominal_kminus1), one_hot = FALSE) %>%
    
    # ordinal score encoding (creates numeric columns)
    # canopy_coverage -> canopy_coverage_score (0/1/2)
    step_ordinalscore(all_of(ordinal_cov))
  
  # ---- preview baked columns to confirm encoding ----
  rec_prep <- prep(rec, training = d, verbose = FALSE)
  baked <- bake(rec_prep, new_data = NULL)
  
  # dummy columns created for nominal vars
  dummy_cols <- unlist(lapply(nominal_kminus1, function(v) {
    grep(paste0("^", v, "_"), names(baked), value = TRUE)
  }))
  
  # ordinal score column names created by step_ordinalscore()
  # usually adds "_score"
  score_cols <- paste0(ordinal_cov, "_score")
  
  cat("\n============================\n")
  cat("Dummy columns created (k-1; reference level dropped):\n")
  print(dummy_cols)
  
  cat("\nOrdinal score columns created (0/1/2):\n")
  print(score_cols)
  
  # -------------------------
  # 5) Random Forest spec (regression) + workflow
  # -------------------------
  rf_spec <- rand_forest(
    trees = 1000,
    mtry  = tune(),
    min_n = tune()
  ) %>%
    set_engine("ranger", importance = "permutation") %>%
    set_mode("regression")
  
  wf <- workflow() %>% add_recipe(rec) %>% add_model(rf_spec)
  
  # -------------------------
  # 6) Bootstrap resampling (no train/test split)
  # -------------------------
  rs <- bootstraps(d, times = boots)
  
  # Determine range for mtry using number of predictors after encoding
  p <- ncol(baked) - 1  # minus outcome column
  
  grid <- grid_regular(
    mtry(range = c(2, max(2, floor(sqrt(p)) * 2))),
    min_n(range = c(2, 25)),
    levels = 5
  )
  
  # -------------------------
  # 7) Tune + select best by R^2
  # -------------------------
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
  
  # -------------------------
  # 8) Bootstrap performance (R^2 + RMSE)
  # -------------------------
  final_rs <- fit_resamples(
    final_wf,
    resamples = rs,
    metrics = metric_set(rsq, rmse),
    control = control_resamples(save_pred = FALSE)
  )
  
  perf <- collect_metrics(final_rs)
  
  # -------------------------
  # 9) Fit on ALL data for permutation importance table
  # -------------------------
  final_fit <- fit(final_wf, data = d)
  ranger_fit <- extract_fit_parsnip(final_fit)$fit
  
  vip_tbl <- tibble(
    variable = names(ranger_fit$variable.importance),
    importance = as.numeric(ranger_fit$variable.importance)
  ) %>% arrange(desc(importance))
  
  list(
    outcome = outcome,
    best = best,
    bootstrap_metrics = perf,
    vip = vip_tbl,
    baked_preview = baked %>% select(any_of(dummy_cols)) %>% head(preview_n),
    dummy_columns = dummy_cols,
    final_workflow = final_fit
  )
}

# =========================
# 6) Run the 3 models
# =========================
results <- purrr::map(outcomes, ~fit_rf_bootstrap(
  data = dat_complete,
  outcome = .x,
  num_kept = num_kept,
  boots = 200
))

# =========================
# 7) Print summary R^2 (bootstrap estimate)
# =========================
summary_tbl <- map_dfr(results, function(res) {
  res$bootstrap_metrics %>%
    mutate(outcome = res$outcome) %>%
    select(outcome, .metric, mean, std_err)
})

cat("\n============================\n")
cat("Bootstrap performance summary (mean ± std_err):\n")
print(summary_tbl)
cat("============================\n\n")

# =========================
# 8) Variable importance outputs
# =========================
for (res in results) {
  cat("\n============================\n")
  cat("Top 20 variable importances for:", res$outcome, "\n")
  print(res$vip %>% slice_head(n = 20))
  cat("============================\n\n")
  
  p <- res$vip %>%
    slice_head(n = 25) %>%
    ggplot(aes(x = reorder(variable, importance), y = importance)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste("Permutation importance:", res$outcome),
      x = NULL,
      y = "Importance"
    )
  print(p)
}
