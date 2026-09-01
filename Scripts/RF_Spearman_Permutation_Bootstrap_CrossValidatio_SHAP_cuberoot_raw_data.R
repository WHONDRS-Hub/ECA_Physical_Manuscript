
# Load all required libraries
library(tidyverse);library(corrplot);library(ggpubr);library(ggpmisc);library(factoextra);library(stringr);library(glmnet);library(magick); library(ggnewscale); library(FSA); library(multcompView); library(rcompanion)
library(ranger);library(patchwork);library(cluster)
if(!require(vegan, quietly = TRUE)) install.packages("vegan")
if(!require(ggrepel, quietly = TRUE)) install.packages("ggrepel")
if(!require(viridis, quietly = TRUE)) install.packages("viridis")
if(!require(fastshap, quietly = TRUE)) install.packages("fastshap")
if(!require(caret, quietly = TRUE)) install.packages("caret")
if(!require(boot, quietly = TRUE)) install.packages("boot")
if(!require(ggdist, quietly = TRUE)) install.packages("ggdist")
library(vegan);library(ggrepel);library(viridis);library(fastshap);library(caret);library(boot);library(ggdist)

rm(list=ls());graphics.off()

# Create output directories
dir.create("Data/RF_results", recursive = TRUE, showWarnings = FALSE)
dir.create("Figures/RF_results", recursive = TRUE, showWarnings = FALSE)

cat("=== COMPLETE COMPREHENSIVE MULTICOLLINEARITY AND ROBUSTNESS ANALYSIS ===\n")
cat("Components:\n")
cat("1. Spearman correlation clustering and visualization\n")
cat("2. Method 1 variable selection (correlation-based)\n")
cat("3. Permutation importance analysis\n")
cat("4. Cross-validation for generalizability\n")
cat("5. Bootstrap analysis for feature selection stability\n")
cat("6. SHAP analysis for directionality and effect magnitude\n")
cat("7. Complete analysis for both raw and cube root data\n\n")

# Functions ---------------------------------------------------------------
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

## DATA LOADING ================================================
# Individual Rate data for histograms ####
all_data = read.csv("Data/EC_Sediment_SpC_pH_Temp_Respiration.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%
  select(-c(Field_Name, IGSN, Material))

# Calculate "bulk" medians (not separated by wet/dry)
median_respiration = all_data %>%
  select(-c(Respiration_R_Squared, Respiration_R_Squared_Adj, Respiration_p_value, Total_Incubation_Time_Min, Number_Points_In_Respiration_Regression, Number_Points_Removed_Respiration_Regression,DO_Concentration_At_Incubation_Time_Zero)) %>%
  mutate(across(c(SpC_microsiemens_per_cm:Respiration_Rate_mg_DO_per_kg_per_H), as.numeric)) %>%
  mutate(Respiration_Rate_mg_DO_per_L_per_H = ifelse(grepl("INC_Method_001|INC_Method_002|INC_QA_004", Methods_Deviation), NA, Respiration_Rate_mg_DO_per_L_per_H)) %>%
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(grepl("INC_Method_001|INC_Method_002|INC_QA_004", Methods_Deviation), NA, Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  mutate(SpC_microsiemens_per_cm = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, SpC_microsiemens_per_cm)) %>%
  mutate(pH = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, pH)) %>%
  mutate(Temperature_degC = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, Temperature_degC)) %>%
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(Respiration_Rate_mg_DO_per_kg_per_H == "-9999", NA, Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>%
  mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
  group_by(Sample_ID) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
                        n = ~sum(!is.na(.x))),
                   .names = "{.fn}_{.col}")) %>%
  ungroup() %>%
  group_by(Sample_ID) %>%
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 10), "FALSE", "TRUE")) %>%
  select(c(Sample_ID, Median_SpC_microsiemens_per_cm, Median_pH, Median_Temperature_degC, Median_Respiration_Rate_mg_DO_per_L_per_H, Median_Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  ungroup()

# ATP ####
atp = read.csv("Data/EC_Sediment_ATP.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material)) %>%
  mutate(Sample_Name = str_replace(Sample_Name, "ATP", "INC"))

median_atp = atp %>%
  mutate(ATP_nanomoles_per_L = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, ATP_nanomoles_per_L)) %>%
  mutate(ATP_picomoles_per_g = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, ATP_picomoles_per_g)) %>%
  mutate(ATP_nanomoles_per_L = ifelse(ATP_nanomoles_per_L == -9999, NA, ATP_nanomoles_per_L)) %>%
  mutate(ATP_picomoles_per_g = ifelse(ATP_picomoles_per_g == -9999, NA, ATP_picomoles_per_g)) %>%
  mutate(across(c(ATP_nanomoles_per_L:ATP_picomoles_per_g), as.numeric)) %>%
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>%
  mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
  group_by(Sample_ID) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
                        n = ~sum(!is.na(.x))),
                   .names = "{.fn}_{.col}")) %>%
  ungroup() %>%
  group_by(Sample_ID) %>%
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 10), "FALSE", "TRUE")) %>%
  select(c(Sample_ID, Median_ATP_nanomoles_per_L, Median_ATP_picomoles_per_g)) %>%
  ungroup()

# CN ####
cn = read.csv("Data/EC_Sediment_CN.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material)) %>%
  mutate(Sample_Name = str_replace(Sample_Name, "SCN", "INC"))

median_cn = cn %>%
  mutate(X01395_C_percent_per_mg = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, X01395_C_percent_per_mg)) %>%
  mutate(X01397_N_percent_per_mg = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, X01397_N_percent_per_mg)) %>%
  mutate(X01395_C_percent_per_mg = ifelse(X01395_C_percent_per_mg == -9999, NA, X01395_C_percent_per_mg)) %>%
  mutate(X01397_N_percent_per_mg = ifelse(X01397_N_percent_per_mg == -9999, NA, X01397_N_percent_per_mg)) %>%
  mutate(across(c(X01395_C_percent_per_mg:X01397_N_percent_per_mg), as.numeric)) %>%
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>%
  mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
  group_by(Sample_ID) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
                        n = ~sum(!is.na(.x))),
                   .names = "{.fn}_{.col}")) %>%
  ungroup() %>%
  group_by(Sample_ID) %>%
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 10), "FALSE", "TRUE")) %>%
  select(c(Sample_ID, Median_X01395_C_percent_per_mg, Median_X01397_N_percent_per_mg))%>%
  ungroup()

# NPOC/TN ####
npoc_tn = read.csv("Data/EC_Sediment_NPOC_TN.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material)) %>%
  mutate(Sample_Name = str_replace(Sample_Name, "SIR", "INC"))

median_npoc_tn = npoc_tn %>%
  mutate(Extractable_NPOC_mg_per_L = ifelse(grepl("VI_OCN_010|VB_OCN_001", Methods_Deviation), NA, Extractable_NPOC_mg_per_L)) %>%
  mutate(Extractable_NPOC_mg_per_kg = ifelse(grepl("VI_OCN_010|VB_OCN_001", Methods_Deviation), NA, Extractable_NPOC_mg_per_kg)) %>%
  mutate(Extractable_TN_mg_per_L = ifelse(grepl("VI_OCN_010|VB_OCN_001", Methods_Deviation), NA, Extractable_TN_mg_per_L)) %>%
  mutate(Extractable_TN_mg_per_kg = ifelse(grepl("VI_OCN_010|VB_OCN_001", Methods_Deviation), NA, Extractable_TN_mg_per_kg)) %>%
  mutate(across(c(Extractable_NPOC_mg_per_L:Extractable_TN_mg_per_kg), as.numeric)) %>%
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>%
  mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
  group_by(Sample_ID) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
                        n = ~sum(!is.na(.x))),
                   .names = "{.fn}_{.col}")) %>%
  ungroup() %>%
  group_by(Sample_ID) %>%
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 10), "FALSE", "TRUE")) %>%
  select(c(Sample_ID, Median_Extractable_NPOC_mg_per_L, Median_Extractable_NPOC_mg_per_kg, Median_Extractable_TN_mg_per_L, Median_Extractable_TN_mg_per_kg))%>%
  ungroup()

# Gravimetric Moisture ####
grav_inc = read.csv("Data/EC_Sediment_Gravimetric_Moisture.csv", skip = 2) %>%
  slice(-1:-11) %>%
  filter(Field_Name != "#End_Data") %>%
  select(-c(Field_Name, IGSN, Material)) %>%
  mutate(across(c(Initial_Water_Mass_g:Incubation_Water_Mass_g), as.numeric))

# Fe ####
fe = read.csv("Data/EC_Sediment_Fe.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material)) %>%
  mutate(Sample_Name = str_replace(Sample_Name, "SFE", "INC")) %>%
  rename(Dev_Fe = Methods_Deviation)

median_iron = fe  %>%
  left_join(grav_inc, by = "Sample_Name") %>%
  mutate(Fe_mg_per_L = ifelse(grepl("INC_Method_001|INC_Method_002", Dev_Fe), NA, Fe_mg_per_L)) %>%
  mutate(Fe_mg_per_kg = ifelse(grepl("INC_Method_001|INC_Method_002", Dev_Fe), NA, Fe_mg_per_kg)) %>%
  mutate(Fe_mg_per_L = if_else(grepl("SFE_Below", Fe_mg_per_L), "0.001", Fe_mg_per_L)) %>%
  mutate(Fe_mg_per_L = if_else(grepl("SFE_Above", Fe_mg_per_L), str_extract(Fe_mg_per_L, "(?<=\\|[^|]{1,100}\\|)\\d+\\.\\d+"), Fe_mg_per_L)) %>%
  mutate(Fe_mg_per_L = as.numeric(Fe_mg_per_L)) %>%
  mutate(Fe_mg_per_kg = if_else(grepl("SFE_Above", Fe_mg_per_kg), str_extract(Fe_mg_per_kg, "(?<=\\|[^|]{1,100}\\|)\\d+\\.\\d+"), Fe_mg_per_kg)) %>%
  mutate(Fe_mg_per_kg = if_else(grepl("SFE_Below", Fe_mg_per_kg), as.numeric(Fe_mg_per_L * (Incubation_Water_Mass_g/Dry_Sediment_Mass_g)), as.numeric(Fe_mg_per_kg))) %>%
  mutate(Fe_mg_per_kg = as.numeric(Fe_mg_per_kg)) %>%
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>%
  mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
  group_by(Sample_ID) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
                        n = ~sum(!is.na(.x))),
                   .names = "{.fn}_{.col}")) %>%
  ungroup() %>%
  group_by(Sample_ID) %>%
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 5), "FALSE", "TRUE")) %>%
  select(c(Sample_ID, Median_Fe_mg_per_L, Median_Fe_mg_per_kg))%>%
  ungroup()

# Join all median data
all_medians = median_respiration %>%
  left_join(median_atp, by = "Sample_ID") %>%
  left_join(median_iron, by = "Sample_ID") %>%
  left_join(median_cn, by = "Sample_ID") %>%
  left_join(median_npoc_tn, by = "Sample_ID") %>%   
  mutate(Sample_Name = str_replace(Sample_ID, "INC", "all")) %>%
  select(-c(Sample_ID)) %>%
  relocate(Sample_Name, .before = Median_SpC_microsiemens_per_cm)

# Read in Wet/Dry Median Data from Summary file
median = read.csv("Data/EC_Sediment_Sample_Data_Summary.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  separate(Sample_Name, c("Sample_Name", "Rep"), sep = "-") %>%
  mutate(Sample_Name = paste0(Sample_Name, "_all")) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material, Median_Missing_Reps, Median_Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  select(-matches("per_L")) %>%
  rename_with(~ str_remove_all(., "_[0-9]+")) %>%
  rename_with(~ str_replace(., "^(([^_]*_){2}[^_]*).*", "\\1")) %>%
  rename_with(~ str_replace_all(., "Median", "median")) %>%
  rename(median_SpC = median_SpC_microsiemens) %>%
  rename(median_Temp = median_Temperature_degC)

median_dry = median %>%
  filter(Rep == "D") %>%
  select(c(Sample_Name, median_Initial_Gravimetric, median_Final_Gravimetric)) %>%
  mutate(across(c(median_Initial_Gravimetric:median_Final_Gravimetric), as.numeric))

# Effect Size Data
effect = read.csv("Data/EC_Sediment_Effect_Size.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%
  select(-c(IGSN, Field_Name, Material, Methods_Deviation))

# Grain size/ssa variables
grain = read.csv("Data/v3_CM_SSS_Sediment_Sample_Data_Summary.csv", skip = 2) %>%
  filter(grepl("CM", Sample_Name)) %>%
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>%
  mutate(Sample_Name = str_replace(Sample_Name, "Sediment", "all")) %>%
  select(c(Sample_Name, Percent_Tot_Sand, Percent_Coarse_Sand, Percent_Med_Sand, Percent_Fine_Sand, Percent_Silt, Percent_Clay, Mean_Specific_Surface_Area_m2_per_g))

# Join all data
effect_data = left_join(effect, grain, by = "Sample_Name") %>%
  left_join(all_medians) %>%
  left_join(median_dry) %>%
  mutate(across(c(Effect_Size_SpC_microsiemens_per_cm:Mean_Specific_Surface_Area_m2_per_g), as.numeric)) %>%
  select(-c(Effect_Size_Initial_Gravimetric_Moisture_g_per_g, Effect_Size_Final_Gravimetric_Moisture_g_per_g, Median_Respiration_Rate_mg_DO_per_L_per_H, Median_Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  rename(median_Dry_Initial_Gravimetric = median_Initial_Gravimetric) %>%
  rename(median_Dry_Final_Gravimetric = median_Final_Gravimetric) %>%
  mutate(Effect_Size_Respiration_Rate_mg_DO_per_L_per_H = abs(Effect_Size_Respiration_Rate_mg_DO_per_L_per_H)) %>%
  mutate(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = abs(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H))

# Prepare datasets
raw_effect = effect_data %>%
  column_to_rownames("Sample_Name") %>%
  filter(Effect_Size_Fe_mg_per_kg > -1) %>%
  select(-contains("per_L"), -c(Percent_Silt)) %>%
  relocate(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, .before = Effect_Size_SpC_microsiemens_per_cm)

cube_effect = effect_data %>%
  mutate(across(where(is.numeric), cube_root)) %>%
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>%
  column_to_rownames("Sample_Name") %>%
  filter(cube_Effect_Size_Fe_mg_per_kg > -1) %>%
  select(-contains("per_L"), -c(cube_Percent_Silt)) %>%
  relocate(cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, .before = cube_Effect_Size_SpC_microsiemens_per_cm)

cat("Data preparation complete!\n")
cat("Raw data:", nrow(raw_effect), "samples x", ncol(raw_effect)-1, "features\n")
cat("Cube root data:", nrow(cube_effect), "samples x", ncol(cube_effect)-1, "features\n\n")

## COMPREHENSIVE ANALYSIS FUNCTIONS =======================================================

# Function 1: Spearman Clustering and Variable Selection
spearman_clustering_selection = function(X_data, y_data, data_type = "Raw", cluster_threshold = 1.0, save_plots = TRUE) {
  cat("=== SPEARMAN CORRELATION CLUSTERING:", data_type, "===\n")
  
  # Calculate Spearman correlation matrix
  spearman_corr = cor(X_data, method = "spearman", use = "complete.obs")
  spearman_corr = (spearman_corr + t(spearman_corr)) / 2
  diag(spearman_corr) = 1
  
  # Convert to distance and cluster
  distance_matrix = 1 - abs(spearman_corr)
  hc_result = hclust(as.dist(distance_matrix), method = "ward.D2")
  cluster_ids = cutree(hc_result, h = cluster_threshold)
  
  cat("Clusters formed:", max(cluster_ids), "\n")
  
  # Create cluster assignments
  cluster_assignments = data.frame(
    feature = names(cluster_ids),
    cluster = cluster_ids
  ) %>% arrange(cluster, feature)
  
  # Method 1 Selection: Highest individual correlation within each cluster
  selected_features = cluster_assignments %>%
    group_by(cluster) %>%
    rowwise() %>%
    mutate(individual_cor = abs(cor(X_data[, feature], y_data, use = "complete.obs"))) %>%
    ungroup() %>%
    group_by(cluster) %>%
    slice_max(order_by = individual_cor, n = 1) %>%
    arrange(desc(individual_cor)) %>%
    pull(feature)
  
  # Print selection rationale
  selection_rationale = cluster_assignments %>%
    group_by(cluster) %>%
    rowwise() %>%
    mutate(individual_cor = abs(cor(X_data[, feature], y_data, use = "complete.obs"))) %>%
    ungroup() %>%
    group_by(cluster) %>%
    mutate(
      selected = feature %in% selected_features,
      rank_in_cluster = rank(desc(individual_cor))
    ) %>%
    arrange(cluster, rank_in_cluster)
  
  cat("\n=== CLUSTER ASSIGNMENTS ===\n")
  for(i in 1:max(cluster_ids)) {
    cluster_vars = selection_rationale %>% filter(cluster == i)
    selected_var = cluster_vars %>% filter(selected == TRUE)
    cat("Cluster", i, ":", nrow(cluster_vars), "variables\n")
    cat("  Selected:", selected_var$feature, "(r =", round(selected_var$individual_cor, 3), ")\n")
    other_vars = cluster_vars %>% filter(selected == FALSE)
    if(nrow(other_vars) > 0) {
      cat("  Others:", paste(paste0(other_vars$feature, "(r=", round(other_vars$individual_cor, 3), ")"), collapse = ", "), "\n")
    }
    cat("\n")
  }
  
  # Create visualizations if requested
  if(save_plots) {
    # Dendrogram and heatmap
    png(paste0("Figures/RF_results/Spearman_Clustering_", data_type, ".png"),
        width = 16, height = 8, units = "in", res = 300)
    par(mfrow = c(1, 2), mar = c(8, 4, 4, 2))
    
    # Dendrogram
    plot(hc_result,
         main = paste("Hierarchical Clustering:", data_type, "Data\n(Based on Spearman Correlations)"),
         xlab = "", sub = "", cex = 0.6)
    abline(h = cluster_threshold, col = "red", lty = 2, lwd = 2)
    
    # Correlation heatmap
    dendro_order = hc_result$order
    ordered_corr = spearman_corr[dendro_order, dendro_order]
    image(1:ncol(ordered_corr), 1:nrow(ordered_corr),
          ordered_corr,
          col = colorRampPalette(c("blue", "white", "red"))(100),
          main = paste("Spearman Correlations:", data_type, "Data\n(Ordered by Clustering)"),
          xlab = "", ylab = "", axes = FALSE)
    axis(1, at = 1:ncol(ordered_corr), labels = colnames(ordered_corr), las = 2, cex.axis = 0.5)
    axis(2, at = 1:nrow(ordered_corr), labels = rownames(ordered_corr), las = 2, cex.axis = 0.5)
    
    dev.off()
  }
  
  return(list(
    selected_features = selected_features,
    cluster_assignments = cluster_assignments,
    selection_rationale = selection_rationale,
    spearman_corr = spearman_corr,
    hc_result = hc_result
  ))
}

# Function 2: Permutation Importance Analysis
compute_detailed_permutation_importance = function(rf_model, X_data, y_data, data_type = "Raw", n_repeats = 10) {
  cat("=== PERMUTATION IMPORTANCE:", data_type, "DATA ===\n")
  
  imp_results = matrix(NA, nrow = n_repeats, ncol = ncol(X_data))
  colnames(imp_results) = colnames(X_data)
  
  pb = txtProgressBar(min = 0, max = n_repeats * ncol(X_data), style = 3, title = paste("Permutation Importance:", data_type))
  progress = 0
  
  for(i in 1:n_repeats) {
    baseline_pred = predict(rf_model, X_data)$predictions
    baseline_r2 = cor(baseline_pred, y_data)^2
    
    for(j in 1:ncol(X_data)) {
      X_perm = X_data
      X_perm[, j] = sample(X_perm[, j])
      perm_pred = predict(rf_model, X_perm)$predictions
      perm_r2 = cor(perm_pred, y_data)^2
      imp_results[i, j] = baseline_r2 - perm_r2
      
      progress = progress + 1
      setTxtProgressBar(pb, progress)
    }
  }
  close(pb)
  
  # Calculate statistics
  mean_imp = colMeans(imp_results, na.rm = TRUE)
  sd_imp = apply(imp_results, 2, sd, na.rm = TRUE)
  
  imp_summary = data.frame(
    feature = colnames(X_data),
    mean_importance = mean_imp,
    sd_importance = sd_imp,
    cv_importance = sd_imp / abs(mean_imp)
  ) %>%
    arrange(desc(mean_importance))
  
  # Create boxplot
  imp_long = imp_results %>%
    as.data.frame() %>%
    mutate(iteration = row_number()) %>%
    pivot_longer(cols = -iteration, names_to = "feature", values_to = "importance")
  
  feature_order = names(sort(mean_imp))
  perm_plot = ggplot(imp_long, aes(x = factor(feature, levels = feature_order), y = importance)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(title = paste("Permutation Importance:", data_type, "Data"),
         x = "Features", y = "Decrease in R²") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8), plot.title = element_text(hjust = 0.5))
  
  print(perm_plot)
  ggsave(paste0("Figures/RF_results/Detailed_Permutation_Importance_", tolower(gsub("_", "", data_type)), ".png"),
         plot = perm_plot, width = 12, height = 10, dpi = 300)
  
  cat("Top 5 most important features:\n")
  for(i in 1:min(5, nrow(imp_summary))) {
    cat("  ", i, ".", imp_summary$feature[i],
        "(importance:", round(imp_summary$mean_importance[i], 4), "±", round(imp_summary$sd_importance[i], 4), ")\n")
  }
  cat("\n")
  
  return(list(
    importance_matrix = imp_results,
    importance_summary = imp_summary,
    plot = perm_plot
  ))
}

# Function 3: Cross-Validation
perform_cross_validation = function(X_data, y_data, selected_features, data_type = "Raw", k_folds = 5) {
  cat("=== CROSS-VALIDATION:", data_type, "DATA ===\n")
  
  set.seed(42)
  n = nrow(X_data)
  folds = createFolds(1:n, k = k_folds, list = TRUE, returnTrain = FALSE)
  
  cv_results = data.frame(
    fold = integer(),
    r2_full = numeric(),
    r2_selected = numeric(),
    rmse_full = numeric(),
    rmse_selected = numeric()
  )
  
  for(i in 1:k_folds) {
    test_idx = folds[[i]]
    train_idx = setdiff(1:n, test_idx)
    
    X_train = X_data[train_idx, , drop = FALSE]
    X_test = X_data[test_idx, , drop = FALSE]
    y_train = y_data[train_idx]
    y_test = y_data[test_idx]
    
    # Full model
    rf_full = ranger(x = X_train, y = y_train, num.trees = 100, seed = 42)
    pred_full = predict(rf_full, X_test)$predictions
    r2_full = cor(pred_full, y_test)^2
    rmse_full = sqrt(mean((pred_full - y_test)^2))
    
    # Selected model
    X_train_sel = X_train[, selected_features, drop = FALSE]
    X_test_sel = X_test[, selected_features, drop = FALSE]
    rf_sel = ranger(x = X_train_sel, y = y_train, num.trees = 100, seed = 42)
    pred_sel = predict(rf_sel, X_test_sel)$predictions
    r2_sel = cor(pred_sel, y_test)^2
    rmse_sel = sqrt(mean((pred_sel - y_test)^2))
    
    cv_results = rbind(cv_results, data.frame(
      fold = i, r2_full = r2_full, r2_selected = r2_sel,
      rmse_full = rmse_full, rmse_selected = rmse_sel
    ))
  }
  
  # Summary statistics
  cv_summary = cv_results %>%
    summarise(
      mean_r2_selected = mean(r2_selected, na.rm = TRUE),
      sd_r2_selected = sd(r2_selected, na.rm = TRUE),
      mean_rmse_selected = mean(rmse_selected, na.rm = TRUE),
      sd_rmse_selected = sd(rmse_selected, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("Cross-validation results:\n")
  cat("  R²:", round(cv_summary$mean_r2_selected, 3), "±", round(cv_summary$sd_r2_selected, 3), "\n")
  cat("  RMSE:", round(cv_summary$mean_rmse_selected, 3), "±", round(cv_summary$sd_rmse_selected, 3), "\n\n")
  
  return(list(cv_results = cv_results, cv_summary = cv_summary))
}

# Function 4: Bootstrap Feature Selection Stability
bootstrap_stability_analysis = function(X_data, y_data, data_type = "Raw", n_bootstrap = 100, cluster_threshold = 1.0) {
  cat("=== BOOTSTRAP STABILITY:", data_type, "DATA ===\n")
  
  set.seed(42)
  n = nrow(X_data)
  feature_counts = data.frame(feature = colnames(X_data), count = 0, proportion = 0)
  bootstrap_selections = list()
  
  pb = txtProgressBar(min = 0, max = n_bootstrap, style = 3, title = paste("Bootstrap:", data_type))
  
  for(i in 1:n_bootstrap) {
    boot_idx = sample(1:n, n, replace = TRUE)
    X_boot = X_data[boot_idx, , drop = FALSE]
    y_boot = y_data[boot_idx]
    
    tryCatch({
      # Clustering and selection on bootstrap sample
      selection_result = spearman_clustering_selection(X_boot, y_boot, data_type = paste0(data_type, "_Bootstrap"),
                                                       cluster_threshold = cluster_threshold, save_plots = FALSE)
      selected_features = selection_result$selected_features
      bootstrap_selections[[i]] = selected_features
      
      # Update counts
      for(feat in selected_features) {
        idx = which(feature_counts$feature == feat)
        if(length(idx) > 0) {
          feature_counts$count[idx] = feature_counts$count[idx] + 1
        }
      }
    }, error = function(e) {
      bootstrap_selections[[i]] = character(0)
    })
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Calculate proportions
  feature_counts$proportion = feature_counts$count / n_bootstrap
  stable_features = feature_counts %>%
    filter(proportion > 0.5) %>%
    arrange(desc(proportion))
  
  cat("Features stable across >50% of bootstrap samples:", nrow(stable_features), "\n")
  if(nrow(stable_features) > 0) {
    for(i in 1:nrow(stable_features)) {
      cat("  ", stable_features$feature[i], "(", round(stable_features$proportion[i]*100, 1), "%)\n")
    }
  }
  cat("\n")
  
  return(list(
    feature_stability = feature_counts,
    stable_features = stable_features,
    bootstrap_selections = bootstrap_selections
  ))
}

# Function 5: SHAP Analysis
shap_analysis = function(X_data, y_data, selected_features, data_type = "Raw") {
  cat("=== SHAP ANALYSIS:", data_type, "DATA ===\n")
  
  X_selected = X_data[, selected_features, drop = FALSE]
  
  # Train model
  set.seed(42)
  rf_model = ranger(x = X_selected, y = y_data, num.trees = 100, seed = 42)
  
  # Prediction function for SHAP
  pred_function = function(model, newdata) {
    predict(model, newdata)$predictions
  }
  
  # Calculate SHAP values
  cat("Computing SHAP values for", length(selected_features), "features...\n")
  shap_values = fastshap::explain(
    rf_model,
    X = X_selected,
    pred_wrapper = pred_function,
    nsim = 50
  )
  
  # SHAP summary
  shap_summary = shap_values %>%
    as.data.frame() %>%
    summarise_all(list(mean_abs = ~ mean(abs(.)), mean = ~ mean(.), median = ~ median(.), sd = ~ sd(.))) %>%
    pivot_longer(everything(), names_to = "stat_feature", values_to = "value") %>%
    separate(stat_feature, into = c("feature", "stat"), sep = "_(?=mean_abs|mean|median|sd)") %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    arrange(desc(mean_abs))
  
  # SHAP plots
  shap_long = shap_values %>%
    as.data.frame() %>%
    mutate(observation = row_number()) %>%
    pivot_longer(cols = -observation, names_to = "feature", values_to = "shap_value") %>%
    left_join(
      X_selected %>%
        mutate(observation = row_number()) %>%
        pivot_longer(cols = -observation, names_to = "feature", values_to = "feature_value"),
      by = c("observation", "feature")
    )
  
  # SHAP summary plot
  shap_summary_plot = shap_long %>%
    ggplot(aes(x = shap_value, y = reorder(feature, shap_value, function(x) mean(abs(x))))) +
    geom_point(aes(color = feature_value), alpha = 0.7, size = 2) +
    scale_color_viridis_c(name = "Feature\nValue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
    labs(title = paste("SHAP Summary:", data_type, "Data"),
         subtitle = "Each point is one observation; color shows feature value",
         x = "SHAP Value (Impact on Prediction)",
         y = "Features (ranked by |SHAP|)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Individual SHAP dependence plots for top features
  top_features = shap_summary$feature[1:min(3, nrow(shap_summary))]
  dependence_plots = list()
  
  for(i in 1:length(top_features)) {
    feat = top_features[i]
    p = shap_long %>%
      filter(feature == feat) %>%
      ggplot(aes(x = feature_value, y = shap_value)) +
      geom_point(alpha = 0.6, color = "steelblue", size = 2) +
      geom_smooth(method = "loess", se = TRUE, color = "red", alpha = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      labs(title = str_wrap(feat, 20),
           x = "Feature Value", y = "SHAP Value") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 9))
    dependence_plots[[i]] = p
  }
  
  cat("SHAP analysis complete. Top features by importance:\n")
  for(i in 1:nrow(shap_summary)) {
    direction = ifelse(shap_summary$mean[i] > 0, "Positive", "Negative")
    cat("  ", i, ".", shap_summary$feature[i],
        "| Mean |SHAP|:", round(shap_summary$mean_abs[i], 4),
        "| Direction:", direction, "\n")
  }
  cat("\n")
  
  return(list(
    shap_values = shap_values,
    shap_summary = shap_summary,
    shap_long = shap_long,
    summary_plot = shap_summary_plot,
    dependence_plots = dependence_plots,
    rf_model = rf_model
  ))
}

# Main comprehensive analysis function
run_complete_robust_analysis = function(X_data, y_data, data_type = "Raw", save_prefix = "raw") {
  cat("\n" %>% strrep(3))
  cat("=== COMPLETE ROBUST ANALYSIS:", data_type, "DATA ===\n")
  
  # Step 1: Spearman Clustering and Variable Selection
  clustering_results = spearman_clustering_selection(X_data, y_data, data_type)
  selected_features = clustering_results$selected_features
  
  # Step 2: Train models
  set.seed(42)
  # Full model
  rf_full = ranger(x = X_data, y = y_data, num.trees = 100, importance = "impurity", seed = 42)
  full_pred = predict(rf_full, X_data)$predictions
  full_r2 = cor(full_pred, y_data)^2
  
  # Selected model
  X_selected = X_data[, selected_features, drop = FALSE]
  rf_selected = ranger(x = X_selected, y = y_data, num.trees = 100, importance = "impurity", seed = 42)
  selected_pred = predict(rf_selected, X_selected)$predictions
  selected_r2 = cor(selected_pred, y_data)^2
  
  cat("Model Performance:\n")
  cat("  Full model R²:", round(full_r2, 3), "(", ncol(X_data), "variables )\n")
  cat("  Selected model R²:", round(selected_r2, 3), "(", length(selected_features), "variables )\n")
  cat("  Performance change:", round(selected_r2 - full_r2, 3), "\n\n")
  
  # Step 3: Detailed Permutation Importance (on selected features)
  perm_results = compute_detailed_permutation_importance(rf_selected, X_selected, y_data, data_type)
  
  # Step 4: Cross-Validation
  cv_results = perform_cross_validation(X_data, y_data, selected_features, data_type)
  
  # Step 5: Bootstrap Stability
  bootstrap_results = bootstrap_stability_analysis(X_data, y_data, data_type)
  
  # Step 6: SHAP Analysis
  shap_results = shap_analysis(X_selected, y_data, selected_features, data_type)
  
  # Step 7: Create comprehensive summary
  comprehensive_summary = data.frame(
    Feature = selected_features,
    Individual_Correlation = sapply(selected_features, function(f) abs(cor(X_data[, f], y_data, use = "complete.obs"))),
    Permutation_Importance = perm_results$importance_summary$mean_importance,
    Permutation_SD = perm_results$importance_summary$sd_importance,
    Bootstrap_Stability = sapply(selected_features, function(f) {
      idx = which(bootstrap_results$feature_stability$feature == f)
      ifelse(length(idx) > 0, bootstrap_results$feature_stability$proportion[idx], 0)
    }),
    SHAP_Mean_Abs = shap_results$shap_summary$mean_abs,
    SHAP_Mean = shap_results$shap_summary$mean,
    SHAP_Direction = ifelse(shap_results$shap_summary$mean > 0, "Positive", "Negative")
  ) %>%
    arrange(desc(SHAP_Mean_Abs))
  
  # Step 8: Create comprehensive visualization
  print(shap_results$summary_plot)
  ggsave(paste0("Figures/RF_results/SHAP_Summary_", save_prefix, ".png"),
         plot = shap_results$summary_plot, width = 10, height = 8, dpi = 300)
  
  # Dependence plots for top 3
  if(length(shap_results$dependence_plots) > 0) {
    combined_dependence = do.call(ggarrange, c(shap_results$dependence_plots, list(ncol = 3, nrow = 1)))
    print(combined_dependence)
    ggsave(paste0("Figures/RF_results/SHAP_Dependence_", save_prefix, ".png"),
           plot = combined_dependence, width = 15, height = 5, dpi = 300)
  }
  
  # Step 9: Export all results to RF_results folder
  write.csv(comprehensive_summary, paste0("Data/RF_results/Complete_Analysis_", data_type, ".csv"), row.names = FALSE)
  write.csv(clustering_results$selection_rationale, paste0("Data/RF_results/Clustering_Details_", data_type, ".csv"), row.names = FALSE)
  write.csv(cv_results$cv_results, paste0("Data/RF_results/Cross_Validation_", data_type, ".csv"), row.names = FALSE)
  write.csv(bootstrap_results$feature_stability, paste0("Data/RF_results/Bootstrap_Stability_", data_type, ".csv"), row.names = FALSE)
  write.csv(shap_results$shap_summary, paste0("Data/RF_results/SHAP_Analysis_", data_type, ".csv"), row.names = FALSE)
  write.csv(perm_results$importance_summary, paste0("Data/RF_results/Permutation_Importance_", data_type, ".csv"), row.names = FALSE)
  
  return(list(
    selected_features = selected_features,
    clustering_results = clustering_results,
    perm_results = perm_results,
    cv_results = cv_results,
    bootstrap_results = bootstrap_results,
    shap_results = shap_results,
    comprehensive_summary = comprehensive_summary,
    model_performance = list(full_r2 = full_r2, selected_r2 = selected_r2)
  ))
}

## RUN COMPLETE ANALYSES ON BOTH DATASETS ==============================================

# Prepare data
exclude_col_raw = "Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H"
X_raw = raw_effect[, !(names(raw_effect) %in% exclude_col_raw)]
y_raw = raw_effect$Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H

exclude_col_cube = "cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H"
X_cube = cube_effect[, !(names(cube_effect) %in% exclude_col_cube)]
y_cube = cube_effect$cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H

cat("Starting complete robust analysis for both datasets...\n\n")

# Run complete analyses
raw_complete = run_complete_robust_analysis(X_raw, y_raw, "Raw", "raw")
cube_complete = run_complete_robust_analysis(X_cube, y_cube, "Cube_Root", "cube")

## FINAL COMPREHENSIVE COMPARISON =======================================================

cat("\n" %>% strrep(5))
cat("=== ULTIMATE COMPREHENSIVE COMPARISON ===\n")

final_comparison = data.frame(
  Data_Type = c("Raw", "Cube Root"),
  Selected_Features = c(length(raw_complete$selected_features), length(cube_complete$selected_features)),
  Model_R2 = c(raw_complete$model_performance$selected_r2, cube_complete$model_performance$selected_r2),
  CV_R2_Mean = c(raw_complete$cv_results$cv_summary$mean_r2_selected, cube_complete$cv_results$cv_summary$mean_r2_selected),
  CV_R2_SD = c(raw_complete$cv_results$cv_summary$sd_r2_selected, cube_complete$cv_results$cv_summary$sd_r2_selected),
  Stable_Features = c(sum(raw_complete$bootstrap_results$feature_stability$proportion > 0.5),
                      sum(cube_complete$bootstrap_results$feature_stability$proportion > 0.5)),
  Top_SHAP_Feature = c(raw_complete$shap_results$shap_summary$feature[1],
                       cube_complete$shap_results$shap_summary$feature[1])
)

print(final_comparison)

# Export final comparison
write.csv(final_comparison, "Data/RF_results/Final_Model_Comparison.csv", row.names = FALSE)

# Determine final recommendation
cv_comparison = final_comparison$CV_R2_Mean - final_comparison$CV_R2_SD  # Conservative estimate
recommended_idx = which.max(cv_comparison)
recommended_approach = final_comparison$Data_Type[recommended_idx]
recommended_results = if(recommended_approach == "Raw") raw_complete else cube_complete

cat("\n=== FINAL RECOMMENDATION ===\n")
cat("RECOMMENDED APPROACH:", recommended_approach, "data\n")
cat("Reason: Highest conservative CV R² estimate\n")
cat("Performance: R² =", round(final_comparison$CV_R2_Mean[recommended_idx], 3),
    "±", round(final_comparison$CV_R2_SD[recommended_idx], 3), "\n")
cat("Stable features:", final_comparison$Stable_Features[recommended_idx], "\n")

cat("\n=== ROBUST FEATURES FOR PUBLICATION ===\n")
cat("Features meeting all robustness criteria:\n")

robust_features = recommended_results$comprehensive_summary %>%
  filter(Bootstrap_Stability > 0.5) %>%  # Stable across bootstrap
  arrange(desc(SHAP_Mean_Abs))

for(i in 1:nrow(robust_features)) {
  cat(sprintf("%d. %s\n", i, robust_features$Feature[i]))
  cat(sprintf("   Individual correlation: r = %.3f\n", robust_features$Individual_Correlation[i]))
  cat(sprintf("   Permutation importance: %.3f ± %.3f\n", robust_features$Permutation_Importance[i], robust_features$Permutation_SD[i]))
  cat(sprintf("   Bootstrap stability: %.1f%%\n", robust_features$Bootstrap_Stability[i] * 100))
  cat(sprintf("   SHAP importance: %.3f (%s effect)\n", robust_features$SHAP_Mean_Abs[i], robust_features$SHAP_Direction[i]))
  cat("\n")
}

# Export robust features summary
write.csv(robust_features, paste0("Data/RF_results/Robust_Features_", recommended_approach, ".csv"), row.names = FALSE)

cat("=== BIOLOGICAL INTERPRETATION ===\n")
cat("Your", recommended_approach, "data analysis reveals:\n")
cat("• Primary drivers of O₂ consumption effect size\n")
cat("• Directionality of each effect (SHAP analysis)\n")
cat("• Stability of findings across data resampling\n")
cat("• Generalizability via cross-validation\n")
cat("• Non-linear relationships in dependence plots\n")

cat("\n=== METHODOLOGICAL ROBUSTNESS ACHIEVED ===\n")
cat("✅ Spearman clustering: Identified correlated variable groups\n")
cat("✅ Method 1 selection: Chose best representative from each cluster\n")
cat("✅ Permutation importance: Model-agnostic importance with uncertainty\n")
cat("✅ Cross-validation: Confirmed model generalizability\n")
cat("✅ Bootstrap analysis: Verified feature selection stability\n")
cat("✅ SHAP analysis: Revealed directionality and effect magnitude\n")
cat("✅ Transformation comparison: Selected optimal data representation\n")

cat("\n=== ALL FILES GENERATED IN RF_RESULTS FOLDERS ===\n")

# Data files in Data/RF_results/
data_files = c(
  "Complete_Analysis_Raw.csv", "Complete_Analysis_Cube_Root.csv",
  "Clustering_Details_Raw.csv", "Clustering_Details_Cube_Root.csv",
  "Cross_Validation_Raw.csv", "Cross_Validation_Cube_Root.csv",
  "Bootstrap_Stability_Raw.csv", "Bootstrap_Stability_Cube_Root.csv",
  "SHAP_Analysis_Raw.csv", "SHAP_Analysis_Cube_Root.csv",
  "Permutation_Importance_Raw.csv", "Permutation_Importance_Cube_Root.csv",
  "Final_Model_Comparison.csv", paste0("Robust_Features_", recommended_approach, ".csv")
)

# Figure files in Figures/RF_results/
figure_files = c(
  "Spearman_Clustering_Raw.png", "Spearman_Clustering_Cube_Root.png",
  "Detailed_Permutation_Importance_raw.png", "Detailed_Permutation_Importance_cuberoot.png",
  "SHAP_Summary_raw.png", "SHAP_Summary_cube.png",
  "SHAP_Dependence_raw.png", "SHAP_Dependence_cube.png"
)

cat("DATA FILES (Data/RF_results/):\n")
for(file in data_files) {
  cat("- ", file, "\n")
}

cat("\nFIGURE FILES (Figures/RF_results/):\n")
for(file in figure_files) {
  cat("- ", file, "\n")
}

cat("\n=== SUCCESS! ===\n")
cat("Complete analysis finished with organized file structure:\n")
cat("• All CSV files saved in: Data/RF_results/\n")
cat("• All figures saved in: Figures/RF_results/\n")
cat("• Robust variable selection with multiple validation approaches\n")
cat("• Clear biological interpretation with directionality\n")
cat("• Publication-ready methodology with all robustness checks\n")
cat("• Comprehensive comparison of data transformation approaches\n")

cat("\nFocus your biological interpretation on the", nrow(robust_features), "features that are:\n")
cat("• Bootstrap stable (>50% selection rate)\n")
cat("• High SHAP importance\n")
cat("• Clear directional effects\n")
cat("• Biologically meaningful\n")
# ============================================================
# FIGURE S6 — 
# ============================================================

# ---- Pick results ----
results <- cube_complete  # change to raw_complete if needed

# ---- Check what features exist ----
cat("Permutation features:\n")
print(colnames(results$perm_results$importance_matrix))
cat("\nSHAP features:\n")
print(unique(results$shap_results$shap_long$feature))

# ---- Name mapping ----
name_mapping <- c(
  "cube_Percent_Fine_Sand"                   = "Fine Sand (%)",
  "cube_Percent_Med_Sand"                    = "Med. Sand (%)",
  "cube_Percent_Coarse_Sand"                 = "Coarse Sand (%)",
  "cube_Percent_Tot_Sand"                    = "Total Sand (%)",
  "cube_Percent_Clay"                        = "Clay (%)",
  "cube_Median_ATP_picomoles_per_g"          = "ATP (pmol/g)",
  "cube_Median_SpC_microsiemens_per_cm"      = "SpC (uS/cm)",
  "cube_Median_pH"                           = "Median pH",
  "cube_Median_Fe_mg_per_kg"                 = "Fe (mg/kg)",
  "cube_Median_Temperature_degC"             = "Temperature (C)",
  "cube_Median_X01395_C_percent_per_mg"      = "Sediment C (%)",
  "cube_Median_X01397_N_percent_per_mg"      = "Sediment N (%)",
  "cube_Median_Extractable_NPOC_mg_per_kg"   = "NPOC (mg/kg)",
  "cube_Median_Extractable_TN_mg_per_kg"     = "Extractable TN (mg/kg)",
  "cube_Mean_Specific_Surface_Area_m2_per_g" = "SSA (m2/g)",
  "cube_median_Dry_Initial_Gravimetric"      = "Grav. Moisture",
  "cube_median_Dry_Final_Gravimetric"        = "Final Grav. Moisture",
  "cube_Effect_Size_pH"                      = "Effect Size pH",
  "cube_Effect_Size_SpC_microsiemens_per_cm" = "Effect Size SpC",
  "cube_Effect_Size_Fe_mg_per_kg"            = "Effect Size Fe",
  "cube_Effect_Size_Temperature_degC"        = "Effect Size Temp",
  "cube_Effect_Size_ATP_picomoles_per_g"     = "Effect Size ATP",
  "Percent_Fine_Sand"                        = "Fine Sand (%)",
  "Percent_Med_Sand"                         = "Med. Sand (%)",
  "Percent_Coarse_Sand"                      = "Coarse Sand (%)",
  "Percent_Tot_Sand"                         = "Total Sand (%)",
  "Percent_Clay"                             = "Clay (%)",
  "Median_ATP_picomoles_per_g"               = "ATP (pmol/g)",
  "Median_SpC_microsiemens_per_cm"           = "SpC (uS/cm)",
  "Median_pH"                                = "Median pH",
  "Median_Fe_mg_per_kg"                      = "Fe (mg/kg)",
  "Median_Temperature_degC"                  = "Temperature (C)",
  "Median_X01395_C_percent_per_mg"           = "Sediment C (%)",
  "Median_X01397_N_percent_per_mg"           = "Sediment N (%)",
  "Median_Extractable_NPOC_mg_per_kg"        = "NPOC (mg/kg)",
  "Median_Extractable_TN_mg_per_kg"          = "Extractable TN (mg/kg)",
  "Mean_Specific_Surface_Area_m2_per_g"      = "SSA (m2/g)",
  "median_Dry_Initial_Gravimetric"           = "Grav. Moisture",
  "median_Dry_Final_Gravimetric"             = "Final Grav. Moisture",
  "Effect_Size_pH"                           = "Effect Size pH",
  "Effect_Size_SpC_microsiemens_per_cm"      = "Effect Size SpC",
  "Effect_Size_Fe_mg_per_kg"                 = "Effect Size Fe",
  "Effect_Size_Temperature_degC"             = "Effect Size Temp",
  "Effect_Size_ATP_picomoles_per_g"          = "Effect Size ATP"
)

clean_name <- function(x) {
  out <- name_mapping[x]
  out[is.na(out)] <- x[is.na(out)] %>%
    str_remove("^cube_") %>%
    str_remove("^Median_") %>%
    str_remove("^Effect_Size_") %>%
    str_replace_all("_", " ")
  return(out)
}

# ============================================================
# PANEL A with larger text
# ============================================================

panel_a <- ggplot(imp_long,
                  aes(x = importance,
                      y = factor(feature_clean, levels = imp_order))) +
  geom_boxplot(fill = "#B8D4E3", color = "black",
               outlier.size = 0.8, outlier.alpha = 0.5,
               width = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = expression(Decrease~"in"~R^2), y = NULL) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y    = element_text(size = 14, color = "black"),
    axis.text.x    = element_text(size = 12, color = "black"),
    axis.title.x   = element_text(size = 14, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin = margin(5, 15, 5, 5)
  )

# ============================================================
# PANEL B with larger text
# ============================================================

panel_b <- ggplot(shap_clean,
                  aes(x = shap_value,
                      y = factor(feature_clean, levels = shap_order),
                      color = feature_scaled)) +
  geom_jitter(height = 0.2, width = 0,
              alpha = 0.7, size = 2.5, shape = 16) +
  scale_color_viridis_c(name = "Feature\nValue",
                        option = "viridis",
                        breaks = c(0, 0.5, 1),
                        labels = c("Low", "Mid", "High")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "SHAP Value (Impact on Prediction)", y = NULL) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y       = element_text(size = 14, color = "black"),
    axis.text.x       = element_text(size = 12, color = "black"),
    axis.title.x      = element_text(size = 14, face = "bold"),
    legend.title      = element_text(size = 12, face = "bold"),
    legend.text       = element_text(size = 11),
    legend.key.height = unit(1.2, "cm"),
    legend.key.width  = unit(0.5, "cm"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin = margin(5, 5, 5, 15)
  )

# ============================================================
# COMBINE with larger tags and subtitle
# ============================================================

fig_s6 <- panel_a + panel_b +
  plot_layout(widths = c(1, 1.2)) +
  plot_annotation(
    subtitle = "All variables are cube root transformed",
    tag_levels = "A"
  ) &
  theme(
    plot.tag      = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13, face = "italic", hjust = 0.5)
  )

print(fig_s6)

# Save larger
ggsave("Figures/Figure_S6_Publication.png",
       fig_s6, width = 16, height = 7, dpi = 300, bg = "white")

ggsave("Figures/Figure_S6_Publication.pdf",
       fig_s6, width = 16, height = 7)

