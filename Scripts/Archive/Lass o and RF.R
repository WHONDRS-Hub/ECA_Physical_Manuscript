#### Sensitivity Analysis For ECA removals with LASSO and Random Forest ####
# This script makes figures for ECA physical manuscript and calculates LASSO regression
# followed by Random Forest with permutation importance and stability assessment

library(tidyverse);library(corrplot);library(ggpubr);library(ggpmisc);library(factoextra);library(stringr);library(glmnet);library(magick); library(ggnewscale); library(FSA); library(multcompView); library(rcompanion)

# Load additional libraries for RF
if(!require(ranger, quietly = TRUE)) install.packages("ranger")
if(!require(patchwork, quietly = TRUE)) install.packages("patchwork")
if(!require(ggrepel, quietly = TRUE)) install.packages("ggrepel")
if(!require(viridis, quietly = TRUE)) install.packages("viridis")
library(ranger); library(patchwork); library(ggrepel); library(viridis)

rm(list=ls());graphics.off()

# Functions ---------------------------------------------------------------
# Transformation for normalization is cube root - have to cube root then add sign back to value to make it positive or negative
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

## Read in/Merge Data  ------------------------------------------------------------
# Individual Rate data for histograms ####
all_data = read.csv("Data/EC_Sediment_SpC_pH_Temp_Respiration.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  # remove samples with too much water (EC_011, EC_012), sample with no mg/kg (EC_023), duplicated NEON sites (EC_052, EC_053, EC_057)
  select(-c(Field_Name, IGSN, Material))

# Take Cube root of all respiration rates for figure
cube_respiration = all_data %>%
  select(c(Sample_Name, Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999) %>% # removes overexposed samples, missing replicates
  mutate(cube_Respiration_mg_kg = cube_root(abs(as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)))) %>% # make respiration positive and cube root
  mutate(Treat = if_else(grepl("D", Sample_Name), "Dry", "Wet"))

# Calculate "bulk" medians (not separated by wet/dry)
median_respiration = all_data %>%
  select(-c(Respiration_R_Squared, Respiration_R_Squared_Adj, Respiration_p_value, Total_Incubation_Time_Min, Number_Points_In_Respiration_Regression, Number_Points_Removed_Respiration_Regression,DO_Concentration_At_Incubation_Time_Zero)) %>%
  mutate(across(c(SpC_microsiemens_per_cm:Respiration_Rate_mg_DO_per_kg_per_H), as.numeric)) %>%
  mutate(Respiration_Rate_mg_DO_per_L_per_H = ifelse(grepl("INC_Method_001|INC_Method_002|INC_QA_004", Methods_Deviation), NA, Respiration_Rate_mg_DO_per_L_per_H)) %>%
  #missing replicates (EC_072-W5/D5),  overexposed samples (EC_027, EC_013, EC_014), less sediment in sample (EC_012-D5)
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(grepl("INC_Method_001|INC_Method_002|INC_QA_004", Methods_Deviation), NA, Respiration_Rate_mg_DO_per_kg_per_H)) %>% #missing replicates (EC_072-W5/D5),  overexposed samples (EC_027, EC_013, EC_014), less sediment in sample (EC_012-D5)
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
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 10), "FALSE", "TRUE")) %>% ## Check CV's, then remove
  select(c(Sample_ID, Median_SpC_microsiemens_per_cm, Median_pH, Median_Temperature_degC, Median_Respiration_Rate_mg_DO_per_L_per_H, Median_Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  ungroup()

## Calculate "bulk" medians (not separated by wet/dry) ---------------------
# ATP ####
atp = read.csv("Data/EC_Sediment_ATP.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material)) %>%
  mutate(Sample_Name = str_replace(Sample_Name, "ATP", "INC"))

median_atp = atp %>%
  mutate(ATP_nanomoles_per_L = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, ATP_nanomoles_per_L)) %>%
  # ATP_002 (don't have sample), INC_Method_001
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
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 10), "FALSE", "TRUE")) %>% ## Check CV's, then remove
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
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 10), "FALSE", "TRUE")) %>% ## Check CV's, then remove
  select(c(Sample_ID, Median_X01395_C_percent_per_mg, Median_X01397_N_percent_per_mg))%>%
  ungroup()

# NPOC/TN ####
npoc_tn = read.csv("Data/EC_Sediment_NPOC_TN.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material)) %>%
  mutate(Sample_Name = str_replace(Sample_Name, "SIR", "INC"))

# remove broken/missing samples
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
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 10), "FALSE", "TRUE")) %>% ## Check CV's, then remove
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

# Calculate "bulk" medians
# if samples flagged as Below LOD/lowest standard - give value of 1/2 of lowest LOD from all ECA Fe analysis (0.002/2 = 0.001)
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
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 5), "FALSE", "TRUE")) %>% ## Check CV's, then remove
  select(c(Sample_ID, Median_Fe_mg_per_L, Median_Fe_mg_per_kg))%>%
  ungroup()

# Join all median data ----------------------------------------------------
## All Bulk Medians
all_medians = median_respiration %>%
  left_join(median_atp, by = "Sample_ID") %>%
  left_join(median_iron, by = "Sample_ID") %>%
  left_join(median_cn, by = "Sample_ID") %>%
  left_join(median_npoc_tn, by = "Sample_ID") %>%   mutate(Sample_Name = str_replace(Sample_ID, "INC", "all")) %>%
  select(-c(Sample_ID)) %>%
  relocate(Sample_Name, .before = Median_SpC_microsiemens_per_cm)

## Read in Wet/Dry Median Data from Summary file to get Dry moisture values
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

# Median Dry Initial Gravimetric moisture used in analysis
median_dry = median %>%
  filter(Rep == "D") %>%
  select(c(Sample_Name, median_Initial_Gravimetric, median_Final_Gravimetric)) %>%
  mutate(across(c(median_Initial_Gravimetric:median_Final_Gravimetric), as.numeric))

## Effect Size Data ####
effect = read.csv("Data/EC_Sediment_Effect_Size.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%
  select(-c(IGSN, Field_Name, Material, Methods_Deviation))

## Read in grain size/ssa variables ####
grain = read.csv("Data/v3_CM_SSS_Sediment_Sample_Data_Summary.csv", skip = 2) %>%
  filter(grepl("CM", Sample_Name)) %>%
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>%
  mutate(Sample_Name = str_replace(Sample_Name, "Sediment", "all")) %>%
  select(c(Sample_Name, Percent_Tot_Sand, Percent_Coarse_Sand, Percent_Med_Sand, Percent_Fine_Sand, Percent_Silt, Percent_Clay, Mean_Specific_Surface_Area_m2_per_g))

## Join all data with effect size ####
effect_data = left_join(effect, grain, by = "Sample_Name") %>%
  left_join(all_medians) %>%
  left_join(median_dry) %>%
  mutate(across(c(Effect_Size_SpC_microsiemens_per_cm:Mean_Specific_Surface_Area_m2_per_g), as.numeric)) %>%  # make data numeric
  select(-c(Effect_Size_Initial_Gravimetric_Moisture_g_per_g, Effect_Size_Final_Gravimetric_Moisture_g_per_g, Median_Respiration_Rate_mg_DO_per_L_per_H, Median_Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  rename(median_Dry_Initial_Gravimetric = median_Initial_Gravimetric) %>%
  rename(median_Dry_Final_Gravimetric = median_Final_Gravimetric) %>%
  mutate(Effect_Size_Respiration_Rate_mg_DO_per_L_per_H = abs(Effect_Size_Respiration_Rate_mg_DO_per_L_per_H)) %>%
  mutate(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = abs(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H))

# Transform Data ----------------------------------------------------------
# Fe outlier not included in analysis - remove from DF
cube_effect = effect_data %>%
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>%
  column_to_rownames("Sample_Name") %>%
  filter(cube_Effect_Size_Fe_mg_per_kg > -1) %>%  # remove Fe outlier for analysis
  select(-contains("per_L")) %>% # remove per_L data, analysis ran on _per_kg
  relocate(cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, .before = cube_Effect_Size_SpC_microsiemens_per_cm)

# Data frame that includes only Fe outlier - this used for Scatter Plot
fe_cube_effect = effect_data %>%
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>%
  column_to_rownames("Sample_Name") %>%
  select(-contains("per_L")) %>%
  filter(cube_Effect_Size_Fe_mg_per_kg < -1) %>%
  relocate(cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, .before = cube_Effect_Size_SpC_microsiemens_per_cm)

# Histograms --------------------------------------------------------------
## Histogram of all Rates
all_hist = all_data %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999) %>% # removes overexposed samples, missing replicates
  mutate(Treat = if_else(grepl("D", Sample_Name), "Dry", "Wet")) %>%
  ggplot(aes(x = abs(as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)))) +
  geom_histogram(position = "identity", alpha = 0.8, aes(fill = Treat))+
  scale_fill_manual(values = c("#D55E00","#0072B2"))  +
  theme(strip.text = element_text(size = 4))+
  theme_bw() +
  theme(legend.position = c(0.85, 0.8),
        legend.key.size = unit(0.15, "in"),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 10)) +
  guides(fill = guide_legend(title="Treatment")) +
  xlab(expression(atop("\n Respiration Rate", "(mg O"[2]*" kg"^-1*" h"^-1*")"))) +
  ylab("Count")

effect_limits <- c(-1400, 1400)
effect_hist = ggplot(effect_data, aes(x = Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H))+
  geom_histogram(binwidth = 25, aes(fill = after_stat(x))) +
  scale_fill_gradient2(name = "Effect Size", limits = effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(effect_limits)+min(effect_limits))/2) +
  theme_bw()+
  xlim(c(-1500, 1500))+
  ylab("Count\n")+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 10)) +
  xlab(expression(atop("\n Effect Size Respiration Rate","(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" h"^-1*")")))

combined_SI_hist = ggarrange(all_hist, effect_hist, ncol = 2, labels = c("A", "B"), hjust = -1, vjust = 2.5)
ggsave("Figures/Combined_SI_Histograms.png", plot = combined_SI_hist, width = 10, height = 5, dpi = 300)

## Cube Root Transformed
all_cube_hist = ggplot(cube_respiration, aes(x = cube_Respiration_mg_kg)) +
  geom_histogram(position = "identity", alpha = 0.8, aes(fill = Treat))+
  scale_fill_manual(values = c("#D55E00","#0072B2"))  +
  theme(strip.text = element_text(size = 4))+
  ylim(0, 87.5)+
  theme_bw() +
  theme(legend.position = c(0.85, 0.8),
        legend.key.size = unit(0.15, "in"),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 10)) +
  guides(fill = guide_legend(title="Treatment")) +
  xlab(expression(atop("\n Respiration Rate" ^(1/3)*"", "(mg O"[2]*" kg"^-1*" h"^-1*")"))) +
  ylab("Count")

cube_effect_limits <- c(-12, 12)
cube_effect_hist = ggplot(cube_effect, aes(x = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H))+
  geom_histogram(binwidth = 0.5, aes(fill = after_stat(x))) +
  scale_fill_gradient2(name = "Cubed Root Effect Size", limits = cube_effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(cube_effect_limits)+min(cube_effect_limits))/2) +
  theme_bw()+
  xlim(c(-12, 12))+
  ylab("Count\n")+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 10)) +
  xlab(expression(atop("\n Effect Size Respiration Rate" ^(1/3)*"","(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" h"^-1*")")))

combined_hist = ggarrange(all_cube_hist, cube_effect_hist, ncol = 2, labels = c("A", "B"), hjust = -5, vjust = 2.5)

## Pearson Correlation Matrix ####
# scale data before it goes into correlation matrix
scale_cube_effect = as.data.frame(scale(cube_effect, center = T, scale = T))
scale_cube_effect_pearson <- cor(scale_cube_effect, method = "pearson", use = "complete.obs")

# make one line correlation matrix with just effect size
corr_effect = matrix(scale_cube_effect_pearson[1, ], nrow = 1)
colnames(corr_effect) = colnames(scale_cube_effect_pearson)
rownames(corr_effect) = rownames(scale_cube_effect_pearson)[1]

# Make dataframe to plot pearson and LASSO together as two lines
corr_effect_df = as.data.frame(corr_effect) %>%
  reshape2::melt() %>%
  rename(Coefficients = value) %>%
  filter(Coefficients != 1) %>%
  mutate(y = "Pearson")

## Make LASSO DF (remove silt)
all_cube_variables = cube_effect %>%
  select(-c(cube_Percent_Silt))

# LASSO Analysis (Option 2: Single Best Model) ------------------------------------------
# Loop through LASSO to get average over 100 seeds  -----------------------

num_seeds = 100 # set number of seeds
seeds = sample(1:500, num_seeds) # choose 100 random seeds

## Set response variable (Cube_Effect_Size) and scale
yvar <- data.matrix(scale(all_cube_variables$cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, center = TRUE, scale = TRUE))
mean(yvar) # mean should be 0
sd(yvar) # sd should be 1

# this for non-cube root transformed data
#yvar <- data.matrix(scale(all_variables$Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, center = TRUE, scale = TRUE))

# list for storing LASSO iterations
norm_coeffs = list()
lasso_coefs_pull = list()
r2_scores = numeric(num_seeds)

## Set predictor variables and scale
exclude_col = "cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H"
#exclude_col = "Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H"

all_x_cube_variables = as.data.frame(scale(all_cube_variables[, !(names(all_cube_variables) %in% exclude_col)], center = T, scale = T))
#mean(all_x_cube_variables$Effect_Size_SpC_microsiemens_per_cm)
#sd(all_x_cube_variables$Effect_Size_SpC_microsiemens_per_cm)

# this for non-cube root transformed data
#all_x_cube_variables = as.data.frame(scale(all_variables[, !(names(all_variables) %in% exclude_col)], center = T, scale = T))

xvars <- data.matrix(all_x_cube_variables)

## Loop through LASSO seeds
for (i in 1:num_seeds) {
  
  seed = seeds[i]
  set.seed(seed)
  
  #perform cross-validation of LASSO regression
  lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                    standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                    #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                    # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
  ) # LASSO has various options for standardizing, we chose to scale earlier ourselves
  
  # this sets the penalty parameter (lambda) to the minimum value  
  best_lambda <- lasso$lambda.min
  #best_lambda <- lasso$lambda.1se # can also be set to 1se, makes a simpler model (less parameters)
  
  #best_lambda
  #plot(lasso) # check plot, should be v-shaped to ensure that the lambda parameter is really the best
  
  # run LASSO regression with best lambda
  best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
                             standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                             #  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                             #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
  )
  
  # normalize LASSO coefficients to top coefficient
  lasso_coefs = as.matrix(coef(best_lasso_model, s = best_lambda))
  
  lasso_coefs_pull[[as.character(seed)]] = lasso_coefs[-1, , drop = FALSE]
  
  norm_coeffs_scale = lasso_coefs/max(abs(lasso_coefs[-1]))
  
  norm_coeffs[[as.character(seed)]] = norm_coeffs_scale[-1, , drop = FALSE]
  
  y_pred = predict(best_lasso_model, newx = xvars, s = best_lambda)
  
  #calculate R2
  sst = sum((yvar - mean(yvar))^2)
  sse = sum((y_pred - yvar)^2)
  r2_scores[i] = 1 - (sse / sst)
  
}

# matrix of non-normalized LASSO coefficients
lasso_coef_mat = as.data.frame(do.call(cbind, lasso_coefs_pull)) 

colnames(lasso_coef_mat) = make.names(colnames(lasso_coef_mat), unique = T)

# Make DF of all LASSO results with mean and std. dev  
lasso_coef_means = lasso_coef_mat %>% 
  dplyr::mutate(variable = rownames(lasso_coef_mat)) %>% 
  rowwise() %>% 
  dplyr::mutate(mean = mean(c_across(dplyr::contains("s1"))), 
         sd = sd(c_across(dplyr::contains("s1")))) %>% 
  dplyr::relocate(mean, .before = s1.1) %>% 
  dplyr::relocate(sd, .before = s1.1) %>% 
  dplyr::relocate(variable, .before = mean)


# Bind all normalized LASSO results from 100 iterations
norm_coeffs_matrix = do.call(cbind, norm_coeffs)

mean_coeffs = as.data.frame(norm_coeffs_matrix, row.names = rownames(norm_coeffs_matrix))

colnames(mean_coeffs) = make.names(colnames(mean_coeffs), unique = T)

# Make DF of all LASSO results with mean and std. dev  
mean_coeffs_df = mean_coeffs %>%
  dplyr::select_if(~all(!is.nan(.))) %>% 
  dplyr::mutate(variable = rownames(mean_coeffs)) %>% 
  rowwise() %>% 
  dplyr::mutate(mean = mean(c_across(dplyr::contains("s1"))), 
         sd = sd(c_across(dplyr::contains("s1")))) %>% 
  dplyr::relocate(mean, .before = s1.1) %>% 
  dplyr::relocate(sd, .before = s1.1) %>% 
  dplyr::relocate(variable, .before = mean)

# calculate mean and SD R2
results_r2 = as.data.frame(r2_scores) 
mean(results_r2$r2_scores)
sd(results_r2$r2_scores)

# Make data frame for one-line pearson/lasso correlation matrix 

#non-normalized
ds_lasso_df = lasso_coef_means %>% 
  select(c(variable, mean)) %>% 
  rename(Coefficients = mean) %>% 
  mutate(y = "LASSO")

#normalized
norm_lasso_df = mean_coeffs_df = mean_coeffs_df %>% 
  select(c(variable, mean)) %>% 
  rename(Coefficients = mean) %>% 
  mutate(y = "LASSO")


# RANDOM FOREST ANALYSIS ------------------------------------------------

# STEP 1: Extract LASSO-Selected Variables
lasso_selected_vars = names(lasso_coeffs_clean)[lasso_coeffs_clean != 0]

cat("\n=== LASSO-RF INTEGRATION ===\n")
cat("LASSO selected", length(lasso_selected_vars), "variables:\n")
cat(paste(lasso_selected_vars, collapse = ", "), "\n\n")

# Prepare data for Random Forest
# Use the original cube-root transformed data (not scaled) for RF
rf_data = all_cube_variables %>%
  select(all_of(c("cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H", lasso_selected_vars)))

# Separate response and predictors
y_rf = rf_data$cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H
X_rf = rf_data %>% select(-cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)

cat("Random Forest data prepared:\n")
cat("  Response variable: cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H\n")
cat("  Number of predictors:", ncol(X_rf), "\n")
cat("  Number of samples:", nrow(X_rf), "\n\n")

# STEP 2: Single Random Forest with Permutation Importance
set.seed(123)  # for reproducibility

# Tune mtry for small sample size (more conservative)
mtry_value = max(1, floor(sqrt(ncol(X_rf))))

cat("Random Forest parameters:\n")
cat("  mtry:", mtry_value, "\n")
cat("  min.node.size: 3 (conservative for small n)\n")
cat("  num.trees: 1000\n\n")

# Fit single RF model
rf_single = ranger(
  x = X_rf,
  y = y_rf,
  num.trees = 1000,
  mtry = mtry_value,
  min.node.size = 3,  # Conservative for small sample size
  importance = "permutation",
  respect.unordered.factors = "order",
  seed = 123,
  keep.inbag = TRUE  # Keep for OOB calculations
)

# Print model performance
cat("=== SINGLE RF MODEL PERFORMANCE ===\n")
cat("R²:", round(rf_single$r.squared, 4), "\n")
cat("OOB Prediction Error (MSE):", round(rf_single$prediction.error, 4), "\n")
cat("RMSE:", round(sqrt(rf_single$prediction.error), 4), "\n")

# Compare to LASSO performance
cat("\nComparison to LASSO:\n")
cat("LASSO R²:", round(r2_final, 4), "\n")
cat("RF R²:   ", round(rf_single$r.squared, 4), "\n")
cat("RF improvement:", round(rf_single$r.squared - r2_final, 4), "\n\n")

# Extract single RF importance
imp_single = rf_single$variable.importance
imp_single_df = data.frame(
  variable = names(imp_single),
  importance = as.numeric(imp_single),
  row.names = NULL
) %>%
  arrange(desc(importance))

print("Single RF Permutation Importance:")
print(imp_single_df)

# STEP 3: Multiple RF Runs for Stability Assessment
cat("\n=== RUNNING MULTIPLE RF ITERATIONS FOR STABILITY ===\n")

set.seed(2025)
n_rf_runs = 100  # Match your LASSO approach
p_vars = length(lasso_selected_vars)

# Storage for results
imp_matrix = matrix(NA_real_, nrow = n_rf_runs, ncol = p_vars)
colnames(imp_matrix) = lasso_selected_vars
r2_scores_rf = numeric(n_rf_runs)
oob_errors = numeric(n_rf_runs)

# Progress tracking
pb = txtProgressBar(min = 0, max = n_rf_runs, style = 3)

for (i in seq_len(n_rf_runs)) {
  # Set different seed for each run
  set.seed(i + 1000)  # Offset to avoid overlap with LASSO seeds
  
  rf_i = ranger(
    x = X_rf,
    y = y_rf,
    num.trees = 1000,
    mtry = mtry_value,
    min.node.size = 3,
    importance = "permutation",
    respect.unordered.factors = "order",
    seed = i + 1000
  )
  
  # Store results
  imp_matrix[i, ] = rf_i$variable.importance
  r2_scores_rf[i] = rf_i$r.squared
  oob_errors[i] = rf_i$prediction.error
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Summarize importance across runs
imp_stability = data.frame(
  variable = lasso_selected_vars,
  mean_imp = apply(imp_matrix, 2, mean, na.rm = TRUE),
  sd_imp = apply(imp_matrix, 2, sd, na.rm = TRUE),
  median_imp = apply(imp_matrix, 2, median, na.rm = TRUE),
  min_imp = apply(imp_matrix, 2, min, na.rm = TRUE),
  max_imp = apply(imp_matrix, 2, max, na.rm = TRUE)
) %>%
  mutate(
    cv_imp = sd_imp / ifelse(abs(mean_imp) < 1e-8, NA, abs(mean_imp)),  # coefficient of variation
    stability = ifelse(is.na(cv_imp), "Stable", 
                       ifelse(cv_imp < 0.5, "Stable", 
                              ifelse(cv_imp < 1.0, "Moderate", "Unstable")))
  ) %>%
  arrange(desc(mean_imp))

# RF model performance summary
rf_performance = data.frame(
  metric = c("Mean R²", "SD R²", "Mean RMSE", "SD RMSE"),
  value = c(mean(r2_scores_rf), sd(r2_scores_rf), 
            mean(sqrt(oob_errors)), sd(sqrt(oob_errors)))
)

cat("\n=== RF STABILITY ANALYSIS RESULTS ===\n")
print("RF Performance across 100 runs:")
print(rf_performance)
cat("\nVariable Importance Stability:\n")
print(imp_stability)

# STEP 4: Comprehensive Visualizations

## 4.1 Single RF Importance Plot
single_rf_plot = ggplot(imp_single_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  coord_flip() +
  labs(title = "Random Forest Permutation Importance (Single Model)",
       subtitle = paste("R² =", round(rf_single$r.squared, 3), "| OOB RMSE =", round(sqrt(rf_single$prediction.error), 3)),
       x = "Variables", y = "Permutation Importance") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 9))

## 4.2 Stability of RF Importance  
stability_plot = ggplot(imp_stability, aes(x = reorder(variable, mean_imp), y = mean_imp)) +
  geom_point(aes(color = stability), size = 3) +
  geom_errorbar(aes(ymin = mean_imp - sd_imp, ymax = mean_imp + sd_imp, color = stability), 
                width = 0.2, alpha = 0.7) +
  coord_flip() +
  scale_color_manual(values = c("Stable" = "#2166AC", "Moderate" = "#FEB24C", "Unstable" = "#B2182B"),
                     name = "Stability") +
  labs(title = "RF Variable Importance Stability",
       subtitle = "Mean ± SD across 100 RF iterations",
       x = "Variables", y = "Mean Permutation Importance") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 9))

## 4.3 Comparison: LASSO vs RF Importance
# Merge LASSO and RF importance for comparison
lasso_rf_comparison = merge(
  data.frame(variable = names(lasso_coeffs_clean), 
             lasso_coeff = abs(as.numeric(lasso_coeffs_clean))),  # Use absolute values
  imp_stability %>% select(variable, rf_importance = mean_imp),
  by = "variable", all = TRUE
) %>%
  mutate(
    lasso_coeff = replace_na(lasso_coeff, 0),
    rf_importance = replace_na(rf_importance, 0)
  ) %>%
  filter(lasso_coeff > 0 | rf_importance > 0) %>%  # Keep if important in either method
  mutate(
    # Normalize to 0-1 scale for comparison
    lasso_norm = lasso_coeff / max(lasso_coeff, na.rm = TRUE),
    rf_norm = rf_importance / max(rf_importance, na.rm = TRUE)
  )

comparison_plot = ggplot(lasso_rf_comparison, aes(x = lasso_norm, y = rf_norm)) +
  geom_point(size = 3, alpha = 0.7, color = "darkgreen") +
  geom_text_repel(aes(label = variable), size = 2.5, max.overlaps = 15) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  labs(title = "LASSO vs Random Forest Variable Importance",
       subtitle = "Normalized importance scores (0-1 scale)",
       x = "LASSO Coefficient (normalized)", y = "RF Permutation Importance (normalized)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

## 4.4 Model Performance Comparison
perf_comparison = data.frame(
  Model = c("LASSO", "RF (Single)", "RF (Mean)"),
  R_squared = c(r2_final, rf_single$r.squared, mean(r2_scores_rf)),
  SD_R_squared = c(NA, NA, sd(r2_scores_rf))
)

perf_plot = ggplot(perf_comparison, aes(x = Model, y = R_squared)) +
  geom_col(fill = c("#E31A1C", "#1F78B4", "#33A02C"), alpha = 0.7) +
  geom_errorbar(aes(ymin = R_squared - SD_R_squared, ymax = R_squared + SD_R_squared), 
                width = 0.2, na.rm = TRUE) +
  geom_text(aes(label = round(R_squared, 3)), vjust = -0.5, fontface = "bold") +
  labs(title = "Model Performance Comparison", 
       y = "R² Score", x = "Model Type") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

## 4.5 Variable Importance Heatmap (Top Variables)
top_vars = imp_stability %>% 
  slice_head(n = min(10, nrow(imp_stability))) %>% 
  pull(variable)

# Create matrix for heatmap
heatmap_data = imp_matrix[, top_vars, drop = FALSE] %>%
  as.data.frame() %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(cols = -iteration, names_to = "variable", values_to = "importance")

heatmap_plot = ggplot(heatmap_data, aes(x = iteration, y = variable, fill = importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Importance") +
  labs(title = "Variable Importance Across RF Iterations",
       subtitle = paste("Top", length(top_vars), "variables by mean importance"),
       x = "RF Iteration", y = "Variables") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 8))

## 4.6 Coefficient of Variation Analysis
cv_plot = ggplot(imp_stability, aes(x = reorder(variable, cv_imp), y = cv_imp)) +
  geom_col(aes(fill = stability), alpha = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Stable" = "#2166AC", "Moderate" = "#FEB24C", "Unstable" = "#B2182B")) +
  geom_hline(yintercept = c(0.5, 1.0), linetype = "dashed", alpha = 0.5) +
  labs(title = "Variable Importance Stability (Coefficient of Variation)",
       subtitle = "Lower values indicate more stable importance",
       x = "Variables", y = "Coefficient of Variation (SD/Mean)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 9))

# STEP 5: Create Combined Visualization

# Arrange plots
top_row = (single_rf_plot | stability_plot)
middle_row = (comparison_plot | perf_plot)  
bottom_row = (cv_plot | heatmap_plot)

combined_rf_analysis = top_row / middle_row / bottom_row +
  plot_annotation(title = "Comprehensive Random Forest Analysis Following LASSO Selection",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))

print(combined_rf_analysis)
ggsave("Figures/RF_Comprehensive_Analysis.png", plot = combined_rf_analysis, 
       width = 16, height = 20, dpi = 300)

# Save individual plots
ggsave("Figures/RF_Single_Importance.png", plot = single_rf_plot, width = 10, height = 8, dpi = 300)
ggsave("Figures/RF_Stability_Analysis.png", plot = stability_plot, width = 10, height = 8, dpi = 300)
ggsave("Figures/LASSO_RF_Comparison.png", plot = comparison_plot, width = 10, height = 8, dpi = 300)

# STEP 6: Export Results

# Create comprehensive results dataframe
rf_results = list(
  lasso_selected_variables = lasso_selected_vars,
  single_rf_performance = list(
    r_squared = rf_single$r.squared,
    oob_error = rf_single$prediction.error,
    rmse = sqrt(rf_single$prediction.error)
  ),
  single_rf_importance = imp_single_df,
  stability_analysis = imp_stability,
  performance_summary = rf_performance,
  model_comparison = perf_comparison
)

# Export to CSV files
write.csv(imp_stability, "Data/RF_Variable_Importance_Stability.csv", row.names = FALSE)
write.csv(lasso_rf_comparison, "Data/LASSO_RF_Importance_Comparison.csv", row.names = FALSE)
write.csv(perf_comparison, "Data/Model_Performance_Comparison.csv", row.names = FALSE)

# Export importance matrix for further analysis
write.csv(imp_matrix, "Data/RF_Importance_Matrix_All_Iterations.csv", row.names = FALSE)

# STEP 7: LASSO Visualizations (for comparison)

# Make data frame for correlation matrix comparison
ds_lasso_df = data.frame(
  variable = names(lasso_coeffs_clean),
  mean = as.numeric(lasso_coeffs_clean)
) %>%
  rename(Coefficients = mean) %>%
  mutate(y = "LASSO")

norm_lasso_df = data.frame(
  variable = names(norm_lasso_coeffs),
  mean = as.numeric(norm_lasso_coeffs)
) %>%
  rename(Coefficients = mean) %>%
  mutate(y = "LASSO")

# Scatter Plots for manuscript
fs = ggplot(cube_effect, aes(y = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, x = cube_Percent_Fine_Sand)) +
  geom_point(size = 2, shape = 1) +
  theme_bw() +
  stat_cor(data = cube_effect, label.x = 1.1, label.y = 11, size = 3.5, digits = 2, cor.coef.name = "r", aes(label = paste(..r.label..)))+
  stat_poly_line(data = cube_effect, se = FALSE, linetype = 'dashed') +
  ylab("")+
  xlab(expression("Fine Sand (%)"^(1/3))) +
  theme(legend.position  = "none", aspect.ratio = 1, axis.title.x = element_text(size = 10))

atp = ggplot(cube_effect, aes(y = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, x = cube_Median_ATP_picomoles_per_g)) +
  geom_point(size = 2, shape  = 1) +
  theme_bw() +
  stat_cor(data = cube_effect, label.x = 1.1, label.y = 11, size = 3.5, digits = 2, cor.coef.name = "r", aes(label = paste(..r.label..)))+
  stat_poly_line(data = cube_effect, se = FALSE, linetype = 'dashed') +
  ylab("")+
  xlab(expression("Median ATP (pmol g"^-1*")"^(1/3)))+
  theme(legend.position  = "none", aspect.ratio = 1, axis.title.x = element_text(size = 10))

# Pull together scatter plots
col_scatter = ggarrange(fs, atp, ncol = 2, nrow = 1, common.legend =  T, legend = "right", align = "hv")

print(col_scatter)
ggsave("Figures/Example_Scatter_Plots.png", plot = col_scatter, width = 10, height = 5, dpi = 300)

# STEP 8: Summary Report

cat("\n" %>% strrep(3))
cat("=== COMPREHENSIVE ANALYSIS SUMMARY ===\n")
cat("LASSO Variable Selection:\n")
cat("  Variables selected:", length(lasso_selected_vars), "out of", ncol(all_x_cube_variables), "\n")
cat("  LASSO R²:", round(r2_final, 4), "\n")

cat("\nRandom Forest Analysis:\n")
cat("  Single RF R²:", round(rf_single$r.squared, 4), "\n")
cat("  Mean RF R² (100 runs):", round(mean(r2_scores_rf), 4), "±", round(sd(r2_scores_rf), 4), "\n")
cat("  RF improvement over LASSO:", round(mean(r2_scores_rf) - r2_final, 4), "\n")

cat("\nVariable Importance Stability:\n")
stable_vars = sum(imp_stability$stability == "Stable")
moderate_vars = sum(imp_stability$stability == "Moderate") 
unstable_vars = sum(imp_stability$stability == "Unstable")
cat("  Stable variables:", stable_vars, "\n")
cat("  Moderately stable:", moderate_vars, "\n")
cat("  Unstable variables:", unstable_vars, "\n")

cat("\nTop 5 Most Important Variables (by mean RF importance):\n")
top5 = imp_stability %>% slice_head(n = 5)
for(i in 1:5) {
  cat("  ", i, ".", top5$variable[i], "(", round(top5$mean_imp[i], 4), "±", round(top5$sd_imp[i], 4), ")\n")
}

cat("\nFiles Generated:\n")
cat("  - LASSO_CV_Curve.png\n")
cat("  - RF_Comprehensive_Analysis.png\n")
cat("  - RF_Single_Importance.png\n") 
cat("  - RF_Stability_Analysis.png\n")
cat("  - LASSO_RF_Comparison.png\n")
cat("  - RF_Variable_Importance_Stability.csv\n")
cat("  - LASSO_RF_Importance_Comparison.csv\n")
cat("  - Model_Performance_Comparison.csv\n")
cat("  - RF_Importance_Matrix_All_Iterations.csv\n")

cat("\n=== ANALYSIS COMPLETE ===\n")