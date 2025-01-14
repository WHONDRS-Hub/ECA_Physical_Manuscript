#### Sensitivity Analysis For ECA removals ####

# This script makes figures for ECA physical manuscript and performs PCA Analysis and LASSO regression

library(tidyverse);library(corrplot);library(ggpubr);library(ggpmisc);library(factoextra);library(stringr);library(glmnet);library(magick); library(ggnewscale); library(FSA); library(multcompView); library(rcompanion)

rm(list=ls());graphics.off()

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
getwd()

# Functions ---------------------------------------------------------------

# Transformation for normalization is cube root - have to cube root then add sign back to value to make it positive or negative
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

## Read in/Merge Data  ------------------------------------------------------------

# Individual Rate data for histograms ####
all_data = read.csv("./Data/EC_Sediment_SpC_pH_Temp_Respiration.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  # remove samples with too much water (EC_011, EC_012), sample with no mg/kg (EC_023), duplicated NEON sites (EC_052, EC_053, EC_057)
  select(-c(Field_Name, IGSN, Material))

# Take Cube root of all respiration rates
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
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 10), "FALSE", "TRUE")) %>% ## Check CV's, then remove 
  select(c(Sample_ID, Median_SpC_microsiemens_per_cm, Median_pH, Median_Temperature_degC, Median_Respiration_Rate_mg_DO_per_L_per_H, Median_Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  ungroup()

## Calculate "bulk" medians ####

# ATP ####

atp = read.csv("./Data/EC_Sediment_ATP.csv", skip = 2) %>% 
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

cn = read.csv("./Data/EC_Sediment_CN.csv", skip = 2) %>% 
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

npoc_tn = read.csv("./Data/EC_Sediment_NPOC_TN.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "SIR", "INC"))

# "Bulk" Medians - remove broken/missing samples
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

grav_inc = read.csv("./Data/EC_Sediment_Gravimetric_Moisture.csv", skip = 2) %>% 
  slice(-1:-11) %>% 
  filter(Field_Name != "#End_Data") %>% 
  select(-c(Field_Name, IGSN, Material)) %>% 
  mutate(across(c(Initial_Water_Mass_g:Incubation_Water_Mass_g), as.numeric))

#Some dry reps have high CV: EC_057 (low moisture), EC_081 (one lower sample), EC_063 (low moisture), EC_088 (one lower sample), EC_076 (low moisture), EC_071 (one lower sample), EC_056 (low moisture), EC_069 (one slightly lower) 

median_grav = grav_inc %>%
  mutate(X62948_Initial_Gravimetric_Moisture_g_per_g = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, X62948_Initial_Gravimetric_Moisture_g_per_g)) %>%
  mutate(X62948_Final_Gravimetric_Moisture_g_per_g = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, X62948_Final_Gravimetric_Moisture_g_per_g)) %>%
  #missing replicates (EC_072-W5/D5), less sediment in sample (EC_012-D5)
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>%
  mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
  group_by(Sample_ID, Rep) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm =TRUE),
                        n = ~sum(!is.na(.x))),
                   .names = "{.fn}_{.col}")) %>%
  ungroup() %>%
  group_by(Sample_ID, Rep) %>%
  mutate(Remove = ifelse(all(c_across(starts_with("n_")) == 5), "FALSE", "TRUE")) %>% ## Check CV's, then remove
  select(c(Sample_ID, Rep, Median_X62948_Initial_Gravimetric_Moisture_g_per_g, Median_X62948_Final_Gravimetric_Moisture_g_per_g))

# Fe #### 

fe = read.csv("./Data/EC_Sediment_Fe.csv", skip = 2) %>% 
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


# Join All Data -----------------------------------------------------------

# remove missing/broken samples
# set Fe LOD again
all_join = all_data %>% 
  left_join(atp, by = "Sample_Name") %>% 
  left_join(fe, by = "Sample_Name") %>% 
  left_join(cn, by = "Sample_Name") %>% 
  left_join(npoc_tn, by = "Sample_Name") %>% 
  left_join(grav_inc, by = "Sample_Name") %>% 
  mutate(Fe_mg_per_L = ifelse(grepl("INC_Method_001|INC_Method_002", Dev_Fe), NA, Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_kg = ifelse(grepl("INC_Method_001|INC_Method_002", Dev_Fe), NA, Fe_mg_per_kg)) %>% 
  mutate(Fe_mg_per_L = if_else(grepl("SFE_Below", Fe_mg_per_L), "0.001", Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_L = if_else(grepl("SFE_Above", Fe_mg_per_L), str_extract(Fe_mg_per_L, "(?<=\\|[^|]{1,100}\\|)\\d+\\.\\d+"), Fe_mg_per_L)) %>% 
  mutate(Fe_mg_per_L = as.numeric(Fe_mg_per_L)) %>%
  mutate(Fe_mg_per_kg = if_else(grepl("SFE_Above", Fe_mg_per_kg), str_extract(Fe_mg_per_kg, "(?<=\\|[^|]{1,100}\\|)\\d+\\.\\d+"), Fe_mg_per_kg)) %>% 
  mutate(Fe_mg_per_kg = if_else(grepl("SFE_Below", Fe_mg_per_kg), as.numeric(Fe_mg_per_L * (Incubation_Water_Mass_g/Dry_Sediment_Mass_g)), as.numeric(Fe_mg_per_kg))) %>%
  mutate(Fe_mg_per_kg = as.numeric(Fe_mg_per_kg)) %>% 
  mutate(ATP_nanomoles_per_L = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, ATP_nanomoles_per_L)) %>% 
  # ATP_002 (don't have sample), INC_Method_001
  mutate(ATP_picomoles_per_g = ifelse(grepl("INC_Method_001", Methods_Deviation), NA, ATP_picomoles_per_g)) %>%
  mutate(ATP_nanomoles_per_L = ifelse(ATP_nanomoles_per_L == -9999, NA, ATP_nanomoles_per_L)) %>% 
  mutate(ATP_picomoles_per_g = ifelse(ATP_picomoles_per_g == -9999, NA, ATP_picomoles_per_g)) %>% 
  mutate(across(c(ATP_nanomoles_per_L:ATP_picomoles_per_g), as.numeric)) %>%
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(grepl("INC_Method_001|INC_Method_002|INC_QA_004", Methods_Deviation), NA, Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(Respiration_Rate_mg_DO_per_kg_per_H == "-9999", NA, Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  mutate(X62948_Initial_Gravimetric_Moisture_g_per_g = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, X62948_Initial_Gravimetric_Moisture_g_per_g)) %>%
  mutate(X62948_Final_Gravimetric_Moisture_g_per_g = ifelse(grepl("INC_Method_001|INC_Method_002", Methods_Deviation), NA, X62948_Final_Gravimetric_Moisture_g_per_g)) %>%
  select(c(Sample_Name, Respiration_Rate_mg_DO_per_kg_per_H, ATP_picomoles_per_g, Fe_mg_per_kg, X01395_C_percent_per_mg, X01397_N_percent_per_mg, Extractable_NPOC_mg_per_kg, Extractable_TN_mg_per_kg, SpC_microsiemens_per_cm, X62948_Final_Gravimetric_Moisture_g_per_g)) 

## Median Data ####

## Bulk Medians

all_medians = median_respiration %>% 
  left_join(median_atp, by = "Sample_ID") %>% 
  left_join(median_iron, by = "Sample_ID") %>% 
  left_join(median_cn, by = "Sample_ID") %>% 
  left_join(median_npoc_tn, by = "Sample_ID") %>%   mutate(Sample_Name = str_replace(Sample_ID, "INC", "all")) %>% 
  select(-c(Sample_ID)) %>% 
  relocate(Sample_Name, .before = Median_SpC_microsiemens_per_cm)

## Read in Wet/Dry Median Data from Summary file to get Dry moisture values

median = read.csv("./Data/EC_Sediment_Sample_Data_Summary.csv", skip = 2) %>% 
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

effect = read.csv("./Data/EC_Sediment_Effect_Size.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%
 select(-c(IGSN, Field_Name, Material, Methods_Deviation)) 

## Read in grain size/ssa variables ####

grain = read.csv("./Data/v3_CM_SSS_Sediment_Sample_Data_Summary.csv", skip = 2) %>% 
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
  mutate(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = abs(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  mutate(across(contains("percent"), as.numeric))
  
# Transform Data ----------------------------------------------------------

# Fe outlier not included in analysis - remove from DF
cube_effect = effect_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>% 
  column_to_rownames("Sample_Name") %>%
  filter(cube_Effect_Size_Fe_mg_per_kg > -1) %>%  # remove Fe outlier for analysis 
  select(-contains("per_L")) %>% 
  relocate(cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, .before = cube_Effect_Size_SpC_microsiemens_per_cm)

# Data frame that includes Fe outlier - this used for Scatter Plot
fe_cube_effect = effect_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>% 
  column_to_rownames("Sample_Name") %>%
  select(-contains("per_L")) %>% 
  filter(cube_Effect_Size_Fe_mg_per_kg < -1) %>%
  relocate(cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, .before = cube_Effect_Size_SpC_microsiemens_per_cm)


# Histograms --------------------------------------------------------------

## Histogram of all Rates
# These used with Nathan Johnson annotations in manuscript

all_cube_hist = ggplot(cube_respiration, aes(x = cube_Respiration_mg_kg)) +
    geom_histogram(position = "identity", alpha = 0.8, aes(fill = Treat))+
    scale_fill_manual(values = c("#D55E00","#0072B2"))  +
    theme(strip.text = element_text(
      size = 4))+
    ylim(0, 87.5)+
    theme_bw() + 
    theme(legend.position = c(0.85, 0.8), 
          legend.key.size = unit(0.15, "in"), 
          legend.title = element_text(size = 8),
          axis.title.x = element_text(size = 10)) +
    guides(fill = guide_legend(title="Treatment")) + 
    xlab(expression(atop("\n Respiration Rate" ^(1/3)*"", "(mg O"[2]*" kg"^-1*" H"^-1*")"))) +
    ylab("Count")
  
  all_cube_hist

  # Histogram of effect size

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
    xlab(expression(atop("\n Effect Size Respiration Rate" ^(1/3)*"","(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" H"^-1*")")))
  
cube_effect_hist

combined_hist = ggarrange(all_cube_hist, cube_effect_hist, ncol = 2, labels = c("A", "B"), hjust = -5, vjust = 2.5)

ggsave("./Physical_Manuscript_Figures/Combined_Cube_Histograms.pdf", plot = combined_hist, width = 10, height = 5)

## Pearson Correlation Matrix ####

## All Medians
# scale data before it goes into correlation matrix

scale_cube_effect = as.data.frame(scale(cube_effect, center = T, scale = T))

scale_cube_effect_pearson <- cor(scale_cube_effect, method = "pearson", use = "complete.obs")

# png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_Pearson_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)
#   
# corrplot(scale_cube_effect_pearson,type = "upper", method = "number", tl.col = "black", tl.cex = 0.8, cl.cex = 0.5, number.cex = 0.5,  title = "Effect Samples Pearson Correlation")
#   
# 
# dev.off()

# make one line correlation matrix with just effect size

corr_effect = matrix(scale_cube_effect_pearson[1, ], nrow = 1)

colnames(corr_effect) = colnames(scale_cube_effect_pearson)

rownames(corr_effect) = rownames(scale_cube_effect_pearson)[1]

# Try to plot pearson and LASSO together as two lines

corr_effect_df = as.data.frame(corr_effect) %>% 
  reshape2::melt() %>% 
  rename(Coefficients = value) %>% 
  filter(Coefficients != 1) %>% 
  mutate(y = "Pearson")

## Downselected LASSO - Loop through coefficients to choose for LASSO ####

## All Medians

all_pearson_df = as.data.frame(scale_cube_effect_pearson)

row_names_pearson <- rownames(all_pearson_df)

all_pearson_df$Variable <- row_names_pearson

# Melt the dataframe for plotting
all_pearson_melted <- reshape2::melt(all_pearson_df, id.vars = "Variable") %>% 
  filter(value != 1) %>% 
  mutate(value = abs(value)) %>% 
  filter(!grepl("Respiration", Variable))

all_effect_melted <- all_pearson_melted %>% 
  filter(grepl("Effect_Size_Respiration", variable)) %>%
  filter(!grepl("Silt", Variable)) # remove silt from variables, 0 values so not using

all_choose_melted <- all_pearson_melted %>% 
  filter(!grepl("Respiration", variable)) %>%
  filter(!grepl("Silt", variable)) %>%
  filter(!grepl("Silt", Variable)) %>% #try removing silt (has 0 values)
  #distinct(value, .keep_all = TRUE) %>% 
  left_join(all_effect_melted, by = "Variable") %>% 
  rename(Variable_1 = Variable) %>% 
  rename(Variable_2 = variable.x) %>% 
  rename(Correlation = value.x) %>% 
  rename(Variable_1_Effect_Correlation = value.y) %>% 
  select(-c(variable.y)) %>% 
  left_join(all_effect_melted, by = c("Variable_2" = "Variable")) %>% 
  rename(Variable_2_Effect_Correlation = value) %>% 
  select(-c(variable))

all_loop_melt = all_choose_melted %>% 
  arrange(desc(Correlation))

# Pearson correlation coefficient to remove above
correlation = 0.7

## Start loop to remove highly correlated (> 0.7)
all_effect_filter = function(all_loop_melt) {
  
  rows_to_keep = rep(TRUE, nrow(all_loop_melt))
  
  for (i in seq_len(nrow(all_loop_melt))) {
    
    if (!rows_to_keep[i]) next
    
    row = all_loop_melt[i, ]
    
    if (row$Correlation < correlation) next
    
    if(row$Variable_1_Effect_Correlation >= row$Variable_2_Effect_Correlation) {
      
      var_to_keep = row$Variable_1
      var_to_remove = row$Variable_2
      
    } else {
      
      var_to_keep = row$Variable_2
      var_to_remove = row$Variable_1
      
    }
    
    all_loop_melt$Variable_to_Keep[i] = var_to_keep
    all_loop_melt$Variable_to_Remove[i] = var_to_remove
    
    for (j in seq(i + 1, nrow(all_loop_melt))) {
      
      if(all_loop_melt$Variable_1[j] == var_to_remove || all_loop_melt$Variable_2[j] == var_to_remove) {
        
        rows_to_keep[j] = FALSE
        
      }
      
    }
    
    
  }
  
  return(all_loop_melt[rows_to_keep, ])
  
}

all_filtered_data = all_effect_filter(all_loop_melt) 

# pull out variables to remove
all_removed_variables = all_filtered_data %>% 
  distinct(Variable_to_Remove)

# pull out all variables 
all_variables = all_effect_melted %>% 
  select(c(Variable))

# remove variables from all variables to get variables to keep for LASSO 
all_kept_variables = all_effect_melted[!(all_effect_melted$Variable %in% all_removed_variables$Variable_to_Remove), ] #keeps SpC, Temp, pH, ATP, NPOC, TN (ext), TOC, TN (solid), med sand, silt

# if silt is removed, keeps SpC, Temp, pH, ATP, NPOC, TOC, TN (solid), med sand, fine sand, lost grav. moisture

## LASSO VARIABLES ####

# Keep variables selected from down-selected correlation matrix and add Cube_Effect_Size
col_to_keep = unique(all_kept_variables$Variable)
col_to_keep = c(col_to_keep, "cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H")

all_cube_variables = cube_effect[, col_to_keep, drop = FALSE]

## Loop through LASSO to get average over a lot of seeds ####

num_seeds = 100
seeds = sample(1:500, num_seeds)

## Set response variable (Cube_Effect_Size) and scale
yvar <- data.matrix(scale(all_cube_variables$cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, center = TRUE, scale = TRUE))
mean(yvar)
sd(yvar)

# list for storing LASSO iterations
norm_coeffs = list()
lasso_coefs_pull = list()
r2_scores = numeric(num_seeds)

## Set predictor variables and scale
exclude_col = "cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H"

all_x_cube_variables = as.data.frame(scale(all_cube_variables[, !(names(all_cube_variables) %in% exclude_col)], center = T, scale = T))
#mean(x_cube_variables$Cube_SpC_Diff)
#sd(x_cube_variables$Cube_SpC_Diff)

xvars <- data.matrix(all_x_cube_variables)

for (i in 1:num_seeds) {
  
  seed = seeds[i]
  set.seed(seed)
  
  lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                    standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                    #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                    # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
  )
  
  best_lambda <- lasso$lambda.min
  #best_lambda
  #plot(lasso)
  
  best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
                             standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                             #  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                             #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
  )
  
  lasso_coefs = as.matrix(coef(best_lasso_model, s = best_lambda))
  
  lasso_coefs_pull[[as.character(seed)]] = lasso_coefs[-1, , drop = FALSE]
  
  norm_coeffs_scale = lasso_coefs/max(abs(lasso_coefs[-1]))
  
  norm_coeffs[[as.character(seed)]] = norm_coeffs_scale[-1, , drop = FALSE]
  
  y_pred = predict(best_lasso_model, newx = xvars, s = best_lambda)
  
  sst = sum((yvar - mean(yvar))^2)
  sse = sum((y_pred - yvar)^2)
  r2_scores[i] = 1 - (sse / sst)
  
}

lasso_coef_mat = as.data.frame(do.call(cbind, lasso_coefs_pull)) 

colnames(lasso_coef_mat) = make.names(colnames(lasso_coef_mat), unique = T)

# Make DF of all LASSO results with mean and std. dev  
lasso_coef_means = lasso_coef_mat %>% 
  mutate(variable = rownames(lasso_coef_mat)) %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(contains("s1"))), 
         sd = sd(c_across(contains("s1")))) %>% 
  relocate(mean, .before = s1) %>% 
  relocate(sd, .before = s1) %>% 
  relocate(variable, .before = mean)


# Bind all LASSO results from 100 iterations
norm_coeffs_matrix = do.call(cbind, norm_coeffs)

mean_coeffs = as.data.frame(norm_coeffs_matrix, row.names = rownames(norm_coeffs_matrix))

colnames(mean_coeffs) = make.names(colnames(mean_coeffs), unique = T)

# Make DF of all LASSO results with mean and std. dev  
mean_coeffs_df = mean_coeffs %>% 
  mutate(variable = rownames(mean_coeffs)) %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(contains("s1"))), 
         sd = sd(c_across(contains("s1")))) %>% 
  relocate(mean, .before = s1) %>% 
  relocate(sd, .before = s1) %>% 
  relocate(variable, .before = mean)

results_r2 = as.data.frame(r2_scores) 
mean(results_r2$r2_scores)
sd(results_r2$r2_scores)

ds_lasso_df = lasso_coef_means %>% 
  select(c(variable, mean)) %>% 
  rename(Coefficients = mean) %>% 
  mutate(y = "LASSO")

### Make Figures ---------------------------

# Scatter Plots #### 
fs = ggplot(cube_effect, aes(y = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, x = cube_Percent_Fine_Sand)) +
  geom_point(size = 2, shape = 1) +
  theme_bw() +
  #stat_cor(data = fe_cube_out, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_cor(data = cube_effect, label.x = 1.1, label.y = 11, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_effect, label.x = 1.1, label.y = 10.25, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_effect, se = FALSE, linetype = 'dashed') + 
  #ylab(expression("Effect Size Respiration Rate (mg kg"^-1*")"^(1/3))) +
  ylab("")+
  xlab(expression("Fine Sand (%)"^(1/3))) + 
  theme(legend.position  = "none", aspect.ratio = 1, axis.title.x = element_text(size = 10))

atp = ggplot(cube_effect, aes(y = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, x = cube_Median_ATP_picomoles_per_g)) +
  geom_point(size = 2, shape  = 1) +
  theme_bw() +
  #stat_cor(data = fe_cube_out, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_cor(data = cube_effect, label.x = 1.1, label.y = 11, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_effect, label.x = 1.1, label.y = 10.25, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_effect, se = FALSE, linetype = 'dashed') + 
  #ylab(expression("Effect Size Respiration Rate (mg kg"^-1*")"^(1/3))) +
  ylab("")+
  xlab(expression("Median ATP (pmol g"^-1*")"^(1/3)))+ 
  theme(legend.position  = "none", aspect.ratio = 1, axis.title.x = element_text(size = 10))

toc = ggplot(cube_effect, aes(y = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, x = cube_Median_X01395_C_percent_per_mg)) +
  geom_point(shape = 1, size = 2) +
  theme_bw() +
  #stat_cor(data = fe_cube_out, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_cor(data = cube_effect, label.x = 0.55, label.y = 11, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_effect, label.x = 0.55, label.y = 10.25, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_effect, se = FALSE, linetype = 'dashed') + 
  #ylab(expression("Effect Size Respiration Rate (mg kg"^-1*")"^(1/3))) +
  ylab("")+
  xlab(expression("Median TOC (%)"^(1/3)))+ 
  theme(legend.position  = "none", aspect.ratio = 1, axis.title.x = element_text(size = 10))

tn = ggplot(cube_effect, aes(y = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, x = cube_Median_X01397_N_percent_per_mg)) +
  geom_point(shape =1, size = 2) +
  theme_bw() +
  #stat_cor(data = fe_cube_out, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_cor(data = cube_effect, label.x = 0.225, label.y = 11, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_effect, label.x = 0.225, label.y = 10.25, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_effect, se = FALSE, linetype = 'dashed')+ 
 # ylab("Effect Size Respiration Rate (mg/kg)") +
  ylab("")+
  xlab(expression("Median TN (%)"^(1/3)))+ 
  theme(legend.position  = "none", aspect.ratio = 1, axis.title.x = element_text(size = 10))

spc = ggplot(cube_effect, aes(y = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, x = cube_Median_SpC_microsiemens_per_cm)) +
  geom_point(shape = 1, size = 2) +
  theme_bw() +
  #stat_cor(data = fe_cube_out, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_cor(data = cube_effect, label.x = 2.5, label.y = 11, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_effect, label.x = 2.5, label.y = 10.25, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  #stat_poly_line(data = cube_effect, se = FALSE, linetype = 'dashed')+ 
  #ylab("Effect Size Respiration Rate (mg/kg)") +
  ylab("")+
  xlab(expression("Median SpC (\u03BCS cm"^-1*")"^(1/3)))+ 
  theme(legend.position  = "none", aspect.ratio = 1, axis.title.x = element_text(size = 10))

fe = ggplot(cube_effect, aes(y = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, x = cube_Effect_Size_Fe_mg_per_kg)) +
  geom_point(shape = 1, size = 2) +
  geom_point(fe_cube_effect, mapping = aes(y = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, x = cube_Effect_Size_Fe_mg_per_kg), color = "red", size = 2) +
  theme_bw() +
  #stat_cor(data = fe_cube_out, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
  stat_cor(data = cube_effect, label.x = -2.5, label.y = 11, size = 3, digits = 2, aes(label = paste(..rr.label..)))+
  stat_cor(data = cube_effect, label.x = -2.5, label.y = 10.25, size = 3, digits = 2, aes(label = paste(..p.label..)))+
  stat_poly_line(data = cube_effect, se = FALSE, linetype = 'dashed')+ 
  #ylab(expression("Effect Size Respiration Rate (mg kg"^-1*")"^(1/3))) +
  ylab("")+
  xlab(expression("Effect Size Fe (II) (mg kg"^-1*")"^(1/3)))+ 
  theme(legend.position  = "none", aspect.ratio = 1, axis.title.x = element_text(size = 10))

## Combine LASSO and Pearson Coefficients into Heat Map ####

# Update names for figure
lasso_pear_df = bind_rows(ds_lasso_df, corr_effect_df) %>% 
  rename(type = y) %>% 
  filter(!grepl("Silt", variable)) %>% 
  mutate(variable = ifelse(variable == "cube_Effect_Size_SpC_microsiemens_per_cm", "Effect Size Specific Conductivity (\u03BCS cm\u207b\u00b9)", ifelse(variable == "cube_Effect_Size_pH", "Effect Size pH", ifelse(variable == "cube_Effect_Size_Temperature_degC", " Effect Size Temperature (\u00B0C)",  ifelse(variable == "cube_Effect_Size_Fe_mg_per_kg", "Effect Size Fe (II) (mg kg\u207b\u00b9)", ifelse(variable == "cube_Effect_Size_ATP_picomoles_per_g", "Effect Size ATP (pmol g\u207b\u00b9)", ifelse(variable == "cube_Effect_Size_Extractable_NPOC_mg_per_kg", "Effect Size Extractable NPOC (mg kg\u207b\u00b9)", ifelse(variable == "cube_Effect_Size_Extractable_TN_mg_per_kg", "Effect Size Extractable TN (mg kg\u207b\u00b9)", ifelse(variable == "cube_Effect_Size_C_percent_per_mg", "Effect Size TOC (%)", ifelse(variable == "cube_Effect_Size_N_percent_per_mg", "Effect Size TN (%)", ifelse(variable == "cube_Percent_Tot_Sand", "Total Sand (%)", ifelse(variable == "cube_Percent_Med_Sand", "Medium Sand (%)", ifelse(variable == "cube_Percent_Fine_Sand", "Fine Sand (%)", ifelse(variable == "cube_Median_SpC_microsiemens_per_cm", "Median Specific Conductivity (\u03BCS cm\u207b\u00b9)", ifelse(variable == "cube_Median_pH", "Median pH", ifelse(variable == "cube_Median_Temperature_degC", "Median Temperature (\u00B0C)", ifelse(variable == "cube_Median_ATP_picomoles_per_g", "Median ATP (pmol g\u207b\u00b9)",  ifelse(variable == "cube_Median_Fe_mg_per_kg", "Median Fe (II) (mg kg\u207b\u00b9)", ifelse(variable == "cube_Median_X01395_C_percent_per_mg", "Median TOC (%)",   ifelse(variable ==  "cube_Median_X01397_N_percent_per_mg", "Median TN (%)", ifelse(variable == "cube_Median_Extractable_NPOC_mg_per_kg", "Median Extractable NPOC (mg kg\u207b\u00b9)", ifelse(variable == "cube_Median_Extractable_TN_mg_per_kg", "Median Extractable TN (mg kg\u207b\u00b9)", ifelse(variable == "cube_median_Dry_Initial_Gravimetric", "Median Initial Dry Gravimetric Moisture (g g\u207b\u00b9)", ifelse(variable == "cube_Percent_Coarse_Sand", "Coarse Sand (%)", ifelse(variable == "cube_Percent_Clay", "Clay (%)", ifelse(variable == "cube_Mean_Specific_Surface_Area_m2_per_g", "Specific Surface Area (m\u00B2 g\u207b\u00b9)", ifelse(variable == "cube_median_Dry_Final_Gravimetric", "Median Final Dry Gravimetric Moisture (g g\u207b\u00b9)",
 variable)))))))))))))))))))))))))))

# Pull out max and min LASSO coeffiecients
lasso_df = lasso_pear_df %>% 
  filter(type == "LASSO")%>% 
  arrange(Coefficients)

max_l = max(lasso_df$Coefficients)
min_l = min(lasso_df$Coefficients)

##Pull out max and min Pearson coefficients, and reorder based on variable
pearson_df = lasso_pear_df %>% 
  filter(type == "Pearson") %>% 
  arrange(Coefficients) %>% 
  mutate(variable = factor(variable, levels = variable[order(Coefficients)]))

max_p = max(pearson_df$Coefficients)
min_p = min(pearson_df$Coefficients)

# Set color scales for heat maps
color_palette_p = colorRampPalette(c("#B2182B", "#F7F7F7", "#2166AC"))(200)

color_palette_l = colorRampPalette(c("#F7F7F7", "#2166AC"))(200)

# Make heat map plot
combined_matrix = ggplot() +
  geom_tile(pearson_df, fill = "white", color = "black", mapping = aes(variable, type)) +
  geom_text(pearson_df, mapping = aes(x = variable, y = type, label = round(Coefficients, 2), color = Coefficients), size = 3, fontface = "bold") + 
  scale_color_gradientn(colors = color_palette_p, 
                        limit = c(min_p, max_p),
                        guide = "none") +
  new_scale_color() +
  geom_tile(lasso_df, fill = "white", color = "black", mapping = aes(variable, type)) +
  geom_text(lasso_df, mapping = aes(x = variable, y = type, label = round(Coefficients, 2), color = Coefficients), size = 3, fontface = "bold") + 
  scale_color_gradientn(colors = color_palette_l, 
                        limit = c(min_l, max_l),
                        guide = "none") +
  theme_bw() + 
  theme(aspect.ratio = 0.1, 
        axis.text.x = element_text(angle = 90, hjust = 0, size = 10), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())+#, 
  #axis.text.y = element_blank()) +
  scale_x_discrete(position = "top")
 
combined_matrix
 
ggsave("./Physical_Manuscript_Figures/Combined_Pearson_Lasso_Matrix.png", plot = combined_matrix, width = 9, height = 4)

# Merge Matrices + Scatter Plots ####

 # Pull together scatter plots
col_scatter = ggarrange(fs, atp, tn, fe, toc, spc, ncol = 3, nrow = 2, common.legend =  T, legend = "right",  labels = c("B", "C", "D", "E", "F", "G", "H"), label.x = 0.875, label.y = 0.28, align = "hv", heights = c(1,1), font.label = list(size = 12))

 # Annotate Figure by adding common "Effect Size" y-axis
col_scatter_ann = annotate_figure(col_scatter, left = text_grob(expression("Effect Size Respiration Rate mg kg"^-1*")"^(1/3)), rot = 90, size = 12))

col_scatter_ann

ggsave("./Physical_Manuscript_Figures/Scatter_Plots.png", plot = col_scatter_ann, width = 9, height = 6)


## Make sure to push figures - it was easier to combine/scale with saved figures. You will pull them back in and combine here

combine_hm_image = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/Combined_Pearson_Lasso_Matrix.png")

combine_label_image = image_annotate(combine_hm_image, "A", size = 100, location = "+25+25", color = "black")

scatter_plot = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/Scatter_Plots.png")

scatter_scale = image_scale(scatter_plot, "100%")

# Get the original dimensions of scatter_scale
scatter_width <- image_info(scatter_scale)$width
scatter_height <- image_info(scatter_scale)$height

# Add padding to the left by increasing the canvas width
scatter_scale_padded <- image_extent(
  scatter_scale,
  geometry = paste0(scatter_width, "x", scatter_height), # Add 200 pixels to the width
  gravity = "center" # Content stays on the left, padding is added to the left
)

whole_image = image_append(c(combine_label_image, scatter_scale_padded), stack = TRUE)

image_write(whole_image, path = "./Physical_Manuscript_Figures/Scatter_Heat_Map.png")

## Conceptual Model ####

d50 = read.csv("./Data/D50_Calculations.csv")
  
effect_d50 = left_join(effect_data, d50, by = "Sample_Name")
 
dummy_row = data.frame(d50 = 0.053)

d50_plot = effect_d50 %>%  
  bind_rows(dummy_row) %>% 
  ggplot(aes(x = round(d50, 2), y = Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) + 
   geom_bar(width = 0.005, stat = "identity", aes(fill = Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) +
  scale_fill_gradient2(name = "Effect Size", limits = c(-1400, 1400), low = "firebrick2", mid = "goldenrod2",
                             high = "dodgerblue2", midpoint = (max(1400)+min(-1400))/2) +
   geom_vline(xintercept = 0.053, linetype = 2) + 
   geom_vline(xintercept = 0.25, linetype = 2)+
  ylab("Effect Size Respiration (mg/kg)") +
  xlab("D50") +
   theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

fine_sand = effect_d50 %>% 
  rename(Sample = Sample_Name) %>% 
  select(c(Sample, Percent_Fine_Sand, Mean_Specific_Surface_Area_m2_per_g))
 
all_d50 = all_join %>% 
  separate(Sample_Name, c("Sample", "Rep"), sep = "-", remove = FALSE) %>% 
  mutate(Sample = str_replace(Sample, "INC", "all")) %>% 
  left_join(d50, by = c("Sample" = "Sample_Name")) %>% 
  mutate(Treat = ifelse(grepl("W", Rep), "Wet", "Dry")) %>% 
  mutate(category = cut(d50, breaks = c(0, 0.053, 0.25, 2), 
                        labels = c("Clay/Silt", "Fine Sand", "Med/Coarse Sand"), 
                        include.lowest = T, right = F)) %>% 
  mutate(TN = as.numeric(X01397_N_percent_per_mg)) %>% 
  select(-c(X01397_N_percent_per_mg)) %>% 
  filter(TN != -9999) %>% 
  left_join(fine_sand) %>% 
  mutate(SpC_microsiemens_per_cm = as.numeric(SpC_microsiemens_per_cm)) %>% 
  mutate(TOC = as.numeric(X01395_C_percent_per_mg)) %>% 
  select(-c(X01395_C_percent_per_mg)) %>% 
  rename(Final_Grav = X62948_Final_Gravimetric_Moisture_g_per_g) %>% 
  mutate(Grav_by_SSA = Final_Grav/Mean_Specific_Surface_Area_m2_per_g)


effect_analysis = effect_d50 %>% 
  select(c(Sample_Name, Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, Effect_Size_Fe_mg_per_kg, d50, Percent_Fine_Sand,  Median_ATP_picomoles_per_g, Median_X01397_N_percent_per_mg, Median_X01395_C_percent_per_mg, Median_SpC_microsiemens_per_cm, median_Dry_Final_Gravimetric)) %>% 
  mutate(category = cut(d50, breaks = c(0, 0.053, 0.25, 2), 
                        labels = c("Clay/Silt", "Fine Sand", "Med/Coarse Sand"), 
                        include.lowest = T, right = F)) %>% 
  mutate(Effect_Size_Fe_mg_per_kg = ifelse(Effect_Size_Fe_mg_per_kg < -10, NA, Effect_Size_Fe_mg_per_kg))


## ANOVA + Pairwise Tests for CLD Boxplots ####
cube_effect_analysis = effect_analysis %>% 
  mutate(cube_Effect = cube_root(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  mutate(cube_FS = cube_root(Percent_Fine_Sand)) %>% 
  mutate(cube_Fe_Effect = cube_root(Effect_Size_Fe_mg_per_kg)) %>% 
  mutate(cube_Median_ATP = cube_root(Median_ATP_picomoles_per_g)) %>% 
  mutate(cube_TN = cube_root(Median_X01397_N_percent_per_mg)) %>% 
  mutate(cube_TOC = cube_root(Median_X01395_C_percent_per_mg)) %>% 
  mutate(cube_Median_SpC = cube_root(Median_SpC_microsiemens_per_cm)) %>% 
  mutate(cube_Final_Dry_Grav = cube_root(median_Dry_Final_Gravimetric))

## Effect Size 
cube_aov = aov(cube_Effect ~ category, data = cube_effect_analysis)
summary(cube_aov)

# Pairwise - same results 
#Clay/Silt vs. Fine sand - not sig. (0.36)
#Clay/Silt vs. Med/Coarse sand - not sig. (0.17)
#Fine sand vs. Med/Coarse sand - p = 0.000007
tukey_eff = TukeyHSD(cube_aov, "category")
cld_eff = multcompLetters4(cube_aov, tukey_eff)
cld_eff = as.data.frame.list(cld_eff$category)

effect_analysis$cld_eff = cld_eff$Letters[match(effect_analysis$category, rownames(cld_eff))]


## Fine Sand
fs_aov= aov(cube_FS ~ category, data = cube_effect_analysis)
summary(fs_aov)

# Pairwise - same results 
#Clay/Silt vs. Fine sand - (0.005)
#Clay/Silt vs. Med/Coarse sand - not sig. (0.22)
#Fine sand vs. Med/Coarse sand - p = 0.00000
tukey_fs = TukeyHSD(fs_aov, "category")
cld_fs = multcompLetters4(fs_aov, tukey_fs)
cld_fs = as.data.frame.list(cld_fs$category)

effect_analysis$cld_fs = cld_fs$Letters[match(effect_analysis$category, rownames(cld_fs))]

## ATP
atp_aov= aov(cube_Median_ATP ~ category, data = cube_effect_analysis)
summary(atp_aov)

# Pairwise 
#Clay/Silt vs. Fine sand - (0.028)
#Clay/Silt vs. Med/Coarse sand - not sig. (0.95)
#Fine sand vs. Med/Coarse sand - p = 0.0003
tukey_atp = TukeyHSD(atp_aov, "category")
cld_atp = multcompLetters4(atp_aov, tukey_atp)
cld_atp = as.data.frame.list(cld_atp$category)

effect_analysis$cld_atp = cld_atp$Letters[match(effect_analysis$category, rownames(cld_atp))]

## Parametric Boxplots

fs_cat = ggplot(effect_analysis, aes(x = category, y = Percent_Fine_Sand,  fill = category)) +
  geom_boxplot() +
  geom_text(aes(label = cld_fs, y = 85)) + 
  theme_bw() +
  ylab("Fine Sand (%)") + 
  xlab("") +
  theme(legend.title = element_blank())

atp_cat = ggplot(effect_analysis, aes(x = category, y = Median_ATP_picomoles_per_g,  fill = category)) +
  geom_boxplot() +
  geom_text(aes(label = cld_atp, y = 325))+ 
  theme_bw() +
  ylab(expression("ATP (pmol g"^-1*")")) + 
  xlab("") +
  theme(legend.title = element_blank())

# Only ATP and FS for final model, other pairwise tests below

## Final Conceptual Figure?? ####
d50_plot = d50_plot +
  ylab(expression("Effect Size Respiration (mg kg"^-1*")")) +
  theme(axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10), 
        legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8))

fs_cat = fs_cat + 
  theme(legend.position = "none", aspect.ratio = 1,
        plot.margin = unit(c(0, -1, 0, 0), "cm")) 

atp_cat_save = atp_cat +
  theme(legend.position = "none", aspect.ratio = 1,
        plot.margin = unit(c(0, 0, 0, -1), "cm"))


d50_box = ggarrange(d50_plot, labels = c("A"), nrow = 2, ggarrange(fs_cat, atp_cat_save, ncol = 2, widths = c(3,3), labels = c("B", "C"), hjust = -5, align = "h"))

ggsave("./Physical_Manuscript_Figures/d50_boxes.png", width = 12, height = 9, plot = d50_box)

## Continue Pairwise Comparisons ####

## Effect Fe 
fe_aov= aov(cube_Fe_Effect ~ category, data = cube_effect_analysis)
summary(fe_aov)

#check variance
res = bartlett.test(cube_Fe_Effect ~ category, data = cube_effect_analysis)
res

# Pairwise 
#Clay/Silt vs. Fine sand - (0.18)
#Clay/Silt vs. Med/Coarse sand - not sig. (0.0008)
#Fine sand vs. Med/Coarse sand - p = 0.0014
tukey_fe = TukeyHSD(fe_aov, "category")
cld_fe = multcompLetters4(fe_aov, tukey_fe)
cld_fe = as.data.frame.list(cld_fe$category)

effect_analysis$cld_fe = cld_fe$Letters[match(effect_analysis$category, rownames(cld_fe))]

## TN
tn_aov= aov(cube_TN ~ category, data = cube_effect_analysis)
summary(tn_aov)

# Pairwise 
#Clay/Silt vs. Fine sand - (0.84)
#Clay/Silt vs. Med/Coarse sand - not sig. (0.11)
#Fine sand vs. Med/Coarse sand - p = 0.06
tukey_tn = TukeyHSD(tn_aov, "category")
cld_tn = multcompLetters4(tn_aov, tukey_tn)
cld_tn = as.data.frame.list(cld_tn$category)

effect_analysis$cld_tn = cld_tn$Letters[match(effect_analysis$category, rownames(cld_tn))]

## TOC (non-equal variances)
toc_aov= aov(cube_TOC ~ category, data = cube_effect_analysis)
summary(toc_aov)

# Pairwise 
#Clay/Silt vs. Fine sand - (0.88)
#Clay/Silt vs. Med/Coarse sand - not sig. (0.004)
#Fine sand vs. Med/Coarse sand - p = 0.00007
tukey_toc = TukeyHSD(toc_aov, "category")
cld_toc = multcompLetters4(toc_aov, tukey_toc)
cld_toc = as.data.frame.list(cld_toc$category)

effect_analysis$cld_toc = cld_toc$Letters[match(effect_analysis$category, rownames(cld_toc))]

## SpC
spc_aov= aov(cube_Median_SpC ~ category, data = cube_effect_analysis)
summary(spc_aov)

# Pairwise 
#Clay/Silt vs. Fine sand - (0.36)
#Clay/Silt vs. Med/Coarse sand - not sig. (0.06)
#Fine sand vs. Med/Coarse sand - p = 0.38
tukey_spc = TukeyHSD(spc_aov, "category")
cld_spc = multcompLetters4(spc_aov, tukey_spc)
cld_spc = as.data.frame.list(cld_spc$category)

effect_analysis$cld_spc = cld_spc$Letters[match(effect_analysis$category, rownames(cld_spc))]

## Dry Grav
grav_aov= aov(cube_Final_Dry_Grav ~ category, data = cube_effect_analysis)
summary(grav_aov)

# Pairwise 
#Clay/Silt vs. Fine sand - (0.6)
#Clay/Silt vs. Med/Coarse sand - not sig. (0.1)
#Fine sand vs. Med/Coarse sand - p = 0.2
tukey_grav = TukeyHSD(grav_aov, "category")
cld_grav = multcompLetters4(grav_aov, tukey_grav)
cld_grav = as.data.frame.list(cld_grav$category)

effect_analysis$cld_grav = cld_grav$Letters[match(effect_analysis$category, rownames(cld_grav))]
