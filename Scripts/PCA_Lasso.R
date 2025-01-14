#### Sensitivity Analysis For ECA removals ####

# This script makes figures for ECA physical manuscript and performs PCA Analysis and LASSO regression

library(tidyverse)
library(corrplot)
library(ggpubr)
library(ggpmisc)
library(factoextra)
library(stringr)
library(glmnet)
library(magick)
library(tidytext)

rm(list=ls());graphics.off()

## Set image export

print.images = T

# Functions ---------------------------------------------------------------

# Transformation for normalization is cube root - have to cube root then add sign back to value to make it positive or negative
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

# Read in/Merge Data ------------------------------------------------------------

## Individual Rate data for histograms
all_data = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/EC_Sediment_SpC_pH_Temp_Respiration.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material))

cube_respiration = all_data %>% 
  select(c(Sample_Name, Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999) %>% 
  mutate(cube_Respiration_mg_kg = cube_root(abs(as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)))) %>% 
  mutate(Treat = if_else(grepl("D", Sample_Name), "Dry", "Wet"))


## Effect Size Data

#change this to published data when ready

effect = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/EC_Sediment_Effect_Size.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material, Methods_Deviation, Effect_Size_Initial_Gravimetric_Moisture_g_per_g,Effect_Size_Final_Gravimetric_Moisture_g_per_g))
  
## Read in Median Data to get Dry moisture values

grav = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/EC_Sediment_Sample_Data_Summary.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(grepl("D", Sample_Name)) %>%  
  mutate(Sample_Name = str_replace(Sample_Name, "-D", "_all")) %>% 
  mutate(Dry_Initial_Grav = as.numeric( Median_62948_Initial_Gravimetric_Moisture_g_per_g)) %>% 
  mutate(Dry_Final_Grav = as.numeric( Median_62948_Final_Gravimetric_Moisture_g_per_g)) %>% 
  mutate(Dry_Lost_Grav = Dry_Initial_Grav - Dry_Final_Grav) %>% 
  select(c(Sample_Name, Dry_Initial_Grav, Dry_Final_Grav, Dry_Lost_Grav)) 

## Read in grain size/ssa variables

grain = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/v3_CM_SSS_Sediment_Sample_Data_Summary.csv", skip = 2) %>% 
  filter(grepl("CM", Sample_Name)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "Sediment", "all")) %>% 
  select(c(Sample_Name, Percent_Tot_Sand, Percent_Coarse_Sand, Percent_Med_Sand, Percent_Fine_Sand, Percent_Silt, Percent_Clay, Mean_Specific_Surface_Area_m2_per_g))

## Join all data

effect_data = left_join(effect, grain, by = "Sample_Name") %>% 
  left_join(grav, by = "Sample_Name") %>% 
  mutate_at(vars(Effect_Size_SpC_microsiemens_per_cm:Dry_Lost_Grav), as.numeric)  # make data numeric 
  
# Transform Data ----------------------------------------------------------

## Cube PCA for LASSO####
# Fe outlier not in analysis

cube_effect_data = effect_data %>% 
  mutate(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = abs(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) %>% # made these abs. in original data analysis
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>% 
  column_to_rownames("Sample_Name") %>% 
  select(-matches("per_L")) %>% # remove samples with per_L in sample name, we're using mg/kg values
  rename(cube_Effect_SpC = cube_Effect_Size_SpC_microsiemens_per_cm) %>% 
  rename(cube_Effect_pH = cube_Effect_Size_pH) %>%
  rename(cube_Effect_Temp = cube_Effect_Size_Temperature_degC) %>%
  rename(cube_Effect_Respiration_mg_kg = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H) %>%
  rename(cube_Effect_Fe_mg_kg= cube_Effect_Size_Fe_mg_per_kg) %>%
  rename(cube_Effect_ATP_pmol_g = cube_Effect_Size_ATP_picomoles_per_g) %>%
  rename(cube_Effect_NPOC_mg_kg = cube_Effect_Size_Extractable_NPOC_mg_per_kg) %>% 
  rename(cube_Effect_TN_mg_kg = cube_Effect_Size_Extractable_TN_mg_per_kg) %>% 
  rename(cube_Effect_TOC_percent = cube_Effect_Size_C_percent_per_mg) %>%
  rename(cube_Effect_TN_percent = cube_Effect_Size_N_percent_per_mg) %>% 
  rename(cube_SSA = cube_Mean_Specific_Surface_Area_m2_per_g) %>% 
  relocate(cube_Effect_Respiration_mg_kg, .before = cube_Effect_SpC) %>% 
  filter(cube_Effect_Fe_mg_kg > -1) # remove Fe outlier for analysis


# LASSO from PCA Loadings-----------------------------------------------------------------

## Run PCA without effect size 
cube_data_pca = cube_effect_data %>% 
  select(-c(cube_Effect_Respiration_mg_kg))

## Scale and Center Data
cube_effect_pca <- prcomp(cube_data_pca, scale = TRUE,
                          center = TRUE, retx = T)

summary(cube_effect_pca)

# Plot cumulative variation
cumulative_var = cumsum(cube_effect_pca$sdev^2)/sum(cube_effect_pca$sdev^2)

pca_plot = data.frame(PC = 1:length(cube_effect_pca$sdev), Cumulative_Variance = cumulative_var)

if (print.images == T){
  
  png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cubed_Root_PCA_Axes.png"), width = 8, height = 8, units = "in", res = 300)

ax = ggplot(pca_plot, aes(x = PC, y = Cumulative_Variance)) +
  geom_col()

ax

}
dev.off()

# Look at PCA loadings
loadings = cube_effect_pca$rotation

print(loadings)

loadings_df <- as.data.frame(loadings) %>% 
  select(c(PC1, PC2, PC3, PC4, PC5)) # decided to use PC1 - 5, ~70% of explanatory power

# Generate heatmap
row_names <- rownames(loadings_df)
loadings_df$Variable <- row_names

# Melt the dataframe for plotting
loadings_melted <- reshape2::melt(loadings_df, id.vars = "Variable")

if (print.images == T) {

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cubed_Root_PCA_Heatmap.png"), width = 8, height = 8, units = "in", res = 300)

# Plotting the heatmap
ggplot(loadings_melted, aes(x = variable, y = Variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Principal Component", y = "Variable", fill = "Loadings") +
  ggtitle("PCA Loadings Heatmap")+
  geom_text(aes(label = round(value, 2)), color = "black", size = 3)

}

dev.off()

## LASSO with PCA Loadings ####
set.seed(42)
## Set response variable 
yvar <- data.matrix(scale(cube_effect_data$cube_Effect_Respiration_mg_kg, center = TRUE, scale = TRUE))
#yvar <- data.matrix(cube_effect_data_corr$Cube_Effect_Size) #tried looking at differences without scaling
mean(yvar) #check scaling
sd(yvar)

## Set predictor variables (PCA Loadings)
pca_scores = as.data.frame(scale(cube_effect_pca$x, center = TRUE, scale = TRUE))
#pca_scores = as.data.frame(cube_effect_pca$x) #tried looking at differences without scaling
mean(pca_scores$PC1)
sd(pca_scores$PC1)

# only using first 5 PC's, which explain ~70% of data
xvars <- data.matrix(pca_scores[, c("PC1",  "PC2", "PC3", "PC4", "PC5")])

lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                  ,standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                  #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                 # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
                  )

best_lambda <- lasso$lambda.min
best_lambda

plot(lasso)

best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
   , standardize = FALSE, standardize.response = FALSE, intercept = FALSE
#  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
  #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

lasso_coefs = coef(best_lasso_model)

yvar_predict <- predict(best_lasso_model, s = best_lambda, newx = xvars)

sst <- sum((yvar - mean(yvar))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst

rsq

lasso_df = as.data.frame(as.matrix(lasso_coefs))

colnames(lasso_df) = c("Coefficients")

lasso_df = lasso_df %>% 
  rownames_to_column(var = "variable") %>% 
  slice(-1) 

lasso_df$y = 0

## Make heatmap of LASSO coefficients

if (print.images == T) {
  
  png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_PCA_LASSO_Heat_Matrix.png"), width = 8, height = 8, units = "in", res = 300)
  
ggplot(lasso_df, aes(variable, y)) +
  geom_tile(fill = "white", color = "black") +
  geom_text(aes(label = round(Coefficients, 2), color = Coefficients), size = 3, fontface = "bold") + 
  scale_color_gradient2(high = "#2166ac", low = "#b2182b", mid = "white", 
                       midpoint = 0, limit = c(-1, 1),
                       space = "Lab", name="Coefficient") +
  theme_minimal() + 
  theme(aspect.ratio = 0.25, 
        axis.text.x = element_text(angle = 90, 
                                   hjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 0), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank()) +
  scale_x_discrete(position = "top") +
  labs(y = "LASSO using PCA Loadings")

}

dev.off()


weighted_pc = merge(lasso_df, loadings_melted, by = "variable") %>% 
  mutate(pc_weight = value * Coefficients)

weighted_total = weighted_pc %>% 
  group_by(Variable) %>% 
  mutate(total_weight = sum(pc_weight))

abs_total = weighted_total %>% 
  mutate(across(where(is.numeric), abs)) # abs values

if (print.images == T) {
  
  png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_PCA_LASSO_parsed_weights.png"), width = 8, height = 8, units = "in", res = 300)
  
  ggplot(weighted_total, aes(x = fct_reorder(Variable, total_weight))) +
    geom_col(aes(y = pc_weight, fill = variable)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  }

dev.off()

if (print.images == T) {
  
  png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_PCA_LASSO_total_weighted.png"), width = 8, height = 8, units = "in", res = 300)
  
  ggplot(weighted_total, aes(x = fct_reorder(Variable, total_weight), y = total_weight)) + 
    geom_col() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
}

dev.off()


if (print.images == T) {
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_PCA_LASSO_abs_total_weighted.png"), width = 8, height = 8, units = "in", res = 300)
  
  ggplot(abs_total, aes(x = fct_reorder(Variable, total_weight), y = total_weight)) + 
    geom_col() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
}

dev.off()

cube_limits = c(-12, 12)

if (print.images == T){

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_PCA.png"), width = 8, height = 8, units = "in", res = 300)

cube_pca = fviz_pca_biplot(cube_effect_pca, col.var = "black",geom = "point"
) +
  geom_point(aes(color = cube_effect_data$cube_Effect_Respiration_mg_kg), size = 3.5) +
  scale_color_gradient2(limits = cube_limits, low = "firebrick2", mid = "goldenrod2",
                        high = "dodgerblue2", midpoint = (max(cube_limits)+min(cube_limits))/2) +
  labs(color = paste0("Cubed Root Wet - Dry Rate"))+ theme(legend.position = c(0.2, 0.85), 
                                                           legend.key.size = unit(0.25, "in"), 
                                                           legend.title = element_text(size = 8),
                                                           axis.title.x = element_text(size = 10))

cube_pca

}

dev.off()

# LASSO from downselected Pearson Correlation Variables  ------------------

## Pearson Correlation Matrix ####

# scale data before it goes into correlation matrix
scale_cube_effect = as.data.frame(scale(cube_effect_data))

scale_cube_effect_pearson <- cor(scale_cube_effect, method = "pearson")

if (print.images == T){

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_Pearson_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

corrplot(scale_cube_effect_pearson,type = "upper", method = "number", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "Effect Samples Pearson Correlation")
 
}

dev.off()

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

color_palette <- colorRampPalette(c("#B2182B", "#F7F7F7", "#2166AC"))(200)

if (print.images == T){

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_Pearson_Correlation_Matrix_One_Line.png"), width = 12, height = 5, units = "in", res = 300)
 

corrplot(corr_effect, type = "upper", method = "number", tl.col = "black", tl.cex = 1.6, cl.cex = 1, diag = FALSE, is.corr = FALSE, cl.pos = 'n', col = color_palette)

}

dev.off()

## Downselected LASSO - Loop through coefficients to choose for LASSO ####

# 1) Pivot data frame and sort highest to lowest

pearson_df <- as.data.frame(scale_cube_effect_pearson)

row_names_pearson <- rownames(pearson_df)

pearson_df$Variable <- row_names_pearson

# Melt the dataframe for plotting
pearson_melted <- reshape2::melt(pearson_df, id.vars = "Variable") %>% 
  filter(value != 1) %>% 
  mutate(value = abs(value)) %>% 
  filter(!grepl("Respiration", Variable))

effect_melted <- pearson_melted %>% 
  filter(grepl("Respiration", variable)) %>%
  filter(!grepl("Silt", Variable)) # remove silt from variables, lots of 0 values so not using

choose_melted <- pearson_melted %>% 
  filter(!grepl("Respiration", variable)) %>%
 filter(!grepl("Silt", variable)) %>%
 filter(!grepl("Silt", Variable)) %>% #try removing silt (0 values)
  #distinct(value, .keep_all = TRUE) %>% 
  left_join(effect_melted, by = "Variable") %>% 
  rename(Variable_1 = Variable) %>% 
  rename(Variable_2 = variable.x) %>% 
  rename(Correlation = value.x) %>% 
  rename(Variable_1_Effect_Correlation = value.y) %>% 
  select(-c(variable.y)) %>% 
  left_join(effect_melted, by = c("Variable_2" = "Variable")) %>% 
  rename(Variable_2_Effect_Correlation = value) %>% 
  select(-c(variable))

loop_melt = choose_melted %>% 
  arrange(desc(Correlation))

# Pearson correlation coefficient to remove above
correlation = 0.7

## Start loop to remove highly correlated (> 0.5)
effect_filter = function(loop_melt) {
  
  rows_to_keep = rep(TRUE, nrow(loop_melt))
  
  for (i in seq_len(nrow(loop_melt))) {
    
    if (!rows_to_keep[i]) next
    
    row = loop_melt[i, ]
    
    if (row$Correlation < correlation) next
    
    if(row$Variable_1_Effect_Correlation >= row$Variable_2_Effect_Correlation) {
      
      var_to_keep = row$Variable_1
      var_to_remove = row$Variable_2
      
    } else {
      
      var_to_keep = row$Variable_2
      var_to_remove = row$Variable_1
      
    }
    
    loop_melt$Variable_to_Keep[i] = var_to_keep
    loop_melt$Variable_to_Remove[i] = var_to_remove
   
    for (j in seq(i + 1, nrow(loop_melt))) {
      
      if(loop_melt$Variable_1[j] == var_to_remove || loop_melt$Variable_2[j] == var_to_remove) {
        
        rows_to_keep[j] = FALSE
        
      }
      
    }
     
    
  }
  
  return(loop_melt[rows_to_keep, ])
  
}
  
filtered_data = effect_filter(loop_melt) 

# pull out variables to remove
removed_variables = filtered_data %>% 
  distinct(Variable_to_Remove)

# pull out all variables 
all_variables = effect_melted %>% 
  select(c(Variable))

# remove variables from all variables to get variables to keep for LASSO 
kept_variables = effect_melted[!(effect_melted$Variable %in% removed_variables$Variable_to_Remove), ] #keeps SpC, Temp, pH, ATP, NPOC, TN (ext), TOC, TN (solid), med sand, silt

# if silt is removed, keeps SpC, Temp, pH, ATP, NPOC, TOC, TN (solid), med sand, fine sand, lost grav. moisture

## LASSO VARIABLES ####

# Keep variables selected from down-selected correlation matrix and add Cube_Effect_Size
col_to_keep = unique(kept_variables$Variable)
col_to_keep = c(col_to_keep, "cube_Effect_Respiration_mg_kg")

cube_variables = cube_effect_data[, col_to_keep, drop = FALSE]

## LASSO with Correlation Matrix Selected Variables ####
set.seed(42)
## Set response variable (Cube_Effect_Size) and scale
yvar <- data.matrix(scale(cube_variables$cube_Effect_Respiration_mg_kg, center = TRUE, scale = TRUE))
mean(yvar)
sd(yvar)

## Set predictor variables and scale
exclude_col = "cube_Effect_Respiration_mg_kg"

x_cube_variables = as.data.frame(scale(cube_variables[, !(names(cube_variables) %in% exclude_col)], center = T, scale = T))
#mean(x_cube_variables$Cube_SpC_Diff)
#sd(x_cube_variables$Cube_SpC_Diff)

xvars <- data.matrix(x_cube_variables)

lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                  ,standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                  #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                  # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

best_lambda <- lasso$lambda.min
best_lambda

plot(lasso)

best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
                           , standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                           #  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                           #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

#lasso_coefs = coef(best_lasso_model)
ds_lasso_coefs = coef(best_lasso_model)

yvar_predict <- predict(best_lasso_model, s = best_lambda, newx = xvars)

sst <- sum((yvar - mean(yvar))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst

rsq

ds_lasso_df = as.data.frame(as.matrix(ds_lasso_coefs))

colnames(ds_lasso_df) = c("Coefficients")

ds_lasso_df = ds_lasso_df %>% 
  rownames_to_column(var = "variable") %>% 
  slice(-1) 

ds_lasso_df$y = "LASSO"

color_palette <- colorRampPalette(c("#B2182B", "#F7F7F7", "#2166AC"))(200)

if (print.images == T) {
  
  png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_LASSO_Heat_Matrix.png"), width = 12, height = 4, units = "in", res = 300)
  
  ggplot(ds_lasso_df, aes(variable, y)) +
    geom_tile(fill = "white", color = "black") +
    geom_text(aes(label = round(Coefficients, 2), color = Coefficients), size = 8, fontface = "bold") + 
    scale_color_gradientn(colors = color_palette, 
      limit = c(-1, 1),
      guide = "none") +
    theme_minimal() + 
    theme(aspect.ratio = 0.1, 
          axis.text.x = element_text(angle = 90, hjust = 0, face = "bold", size = 13), 
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0), 
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank()) +
    scale_x_discrete(position = "top") +
    labs(y = "LASSO Coefficients")
  
}

dev.off()

#### Combined Heat maps ####

lasso_pear_df = bind_rows(ds_lasso_df, corr_effect_df) %>% 
  rename(type = y)

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Combined_Heat_Matrix.png"), width = 12, height = 4, units = "in", res = 300)

ggplot(lasso_pear_df, aes(variable, type)) +
  geom_tile(fill = "white", color = "black") +
  geom_text(aes(label = round(Coefficients, 2), color = Coefficients), size = 6, fontface = "bold") + 
  scale_color_gradientn(colors = color_palette, 
                        limit = c(-1, 1),
                        guide = "none") +
  theme_bw() + 
  theme(aspect.ratio = 0.1, 
        axis.text.x = element_text(angle = 90, hjust = 0, size = 13), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())+#, 
        #axis.text.y = element_blank()) +
  scale_x_discrete(position = "top") #+
 # labs(y = "LASSO Coefficients")

dev.off()

## LASSO with all variables to check for collinearity effects ####

## Set response variable (Cube_Effect_Size) and scale
yvar <- data.matrix(scale(cube_effect_data$cube_Effect_Respiration_mg_kg, center = TRUE, scale = TRUE))
mean(yvar)
sd(yvar)

## Set predictor variables and scale
x_cube_variables =  as.data.frame(scale(cube_effect_data[, !(names(cube_effect_data) %in% exclude_col)], center = T, scale = T))
#mean(x_cube_variables$Cube_SpC_Diff)
#sd(x_cube_variables$Cube_SpC_Diff)

xvars <- data.matrix(x_cube_variables)

lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                  ,standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                  #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                  # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

best_lambda <- lasso$lambda.min
best_lambda

plot(lasso)

best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
                           , standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                           #  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                           #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
)

lasso_coefs = coef(best_lasso_model)

yvar_predict <- predict(best_lasso_model, s = best_lambda, newx = xvars)

sst <- sum((yvar - mean(yvar))^2)
sse <- sum((yvar_predict - yvar)^2)

rsq = 1 - sse/sst

rsq


## calculate VIF

model = lm(Cube_Effect_Size ~ .,  data = cube_effect_data_corr)

vif_values = vif(model)
vif_values

model_corr = lm(Cube_Effect_Size ~ Cube_SpC_Diff + Cube_pH_Diff + Cube_Temp_Diff + Cube_Fe_mg_kg_Diff + Cube_ATP_pmol_g_Diff + Cube_NPOC_mg_C_per_kg_Diff + Cube_TN_mg_N_per_kg_Diff + Cube_TOC_Percent_Diff + Cube_TN_Percent_Diff + Cube_Fine_Sand + Cube_Dry_LostGravMoi, data = cube_effect_data_corr)

vif_corr = vif(model_corr)
vif_corr

# Make Figures ------------------------------------

## Histogram of all Rates

if (print.images == T) {
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_All_Rates_Cubed_Histogram.png"), width = 4, height = 4, units = "in", res = 300)

all_cube_hist = ggplot(cube_respiration, aes(x = cube_Respiration_mg_kg)) +
  geom_histogram(position = "identity", alpha = 0.8, aes(fill = Treat))+
  scale_fill_manual(values = c("#D55E00","#0072B2"))  +
  #ggtitle("Wet Rates")+
  theme(strip.text = element_text(
    size = 4))+
  ylim(0, 87.5)+
  theme_bw() + 
  theme(legend.position = c(0.85, 0.8), 
        legend.key.size = unit(0.15, "in"), 
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 10)) +
  guides(fill = guide_legend(title="Treatment")) + 
  xlab(expression(atop("\n Cubed Root Respiration Rate", "(mg O"[2]*" kg"^-1*" H"^-1*")"))) +
  ylab("Count")

all_cube_hist

}
dev.off()

# Histogram of effect size

cube_effect_limits <- c(-12, 12)

if (print.images == T) {
png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_Histogram.png"), width = 4, height = 4, units = "in", res = 300)

cube_effect_hist = ggplot(cube_effect_data, aes(x = cube_Effect_Respiration_mg_kg))+
  # geom_histogram(binwidth = 0.15, fill = "#009E73")+
  geom_histogram(binwidth = 0.5, aes(fill = after_stat(x))) +
  scale_fill_gradient2(name = "Cubed Root Effect Size", limits = cube_effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(cube_effect_limits)+min(cube_effect_limits))/2) +
  theme_bw()+
  #theme(axis.title.x = element_text(size = 4),
  #  axis.title.y = element_text(size = 4),
  #  axis.text.x = element_text(size = 4),
  #  axis.text.y = element_text(size =4))+
  xlim(c(-12, 12))+
  ylab("Count\n")+
  theme(legend.position = c(0.8, 0.8),
        legend.key.size = unit(0.15, "in"), 
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 10)) + 
  xlab(expression(atop("\n Cubed Root Effect Size", "(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" H"^-1*")"))) 


cube_effect_hist
}

dev.off()


# Combined Map/Rates/Effect Size histograms - Map being made own figure
#map_path = "C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/Map/Ecoregion_EffectSize_Map.pdf"

#map_image = image_read_pdf(map_path, density = 300)

#scale_map_image = image_scale(map_image, "80%")

#scale_map_label_image = image_annotate(scale_map_image, "C", size = 15, location = "+10+80", color = "black")

# Read in Rate histograms figure .png

rate_image = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-08-29_All_Rates_Cubed_Histogram.png")

rate_label_image = image_annotate(rate_image, "A", size = 65, location = "+30+20", color = "black")

# Read in Effect Size Figure 

effect_image = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-08-29_Cube_Median_Effect_Histogram.png")

effect_label_image = image_annotate(effect_image, "B", size = 65, location = "+30+20", color = "black")

com_image = image_append(c(rate_label_image, effect_label_image))

#com_map_image = image_append(c(com_image, scale_map_label_image), stack = TRUE)

image_write(com_image, path = "C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-08-29_Combined_Rates.png")


#Fe, Fine Sand most important

# plot with Fe outlier, but don't perform stats with Fe outlier
cube_effect_data_fe_inc = effect_data %>% 
  mutate(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = abs(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) %>% # made these abs. in original data analysis
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>% 
  column_to_rownames("Sample_Name") %>% 
  select(-matches("per_L")) %>% # remove samples with per_L in sample name, we're using mg/kg values
  rename(cube_Effect_SpC = cube_Effect_Size_SpC_microsiemens_per_cm) %>% 
  rename(cube_Effect_pH = cube_Effect_Size_pH) %>%
  rename(cube_Effect_Temp = cube_Effect_Size_Temperature_degC) %>%
  rename(cube_Effect_Respiration_mg_kg = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H) %>%
  rename(cube_Effect_Fe_mg_kg= cube_Effect_Size_Fe_mg_per_kg) %>%
  rename(cube_Effect_ATP_pmol_g = cube_Effect_Size_ATP_picomoles_per_g) %>%
  rename(cube_Effect_NPOC_mg_kg = cube_Effect_Size_Extractable_NPOC_mg_per_kg) %>% 
  rename(cube_Effect_TN_mg_kg = cube_Effect_Size_Extractable_TN_mg_per_kg) %>% 
  rename(cube_Effect_TOC_percent = cube_Effect_Size_C_percent_per_mg) %>%
  rename(cube_Effect_TN_percent = cube_Effect_Size_N_percent_per_mg) %>% 
  rename(cube_SSA = cube_Mean_Specific_Surface_Area_m2_per_g) %>% 
  relocate(cube_Effect_Respiration_mg_kg, .before = cube_Effect_SpC)


if (print.images == T) {
  
  png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_vs_Fe_Scatter.png"), width = 6, height = 6, units = "in", res = 300)
  
  fe_cube = ggplot(cube_effect_data_fe_inc, aes(x = cube_Effect_Fe_mg_kg, y = cube_Effect_Respiration_mg_kg)) +
    geom_point() +
    theme_bw() +
    #stat_cor(data = fe_cube_out, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")))+
    stat_cor(data = cube_effect_data, label.x = -2.5, label.y = 11, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
    stat_cor(data = cube_effect_data, label.x = -2.5, label.y = 10.25, size = 4, digits = 2, aes(label = paste(..p.label..)))+
    stat_poly_line(data = cube_effect_data, se = FALSE)+
    scale_y_continuous(name="Cubed Root Effect Size (mg/kg) (Wet - Dry)", limits=c(0, 11.25))+
    xlab("Cubed Root Fe (II) (mg/kg) Difference (Wet - Dry)") +
    #ylab("Cubed Root Effect Size (mg/kg) (Wet - Dry)")+ 
    theme(text = element_text(size = 12)) 
  
  fe_cube
  
}
dev.off()

if (print.images == T){ 
  png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Cube_Median_Effect_vs_Fine_Sand_Scatter.png"), width = 6, height = 6, units = "in", res = 300)
  
  fs_cube = ggplot(cube_effect_data, aes(x = cube_Percent_Fine_Sand, y = cube_Effect_Respiration_mg_kg)) +
    geom_point() +
    theme_bw() +
    #stat_cor(data = cube_effect_data_corr, size = 5, digits = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`\n`~")))+ #sep = "~`;`~"
    stat_cor(data = cube_effect_data, label.x = 0.9, label.y = 11.25, size = 4, digits = 2, aes(label = paste(..rr.label..)))+
    stat_cor(data = cube_effect_data, label.x = 0.9, label.y = 10.5, size = 4, digits = 2, aes(label = paste(..p.label..)))+
    stat_poly_line(data = cube_effect_data, se = FALSE)+
    scale_y_continuous(name="Cubed Root Effect Size (mg/kg) (Wet - Dry)", limits=c(0, 11.25))+
    xlab("Cubed Root Fine Sand (%)") +
    #ylab("Cubed Root Effect Size (mg/kg) (Wet - Dry)")+ 
    theme(text = element_text(size = 12)) 
  
  fs_cube
}

dev.off()

## Make combined one line heat maps (Pearson + Downselected LASSO)

combine_hm_image = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-08-30_Combined_Heat_Matrix.png")

combine_label_image = image_annotate(combine_hm_image, "A", size = 100, location = "+25+50", color = "black")

#lasso_hm_image = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-08-29_LASSO_Heat_Matrix.png")

#lasso_label_image = image_annotate(lasso_hm_image, "B", size = 100, location = "+100+100", color = "black")

#hm_combine_image = image_append(c(pearson_label_image, lasso_label_image), stack = TRUE)

fine_sand_image = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-08-30_Cube_Median_Effect_vs_Fine_Sand_Scatter.png")

fine_sand_label_image = image_annotate(fine_sand_image, "B", size = 65, location = "+30+20", color = "black")

fine_sand_scale_image = image_scale(fine_sand_label_image, "100%")

fe_image = image_read("C:/GitHub/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-08-30_Cube_Median_Effect_vs_Fe_Scatter.png")

fe_label_image = image_annotate(fe_image, "C", size = 65, location = "+30+20", color = "black")

fe_scale_image = image_scale(fe_label_image, "100%")

scatter_image = image_append(c(fine_sand_scale_image, fe_scale_image))

whole_image = image_append(c(combine_label_image, scatter_image), stack = TRUE)

image_write(whole_image, path = "C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/2024-08-30_Combined_Heat_Map.png")

## Bar Plots of things in LASSO colored by effect size

medians = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/EC_Sediment_Sample_Data_Summary.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material, Median_Missing_Reps, Median_Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  select(-matches("per_L")) %>% 
  mutate(Initial_Grav = as.numeric( Median_62948_Initial_Gravimetric_Moisture_g_per_g)) %>% 
  mutate(Final_Grav = as.numeric( Median_62948_Final_Gravimetric_Moisture_g_per_g)) %>% 
  mutate(Lost_Grav = Initial_Grav - Final_Grav) %>% 
  separate(Sample_Name, c("Sample_Name", "Rep"), sep = "-") %>% 
  mutate(Sample_Name = paste0(Sample_Name, "_all")) 

effect_size = effect_data %>%
  select(c(Sample_Name, Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) %>% 
  rename(Effect_Size = Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H) %>% 
  mutate(Effect_Size = abs(Effect_Size))

medians_effect = medians %>% 
  left_join(grain) %>% 
  left_join(effect_size) %>% 
  reshape2::melt(id.vars = c("Sample_Name",  "Rep", "Effect_Size"))

effect_limits  = c(-1400, 1400)

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Wet_Medians_by_Effect.png"), width = 12, height = 12, units = "in", res = 300)

medians_effect %>% 
  filter(Rep == "W") %>% 
  filter(!grepl("Lost_Grav|Initial_Grav|Final_Grav", variable)) %>% 
  mutate(value = as.numeric(value)) %>% 
  ggplot(aes(x = reorder_within(Sample_Name, value, variable), y = value, fill = Effect_Size)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +  
  scale_fill_gradient2(name = "Effect Size", limits = effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(effect_limits)+min(effect_limits))/2) +
  theme_bw()

dev.off()

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Dry_Medians_by_Effect.png"), width = 12, height = 12, units = "in", res = 300)

medians_effect %>% 
  filter(Rep == "D") %>% 
  #filter(!grepl("Lost_Grav|Initial_Grav|Final_Grav", variable)) %>% 
  mutate(value = as.numeric(value)) %>% 
  ggplot(aes(x = reorder_within(Sample_Name, value, variable), y = value, fill = Effect_Size)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +  
  scale_fill_gradient2(name = "Effect Size", limits = effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(effect_limits)+min(effect_limits))/2) +
  theme_bw() + 
  ggtitle("Dry Treatment Medians")

dev.off()

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_Medians_by_Effect.png"), width = 12, height = 12, units = "in", res = 300)

medians_effect %>% 
  #filter(!grepl("Lost_Grav|Initial_Grav|Final_Grav", variable)) %>% 
  mutate(value = as.numeric(value)) %>% 
  ggplot(aes(x = reorder_within(Sample_Name, value, variable), y = value, fill = Effect_Size)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +  
  scale_fill_gradient2(name = "Effect Size", limits = effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(effect_limits)+min(effect_limits))/2) +
  theme_bw() + 
  ggtitle("All Treatment Medians")

dev.off()

