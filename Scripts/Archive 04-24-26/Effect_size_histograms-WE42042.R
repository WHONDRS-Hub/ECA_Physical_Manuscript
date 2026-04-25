# This script makes effect size histograms from ECA and EV

library(tidyverse);library(corrplot);library(ggpubr);library(ggpmisc);library(factoextra);library(stringr);library(glmnet);library(magick); library(ggnewscale); library(FSA); library(multcompView); library(rcompanion)

rm(list=ls());graphics.off()

# current_path <- rstudioapi::getActiveDocumentContext()$path
# setwd(dirname(current_path))
# setwd("./..")
# getwd()

# Functions ---------------------------------------------------------------

# Transformation for normalization is cube root - have to cube root then add sign back to value to make it positive or negative
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

## Read in/Merge Data  ------------------------------------------------------------

## Effect Size Data for ECA ####

effect = read.csv("Data/EC_Sediment_Effect_Size.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%
 select(-c(IGSN, Field_Name, Material, Methods_Deviation)) %>%
  dplyr::select(Sample_Name,Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)

effect_ev = read.csv('Data/WHONDRS_EV_Data_Package/Sample_Data/WHONDRS_EV_Sediment_Effect_Size.csv', skip =2 ) %>%
  filter(grepl("EV", Sample_Name)) %>% 
  select(-c(IGSN, Field_Name, Material, Methods_Deviation)) %>%
  dplyr::select(Sample_Name,Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = Effect_Size_Median_Respiration_Rate_mg_DO_per_kg_per_H)

## Join all data with effect size ####
effect_data = rbind(effect,effect_ev)
effect_data$Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = as.numeric(effect_data$Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)  

# All data for respiration rates: 
all_data = read.csv("./Data/EC_Sediment_SpC_pH_Temp_Respiration.csv", skip = 2) %>% 
  filter(grepl("EC", Sample_Name)) %>% 
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  # remove samples with too much water (EC_011, EC_012), sample with no mg/kg (EC_023), duplicated NEON sites (EC_052, EC_053, EC_057)
  select(-c(Field_Name, IGSN, Material)) %>%
  select(Sample_Name,Respiration_Rate_mg_DO_per_kg_per_H) %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999)


ev_all = read.csv('Data/WHONDRS_EV_Data_Package/Sample_Data/WHONDRS_EV_Sediment_SpC_pH_Temp_Respiration.csv', skip = 2) %>%
  filter(grepl("EV", Sample_Name)) %>%
  select(Sample_Name,Respiration_Rate_mg_DO_per_kg_per_H) %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999)
  
all_data = rbind(all_data,ev_all)
all_data$Respiration_Rate_mg_DO_per_kg_per_H = as.numeric(all_data$Respiration_Rate_mg_DO_per_kg_per_H)

# Transform Data --------------------------
--------------------------------

cube_effect = effect_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  dplyr::rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>% 
  column_to_rownames("Sample_Name") 

cube_respiration = all_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% # cube root transform data
  dplyr::rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) %>% 
  column_to_rownames("Sample_Name") 

# Histograms --------------------------------------------------------------
all_hist = all_data %>% 
  filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999) %>% # removes overexposed samples, missing replicates
  mutate(Treat = if_else(grepl("D", Sample_Name), "Dry", "Wet")) %>% 
  ggplot(aes(x = (as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)))) +
  geom_histogram(position = "identity", alpha = 0.8, aes(fill = Treat))+
  scale_fill_manual(values = c("#D55E00","#0072B2"))  +
  theme(strip.text = element_text(
    size = 4))+
  #ylim(0, 87.5)+
  theme_bw() + 
  theme(legend.position = c(0.85, 0.8), 
        legend.key.size = unit(0.15, "in"), 
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 10)) +
  guides(fill = guide_legend(title="Treatment")) + 
  xlab(expression(atop("\n Respiration Rate", "(mg O"[2]*" kg"^-1*" h"^-1*")"))) +
  ylab("Count")

all_hist

effect_limits <- c(-1400, 1400)

effect_hist = ggplot(effect_data, aes(x = (Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)))+
  geom_histogram(binwidth = 25, aes(fill = after_stat(x))) +
  scale_fill_gradient2(name = "Effect Size", limits = effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(effect_limits)+min(effect_limits))/2) +
  theme_bw()+
  xlim(c(-1500, 1500))+
  ylab("Count\n")+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 10)) + 
  xlab(expression(atop("\n Effect Size Respiration Rate","(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" h"^-1*")")))

effect_hist

combined_SI_hist = ggarrange(all_hist, effect_hist, ncol = 2, labels = c("A", "B"), hjust = -1, vjust = 2.5)

#ggsave("./Figures/Combined_SI_Histograms.png", plot = combined_SI_hist, width = 10, height = 5, dpi = 300)

## Cube Root Transformed
# These used with Nathan Johnson annotations in manuscript

cube_respiration = cube_respiration %>%
  mutate(Treat = if_else(grepl("D", Sample_Name), "Dry", "Wet")) 
  
all_cube_hist = ggplot(cube_respiration, aes(x = (cube_Respiration_Rate_mg_DO_per_kg_per_H))) +
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
    xlab(expression(atop("\n Respiration Rate" ^(1/3)*"", "(mg O"[2]*" kg"^-1*" h"^-1*")"))) +
    ylab("Count")
  
  all_cube_hist

  # Histogram of effect size

cube_effect_limits <- c(-12, 12)

cube_effect_hist = ggplot(cube_effect, aes(x = (cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)))+
    geom_histogram(binwidth = 0.5, aes(fill = after_stat(x))) +
    scale_fill_gradient2(name = "Cubed Root Effect Size", limits = cube_effect_limits, low = "firebrick2", mid = "goldenrod2",
                         high = "dodgerblue2", midpoint = (max(cube_effect_limits)+min(cube_effect_limits))/2) +
    theme_bw()+
    xlim(c(-12, 12))+
    ylab("Count\n")+
    theme(legend.position = "none",
          axis.title.x = element_text(size = 10)) + 
    xlab(expression(atop("\n Effect Size Respiration Rate" ^(1/3)*"","(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" h"^-1*")")))
  
cube_effect_hist

combined_hist = ggarrange(all_cube_hist, cube_effect_hist, ncol = 2, labels = c("A", "B"), hjust = -5, vjust = 2.5)

#ggsave("./Physical_Manuscript_Figures/Combined_Cube_Histograms.pdf", plot = combined_hist, width = 10, height = 5, dpi = 300)

