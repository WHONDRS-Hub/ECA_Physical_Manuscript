# This script makes scatterplot figures for ECA physical manuscript
library(tidyverse); library(ggpubr)
rm(list=ls()); graphics.off()

# Functions ---------------------------------------------------------------
# Transformation for normalization is cube root - have to cube root then add sign back to value to make it positive or negative
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

## Read in/Merge Data for Histograms -------------------------------------
# Individual Rate data for histograms ####
data = read.csv("Data/EC_Sediment_Sample_Data_Summary.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material)) %>%
  select(Sample_Name, Median_Respiration_Rate_mg_DO_per_kg_per_H)
  

ev_data = read.csv("Data/WHONDRS_EV_Data_Package/Sample_Data/WHONDRS_EV_Sediment_Sample_Data_Summary.csv", skip = 2) %>%
  filter(grepl("EV", Sample_Name)) %>%
  select(-c(Field_Name, IGSN, Material)) %>%
  select(Sample_Name, Median_Respiration_Rate_mg_DO_per_kg_per_H)
  

all_data = rbind(data, ev_data) 

resp <- all_data %>%
  mutate(
    rate_mgkg_h = as.numeric(Median_Respiration_Rate_mg_DO_per_kg_per_H),
    treatment = str_extract(Sample_Name, "(?<=-)[WD]$"),
    site = str_remove(Sample_Name, "-[WD]$")
  ) %>%
  filter(!is.na(site), !is.na(treatment), !is.na(rate_mgkg_h))

# --- pair wet and dry samples by site ---
paired <- resp %>%
  select(site, treatment, rate_mgkg_h) %>%
  pivot_wider(
    names_from = treatment,
    values_from = rate_mgkg_h
  ) %>%
  filter(!is.na(W), !is.na(D)) %>%
  mutate(
    abs_diff = abs(W - D),
    slower_rate = pmin(abs(W),abs(D)),
    
    slower_treatment = case_when(
      abs(D) < abs(W) ~ "Dry",
      abs(W) < abs(D) ~ "Wet",
      TRUE ~ "Tie"
    ),
    
    rel_diff = abs_diff / (slower_rate)
  ) %>%
  filter(slower_treatment != "Tie")

paired2 = paired %>%
  filter(D< 0)
# --- plot ---
ggplot(paired2, aes(x = abs_diff, y = rel_diff, color = slower_treatment)) +
  geom_point(size = 3, alpha = 0.5, stroke = 1) +  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.8) +
  labs(
    x = "Absolute difference between wet and dry median respiration rates (mg/kg/h)",
    y = "Absolute difference / slower rate",
    color = "Slower rate from"
  ) +
  scale_color_manual(values = c("Wet" = "red", "Dry" = "dodgerblue2")) +
  theme_bw(base_size = 12)


library(ggrepel)

ggplot(paired, aes(x = abs_diff, y = rel_diff, color = slower_treatment, label = site)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(size = 3.5, max.overlaps = Inf) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.8) +
  labs(
    x = "Absolute difference between wet and dry median respiration rates (mg/kg/h)",
    y = "Absolute difference / slower rate",
    color = "Slower rate from",
    title = "Wet vs Dry respiration differences by site"
  ) +
  theme_bw(base_size = 12)

