# This script makes scatterplot figures for ECA physical manuscript
library(tidyverse)
library(ggpubr)
library(ggrepel)

rm(list = ls()); graphics.off()

# Functions ---------------------------------------------------------------
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

## Read in / merge sediment data ------------------------------------------
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

## Read in field metadata to get perennial/intermittent -------------------
ec_meta = read.csv("Data/EC_Field_Metadata.csv", fileEncoding = "latin1") %>%
  select(
    Parent_ID,
    Intermittent_or_Perennial
  ) %>%
  rename(site = Parent_ID)

ev_meta = read.csv("Data/WHONDRS_EV_Data_Package/WHONDRS_EV_Field_Metadata.csv") %>%
  select(
    Parent_ID,
    Intermittent_or_Perennial
  ) %>%
  rename(site = Parent_ID)

site_meta = bind_rows(ec_meta, ev_meta) %>%
  mutate(
    flow_regime = case_when(
      str_detect(Intermittent_or_Perennial, regex("perennial", ignore_case = TRUE)) ~ "Perennial",
      str_detect(Intermittent_or_Perennial, regex("intermittent", ignore_case = TRUE)) ~ "Non-perennial",
      TRUE ~ NA_character_
    )
  ) %>%
  distinct(site, flow_regime)

## Prepare paired wet/dry data --------------------------------------------
resp <- all_data %>%
  mutate(
    rate_mgkg_h = as.numeric(Median_Respiration_Rate_mg_DO_per_kg_per_H),
    treatment = str_extract(Sample_Name, "(?<=-)[WD]$"),
    site = str_remove(Sample_Name, "-[WD]$")
  ) %>%
  filter(!is.na(site), !is.na(treatment), !is.na(rate_mgkg_h))

paired <- resp %>%
  select(site, treatment, rate_mgkg_h) %>%
  pivot_wider(
    names_from = treatment,
    values_from = rate_mgkg_h
  ) %>%
  filter(!is.na(W), !is.na(D)) %>%
  mutate(
    abs_diff = abs(W - D),
    wet_rate = abs(W),
    slower_treatment = case_when(
      abs(D) < abs(W) ~ "Dry",
      abs(W) < abs(D) ~ "Wet",
      TRUE ~ "Tie"
    ),
    rel_diff = abs_diff / wet_rate
  ) %>%
  filter(slower_treatment != "Tie") %>%
  left_join(site_meta, by = "site")

paired2 = paired %>%
  filter(D < 0)

## Plot 1 -----------------------------------------------------------------
ggplot(
  paired,
  aes(
    x = abs_diff,
    y = rel_diff,
    color = slower_treatment,
    shape = flow_regime
  )
) +
  geom_point(size = 3, alpha = 0.7, stroke = 1) +
 # geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.8) +
  labs(
    x = "Absolute difference between wet and dry median respiration rates (mg/kg/h)",
    y = "Absolute difference / wet absolute rate",
    color = "Slower rate from",
    shape = "Flow regime"
  ) +
  scale_color_manual(values = c("Wet" = "red", "Dry" = "dodgerblue2")) +
  scale_shape_manual(values = c("Perennial" = 16, "Non-perennial" = 17)) +
  theme_bw(base_size = 12)

#

max_of_wet <- max(abs(paired$W), na.rm = TRUE)

ggplot(
  paired,
  aes(
    x = abs_diff,
    y = rel_diff,
    color = slower_treatment,
    shape = flow_regime
  )
) +
  geom_point(size = 3, alpha = 0.7, stroke = 1) +
  #geom_point(x = 1000, y = (1000/3000), size = 3)+
  geom_abline(intercept = 0, slope = 1 / abs(max_of_wet), linetype = "dashed") +
  labs(
    x = "Absolute difference between wet and dry median respiration rates (mg/kg/h)",
    y = "Absolute difference / wet absolute rate",
    color = "Slower rate from",
    shape = "Flow regime"
  ) +
  scale_color_manual(values = c("Wet" = "red", "Dry" = "dodgerblue2")) +
  scale_shape_manual(values = c("Perennial" = 16, "Non-perennial" = 17)) +
  theme_bw(base_size = 12)
