# This script makes histogram figures for ECA physical manuscript
library(tidyverse); library(ggpubr)
rm(list=ls()); graphics.off()

# Functions ---------------------------------------------------------------
# Transformation for normalization is cube root - have to cube root then add sign back to value to make it positive or negative
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

## Read in/Merge Data for Histograms -------------------------------------
# Individual Rate data for histograms ####
all_data = read.csv("Data/EC_Sediment_SpC_pH_Temp_Respiration.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%  
  select(-c(Field_Name, IGSN, Material))
all_data = all_data %>%
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = (as.numeric(Respiration_Rate_mg_DO_per_kg_per_H))*-1)

# Reading in metadata from both files to get river type
metadata = read.csv("Data/EC_Field_Metadata.csv") %>%
  select(c(Parent_ID, Intermittent_or_Perennial))

all_data_with_metadata = all_data %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^[^_]+_[^_]+")) %>%  # Extract "EC_005" from "EC_005_INC-D1"
  left_join(metadata, by = "Parent_ID")


## Effect Size Data for histograms ####
effect = read.csv("Data/EC_Sediment_Effect_Size.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%
  select(c(Sample_Name, Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  mutate(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = (as.numeric(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H))*-1)

effect_with_metadata = effect %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^[^_]+_[^_]+")) %>%  # Extract "EC_005" from "EC_005_all"
  left_join(metadata, by = "Parent_ID")

all_data_filtered = all_data_with_metadata

# Get O2 consumption rate values for intermittent sites (from filtered data)
intermittent_o2_values = all_data_filtered %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != 9999) %>%
  filter(Intermittent_or_Perennial == "Intermittent", !is.na(Intermittent_or_Perennial)) %>%
  pull(Respiration_Rate_mg_DO_per_kg_per_H) %>%
  as.numeric() 

print("Intermittent O2 consumption rate values (filtered):")
print(intermittent_o2_values)

# Take Cube root of all respiration rates for figure (from filtered data)
cube_respiration = all_data_filtered %>%
  select(c(Sample_Name, Respiration_Rate_mg_DO_per_kg_per_H, Parent_ID, Intermittent_or_Perennial)) %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != 9999) %>%
  mutate(cube_Respiration_mg_kg = cube_root((as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)))) %>%
  mutate(Treat = if_else(grepl("D", Sample_Name), "Dry", "Wet"))

# Get cube root O2 consumption rate values for intermittent sites BY TREATMENT
intermittent_cube_o2_by_treatment = cube_respiration %>%
  filter(Intermittent_or_Perennial == "Intermittent", !is.na(Intermittent_or_Perennial)) %>%
  select(Parent_ID, cube_Respiration_mg_kg, Treat)

print("Intermittent cube root O2 consumption rate values by treatment:")
print(intermittent_cube_o2_by_treatment)


####################################################################
# Two-panel figure: (A) untransformed rates, (B) cube-root transformed rates
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# --- Build site-level paired (median) treatment rates ---
site_treat_rates <- all_data_filtered %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != 9999) %>%
  mutate(
    Treat = if_else(grepl("D", Sample_Name), "Dry", "Wet"),
    Rate  = as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)
  ) %>%
  group_by(Parent_ID, Intermittent_or_Perennial, Treat) %>%
  summarise(Treat_Rate = median(Rate, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Treat, values_from = Treat_Rate) %>%
  filter(!is.na(Dry) & !is.na(Wet)) %>%
  mutate(
    Dry_cube = cube_root(Dry),
    Wet_cube = cube_root(Wet)
  )

# --- One effect size per site (for color) ---
site_effects <- effect_with_metadata %>%
  mutate(
    Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H =
      as.numeric(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)
  ) %>%
  group_by(Parent_ID, Intermittent_or_Perennial) %>%
  summarise(
    Effect_Size = median(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Cube_Effect_Size = cube_root(Effect_Size),
    Abs_Cube_Effect_Size = abs(Cube_Effect_Size)   # == cube_root(abs(Effect_Size))
  )

scatter_df <- site_treat_rates %>%
  left_join(site_effects, by = c("Parent_ID", "Intermittent_or_Perennial")) %>%
  filter(!is.na(Effect_Size))


# --- Color scale to MATCH Figure 1 histogram palette ---
# Histogram used: low="dodgerblue2", mid="goldenrod2", high="firebrick2"
col_scale_abs <- scale_color_gradient(
  name = "cube-root\nEffect Size",
  low = "goldenrod2",      # 0
  high = "dodgerblue2",    # 10
  limits = c(0, 10),
  oob = squish
)

legend_topleft <- theme(
  legend.position = c(0.15, 0.85),
  legend.justification = c(0, 1),
  legend.background = element_rect(fill = "white", color = NA),
  legend.key = element_rect(fill = "white", color = NA)
)

base_theme <- theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9)
  )

# --- Panel A: untransformed axes ---
pA <- ggplot(scatter_df, aes(x = Dry, y = Wet, color = Abs_Cube_Effect_Size)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_point(size = 3, alpha = 0.9) +
  col_scale_abs +
  base_theme +
  legend_topleft +
  xlab(expression(atop("Dry treatment median O"[2]*" consumption rate",
                       "(mg O"[2]*" kg"^-1*" h"^-1*")"))) +
  ylab(expression(atop("Wet treatment median O"[2]*" consumption rate",
                       "(mg O"[2]*" kg"^-1*" h"^-1*")")))

# --- Panel B: cube-root transformed axes ---
pB <- ggplot(scatter_df, aes(x = Dry_cube, y = Wet_cube, color = Abs_Cube_Effect_Size)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_point(size = 3, alpha = 0.9) +
  col_scale_abs +
  base_theme +
  theme(legend.position = "none") +
  xlab(expression(atop("Dry treatment median O"[2]*" consumption rate"^(1/3),
                       "(mg O"[2]*" kg"^-1*" h"^-1*")"))) +
  ylab(expression(atop("Wet treatment median O"[2]*" consumption rate"^(1/3),
                       "(mg O"[2]*" kg"^-1*" h"^-1*")")))

# --- Combine into 2-panel figure labeled A and B ---
two_panel_scatter <- ggarrange(
  pA, pB,
  ncol = 2,
  labels = c("A", "B"),
  hjust = -1,
  vjust = 1.5,
  common.legend = FALSE
)

ggsave("Figures/SI_TwoPanel_Scatter_Untransformed_vs_CubeRates_EffectColor_Fig1Colors_0yellow_10blue.png",
       plot = two_panel_scatter, width = 12, height = 5.5, dpi = 300)
ggsave("Figures/SI_TwoPanel_Scatter_Untransformed_vs_CubeRates_EffectColor_Fig1Colors_0yellow_10blue.pdf",
       plot = two_panel_scatter, width = 12, height = 5.5)
