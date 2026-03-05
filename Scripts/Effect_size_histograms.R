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
ev_data = read.csv("Data/WHONDRS_EV_Data_Package/Sample_Data/WHONDRS_EV_Sediment_SpC_pH_Temp_Respiration.csv", skip = 2) %>%
  filter(grepl("EV", Sample_Name)) %>%
  select(-c(Field_Name, IGSN, Material))
all_data = rbind(all_data, ev_data) 
all_data = all_data %>%
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = (as.numeric(Respiration_Rate_mg_DO_per_kg_per_H))*-1)

# Reading in metadata from both files to get river type
metadata = read.csv("Data/EC_Field_Metadata.csv") %>%
  select(c(Parent_ID, Intermittent_or_Perennial))
metadata2 = read.csv("Data/WHONDRS_EV_Data_Package/WHONDRS_EV_Field_Metadata.csv") %>%
  select(c(Parent_ID, Intermittent_or_Perennial))
metadata = rbind(metadata, metadata2)

all_data_with_metadata = all_data %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^[^_]+_[^_]+")) %>%  # Extract "EC_005" from "EC_005_INC-D1"
  left_join(metadata, by = "Parent_ID")

# Check merge
print("All data with metadata:")
print(head(all_data_with_metadata %>% select(Sample_Name, Parent_ID, Intermittent_or_Perennial)))
print(paste("Number of intermittent samples in all_data:", sum(all_data_with_metadata$Intermittent_or_Perennial == "Intermittent", na.rm = TRUE)))

## Effect Size Data for histograms ####
effect = read.csv("Data/EC_Sediment_Effect_Size.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%
  select(c(Sample_Name, Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  mutate(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = (as.numeric(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H))*-1)

effect_ev = read.csv("Data/WHONDRS_EV_Data_Package/Sample_Data/WHONDRS_EV_Sediment_Effect_Size.csv", skip = 2) %>%
  filter(grepl("EV", Sample_Name)) %>%
  select(c(Sample_Name, Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = Effect_Size_Median_Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  mutate(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = (as.numeric(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H))*-1)

effect = rbind(effect, effect_ev)

effect_with_metadata = effect %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^[^_]+_[^_]+")) %>%  # Extract "EC_005" from "EC_005_all"
  left_join(metadata, by = "Parent_ID")

# Check merge
print("Effect size data with metadata:")
print(head(effect_with_metadata %>% select(Sample_Name, Parent_ID, Intermittent_or_Perennial)))
print(paste("Number of intermittent sites in effect data:", sum(effect_with_metadata$Intermittent_or_Perennial == "Intermittent", na.rm = TRUE)))

#  Filter to only keep sites with positive effect size (excluding negative effect sizes due to uncertainty in the calculation of the rates)
sites_with_negative_effect = effect_with_metadata %>%
  filter(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H >= 0) %>%
  pull(Parent_ID)

print(paste("Number of sites with negative effect size:", length(sites_with_negative_effect)))
print("Sites with negative effect size:")
print(sites_with_negative_effect)

all_data_filtered = all_data_with_metadata %>%
  filter(Parent_ID %in% sites_with_negative_effect)

print(paste("Number of samples after filtering for negative effect size sites:", nrow(all_data_filtered)))
print(paste("Number of samples BEFORE filtering:", nrow(all_data_with_metadata)))

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

# Create data frame for vertical lines with treatment information**
vlines_cube_o2_by_treatment = intermittent_cube_o2_by_treatment %>%
  rename(xintercept = cube_Respiration_mg_kg) %>%
  mutate(line_type = paste("Intermittent:", Treat))

# Get effect size values for intermittent sites (already filtered to negative)
intermittent_effect_values = effect_with_metadata %>%
  filter(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H >= 0) %>%  # Changed to => to match filtering above
  filter(Intermittent_or_Perennial == "Intermittent", !is.na(Intermittent_or_Perennial)) %>%
  pull(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)

print("Intermittent effect size values:")
print(intermittent_effect_values)

# Cube root transform effect size data 
cube_effect = effect_with_metadata %>%
  filter(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H >= 0) %>%  # Changed the sign here
  mutate(cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = cube_root(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H))

# Get cube root effect size values for intermittent sites
intermittent_cube_effect_values = cube_effect %>%
  filter(Intermittent_or_Perennial == "Intermittent", !is.na(Intermittent_or_Perennial)) %>%
  pull(cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)

print("Intermittent cube root effect size values:")
print(intermittent_cube_effect_values)

# Create data frames for vertical lines (for legend purposes)
vlines_o2 = data.frame(
  xintercept = intermittent_o2_values,
  line_type = "Intermittent sites"
)

vlines_effect = data.frame(
  xintercept = (intermittent_effect_values),
  line_type = "Intermittent sites"
)

vlines_cube_effect = data.frame(
  xintercept = intermittent_cube_effect_values,
  line_type = "Intermittent sites"
)

# Histograms --------------------------------------------------------------
## Histogram of all O2 consumption rates WITH VERTICAL LINES FOR INTERMITTENT SITES
all_hist = all_data_filtered %>%  # **CHANGED: Use filtered data**
  filter(Respiration_Rate_mg_DO_per_kg_per_H != 9999) %>%
  mutate(Treat = if_else(grepl("D", Sample_Name), "Dry", "Wet")) %>%
  ggplot(aes(x = (as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)))) +
  geom_histogram(position = "identity", alpha = 0.8, aes(fill = Treat)) +
  scale_fill_manual(values = c("#D55E00","#0072B2")) +
  scale_linetype_manual(name = "", values = c("Intermittent sites" = "dashed")) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.85),
        legend.key.size = unit(0.15, "in"),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 14),  # **CHANGED: Increased from 10 to 14**
        axis.title.y = element_text(size = 14),  # **NEW: Added y-axis title size**
        axis.text.x = element_text(size = 12),   # **NEW: Increased axis text size**
        axis.text.y = element_text(size = 12),   # **NEW: Increased axis text size**
        legend.box = "vertical") +
  guides(fill = guide_legend(title = "Treatment", order = 1),
         linetype = guide_legend(title = "", order = 2)) +
  xlab(expression(atop("\n O"[2]*" Consumption Rate", "(mg O"[2]*" kg"^-1*" h"^-1*")"))) +
  ylab("Count")

## Effect size histogram WITH VERTICAL LINES FOR INTERMITTENT SITES
effect_limits <- c(0, 1400)  # Changed to start at 0 
effect_hist = ggplot(cube_effect, aes(x = (Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H))) +
  geom_histogram(binwidth = 25, aes(fill = after_stat(x))) +
  # Add vertical lines for each intermittent site's effect size value
  {if(length(intermittent_effect_values) > 0)
    geom_vline(data = vlines_effect, aes(xintercept = xintercept, linetype = line_type),
               color = "red", linewidth = 1)} +
  scale_fill_gradient2(name = "Effect Size", limits = effect_limits, low = "dodgerblue2", mid = "goldenrod2",
                       high = "firebrick2", midpoint = (max(effect_limits)+min(effect_limits))/2) +
  scale_linetype_manual(name = "", values = c("Intermittent sites" = "dashed")) +
  theme_bw() +
  xlim(c(0, 1500)) +  # Changed to start at 0
  ylab("Count\n") +
  theme(legend.position = c(0.85, 0.85),
        axis.title.x = element_text(size = 14),  # **CHANGED: Increased from 10 to 14**
        axis.title.y = element_text(size = 14),  # **NEW: Added y-axis title size**
        axis.text.x = element_text(size = 12),   # **NEW: Increased axis text size**
        axis.text.y = element_text(size = 12),   # **NEW: Increased axis text size**
        legend.key.width = unit(0.4, "in")) +  # **NEW: Wider legend key for clearer dashed line**
  guides(linetype = guide_legend(title = "", override.aes = list(linewidth = 1))) +  # **NEW: Make dashed line thicker in legend**
  xlab(expression(atop("\n |Effect Size O"[2]*" Consumption Rate|","(|Median Wet - Median Dry Rate|; mg O"[2]*" kg"^-1*" h"^-1*")")))

## Combine SI histograms
combined_SI_hist = ggarrange(all_hist, effect_hist, ncol = 2, labels = c("A", "B"), hjust = -1, vjust = 2.5)
ggsave("Figures/Combined_SI_Histograms.png", plot = combined_SI_hist, width = 10, height = 5, dpi = 300)

all_cube_hist = ggplot(cube_respiration, aes(x = cube_Respiration_mg_kg)) +
  geom_histogram(position = "identity", alpha = 0.8, aes(fill = Treat)) +
  # Add vertical lines for each intermittent site's cube root O2 consumption rate BY TREATMENT
  {if(nrow(vlines_cube_o2_by_treatment) > 0)
    geom_vline(data = vlines_cube_o2_by_treatment,
               aes(xintercept = xintercept, color = Treat, linetype = "Intermittent sites"),
               linewidth = 1)} +
  scale_fill_manual(values = c("#D55E00","#0072B2")) +
  scale_color_manual(values = c("Dry" = "#D55E00", "Wet" = "#0072B2"), guide = "none") +  # Match fill colors
  scale_linetype_manual(name = "", values = c("Intermittent sites" = "dashed")) +
  ylim(0, 87.5) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.85),
        legend.key.size = unit(0.15, "in"),
        legend.key.width = unit(0.4, "in"),  # **NEW: Wider legend key for clearer dashed line**
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 14),  # **CHANGED: Increased from 10 to 14**
        axis.title.y = element_text(size = 14),  # **NEW: Added y-axis title size**
        axis.text.x = element_text(size = 12),   # **NEW: Increased axis text size**
        axis.text.y = element_text(size = 12),   # **NEW: Increased axis text size**
        legend.box = "vertical") +
  guides(fill = guide_legend(title = "Treatment", order = 1),
         linetype = guide_legend(title = "", order = 2, override.aes = list(linewidth = 1))) +  # **NEW: Make dashed line thicker in legend**
  xlab(expression(atop("\n O"[2]*" Consumption Rate" ^(1/3)*"", "(mg O"[2]*" kg"^-1*" h"^-1*")"))) +
  ylab("Count")

# Cube root effect size histogram WITH VERTICAL LINES FOR INTERMITTENT SITES
cube_effect_limits <- c(0, 10)  
cube_effect_hist = ggplot(cube_effect, aes(x = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) +  # **Plot absolute values**
  geom_histogram(binwidth = 0.5, aes(fill = after_stat(x))) +
  # Add vertical lines for each intermittent site's cube root effect size value
  {if(length(intermittent_cube_effect_values) > 0)
    geom_vline(data = vlines_cube_effect,
               aes(xintercept = xintercept, linetype = line_type),
               color = "red", linewidth = 1)} +
  scale_fill_gradient2(name = "Cubed Root\nEffect Size", limits = cube_effect_limits, low = "dodgerblue2", mid = "goldenrod2",
                       high = "firebrick2", midpoint = 5) +  # **Midpoint at 5 (middle of 0-10)**
  scale_linetype_manual(name = "", values = c("Intermittent sites" = "dashed")) +
  theme_bw() +
  xlim(c(-10, 10)) +  # **CHANGED: X-axis from -10 to 10 (but data only on 0-10 side)**
  ylab("Count\n") +
  theme(legend.position = c(0.85, 0.85),
        axis.title.x = element_text(size = 14),  # **CHANGED: Increased from 10 to 14**
        axis.title.y = element_text(size = 14),  # **NEW: Added y-axis title size**
        axis.text.x = element_text(size = 12),   # **NEW: Increased axis text size**
        axis.text.y = element_text(size = 12),   # **NEW: Increased axis text size**
        legend.key.width = unit(0.4, "in")) +  # **NEW: Wider legend key for clearer dashed line**
  guides(linetype = guide_legend(title = "", override.aes = list(linewidth = 1))) +  # **NEW: Make dashed line thicker in legend**
  xlab(expression(atop("\n Effect Size O"[2]*" Consumption Rate" ^(1/3),"(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" h"^-1*")")))

## Combine cube root histograms
combined_hist = ggarrange(all_cube_hist, cube_effect_hist, ncol = 2, labels = c("A", "B"), hjust = -5, vjust = 2.5)
ggsave("Figures/Combined_Cube_Histograms.png", plot = combined_hist, width = 10, height = 5, dpi = 300)
ggsave("Figures/Combined_Cube_Histograms.pdf", plot = combined_hist, width = 10, height = 5, dpi = 300)


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
