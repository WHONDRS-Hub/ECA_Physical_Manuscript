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
  filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999) %>% # removes overexposed samples, missing replicates
  select(-c(Field_Name, IGSN, Material))

ev_data = read.csv("Data/WHONDRS_EV_Data_Package/Sample_Data/WHONDRS_EV_Sediment_SpC_pH_Temp_Respiration.csv", skip = 2) %>%
  filter(grepl("EV", Sample_Name)) %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999) %>% # removes overexposed samples, missing replicates
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

# Diagnostics only: these do NOT filter the data
sites_with_negative_effect = effect_with_metadata %>%
  filter(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H < 0)

print(paste("Number of sites with negative effect size:", nrow(sites_with_negative_effect)))
print("Sites with negative effect size:")
print(sites_with_negative_effect)

all_data_filtered = all_data_with_metadata

abs(median(all_data_filtered$Respiration_Rate_mg_DO_per_kg_per_H[grep('-D', all_data_filtered$Sample_Name )], na.rm = TRUE))
abs(median(all_data_filtered$Respiration_Rate_mg_DO_per_kg_per_H[grep('-W', all_data_filtered$Sample_Name )], na.rm = TRUE))

print(paste("Number of samples after filtering for data quality:", nrow(all_data_filtered)))
print(paste("Number of samples BEFORE filtering:", nrow(all_data_with_metadata)))

# Get O2 consumption rate values for intermittent sites
intermittent_o2_values = all_data_filtered %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != 9999) %>%
  filter(Intermittent_or_Perennial == "Intermittent", !is.na(Intermittent_or_Perennial)) %>%
  pull(Respiration_Rate_mg_DO_per_kg_per_H) %>%
  as.numeric() 

print("Intermittent O2 consumption rate values:")
print(intermittent_o2_values)

# Take Cube root of all respiration rates for figure
cube_respiration = all_data_filtered %>%
  select(c(Sample_Name, Respiration_Rate_mg_DO_per_kg_per_H, Parent_ID, Intermittent_or_Perennial)) %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != 9999) %>%
  mutate(cube_Respiration_mg_kg = cube_root((as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)))) %>%
  mutate(Treat = if_else(grepl("-D", Sample_Name), "Dry", "Wet"))

# Get cube root O2 consumption rate values for intermittent sites BY TREATMENT
intermittent_cube_o2_by_treatment = cube_respiration %>%
  filter(Intermittent_or_Perennial == "Intermittent", !is.na(Intermittent_or_Perennial)) %>%
  select(Parent_ID, cube_Respiration_mg_kg, Treat)

print("Intermittent cube root O2 consumption rate values by treatment:")
print(intermittent_cube_o2_by_treatment)

# Create data frame for vertical lines with treatment information
vlines_cube_o2_by_treatment = intermittent_cube_o2_by_treatment %>%
  rename(xintercept = cube_Respiration_mg_kg) %>%
  mutate(line_type = paste("Intermittent:", Treat))

# Get effect size values for intermittent sites
# IMPORTANT: no filtering by positive/negative effect size here
intermittent_effect_values = effect_with_metadata %>%
  filter(Intermittent_or_Perennial == "Intermittent", !is.na(Intermittent_or_Perennial)) %>%
  pull(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)

print("Intermittent effect size values:")
print(intermittent_effect_values)

# Cube root transform effect size data 
# IMPORTANT: no filtering by positive/negative effect size here
cube_effect = effect_with_metadata %>%
  mutate(cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = cube_root(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H))

# Get cube root effect size values for intermittent sites
intermittent_cube_effect_values = cube_effect %>%
  filter(Intermittent_or_Perennial == "Intermittent", !is.na(Intermittent_or_Perennial)) %>%
  pull(cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)

print("Intermittent cube root effect size values:")
print(intermittent_cube_effect_values)

# Create data frames for vertical lines
# For O2, include Treat so dashed lines can be colored by treatment
vlines_o2 = all_data_filtered %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != 9999) %>%
  filter(Intermittent_or_Perennial == "Intermittent", !is.na(Intermittent_or_Perennial)) %>%
  mutate(Treat = if_else(grepl("-D", Sample_Name), "Dry", "Wet")) %>%
  select(
    xintercept = Respiration_Rate_mg_DO_per_kg_per_H,
    Treat
  )

vlines_effect = data.frame(
  xintercept = intermittent_effect_values,
  line_type = "Intermittent sites"
)

vlines_cube_effect = data.frame(
  xintercept = intermittent_cube_effect_values,
  line_type = "Intermittent sites"
)

# Histograms --------------------------------------------------------------
## Histogram of all O2 consumption rates WITH VERTICAL LINES FOR INTERMITTENT SITES
all_hist = all_data_filtered %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != 9999) %>%
  mutate(Treat = if_else(grepl("-D", Sample_Name), "Dry", "Wet")) %>%
  ggplot(aes(x = (as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)))) +
  geom_histogram(
    position = "identity",
    alpha = 0.8,
    aes(fill = Treat)
  ) +
  
  # Add vertical lines for each intermittent site's O2 consumption rate
  # Dashed lines are colored by the treatment that each value comes from
  {if(nrow(vlines_o2) > 0)
    geom_vline(
      data = vlines_o2,
      aes(xintercept = xintercept, color = Treat),
      linetype = "dashed",
      linewidth = 1,
      show.legend = FALSE
    )
  } +
  
  scale_fill_manual(
    name = "Treatment",
    values = c("Dry" = "#D55E00", "Wet" = "#0072B2")
  ) +
  scale_color_manual(
    values = c("Dry" = "#D55E00", "Wet" = "#0072B2"),
    guide = "none"
  ) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.85),
        legend.key.size = unit(0.15, "in"),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.box = "vertical") +
  guides(fill = guide_legend(title = "Treatment", order = 1)) +
  xlab(expression(atop("\n O"[2]*" Consumption Rate", "(mg O"[2]*" kg"^-1*" h"^-1*")"))) +
  ylab("Count")

## Effect size histogram WITH VERTICAL LINES FOR INTERMITTENT SITES
effect_limits <- c(-1500, 1500)

effect_hist <- ggplot(cube_effect, aes(x = Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_histogram(
    binwidth = 25,
    boundary = 0,
    aes(fill = after_stat(ifelse(x < 0, NA, x)))
  ) +
  
  # Add vertical lines for each intermittent site's effect size value
  {if(length(intermittent_effect_values) > 0)
    geom_vline(
      data = vlines_effect,
      aes(xintercept = xintercept),
      color = "red",
      linetype = "dashed",
      linewidth = 1,
      show.legend = FALSE
    )
  } +
  
  scale_fill_gradient2(
    name = "Effect Size",
    limits = c(0, 1500),
    low = "dodgerblue2",
    mid = "goldenrod2",
    high = "firebrick2",
    midpoint = 750,
    na.value = "grey70",
    guide = "none"
  ) +
  
  scale_x_continuous(
    breaks = c(-1500, -1000, -500, 0, 500, 1000, 1500)
  ) +
  
  coord_cartesian(xlim = effect_limits) +
  
  theme_bw() +
  ylab("Count\n") +
  theme(
    legend.position = c(0.85, 0.85),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.key.width = unit(0.4, "in")
  ) +
  xlab(expression(atop(
    "\n Effect Size O"[2]*" Consumption Rate",
    "(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" h"^-1*")"
  )))

## Combine SI histograms
combined_SI_hist = ggarrange(all_hist, effect_hist, ncol = 2, labels = c("A", "B"), hjust = -1, vjust = 2.5)
ggsave("Figures/Figure_S1_Combined_SI_Histograms.png", plot = combined_SI_hist, width = 10, height = 5, dpi = 300)
ggsave("Figures/Figure_S1_Combined_SI_Histograms.pdf", plot = combined_SI_hist, width = 10, height = 5, dpi = 300)

all_cube_hist = ggplot(cube_respiration, aes(x = cube_Respiration_mg_kg)) +
  geom_histogram(position = "identity", alpha = 0.8, aes(fill = Treat)) +
  # Add vertical lines for each intermittent site's cube root O2 consumption rate BY TREATMENT
  {if(nrow(vlines_cube_o2_by_treatment) > 0)
    geom_vline(data = vlines_cube_o2_by_treatment,
               aes(xintercept = xintercept, color = Treat),
               linetype = "dashed",
               linewidth = 1,
               show.legend = FALSE)
  } +
  scale_fill_manual(
    name = "Treatment",
    values = c("Dry" = "#D55E00", "Wet" = "#0072B2")
  ) +
  scale_color_manual(values = c("Dry" = "#D55E00", "Wet" = "#0072B2"), guide = "none") +
  coord_cartesian(ylim = c(0, 87.5)) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.85),
        legend.key.size = unit(0.15, "in"),
        legend.key.width = unit(0.4, "in"),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.box = "vertical") +
  guides(fill = guide_legend(title = "Treatment", order = 1)) +
  xlab(expression(atop("\n O"[2]*" Consumption Rate" ^(1/3)*"", "(mg O"[2]*" kg"^-1*" h"^-1*")"))) +
  ylab("Count")

# Cube root effect size histogram WITH VERTICAL LINES FOR INTERMITTENT SITES
cube_effect_limits <- c(-12, 12)

cube_effect_hist = ggplot(cube_effect, aes(x = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_histogram(
    binwidth = 0.5,
    boundary = 0,
    aes(fill = after_stat(ifelse(x < 0, NA, x)))
  ) +
  # Add vertical lines for each intermittent site's cube root effect size value
  {if(length(intermittent_cube_effect_values) > 0)
    geom_vline(data = vlines_cube_effect,
               aes(xintercept = xintercept),
               color = "red",
               linetype = "dashed",
               linewidth = 1,
               show.legend = FALSE)
  } +
  scale_fill_gradient2(name = "Cubed Root\nEffect Size",
                       limits = c(0, 12),
                       low = "dodgerblue2",
                       mid = "goldenrod2",
                       high = "firebrick2",
                       midpoint = 6,
                       na.value = "grey70",
                       guide = "none") +
  scale_x_continuous(
    breaks = c(-12, -8, -4, 0, 4, 8, 12)
  ) +
  theme_bw() +
  coord_cartesian(xlim = cube_effect_limits) +
  ylab("Count\n") +
  theme(legend.position = c(0.85, 0.85),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.key.width = unit(0.4, "in")) +
  xlab(expression(atop("\n Effect Size O"[2]*" Consumption Rate" ^(1/3),"(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" h"^-1*")")))

## Combine cube root histograms
combined_hist = ggarrange(all_cube_hist, cube_effect_hist, ncol = 2, labels = c("A", "B"), hjust = -5, vjust = 2.5)
ggsave("Figures/Figure_2_Combined_Cube_Histograms.png", plot = combined_hist, width = 10, height = 5, dpi = 300)
ggsave("Figures/Figure_2_Combined_Cube_Histograms.pdf", plot = combined_hist, width = 10, height = 5, dpi = 300)
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
    Treat = case_when(
      grepl("-D", Sample_Name) ~ "Dry",
      grepl("-W", Sample_Name) ~ "Wet",
      TRUE ~ NA_character_
    ),
    Rate  = as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)
  ) %>%
  filter(!is.na(Treat)) %>%
  group_by(Parent_ID, Intermittent_or_Perennial, Treat) %>%
  summarise(Treat_Rate = median(Rate, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Treat, values_from = Treat_Rate)

# --- One effect size per site (for color) ---
# IMPORTANT: this keeps the sign of the effect size
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
    
    # Match histogram color logic:
    # negative effect sizes become NA and are colored grey using na.value
    # positive effect sizes use the blue-gold-red gradient
    Cube_Effect_Size_for_Color = if_else(Cube_Effect_Size < 0, NA_real_, Cube_Effect_Size)
  )

scatter_df <- site_treat_rates %>%
  left_join(site_effects, by = c("Parent_ID", "Intermittent_or_Perennial"))

# --- Color scale to MATCH histogram palette and show negative values ---
# Histogram logic:
# negative effect sizes = grey
# positive effect sizes = dodgerblue2 to goldenrod2 to firebrick2
col_scale_effect <- scale_color_gradientn(
  name = "cube-root\nEffect Size",
  colors = c("grey70", "grey70", "dodgerblue2", "goldenrod2", "firebrick2"),
  values = scales::rescale(c(-12, -0.001, 0, 6, 12), from = c(-12, 12)),
  limits = c(-12, 12),
  breaks = c(-12, -8, -4, 0, 4, 8, 12),
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
pA <- ggplot(scatter_df, aes(x = Dry, y = Wet, color = Cube_Effect_Size)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_point(size = 3, alpha = 0.9) +
  col_scale_effect +
  base_theme +
  legend_topleft +
  xlab(expression(atop("Dry treatment median O"[2]*" consumption rate",
                       "(mg O"[2]*" kg"^-1*" h"^-1*")"))) +
  ylab(expression(atop("Wet treatment median O"[2]*" consumption rate",
                       "(mg O"[2]*" kg"^-1*" h"^-1*")")))

# --- Panel B: cube-root transformed axes ---
pB <- ggplot(scatter_df, aes(x = Dry_cube, y = Wet_cube, color = Cube_Effect_Size)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_point(size = 3, alpha = 0.9) +
  col_scale_effect +
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

ggsave("Figures/Figure_S2_TwoPanel_Scatter_Untransformed_vs_CubeRates_EffectColor_Fig1Colors_0yellow_10blue.png",
       plot = two_panel_scatter, width = 12, height = 5.5, dpi = 300)
ggsave("Figures/Figure_S2_TwoPanel_Scatter_Untransformed_vs_CubeRates_EffectColor_Fig1Colors_0yellow_10blue.pdf",
       plot = two_panel_scatter, width = 12, height = 5.5)
