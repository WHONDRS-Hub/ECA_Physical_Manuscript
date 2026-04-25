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
  mutate(Parent_ID = str_extract(Sample_Name, "^[^_]+_[^_]+")) %>%  
  left_join(metadata, by = "Parent_ID")

# Check merge
print("All data with metadata:")
print(head(all_data_with_metadata %>% select(Sample_Name, Parent_ID, Intermittent_or_Perennial)))
print(paste("Number of intermittent samples in all_data:", sum(all_data_with_metadata$Intermittent_or_Perennial == "Intermittent", na.rm = TRUE)))
print(paste("Total EC samples:", sum(grepl("^EC_", all_data_with_metadata$Sample_Name))))
print(paste("Total EV samples:", sum(grepl("^EV_", all_data_with_metadata$Sample_Name))))

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
  mutate(Parent_ID = str_extract(Sample_Name, "^[^_]+_[^_]+")) %>%  
  left_join(metadata, by = "Parent_ID")

# Check merge
print("Effect size data with metadata:")
print(head(effect_with_metadata %>% select(Sample_Name, Parent_ID, Intermittent_or_Perennial)))
print(paste("Number of intermittent sites in effect data:", sum(effect_with_metadata$Intermittent_or_Perennial == "Intermittent", na.rm = TRUE)))
print(paste("Total EC effect sites:", sum(grepl("^EC_", effect_with_metadata$Sample_Name))))
print(paste("Total EV effect sites:", sum(grepl("^EV_", effect_with_metadata$Sample_Name))))

# NO FILTERING - keep all sites
all_data_filtered = all_data_with_metadata

print(paste("Number of samples (all included):", nrow(all_data_filtered)))

# Get O2 consumption rate values for NON-PERENNIAL sites (changed from intermittent)
non_perennial_o2_values = all_data_filtered %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != 9999) %>%
  filter(Intermittent_or_Perennial != "Perennial", !is.na(Intermittent_or_Perennial)) %>%
  pull(Respiration_Rate_mg_DO_per_kg_per_H) %>%
  as.numeric()

print("Non-perennial O2 consumption rate values:")
print(non_perennial_o2_values)

# Take Cube root of all respiration rates for figure
cube_respiration = all_data_filtered %>%
  select(c(Sample_Name, Respiration_Rate_mg_DO_per_kg_per_H, Parent_ID, Intermittent_or_Perennial)) %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != 9999) %>%
  mutate(cube_Respiration_mg_kg = cube_root((as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)))) %>%
  mutate(Treat = if_else(grepl("D", Sample_Name), "Dry", "Wet"))

# Get cube root O2 consumption rate values for NON-PERENNIAL sites BY TREATMENT
non_perennial_cube_o2_by_treatment = cube_respiration %>%
  filter(Intermittent_or_Perennial != "Perennial", !is.na(Intermittent_or_Perennial)) %>%
  select(Parent_ID, cube_Respiration_mg_kg, Treat)

print("Non-perennial cube root O2 consumption rate values by treatment:")
print(non_perennial_cube_o2_by_treatment)

# Create data frame for vertical lines with treatment information
vlines_cube_o2_by_treatment = non_perennial_cube_o2_by_treatment %>%
  rename(xintercept = cube_Respiration_mg_kg) %>%
  mutate(line_type = paste("Non-perennial:", Treat))

# Get effect size values for NON-PERENNIAL sites
non_perennial_effect_values = effect_with_metadata %>%
  filter(Intermittent_or_Perennial != "Perennial", !is.na(Intermittent_or_Perennial)) %>%
  pull(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)

print("Non-perennial effect size values:")
print(non_perennial_effect_values)

# Cube root transform effect size data (NO FILTERING)
cube_effect = effect_with_metadata %>%
  mutate(cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = cube_root(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H))

# Get cube root effect size values for NON-PERENNIAL sites
non_perennial_cube_effect_values = cube_effect %>%
  filter(Intermittent_or_Perennial != "Perennial", !is.na(Intermittent_or_Perennial)) %>%
  pull(cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)

print("Non-perennial cube root effect size values:")
print(non_perennial_cube_effect_values)

# Create data frames for vertical lines (for legend purposes)
vlines_o2 = data.frame(
  xintercept = non_perennial_o2_values,
  line_type = "Non-perennial sites"
)

vlines_effect = data.frame(
  xintercept = non_perennial_effect_values,
  line_type = "Non-perennial sites"
)

vlines_cube_effect = data.frame(
  xintercept = non_perennial_cube_effect_values,
  line_type = "Non-perennial sites"
)

# Histograms --------------------------------------------------------------
## Histogram of all O2 consumption rates WITH VERTICAL LINES FOR NON-PERENNIAL SITES
all_hist = all_data_filtered %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != 9999) %>%
  mutate(Treat = if_else(grepl("D", Sample_Name), "Dry", "Wet")) %>%
  ggplot(aes(x = (as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)))) +
  geom_histogram(position = "identity", alpha = 0.8, aes(fill = Treat)) +
  {if(length(non_perennial_o2_values) > 0)
    geom_vline(data = vlines_o2, aes(xintercept = xintercept, linetype = line_type),
               color = "red", linewidth = 1)} +
  scale_fill_manual(values = c("#D55E00","#0072B2")) +
  scale_linetype_manual(name = "", values = c("Non-perennial sites" = "dashed")) +
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
  guides(fill = guide_legend(title = "Treatment", order = 1),
         linetype = guide_legend(title = "", order = 2, override.aes = list(linewidth = 1))) +
  xlab(expression(atop("\n O"[2]*" Consumption Rate", "(mg O"[2]*" kg"^-1*" h"^-1*")"))) +
  ylab("Count")

## Effect size histogram WITH VERTICAL LINES FOR NON-PERENNIAL SITES
effect_limits <- c(min(effect_with_metadata$Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, na.rm = TRUE), 
                   max(effect_with_metadata$Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H, na.rm = TRUE))

effect_hist = ggplot(effect_with_metadata, aes(x = (Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H))) +
  geom_histogram(binwidth = 25, aes(fill = after_stat(x))) +
  {if(length(non_perennial_effect_values) > 0)
    geom_vline(data = vlines_effect, aes(xintercept = xintercept, linetype = line_type),
               color = "red", linewidth = 1)} +
  scale_fill_gradient2(name = "Effect Size", limits = effect_limits, low = "dodgerblue2", mid = "goldenrod2",
                       high = "firebrick2", midpoint = (max(effect_limits)+min(effect_limits))/2) +
  scale_linetype_manual(name = "", values = c("Non-perennial sites" = "dashed")) +
  theme_bw() +
  ylab("Count\n") +
  theme(legend.position = c(0.85, 0.85),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.key.width = unit(0.4, "in")) +
  guides(linetype = guide_legend(title = "", override.aes = list(linewidth = 1))) +
  xlab(expression(atop("\n Effect Size O"[2]*" Consumption Rate","(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" h"^-1*")")))

## Combine SI histograms
combined_SI_hist = ggarrange(all_hist, effect_hist, ncol = 2, labels = c("A", "B"), hjust = -1, vjust = 2.5)
ggsave("Figures/Combined_SI_Histograms.png", plot = combined_SI_hist, width = 10, height = 5, dpi = 300)

## Cube root histogram WITH VERTICAL LINES FOR NON-PERENNIAL SITES BY TREATMENT
all_cube_hist = ggplot(cube_respiration, aes(x = cube_Respiration_mg_kg)) +
  geom_histogram(position = "identity", alpha = 0.8, aes(fill = Treat)) +
  {if(nrow(vlines_cube_o2_by_treatment) > 0)
    geom_vline(data = vlines_cube_o2_by_treatment,
               aes(xintercept = xintercept, color = Treat, linetype = "Non-perennial sites"),
               linewidth = 1)} +
  scale_fill_manual(values = c("#D55E00","#0072B2")) +
  scale_color_manual(values = c("Dry" = "#D55E00", "Wet" = "#0072B2"), guide = "none") +
  scale_linetype_manual(name = "", values = c("Non-perennial sites" = "dashed")) +
  ylim(0, 87.5) +
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
  guides(fill = guide_legend(title = "Treatment", order = 1),
         linetype = guide_legend(title = "", order = 2, override.aes = list(linewidth = 1))) +
  xlab(expression(atop("\n O"[2]*" Consumption Rate" ^(1/3)*"", "(mg O"[2]*" kg"^-1*" h"^-1*")"))) +
  ylab("Count")

# Cube root effect size histogram WITH VERTICAL LINES FOR NON-PERENNIAL SITES
c# Cube root effect size histogram WITH VERTICAL LINES FOR NON-PERENNIAL SITES
cube_effect_limits <- c(-12, 12)  # Fixed limits for color scale

cube_effect_limits <- c(-12, 12)  # Match the x-axis limits

cube_effect_hist = ggplot(cube_effect, aes(x = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_histogram(binwidth = 0.5, aes(fill = after_stat(x))) +
  {if(length(non_perennial_cube_effect_values) > 0)
    geom_vline(data = vlines_cube_effect,
               aes(xintercept = xintercept, linetype = line_type),
               color = "red", linewidth = 1)} +
  scale_fill_gradient2(
    name = "Cubed Root\nEffect Size", 
    low = "firebrick2",       # Negative values = red
    mid = "goldenrod2",       # Zero = yellow
    high = "dodgerblue2",     # Positive values = blue
    midpoint = 0,             # Yellow at zero
    limits = cube_effect_limits,
    oob = scales::squish      # ADD THIS - handles out-of-bounds colors
  ) +
  scale_linetype_manual(name = "", values = c("Non-perennial sites" = "dashed")) +
  theme_bw() +
  coord_cartesian(xlim = c(-12, 12)) +  # USE THIS instead of xlim()
  ylab("Count\n") +
  theme(legend.position = c(0.85, 0.85),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.key.width = unit(0.4, "in")) +
  guides(linetype = guide_legend(title = "", override.aes = list(linewidth = 1))) +
  xlab(expression(atop("\n Effect Size O"[2]*" Consumption Rate" ^(1/3),"(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" h"^-1*")")))
## Combine cube root histograms
combined_hist = ggarrange(all_cube_hist, cube_effect_hist, ncol = 2, labels = c("A", "B"), hjust = -5, vjust = 2.5)
ggsave("Figures/FigS3Combined_Cube_Histograms.png", plot = combined_hist, width = 10, height = 5, dpi = 300)
ggsave("Figures/Figs3Combined_Cube_Histograms.pdf", plot = combined_hist, width = 10, height = 5, dpi = 300)

