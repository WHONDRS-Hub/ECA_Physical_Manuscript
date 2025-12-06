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

# Reading in metadata from both files to get river type
metadata = read.csv("Data/EC_Field_Metadata.csv") %>%
  select(c(Parent_ID, Intermittent_or_Perennial))
metadata2 = read.csv("Data/WHONDRS_EV_Data_Package/WHONDRS_EV_Field_Metadata.csv") %>%
  select(c(Parent_ID, Intermittent_or_Perennial))
metadata = rbind(metadata, metadata2)

# FIXED: MERGE ALL_DATA WITH METADATA TO IDENTIFY INTERMITTENT SITES
all_data_with_metadata = all_data %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^[^_]+_[^_]+")) %>%  # Extract "EC_005" from "EC_005_INC-D1"
  left_join(metadata, by = "Parent_ID")

# Check merge
print("All data with metadata:")
print(head(all_data_with_metadata %>% select(Sample_Name, Parent_ID, Intermittent_or_Perennial)))
print(paste("Number of intermittent samples in all_data:", sum(all_data_with_metadata$Intermittent_or_Perennial == "Intermittent", na.rm = TRUE)))

# Get O2 consumption rate values for intermittent sites
intermittent_o2_values = all_data_with_metadata %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999) %>%
  filter(Intermittent_or_Perennial == "Intermittent", !is.na(Intermittent_or_Perennial)) %>%
  pull(Respiration_Rate_mg_DO_per_kg_per_H) %>%
  as.numeric() %>%
  abs() 

print("Intermittent O2 consumption rate values:")
print(intermittent_o2_values)

# Take Cube root of all respiration rates for figure
cube_respiration = all_data_with_metadata %>%
  select(c(Sample_Name, Respiration_Rate_mg_DO_per_kg_per_H, Parent_ID, Intermittent_or_Perennial)) %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999) %>% 
  mutate(cube_Respiration_mg_kg = cube_root(abs(as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)))) %>% 
  mutate(Treat = if_else(grepl("D", Sample_Name), "Dry", "Wet"))

# Get cube root O2 consumption rate values for intermittent sites
intermittent_cube_o2_values = cube_respiration %>%
  filter(Intermittent_or_Perennial == "Intermittent", !is.na(Intermittent_or_Perennial)) %>%
  pull(cube_Respiration_mg_kg)

print("Intermittent cube root O2 consumption rate values:")
print(intermittent_cube_o2_values)

# Calculate medians for effect size calculation
median_respiration = all_data %>%
  select(-c(Respiration_R_Squared, Respiration_R_Squared_Adj, Respiration_p_value, Total_Incubation_Time_Min, Number_Points_In_Respiration_Regression, Number_Points_Removed_Respiration_Regression, DO_Concentration_At_Incubation_Time_Zero)) %>%
  mutate(across(c(SpC_microsiemens_per_cm:Respiration_Rate_mg_DO_per_kg_per_H), as.numeric)) %>%
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(grepl("INC_Method_001|INC_Method_002|INC_QA_004", Methods_Deviation), NA, Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  mutate(Respiration_Rate_mg_DO_per_kg_per_H = ifelse(Respiration_Rate_mg_DO_per_kg_per_H == "-9999", NA, Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  separate(Sample_Name, c("Sample_ID", "Rep"), sep = "-") %>%
  mutate(Rep = if_else(grepl("D", Rep), "D", "W")) %>%
  group_by(Sample_ID) %>%
  summarise(across(where(is.numeric),
                   list(Median = ~median(.x, na.rm = TRUE),
                        cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm = TRUE),
                        n = ~sum(!is.na(.x))),
                   .names = "{.fn}_{.col}")) %>%
  ungroup() %>%
  select(c(Sample_ID, Median_Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  mutate(Sample_Name = str_replace(Sample_ID, "INC", "all")) %>%
  select(-Sample_ID)

## Effect Size Data for histograms ####
effect = read.csv("Data/EC_Sediment_Effect_Size.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%
  select(c(Sample_Name, Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  mutate(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = as.numeric(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H))

# FIXED: Corrected the filter issue in effect_ev
effect_ev = read.csv("Data/WHONDRS_EV_Data_Package/Sample_Data/WHONDRS_EV_Sediment_Effect_Size.csv", skip = 2) %>%
  filter(grepl("EV", Sample_Name)) %>%
  select(c(Sample_Name, Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = Effect_Size_Median_Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  mutate(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = as.numeric(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) %>%
  filter(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H < 0)  # Fixed: moved filter outside mutate

effect = rbind(effect, effect_ev)

# FIXED: MERGE EFFECT SIZE DATA WITH METADATA TO IDENTIFY INTERMITTENT SITES
effect_with_metadata = effect %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^[^_]+_[^_]+")) %>%  # Extract "EC_005" from "EC_005_all" 
  left_join(metadata, by = "Parent_ID")

# Check merge
print("Effect size data with metadata:")
print(head(effect_with_metadata %>% select(Sample_Name, Parent_ID, Intermittent_or_Perennial)))
print(paste("Number of intermittent sites in effect data:", sum(effect_with_metadata$Intermittent_or_Perennial == "Intermittent", na.rm = TRUE)))

# Get effect size values for intermittent sites
intermittent_effect_values = effect_with_metadata %>%
  filter(Intermittent_or_Perennial == "Intermittent", !is.na(Intermittent_or_Perennial)) %>%
  pull(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)

print("Intermittent effect size values:")
print(intermittent_effect_values)

# Cube root transform effect size data
cube_effect = effect_with_metadata %>%
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

vlines_cube_o2 = data.frame(
  xintercept = intermittent_cube_o2_values,
  line_type = "Intermittent sites"
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
all_hist = all_data_with_metadata %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999) %>% 
  mutate(Treat = if_else(grepl("D", Sample_Name), "Dry", "Wet")) %>%
  ggplot(aes(x = abs(as.numeric(Respiration_Rate_mg_DO_per_kg_per_H)))) +
  geom_histogram(position = "identity", alpha = 0.8, aes(fill = Treat)) +
  # # Add vertical lines for each intermittent site's O2 consumption rate
  # {if(length(intermittent_o2_values) > 0) 
  #   geom_vline(data = vlines_o2, aes(xintercept = xintercept, linetype = line_type), 
  #              color = "red", linewidth = 1)} +
  scale_fill_manual(values = c("#D55E00","#0072B2")) +
  scale_linetype_manual(name = "", values = c("Intermittent sites" = "dashed")) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.85),
        legend.key.size = unit(0.15, "in"),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        legend.box = "vertical") +
  guides(fill = guide_legend(title = "Treatment", order = 1),
         linetype = guide_legend(title = "", order = 2)) +
  xlab(expression(atop("\n O"[2]*" Consumption Rate", "(mg O"[2]*" kg"^-1*" h"^-1*")"))) +
  ylab("Count")

## Effect size histogram WITH VERTICAL LINES FOR INTERMITTENT SITES
effect_limits <- c(-1400, 1400)
effect_hist = ggplot(effect_with_metadata, aes(x = Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_histogram(binwidth = 25, aes(fill = after_stat(x))) +
  # Add vertical lines for each intermittent site's effect size value
  {if(length(intermittent_effect_values) > 0) 
    geom_vline(data = vlines_effect, aes(xintercept = xintercept, linetype = line_type), 
               color = "red", linewidth = 1)} +
  scale_fill_gradient2(name = "Effect Size", limits = effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(effect_limits)+min(effect_limits))/2) +
  scale_linetype_manual(name = "", values = c("Intermittent sites" = "dashed")) +
  theme_bw() +
  xlim(c(-1500, 1500)) +
  ylab("Count\n") +
  theme(legend.position = c(0.15, 0.85),
        axis.title.x = element_text(size = 10)) +
  guides(linetype = guide_legend(title = "")) +
  xlab(expression(atop("\n Effect Size O"[2]*" Consumption Rate","(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" h"^-1*")")))

## Combine SI histograms
combined_SI_hist = ggarrange(all_hist, effect_hist, ncol = 2, labels = c("A", "B"), hjust = -1, vjust = 2.5)
ggsave("Figures/Combined_SI_Histograms.png", plot = combined_SI_hist, width = 10, height = 5, dpi = 300)

## Cube Root Transformed histograms WITH VERTICAL LINES FOR INTERMITTENT SITES
all_cube_hist = ggplot(cube_respiration, aes(x = cube_Respiration_mg_kg)) +
  geom_histogram(position = "identity", alpha = 0.8, aes(fill = Treat)) +
  # Add vertical lines for each intermittent site's cube root O2 consumption rate
  # {if(length(intermittent_cube_o2_values) > 0) 
  #   geom_vline(data = vlines_cube_o2, aes(xintercept = xintercept, linetype = line_type), 
  #              color = "red", linewidth = 1)} +
  scale_fill_manual(values = c("#D55E00","#0072B2")) +
  scale_linetype_manual(name = "", values = c("Intermittent sites" = "dashed")) +
  ylim(0, 87.5) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.85),
        legend.key.size = unit(0.15, "in"),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        legend.box = "vertical") +
  guides(fill = guide_legend(title = "Treatment", order = 1),
         linetype = guide_legend(title = "", order = 2)) +
  xlab(expression(atop("\n O"[2]*" Consumption Rate" ^(1/3)*"", "(mg O"[2]*" kg"^-1*" h"^-1*")"))) +
  ylab("Count")

# Cube root effect size histogram WITH VERTICAL LINES FOR INTERMITTENT SITES
cube_effect_limits <- c(-12, 12)
cube_effect_hist = ggplot(cube_effect, aes(x = cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)) +
  geom_histogram(binwidth = 0.5, aes(fill = after_stat(x))) +
  # Add vertical lines for each intermittent site's cube root effect size value
  {if(length(intermittent_cube_effect_values) > 0) 
    geom_vline(data = vlines_cube_effect, aes(xintercept = xintercept, linetype = line_type), 
               color = "red", linewidth = 1)} +
  scale_fill_gradient2(name = "Cubed Root Effect Size", limits = cube_effect_limits, low = "firebrick2", mid = "goldenrod2",
                       high = "dodgerblue2", midpoint = (max(cube_effect_limits)+min(cube_effect_limits))/2) +
  scale_linetype_manual(name = "", values = c("Intermittent sites" = "dashed")) +
  theme_bw() +
  xlim(c(-12, 12)) +
  ylab("Count\n") +
  theme(legend.position = c(0.15, 0.85),
        axis.title.x = element_text(size = 10)) +
  guides(linetype = guide_legend(title = "")) +
  xlab(expression(atop("\n Effect Size O"[2]*" Consumption Rate" ^(1/3)*"","(Median Wet - Median Dry Rate; mg O"[2]*" kg"^-1*" h"^-1*")")))

## Combine cube root histograms
combined_hist = ggarrange(all_cube_hist, cube_effect_hist, ncol = 2, labels = c("A", "B"), hjust = -5, vjust = 2.5)
ggsave("Figures/Combined_Cube_Histograms.png", plot = combined_hist, width = 10, height = 5, dpi = 300)
ggsave("Figures/Combined_Cube_Histograms.pdf", plot = combined_hist, width = 10, height = 5, dpi = 300)
