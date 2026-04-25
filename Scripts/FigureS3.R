# This script makes scatterplot figures for ECA physical manuscript
library(tidyverse)
library(ggpubr)

rm(list = ls()); graphics.off()

# Functions ---------------------------------------------------------------
cube_root <- function(x) sign(x) * (abs(x))^(1/3)

## Read in / merge sediment data ------------------------------------------
data <- read.csv("Data/EC_Sediment_Sample_Data_Summary.csv", skip = 2) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%
  select(-c(Field_Name, IGSN, Material)) %>%
  select(Sample_Name, Median_Respiration_Rate_mg_DO_per_kg_per_H)

ev_data <- read.csv("Data/WHONDRS_EV_Data_Package/Sample_Data/WHONDRS_EV_Sediment_Sample_Data_Summary.csv", skip = 2) %>%
  filter(grepl("EV", Sample_Name)) %>%
  select(-c(Field_Name, IGSN, Material)) %>%
  select(Sample_Name, Median_Respiration_Rate_mg_DO_per_kg_per_H)

all_data <- rbind(data, ev_data)

## Read in field metadata to get perennial/intermittent -------------------
ec_meta <- read.csv("Data/EC_Field_Metadata.csv", fileEncoding = "latin1") %>%
  select(Parent_ID, Intermittent_or_Perennial) %>%
  rename(site = Parent_ID)

ev_meta <- read.csv("Data/WHONDRS_EV_Data_Package/WHONDRS_EV_Field_Metadata.csv") %>%
  select(Parent_ID, Intermittent_or_Perennial) %>%
  rename(site = Parent_ID)

site_meta <- bind_rows(ec_meta, ev_meta) %>%
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
    treatment   = str_extract(Sample_Name, "(?<=-)[WD]$"),
    site        = str_remove(Sample_Name, "-[WD]$")
  ) %>%
  filter(!is.na(site), !is.na(treatment), !is.na(rate_mgkg_h))

paired <- resp %>%
  select(site, treatment, rate_mgkg_h) %>%
  pivot_wider(names_from = treatment, values_from = rate_mgkg_h) %>%
  filter(!is.na(W), !is.na(D)) %>%
  mutate(
    abs_diff = abs(W - D),
    wet_rate = abs(W),
    slower_treatment = case_when(
      abs(D) < abs(W) ~ "Dry",
      abs(W) < abs(D) ~ "Wet",
      TRUE             ~ "Tie"
    ),
    rel_diff = abs_diff / wet_rate
  ) %>%
  filter(slower_treatment != "Tie") %>%
  left_join(site_meta, by = "site")

## Publication-ready figure with programmatic arrows ----------------------
max_of_wet <- max(abs(paired$W), na.rm = TRUE)
x_max      <- max(paired$abs_diff, na.rm = TRUE)
arrow_col  <- "#0B6E8E"

p <- ggplot(paired, aes(x     = abs_diff,
                        y     = rel_diff,
                        color = slower_treatment,
                        shape = flow_regime)) +
  
  geom_point(size = 4, alpha = 0.75, stroke = 1.1) +
  
  geom_abline(intercept = 0,
              slope     = 1 / max_of_wet,
              linetype  = "dashed",
              linewidth = 0.5,
              color     = "grey40") +
  
  # --- Scales ---
  scale_color_manual(
    values = c("Dry" = "#4A90D9", "Wet" = "#E41A1C"),
    name   = "Slower rate from"
  ) +
  scale_shape_manual(
    values = c("Perennial" = 16, "Non-perennial" = 17),
    name   = "Flow regime"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  # --- Labels with proper units ---
  labs(
    x = expression(
      "Absolute difference between wet and dry median respiration rates (mg DO "
      * kg^{-1} ~ h^{-1} * ")"
    ),
    y = "Absolute difference / wet absolute rate"
  ) +
  
  # --- Theme ---
  theme_bw(base_size = 16) +
  theme(
    # Axis titles
    axis.title.x    = element_text(size = 15, face = "bold",
                                   margin = margin(t = 12)),
    axis.title.y    = element_text(size = 15, face = "bold",
                                   margin = margin(r = 12)),
    # Axis tick labels
    axis.text       = element_text(size = 13, color = "black"),
    
    # Legend
    legend.title      = element_text(size = 13, face = "bold"),
    legend.text       = element_text(size = 12),
    legend.position   = c(0.73, 0.45),
    # For ggplot2 >= 3.5 use instead:
    # legend.position         = "inside",
    # legend.position.inside  = c(0.73, 0.45),
    legend.background = element_rect(fill   = "white",
                                     color  = "grey60",
                                     linewidth = 0.4),
    legend.key.size   = unit(1.2, "lines"),
    legend.spacing.y  = unit(4, "pt"),
    
    # Panel
    panel.grid.minor = element_blank(),
    
    # Extra margins for the arrows
    plot.margin = margin(t = 55, r = 65, b = 12, l = 12, unit = "pt")
  ) +
  
  # Allow drawing outside the panel
  
  coord_cartesian(clip = "off") +
  
  # ── Horizontal arrow at top ──
  annotate("segment",
           x = 0,   xend = x_max * 0.95,
           y = 1.12, yend = 1.12,
           arrow     = arrow(length = unit(0.30, "cm"), type = "closed"),
           color     = arrow_col,
           linewidth = 1.5) +
  annotate("text",
           x = x_max * 0.475, y = 1.19,
           label    = "Increasingly meaningful absolute effect size",
           size     = 5,
           fontface = "bold",
           color    = arrow_col) +
  
  # ── Vertical arrow on right ──
  annotate("segment",
           x = x_max * 1.07, xend = x_max * 1.07,
           y = 0.0,          yend = 0.98,
           arrow     = arrow(length = unit(0.30, "cm"), type = "closed"),
           color     = arrow_col,
           linewidth = 1.5) +
  annotate("text",
           x = x_max * 1.15, y = 0.50,
           label    = "Increasingly meaningful relative effect",
           size     = 5,
           fontface = "bold",
           color    = arrow_col,
           angle    = 270)

print(p)

# Save at publication quality
ggsave("Figures/FigureS3_effect_size.png", p,
       width = 10, height = 7, dpi = 300, bg = "white")
ggsave("Figures/FigureS3_effect_size.pdf", p,
       width = 10, height = 7)
