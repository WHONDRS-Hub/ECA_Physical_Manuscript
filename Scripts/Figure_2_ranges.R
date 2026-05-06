# =========================================================
# Figure 2 cube-root histograms with RANGE BARS
# instead of individual vertical lines
# =========================================================

library(tidyverse)
library(ggpubr)

rm(list = ls())
graphics.off()

# ---------------------------------------------------------
# Functions
# ---------------------------------------------------------

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

# ---------------------------------------------------------
# Read in respiration data
# ---------------------------------------------------------

all_data = read.csv(
  "Data/EC_Sediment_SpC_pH_Temp_Respiration.csv",
  skip = 2
) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999) %>%
  select(-c(Field_Name, IGSN, Material))

ev_data = read.csv(
  "Data/WHONDRS_EV_Data_Package/Sample_Data/WHONDRS_EV_Sediment_SpC_pH_Temp_Respiration.csv",
  skip = 2
) %>%
  filter(grepl("EV", Sample_Name)) %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != -9999) %>%
  select(-c(Field_Name, IGSN, Material))

all_data = rbind(all_data, ev_data)

all_data = all_data %>%
  mutate(
    Respiration_Rate_mg_DO_per_kg_per_H =
      as.numeric(Respiration_Rate_mg_DO_per_kg_per_H) * -1
  )

# ---------------------------------------------------------
# Metadata
# ---------------------------------------------------------

metadata = read.csv("Data/EC_Field_Metadata.csv") %>%
  select(c(Parent_ID, Intermittent_or_Perennial))

metadata2 = read.csv(
  "Data/WHONDRS_EV_Data_Package/WHONDRS_EV_Field_Metadata.csv"
) %>%
  select(c(Parent_ID, Intermittent_or_Perennial))

metadata = rbind(metadata, metadata2)

all_data_with_metadata = all_data %>%
  mutate(
    Parent_ID = str_extract(Sample_Name, "^[^_]+_[^_]+")
  ) %>%
  left_join(metadata, by = "Parent_ID")

# ---------------------------------------------------------
# Effect size data
# ---------------------------------------------------------

effect = read.csv(
  "Data/EC_Sediment_Effect_Size.csv",
  skip = 2
) %>%
  filter(grepl("EC", Sample_Name)) %>%
  filter(!grepl("EC_011|EC_012|EC_023|EC_052|EC_053|EC_057", Sample_Name)) %>%
  select(
    Sample_Name,
    Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H
  ) %>%
  mutate(
    Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H =
      as.numeric(
        Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H
      ) * -1
  )

effect_ev = read.csv(
  "Data/WHONDRS_EV_Data_Package/Sample_Data/WHONDRS_EV_Sediment_Effect_Size.csv",
  skip = 2
) %>%
  filter(grepl("EV", Sample_Name)) %>%
  select(
    Sample_Name,
    Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H =
      Effect_Size_Median_Respiration_Rate_mg_DO_per_kg_per_H
  ) %>%
  mutate(
    Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H =
      as.numeric(
        Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H
      ) * -1
  )

effect = rbind(effect, effect_ev)

effect_with_metadata = effect %>%
  mutate(
    Parent_ID = str_extract(Sample_Name, "^[^_]+_[^_]+")
  ) %>%
  left_join(metadata, by = "Parent_ID")

# ---------------------------------------------------------
# Cube-root transformed respiration data
# ---------------------------------------------------------

cube_respiration = all_data_with_metadata %>%
  filter(Respiration_Rate_mg_DO_per_kg_per_H != 9999) %>%
  mutate(
    cube_Respiration_mg_kg =
      cube_root(
        as.numeric(
          Respiration_Rate_mg_DO_per_kg_per_H
        )
      ),
    Treat = if_else(
      grepl("-D", Sample_Name),
      "Dry",
      "Wet"
    )
  )

# ---------------------------------------------------------
# Cube-root transformed effect size data
# ---------------------------------------------------------

cube_effect = effect_with_metadata %>%
  mutate(
    cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H =
      cube_root(
        Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H
      )
  )

# ---------------------------------------------------------
# RANGE DATA FOR NON-PERENNIAL STREAMS
# ---------------------------------------------------------

# Panel A: wet and dry ranges

intermittent_cube_o2_by_treatment =
  cube_respiration %>%
  filter(
    Intermittent_or_Perennial == "Intermittent",
    !is.na(Intermittent_or_Perennial)
  ) %>%
  select(
    Parent_ID,
    cube_Respiration_mg_kg,
    Treat
  )

range_cube_o2_by_treatment =
  intermittent_cube_o2_by_treatment %>%
  group_by(Treat) %>%
  summarise(
    xmin = min(cube_Respiration_mg_kg, na.rm = TRUE),
    xmax = max(cube_Respiration_mg_kg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y = if_else(Treat == "Wet", 84, 78)
  )

# Panel B: effect size range

intermittent_cube_effect_values =
  cube_effect %>%
  filter(
    Intermittent_or_Perennial == "Intermittent",
    !is.na(Intermittent_or_Perennial)
  ) %>%
  pull(
    cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H
  )

range_cube_effect = data.frame(
  xmin = min(intermittent_cube_effect_values, na.rm = TRUE),
  xmax = max(intermittent_cube_effect_values, na.rm = TRUE),
  y = 7
)

# ---------------------------------------------------------
# PANEL A
# ---------------------------------------------------------

all_cube_hist =
  ggplot(
    cube_respiration,
    aes(x = cube_Respiration_mg_kg)
  ) +
  
  geom_histogram(
    position = "identity",
    alpha = 0.8,
    aes(fill = Treat)
  ) +
  
  # RANGE BARS
  
  geom_segment(
    data = range_cube_o2_by_treatment,
    aes(
      x = xmin,
      xend = xmax,
      y = y,
      yend = y,
      color = Treat
    ),
    linewidth = 2.5,
    lineend = "round",
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  
  annotate(
    "text",
    x = -Inf,
    y = 84,
    label = "Wet non-perennial range",
    hjust = -0.05,
    vjust = -0.4,
    size = 3.2
  ) +
  
  annotate(
    "text",
    x = -Inf,
    y = 78,
    label = "Dry non-perennial range",
    hjust = -0.05,
    vjust = -0.4,
    size = 3.2
  ) +
  
  scale_fill_manual(
    name = "Treatment",
    values = c(
      "Dry" = "#D55E00",
      "Wet" = "#0072B2"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "Dry" = "#D55E00",
      "Wet" = "#0072B2"
    ),
    guide = "none"
  ) +
  
  coord_cartesian(ylim = c(0, 87.5)) +
  
  theme_bw() +
  
  theme(
    legend.position = c(0.75, 0.85),
    legend.key.size = unit(0.15, "in"),
    legend.key.width = unit(0.4, "in"),
    legend.title = element_text(size = 8),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.box = "vertical"
  ) +
  
  xlab(
    expression(
      atop(
        "\n O"[2] * " Consumption Rate"^(1/3),
        "(mg O"[2] * " kg"^-1 * " h"^-1 * ")"
      )
    )
  ) +
  
  ylab("Count")

# ---------------------------------------------------------
# PANEL B
# ---------------------------------------------------------

cube_effect_limits <- c(-12, 12)

cube_effect_hist =
  ggplot(
    cube_effect,
    aes(
      x =
        cube_Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H
    )
  ) +
  
  geom_histogram(
    binwidth = 0.5,
    boundary = 0,
    aes(fill = after_stat(ifelse(x < 0, NA, x)))
  ) +
  
  # RANGE BAR
  
  geom_segment(
    data = range_cube_effect,
    aes(
      x = xmin,
      xend = xmax,
      y = y,
      yend = y
    ),
    color = "black",
    linewidth = 2.5,
    lineend = "round",
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  
  annotate(
    "text",
    x = -12,
    y = 8,
    label = "Non-perennial range",
    hjust = 0,
    vjust = -0.5,
    size = 3.2,
    color = "black"
  ) +
  
  scale_fill_gradient2(
    name = "Cubed Root\nEffect Size",
    limits = c(0, 12),
    low = "dodgerblue2",
    mid = "goldenrod2",
    high = "firebrick2",
    midpoint = 6,
    na.value = "grey70",
    guide = "none"
  ) +
  
  scale_x_continuous(
    breaks = c(-12, -8, -4, 0, 4, 8, 12)
  ) +
  
  coord_cartesian(
    xlim = cube_effect_limits
  ) +
  
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
  
  xlab(
    expression(
      atop(
        "\n Effect Size O"[2] *
          " Consumption Rate"^(1/3),
        "(Median Wet - Median Dry Rate; mg O"[2] *
          " kg"^-1 * " h"^-1 * ")"
      )
    )
  )
# ---------------------------------------------------------
# Combine panels
# ---------------------------------------------------------

combined_hist =
  ggarrange(
    all_cube_hist,
    cube_effect_hist,
    ncol = 2,
    labels = c("A", "B"),
    hjust = -5,
    vjust = 2.5
  )

# ggsave(
#   "Figures/Figure_2_Combined_Cube_Histograms.png",
#   plot = combined_hist,
#   width = 10,
#   height = 5,
#   dpi = 300
# )

ggsave(
  "Figures/Figure_2_Combined_Cube_Histograms_range.pdf",
  plot = combined_hist,
  width = 10,
  height = 5,
  dpi = 300
)

# ---------------------------------------------------------
# Counts for Figure 2 caption
# ---------------------------------------------------------

# Panel A: rate estimates by treatment and stream type
panelA_counts <- cube_respiration %>%
  filter(!is.na(Intermittent_or_Perennial)) %>%
  count(Treat, Intermittent_or_Perennial, name = "n_rate_estimates")

print("Panel A counts: rate estimates by treatment and stream type")
print(panelA_counts)

# Panel B: effect-size estimates by stream type
panelB_counts <- cube_effect %>%
  filter(!is.na(Intermittent_or_Perennial)) %>%
  count(Intermittent_or_Perennial, name = "n_effect_size_estimates")

print("Panel B counts: effect-size estimates by stream type")
print(panelB_counts)
