library(tidyverse)

grain = read.csv("C:/GitHub/ECA_Multireactor_Incubations/Data/v3_CM_SSS_Sediment_Sample_Data_Summary.csv", skip = 2) %>% 
  filter(grepl("CM", Sample_Name)) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "CM", "EC")) %>% 
  mutate(Sample_Name = str_replace(Sample_Name, "Sediment", "all")) %>% 
  select(c(Sample_Name, Percent_Tot_Sand, Percent_Coarse_Sand, Percent_Med_Sand, Percent_Fine_Sand, Percent_Silt, Percent_Clay, Mean_Specific_Surface_Area_m2_per_g)) %>% 
  filter(!grepl("025|026|029|030", Sample_Name))

grain_range = grain %>%
  select(-c(Percent_Tot_Sand, Mean_Specific_Surface_Area_m2_per_g)) %>%
  pivot_longer(!Sample_Name, names_to = "Fraction", values_to = "Percent") %>%
  mutate(Percent = as.numeric(Percent)) %>% 
  mutate(Size = ifelse(Fraction == "Percent_Coarse_Sand", 2, ifelse(Fraction == "Percent_Med_Sand", 0.5, ifelse(Fraction == "Percent_Fine_Sand", 0.25, ifelse(Fraction == "Percent_Silt", 0.053, 0.002))))) %>% 
  arrange(Size) %>% 
  group_by(Sample_Name) %>%
  mutate(percent_finer = cumsum(Percent)) %>%
  mutate(below_50 = ifelse(percent_finer < 50, TRUE, FALSE)) %>% 
  mutate(above_50 = ifelse(percent_finer > 50, TRUE, FALSE)) %>% 
  ungroup()



d50 = grain_range %>% 
  group_by(Sample_Name) %>% 
  mutate(below_50 = ifelse(row_number() == max(row_number()[below_50 == TRUE]), 
                       TRUE, 
                       FALSE)) %>% 
  mutate(keep = ifelse(below_50 == TRUE, TRUE, FALSE)) %>% 
  mutate(keep = ifelse(keep == TRUE, TRUE, ifelse(row_number() == min(row_number()[above_50 == TRUE]), TRUE, keep))) %>% 
  filter(keep == TRUE) %>% 
  mutate(d50 = first(Size[below_50]) +
            (50 - first(percent_finer[below_50])) *
            (first(Size[above_50]) - first(Size[below_50])) /
           (first(percent_finer[above_50]) - first(percent_finer[below_50]))) %>% 
  na.omit() %>% 
  select(c(Sample_Name, d50)) %>% 
  distinct(Sample_Name, .keep_all = TRUE)

png(file = paste0("C:/Github/ECA_Multireactor_Incubations/Physical_Manuscript_Figures/", as.character(Sys.Date()),"_D50_Range.png"), width = 12, height = 12, units = "in", res = 300)

ggplot(grain_range, mapping=aes(y = percent_finer, x = log10(Size))) + 
  geom_point() +
  geom_line() + 
  facet_wrap(~ Sample_Name) +
  geom_point(d50, color = "red", mapping = aes(y = 50, x = log10(d50)))

dev.off()

write.csv(d50, "C:/Github/ECA_Multireactor_Incubations/Data/D50_Calculations.csv", row.names = FALSE)