# ==============================================================================
#
# Create input file to map
# 
# Status: COMPLETE
#
# ==============================================================================
#
# Author: Brieanne Forbes (brieanne.forbes@pnnl.gov)
# 13 September 2024

# remove all files
rm(list=ls(all=TRUE))

#set working directory to parent github folder
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
setwd("./../..")

# ==============================================================================

library(tidyverse)

# ================================= User inputs ================================

data_dir <- './Data'

# ============================ find and read files =============================

metadata <- list.files(data_dir, 'Metadata', full.names = T) %>%
  read_csv() %>%
  select(Parent_ID, Sample_Latitude, Sample_Longitude) %>% # select columns with coordinates 
  rename(Latitude = Sample_Latitude, # rename for simplicity
         Longitude = Sample_Longitude)

data <- list.files(data_dir, 'Effect_Size', full.names = T) %>%
  read_csv(skip = 2, na =c('', NA, 'N/A', '-9999')) %>% # skip metadata rows and replace all different missing value codes with NA
  filter(!is.na(Sample_Name)) %>% # filter out metadata rows
  mutate(Parent_ID = str_remove(Sample_Name, '_all'), # create Parent_ID column to match to metadata
         Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H = as.numeric(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H) * -1) %>%  #make values positive
select(Parent_ID, Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H) #select column used for map

input_file <- metadata %>%
  full_join(data) %>%
  filter(!Parent_ID %in% c('EC_011','EC_012', 'EC_023','EC_052','EC_053','EC_057')) # remove Parent_IDs not used in the analysis

write_csv(input_file, './Physical_Manuscript_Figures/Map/Map_Input_File.csv')

# This input file will not be included in the data package. Rerun this script in order to 
# reproduce the map. Use the "Repair Data Source" function to select the file created by this script. 
