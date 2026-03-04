rm(list=ls(all=T))
library(tidyverse)
library(dplyr) # For reorganization
library(stringr) # For string manipulation
# ==== Defining paths and working directories ======
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Geospatial_variables/'
# ====== Read in data ======
comids =  read.csv("C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Geospatial_variables/Example_Code/v4_RCSFA_Geospatial_Data_Package/v4_RCSFA_Geospatial_Site_Information.csv")%>%
dplyr::select(site = Site_ID, comid = COMID)

sample_data = read_csv(paste0('EC_Data_Package/EC_Field_Metadata.csv'),comment = '#', na = c('N/A', -9999)) %>%
 dplyr::select(-CM_Parent_ID,-IGSN_Sample_Name,-IGSN,-Sample_Date,-Sample_Latitude,
               -Sample_Longitude, -miniDOT_Date, -miniDOT_Latitude, -miniDOT_Longitude,
               -miniDOT_Start_Time, -miniDOT_End_Time, -Time_Zone, -miniDOT_SN,
               -Contact_First_Name, -Contact_Last_Name, -Organization, -miniDOT_Notes,               
               -Additional_Sampling_Notes,-Hydrograph_Online, -Site_Name)
  
  
names(sample_data)[3] = 'site'

geospatial = read.csv(paste0(github,'Archived_versions/v4_RCSFA_Extracted_Geospatial_Data_2025-01-31.csv'))%>%
  mutate(PctFst = pctmxfst2019ws + pctdecid2019ws + pctconif2019ws,
         PctAg = pctcrop2019ws + pcthay2019ws) %>%
  dplyr::select(site,stream_order = streamorde,total_drainage_area_sq_km = totdasqkm, slope,elevation = elevws,AridityWs,Pct_shrub= pctshrb2019ws,precipitation = PrecipSite,PctFst,PctAg) %>%
  distinct(site, .keep_all = TRUE)



# ==== Merge data  ===
data_merged = sample_data %>%
  left_join(geospatial, by='site')
# ==== Save data ====
write.csv(data_merged, 'data/EC_Field_metadata_and_Geospatial.csv', row.names = FALSE)
