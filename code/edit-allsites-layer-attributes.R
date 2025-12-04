## =================================================
##
## Title: edit-allsites-layer-attributes
## Author: Jasmine Williamson 
## Date Created: 11/12/2025
##
## Description: editing attribute table for allsites_2024_joined shapefile so that
## each site has a treatment designation
##
## =================================================

## some sites have multiple trt designations for 2023-2024, will just use the most recent one used

## Settings, load
  library(dplyr)

  att_table_all <- read.csv("/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/oss-occu/data/allsites_2024_joined_layer_attribute_table.csv")
  site_level_matrix <- read.csv("/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/oss-occu/data/covariate matrices/site_level_matrix.csv")
  
  
## edit, subset existing attribute table cols  
  att_table_all <- att_table_all %>%
    mutate(Srvy_13_19 = ifelse(
      rowSums(select(., Srvy_2013:Srvy_2019), na.rm = TRUE) > 0, 
      1, 
      0
    ))
  
  att_table_all <- att_table_all %>%
    mutate(Srvy_23 = ifelse(X2023_Surve == "y", 1,0))
  
  att_table <- att_table_all[, -c(4:9,28,30:37)]
  names(att_table)[21] <- "Site_ID" 
  head(att_table)
  
  
## add col saying what sites were surveyed in 2024
  
  # Get list of stands surveyed in 2024
  stands_2024 <- site_level_matrix %>%
    filter(year == 2024) %>%
    pull(stand) %>%
    unique()
  
  # Add Srvy_24 column to att_table
  att_table <- att_table %>%
    mutate(Srvy_24 = ifelse(Site_ID %in% stands_2024, 1, 0))
  


  write.csv(att_table, "att-table-edited-111225.csv", row.names = FALSE)  
  
  
  
  
    
