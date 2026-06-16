## =================================================
##
## Title: occupancy-matrices-updated2025
## Author: Jasmine Williamson
## Date Created: 5/28/2025
##
## Description: Updates code for occupancy matrices for multi scale model
## These have rows at the subplot level and columns at the repeat survey (1-3) level
##
## =================================================


## settings -----------------------------------------------------------------------------------------------

 # rm(list=ls())
  #setwd("/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/oss-occu/data")
  
  library(unmarked)
  library(ggplot2)
  library(stats)
  library(MASS)
  library(tidyverse)
  library(dplyr)


## load data ----------------------------------------------------------------------------------------------
  
  subplot <- read.csv("data/covariate matrices/subplot.complete.csv")
  sals <- read.csv("data/occupancy/sals.complete.csv", 
                   colClasses = c(landowner="factor", stand="character", trt="factor",
                                  obs="factor", subplot="factor", recap="factor",
                                  pass="factor", spp="factor", cover_obj="factor", 
                                  substrate="factor", age_class="factor"))
  
  
## df for merge ------------------------------------------------------------------------------------------
  
  #creating occupancy df with non-detections
  df.new <- subplot[,c(1,3,5,6,9)] # siteID, stand, subplots, year, treatment
  df.new$subplot <- as.factor(df.new$subplot)
  df.new <- df.new[order(df.new$site_id, df.new$subplot),] # reorder
  
  # add 3 passes for each subplot
  df.new <- df.new %>%
    slice(rep(1:n(), each = 3)) %>%  # duplicate each row 3 times
    group_by(site_id, subplot) %>%
    mutate(pass = 1:3) %>%
    ungroup()
  df.new$pass <- as.factor(df.new$pass)

    
# oss matrix ---------------------------------------------------------------------------------------
  
  #remove recaptures
  #sals <- sals[sals$recap != 1, ]
  
  #subset of oss with site, subplot, spp, pass
  sals.oss <- subset(sals, spp=="OSS")
  sals.oss <- sals.oss[,c(1,12,13,14)] 
  sals.oss$detect <- 1

# merge
  
  #add in sites with no detections by merging with df that has all site/subplot combos listed
  df.merge.o <- full_join(df.new,sals.oss,by=c("site_id","subplot","pass"))
  df.merge.o$detect <- ifelse(is.na(df.merge.o$detect), 0, df.merge.o$detect) #make NA's = 0
  df.merge.o <- df.merge.o[-7]  
  df.merge.o <- as.data.frame(df.merge.o)
  
  # summarize count and detections per site/subplot/pass so the rows are unique (no repliate rows)
  df.sum.o <- df.merge.o %>%
    group_by(site_id, stand, trt, year, subplot, pass) %>%
    summarise(
      count = sum(detect, na.rm = TRUE),
      detect = as.integer(sum(detect, na.rm = TRUE) > 0),
      .groups = "drop"
    )
  
  # reshape to wide format for count and detections
  # df.wide.o <- df.sum.o %>%
  #   pivot_wider(
  #     id_cols = c(site_id, subplot),
  #     names_from = pass,
  #     values_from = c(count, detect),
  #     names_glue = "{.value}_pass{pass}"
  #   )
  
  # reshape for only detections 
  dets.o <-  df.sum.o %>%
    group_by(site_id,  stand, trt, year, subplot, pass) %>%
    summarise(detect = max(detect), .groups = "drop") %>%
    pivot_wider(
      names_from = pass,
      values_from = detect,
      names_prefix = "V"
    )
  
  dets.o <- as.data.frame(dets.o)
  
  write.csv(dets.o, "data/occupancy/dets.o.post.csv", row.names = FALSE)

# enes matrix ---------------------------------------------------------------------------------------
  
  #subset of enes with site, subplot, spp, pass
  sals.enes <- subset(sals, spp=="ENES")
  sals.enes <- sals.enes[,c(1,12,13,14)] 
  sals.enes$detect <- 1
  

# merge
  
  #add in sites with no detections by merging with df that has all site/subplot combos listed
  df.merge.e <- full_join(df.new,sals.enes,by=c("site_id","subplot","pass"))
  df.merge.e$detect <- ifelse(is.na(df.merge.e$detect), 0, df.merge.e$detect) #make NA's = 0
  df.merge.e <- df.merge.e[-7]  
  df.merge.e <- as.data.frame(df.merge.e)
  
  # summarize count and detections per site/subplot/pass so the rows are unique (no repliate rows)
  df.sum.e <- df.merge.e %>%
    group_by(site_id,  stand, trt, year, subplot, pass) %>%
    summarise(
      count = sum(detect, na.rm = TRUE),
      detect = as.integer(sum(detect, na.rm = TRUE) > 0),
      .groups = "drop"
    )
  
  # reshape for only detections 
  dets.e <-  df.sum.e %>%
    group_by(site_id,  stand, trt, year, subplot, pass) %>%
    summarise(detect = max(detect), .groups = "drop") %>%
    pivot_wider(
      names_from = pass,
      values_from = detect,
      names_prefix = "V"
    )
  
  dets.e <- as.data.frame(dets.e)
  
  write.csv(dets.e, "data/occupancy/dets.e.post.csv", row.names = FALSE)
  
  