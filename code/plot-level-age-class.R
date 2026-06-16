## =================================================
##
## Title: plot-lvl-counts-classes
## Author: Jasmine Williamson
## Date Created: 6/16/2026
##
## Description: Look at raw counts for age classes per plot
## 
##
## =================================================


## settings -----------------------------------------------------------------------------------------------

  # rm(list=ls())
  setwd("/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/occupancy-central/data")
  
  # library(unmarked)
  # library(ggplot2)
  # library(stats)
  # library(MASS)
  # library(tidyverse)
  library(dplyr)


## load data ----------------------------------------------------------------------------------------

  subplot <- read.csv("data/covariate matrices/subplot.complete.csv")
  sals <- read.csv("data/occupancy/sals.complete.csv", 
                   colClasses = c(landowner="factor", stand="character", trt="factor",
                                  obs="factor", subplot="factor", recap="factor",
                                  pass="factor", spp="factor", cover_obj="factor", 
                                  substrate="factor", age_class="factor"))


## df with all plots/occasions ----------------------------------------------------------------------

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
  df.new <- as.data.frame(df.new)


# merge sals with full subplot/pass df -------------------------------------------------------------
  
  sals.new <- sals[,c(1,12:18)] 
  sals.new <- filter(sals.new, spp %in% c("OSS", "ENES"))
  sals.new$detect = 1  
  
  # merge
  
  #add in sites with no detections by merging with df that has all site/subplot combos listed
  df.merge <- full_join(df.new,sals.new,by=c("site_id","subplot","pass"))
  df.merge$detect <- ifelse(is.na(df.merge$detect), 0, df.merge$detect) #make NA's = 0

  
# counts of sals -----------------------------------------------------------------------------------
  
  sal_counts <- df.merge %>%
    filter(detect == 1) %>%
    group_by(spp, age_class) %>%
    summarise(count = n(), .groups = "drop")

  sal_counts_wide <- sal_counts %>%
    unite("spp_age", spp, age_class, sep = "_") %>%
    pivot_wider(names_from = spp_age, values_from = count, values_fill = 0)
  sal_counts_wide <- sal_counts_wide[-6] #removing enes unknown age class
  sal_counts_wide <- sal_counts_wide[-7] #removing oss unknown age class
  
  sum(df.merge$detect[df.merge$spp == "OSS"], na.rm = TRUE)
  sum(df.merge$detect[df.merge$spp == "ENES"], na.rm = TRUE)
  
  
#### counts dataframe --------------------------------------------------------------------------
  
  # oss
  
  sals.oss <- subset(sals.new, spp=="OSS")
  df.merge.o <- full_join(df.new,sals.oss,by=c("site_id","subplot","pass"))
  df.merge.o$detect <- ifelse(is.na(df.merge.o$detect), 0, df.merge.o$detect) #make NA's = 0

  
  # summarize count and detections per site/subplot/pass so the rows are unique (no repliate rows)
  df.counts.o <- df.merge.o %>%
    group_by(site_id, stand, trt, year, subplot, pass) %>%
    summarise(
      count = sum(detect, na.rm = TRUE),
      detect = as.integer(sum(detect, na.rm = TRUE) > 0),
      .groups = "drop"
    ) 

  # reshape for only counts 
  counts.wide.o <-  df.counts.o %>%
    group_by(site_id,  stand, trt, year, subplot, pass) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    pivot_wider(
      names_from = pass,
      values_from = count,
      names_prefix = "V"
    )  
  
  
  
    write.csv(counts.wide.o, "data/abundance/plot-lvl-counts-o.csv", row.names = FALSE)
  
    
  
  # enes
  
  
  sals.enes <- subset(sals.new, spp=="ENES")
  df.merge.e <- full_join(df.new,sals.enes,by=c("site_id","subplot","pass"))
  df.merge.e$detect <- ifelse(is.na(df.merge.e$detect), 0, df.merge.e$detect) #make NA's = 0
  
  
  # summarize count and detections per site/subplot/pass so the rows are unique (no repliate rows)
  df.counts.e <- df.merge.e %>%
    group_by(site_id, stand, trt, year, subplot, pass) %>%
    summarise(
      count = sum(detect, na.rm = TRUE),
      detect = as.integer(sum(detect, na.rm = TRUE) > 0),
      .groups = "drop"
    ) 
  
  # reshape for only counts 
  counts.wide.e <-  df.counts.e %>%
    group_by(site_id,  stand, trt, year, subplot, pass) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    pivot_wider(
      names_from = pass,
      values_from = count,
      names_prefix = "V"
    )  
  
  
  
    write.csv(counts.wide.e, "data/abundance/plot-lvl-counts-e.csv", row.names = FALSE)
  
  
  