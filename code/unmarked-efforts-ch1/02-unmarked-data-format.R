## -------------------------------------------------------------------------------------------------------
##
## 02-unmarked-data-format.R 
##
## Rerun the occupancy models using package Unmarked from 2023 analysis
##
## Jasmine Williamson
## Date Created: 07-09-2024
##

## goals
# using the unmarked code that i started last august, rework to use new covariate matrix

## insights
# 02/27/25 tried formatting for use of csvToUMF but couldnt get function to work, gave up on it 
# 03/18/25 updated matrix for UMF using all variables for global model


## settings -----------------------------------------------------------------------------------------------

    rm(list=ls())
    #setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/oss-occu/data")
    
    library(stats)
    library(MASS)
    library(tidyverse)

    

## load data----------------------------------------------------------------------------------------------

    site <- read.csv("data/covariate matrices/env_subset_corr.csv", row.names = 1)
    trt <- read.csv("data/site.complete.csv", row.names = 1)
    subplot <- read.csv("data/subplot.complete.csv")
    
    # oss.dets <- read.csv("data/occupancy/oss.occu.wide.csv")
    # enes.dets <- read.csv("data/occupancy/enes.occu.wide.csv")
    forcounts <- readRDS("data/covariate matrices/habitat.occu.complete.rds")

  
    

## Detection covs  ------------------------------------------------------------------------------------------

#format subplot/det covs from long to wide dataframes

# soil moisture 

    #subset df with only cols i need (cant use site id here bc its a character)
    df <- subset(subplot[,c('subplot','soil_moist_avg')])
    
    #change to long format
    df_long <- df %>%
      pivot_longer(cols = -subplot, names_to = "variable", values_to = "value")
    #add site id back in
    df_long <- cbind(df_long, subplot$site_id)
    
    # Reshape the data back to wide format, flipping the dimensions
    df_wide <- df_long %>%
      pivot_wider(names_from = subplot, values_from = value)
    
    df_wide <- df_wide[,-1]
    df_wide <- as.data.frame(df_wide)
    rownames(df_wide) <- df_wide[,1]
    soil_moist <- df_wide[,-1]
    soil_moist <- round(soil_moist, 2)
    names(soil_moist) <- c("soilmoist-1","soilmoist-2","soilmoist-3","soilmoist-4","soilmoist-5","soilmoist-6","soilmoist-7")
    
    
# dwd cover
    
    df <- subset(subplot[,c('subplot','dwd_cov')])
    
    #change to long format
    df_long <- df %>%
      pivot_longer(cols = -subplot, names_to = "variable", values_to = "value")
    #add site id back in
    df_long <- cbind(df_long, subplot$site_id)
    
    # Reshape the data back to wide format, flipping the dimensions
    df_wide <- df_long %>%
      pivot_wider(names_from = subplot, values_from = value)
    
    df_wide <- df_wide[,-1]
    df_wide <- as.data.frame(df_wide)
    rownames(df_wide) <- df_wide[,1]
    dwd_cov <- df_wide[,-1]
    names(dwd_cov) <- c("dwdcov-1","dwdcov-2","dwdcov-3","dwdcov-4","dwdcov-5","dwdcov-6","dwdcov-7")
    
# time since rain
    
    df <- site[,"days_since_rain", drop=FALSE]
    rain <- cbind(df, rep(df, each = 6))
    names(rain) <- c("rain-1","rain-2","rain-3","rain-4","rain-5","rain-6","rain-7")

# temp
    
    df <- site[,"temp", drop=FALSE]
    temp <- cbind(df, rep(df, each = 6))
    names(temp) <- c("temp-1","temp-2","temp-3","temp-4","temp-5","temp-6","temp-7")
    
# treatment
    
    df <- trt[,"trt", drop=FALSE]
    treatment <- cbind(df, rep(df, each = 6))
    names(treatment) <- c("trt-1","trt-2","trt-3","trt-4","trt-5","trt-6","trt-7")
    
# year
    
    df <- trt[,"year", drop=FALSE]
    year <- cbind(df, rep(df, each = 6))
    names(year) <- c("yr-1","yr-2","yr-3","yr-4","yr-5","yr-6","yr-7")
    
    
## Site/occu covs ----------------------------------------------------------------------------------
    
    site_subset <- subset(site[,c('veg_cov', 'canopy_cov', 'soil_moist', 'fwd_cov', 'dwd_count', 'avg_volume')])
    site_subset$trt <- trt$trt
    site_subset$trt <- factor(trt$trt, 
                                  levels = c("UU", "BU", "HB", "HU", "BS"))
    site_subset$year <- trt$year
    
    # for this method, scale in the formula, not in the input data
    
    
## Sal counts/detections ---------------------------------------------------------------------
    
    counts <- forcounts[,c("site_id","subplot","OSS","ENES")]
    
    oss_counts <- with(counts, tapply(OSS, list(site_id, subplot), sum, na.rm=TRUE))
    colnames(oss_counts) <- c("OC1","OC2","OC3","OC4","OC5","OC6","OC7")    
    enes_counts <- with(counts, tapply(ENES, list(site_id, subplot), sum, na.rm=TRUE))
    colnames(enes_counts) <- c("EC1","EC2","EC3","EC4","EC5","EC6","EC7")    
    
    oss_dets <- ifelse(oss_counts > 0, 1, 0)
    colnames(oss_dets) <- c("OD1","OD2","OD3","OD4","OD5","OD6","OD7")    
    enes_dets <- ifelse(enes_counts > 0, 1, 0)
    colnames(enes_dets) <- c("ED1","ED2","ED3","ED4","ED5","ED6","ED7")    
    
    
## merge and format --------------------------------------------------------------------------------------
 
    
    # merge sal data
    df1 <- cbind(oss_counts, enes_counts)
    df1 <- cbind(df1, oss_dets)
    df1 <- cbind(df1, enes_dets)
    
    
    # merge covariates
    df2 <- cbind(site_subset, soil_moist)
    df2 <- cbind(df2, temp)
    df2 <- cbind(df2, rain)
    df2 <- cbind(df2, dwd_cov)
    df2 <- cbind(df2, treatment)
    df2 <- cbind(df2, year)
    
    
    # merge sals with covariates
    df3 <- cbind(df1, df2)
    df3$site_id <- rownames(df3)
    
    
    write.csv(df3, "data/occupancy/data-UMF-1.csv", row.names = FALSE)
    

    
    
    

    



