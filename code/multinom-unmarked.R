## =================================================
##
## Title: multinom-unmarked
## Author: Jasmine Williamson
## Date Created: 6/17/2026
##
## Description: Plot level multinomial n-mixture model in unmarked
## 
##
## =================================================


## settings -----------------------------------------------------------------------------------------------

  rm(list=ls())
  setwd("/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/occupancy-central")
  
  library(unmarked)
  library(dplyr)


## load data ----------------------------------------------------------------------------------------

  subplot <- read.csv("data/covariate matrices/subplot.complete.new.csv") # site covs from here
  enes.df <- read.csv("data/abundance/pass-level-counts-e.csv") # counts at pass level
  oss.df <- read.csv("data/abundance/pass-level-counts-o.csv") # counts at pass level
  weather <- read.csv("data/covariate matrices/site_aspect_precip_all_vars.csv") # precip and days since rain

  
## format loaded data 
  
  # subset
  subplot <- subplot[, c("site_id", "stand", "trt", "year", "subplot", "lat", "long", "elev", "soil_moist_avg", 
                       "jul_date", "date", "dwd_count", "decay_cl", "char_cl", "temp")] 
  
  # set decay cl and char cl NA's to zero
  subplot <- subplot %>%
    mutate(
      decay_cl = as.numeric(ifelse(is.na(decay_cl), 0, decay_cl)),
      char_cl = as.numeric(ifelse(is.na(char_cl), 0, char_cl))
    )

  # subset
  weather <- weather[, c("site_id", "stand", "year", "precip_mm", "days_since_rain")]
  
  # expand from site to plot lvl
  weather <- weather %>%
    slice(rep(1:n(), each = 7)) %>%
    group_by(site_id) %>%
    mutate(subplot = 1:7) %>%
    ungroup()
  
  
## format model data ---------------------------------------------------------------------------------

  # make counts wide format
  counts.wide.o <-  oss.df %>%
    group_by(site_id,  stand, trt, year, subplot, pass) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    pivot_wider(
      names_from = pass,
      values_from = count,
      names_prefix = "C"
    )  

  counts.wide.e <-  enes.df %>%
    group_by(site_id,  stand, trt, year, subplot, pass) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    pivot_wider(
      names_from = pass,
      values_from = count,
      names_prefix = "C"
    ) 
  
  
  # make obs wide format
  obs.wide.o <- oss.df %>%
    distinct(site_id, stand, year, subplot, pass, obs) %>%
    pivot_wider(
      names_from = pass,
      values_from = obs,
      names_prefix = "O"
    )
  
  obs.wide.e <- enes.df %>%
    distinct(site_id, stand, year, subplot, pass, obs) %>%
    pivot_wider(
      names_from = pass,
      values_from = obs,
      names_prefix = "O"
    )
  
  
  # join counts + obs + site covariates into one aligned dataframe for each spp
  model_data_o <- counts.wide.o %>%
    left_join(obs.wide.o, by = c("site_id", "stand", "year", "subplot")) %>%
    left_join(subplot, by = c("site_id", "subplot", "year", "stand", "trt")) %>%
    left_join(weather, by = c("site_id", "stand", "subplot", "year"))
  
  model_data_e <- counts.wide.e %>%
    left_join(obs.wide.e, by = c("site_id", "stand", "year", "subplot")) %>%
    left_join(subplot, by = c("site_id", "subplot", "year", "stand", "trt")) %>%
    left_join(weather, by = c("site_id", "stand", "subplot", "year"))
  
  
  # scale continuous covs
  model_data_o <- model_data_o %>%
    mutate(across(c(soil_moist_avg, dwd_count, decay_cl, temp, precip_mm, days_since_rain),
                  ~ as.numeric(scale(.)), .names = "{.col}_z"))
  
  model_data_e <- model_data_e %>%
    mutate(across(c(soil_moist_avg, dwd_count, decay_cl, temp, precip_mm, days_since_rain),
                  ~ as.numeric(scale(.)), .names = "{.col}_z"))
  
  
  # sanity check 
  nrow(counts.wide.o) # 889
  nrow(model_data_o) # 889
  sum(is.na(model_data_o$O1))  # 0
  nrow(counts.wide.e) # 889
  nrow(model_data_e) # 889
  sum(is.na(model_data_e$O1))  # 0
  
  
#### define model components -----------------------------------------------------  
  
  # abundance
  y_oss <- as.matrix(model_data_o[, c("C1", "C2", "C3")])
  y_enes <- as.matrix(model_data_e[, c("C1", "C2", "C3")])

  # site covs
  site_covs_o <- model_data_o %>%
    select(trt, soil_moist_avg_z, dwd_count_z, decay_cl_z, temp_z, precip_mm_z, days_since_rain_z)
  site_covs_e <- model_data_e %>%
    select(trt, soil_moist_avg_z, dwd_count_z, decay_cl_z, temp_z, precip_mm_z, days_since_rain_z)
  
  site_covs_o$trt <- relevel(as.factor(site_covs_o$trt), ref = "UU")
  site_covs_e$trt <- relevel(as.factor(site_covs_e$trt), ref = "UU")
  

  # obs covs
  
  # take out CH and LS obs, too few observations by them messing up model
  model_data_o <- model_data_o %>%
    mutate(across(c(O1, O2, O3), ~ ifelse(. %in% c("CH", "LS"), "Other", .)))
  
  model_data_e <- model_data_e %>%
    mutate(across(c(O1, O2, O3), ~ ifelse(. %in% c("CH", "LS"), "Other", .)))
  
  # rebuild obs_covs_o and umf_o with the updated observer matrix
  obs_covs_o <- list(observer = as.matrix(model_data_o[, c("O1", "O2", "O3")]))
  obs_covs_e <- list(observer = as.matrix(model_data_e[, c("O1", "O2", "O3")]))

  
  
#### make umf objects ------------------------------------------------------------  
  
  # oss
  
  umf_o <- unmarkedFrameMPois(
    y = y_oss,
    siteCovs = site_covs_o,
    obsCovs = obs_covs_o,
    type = "removal"
  )
  
  # set JW as obs reference
  umf_o@obsCovs$observer <- relevel(umf_o@obsCovs$observer, ref = "JW")
  
  summary(umf_o)
  
  
  # enes
  
  umf_e <- unmarkedFrameMPois(
    y = y_enes,
    siteCovs = site_covs_e,
    obsCovs = obs_covs_e,
    type = "removal"
  )
  
  # set JW as obs reference
  umf_e@obsCovs$observer <- relevel(umf_e@obsCovs$observer, ref = "JW")
  
  summary(umf_e)
  
  
  
#### oss model ---------------------------------------------------------------------------
  
  # output_name <- multinomPois(~detection formula ~abundance formula, umf_o/e)
  
  # null model - intercept only
  fm_null_o <- multinomPois(~1 ~1, umf_o)
  summary(fm_null_o)
  
  exp(-0.359)      # mean abundance per subplot 0.70
  plogis(-1.66)    # detection probability per pass 16%
  
  # no errors with running null
  
  
  
  # full model with hypothesized covs
  fm_full_o <- multinomPois(~temp_z + precip_mm_z + days_since_rain_z + observer # detection formula
                          ~trt + soil_moist_avg_z + dwd_count_z + decay_cl_z,  # abundance formula
                          umf_o)
  summary(fm_full_o)
  
  
  
#### enes model --------------------------------------------------------------------------
  
  # null model - intercept only
  fm_null_e <- multinomPois(~1 ~1, umf_e)
  summary(fm_null_e)
  
  exp(-1.17)      # mean abundance per subplot 0.31
  plogis(-1.35)    # detection probability per pass 20%
  
  # no errors with running null
  
  
  # full model - all hypothesized covs
  fm_full_e <- multinomPois(~temp_z + precip_mm_z + days_since_rain_z + observer # detection formula
                            ~trt + soil_moist_avg_z + dwd_count_z + decay_cl_z,  # abundance formula
                            umf_e)
  summary(fm_full_e)
  
  
  # getting warning message: large or missing SE values, be very cautious using these results
  # enes sample size is much smaller than oss
  # instead of running ideal full model, may need to start with null and ramp up the covs sequentially
  
  
  
  
  
  
  
  
  
  
  
  
  
  