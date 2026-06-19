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

  # rm(list=ls())
  setwd("/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/occupancy-central")
  
  library(unmarked)
  library(dplyr)


## load data ----------------------------------------------------------------------------------------

  subplot <- read.csv("data/covariate matrices/subplot.complete.csv") # site covs from here
  enes.df <- read.csv("data/abundance/pass-level-counts-e.csv") # counts at pass level
  oss.df <- read.csv("data/abundance/pass-level-counts-o.csv") # counts at pass level
  
  
## format data for model -----------------------------------------------------------------------------

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
  
  
  
  # join counts + obs + site covariates into one aligned dataframe
  model_data_o <- counts.wide.o %>%
    left_join(obs.wide.o, by = c("site_id", "stand", "year", "subplot")) %>%
    left_join(subplot, by = c("site_id", "subplot", "year", "stand", "trt"))
  
  model_data_e <- counts.wide.e %>%
    left_join(obs.wide.e, by = c("site_id", "stand", "year", "subplot")) %>%
    left_join(subplot, by = c("site_id", "subplot", "year", "stand", "trt"))
  
  
  
  
  
  