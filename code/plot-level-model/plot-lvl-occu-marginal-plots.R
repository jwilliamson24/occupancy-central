## =================================================
##
## Title: plot-lvl-occu-marginal-plots
## Author: Jasmine Williamson 
## Date Created: 07/25/2025
##
## Description: Visualize and interpret outputs for 
## plot-level occupancy model, both spp, 2023-2024 dataset
##
## =================================================



# Load packages and data
  library(ggplot2)
  library(dplyr)
  
  load("data/msc model rdata/msc-enes-data-workspace.RData")
  load("data/msc model rdata/multiscale_output_and_data_072125_enes_small.RData")
  E = a2
  load("data/msc model rdata/multiscale_output_072525_oss_small.RData")
  O = a2
  
# Estimates
  summary(E)
  summary(O)

