## -------------------------------------------------------------------------------------------------------
##
## 04-unmarked-output-plot.R 
##
## Prediction outputs and plots from unmarked models using 2023-2024 data
##
## Jasmine Williamson
## Date Created: 02-27-2025
##
## 
## goals
## 

## insights
## 


## settings -----------------------------------------------------------------------------------------------

    rm(list=ls())
    #setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/oss-occu/data")
    
    library(unmarked)
    library(ggplot2)
    library(stats)
    library(MASS)
    library(tidyverse)


## load data----------------------------------------------------------------------------------------------

    df_oss <- read.csv("data/occupancy/oss_forUMF.csv")


## Extract treatment occupancy predictions  ------------------------------------------------------------

    # Extract psi predictions from model
    preds <- predict(m1, type="state")
    sites_trt <- sitecovs[, c("site_id","trt")]
    site_preds <- cbind(sites_trt, preds)
    
    confint(m1, type='det', method = 'normal')
    
    
    
    
 ## plot --------------------------------------------------------------------------------------------------
    
    # Barplot of predictions for each treatment
    ggplot(site_preds, aes(x = trt, y = Predicted, ymin = lower, ymax = upper, fill = trt)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.25, position = position_dodge(width = 0.7), linewidth = 1) +
      labs(title = "OSS: Treatment vs. Predicted Occupancy",
           x = "Treatment", y = "Predicted Value") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_brewer(palette = "Set1")
    
    #ggsave("unmarked_occu_barplot_oss.png", 
    #  path = "C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/Occupancy-2023/figures/Unmarked_occu_TWS_2023")
    
    
    # Boxplot of preds with confidence intervals
    ggplot(site_preds, aes(x = trt, y = Predicted, fill = trt)) +
      geom_boxplot() +
      geom_point(position = position_dodge(width = 0.75), size = 3) +
      geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.25, position = position_dodge(width = 0.75), size = 1) +
      labs(title = "OSS: Treatment vs. Predicted Occupancy",
           x = "Treatment", y = "Predicted Value") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_brewer(palette = "Set1")
    
    #ggsave("unmarked_occu_pt_est_oss.png", 
    # path = "C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/Occupancy-2023/figures/Unmarked_occu_TWS_2023")
    
    
    



