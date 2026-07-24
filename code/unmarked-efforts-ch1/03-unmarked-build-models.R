## -------------------------------------------------------------------------------------------------------
##
## 03-unmarked-build-models.R 
##
## Fit models in unmarked using dataframe built in previous script for 2023-2024 data
##
## Jasmine Williamson
## Date Created: 02-27-2025
##
## 
## goals
## run model selection for different dwd metrics
## use chosen dwd metric and run model selection for different detection covariates
## use chosen dwd and detection covs to run large model selection

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

    df_oss <- read.csv("occupancy/oss_forUMF_scaled.csv")
    df_oss <- read.csv("occupancy/oss_forUMF_unscaled.csv")


## Build OSS unmarkedFrameOccu Object -----------------------------------------------------------------------

    # scaled
    UMF.oss.1 <- unmarkedFrameOccu(
      y = df_oss[, grep("^X", names(df_oss))],  # selects cols starting with X
      siteCovs = df_oss[, c("veg_cov","dwd_cov","fwd_cov","dwd_count","decay_cl","avg_volume","trt")],
      obsCovs = list(
        temp = df_oss[, grep("temp.", names(df_oss))], # selects cols including temp-
        soilmoist = df_oss[, grep("soilmoist.", names(df_oss))],
        rain = df_oss[, grep("rain.", names(df_oss))]
      )
    )


    # unscaled
    UMF.oss.2 <- unmarkedFrameOccu(
      y = df_oss_unsc[, grep("^X", names(df_oss_unsc))],  # selects cols starting with X
      siteCovs = df_oss_unsc[, c("veg_cov","dwd_cov","fwd_cov","dwd_count","decay_cl","avg_volume","trt")],
      obsCovs = list(
        temp = df_oss_unsc[, grep("temp.", names(df_oss_unsc))], # selects cols including temp-
        soilmoist = df_oss_unsc[, grep("soilmoist.", names(df_oss_unsc))],
        rain = df_oss_unsc[, grep("rain.", names(df_oss_unsc))]
      )
    )

    
    # original model object code:
    # UMF.oss.3 <- unmarkedFrameOccu(
    #   y = oss.dets,
    #   obsCovs = list(soilmoist = scaled_soilmoist),
    #   siteCovs = scaled_sitecovs)
    # 
    # 
    
## OSS dwd cov exploration -------------------------------------------------------------------------------
    
    # Run several models to determine which dwd covariate is most important for occupancy
    # formula: occu(~det covs ~occu covs, data=, se=TRUE)

    #1: psi(.) p(.)
    m1 <- occu(~1 ~1, data=UMF.oss.1)
    
    #2: psi(dwd_cov) p(.)
    m2 <- occu( ~dwd_cov ~trt +dwd_cov , data=UMF.oss.1)
   
    #3: psi(dwd_count) p(.)
    m3 <- occu(~dwd_count ~trt+dwd_count, data=UMF.oss.1)
    
    #4: psi(avg_volume) p(.)
    m4 <- occu(~avg_volume ~trt+avg_volume, data=UMF.oss.1)
    
    #scale in the model with unmarked
    
    # dwd_count is the best option here, with the lowest AIC and p-value
    # all dwd variables are showing up as negatively correlated with occupancy, which is really troubling
    # using the unscaled data gives me NA's and warnings for the dwd_count variable
    
    
## OSS det cov exploration ----------------------------------------------------------------------------------
    
    # Run several models to determine which det covariate is most important
    
    #5: psi(dwd_count) p(soimoist)
    m5 <- occu(~soilmoist ~dwd_count, data=UMF.oss.2)

    #6: psi(dwd_count) p(temp)
    m6 <- occu( ~temp ~dwd_count, data=UMF.oss.2)

    #7: psi(dwd_count) p(rain)
    m7 <- occu(~rain ~dwd_count, data=UMF.oss.2)

    
    # none of the models with det covs have higher AIC or significant variables
    # according to my random forest, time since rain is pretty important. this is confusing
    
    
     
## OSS top model exploration -------------------------------------------------------------------------------
    

    









