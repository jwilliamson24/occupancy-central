## =================================================
##
## Title: vif_correlation_matrices
## Author: Jasmine Williamson
## Date Created: 07/01/2025
##
## Description: Checking for collinearity for occupancy models
##
## VIF is below 5 for all tested covs, so i should be good to include any of them.
##
## =================================================

## load packages
  library(corrplot)
  library(car)
  library(dplyr)

## load data -------------------------------------------------------------------

  env <- read.csv("data/covariate matrices/env_subset_corr2.csv", row.names = 1) # already subsetted from earlier correlations script
  env$site_id <- rownames(env)
  
  #site-level data
  dat <- read.csv("data/covariate matrices/site_aspect_precip_all_vars.csv", row.names = 1) 
  env2 <- dat[c(1:25,28:31)]
  
  #extra dwd metrics
  dwd <- read.csv("data/covariate matrices/dwd.extra.metrics.csv")
  
  #add dwd dens and avg volume from dwd to env df
  env_merge <- cbind(env2, dwd[,c("dwd_dens","avg_volume")])
  
  
  
## Pearson Correlation matrix -----------------------------------------------------
  
  # did this the long way in earlier correlations script,
  # just running these two to confirm my choices for covariate removal 

  #par(mfrow = c(1, 1))  
  
  # all vars
  num_vars <- env_merge[, sapply(env_merge, is.numeric)] # numeric only
  cor_matrix <- cor(num_vars, use = "complete.obs", method = "pearson") # (use = "complete.obs" handles NAs)
  corrplot(cor_matrix, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
  
  # correlated pairs 
  threshold <- 0.6
  cor_pairs <- which(abs(cor_matrix) > threshold & abs(cor_matrix) < 1, arr.ind = TRUE)
  cor_pairs <- cor_pairs[cor_pairs[, 1] < cor_pairs[, 2], ]
  
  high_corrs <- data.frame(
    Var1 = rownames(cor_matrix)[cor_pairs[, 1]],
    Var2 = colnames(cor_matrix)[cor_pairs[, 2]],
    Correlation = cor_matrix[cor_pairs]
  )
  
  
  
  # removed correlated vars
  num_vars2 <- env[sapply(env, is.numeric)] # numeric only
  cor_matrix2 <- cor(num_vars2, use = "complete.obs", method = "pearson")
  corrplot(cor_matrix2, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
 
  # correlated pairs 
  threshold <- 0.5
  cor_pairs2 <- which(abs(cor_matrix2) > threshold & abs(cor_matrix2) < 1, arr.ind = TRUE)
  cor_pairs2 <- cor_pairs2[cor_pairs2[, 1] < cor_pairs2[, 2], ]
  
  low_corrs <- data.frame(
    Var1 = rownames(cor_matrix2)[cor_pairs2[, 1]],
    Var2 = colnames(cor_matrix2)[cor_pairs2[, 2]],
    Correlation = cor_matrix2[cor_pairs2]
  )
  
  
  
## Variance Inflation Factor -----------------------------------------------------
  
  # write linear model with all variables of interest

  # data for lm
  env_lm <- merge(env, dat[, c("site_id","trt", "oss", "landowner")], by = "site_id", all.x = TRUE)
  env_lm$oss[env_lm$oss > 1] <- 1 # oss from counts to 0/1
  
  # landowner to mgmt type                   
  env_lm <- env_lm %>%
    mutate(mgmt_type = case_when(
      landowner == "PB" ~ "0", # private
      landowner == "WY" ~ "0",  # private
      landowner == "BLM" ~ "1",  # public
      landowner == "ODF" ~ "1",  # public
      TRUE ~ "default_value"  
    ))
  

  
  glm_vif <- glm(oss ~ trt + lat + long + elev + temp + canopy_cov + veg_cov + fwd_cov +
                   soil_moist + dwd_count + decay_cl + precip_mm + days_since_rain + mgmt_type,
                 data = env_lm, family = "binomial")
  
  
  # Calculate VIFs
  # use GVIF^(1/(2*Df)) 
  vif(glm_vif)
  
    #                      GVIF  Df    GVIF^(1/(2*Df))
    # trt             66.225707  4        1.688995
    # lat              4.118525  1        2.029415
    # long             3.906639  1        1.976522
    # elev             2.448984  1        1.564923
    # temp             2.029414  1        1.424575
    # canopy_cov      11.886901  1        3.447739
    # veg_cov          1.773410  1        1.331694
    # fwd_cov          2.104887  1        1.450823
    # soil_moist       2.514022  1        1.585567
    # dwd_count        1.939834  1        1.392779
    # decay_cl         1.876228  1        1.369755
    # precip_mm        1.742795  1        1.320150
    # days_since_rain  2.543149  1        1.594725
    # mgmt_type        5.105111  1        2.259449
  
  
  ## They're all below 5, so i should be good to include any of them.
  
  