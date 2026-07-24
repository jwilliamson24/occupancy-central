## -------------------------------------------------------------------------------------------------------
##
## 02-unmarked-models-trt.R 
##
## Rerun the occupancy models using package Unmarked from 2023 analysis
##
## Jasmine Williamson
## Date Created: 07-09-2024
##
## 


# where i left off: 

# 7/4/24
# i figured out how to format the detection covariates from long format into wide format.
# the code below works for soil moisture and weather
# i also created a site covs data frame that includes many of the covs i will want
# 
# what i still need to do:
# i want to create another detection df for downed wood count per subplot
# i want to add more covs to the site data frame

# 7/9/24
# created occu dfs for oss and enes
# moved them from long to wide format

# 7/9/24
# made unmarked models that run
# but not sure about the results
# cant use temp and hum as det covs because i only took them at the sites and not subplots
# unsure of what to do from here
# nothing showing up as significant



## settings -----------------------------------------------------------------------------------------------

    rm(list=ls())
    setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/oss-occu/data")
    
    library(unmarked)
    library(ggplot2)
    library(stats)
    library(MASS)
    library(tidyverse)



## load data----------------------------------------------------------------------------------------------

    site <- read.csv("site.complete.csv")
    subplot <- read.csv("subplot.complete.csv")
    sals <- read.csv("sals.complete.csv", 
                     colClasses = c(landowner="factor", stand="character", trt="factor",
                                    obs="factor", subplot="factor", recap="factor",
                                    pass="factor", spp="factor", cover_obj="factor", 
                                    substrate="factor", age_class="factor"))
    
    oss.dets <- read.csv("occupancy/oss.occu.wide.csv")
    enes.dets <- read.csv("occupancy/enes.occu.wide.csv")


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
    
    soil_moist_wide <- df_wide[,-1]
    soil_moist_wide <- as.data.frame(soil_moist_wide)
    rownames(soil_moist_wide) <- soil_moist_wide[,1]
    soil_moist <- soil_moist_wide[,-1]


# weather 

    #subset df with only cols i need (cant use site id here bc its a character)
    df <- subset(subplot[,c('subplot','weather')])
    
    #change to long format
    df_long <- df %>%
      pivot_longer(cols = -subplot, names_to = "variable", values_to = "value")
    #add site id back in
    df_long <- cbind(df_long, subplot$site_id)
    
    #reshape the data back to wide format, flipping the dimensions
    df_wide <- df_long %>%
      pivot_wider(names_from = subplot, values_from = value)
    
    weather_wide <- df_wide[,-1]
    weather_wide <- as.data.frame(weather_wide)
    rownames(weather_wide) <- weather_wide[,1]
    weather <- weather_wide[,-1]



## Site/occu covs ----------------------------------------------------------------------------------

    sitecovs <- subset(site[,c('site_id','trt','date_mdy','elev','temp','hum')])
    sitecovs$trt <- factor(sitecovs$trt, 
                           levels = c("UU", "BU", "HB", "HU", "BS"))
    
    # Scale Occupancy Covariates
    sitecovs$site_id <- as.character(sitecovs$site)
    # Identify numeric sitecovs
    numeric_sitecovs <- sitecovs[sapply(sitecovs, is.numeric)]
    # Scale the numeric sitecovs
    scaled_numeric_sitecovs <- as.data.frame(scale(numeric_sitecovs))
    # Combine the scaled numeric sitecovs with non-numeric sitecovs
    scaled_sitecovs <- cbind(scaled_numeric_sitecovs, sitecovs[setdiff(names(sitecovs), names(numeric_sitecovs))])
    
    #write.csv(scaled_sitecovs, "C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/oss-occu/data/scaled.sitecovs.csv",
          #    row.names = FALSE)


## name and scale covs for unmarked object ----------------------------------------------------------

#data frames with site rows and subplots columns
    hum <- sitecovs$hum
    temp <- sitecovs$temp
    weather <- weather
    date <- site$date_mdy
    soilmoist <- soil_moist
    date <-  site$date_mdy
    
    #scale det dovs
    scaled_hum <- as.data.frame(scale(hum))
    scaled_temp <- as.data.frame(scale(temp))
    #scaled_date <- as.data.frame(scale(date))
    scaled_soilmoist <- as.data.frame(scale(soilmoist))



## OSS Top Model Exploration-------------------------------------------------------------------------------

# Build OSS unmarkedFrameOccu Object
    UMF.oss.scaled <- unmarkedFrameOccu(
      y = oss.dets, 
      obsCovs = list(weather = weather, soilmoist = scaled_soilmoist), 
      siteCovs = scaled_sitecovs)

# Run several models to determine Psi and P covariates
    # occu(~det covs ~occu covs, data=, se=TRUE)

    #1: psi(.) p(.)
    m1 <- occu(~1 ~1, data=UMF.oss.scaled)
    #AIC: 819.7476
    
    # 2: psi(treatment) p(.)
    m2 <- occu(~1 ~trt, data=UMF.oss.scaled)
    #AIC: 819.3336 
    
    # 3: psi(.) p(x) , x = date, treatment, or soil moisture                
    m3 <- occu(~soilmoist ~1, data=UMF.oss.scaled)
    #AIC: 820.4457 
    
    # 4: psi(treatment) p(x) , x = date, treatment, or soil moisture
    m4 <- occu(~temp ~trt, data=UMF.oss.scaled)
    #AIC: 821.3062 

    


## Extract treatment occupancy predictions and plot ------------------------------------------------------------

# Extract psi predictions from model
    preds <- predict(m2, type="state")
    sites_trt <- sitecovs[, c("site_id","trt")]
    site_preds <- cbind(sites_trt, preds)
    
    
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


