# -------------------------------------------------------------------------------------------------------
##
## 04-nimble-results-figures.R 
##
## Figures created from the nimble occupancy models
##
## Jasmine Williamson
## Date Created: 02-13-2025
##
## 

## settings -----------------------------------------------------------------------------------------------

    rm(list=ls())
    setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/oss-occu/data")
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    library(boot)
    source('attach.nimble_v2.R')


## load and name data  ------------------------------------------------------------------------------------

# saved model outputs
    # all species model:
    # load("./all.spp_model_2.RData")
    
    # oss model:
    load("./oss_model.RData")
    attach.nimble(mcmc.output.3$samples)
    oss.occu.prob <- inv.logit(TreatmentIntercept) # Inverse logit the trt intercept to get occu probabilities
    oss.det.prob <- inv.logit(DetectionIntercept) # Inverse logit the detection intercept to get det probabilities
    
    # enes model:
    load("./enes_model.RData")
    attach.nimble(mcmc.output.4$samples)
    enes.occu.prob <- inv.logit(TreatmentIntercept)
    enes.det.prob <- inv.logit(DetectionIntercept)
    
# detections data
    oss.dets <- read.csv("occupancy/oss.occu.wide.csv")
    enes.dets <- read.csv("occupancy/enes.occu.wide.csv")
    
# site data
    site <- read.csv("site.complete.csv")
    site_subset <- subset(site[,c(1,5)])

## naive occupancy ------------------------------------------------------------------------------------------
    
# Compute naive occupancy: proportion of sites with at least one detection
    oss_naive_occu <- mean(rowSums(oss.dets) > 0)
    
        # 0.6062992
    
    enes_naive_occu <- mean(rowSums(enes.dets) > 0)
    
        # 0.4409449
    
    # by treatment
        #  1 = BS     
        #  2 = BU     
        #  3 = HB     
        #  4 = HU     
        #  5 = UU
    
# oss
    oss_site_occupancy <- rowSums(oss.dets) > 0 # Identify which sites have at least one detection
    
    # add occu to site and trt
    site_subset$oss_occu <- oss_site_occupancy
    
    oss_naive_occ_by_trt <- aggregate(oss_occu ~ trt, data = site_subset, mean)
    
    # treatment  oss_occu
    # 1         1 0.5000000
    # 2         2 0.7307692
    # 3         3 0.6000000
    # 4         4 0.4444444
    # 5         5 0.7600000
    
# enes
    enes_site_occupancy <- rowSums(enes.dets) > 0 # Identify which sites have at least one detection
    
    # add occu to site and trt
    site_subset$enes_occu <- enes_site_occupancy
    
    enes_naive_occ_by_trt <- aggregate(enes_occu ~ trt, data = site_subset, mean)
    
    # treatment enes_occu
    # 1         1 0.2083333
    # 2         2 0.6538462
    # 3         3 0.2400000
    # 4         4 0.2222222
    # 5         5 0.8800000
    
## interpreting model outputs----------------------------------------------------------------------------------------------------
    
    #this code is not that useful but i kept it in here in case I need it later

    hist(oss.det.prob)
    
    mean(oss.det.prob) # = 0.26108
    mean(oss.det.prob>0)  # = 1
    median(oss.det.prob)  # = 0.2597028
    boxplot(oss.det.prob)
    
    # Median Detection Estimates
    median(enes.det.prob[,1]) # 0.3240405   BS 
    median(enes.det.prob[,2]) # 0.2747833    BU
    median(enes.det.prob[,3]) # 0.2573965    HB
    median(enes.det.prob[,4]) # 0.2087794    HU
    median(enes.det.prob[,5]) # 0.2392232    UU
    
    # occu estimates before adding separate detection intercepts
    median(oss.occu.prob[,1]) # 0.5962445    BS 
    median(oss.occu.prob[,2]) # 0.7885033    BU
    median(oss.occu.prob[,3]) # 0.5864983    HB
    median(oss.occu.prob[,4]) # 0.8241218    HU
    median(oss.occu.prob[,5]) # 0.9940393    UU

    # to calculate confidence intervals from these intercepts
    CI_BS <- quantile(oss.occu.prob[,1], probs = c(0.025, 0.975))  
    CI_BU <- quantile(oss.occu.prob[,2], probs = c(0.025, 0.975))  
    CI_HB <- quantile(oss.occu.prob[,3], probs = c(0.025, 0.975))   
    CI_HU <- quantile(oss.occu.prob[,4], probs = c(0.025, 0.975))   
    CI_UU <- quantile(oss.occu.prob[,5], probs = c(0.025, 0.975))  

    
    
    
## base R boxplot for detection probability -----------------------------------------------------------------------------------------------
    
    new.names <- c("Salvage Logged", "Wildfire", "Harvest, Wildfire", "Harvest", "Control")

# oss
    oss.det.prob <- inv.logit(DetectionIntercept)
    
    colnames(oss.det.prob) <- new.names
    desired.order <- c("Control", "Wildfire", "Harvest, Wildfire", "Harvest", "Salvage Logged")
    
    box.colors <- c('lightgreen','steelblue', 'coral2', '#f9d62e', '#b967ff' )
    
    #boxplot of Treatment Estimates - occu prob for each treatment
    boxplot(oss.det.prob[, match(desired.order, colnames(oss.det.prob))], 
            main = "Treatment Intercepts for OSS", 
            xlab = "Treatment", ylab = "Detection Probability",
            col = box.colors)
    
    
# enes
    enes.det.prob <- inv.logit(DetectionIntercept)
    
    colnames(enes.det.prob) <- new.names
    desired.order <- c("Control", "Wildfire", "Harvest, Wildfire", "Harvest", "Salvage Logged")
    
    box.colors <- c('lightgreen','steelblue', 'coral2', '#f9d62e', '#b967ff' )
    
    #boxplot of Treatment Estimates - occu prob for each treatment
    boxplot(enes.det.prob[, match(desired.order, colnames(enes.det.prob))], 
            main = "Detection Intercepts for ENES", 
            xlab = "Treatment", ylab = "Detection Probability",
            col = box.colors)
    
    
## oss ggplot boxplot treatment effect -----------------------------------------------------------------------------------------------
    

    new.names <- c("Salvage Logged", "Wildfire", "Harvest, Wildfire", "Harvest", "Control")
    colnames(oss.occu.prob) <- new.names
    oss.occu.prob <- as.data.frame(oss.occu.prob)
    
    oss.long <- gather(oss.occu.prob, key = "Treatment", value = "Occupancy_Probability") #reshape data to long format bc ggplot required
    
    oss.long$Treatment <- factor(oss.long$Treatment,#specify treatment order for plot
                                 levels = c("Control", "Wildfire", "Harvest, Wildfire", "Harvest", "Salvage Logged")) # Replace "DesiredOrder", "Treat1", etc., with your actual treatment names in the desired order
    
    
    box.colors <- c('lightgreen','steelblue', 'coral2', '#f9d62e', '#b967ff' )
    
    p <- ggplot(oss.long, aes(x = Treatment, y = Occupancy_Probability, fill = Treatment)) +
      geom_boxplot() +
      theme_classic() +
      scale_fill_manual(values = box.colors) +
      labs(title = "Treatment Intercepts for OSS", x = "Treatment", y = "Occupancy Probability") +
      theme(plot.title = element_text(hjust = 0.5)) # Center the title
    p
    
    ggsave(filename = "oss_trt_occu_prob_0204.png", plot = p, device = "png", 
           path = "~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/oss-occu/figures/06-rework-nimble-models",
           width = 10, height = 6, units = "in", dpi = 300)
    
    
## enes ggplot boxplot treatment effect -----------------------------------------------------------------------------------------------

    new.names <- c("Salvage Logged", "Wildfire", "Harvest, Wildfire", "Harvest", "Control")
    colnames(enes.occu.prob) <- new.names
    enes.occu.prob <- as.data.frame(enes.occu.prob)
    
    enes.long <- gather(enes.occu.prob, key = "Treatment", value = "Occupancy_Probability") #reshape data to long format bc ggplot required
    
    #specify treatment order for plot
    enes.long$Treatment <- factor(enes.long$Treatment, 
                                  levels = c("Control", "Wildfire", "Harvest, Wildfire", "Harvest", "Salvage Logged")) # Replace "DesiredOrder", "Treat1", etc., with your actual treatment names in the desired order
    
    
    box.colors <- c('lightgreen','steelblue', 'coral2', '#f9d62e', '#b967ff' )
    
    p2 <- ggplot(enes.long, aes(x = Treatment, y = Occupancy_Probability, fill = Treatment)) +
      geom_boxplot() +
      theme_classic() +
      scale_fill_manual(values = box.colors) +
      labs(title = "Treatment Intercepts for ENES", x = "Treatment", y = "Occupancy Probability") +
      theme(plot.title = element_text(hjust = 0.5)) # Center the title
    p2
    
    ggsave(filename = "enes_trt_occu_prob.png", plot = p2, device = "png", 
           path = "C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/oss-occu/figures/06-rework-nimble-models",
           width = 10, height = 6, units = "in", dpi = 300)
    
    
## ggplot with both species -----------------------------------------------------------------------------------------------

    # Rename treatments
    new.names <- c("Salvage", "Burn", "Harvest,Burn", "Harvest", "Control")
    
    # Convert OSS occupancy probabilities to a data frame and reshape
    colnames(oss.occu.prob) <- new.names
    oss.occu.long <- gather(as.data.frame(oss.occu.prob), key = "Treatment", value = "Occupancy_Probability")
    oss.occu.long$Species <- "OSS"  # Add species identifier
    
    # Convert ENES occupancy probabilities to a data frame and reshape
    colnames(enes.occu.prob) <- new.names
    enes.occu.long <- gather(as.data.frame(enes.occu.prob), key = "Treatment", value = "Occupancy_Probability")
    enes.occu.long$Species <- "ENES"  # Add species identifier
    
    # Combine both datasets
    combined_data <- bind_rows(oss.occu.long, enes.occu.long)
    
    # Ensure Treatment order is correct
    combined_data$Treatment <- factor(combined_data$Treatment, 
                                      levels = c("Control", "Burn", "Harvest,Burn", "Harvest", "Salvage"))
    
    # Define colors
    box.colors <- c('lightgreen','steelblue', 'coral2', '#f9d62e', '#b967ff')
 
       
# plot
    
    p3 <- ggplot(combined_data, aes(x = Treatment, y = Occupancy_Probability, fill = Treatment)) +
      geom_boxplot() +
      facet_wrap(~Species, scales = "free_y") +  # Separate plots for OSS and ENES
      theme_classic() +
      scale_fill_manual(values = box.colors) +
      labs(title = "Occupancy Probability by Treatment for OSS and ENES", 
           x = "Treatment", y = "Occupancy Probability") +
      theme(plot.title = element_text(hjust = 0.5))  # Center title
    

ggsave(filename = "faceted_spp_occu_prob.png", plot = p3, device = "png", 
           path = "~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/oss-occu/figures/03-nimble-models",
           width = 12, height = 6, units = "in", dpi = 300)
    
    
    
    p4 <- ggplot(combined_data, aes(x = Treatment, y = Occupancy_Probability, fill = Species)) +
      geom_boxplot(position = position_dodge(0.8)) +  # Dodge for side-by-side boxes
      theme_classic() +
      scale_fill_manual(values = c("OSS" = "steelblue", "ENES" = "coral2")) +
      labs(title = "Occupancy Probability by Treatment for OSS and ENES", 
           x = "Treatment", y = "Occupancy Probability") +
      theme(plot.title = element_text(hjust = 0.5))  # Center title

    
ggsave(filename = "both_occu_prob.png", plot = p3, device = "png", 
           path = "~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/oss-occu/figures/06-rework-nimble-models",
           width = 10, height = 6, units = "in", dpi = 300)
    
    
plot(combined_data)
