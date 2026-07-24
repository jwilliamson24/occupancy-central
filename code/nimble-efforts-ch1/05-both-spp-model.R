# -------------------------------------------------------------------------------------------------------
##
## 05-both-spp-model.R 
##
## Task: Create model with both spp together, separate occu prob, with same det prob
##
## Jasmine Williamson
## Date Created: 03/07/2025
##
## 

## settings -----------------------------------------------------------------------------------------------

    rm(list=ls())
    library(tidyr)
    library(nimble)
    library(ggplot2)
    library(data.table)
    library(tidyverse)
    library(mcmcplots)
    library(MCMCvis)
    library(boot)
    source('attach.nimble_v2.R')

## load Data - Github version for Josh ----------------------------------------------------------------------------

# # model output
#     load("./all-model.RData")
# 
#     source('https://raw.githubusercontent.com/jwilliamson24/oss-occu/refs/heads/main/data/attach.nimble_v2.R')
#     attach.nimble(mcmc.output.1$samples)
# 
# # input data
#     site <- read.csv("https://raw.githubusercontent.com/jwilliamson24/oss-occu/refs/heads/main/data/site.complete.csv")
#     subplot <- read.csv("https://raw.githubusercontent.com/jwilliamson24/oss-occu/refs/heads/main/data/subplot.complete.csv")
#     
#     oss.long <- read.csv("https://raw.githubusercontent.com/jwilliamson24/oss-occu/refs/heads/main/data/occupancy/oss.occu.long.csv")
#     enes.long <- read.csv("https://raw.githubusercontent.com/jwilliamson24/oss-occu/refs/heads/main/data/occupancy/enes.occu.long.csv")
#     
#     env_subset_corr <- read.csv("https://raw.githubusercontent.com/jwilliamson24/oss-occu/refs/heads/main/data/covariate%20matrices/env_subset_corr.csv")    # df with precip data
    

## load Data--------------------------------------------------------------------------------------------------

    site <- read.csv("site.complete.csv")
    subplot <- read.csv("subplot.complete.csv")
    sals <- read.csv("sals.complete.csv",
                     colClasses = c(landowner="factor", stand="character", trt="factor",
                                    obs="factor", subplot="factor", recap="factor",
                                    pass="factor", spp="factor", cover_obj="factor",
                                    substrate="factor", age_class="factor"))

    all.long <- read.csv("occupancy/all.occu.long.csv")
    oss.long <- read.csv("occupancy/oss.occu.long.csv")
    enes.long <- read.csv("occupancy/enes.occu.long.csv")

    env_subset_corr <- read.csv("covariate matrices/env_subset_corr.csv")    # df with precip data

    
## format data --------------------------------------------------------------------------------------------------
    
    # we need two parts: one subplot-level matrix, and one site-level matrix
    
    
    # subplot-level matrix:
    # "date" "site" "temp" "humidity" "soil.moist" "all.obs" "oss.obs" "enes.obs" 
    
    # subset and merge the various dataframes
    df1 <- subset(site[,c(1,10,11)]) #subset site vars
    df2 <- subset(subplot[,c(8,1,9,18)]) #subset subplot vars
    df3 <- subset(env_subset_corr[,c(1,14)]) #subset precip: site and days since rain
    colnames(df3)[1] <- "site_id"
    
    #merge site subset with precip subset
    df4 <- full_join(df1,df3,by="site_id")
    
    #merge subplot df and new site df
    dfmerge <- full_join(df2,df4,by="site_id") 
    #joins site vars to subplot matrix, which automatically repeats the values in blank 
    #subplot cells from same site id for both years bc its currently only listed for each site
    
    
    #rename sal columns and rbind
    colnames(oss.long)[3] <- "det"
    oss.long$spp <- "oss"
    colnames(enes.long)[3] <- "det"
    enes.long$spp <- "enes"
    all.spp.long <- rbind(oss.long, enes.long)
    
    colnames(oss.long)[3] <- "oss"
    colnames(enes.long)[3] <- "enes"
    sals <- merge(oss.long, enes.long, by=c("site_id","subplot")) 
    
    det <- sals[, c("oss","enes")]
    # det <- as.data.frame(nrow = nrow(sals), ncol = 2, byrow = TRUE)
    # det[,1] <- sals$oss  # First column for oss
    # det[,2] <- sals$enes # Second column for enes
    
    # 1 = oss
    # 2 = enes
    
    # merge sal columns with env variables
    #merge1 <- merge(dfmerge, all.spp.long, by=c("site_id","subplot"))

    merge4 <- cbind(dfmerge, det)
    
    # site-level matrix: "site" and "treatment"
    df5 <- subset(site[,c(1,5)])
    df5$trt <- as.numeric(as.factor(df5$trt))
    table(df5$trt)
    colnames(df5)[2] <- "treatment"
    
    #  1 = BS     
    #  2 = BU     
    #  3 = HB     
    #  4 = HU     
    #  5 = UU   
    
    merge5 <- merge(merge4, df5, by="site_id")
    colnames(merge5) <- c("site", "date","subplot", "soil.moist","temp",
                          "humidity","days.since.rain","oss","enes","trt")
    
## name data for model --------------------------------------------------------------------------------------------------
    
    data <- merge5
    #data$det <- as.numeric(data$det)
    scaled_temp <- c(scale(data$temp))
    
    #data2 <- df5
    
    
## MODEL  ------------------------------------------------------------------------------------------------
    
# model that includes both species
# enes and oss will be allowed different occu probs
# det prob will be the same for both spp
    
    
    all.model <- nimbleCode ({
      
      # Priors - uninformative/vague 
      
      for(t in 1:n.treatments){
        TreatmentIntercept[t,1] ~ dunif(-10,10) # spp 1 occu int
        TreatmentIntercept[t,2] ~ dunif(-10,10) # spp 2 occu int
        DetectionIntercept[t] ~ dunif(-5,5) # shared app detection intercept (per trt)
      }#t
      
      # env priors
      betaTemp ~ dnorm(-0, 2)
      betaTemp2 ~ dnorm(-0, 2)
      betaRain ~ dnorm(-0, 2)
      
      # Likelihood
      
      # Process/Biological model = Occupancy
      # separate occu probs for each spp
      for(i in 1:n.sites) {
        logit(psi[i,1]) <- TreatmentIntercept[treatment[i],1]
        logit(psi[i,2]) <- TreatmentIntercept[treatment[i],2]
        z[i,1] ~ dbern(psi[i,1])  
        z[i,2] ~ dbern(psi[i,2])
      }#i
      # wrap in sub forloops
      
      # Observation model = Detection
      # need two pieces: one for p(det prob coeff) and one defining Y distribution
      for(j in 1:n.obs) { 
        logit(p[j]) <-  DetectionIntercept[treatment[site[j]]] + #p=detprob for site i, survey j
          #DetectionIntercept +
          betaTemp*temp[j] + 
          betaTemp2*temp[j]^2 +
          betaRain*days.since.rain[j]
        Y[j,1] ~ dbern(p[j] * z[site[j],1]) #Y=my data observations #wrap in subloop
        Y[j,2] ~ dbern(p[j] * z[site[j],2])
      }#j
      
    })
    #data input change - matrix with spp1,2
    #shared det int and relship with temp
    
    #
    
    # Parameters monitored
    parameters <- c("z", "p", "TreatmentIntercept", 
                    "DetectionIntercept","betaTemp", "betaTemp2", "betaRain")
    
    # MCMC Settings
    ni <- 80000
    nt <- 80
    nb <- 40000
    nc <- 3
    
    
    # Data
    nimble.data = list(Y=det,
                       temp=scaled_temp,
                       days.since.rain=data$days.since.rain)
    
    nimble.constants = list(n.sites = length(unique(data$site)),
                            n.treatments = length(unique(data$trt)),
                            treatment=data$trt,
                            site=as.numeric(as.factor(data$site)),
                            n.obs = nrow(data))
    
    mcmc.output.1 <- nimbleMCMC(code = all.model,
                                data = nimble.data,
                                constants=nimble.constants,
                                monitors = parameters,
                                niter = ni,
                                nburnin = nb,
                                nchains = nc,
                                thin=nt,
                                summary=TRUE,
                                samplesAsCodaMCMC = TRUE)
    

    
    attach.nimble(mcmc.output.1$samples)
    save(mcmc.output.1, file="./all-model.RData")
    # load("./enes_model.RData")
    
    
    
    
## assessing convergence-------------------------------------------------------------------------------------

# Gelman-Rubin diagnostic 
    z <- mcmc.output.1$samples
    g <- matrix(NA, nrow=nvar(z), ncol=2)
    for (v in 1:nvar(z)) { g[v,] <- gelman.diag(z[,v])$psrf }
    PSRF <- bind_cols(colnames(z$chain1),g) %>% rename(Parameter = ...1 ,PSRF = ...2 ,PSRFUpperCI = ...3) 
    
    PSRF 
    # most values are < 1.1 = good
    
    
# Trace plots and parameter estimates
    MCMCtrace(object = mcmc.output.1$samples,
              pdf = FALSE, # no export to PDF
              ind = TRUE, # separate density lines per chain
              params = c("DetectionIntercept", "betaTemp", "TreatmentIntercept"))
    
    plot(mcmc.output.1$samples[, 1:5]) 
    par(mfrow=c(1,1))
    
    # these look okay, some of the parameter estimates are wobbly/ not entirely overlapping
    
    
## plots, parameter summaries    -------------------------------------------------------------------------------------
   
    summary(mcmc.output.1$samples)
    
    hist(TreatmentIntercept)
    hist(DetectionIntercept)
    
    det.probs.inv <- inv.logit(DetectionIntercept)
    boxplot(det.probs.inv)
    
    trt.int.inv <- inv.logit(TreatmentIntercept)
    boxplot(trt.int.inv)    
    
    
    
    
# Treatment Effect
    
    get_treatment_intercepts <- function(chain) {
      # Get only the TreatmentIntercept parameters
      treatment_params <- grep("TreatmentIntercept", colnames(chain))
      as.matrix(chain[, treatment_params])
    }
    
    # Extract intercepts from each chain
    chain1_intercepts <- get_treatment_intercepts(mcmc.output.1$samples[[1]])
    chain2_intercepts <- get_treatment_intercepts(mcmc.output.1$samples[[2]])
    chain3_intercepts <- get_treatment_intercepts(mcmc.output.1$samples[[3]])
    
    # Calculate mean and SD for each treatment-species combination
    treatment_means <- apply(rbind(chain1_intercepts, chain2_intercepts, chain3_intercepts), 2, mean)
    treatment_sds <- apply(rbind(chain1_intercepts, chain2_intercepts, chain3_intercepts), 2, sd)
    
    
    treatment_names <- c("BS", "BU", "HB", "HU", "UU")
    
    # Prepare data for plotting
    plot_data <- data.frame(
      Treatment = rep(treatment_names, 2),
      Species = rep(c("Spp1 OSS", "Spp2 ENES"), each = 5),
      Mean = c(treatment_means[1:5], treatment_means[6:10]),
      SD = c(treatment_sds[1:5], treatment_sds[6:10])
    )
    plot_data$Treatment <- factor(plot_data$Treatment, levels = c("UU","BU","HB","HU","BS"))
    
    
    # Create the plot
    ggplot(plot_data, aes(x = factor(Treatment), y = Mean, color = Species)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = Mean - 1.96*SD, ymax = Mean + 1.96*SD), width = 0.2) +
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(title = "Treatment Intercepts by Species",
           x = "Treatment",
           y = "Intercept Value") +
      theme(axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14))

    
    
    
    
# Occupancy probabilities

    
    calculate_occupancy <- function(intercept) {
      # Convert log-odds to probs
      exp(intercept) / (1 + exp(intercept))
    }

    # Calculate mean occupancy prob for each treatment-species combination
    occupancy_means <- apply(rbind(chain1_intercepts, chain2_intercepts, chain3_intercepts),
                             2, function(x) mean(calculate_occupancy(x)))
    
    # Create df for occupancy probs
    occupancy_data <- data.frame(
      Treatment = rep(treatment_names, 2),
      Species = rep(c("Spp1 OSS", "Spp2 ENES"), each = 5),
      Occupancy = c(occupancy_means[1:5], occupancy_means[6:10])
    )
    
    occupancy_data$Treatment <- factor(occupancy_data$Treatment, levels = c("UU","BU","HB","HU","BS"))
    
    # Create occupancy prob plot
    ggplot(occupancy_data, aes(x = factor(Treatment), y = Occupancy, color = Species)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = Occupancy - 0.1, ymax = Occupancy + 0.1), width = 0.2) +
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(title = "Occupancy Probabilities by Treatment and Species",
           x = "Treatment",
           y = "Occupancy Probability") +
      theme(axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14)) +
      coord_cartesian(ylim = c(0, 1))  

    # treatments
    #  1 = BS     
    #  2 = BU     
    #  3 = HB     
    #  4 = HU     
    #  5 = UU  
    
    # species
    # 1 = oss
    # 2 = enes
    
   