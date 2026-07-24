# -------------------------------------------------------------------------------------------------------
##
## 03-nimble-models.R 
##
## Task: Rework the nimble occupancy models from 2023 bayesian class
##
## Jasmine Williamson
## Date Created: 07-10-2024
##
## 

## settings -----------------------------------------------------------------------------------------------

  rm(list=ls())
  #setwd("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/oss-occu/data")
  setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/oss-occu/data")
  library(tidyr)
  library(nimble)
  library(ggplot2)
  library(data.table)
  library(tidyverse)
  library(mcmcplots)
  library(MCMCvis)
  library(boot)
  source('attach.nimble_v2.R')


## saved model outputs -----------------------------------------------------------------------------------

    # all species model:
    # load("./all.spp_model_2.RData")
    
    # oss model:
    # load("./oss_model.RData")
    # attach.nimble(mcmc.output.3$samples)
  
    # enes model:
    # load("./enes_model.RData")
    # attach.nimble(mcmc.output.4$samples)
  

## load Data--------------------------------------------------------------------------------------------------

    site <- read.csv("site.complete.csv")
    subplot <- read.csv("subplot.complete.csv")
    sals <- read.csv("sals.complete.csv", 
                     colClasses = c(landowner="factor", stand="character", trt="factor",
                                    obs="factor", subplot="factor", recap="factor",
                                    pass="factor", spp="factor", cover_obj="factor", 
                                    substrate="factor", age_class="factor"))
    
    oss.dets <- read.csv("oss.occu.wide.csv")
    enes.dets <- read.csv("enes.occu.wide.csv")
    
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
    
    #rename sal columns
    colnames(oss.long)[3] <- "oss.obs"
    colnames(enes.long)[3] <- "enes.obs"
    colnames(all.long)[3] <- "all.obs"
    
    #merge sal columns with env variables
    merge1 <- merge(dfmerge, all.long, by=c("site_id","subplot"))
    merge1 <- merge(merge1, oss.long, by=c("site_id","subplot"))
    complete.merge <- merge(merge1,enes.long,by=c("site_id","subplot"))
    
    colnames(complete.merge) <- c("site","subplot","date","soil.moist","temp", 
                                  "humidity","days.since.rain","all.obs","oss.obs","enes.obs")


# site-level matrix:
    # "site" and "treatment"
    df5 <- subset(site[,c(1,5)])
    df5$trt <- as.numeric(as.factor(df5$trt))
    table(df5$trt)
    colnames(df5)[2] <- "treatment"
    
    #  1 = BS     
    #  2 = BU     
    #  3 = HB     
    #  4 = HU     
    #  5 = UU   


## name data for model --------------------------------------------------------------------------------------------------
    
    data <- complete.merge
    
    data$all.obs <- as.numeric(data$all.obs)
    data$oss.obs <- as.numeric(data$oss.obs)
    data$enes.obs <- as.numeric(data$enes.obs)
    
    scaled_temp <- c(scale(data$temp))

    data2 <- df5
    
#how to treat days since rain as detection covariate? priors?    
    # ggplot(complete.merge, aes(x = days.since.rain, y = oss.obs)) +
    #   geom_point(color = "blue", size = 2) +  
    #   geom_smooth(method = "lm", color = "red", se = TRUE) 
    
    # days since rain and occu have slight negative relationship,
    # set prior as normal distribution with a negative mean (mu = -1) and a moderate stdev (sd = 1).
    


## MODEL 1 ------------------------------------------------------------------------------------------------
# separate det intercepts per treatment
    
    
    oss.model <- nimbleCode ({
      
      # Priors 
      # uninformative, vague priors
      
      for(t in 1:n.treatments){
        TreatmentIntercept[t] ~ dunif(-10,10)
        DetectionIntercept[t] ~ dunif(-5,5) # separate det int per treatment
      }#t
      
      #DetectionIntercept ~ dunif(-5,5) # single det int
      betaTemp ~ dnorm(-0, 2) 
      betaTemp2 ~ dnorm(-0, 2)
      betaRain ~ dnorm(-0, 2) 
      
      # Likelihood
      
      # Process/Biological model = Occupancy
      # need two pieces: one defining psi(occu prob coeff) and covariates, and one defining z dist
      for(i in 1:n.sites) {
        logit(psi[i]) <- TreatmentIntercept[treatment[i]]  #psi=occupancy probability
        z[i] ~ dbern(psi[i])  # z=1 if occupied, z=latent true occupancy
      }#i
      
      # Observation model = Detection
      # need two pieces: one for p(det prob coeff) and one defining Y distribution
      for(j in 1:n.obs) {
        logit(p[j]) <- DetectionIntercept[treatment[site[j]]] + 
                       #prob of detecting a sal in observation j depends on the trt assigned to the site that obs belongs to
                       #DetectionIntercept +
                       #DetectionIntercept[treatment[j]] +
                       betaTemp*temp[j] + 
                       betaTemp2*temp[j]^2 + #using temp as covariate with a quadratic relationship
                       betaRain*days.since.rain[j]
        
        #p=detection probability for site i and survey j
        Y[j] ~ dbern(p[j] * z[site[j]]) #Y=my actual data observations
        #z=1 or 0, turns this on or off
      }#j
      
    })
    
    # Parameters monitored
    parameters <- c("z","p","TreatmentIntercept","DetectionIntercept","betaTemp", "betaTemp2", "betaRain")
    
    # MCMC Settings
    ni <- 40000
    nt <- 40
    nb <- 20000
    nc <- 3
  
    
## For OSS Model 1 ---------------------------------------------------------------------------------
    
    # Data
    nimble.data = list(Y=data$oss.obs,
                       temp=scaled_temp,
                       days.since.rain=data$days.since.rain)
    
    nimble.constants = list(n.sites = length(unique(data$site)),
                            n.treatments = length(unique(data2$treatment)),
                            treatment=data2$treatment,
                            site=as.numeric(as.factor(data$site)),
                            n.obs = length(data$all.obs))
    
    mcmc.output.3 <- nimbleMCMC(code = oss.model,
                                data = nimble.data,
                                constants=nimble.constants,
                                monitors = parameters,
                                niter = ni,
                                nburnin = nb,
                                nchains = nc,
                                thin=nt,
                                summary=TRUE,
                                samplesAsCodaMCMC = TRUE)

    attach.nimble(mcmc.output.3$samples)
    save(mcmc.output.3, file="./oss_model.RData")
    # load("./oss_model.RData")

    
## For ENES Model 1 ---------------------------------------------------------------------------------

    # Data
    nimble.data = list(Y=data$enes.obs,
                       temp=scaled_temp,
                       days.since.rain=data$days.since.rain)
    
    nimble.constants = list(n.sites = length(unique(data$site)),
                            n.treatments = length(unique(data2$treatment)),
                            treatment=data2$treatment,
                            site=as.numeric(as.factor(data$site)),
                            n.obs = length(data$all.obs))
    
    mcmc.output.4 <- nimbleMCMC(code = oss.model,
                                data = nimble.data,
                                constants=nimble.constants,
                                monitors = parameters,
                                niter = ni,
                                nburnin = nb,
                                nchains = nc,
                                thin=nt,
                                summary=TRUE,
                                samplesAsCodaMCMC = TRUE)
    
    #######     warning: logProb of data node Y[724]: logProb is -Inf. what does this mean??
    
    
    attach.nimble(mcmc.output.4$samples)
    save(mcmc.output.4, file="./enes_model.RData")
    # load("./enes_model.RData")

    
## assessing convergence------------------------------------------------------------------------------------------------------

#mcmcplot(mcmc.output.1$samples)

# Gelman-Rubin diagnostic (AKA RHat or PSRF)
    z <- mcmc.output.3$samples
    g <- matrix(NA, nrow=nvar(z), ncol=2)
    for (v in 1:nvar(z)) { g[v,] <- gelman.diag(z[,v])$psrf }
    PSRF <- bind_cols(colnames(z$chain1),g) %>% rename(Parameter = ...1 ,PSRF = ...2 ,PSRFUpperCI = ...3) 
    
    PSRF # Values are mostly below 1.05 (good), but some are over
    
# Looking at trace plots and parameter estimates
    MCMCtrace(object = mcmc.output.3$samples,
              pdf = FALSE, # no export to PDF
              ind = TRUE, # separate density lines per chain
              params = c("DetectionIntercept", "betaTemp", "TreatmentIntercept"))
    #this looks good:
    #caterpillars on traceplot are mostly lined up/overlapping
    #parameter estimate lines on density plot are mostly overlapping
    
    #unsure about trtint 1 & 3 trace plots
    

## MODEL 2 ------------------------------------------------------------------------------------------------
# same det intercept for all treatments
    
    enes.model <- nimbleCode ({
  
      # Priors 
      # uninformative, vague priors
      
      for(t in 1:n.treatments){
        TreatmentIntercept[t,1] ~ dunif(-10,10) 
      }#t
      
      DetectionIntercept ~ dunif(-5,5) # single det intercept
      betaTemp ~ dnorm(-0, 2)
      betaTemp2 ~ dnorm(-0, 2)
      betaRain ~ dnorm(-0, 2)
      
      # Likelihood
      
      # Process/Biological model = Occupancy
      # need two pieces: one defining psi(occu prob coeff) and covariates, and one defining z dist
      for(i in 1:n.sites) {
        logit(psi[i,s]) <- TreatmentIntercept[treatment[i]]  #psi=occupancy probability
        z[i,s] ~ dbern(psi[i,s])  # z=1 if occupied, z=latent true occupancy
      }#i
      # wrap in sub forloops
      
      # Observation model = Detection
      # need two pieces: one for p(det prob coeff) and one defining Y distribution
      for(j in 1:n.obs) { 
        logit(p[j]) <-  DetectionIntercept +
                        betaTemp*temp[j] + 
                        betaTemp2*temp[j]^2 +
                        betaRain*days.since.rain[j]
        #using temp as covariate with a quadratic relationship
        #p=detection probability for site i and survey j
        Y[j,s] ~ dbern(p[j] * z[site[j] ,s]) #Y=my actual data observations #wrap in subloop 
        #z=1 or 0, turns this on or off
      }#j
      
    })
    #data input change - matrix with spp1,2
    #shared det int and relship with temp
    
    #
    
    # Parameters monitored
    parameters <- c("z","p","TreatmentIntercept","DetectionIntercept","betaTemp", "betaTemp2", "betaRain")
    
    # MCMC Settings
    ni <- 80000
    nt <- 80
    nb <- 40000
    nc <- 3

    
    # Data
    nimble.data = list(Y=data$enes.obs,
                       temp=scaled_temp,
                       days.since.rain=data$days.since.rain)
    
    nimble.constants = list(n.sites = length(unique(data$site)),
                            n.treatments = length(unique(data2$treatment)),
                            treatment=data2$treatment,
                            site=as.numeric(as.factor(data$site)),
                            n.obs = length(data$all.obs))
    
    mcmc.output.4 <- nimbleMCMC(code = enes.model,
                                data = nimble.data,
                                constants=nimble.constants,
                                monitors = parameters,
                                niter = ni,
                                nburnin = nb,
                                nchains = nc,
                                thin=nt,
                                summary=TRUE,
                                samplesAsCodaMCMC = TRUE)

#######     warning: logProb of data node Y[724]: logProb is -Inf. what does this mean??

    
    attach.nimble(mcmc.output.4$samples)
    save(mcmc.output.4, file="./enes_model.RData")
    # load("./enes_model.RData")


## assessing convergence------------------------------------------------------------------------------------------------------

# mcmcplot(mcmc.output.4$samples)

# Gelman-Rubin diagnostic (AKA RHat or PSRF)
    z <- mcmc.output.4$samples
    g <- matrix(NA, nrow=nvar(z), ncol=2)
    for (v in 1:nvar(z)) { g[v,] <- gelman.diag(z[,v])$psrf }
    PSRF <- bind_cols(colnames(z$chain1),g) %>% rename(Parameter = ...1 ,PSRF = ...2 ,PSRFUpperCI = ...3) 
    
    PSRF # some values were 1.05 and above, but i dont know what that means. 
    # fixed that after I changed iterations from 40000 to 80000
    
# ACF plots
    acf(DetectionIntercept)
    acf(TreatmentIntercept)
    #z
    
    
# Looking at trace plots and parameter estimates
    MCMCtrace(object = mcmc.output.4$samples,
              pdf = FALSE, # no export to PDF
              ind = TRUE, # separate density lines per chain
              params = c("DetectionIntercept", "betaTemp", "TreatmentIntercept"))
   
     #treatment 2 and 3 caterpillars dont look great
    
    
## josh ran these things after model to see outputs
    
    hist(betaRain)
    mean(betaRain<0)

    plot(TreatmentIntercept[,5],DetectionIntercept[,5])
    plot(TreatmentIntercept[,4],DetectionIntercept[,4])
    plot(TreatmentIntercept[,1],DetectionIntercept[,1])
    plot(TreatmentIntercept[,3],DetectionIntercept[,3])

