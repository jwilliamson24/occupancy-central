## =================================================
##
## Title: plot-level-model
## Author: Jasmine Williamson
## Date Created: 06/30/2025
##
## Description: Occupancy model for 2023-2024 data 
## Took the hierarchical simulation model from JT and reduced it,
## to model 2023-24 data at plot level instead of site level (to increase power)
##
## =================================================
# rm(list=ls())

  load("data/plot lvl occu rdata/e-plot_lvl_output_and_data_072525.rdata")

## load packages
  library(nimble)
  library(coda)
  library(dplyr)

## load data 
  # from occu-matrices-postfire.R
  dets.o <- read.csv("data/occupancy/dets.o.post.csv") # post-fire oss
  dets.e <- read.csv("data/occupancy/dets.e.post.csv") # post-fire enes
  
  # 2023-2024 covariate data
  site.lvl <- read.csv("data/covariate matrices/site_aspect_precip_all_vars.csv") # all site level covariates
  subplot.lvl <- read.csv("data/covariate matrices/subplot.complete.new.csv") # covariates at subplot lvl 

  
## format  --------------------------------------------------------------------
  
# site covariates 
  
  # landowner to management type column
  site.lvl <- site.lvl %>%
    mutate(mgmt_type = case_when(
      landowner == "PB" ~ "0", # private
      landowner == "WY" ~ "0",  # private
      landowner == "BLM" ~ "1",  # public
      landowner == "ODF" ~ "1",  # public
      TRUE ~ "default_value"  
    ))  
  
  # add treatment columns
  treatments <- unique(site.lvl$trt)
  for (trt in treatments) {
    site.lvl[[trt]] <- as.numeric(site.lvl$trt == trt)
  }  
  
  # subset site, date, owner, mgmt, treatments
  site.info <- site.lvl %>%
    select(site_id, jul_date, landowner, mgmt_type, HB,BU,BS,HU,UU, precip_mm, days_since_rain)
  
  site.info$mgmt_type <- as.numeric(site.info$mgmt_type)

  
# subplot covariates
  
  # site, subplot, temp, lat, long, elev
  subplot.info <- subplot.lvl %>%
    select(site_id, subplot, temp, lat, long, elev, veg_cov, soil_moist_avg, logs, decay_cl, canopy_cov)
  
  names(subplot.info)[names(subplot.info) == "logs"] <-  "DW"
  
  # change decay_cl=NA to =0
  subplot.info$decay_cl[is.na(subplot.info$decay_cl)] <- 0
  
  
# add to detection matrices
  
  # oss  
  occu.o <- dets.o %>%
    left_join(site.info, by = "site_id") %>%
    left_join(subplot.info, by = c("site_id", "subplot")) 
  
  # enes
  occu.e <- dets.e %>%
    left_join(site.info, by = "site_id") %>%
    left_join(subplot.info, by = c("site_id", "subplot")) 

  
# scale all numeric covariates
  covs_scaled <- occu.e
  covs_scaled[sapply(covs_scaled, is.numeric)] <- scale(covs_scaled[sapply(covs_scaled, is.numeric)])

  
  
## Global Model written -------------------------------------------------------
  
  # psi(occupancy covariates) p(detection covariates)
  
  # psi( trt + canopy_cov + dwd_count + lat + long + elev + veg_cov + 
  #      soil_moist + decay_cl + mgmt_type )
  
  # p( temp + precip_mm + days_since_rain )
  
  
  
#### Model --------------------------------------------------------------------

  ## Define 
  I <- 889 # I = sites (subplots as sites)
  K <- 3 # K = occasions
  y = occu.e[,c(6:8)]   ################# choose species here
  y2 = occu.e       ################## and here
  covs = covs_scaled
  
# run the chains
  n.chains = 3
  chains = vector("list", n.chains)
  for(chain in 1:n.chains){
    
    #fit model
    constants <- list(I = I, K = K, 
                      HU = y2$HU, HB = y2$HB, BU = y2$BU, BS = y2$BS,
                      canopy = covs$canopy_cov,
                      DW = covs$DW,
                      lat = covs$lat,
                      lon = covs$long,
                      elev = covs$elev,
                      veg = covs$veg_cov,
                      soilmoist = covs$soil_moist_avg,
                      decay = covs$decay_cl,
                      #mgmt = y2$mgmt_type,
                      temp = covs$temp,
                      precip = covs$precip_mm,
                      raindays = covs$days_since_rain
                      )
    
    Nimdata <- list(y=y)
    
    # set initial values
    Niminits <- list( beta0.psi=0, beta1.psi.BU=0, beta2.psi.HB=0, beta3.psi.HU=0, beta4.psi.BS=0,  
                      beta5.psi.lat=0, beta6.psi.lon=0, #beta7.psi.mgmt=0, 
                      beta8.psi.elev=0, beta9.psi.DW=0, 
                      beta10.psi.cc=0, beta11.psi.veg=0, beta12.psi.soil=0, beta13.psi.dec=0, 
                      alpha0=0, alpha1.t=0, alpha2.t=0, alpha3.precip=0, alpha4.rdays=0)
    
    # data summaries
     z.data <- 1*(rowSums(y, na.rm = TRUE)>0)
     z.data[z.data==0] <- NA
    #z.data <- ifelse(rowSums(y, na.rm = TRUE) > 0, 1, 0)  # changed this to try to stop NA's problem
    
    
    # provide data for model
    Nimdata <- list(y=y, z=z.data)
    
    #set parameters to monitor
    parameters <- c("beta0.psi", "beta1.psi.BU", "beta2.psi.HB", "beta3.psi.HU", "beta4.psi.BS", 
                    "beta5.psi.lat", "beta6.psi.lon", #"beta7.psi.mgmt", 
                    "beta8.psi.elev", "beta9.psi.DW", 
                    "beta10.psi.cc", "beta11.psi.veg", "beta12.psi.soil", "beta13.psi.dec", 
                    "alpha0", "alpha1.t", "alpha2.t", "alpha3.precip", "alpha4.rdays", 'zsum')
  
    
NimModel <- nimbleCode({
      
      # priors
      beta0.psi ~ dlogis(0,1) # mean occupancy (intercept)
      beta1.psi.BU ~ dnorm(0, sd = 5)  # slope term for burn effect on site occu
      beta2.psi.HB ~ dnorm(0, sd = 5)  # harvest burn
      beta3.psi.HU ~ dnorm(0, sd = 5)  # harvest
      beta4.psi.BS ~ dnorm(0, sd = 5)  # burn salvage
      
      beta5.psi.lat ~ dnorm(0, sd = 5) # latitude
      beta6.psi.lon ~ dnorm(0, sd = 5) # longitude
      #beta7.psi.mgmt ~ dnorm(0, sd = 5) # management type
      beta8.psi.elev ~ dnorm(0, sd = 5) # elevation
      beta9.psi.DW ~ dnorm(0, sd = 5) # downed wood 
      beta10.psi.cc ~ dnorm(0, sd = 5) # canopy cover
      beta11.psi.veg ~ dnorm(0, sd = 5) # veg cover
      beta12.psi.soil ~ dnorm(0, sd = 5) # soil moisture
      beta13.psi.dec ~ dnorm(0, sd = 5) # decay class
      
      alpha0 ~ dlogis(0,1)  # prior for detection intercept, mean det prob
      alpha1.t ~ dnorm(0, sd = 5) # linear temp effect on detection
      alpha2.t ~ dnorm(0, sd = 5) # quadratic temp
      alpha3.precip ~ dnorm(0, sd = 5) # precipitation amount
      alpha4.rdays ~ dnorm(0, sd = 5) # days since rain
          
      
      # likelihood for state model
      for(i in 1:I){
        
          logit(psi[i]) <- beta0.psi + beta1.psi.BU * BU[i] + beta2.psi.HB * HB[i] + 
            beta3.psi.HU * HU[i] + beta4.psi.BS * BS[i] + beta5.psi.lat * lat[i] + 
            beta6.psi.lon * lon[i] + #beta7.psi.mgmt * mgmt[i] + 
            beta8.psi.elev * elev[i] +
            beta9.psi.DW * DW[i] + beta10.psi.cc * canopy[i] + beta11.psi.veg * veg[i] +
            beta12.psi.soil * soilmoist[i] + beta13.psi.dec * DW[i] * decay[i] #when DW=0, decay term drops
            
        z[i] ~ dbern(psi[i]) 
        
        
      # likelihood for observation model
        for(k in 1:K){
            logit(p[i,k]) <- alpha0 + alpha1.t * temp[i] + alpha2.t * temp[i]^2 +
              alpha3.precip * precip[i] + alpha4.rdays * raindays[i]
            
            y[i,k] ~ dbern(p[i,k]*z[i]) 
          }
        }
      
      zsum <- sum(z[1:I])
      
    })# end model
    
   
  # Custom configurations from JT 
    # Build the model, configure the mcmc, and compileConfigure
    start.time <- Sys.time()
    Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
    conf <- configureMCMC(Rmodel,monitors=parameters, thin=50, useConjugacy=FALSE) # thinning interval
    
    # Build and compile
    Rmcmc <- buildMCMC(conf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    
    # Run the model
    start.time2 <- Sys.time()
    Cmcmc$run(1000000,reset=FALSE) # this sets n.iterations
    end.time <- Sys.time()
    end.time - start.time  # total time for compilation, replacing samplers, and fitting
    end.time - start.time2 # post-compilation run time
    
    #get the chains
    mvSamples = as.matrix(Cmcmc$mvSamples)
    chains[[chain]]=mvSamples
  }
  
  
  n.iter = 20000
  n.burn = 10000
  #combine the chains and burn
  a=mcmc.list(mcmc(chains[[1]][n.burn:n.iter,]),
              mcmc(chains[[2]][n.burn:n.iter,]),
              mcmc(chains[[3]][n.burn:n.iter,]))
  
  # save
  save(a, file = "e-plot_lvl_output_072525.RData")
  save(a, constants, Nimdata, NimModel, Niminits, file = "e-plot_lvl_output_and_data_072525.RData")
  
  
  
## Diagnostics ------------------------------------------------------------
  
# R-hat values 
  gelman.diag(a)
  # looks great
  
# Trace plots
  plot(a)
  # looks great
  
# Estimates
  summary(a)
  # a few estimates close to zero, potentially concerning?
  # DW is really small? cc is really high?
  
# Combine chains
  a=runjags::combine.mcmc(a)
  colnames(mvSamples)
  
  
# Correlation matrix
  cor <- cor(a[, c("alpha0", "alpha1.t", "alpha2.t", "alpha3.precip", "alpha4.rdays",   
                   "beta0.psi", "beta1.psi.BU", "beta10.psi.cc", "beta11.psi.veg", "beta12.psi.soil",
                   "beta13.psi.dec", "beta2.psi.HB", "beta3.psi.HU", "beta4.psi.BS", "beta5.psi.lat",  
                   "beta6.psi.lon", "beta8.psi.elev", "beta9.psi.DW"
                   )])
  
  # some really high correlations here -oopsies ######################
  
  threshold <- 0.6
  high_corr <- abs(cor) > threshold
  diag(high_corr) <- FALSE
  high_corr_pairs <- which(high_corr, arr.ind = TRUE)
  data.frame(
    Var1 = rownames(cor)[high_corr_pairs[, 1]],
    Var2 = colnames(cor)[high_corr_pairs[, 2]],
    Correlation = cor[high_corr]
  )
  
  # cc is moderately correlated to the harvest treatments - should i take it out?
  
  # harvest treatments are highly correlated to each other - should i combine them?
  
  # lat and long are highly correlated, chat gpt says i should remove them and 
  # add a spatial random effect instead?
  
  
# Effective sample size: 
  effectiveSize(a)
  # looks good 
  
  
  
