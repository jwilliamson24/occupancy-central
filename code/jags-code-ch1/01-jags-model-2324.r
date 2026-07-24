# -------------------------------------------------------------------------------------------------------
##
## 01-jags-model-2324
##
## Task: Define OSS/ENES JAGS Occupancy Models 2023-2024
##
## Start Date: 09/23/2024
## Jasmine Williamson
##

## goals
# adapt Jay Jones' model to run with just 2023-2024 data
# define models that will be fit in following script

## insights
# model runs

# -----------------------------------------------------------------------------------------
# 'hierarchical' model (i.e., random effects model)
#
#   occupancy:  stand-level effects:  tree farm, ownership
#               stand-year level:     year, treatment group 
#               plot-level effects:   downed wood count
#               random effects:       random intercept by stand
#   detection:  plot-level effects:   treatment group, AT, AT^2, interactions
#               random effects:       all plot-level effects by year
# -----------------------------------------------------------------------------------------

#rm(list=ls())

# 
model.h1 <- function(){
  
  # priors for occupancy model
  beta0 ~ dnorm(0, 0.333) 
  betaTFCL ~ dnorm(0, 1) # tree farm: CL, NC, (SP = ref)
  betaTFNC ~ dnorm(0, 1)
  betaOWPB ~ dnorm(0, 1) # ownership: PB, ODF, BLM, (WY = ref)
  betaOWODF ~ dnorm(0, 1)
  betaOWBLM ~ dnorm(0, 1)
  betaYr24 ~ dnorm(0, 1) # 2023 = ref
  betaHU ~ dnorm(0, 0.25) # treatments: HU, BU, HB, BS, (control = ref)
  betaBU ~ dnorm(0, 0.25) # unsure about prior values
  betaHB ~ dnorm(0, 0.25)
  betaBS ~ dnorm(0, 0.25)
  betaDW ~ dnorm(0, 0.25) # dwd count
  sd.b0 ~ dgamma(1, 2) # did not change 
  tau.b0 <- 1/(sd.b0 * sd.b0) # did not change
  
  # treatment effect estimators
  SalvEffect <- betaBS - betaBU 
  BurnEffect <- betaHB - betaHU
  
  # priors for detection random effects
  mu.a0 ~ dnorm(0, 0.333) ### mu/sd/tau terms for all my trt betas, including interactions, kept same priors as Jay
  mu.a1 ~ dnorm(0, 1) # harv
  mu.a2 ~ dnorm(0, 1) # AT
  mu.a3 ~ dnorm(0, 1) # AT^2
  mu.a4 ~ dnorm(0, 1) # harv*AT
  mu.a5 ~ dnorm(0, 1) # harv*AT^2
  
  mu.a6 ~ dnorm(0, 1) # burn
  mu.a7 ~ dnorm(0, 1) # burn*AT
  mu.a8 ~ dnorm(0, 1) # burn*AT^2
  
  mu.a9 ~ dnorm(0, 1) # harvburn
  mu.a10 ~ dnorm(0, 1) # harvburn*AT
  mu.a11 ~ dnorm(0, 1) # harvburn*AT^2
  
  mu.a12 ~ dnorm(0, 1) # salv
  mu.a13 ~ dnorm(0, 1) # salv*AT
  mu.a14 ~ dnorm(0, 1) # salv*AT^2
  
  sd.a0 ~ dgamma(1, 2) # used in tau below; one for each alpha term
  sd.a1 ~ dgamma(2, 1) # Jay: hyper-prior for Trt effect; allow for greater prior uncertainty than other effects. what to do with this?
  sd.a2 ~ dgamma(1, 2)
  sd.a3 ~ dgamma(1, 2)
  sd.a4 ~ dgamma(1, 2)
  sd.a5 ~ dgamma(1, 2)
  sd.a6 ~ dgamma(1, 2)
  sd.a7 ~ dgamma(1, 2)
  sd.a8 ~ dgamma(1, 2)
  sd.a9 ~ dgamma(1, 2)
  sd.a10 ~ dgamma(1, 2)
  sd.a11 ~ dgamma(1, 2)
  sd.a12 ~ dgamma(1, 2)
  sd.a13 ~ dgamma(1, 2)
  sd.a14 ~ dgamma(1, 2)
  
  tau.a0 <- 1/(sd.a0 * sd.a0) # one for each alpha term
  tau.a1 <- 1/(sd.a1 * sd.a1)
  tau.a2 <- 1/(sd.a2 * sd.a2)
  tau.a3 <- 1/(sd.a3 * sd.a3)
  tau.a4 <- 1/(sd.a4 * sd.a4)
  tau.a5 <- 1/(sd.a5 * sd.a5)
  tau.a6 <- 1/(sd.a6 * sd.a6)
  tau.a7 <- 1/(sd.a7 * sd.a7)
  tau.a8 <- 1/(sd.a8 * sd.a8)
  tau.a9 <- 1/(sd.a9 * sd.a9)
  tau.a10 <- 1/(sd.a10 * sd.a10)
  tau.a11 <- 1/(sd.a11 * sd.a11)
  tau.a12 <- 1/(sd.a12 * sd.a12)
  tau.a13 <- 1/(sd.a13 * sd.a13)
  tau.a14 <- 1/(sd.a14 * sd.a14)

  # detection model random effects (by year)
  for(i in 1:nyear){ # used all terms from above - did i do this right? took out the "trt" object
    a0[i] ~ dnorm(mu.a0, tau.a0) # each plot level effect by year
    aHU[i] ~ dnorm(mu.a1, tau.a1) # harvest
    aAT[i] ~ dnorm(mu.a2, tau.a2) # ambient temp
    aAT2[i] ~ dnorm(mu.a3, tau.a3) # ambient temp squared (quadratic)
    aHUAT[i] ~ dnorm(mu.a4, tau.a4) # harv/temp interaction
    aHUAT2[i] ~ dnorm(mu.a5, tau.a5) # harv/temp^2 interaction
    aBU[i] ~ dnorm(mu.a6, tau.a6) # burned
    aBUAT[i] ~ dnorm(mu.a7, tau.a7) 
    aBUAT2[i] ~ dnorm(mu.a8, tau.a8) 
    aHB[i] ~ dnorm(mu.a9, tau.a9) # harvested and burned
    aHBAT[i] ~ dnorm(mu.a10, tau.a10) 
    aHBAT2[i] ~ dnorm(mu.a11, tau.a11) 
    aBS[i] ~ dnorm(mu.a12, tau.a12) # salvage logged
    aBSAT[i] ~ dnorm(mu.a13, tau.a13) 
    aBSAT2[i] ~ dnorm(mu.a14, tau.a14) 
  }
  
  for(i in 1:nstand){ 
    # occupancy stand-level effects
    mu1[i] <- beta0 + betaTFCL*TFCL1[i] + betaTFNC*TFNC1[i] + # added multiple tree farm designations
      betaOWPB*OWPB1[i] + betaOWODF*OWODF1[i] + betaOWBLM*OWBLM1[i] # added ownership designations
    mu1i[i] ~ dnorm(mu1[i], tau.b0)
    b0[i] <- mu1i[i] - mu1[i] 
    for(k in 1:nyear){
    # occupancy stand-year level effects
      mu2ik[i,k] <- mu1i[i] + betaYr24*Year2024[k,i] + betaHU*HU2[k,i] + betaBU*BU2[k,i] + 
        betaHB*HB2[k,i] * betaBS*BS2[k,i]  
    } 
  }
  
  for(i in 1:nall){
    # occupancy plot-level effects
    mu3ikj[i] <- mu2ik[Stand3[i], Year3[i]] + betaDW*DW3[i]
    logit(psi[i]) <- mu3ikj[i]
    z[i] ~ dbern(psi[i])
    # full detection model
    logit(p[i]) <- a0[Year3[i]] + aAT[Year3[i]]*AT3[i] + aAT2[Year3[i]]*AT3[i]*AT3[i] + 
      aHU[Year3[i]]*HU3[i] + aHUAT[Year3[i]]*HU3[i]*AT3[i] + aHUAT2[Year3[i]]*HU3[i]*AT3[i]*AT3[i] + #harv
      aBU[Year3[i]]*BU3[i] + aBUAT[Year3[i]]*BU3[i]*AT3[i] + aBUAT2[Year3[i]]*BU3[i]*AT3[i]*AT3[i] + #burn
      aHB[Year3[i]]*HB3[i] + aHBAT[Year3[i]]*HB3[i]*AT3[i] + aHBAT2[Year3[i]]*HB3[i]*AT3[i]*AT3[i] + #harvburn
      aBS[Year3[i]]*BS3[i] + aBSAT[Year3[i]]*BS3[i]*AT3[i] + aBSAT2[Year3[i]]*BS3[i]*AT3[i]*AT3[i]  #salv
     
    p.eff[i] <- z[i] * p[i]                 
    
    for(j in 1:nvisit){
      # likelihood
      y[i,j] ~ dbern(p.eff[i])
    }
  }
  
}



