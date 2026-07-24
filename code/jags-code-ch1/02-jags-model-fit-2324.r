# --------------------------------------------------------------------------------------------
##
## 02-jags-model-fit-2324
##
## Task: Fit OSS/ENES JAGS Occupancy Model to 2023-2024 data
##
## Start Date: 02/19/2025
## Jasmine Williamson
##

## goals
# adapt Jay Jones' code to fit his jags model that i adapted to run 2023-2024 data

## insights
# model runs! 02/20/25
# added detection intercept by treatment, it changes the results a lot
# cant understand why its doing this, pretty sure im doing the code right for that


## settings --------------------------------------------------------------------------------

  library(R2jags)
  #library(beepr)
  
  rm(list=ls())



# --------------------------------
# 1. read in the data
# --------------------------------
  
load("resources/LS-jas data wrangling/output_files/oss_enes_data_packaged_jw.RData")
str(osslist)



# --------------------------------
# 2. load JAGS occu model
# --------------------------------

source("code/jags-code/01-jags-model-2324.r")
  

  
# --------------------------------
# 3. fit the 'hierarchical' model
# --------------------------------

# diagnostic check
myFun <- function(n, y, p){
  out <- matrix(ncol=3, nrow=n)
  miss <- is.na(y)
  for(i in 1:n){
    out[i,1] <- sum(((as.numeric(y) - rep(p[i,], 3))^2)[!miss]) # 253
    out[i,2] <- sum(((rbinom(length(miss), 1, rep(p[i,], 3)) - rep(p[i,], 3))^2)[!miss]) # 252.2
  }
  out[,3] <- ifelse(out[,1] > out[,2], 1, 0)
  out
}



# OSS

  # model setup
  zst.o <- apply(osslist$yo, 1, max, na.rm=T)
  inits.o <- function(){list(z=zst.o)}
  params.o <- c("beta0", "betaTFCL", "betaTFNC", "betaOWPB", "betaOWODF", "betaOWBLM", "betaYr24", "betaHU",
                "betaBU", "betaHB", "betaBS", "betaDW", "SalvEffect", "BurnEffect",
                "mu.a0", "mu.a1", "mu.a2", "mu.a3", "mu.a4", "mu.a5", "mu.a6", "mu.a7", "mu.a8", "mu.a9", "mu.a10",
                "mu.a11", "mu.a12", "mu.a13", "mu.a14", "p.eff", "b0", "sd.b0", "a0")
  
  test.data.o <- with(osslist, list(y=yo, nyear=nyear, nstand=nstand, nall=nall, nvisit=nvisit, 
                                    Year2024=Year2024, TFCL1=TFCL1, TFNC1=TFNC1, 
                                    OWPB1=OWPB1, OWODF1=OWODF1, OWBLM1=OWBLM1,
                                    HB2=HB2, BS2=BS2, HU2=HU2, BU2=BU2, 
                                    DW3=(DW3-3)/2, AT3=(AT3-12)/5.5,  #what is he doing here with these calculations?
                                    HB3=HB3, BS3=BS3, HU3=HU3, BU3=BU3,
                                    Year3=Year3, Stand3=Stand3))

  # fit the model
  system.time(out.oh1 <- jags(test.data.o, inits.o, params.o, model.h1, n.chains=4,
                              #n.thin=10, n.iter=1100, n.burnin=100)) # 22s
                              n.thin=10, n.iter=11000, n.burnin=1000)) # 210s

  # check diagnostics
  hist(out.oh1$BUGSoutput$summary[,8]) # ~all <1.05
  out.oh1$BUGSoutput$summary[rownames(out.oh1$BUGSoutput$summary) %in% 
                               c("betaBU", "betaHB", "betaBS", "betaHU", "betaDW", 
                                 "SalvEffect", "BurnEffect", "mu.a0", "mu.a1"),]
  
  #dont really know what this is doing
  p.eff.oh1 <- out.oh1$BUGSoutput$sims.list$p.eff
  d.oh1 <- myFun(400, y=as.numeric(osslist$yo), p=p.eff.oh1) #9/18 changed to 400 from 1000 bc error: subscript out of bounds
  summary(d.oh1) 

  
  

  
  
  
# ENES
  
  # model setup
  zst.e <- apply(osslist$ye, 1, max, na.rm=T)
  inits.e <- function(){list(z=zst.e)}
  params.e <- c("beta0", "betaTFCL", "betaTFNC", "betaOWPB", "betaOWODF", "betaOWBLM", "betaYr24", "betaHU",
                "betaBU", "betaHB", "betaBS", "betaDW", "SalvEffect", "BurnEffect",
                "mu.a0", "mu.a1", "mu.a2", "mu.a3", "mu.a4", "mu.a5", "mu.a6", "mu.a7", "mu.a8", "mu.a9", "mu.a10",
                "mu.a11", "mu.a12", "mu.a13", "mu.a14", "p.eff", "b0", "sd.b0", "a0")
  
  test.data.e <- with(osslist, list(y=ye, nyear=nyear, nstand=nstand, nall=nall, nvisit=nvisit, 
                                    Year2024=Year2024, TFCL1=TFCL1, TFNC1=TFNC1, 
                                    OWPB1=OWPB1, OWODF1=OWODF1, OWBLM1=OWBLM1,
                                    HB2=HB2, BS2=BS2, HU2=HU2, BU2=BU2, 
                                    DW3=(DW3-3)/2, AT3=(AT3-12)/5.5,  #what is he doing here with these calculations?
                                    HB3=HB3, BS3=BS3, HU3=HU3, BU3=BU3,
                                    Year3=Year3, Stand3=Stand3))
  
  # fit the model
  system.time(out.eh1 <- jags(test.data.e, inits.e, params.e, model.h1, n.chains=4,
                              #n.thin=10, n.iter=1100, n.burnin=100)) # 21s
                              n.thin=10, n.iter=11000, n.burnin=1000)) # 200s
  
  # check diagnostics
  hist(out.eh1$BUGSoutput$summary[,8]) # ~all <1.05
  out.eh1$BUGSoutput$summary[rownames(out.eh1$BUGSoutput$summary) %in% 
                               c("betaBU", "betaHB", "betaBS", "betaHU", "betaDW",
                                 "SalvEffect", "BurnEffect", "mu.a0", "mu.a1"),]
  #dont really know what this is doing
  p.eff.eh1 <- out.eh1$BUGSoutput$sims.list$p.eff
  d.eh1 <- myFun(400, y=as.numeric(osslist$ye), p=p.eff.eh1) 
  summary(d.eh1) 
  
  
  

# --------------------------
# 4. export model objects
# --------------------------
  
  save(out.oh1, out.eh1, file="jags-occu-model-objects.RData")

  






