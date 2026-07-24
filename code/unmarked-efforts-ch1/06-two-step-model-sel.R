## -------------------------------------------------------------------------------------------------------
##
## 06-two-step-model-sel.R 
##
## Model selection for global model in previous script
##
## Jasmine Williamson
## Date Created: 03-18-2025
##

## goals: From Josh T
# two step model selection:

# Create a model where psi is null (~1) and detection has all the covariates of interest 
  # this model might look like ~ rainfall+ temp + avg vol + treatment ~ 1. 
  # If this converges go to step b, if this doesn't remove treatment.

# Conduct model selection on this using the function 'dredge' in MuMIn - take the top model 
  # ^^ watch out for parameter redundancy! (see Arnold, 2010). 

# Create a model using the top detection model, and all covariates of interest on psi 
  # this model might look like ~ top model (whatever came from step 2) ~ treatment + canopy + avg vol

# Conduct model selection this using function 'dredge' take top model <- check for parameter redundancy

## insights
## 


## settings -----------------------------------------------------------------------------------------------

    rm(list=ls())
    
    library(unmarked)
    library(MuMIn)
    library(AICcmodavg)


## load data----------------------------------------------------------------------------------------------

    #df_oss <- read.csv("data/occupancy/oss-for-UMF-1.csv")

    df_oss <- read.csv("https://raw.githubusercontent.com/jwilliamson24/oss-occu/refs/heads/main/data/occupancy/data-UMF-1.csv")
    #this data is unscaled
    

## create UMF object ------------------------------------------------------------------------------------
    
    UMF.oss.1 <- unmarkedFrameOccu(
      y = df_oss[, grep("^OD", names(df_oss))],  # selects cols starting with OD = oss detections
      siteCovs = df_oss[, c("trt","veg_cov","canopy_cov","soil_moist","fwd_cov","dwd_count")],
      obsCovs = list(
        temp = df_oss[, grep("temp.", names(df_oss))], # selects cols including temp-
        soilmoist = df_oss[, grep("soilmoist.", names(df_oss))],
        rain = df_oss[, grep("rain.", names(df_oss))],
        dwd_cov = df_oss[, grep("dwdcov.", names(df_oss))],
        trt = df_oss[, grep("trt.", names(df_oss))]
      )
    )
    
    
## first model ------------------------------------------------------------------------------------
    
    # psi is null (~1) and detection has all the covariates of interest
    # formula: occu(~p ~psi, data=UMF, se=TRUE)
    
    m2 <- occu( ~ scale(temp) + scale(soilmoist) + scale(rain) + dwd_cov + trt
                ~ 1,
                data = UMF.oss.1, se = TRUE)
    
    # i think this technically converges, but some estimates show high SEs compared to their magnitude 

    
    m3 <- occu( ~ scale(temp) + scale(soilmoist) + scale(rain) + dwd_cov
                ~ 1,
                data = UMF.oss.1, se = TRUE)

    # still have some estimates with relatively high SE

    
## first model selection  ---------------------------------------------------------------------------------
    
    # psi(1) and p(global)
    
    ms1 <- dredge(m3, rank = "QAIC", chat = 1.24)
    

    #  Model selection table 
    #     psi(Int)  p(Int) p(dwd_cov) p(scl(ran)) p(scl(slm)) p(scl(tmp)) df   logLik  QAIC delta weight
    # 3    0.8664 -1.0670                -0.2545                          3 -406.577 663.8  0.00  0.133
    # 1    0.7667 -1.0020                                                 2 -407.874 663.9  0.09  0.127
    # 9    0.7930 -1.0200                                       -0.16310  3 -406.761 664.1  0.30  0.115
       
    
## second model  ---------------------------------------------------------------------------------
   
    # p(dwd) psi(global)
    
    m4 <- occu( ~ scale(rain)
                ~ trt + veg_cov + canopy_cov + scale(soil_moist) + fwd_cov + scale(dwd_count),
                data = UMF.oss.1, se = TRUE)
    
     
    # converges, estimates/SE are ok (not great)
    
    
    
## second model selection  ---------------------------------------------------------------------------------
    
    
    ms2 <- dredge(m4, rank = "AIC", fixed = c("psi(trt)", "p(scale(rain))"))
    
    # Model selection table 
    #     psi(Int) psi(cnp_cov) psi(fwd_cov) psi(scl(dwd_cnt)) psi(scl(sol_mst)) psi(trt) psi(veg_cov) p(Int) p(scl(ran)) df
    # 6  -0.44570       1.1680                        0.70950                          +              -1.084     -0.2840  9
    # 14 -0.47100       1.1880                        0.74800            0.3074        +              -1.074     -0.2360 10
    # 8   0.49160       1.2010    -0.379200           0.74460                          +              -1.077     -0.2734 10
    
    
    
    
    