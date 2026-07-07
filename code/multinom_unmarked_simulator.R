## =================================================
##
## Title: multinom_unmarked_simulator.R
##
## Author: Jasmine Williamson
## Date Created: 7/06/2026
##
## Description: Multinomial N-mixture (removal sampling) simulator
##
##  1. Simulates true abundance (N) 
##  2. Simulates removal sampling observations
##  3. Fits multinomPois (unmarked) to the simulated data
##  4. Compares recovered estimates to the true simulated values
##  5. Repeats simulations 50 times and tracks convergence/spread
##  Goal: see if we can recreate the simulated parameter values by
##  modelling the simulated data, and whether ~30 captures over 127
##  sites is enough to do that reliably...
##
##  =================================================


## settings ---------------------------------------------------

  library(unmarked)

## ------------------------------------------------------------
## 1. SIMULATE WORLD
## ------------------------------------------------------------

  # I = sites
  I <- 127
  # K = removal passes (occasions) per site visit
  K <- 3
  
  # assign each site to one of 5 treatments, roughly evenly split
  treatment_levels <- c("control", "burn", "harvest", "salvage", "harvestburn")
  treatment <- sample(rep(treatment_levels, length.out = I))   # shuffles I sites into 5 roughly-equal groups
  table(treatment)
  
  # dummy variables for treatments (control is the reference)
  burn     <- as.numeric(treatment == "burn")
  harvest  <- as.numeric(treatment == "harvest")
  salvage  <- as.numeric(treatment == "salvage")
  harvestburn <- as.numeric(treatment == "harvestburn")
  
  # covariates
  canopycover <- rnorm(I) # simulates canopy cov values at I sites from a normal dist with mean=0 and sd=1 (scaled)
  downedwood <- runif(I, -2, 2) # simulates dwd values at I sites from a uniform dist from -2 to 2
  days_since_rain <- array(scale(rgamma(n = I*K, shape = 2, rate = 0.5)), dim = c(I, K)) # right-skewed, non-neg, scaled, site*rep dimensions

  # set true parameter values
  
  # on log-scale for lambda; use exp(x) to see the value on the abundance scale
  # p uses logit-link, so plogis(alpha0) = per-pass detection probability
  
  beta0.lambda        <- log(0.65)   # mean abundance/site in control sites (30: 0.35), (50: 0.75), (60: 0.9)
  beta.burn.lambda     <- log(1)      # burn effect: no change vs control 
  beta.harvest.lambda  <- log(0.4)    # harvest effect: strong negative 
  beta.salvage.lambda  <- log(0.25)   # salvage effect: strong negative 
  beta.harvestburn.lambda <- log(0.3)    # harv-burn effect: negative
  beta.canopy.lambda   <- 0.3         # positive: more canopy -> more sals
  beta.dwd.lambda      <- 0.2         # positive: more downed wood -> more sals
  alpha0       <- qlogis(0.42)   # mean per-pass detection prob, plogis(alpha0) = 0.4 det prob
  alpha.dsr    <- -0.4           # negative: more days since rain -> drier -> lower detection
  
  
  # calculate lambda (expected abundance) at each site (full covs)
  lambda <- matrix(ncol = I)
  for (i in 1:I) {
    lambda[i] <- exp(beta0.lambda + beta.burn.lambda * burn[i] +
                       beta.harvest.lambda * harvest[i] +
                       beta.salvage.lambda * salvage[i] +
                       beta.harvestburn.lambda * harvestburn[i] +
                       beta.canopy.lambda * canopycover[i] +
                       beta.dwd.lambda * downedwood[i])
  }
  
  
  # intercept-only version:
  # beta0.lambda <- log(0.29)     # mean site abundance -> exp(beta0.lambda) ~ 0.25 individuals/site
  # alpha0       <- qlogis(0.42)   # mean per-pass detection prob -> plogis(alpha0) = 0.4
  # lambda <- matrix(ncol = I)
  # for (i in 1:I) {
  #   lambda[i] <- exp(beta0.lambda)
  # }
  
  
  # calculate true abundance N at each site (Poisson)
  N <- matrix(ncol = I)
  for (i in 1:I) {
    N[i] <- rpois(1, lambda[i])
  }
  
  sum(N)                       # true total abundance
  tapply(N, treatment, sum)    # distribution of true per-site abundance
  
  # calculate detection probability p at each site/pass using a logit-link
  p <- array(0, dim = c(I, K))
  for (i in 1:I) {
    for (k in 1:K) {
      p[i, k] <- plogis(alpha0 + alpha.dsr * days_since_rain[i,k])
    }
  }

  
## ------------------------------------------------------------
## 2. SIMULATE OBSERVATIONS
## ------------------------------------------------------------
  
  # simulate removal-sampling observations
  
  # on pass k, each still-uncaptured individual at site i is independently
  # caught with probability p[i,k]; anyone caught is removed from the pool
  # and can't be caught again on a later pass within the same visit
  
  y <- matrix(0, nrow = I, ncol = K)
  for (i in 1:I) {
    remaining <- N[i]                                     #all N[i] individuals available
    for (k in 1:K) {
      caught       <- rbinom(1, remaining, p[i, k])       # how many remaining individuals get caught this pass
      y[i, k]      <- caught                              # record it
      remaining    <- remaining - caught                  # remove them
    }
  }
  
  sum(y)          # total captures
                  # adjust beta0.lambda / alpha0 above and rerun until sum(y) =~30 (or goal)
  
  y[, 1]          # captures on pass 1 across all sites
  y[, 2]          
  y[, 3]          
  
  str(y)
  tapply(rowSums(y), treatment, sum)
  

## ------------------------------------------------------------
## 3. MODEL
## ------------------------------------------------------------

  # intercept-only version:
  umf <- unmarkedFrameMPois(y = y, type = "removal")
  # fit: ~1 ~1 intercept-only lambda and p
  fit <- multinomPois(~1 ~1, data = umf)
  
  
  # treatment as site-level covariate:
  umf1 <- unmarkedFrameMPois(y = y, type = "removal",
                            siteCovs = data.frame(treatment = factor(treatment, levels = treatment_levels)),
                            obsCovs = 
  )
  
  fit1 <- multinomPois(~1 ~ treatment, data = umf1)
  summary(fit1)
  
  
  # more covs model:
  umf2 <- unmarkedFrameMPois(y = y, type = "removal",
                            siteCovs = data.frame(
                              treatment = factor(treatment, levels = treatment_levels),
                              canopycover = canopycover,
                              downedwood = downedwood
                            ),
                            obsCovs = list(dsr = days_since_rain))
  
  fit2 <- multinomPois(~dsr ~ treatment + canopycover + downedwood, data = umf2)
  summary(fit2)


## ------------------------------------------------------------
## 4. OUTPUT
## ------------------------------------------------------------

  # for full model:
  # recovered estimates, back-transformed from log/logit
  
  coef(fit2, type = "state")    # abundance (log-scale)
  SE(fit2, type = "state")
  
  coef(fit2, type = "det")      # detection (logit-scale)
  SE(fit2, type = "det")
  
  # true values for comparison - order must match coef(fit2, type="state")
  true_betas_fit2 <- c("(Intercept)"          = beta0.lambda,
                       "treatmentburn"        = beta.burn.lambda,
                       "treatmentharvest"     = beta.harvest.lambda,
                       "treatmentsalvage"     = beta.salvage.lambda,
                       "treatmentharvestburn" = beta.harvestburn.lambda,
                       "canopycover"          = beta.canopy.lambda,
                       "downedwood"           = beta.dwd.lambda)
  true_betas_fit2
  
  true_alphas_fit2 <- c("(Intercept)" = alpha0,
                        "dsr"          = alpha.dsr)
  true_alphas_fit2
  
  # predicted abund
  newdat_abund <- data.frame(
    treatment   = factor(treatment_levels, levels = treatment_levels),
    canopycover = 0,
    downedwood  = 0
  )
  predicted_lambda <- predict(fit2, type = "state", newdata = newdat)
  cbind(newdat, predicted_lambda)
  
  
  
  # intercept-only model:
  # backTransform() works here because this model has one estimate per side (no covs)
  lambda.est <- backTransform(fit, type = "state")
  p.est      <- backTransform(fit, type = "det")
  
  lambda.est@estimate   # recovered mean abund/site
  exp(beta0.lambda)     # true mean abund/site - compare these 
  
  p.est@estimate        # recovered det prob
  plogis(alpha0)        # true det prob - compare these 
  
  # even if estimate above is close, huge SE means the CI is too wide to be useful
  SE(lambda.est)
  SE(p.est)
  
  
## ------------------------------------------------------------
## 5. REPEAT SIM
## ------------------------------------------------------------

  # Repeats the whole simulate -> fit sequence n_reps times, storing
  # coefficients each time, so we can see the spread of results across runs.

  n_reps <- 50   # start smaller (10-20) to check it runs, then increase
  
  state_names <- c("(Intercept)", "treatmentburn", "treatmentharvest",
                   "treatmentsalvage", "treatmentharvestburn",
                   "canopycover", "downedwood")
  det_names   <- c("(Intercept)", "dsr")
  
  results_state <- matrix(NA, nrow = n_reps, ncol = length(state_names),
                          dimnames = list(NULL, state_names))
  results_det   <- matrix(NA, nrow = n_reps, ncol = length(det_names),
                          dimnames = list(NULL, det_names))
  converged     <- logical(n_reps)
  
  for (rep in 1:n_reps) {
    
    ## re-simulate the world (section 1)
    treatment <- sample(rep(treatment_levels, length.out = I))
    burn        <- as.numeric(treatment == "burn")
    harvest     <- as.numeric(treatment == "harvest")
    salvage     <- as.numeric(treatment == "salvage")
    harvestburn <- as.numeric(treatment == "harvestburn")
    
    canopycover     <- rnorm(I)
    downedwood      <- runif(I, -2, 2)
    days_since_rain <- array(scale(rgamma(n = I * K, shape = 2, rate = 0.5)), dim = c(I, K))
    
    lambda <- matrix(ncol = I)
    for (i in 1:I) {
      lambda[i] <- exp(beta0.lambda + beta.burn.lambda * burn[i] +
                         beta.harvest.lambda * harvest[i] +
                         beta.salvage.lambda * salvage[i] +
                         beta.harvestburn.lambda * harvestburn[i] +
                         beta.canopy.lambda * canopycover[i] +
                         beta.dwd.lambda * downedwood[i])
    }
    
    N <- matrix(ncol = I)
    for (i in 1:I) {
      N[i] <- rpois(1, lambda[i])
    }
    
    p <- array(0, dim = c(I, K))
    for (i in 1:I) {
      for (k in 1:K) {
        p[i, k] <- plogis(alpha0 + alpha.dsr * days_since_rain[i, k])
      }
    }
    
    ## re-simulate observations (section 2)
    y <- matrix(0, nrow = I, ncol = K)
    for (i in 1:I) {
      remaining <- N[i]
      for (k in 1:K) {
        caught    <- rbinom(1, remaining, p[i, k])
        y[i, k]   <- caught
        remaining <- remaining - caught
      }
    }
    
    ## re-fit the model (same spec as fit2)
    umf_rep <- unmarkedFrameMPois(y = y, type = "removal",
                                  siteCovs = data.frame(
                                    treatment   = factor(treatment, levels = treatment_levels),
                                    canopycover = canopycover,
                                    downedwood  = downedwood
                                  ),
                                  obsCovs = list(dsr = days_since_rain))
    
    fit_rep <- tryCatch(
      multinomPois(~dsr ~ treatment + canopycover + downedwood, data = umf_rep),
      error   = function(e) NULL,
      warning = function(w) NULL
    )
    
    if (is.null(fit_rep)) {
      converged[rep] <- FALSE
    } else {
      converged[rep] <- TRUE
      results_state[rep, ] <- coef(fit_rep, type = "state")
      results_det[rep, ]   <- coef(fit_rep, type = "det")
    }
  }
  
  
  ## summarize spread across reps
  state_summary <- data.frame(
    true_value = c(beta0.lambda, beta.burn.lambda, beta.harvest.lambda,
                   beta.salvage.lambda, beta.harvestburn.lambda,
                   beta.canopy.lambda, beta.dwd.lambda),
    mean_estimate = apply(results_state, 2, mean, na.rm = TRUE),
    sd_across_reps = apply(results_state, 2, sd, na.rm = TRUE)
  )
  rownames(state_summary) <- state_names

  
  det_summary <- data.frame(
    true_value = c(alpha0, alpha.dsr),
    mean_estimate = apply(results_det, 2, mean, na.rm = TRUE),
    sd_across_reps = apply(results_det, 2, sd, na.rm = TRUE)
  )
  rownames(det_summary) <- det_names

  
  
  # boxplots with true value marked in red: want each box centered on red dot (true value)
  boxplot(results_state, main = "Abundance-side coefficients across reps",
          las = 2, ylab = "Estimate (log-scale)")
  points(1:length(state_names), state_summary$true_value, col = "red", pch = 19)
  
  boxplot(results_det, main = "Detection-side coefficients across reps",
          las = 2, ylab = "Estimate (logit-scale)")
  points(1:length(det_names), det_summary$true_value, col = "red", pch = 19)
  
  
  # assess: 
  mean(converged)   # fraction of reps that converged
  print(state_summary) # look at mean est vs true value, want them the same
  print(det_summary) # also look at sd across reps, how much could one rep be off by chance?

  

#### Notes --------------------------------------------------------------------
  
  # Running the sim replicate 50x with 30 captures:
  # some covs are fine, with decent estimates but big-ish SE's
  # treatments that had very few (<5ish) captures had huge SE's and
  # estimates are way off
  
  # with 50 captures:
  # things look better, but start to see bad SE's on the treatment category
  # with low captures; one/some of the reps created a big outlier thats throwing it,
  # similar to the 30 capture sim
  
  # with 62 captures:
  # estimates are close and SE's are reasonable
  
  
  
  
  
  
  
