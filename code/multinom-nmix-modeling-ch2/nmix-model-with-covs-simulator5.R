## =================================================
##
## Title: nmix_model_with_covs_simulator5.R
##
## Author: Jasmine Williamson
## Date Created: 7/15/2026
##
## Description: Age-structured multinomial N-mixture (removal sampling) simulator
## Goal: build model with full suite of hypothesized covariates and run with 
## simulated data similar to real data to see if we can recover estimates
##
##  1. Simulates true site abundance (N), split into 3 age classes
##     via a treatment-specific age-composition vector (pi_age)
##  2. Simulates removal-sampling observations per age class, using
##     age-specific detection probabilities (p_age)
##  3. Fits a Bayesian hierarchical model in NIMBLE to recover pi_age
##     (per treatment) and p_age (per age class) from the simulated data
##  4. Compares recovered estimates to the true simulated values
##
##  =================================================


## settings ---------------------------------------------------

  rm(list=ls())
  library(nimble)
  library(coda)

  ## DEBUG MODE - flip to FALSE for real runs
  ## I/n.iter/n.chains settings
  debug_mode <- FALSE

  if (debug_mode) {
    I        <- 100     # small site count - fast compile, fast run
    n.iter   <- 500      # just enough to see if it runs/errors, not converges
    n.burn   <- 100       # scaled down to match - must be < n.iter
    n.chains <- 2        # 2 = minimum needed for Gelman-Rubin
    
  } else {
    I        <- 889
    n.iter   <- 20000
    n.burn   <- 5000
    n.chains <- 3
  }


## ------------------------------------------------------------
## 1. SIMULATE WORLD
## ------------------------------------------------------------

  # K = removal passes per site visit, 3 age classes (J, SA, A)
  # (I is set above in the DEBUG MODE block, not here)
  K <- 3

  # treatment levels - UU ref lvl (beta.trt[1] <- 0 below)
  treatment_levels <- c("UU", "BS", "BU", "HB", "HU")
  ntrt <- length(treatment_levels)

  # assign sites to treatments, roughly evenly split (~50 sites/trt)
  treatment <- sample(rep(treatment_levels, length.out = I))
  trt <- match(treatment, treatment_levels)   # numeric index 1:ntrt for NIMBLE
  table(treatment)


  ## site-level covariates for abundance model -----
  # canopy cover, dwd count, decay class, char class (z-scored) 
  canopy_cov <- rnorm(I, 0, 1)
  dwd_count  <- rnorm(I, 0, 1)
  decay_cl   <- rnorm(I, 0, 1)
  char_cl    <- rnorm(I, 0, 1)


  ## true abundance model (log-scale, treatment + covariate effects) -----

  # these offsets are set (roughly) so that expected total captures per
  # treatment land near real OSS numbers (UU=71, BU=74, HB=44, HU=36, BS=27)
  # spread across ~178 sites/treatment (889 total / 5 treatments)
  
  # adjust beta0.lambda / beta.trt below and re-run section 1-2 until
  # sum(y) by treatment matches real per-treatment totals

  beta0.lambda      <- log(0.55)         # mean site abundance, UU (reference)
  beta.trt          <- c(UU = 0,
                          BS = log(27/71),
                          BU = log(74/71),
                          HB = log(44/71),
                          HU = log(30/71))

  # true covariate effects on abundance (log scale)
  beta.canopy_true <- 0.30
  beta.dwd_true    <- 0.20
  beta.decay_true  <- -0.20
  beta.char_true   <- -0.15

  lambda <- matrix(ncol = I)
  for (i in 1:I) {
    lambda[i] <- exp(beta0.lambda + beta.trt[trt[i]] +
                      beta.canopy_true * canopy_cov[i] +
                      beta.dwd_true    * dwd_count[i]  +
                      beta.decay_true  * decay_cl[i]   +
                      beta.char_true   * char_cl[i])
  }

  N <- matrix(ncol = I)
  for (i in 1:I) {
    N[i] <- rpois(1, lambda[i])
  }

  sum(N)                     # true total abundance
  tapply(N, treatment, sum)  # true abundance by treatment


  ## true age composition per treatment (pi_age) -----

  # rows = treatment (same order as treatment_levels), columns = J, SA, A
  
  # HYPOTHESIS: more disturbed treatments skew toward adults
  # adults may survive disturbance but be too impacted (body condition, habitat
  # quality) to successfully reproduce
  # UU has the most "normal"/balanced age structure; BS most severe
  
  pi_age_true <- rbind(
    UU = c(J = 0.35, SA = 0.30, A = 0.35),   # least disturbed - balanced structure
    BU = c(J = 0.25, SA = 0.30, A = 0.45),   # burn only
    HU = c(J = 0.20, SA = 0.30, A = 0.50),   # harvest only
    HB = c(J = 0.15, SA = 0.25, A = 0.60),   # harvest + burn 
    BS = c(J = 0.10, SA = 0.25, A = 0.65)    # burn + salvage
  )
  pi_age_true <- pi_age_true[treatment_levels, ]  # keep order aligned w/ trt index
  pi_age_true

  # convert pi_age_true (proportions) into phi0_true (multinomial logit
  # scale, J = reference) so it can combine with covariate effects below 
  phi0_true <- matrix(0, nrow = ntrt, ncol = 3, dimnames = list(treatment_levels, c("J","SA","A")))
  for (t in 1:ntrt) {
    phi0_true[t, ] <- log(pi_age_true[t, ] / pi_age_true[t, "J"])
  }

  # true covariate effects on age composition (multinomial logit scale,
  # relative to J) - dwd_cov added here as new variable, distinct from
  # dwd_count used in the abundance model 
  gam.canopy_true <- c(J = 0, SA = 0.20, A = 0.35)
  gam.dwdcov_true <- c(J = 0, SA = 0.15, A = 0.25)
  gam.soil_true   <- c(J = 0, SA = -0.15, A = -0.30)


## ------------------------------------------------------------
## 2. SIMULATE OBSERVATIONS
## ------------------------------------------------------------

  # detection covariates - site-level, standardized
  temp            <- rnorm(I, 0, 1)
  soil_moist      <- rnorm(I, 0, 1)
  days_since_rain <- rnorm(I, 0, 1)

  # age composition covariate - dwd cover class
  dwd_cov <- rnorm(I, 0, 1)

  # site-level true age composition: treatment baseline (phi0_true) plus
  # canopy/dwd_cov/soil covariate effects, then softmax back to proportions
  phi_true    <- array(0, dim = c(I, 3), dimnames = list(NULL, c("J","SA","A")))
  pi_age_site_true <- array(0, dim = c(I, 3), dimnames = list(NULL, c("J","SA","A")))
  for (i in 1:I) {
    for (a in 1:3) {
      phi_true[i,a] <- phi0_true[trt[i],a] + gam.canopy_true[a]*canopy_cov[i] +
                        gam.dwdcov_true[a]*dwd_cov[i] + gam.soil_true[a]*soil_moist[i]
    }
    pi_age_site_true[i,] <- exp(phi_true[i,]) / sum(exp(phi_true[i,]))
  }

  # observer - fixed effect, not random
  # claude says 5 levels is too few to estimate a stable variance component
  nobs <- 5
  obs <- sample(1:nobs, I, replace = TRUE)   # observer 1 = reference

  # detection prob by age class AND site now (was just age before)
  # baseline age-class detection plus 3 covariate effects plus obs
  # THIS WHOLE COMBINATION IS THE THING WE ARE TESTING
  p_age_true <- c(J = 0.25, SA = 0.35, A = 0.50)   # baseline, at covariates = 0, observer 1
  beta.temp_true <- 0.25
  beta.soil_true <- -0.20
  beta.days_true <- -0.15
  eps.obs_true   <- c(0, 0.20, -0.15, 0.10, -0.25)   # observer 1 = reference (0); EDIT magnitudes

  logit_p_age_true <- matrix(0, nrow = I, ncol = 3, dimnames = list(NULL, c("J","SA","A")))
  p_age_site_true  <- matrix(0, nrow = I, ncol = 3, dimnames = list(NULL, c("J","SA","A")))
  for (i in 1:I) {
    for (a in 1:3) {
      logit_p_age_true[i,a] <- log(p_age_true[a]/(1-p_age_true[a])) +
                                beta.temp_true*temp[i] + beta.soil_true*soil_moist[i] +
                                beta.days_true*days_since_rain[i] +
                                eps.obs_true[obs[i]]
      p_age_site_true[i,a] <- 1/(1+exp(-logit_p_age_true[i,a]))
    }
  }

  # split total N at each site into age classes using the site's
  # true age composition (now site-level: treatment + canopy/dwd_cov/soil)
  N_age <- matrix(0, nrow = I, ncol = 3, dimnames = list(NULL, c("J","SA","A")))
  for (i in 1:I) {
    N_age[i, ] <- rmultinom(1, N[i], pi_age_site_true[i, ])
  }

  colSums(N_age)                       # true total captures-eligible by age
  tapply(rowSums(N_age), treatment, sum)


  # simulate removal-sampling observations, separately per age class
  # p now varies by site too
  y <- array(0, dim = c(I, K, 3), dimnames = list(NULL, NULL, c("J","SA","A")))
  for (i in 1:I) {
    for (a in 1:3) {
      remaining <- N_age[i, a]
      for (k in 1:K) {
        caught       <- rbinom(1, remaining, p_age_site_true[i,a])
        y[i, k, a]   <- caught
        remaining    <- remaining - caught
      }
    }
  }

  sum(y)                                   # total captures
  apply(y, 3, sum)                         # total captures by age class
  tapply(apply(y, 1, sum), treatment, sum)
  
  # site x age matrix
  captures_by_trt_age <- apply(y, c(1,3), sum)                 
  capture_table <- t(sapply(treatment_levels, function(tr) {
    colSums(captures_by_trt_age[treatment == tr, , drop = FALSE])
  }))
  capture_table <- cbind(capture_table, total = rowSums(capture_table))
  capture_table
  
  
## ------------------------------------------------------------
## 3. MODEL
## ------------------------------------------------------------

  # NimbleCode: total abundance (N) with trt effect on lambda,
  # split into age classes via trt-specific pi_age (softmax link),
  # observed via sequential removal-sampling with age-specific detection.
  # Detection is estimated once per age class, shared across treatment -
  # this lets us tell "true age composition differs by 
  # treatment" apart from "detection differs by age"

  NimModel <- nimbleCode({
    
    ## ABUNDANCE MODEL (lambda)
    # estimates abundance per trt, each trt independently, plus 4 site
    # covariate slopes - same tight-prior logic throughout (keeps lambda
    # from wandering to implausible values, same fix used earlier in
    # this script for the ridge/identifiability problem)
    beta0 ~ dnorm(0, sd = 0.5)            # tight prior
    beta.trt[1] <- 0                      # UU = reference treatment
    for (t in 2:ntrt) {                   # trt priors
      beta.trt[t] ~ dnorm(0, sd = 0.5)
    }
    beta.canopy ~ dnorm(0, sd = 0.5)
    beta.dwd    ~ dnorm(0, sd = 0.5)
    beta.decay  ~ dnorm(0, sd = 0.5)
    beta.char   ~ dnorm(0, sd = 0.5)
    
    
    ## AGE COMPOSITION MODEL (pi_age) - PRIORS ONLY
    # pi_age now varies by site too (canopy/dwd_cov/soil moisture are
    # site-level), so the softmax equation moved into the likelihood loop
    # below (Block 4) instead of living here as its own standalone block
    for (t in 1:ntrt) {
      phi0[t,1] <- 0                       # J = reference age class
      phi0[t,2] ~ dnorm(0, sd = 2)
      phi0[t,3] ~ dnorm(0, sd = 2)
    }
    gam.canopy[1] <- 0                     # J = reference age class
    gam.canopy[2] ~ dnorm(0, sd = 1)
    gam.canopy[3] ~ dnorm(0, sd = 1)
    gam.dwdcov[1] <- 0
    gam.dwdcov[2] ~ dnorm(0, sd = 1)
    gam.dwdcov[3] ~ dnorm(0, sd = 1)
    gam.soil[1]   <- 0
    gam.soil[2]   ~ dnorm(0, sd = 1)
    gam.soil[3]   ~ dnorm(0, sd = 1)
    
    
    ## DETECTION MODEL (p_age) - PRIORS ONLY
    # p_age now varies by site too (temp/soil moisture/days-since-rain are
    # site-level), so the actual equation moved into the likelihood loop
    # below (Block 4) instead of living here as its own standalone block -
    # can't compute a site-varying quantity outside the site loop
    mu.p ~ dnorm(-0.6, sd = 0.5)          # logit(0.35) =~ -0.6
    eps.p[1] <- 0                         # J = reference for the offset
    eps.p[2] ~ dnorm(0, sd = 1)
    eps.p[3] ~ dnorm(0, sd = 1)
    beta.temp ~ dnorm(0, sd = 0.5)
    beta.soil ~ dnorm(0, sd = 0.5)
    beta.days ~ dnorm(0, sd = 0.5)
    eps.obs[1] <- 0                       # observer 1 = reference
    for (o in 2:nobs) {
      eps.obs[o] ~ dnorm(0, sd = 1)
    }
    
    
    ## LIKELIHOOD 
    # combines above 3 model pieces to say how many of each age exist per site,
    # and how many are caught across 3 removal passes, using simulated data y
    for (i in 1:I) {
      log(lambda[i]) <- beta0 + beta.trt[trt[i]] +
                         beta.canopy*canopy_cov[i] + beta.dwd*dwd_count[i] +
                         beta.decay*decay_cl[i]    + beta.char*char_cl[i]
      
      # age composition equation - now site AND age specific (was just
      # treatment before). Reuses canopy_cov (from abundance) and
      # soil_moist (from detection) - same covariate, different role here
      for (a in 1:3) {
        phi[i,a] <- phi0[trt[i],a] + gam.canopy[a]*canopy_cov[i] +
                    gam.dwdcov[a]*dwd_cov[i] + gam.soil[a]*soil_moist[i]
        exp.phi[i,a] <- exp(phi[i,a])
      }
      pi_age[i,1:3] <- exp.phi[i,1:3] / sum(exp.phi[i,1:3])
      
      for (a in 1:3) {
        N_age[i,a] ~ dpois(lambda[i] * pi_age[i,a])
        
        # detection equation - now site AND age specific
        logit(p_age[i,a]) <- mu.p + eps.p[a] +
                              beta.temp*temp[i] + beta.soil*soil_moist[i] +
                              beta.days*days_since_rain[i] + eps.obs[obs[i]]
        
        # sequential removal, pass 1-3, same logic as the simulation
        avail[i,1,a] <- N_age[i,a]
        y[i,1,a] ~ dbin(p_age[i,a], avail[i,1,a])
        avail[i,2,a] <- avail[i,1,a] - y[i,1,a]
        y[i,2,a] ~ dbin(p_age[i,a], avail[i,2,a])
        avail[i,3,a] <- avail[i,2,a] - y[i,2,a]
        y[i,3,a] ~ dbin(p_age[i,a], avail[i,3,a])
      }
      N[i] <- sum(N_age[i,1:3])
    }
    
    
    ## DIAGNOSTICS (not part of the model - monitoring only)
    Ntotal <- sum(N[1:I])
    for (t in 1:ntrt) {
      Ntrt[t] <- inprod(N[1:I], trtmat[1:I, t])
    }
    # "baseline" detection per age class at average covariate values (0,
    # since covariates are standardized) - a stand-in for the old p_age[a],
    # useful for a quick sanity check since actual p_age is now 889 x 3
    for (a in 1:3) {
      logit(p_age_baseline[a]) <- mu.p + eps.p[a]
    }
    # "baseline" age composition per treatment at average covariate values
    # (0) - a stand-in for the old pi_age[t,a], since actual pi_age is now
    # 889 x 3
    for (t in 1:ntrt) {
      for (a in 1:3) {
        exp.phi0[t,a] <- exp(phi0[t,a])
      }
      pi_age_baseline[t,1:3] <- exp.phi0[t,1:3] / sum(exp.phi0[t,1:3])
    }
  })
  
  
  ## SETUP: constants, data, inits
  # trtmat: site x treatment 0/1 matrix, feeds the diagnostic Ntrt[] block only
  trtmat <- matrix(0, nrow = I, ncol = ntrt)
  for (t in 1:ntrt) trtmat[trt == t, t] <- 1
  
  # covariates are known/fixed site attributes -> constants, not data
  constants <- list(I = I, K = K, ntrt = ntrt, trt = trt, trtmat = trtmat,
                     canopy_cov = canopy_cov, dwd_count = dwd_count,
                     decay_cl = decay_cl, char_cl = char_cl,
                     temp = temp, soil_moist = soil_moist, days_since_rain = days_since_rain,
                     obs = obs, nobs = nobs, dwd_cov = dwd_cov)
  Nimdata   <- list(y = y)
  
  # N_age inits from observed captures + buffer; N has no inits (derived)
  N_age.init <- apply(y, c(1,3), sum) + 2
  
  Niminits <- list(
    beta0 = 0, phi0 = matrix(0, ntrt, 3), mu.p = 0,
    beta.canopy = 0, beta.dwd = 0, beta.decay = 0, beta.char = 0,
    beta.temp = 0, beta.soil = 0, beta.days = 0,
    N_age = N_age.init
  )
  Niminits$beta.trt   <- rep(NA, ntrt)      # beta.trt[1] fixed in model code
  Niminits$eps.p      <- rep(NA, 3)         # eps.p[1] fixed in model code
  Niminits$eps.obs    <- rep(NA, nobs)      # eps.obs[1] fixed in model code
  Niminits$gam.canopy <- rep(NA, 3)         # gam.canopy[1] fixed in model code (J = ref)
  Niminits$gam.dwdcov <- rep(NA, 3)         # gam.dwdcov[1] fixed in model code
  Niminits$gam.soil   <- rep(NA, 3)         # gam.soil[1] fixed in model code
  
  # Ntotal/Ntrt monitored purely as a convergence check. NOTE: pi_age and
  # p_age are now 889 x 3 (site x age) - too many nodes to monitor
  # directly, so we monitor the _baseline versions instead (at average
  # covariate values / treatment only)
  parameters <- c("pi_age_baseline", "p_age_baseline", "beta0", "beta.trt", "Ntotal", "Ntrt",
                   "beta.canopy", "beta.dwd", "beta.decay", "beta.char",
                   "beta.temp", "beta.soil", "beta.days", "eps.obs",
                   "gam.canopy", "gam.dwdcov", "gam.soil")
  
  
  ## BUILD + COMPILE
  Rmodel <- nimbleModel(code = NimModel, constants = constants,
                        data = Nimdata, inits = Niminits, check = FALSE)
  conf   <- configureMCMC(Rmodel, monitors = parameters, thin = 2, useConjugacy = FALSE)
  Rmcmc  <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc  <- compileNimble(Rmcmc, project = Rmodel)
  
  
  ## FIT: chains from different starting values (n.iter/n.burn/n.chains
  ## set in the DEBUG MODE block at the top of the script)
  # chains agreeing despite different starts = convergence; 
  # chains landing in different places = not converged
  chain_samples <- vector("list", n.chains)
  set.seed(NULL)
  
  for (chain in 1:n.chains) {
    
    N_age.init.chain    <- N_age.init * sample(1:3, 1)   # 1x/2x/3x perturbation
    phi0.init.chain     <- matrix(rnorm(ntrt*3, 0, 1), ntrt, 3)
    phi0.init.chain[,1] <- 0
    
    Cmodel$N_age <- N_age.init.chain   # N derived from this automatically
    Cmodel$beta0 <- rnorm(1, 0, 1)
    Cmodel$phi0  <- phi0.init.chain
    Cmodel$mu.p  <- rnorm(1, 0, 1)
    Cmodel$beta.canopy <- rnorm(1, 0, 1)
    Cmodel$beta.dwd    <- rnorm(1, 0, 1)
    Cmodel$beta.decay  <- rnorm(1, 0, 1)
    Cmodel$beta.char   <- rnorm(1, 0, 1)
    Cmodel$beta.temp   <- rnorm(1, 0, 1)
    Cmodel$beta.soil   <- rnorm(1, 0, 1)
    Cmodel$beta.days   <- rnorm(1, 0, 1)
    Cmodel$eps.obs[2:nobs]    <- rnorm(nobs - 1, 0, 1)   # obs 1 fixed at 0 in model code
    Cmodel$gam.canopy[2:3]    <- rnorm(2, 0, 1)          # index 1 (J) fixed at 0
    Cmodel$gam.dwdcov[2:3]    <- rnorm(2, 0, 1)
    Cmodel$gam.soil[2:3]      <- rnorm(2, 0, 1)
    
    Cmcmc$run(n.iter, reset = TRUE)
    chain_samples[[chain]] <- as.matrix(Cmcmc$mvSamples)
  }
  
  # combine into mcmc.list + Gelman-Rubin (psrf ~1.0 = converged, >1.1 = not)
  # built dynamically so this works whether n.chains is 2 (debug) or 3 (real)
  post.burn <- (n.burn/2 + 1):nrow(chain_samples[[1]])
  a <- as.mcmc.list(lapply(chain_samples, function(cs) mcmc(cs[post.burn, ])))
  

  ## CHECK
  gelman <- gelman.diag(a, multivariate = FALSE)
  print(gelman)   # look for psrf column ~1.0; if > 1.1 = red flag

  # traceplots for the diagnostic Ntotal node + one pi_age/p_age param -
  plot(a[, "Ntotal"])
  plot(a[, "p_age_baseline[3]"])     # adult baseline detection

  # combine chains for point estimates
  samples <- as.matrix(runjags::combine.mcmc(a))

  # true total abundance for comparison (from section 1, before any chains ran)
  cat("true total N:", sum(N), "\n")
  cat("posterior mean Ntotal:", mean(samples[, "Ntotal"]), "\n")
  cat("posterior 95% CI Ntotal:", quantile(samples[, "Ntotal"], c(0.025, 0.975)), "\n\n")

  for (t in 1:ntrt) {
    cat(treatment_levels[t], "- true N:", sum(N[trt == t]),
        " | posterior mean:", round(mean(samples[, paste0("Ntrt[", t, "]")]), 1), "\n")
  }


## ------------------------------------------------------------
## 4. OUTPUT (single run check)
## ------------------------------------------------------------

  # recovered baseline pi_age (at avg covariate values) vs. true pi_age_true,
  # by treatment - pi_age itself is now 889x3 (site-varying), so this
  # compares the treatment-level baseline instead
  for (t in 1:ntrt) {
    cat("\n--", treatment_levels[t], "--\n")
    for (a in 1:3) {
      col <- paste0("pi_age_baseline[", t, ", ", a, "]")
      est <- mean(samples[, col])
      ci  <- quantile(samples[, col], probs = c(0.025, 0.975))
      cat(sprintf("  age %d: true = %.2f, est = %.2f, 95%% CI = [%.2f, %.2f]\n",
                   a, pi_age_true[t, a], est, ci[1], ci[2]))
    }
  }

  # recovered baseline p_age (at covariates = 0) vs. true baseline
  for (a in 1:3) {
    col <- paste0("p_age_baseline[", a, "]")
    est <- mean(samples[, col])
    ci  <- quantile(samples[, col], probs = c(0.025, 0.975))
    cat(sprintf("p_age_baseline[%d]: true = %.2f, est = %.2f, 95%% CI = [%.2f, %.2f]\n",
                 a, p_age_true[a], est, ci[1], ci[2]))
  }

  # recovered detection covariate slopes vs. true
  det_covs_true <- c(beta.temp = beta.temp_true, beta.soil = beta.soil_true,
                      beta.days = beta.days_true)
  for (nm in names(det_covs_true)) {
    est <- mean(samples[, nm])
    ci  <- quantile(samples[, nm], probs = c(0.025, 0.975))
    cat(sprintf("%s: true = %.2f, est = %.2f, 95%% CI = [%.2f, %.2f]\n",
                 nm, det_covs_true[nm], est, ci[1], ci[2]))
  }

  # recovered observer effects vs. true (observer 1 is the fixed reference,
  # so only 2:nobs are actually estimated)
  for (o in 2:nobs) {
    col <- paste0("eps.obs[", o, "]")
    est <- mean(samples[, col])
    ci  <- quantile(samples[, col], probs = c(0.025, 0.975))
    cat(sprintf("eps.obs[%d]: true = %.2f, est = %.2f, 95%% CI = [%.2f, %.2f]\n",
                 o, eps.obs_true[o], est, ci[1], ci[2]))
  }

  # recovered abundance covariate slopes vs. true
  covs_true <- c(beta.canopy = beta.canopy_true, beta.dwd = beta.dwd_true,
                 beta.decay = beta.decay_true, beta.char = beta.char_true)
  for (nm in names(covs_true)) {
    est <- mean(samples[, nm])
    ci  <- quantile(samples[, nm], probs = c(0.025, 0.975))
    cat(sprintf("%s: true = %.2f, est = %.2f, 95%% CI = [%.2f, %.2f]\n",
                 nm, covs_true[nm], est, ci[1], ci[2]))
  }

  # recovered age composition covariate slopes vs. true (age 1 = J = ref,
  # fixed at 0, so only ages 2:3 are actually estimated)
  gam_true_list <- list(gam.canopy = gam.canopy_true, gam.dwdcov = gam.dwdcov_true,
                         gam.soil = gam.soil_true)
  for (nm in names(gam_true_list)) {
    for (a in 2:3) {
      col <- paste0(nm, "[", a, "]")
      est <- mean(samples[, col])
      ci  <- quantile(samples[, col], probs = c(0.025, 0.975))
      cat(sprintf("%s: true = %.2f, est = %.2f, 95%% CI = [%.2f, %.2f]\n",
                   col, gam_true_list[[nm]][a], est, ci[1], ci[2]))
    }
  }



  
  
  
  
  
  
  
  
