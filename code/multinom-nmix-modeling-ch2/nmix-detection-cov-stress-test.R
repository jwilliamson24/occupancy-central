## =================================================
##
## Title: detection_covariate_stress_test.R
##
## Author: Jasmine Williamson
## Date Created: 7/24/2026
##
## Description: Stress-test companion to nmix_model_with_covs_simulator5.R.
##
## That script validated the full covariate model once, at a specific set
## of true detection-covariate effect sizes (temp/soil/days/observer).
## This script asks: how much does age-composition recovery (pi_age_baseline)
## degrade as those detection-covariate effects get bigger?

## This is the covariate-model analog of the detection gap test done on the
## treatment-only model. But now scaling covariate effect magnitude instead of
## a single age-class number, since detection is no longer just 3 fixed numbers.
##
## Abundance and age-composition truth (treatment, covariates, pi_age_true)
## are held identical to the validated model and never scaled - only the
## detection-covariate effect sizes change across scenarios. This isolates
## detection as the thing being stress-tested.
##
## Uses the same model structure as the validated script - only changes
## the true values used to simulate data, and how many times we refit.
##
## =================================================


## settings ---------------------------------------------------

  rm(list=ls())
  library(nimble)
  library(coda)

  ## DEBUG MODE - flip to FALSE for the real stress test
  debug_mode <- FALSE

  if (debug_mode) {
    I      <- 100    # small site count - fast compile, fast run
    n.iter <- 2000    # enough for one rep's posterior mean, not a full check
    n.burn <- 500
    n.reps <- 5       # reps per scenario - just enough to see the code work
  } else {
    I      <- 889
    n.iter <- 8000     # lighter than the 20000 used for the single validated
                       # run - here we need many reps, not one precise fit
    n.burn <- 2000
    n.reps <- 20       # start here; more reps = smoother bias estimates
  }


## ------------------------------------------------------------
## 1. SIMULATE WORLD (identical to the validated model)
## ------------------------------------------------------------

  K <- 3
  treatment_levels <- c("UU", "BS", "BU", "HB", "HU")
  ntrt <- length(treatment_levels)

  treatment <- sample(rep(treatment_levels, length.out = I))
  trt <- match(treatment, treatment_levels)

  ## abundance covariates + true effects 
  canopy_cov <- rnorm(I, 0, 1)
  dwd_count  <- rnorm(I, 0, 1)
  decay_cl   <- rnorm(I, 0, 1)
  char_cl    <- rnorm(I, 0, 1)

  beta0.lambda <- log(0.55)
  beta.trt     <- c(UU = 0, BS = log(27/71), BU = log(74/71),
                     HB = log(44/71), HU = log(30/71))
  beta.canopy_true <- 0.30
  beta.dwd_true    <- 0.20
  beta.decay_true  <- -0.20
  beta.char_true   <- -0.15

  lambda <- matrix(ncol = I)
  for (i in 1:I) {
    lambda[i] <- exp(beta0.lambda + beta.trt[trt[i]] +
                      beta.canopy_true*canopy_cov[i] + beta.dwd_true*dwd_count[i] +
                      beta.decay_true*decay_cl[i]    + beta.char_true*char_cl[i])
  }
  N <- matrix(ncol = I)
  for (i in 1:I) N[i] <- rpois(1, lambda[i])

  ## age composition truth 
  ## fixed across every detection scenario below
  pi_age_true <- rbind(
    UU = c(J = 0.35, SA = 0.30, A = 0.35),
    BU = c(J = 0.25, SA = 0.30, A = 0.45),
    HU = c(J = 0.20, SA = 0.30, A = 0.50),
    HB = c(J = 0.15, SA = 0.25, A = 0.60),
    BS = c(J = 0.10, SA = 0.25, A = 0.65)
  )
  pi_age_true <- pi_age_true[treatment_levels, ]

  phi0_true <- matrix(0, nrow = ntrt, ncol = 3, dimnames = list(treatment_levels, c("J","SA","A")))
  for (t in 1:ntrt) phi0_true[t, ] <- log(pi_age_true[t, ] / pi_age_true[t, "J"])

  gam.canopy_true <- c(J = 0, SA = 0.20, A = 0.35)
  gam.dwdcov_true <- c(J = 0, SA = 0.15, A = 0.25)
  gam.soil_true   <- c(J = 0, SA = -0.15, A = -0.30)

  ## age composition covariate + site-level true composition (fixed, as above)
  dwd_cov    <- rnorm(I, 0, 1)
  soil_moist <- rnorm(I, 0, 1)   # also used in the detection model below

  phi_true         <- array(0, dim = c(I, 3), dimnames = list(NULL, c("J","SA","A")))
  pi_age_site_true <- array(0, dim = c(I, 3), dimnames = list(NULL, c("J","SA","A")))
  for (i in 1:I) {
    for (a in 1:3) {
      phi_true[i,a] <- phi0_true[trt[i],a] + gam.canopy_true[a]*canopy_cov[i] +
                        gam.dwdcov_true[a]*dwd_cov[i] + gam.soil_true[a]*soil_moist[i]
    }
    pi_age_site_true[i,] <- exp(phi_true[i,]) / sum(exp(phi_true[i,]))
  }

  ## detection covariates (values fixed across scenarios; only the
  ## effect sizes on these get scaled per scenario below)
  temp            <- rnorm(I, 0, 1)
  days_since_rain <- rnorm(I, 0, 1)
  nobs <- 5
  obs  <- sample(1:nobs, I, replace = TRUE)

  # baseline age-class detection (unscaled - always the same, only the
  # cov effects on top of this get scaled by scenario)
  p_age_true <- c(J = 0.25, SA = 0.35, A = 0.50)


## ------------------------------------------------------------
## 2. DETECTION-EFFECT-MAGNITUDE SCENARIOS
## ------------------------------------------------------------

  # base (1x) magnitudes - these are exactly what the validated model used.
  # "none" = 0x (no covariate-driven detection variation at all, pure
  # age-class baseline only - closest analog to the old flat p_age[a]).
  # "moderate" = 1x (the validated model's actual values).
  # "large" = 2x (double every detection-covariate effect at once).
  base_beta.temp <- 0.25
  base_beta.soil <- -0.20
  base_beta.days <- -0.15
  base_eps.obs   <- c(0, 0.20, -0.15, 0.10, -0.25)

  detection_scenarios <- list(
    none     = 0,
    moderate = 1,
    large    = 2
  )


## ------------------------------------------------------------
## 3. MODEL (identical to the validated script)
## ------------------------------------------------------------

  NimModel <- nimbleCode({

    ## ABUNDANCE MODEL (lambda)
    beta0 ~ dnorm(0, sd = 0.5)
    beta.trt[1] <- 0
    for (t in 2:ntrt) {
      beta.trt[t] ~ dnorm(0, sd = 0.5)
    }
    beta.canopy ~ dnorm(0, sd = 0.5)
    beta.dwd    ~ dnorm(0, sd = 0.5)
    beta.decay  ~ dnorm(0, sd = 0.5)
    beta.char   ~ dnorm(0, sd = 0.5)

    ## AGE COMPOSITION MODEL (pi_age) - PRIORS ONLY
    for (t in 1:ntrt) {
      phi0[t,1] <- 0
      phi0[t,2] ~ dnorm(0, sd = 2)
      phi0[t,3] ~ dnorm(0, sd = 2)
    }
    gam.canopy[1] <- 0
    gam.canopy[2] ~ dnorm(0, sd = 1)
    gam.canopy[3] ~ dnorm(0, sd = 1)
    gam.dwdcov[1] <- 0
    gam.dwdcov[2] ~ dnorm(0, sd = 1)
    gam.dwdcov[3] ~ dnorm(0, sd = 1)
    gam.soil[1]   <- 0
    gam.soil[2]   ~ dnorm(0, sd = 1)
    gam.soil[3]   ~ dnorm(0, sd = 1)

    ## DETECTION MODEL (p_age) - PRIORS ONLY
    mu.p ~ dnorm(-0.6, sd = 0.5)
    eps.p[1] <- 0
    eps.p[2] ~ dnorm(0, sd = 1)
    eps.p[3] ~ dnorm(0, sd = 1)
    beta.temp ~ dnorm(0, sd = 0.5)
    beta.soil ~ dnorm(0, sd = 0.5)
    beta.days ~ dnorm(0, sd = 0.5)
    eps.obs[1] <- 0
    for (o in 2:nobs) {
      eps.obs[o] ~ dnorm(0, sd = 1)
    }

    ## LIKELIHOOD
    for (i in 1:I) {
      log(lambda[i]) <- beta0 + beta.trt[trt[i]] +
                         beta.canopy*canopy_cov[i] + beta.dwd*dwd_count[i] +
                         beta.decay*decay_cl[i]    + beta.char*char_cl[i]

      for (a in 1:3) {
        phi[i,a] <- phi0[trt[i],a] + gam.canopy[a]*canopy_cov[i] +
                    gam.dwdcov[a]*dwd_cov[i] + gam.soil[a]*soil_moist[i]
        exp.phi[i,a] <- exp(phi[i,a])
      }
      pi_age[i,1:3] <- exp.phi[i,1:3] / sum(exp.phi[i,1:3])

      for (a in 1:3) {
        N_age[i,a] ~ dpois(lambda[i] * pi_age[i,a])

        logit(p_age[i,a]) <- mu.p + eps.p[a] +
                              beta.temp*temp[i] + beta.soil*soil_moist[i] +
                              beta.days*days_since_rain[i] + eps.obs[obs[i]]

        avail[i,1,a] <- N_age[i,a]
        y[i,1,a] ~ dbin(p_age[i,a], avail[i,1,a])
        avail[i,2,a] <- avail[i,1,a] - y[i,1,a]
        y[i,2,a] ~ dbin(p_age[i,a], avail[i,2,a])
        avail[i,3,a] <- avail[i,2,a] - y[i,2,a]
        y[i,3,a] ~ dbin(p_age[i,a], avail[i,3,a])
      }
      N[i] <- sum(N_age[i,1:3])
    }

    ## DIAGNOSTICS
    for (t in 1:ntrt) {
      for (a in 1:3) {
        exp.phi0[t,a] <- exp(phi0[t,a])
      }
      pi_age_baseline[t,1:3] <- exp.phi0[t,1:3] / sum(exp.phi0[t,1:3])
    }
  })


  ## SETUP: constants, inits, monitors - built once, reused across every
  ## scenario and rep below (only `y` and the latent inits get swapped)
  trtmat <- matrix(0, nrow = I, ncol = ntrt)
  for (t in 1:ntrt) trtmat[trt == t, t] <- 1

  constants <- list(I = I, K = K, ntrt = ntrt, trt = trt, trtmat = trtmat,
                     canopy_cov = canopy_cov, dwd_count = dwd_count,
                     decay_cl = decay_cl, char_cl = char_cl,
                     temp = temp, soil_moist = soil_moist, days_since_rain = days_since_rain,
                     obs = obs, nobs = nobs, dwd_cov = dwd_cov)

  # only pi_age_baseline is monitored here (not the full recovery-check
  # list) - this stress test only cares about age-composition bias, which
  # is the thing your advisor's original question was about
  parameters <- c("pi_age_baseline")


## ------------------------------------------------------------
## 4. RUN: refit many times per scenario, scaling detection-covariate
##    effect size each time
## ------------------------------------------------------------

  pi_age_names <- paste0(rep(treatment_levels, each = 3), "_", rep(c("J","SA","A"), ntrt))
  results_all <- vector("list", length(detection_scenarios))
  names(results_all) <- names(detection_scenarios)

  Cmodel <- NULL
  Cmcmc  <- NULL

  for (scen in names(detection_scenarios)) {

    mult <- detection_scenarios[[scen]]
    beta.temp_true <- base_beta.temp * mult
    beta.soil_true <- base_beta.soil * mult
    beta.days_true <- base_beta.days * mult
    eps.obs_true   <- base_eps.obs   * mult

    pi_age_results <- matrix(NA, nrow = n.reps, ncol = length(pi_age_names),
                              dimnames = list(NULL, pi_age_names))

    for (rep in 1:n.reps) {

      ## re-simulate true site-level detection (site x age) at this
      ## scenario's effect-size scale - abundance/age-comp truth (N,
      ## pi_age_site_true) stay untouched, only detection changes
      logit_p <- array(0, dim = c(I, 3))
      p_age_site_true <- array(0, dim = c(I, 3))
      for (i in 1:I) {
        for (a in 1:3) {
          logit_p[i,a] <- log(p_age_true[a]/(1-p_age_true[a])) +
                          beta.temp_true*temp[i] + beta.soil_true*soil_moist[i] +
                          beta.days_true*days_since_rain[i] + eps.obs_true[obs[i]]
          p_age_site_true[i,a] <- 1/(1+exp(-logit_p[i,a]))
        }
      }

      ## re-simulate N_age and y (N and pi_age_site_true stay fixed)
      N_age <- matrix(0, nrow = I, ncol = 3)
      for (i in 1:I) N_age[i, ] <- rmultinom(1, N[i], pi_age_site_true[i, ])

      y <- array(0, dim = c(I, K, 3))
      for (i in 1:I) {
        for (a in 1:3) {
          remaining <- N_age[i, a]
          for (k in 1:K) {
            caught     <- rbinom(1, remaining, p_age_site_true[i,a])
            y[i, k, a] <- caught
            remaining  <- remaining - caught
          }
        }
      }

      N_age.init <- apply(y, c(1,3), sum) + 2

      if (is.null(Cmodel)) {
        # build + compile ONCE on the very first rep of the very first
        # scenario - every subsequent rep/scenario reuses this compiled
        # model via setData(), same pattern as the validated script
        Niminits <- list(
          beta0 = 0, phi0 = matrix(0, ntrt, 3), mu.p = 0,
          beta.canopy = 0, beta.dwd = 0, beta.decay = 0, beta.char = 0,
          beta.temp = 0, beta.soil = 0, beta.days = 0,
          N_age = N_age.init
        )
        Niminits$beta.trt   <- rep(NA, ntrt)
        Niminits$eps.p      <- rep(NA, 3)
        Niminits$eps.obs    <- rep(NA, nobs)
        Niminits$gam.canopy <- rep(NA, 3)
        Niminits$gam.dwdcov <- rep(NA, 3)
        Niminits$gam.soil   <- rep(NA, 3)

        Rmodel <- nimbleModel(code = NimModel, constants = constants,
                               data = list(y = y), inits = Niminits, check = FALSE)
        conf   <- configureMCMC(Rmodel, monitors = parameters, thin = 2, useConjugacy = FALSE)
        Rmcmc  <- buildMCMC(conf)
        Cmodel <- compileNimble(Rmodel)
        Cmcmc  <- compileNimble(Rmcmc, project = Rmodel)
      } else {
        Cmodel$setData(list(y = y))
        Cmodel$N_age <- N_age.init
      }

      fit_ok <- tryCatch({
        Cmcmc$run(n.iter, reset = TRUE)
        TRUE
      }, error = function(e) FALSE)

      if (!fit_ok) next

      samples_rep <- as.matrix(Cmcmc$mvSamples)
      samples_rep <- samples_rep[(n.burn/2 + 1):nrow(samples_rep), ]

      for (t in 1:ntrt) {
        for (a in 1:3) {
          col <- paste0("pi_age_baseline[", t, ", ", a, "]")
          pi_age_results[rep, paste0(treatment_levels[t], "_", c("J","SA","A")[a])] <-
            mean(samples_rep[, col])
        }
      }
    }

    results_all[[scen]] <- pi_age_results
    cat("finished scenario:", scen, "(effect multiplier =", mult, ")\n")
  }


## ------------------------------------------------------------
## 5. COMPARE SCENARIOS
## ------------------------------------------------------------

  # bias = mean(estimate) - truth, per treatment/age, per scenario -
  # same diagnostic as the original none/moderate/large test. Look
  # especially at juvenile (J) bias direction/magnitude as the detection
  # effect size scales up - that's where bias showed up before
  for (scen in names(detection_scenarios)) {
    cat("\n=====", scen, "(", detection_scenarios[[scen]], "x detection effects ) =====\n")
    pi_res <- results_all[[scen]]
    for (t in 1:ntrt) {
      for (a in 1:3) {
        col   <- paste0(treatment_levels[t], "_", c("J","SA","A")[a])
        truth <- pi_age_true[t, a]
        est   <- mean(pi_res[, col], na.rm = TRUE)
        bias  <- est - truth
        cat(sprintf("%s: true = %.2f, mean est = %.2f, bias = %+.3f\n",
                     col, truth, est, bias))
      }
    }
  }

  # boxplots - true value marked in red, one per scenario
  for (scen in names(detection_scenarios)) {
    pi_res <- results_all[[scen]]
    true_vec <- as.vector(t(pi_age_true))
    boxplot(pi_res, main = paste("pi_age_baseline recovery -", scen,
                                  "(", detection_scenarios[[scen]], "x detection effects)"),
            las = 2, ylab = "estimated proportion")
    points(1:length(pi_age_names), true_vec, col = "red", pch = 19)
  }


#### Notes --------------------------------------------------------------------

  # - does juvenile bias increase from none (0x) -> moderate (1x) -> large (2x)?
  #   model more accurately estimates age composition as det covs increase in strength
  
  # - is the pattern/direction consistent with the earlier treatment-only
  #   detection-gap test (juveniles underestimated as detection effects grow)?
  #   no - its kind of the opposite effect
  

