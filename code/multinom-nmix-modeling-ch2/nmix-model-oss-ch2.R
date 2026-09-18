## =================================================
##
## Title: nmix-model-oss-ch2.R
##
## Author: Jasmine Williamson
## Date Created: 9/10/2026
##
## Description: Age-structured multinomial N-mixture (removal sampling) model
## for Oregon slender salamanders (OSS), fitted to real field data.
## Adapted from nmix-model-with-covs-simulator6.R.
##
## Model structure:
##   Abundance: treatment (FE) + stand (RE) + year (FE) + canopy + dwd + decay + char
##   Age composition: treatment baseline + canopy + dwd_cov + soil_moist
##   Detection: age class + temp + soil_moist + days_since_rain + observer (random intercept)
##     Observer levels (6): JW, JB, BZ, RMM, SLG, VO (CH + LS combined)
##     Zero-centered around mu.p; sigma.obs estimated from data
##   Observer varies by pass (obs[I,K] matrix) - each pass has its own searcher
##
## Data sources (all from multinom-data-formatting.R pipeline):
##   pass-level-counts-o.csv  -- site/subplot/pass structure + observer per pass
##   sals.complete.csv        -- individual records with age class
##   subplot.complete.new.csv -- site-level covariates
##   env_subset_corr2.csv     -- days_since_rain (site-level, repeated across subplots)
##
## =================================================


## settings ---------------------------------------------------

  rm(list=ls())
  setwd("/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git")

  library(nimble)
  library(coda)
  library(dplyr)

  ## DATASET — "oss" uses pass-level-counts-o.csv, "enes" uses pass-level-counts-e.csv
  dataset <- "oss"

  ## DEBUG MODE - flip to FALSE for real runs
  debug_mode <- FALSE

  if (debug_mode) {
    n.iter   <- 500
    n.burn   <- 100
    n.chains <- 2
  } else {
    n.iter   <- 20000
    n.burn   <- 5000
    n.chains <- 3
  }

  K <- 3  # removal passes per subplot


## ------------------------------------------------------------
## 1. LOAD DATA
## ------------------------------------------------------------

  counts_file <- switch(dataset,
    oss  = "abundance-ch2/data/abundance/pass-level-counts-o.csv",
    enes = "abundance-ch2/data/abundance/pass-level-counts-e.csv"
  )
  counts.o <- read.csv(counts_file,
                       colClasses = c(stand = "character", obs = "character"))
  counts.o$subplot <- as.integer(counts.o$subplot)
  counts.o$pass    <- as.integer(counts.o$pass)

  sals <- read.csv("occupancy-central/data/occupancy/sals.complete.csv",
                   colClasses = c(stand = "character", subplot = "character",
                                  pass = "character", age_class = "factor",
                                  recap = "character", spp = "factor"))

  subplot.dat <- read.csv("occupancy-central/data/covariate matrices/subplot.complete.new.csv",
                           colClasses = c(stand = "character", subplot = "integer"))


## ------------------------------------------------------------
## 2. PREP SITES AND INDICES
## ------------------------------------------------------------

  ## unique site list — one row per subplot visit (site_id × subplot) -----
  # site_id encodes site + site_rep + year, so each row is one unique visit
  sites <- counts.o %>%
    distinct(site_id, stand, trt, year, subplot) %>%
    arrange(site_id, subplot) %>%
    mutate(site_idx = row_number())
  I <- nrow(sites)
  cat("I (number of sites):", I, "\n")

  # treatment index
  treatment_levels <- c("UU", "BS", "BU", "HB", "HU")
  ntrt <- length(treatment_levels)
  trt  <- match(sites$trt, treatment_levels)
  stopifnot(!any(is.na(trt)))

  # stand index — numeric, 1:nstand
  stand_levels <- sort(unique(sites$stand))
  nstand <- length(stand_levels)
  stand  <- match(sites$stand, stand_levels)
  cat("nstand:", nstand, "\n")

  # year index — 1:nyear (2023 = 1, 2024 = 2, ...)
  year_levels <- sort(unique(sites$year))
  nyear <- length(year_levels)
  year  <- match(sites$year, year_levels)
  cat("years:", year_levels, " -> nyear:", nyear, "\n")

  # trtmat: site x treatment 0/1 matrix for Ntrt diagnostic
  trtmat <- matrix(0, nrow = I, ncol = ntrt)
  for (t in 1:ntrt) trtmat[trt == t, t] <- 1


## ------------------------------------------------------------
## 3. BUILD y ARRAY [I, K, 3]
## ------------------------------------------------------------

  ## age-class counts from sals.complete.csv -----
  # keep OSS only, remove recaps, keep known age classes
  sals.oss <- sals %>%
    filter(spp == "OSS",
           recap != "1",
           age_class %in% c("J", "SA", "A")) %>%
    mutate(subplot = as.integer(subplot),
           pass    = as.integer(pass))

  # count per site_id × subplot × pass × age_class
  sal_counts <- sals.oss %>%
    group_by(site_id, subplot, pass, age_class) %>%
    summarise(count = n(), .groups = "drop")

  # report how many unknown-age individuals are excluded
  n_unknown <- sals %>%
    filter(spp == "OSS", recap != "1", age_class == "U") %>%
    nrow()
  cat("OSS individuals with unknown age class excluded:", n_unknown, "\n")

  # build y array, filling zeros for non-detections
  age_levels <- c("J", "SA", "A")
  y <- array(0L, dim = c(I, K, 3), dimnames = list(NULL, NULL, age_levels))

  for (r in seq_len(nrow(sal_counts))) {
    si <- which(sites$site_id == sal_counts$site_id[r] &
                sites$subplot  == sal_counts$subplot[r])
    k  <- sal_counts$pass[r]
    a  <- match(as.character(sal_counts$age_class[r]), age_levels)
    if (length(si) == 1 && !is.na(k) && k >= 1 && k <= K && !is.na(a)) {
      y[si, k, a] <- sal_counts$count[r]
    }
  }

  cat("Total OSS captures (known age):", sum(y), "\n")
  cat("Captures by age class:\n"); print(apply(y, 3, sum))
  cat("Captures by treatment:\n"); print(tapply(apply(y, 1, sum), sites$trt, sum))


## ------------------------------------------------------------
## 4. BUILD obs MATRIX [I, K]  (6 levels, random effect)
## ------------------------------------------------------------

  # observer varies by pass — each of the 3 searchers per subplot does one pass
  # 6 observer-ID levels, treated as a random intercept in the detection submodel:
  #   1 = JW   2 = JB   3 = BZ   4 = RMM   5 = SLG   6 = VO (CH + LS combined)
  obs_lookup <- c(JW = 1L, JB = 2L, BZ = 3L, RMM = 4L, SLG = 5L, CH = 6L, LS = 6L)

  obs.long <- counts.o %>%
    select(site_id, subplot, pass, obs) %>%
    left_join(sites %>% select(site_id, subplot, site_idx),
              by = c("site_id", "subplot")) %>%
    mutate(obs_idx = obs_lookup[obs])

  # check: every raw initial maps to the right index
  obs.long %>%
    distinct(obs, obs_idx) %>%
    arrange(obs_idx, obs) %>%
    print()

  obs_mat <- matrix(NA_integer_, nrow = I, ncol = K)
  for (r in seq_len(nrow(obs.long))) {
    obs_mat[obs.long$site_idx[r], obs.long$pass[r]] <- obs.long$obs_idx[r]
  }
  stopifnot(!any(is.na(obs_mat)))
  nobs <- 6   # JW, JB, BZ, RMM, SLG, VO
  obs_labels <- c("JW", "JB", "BZ", "RMM", "SLG", "VO")
  cat("Observer pass counts:\n")
  print(table(factor(obs_mat, levels = 1:nobs, labels = obs_labels)))


## ------------------------------------------------------------
## 5. COVARIATES
## ------------------------------------------------------------

  ## join subplot-level covariates to site list -----
  site_covs <- sites %>%
    left_join(subplot.dat %>% select(site_id, subplot, canopy_cov, dwd_count, dwd_cov,
                                     decay_cl, char_cl, fwd_cov, veg_cov,
                                     temp, soil_moist_avg, jul_date),
              by = c("site_id", "subplot"))
  stopifnot(nrow(site_covs) == I)

  ## fill jul_date NAs from subplot.complete.csv (27 subplots with no DWD have NA jul_date
  ## from the DWD-derived source in subplot.complete.new.csv; all subplots have a survey date)
  subplot.base <- read.csv("data/covariate matrices/subplot.complete.csv",
                            colClasses = c(stand = "character")) %>%
    mutate(subplot = as.integer(subplot)) %>%
    select(site_id, subplot, date) %>%
    mutate(jul_date_survey = as.numeric(format(as.Date(date, "%m/%d/%Y"), "%j")))
  site_covs <- site_covs %>%
    left_join(subplot.base %>% select(site_id, subplot, jul_date_survey),
              by = c("site_id", "subplot")) %>%
    mutate(jul_date = dplyr::coalesce(jul_date, jul_date_survey)) %>%
    select(-jul_date_survey)

  ## days_since_rain from env_subset_corr2.csv -----
  # this file is at site level (one row per site_id, 127 sites × 7 subplots = 889)
  # joining on site_id repeats each site's value across its 7 subplots automatically
  env2 <- read.csv("data/covariate matrices/env_subset_corr2.csv")
  site_covs <- site_covs %>%
    left_join(env2 %>% select(site_id, days_since_rain, avg_volume), by = "site_id")

  cat("NAs in covariates:\n"); print(colSums(is.na(site_covs %>% select(-site_idx))))


  ## z-score all continuous covariates -----
  zs <- function(x) as.numeric(scale(x))

  canopy_cov      <- zs(site_covs$canopy_cov)
  dwd_count       <- zs(site_covs$dwd_count)
  decay_cl        <- zs(site_covs$decay_cl)
  char_cl         <- zs(site_covs$char_cl)
  fwd_cov         <- zs(site_covs$fwd_cov)
  veg_cov         <- zs(site_covs$veg_cov)
  avg_volume      <- zs(site_covs$avg_volume)
  dwd_cov         <- zs(site_covs$dwd_cov)
  temp            <- zs(site_covs$temp)
  temp_sq         <- zs(site_covs$temp^2)
  soil_moist      <- zs(site_covs$soil_moist_avg)
  jul_date        <- zs(site_covs$jul_date)
  days_since_rain <- zs(site_covs$days_since_rain)
  
  decay_cl[is.na(decay_cl)] <- 0
  char_cl[is.na(char_cl)]   <- 0

  # check to make sure there are no NA's left
  cov_list <- list(canopy_cov = canopy_cov, dwd_count = dwd_count,
                   decay_cl = decay_cl, char_cl = char_cl, fwd_cov = fwd_cov,
                   veg_cov = veg_cov, avg_volume = avg_volume, dwd_cov = dwd_cov,
                   temp = temp, temp_sq = temp_sq, soil_moist = soil_moist,
                   jul_date = jul_date, days_since_rain = days_since_rain)
  na_counts <- sapply(cov_list, function(x) sum(is.na(x)))
  if (any(na_counts > 0)) {
    cat("WARNING: NAs remain in covariates after imputation:\n")
    print(na_counts[na_counts > 0])
    stop("Fix NA covariates before passing to NIMBLE.")
  } else {
    cat("Covariate NA check passed — no NAs in any z-scored covariate.\n")
  }


## ------------------------------------------------------------
## 6. MODEL
## ------------------------------------------------------------

  NimModel <- nimbleCode({

    ## ABUNDANCE MODEL (lambda) -----
    beta0 ~ dnorm(0, sd = 0.5)
    beta.trt[1] <- 0                      # UU = reference treatment
    for (t in 2:ntrt) {
      beta.trt[t] ~ dnorm(0, sd = 0.5)
    }
    beta.canopy ~ dnorm(0, sd = 0.5)
    beta.dwd.count    ~ dnorm(0, sd = 0.5)
    beta.decay  ~ dnorm(0, sd = 0.5)
    beta.char   ~ dnorm(0, sd = 0.5)
    beta.fwd    ~ dnorm(0, sd = 0.5)
    beta.veg    ~ dnorm(0, sd = 0.5)
    beta.vol    ~ dnorm(0, sd = 0.5)

    # stand random effect
    sigma.stand ~ dexp(1)
    for (s in 1:nstand) {
      alpha.stand[s] ~ dnorm(0, sd = sigma.stand)
    }

    # survey year fixed effect (year 1 = reference)
    beta.year[1] <- 0
    for (yr in 2:nyear) {
      beta.year[yr] ~ dnorm(0, sd = 0.5)
    }


    ## AGE COMPOSITION MODEL (pi_age) — PRIORS ONLY -----
    for (t in 1:ntrt) {
      phi0[t,1] <- 0                       # J = reference age class
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
    gam.fwd[1]    <- 0
    gam.fwd[2]    ~ dnorm(0, sd = 1)
    gam.fwd[3]    ~ dnorm(0, sd = 1)


    ## DETECTION MODEL (p_age) — PRIORS ONLY -----
    # detection varies by site × pass × age via observer tier per pass
    mu.p ~ dnorm(-0.6, sd = 0.5)          # logit(0.35) ≈ -0.6
    eps.p[1] <- 0                          # J = reference age class
    eps.p[2] ~ dnorm(0, sd = 1)
    eps.p[3] ~ dnorm(0, sd = 1)
    beta.temp  ~ dnorm(0, sd = 0.5)
    beta.temp2 ~ dnorm(0, sd = 0.5)       # quadratic temperature term
    beta.soil  ~ dnorm(0, sd = 0.5)
    beta.days  ~ dnorm(0, sd = 0.5)       # days since rain
    beta.jul   ~ dnorm(0, sd = 0.5)       # julian date
    # observer random intercept — 6 levels (JW, JB, BZ, RMM, SLG, VO)
    # zero-centered around mu.p so deviations are relative to the baseline
    # intercept, not to any single reference observer
    sigma.obs ~ dexp(1)
    for (o in 1:nobs) {
      eps.obs[o] ~ dnorm(0, sd = sigma.obs)
    }


    ## LIKELIHOOD -----
    for (i in 1:I) {

      # abundance
      log(lambda[i]) <- beta0 + beta.trt[trt[i]] +
                         alpha.stand[stand[i]] +
                         beta.year[year[i]] +
                         beta.canopy*canopy_cov[i] + beta.dwd.count*dwd_count[i] +
                         beta.decay*decay_cl[i]    + beta.char*char_cl[i] +
                         beta.fwd*fwd_cov[i]       + beta.veg*veg_cov[i] +
                         beta.vol*avg_volume[i]

      # age composition (site-level: treatment baseline + covariates)
      for (a in 1:3) {
        phi[i,a] <- phi0[trt[i],a] + gam.canopy[a]*canopy_cov[i] +
                    gam.dwdcov[a]*dwd_cov[i] + gam.soil[a]*soil_moist[i] +
                    gam.fwd[a]*fwd_cov[i]
        exp.phi[i,a] <- exp(phi[i,a])
      }
      pi_age[i,1:3] <- exp.phi[i,1:3] / sum(exp.phi[i,1:3])

      for (a in 1:3) {
        N_age[i,a] ~ dpois(lambda[i] * pi_age[i,a])

        # detection varies by pass via observer tier (obs[i,k])
        for (k in 1:K) {
          logit(p_age[i,k,a]) <- mu.p + eps.p[a] +
                                  beta.temp*temp[i] + beta.temp2*temp_sq[i] +
                                  beta.soil*soil_moist[i] +
                                  beta.days*days_since_rain[i] + beta.jul*jul_date[i] +
                                  eps.obs[obs[i,k]]
        }

        # sequential removal — pass-specific detection
        avail[i,1,a] <- N_age[i,a]
        y[i,1,a] ~ dbin(p_age[i,1,a], avail[i,1,a])
        avail[i,2,a] <- avail[i,1,a] - y[i,1,a]
        y[i,2,a] ~ dbin(p_age[i,2,a], avail[i,2,a])
        avail[i,3,a] <- avail[i,2,a] - y[i,2,a]
        y[i,3,a] ~ dbin(p_age[i,3,a], avail[i,3,a])
      }
      N[i] <- sum(N_age[i,1:3])
    }


    ## DIAGNOSTICS (monitoring only) -----
    Ntotal <- sum(N[1:I])
    for (t in 1:ntrt) {
      Ntrt[t] <- inprod(N[1:I], trtmat[1:I, t])
    }
    # baseline detection per age class at average covariates (0), observer effect = 0
    # (random effect averages to zero across observers, so this is the population-level baseline)
    for (a in 1:3) {
      logit(p_age_baseline[a]) <- mu.p + eps.p[a]
    }
    # baseline age composition per treatment at average covariates (0)
    for (t in 1:ntrt) {
      for (a in 1:3) {
        exp.phi0[t,a] <- exp(phi0[t,a])
      }
      pi_age_baseline[t,1:3] <- exp.phi0[t,1:3] / sum(exp.phi0[t,1:3])
    }
  })


## ------------------------------------------------------------
## 7. NIMBLE SETUP
## ------------------------------------------------------------

  constants <- list(
    I = I, K = K, ntrt = ntrt, nstand = nstand, nyear = nyear, nobs = nobs,
    trt = trt, trtmat = trtmat, stand = stand, year = year,
    obs = obs_mat,
    canopy_cov = canopy_cov, dwd_count = dwd_count,
    decay_cl = decay_cl, char_cl = char_cl, fwd_cov = fwd_cov, veg_cov = veg_cov,
    avg_volume = avg_volume,
    dwd_cov = dwd_cov,
    temp = temp, temp_sq = temp_sq, soil_moist = soil_moist, jul_date = jul_date,
    days_since_rain = days_since_rain
  )
  Nimdata <- list(y = y)

  # inits: N_age from observed captures + small buffer
  N_age.init <- apply(y, c(1,3), sum) + 2

  Niminits <- list(
    beta0       = 0,
    phi0        = matrix(0, ntrt, 3),
    mu.p        = 0,
    beta.canopy = 0, beta.dwd.count = 0, beta.decay = 0, beta.char = 0,
    beta.fwd    = 0, beta.veg = 0, beta.vol = 0,
    beta.temp   = 0, beta.temp2 = 0, beta.soil = 0, beta.days = 0, beta.jul = 0,
    sigma.stand = 0.5,
    alpha.stand = rep(0, nstand),
    N_age       = N_age.init
  )
  Niminits$beta.trt   <- rep(NA, ntrt)    # [1] fixed in model code
  Niminits$beta.year  <- rep(NA, nyear)   # [1] fixed in model code
  Niminits$eps.p      <- rep(NA, 3)       # [1] fixed in model code
  Niminits$eps.obs    <- rep(0, nobs)     # all free — random effect, no fixed reference
  Niminits$sigma.obs  <- 0.5
  Niminits$gam.canopy <- rep(NA, 3)       # [1] fixed in model code (J = ref)
  Niminits$gam.dwdcov <- rep(NA, 3)       # [1] fixed in model code
  Niminits$gam.soil   <- rep(NA, 3)       # [1] fixed in model code
  Niminits$gam.fwd    <- rep(NA, 3)       # [1] fixed in model code

  parameters <- c(
    "pi_age_baseline", "p_age_baseline",
    "beta0", "beta.trt", "Ntotal", "Ntrt",
    "beta.canopy", "beta.dwd.count", "beta.decay", "beta.char", "beta.fwd", "beta.veg", "beta.vol",
    "sigma.stand",
    "beta.year",
    "beta.temp", "beta.temp2", "beta.soil", "beta.days", "beta.jul",
    "sigma.obs", "eps.obs",
    "gam.canopy", "gam.dwdcov", "gam.soil", "gam.fwd"
  )


## ------------------------------------------------------------
## 8. BUILD + COMPILE
## ------------------------------------------------------------

  Rmodel <- nimbleModel(code = NimModel, constants = constants,
                        data = Nimdata, inits = Niminits, check = FALSE)
  conf   <- configureMCMC(Rmodel, monitors = parameters, thin = 2, useConjugacy = FALSE)
  Rmcmc  <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc  <- compileNimble(Rmcmc, project = Rmodel)


## ------------------------------------------------------------
## 9. FIT
## ------------------------------------------------------------

  chain_samples <- vector("list", n.chains)
  set.seed(NULL)

  for (chain in 1:n.chains) {

    phi0.init.chain     <- matrix(rnorm(ntrt*3, 0, 1), ntrt, 3)
    phi0.init.chain[,1] <- 0

    Cmodel$N_age        <- N_age.init * sample(1:3, 1)
    Cmodel$beta0        <- rnorm(1, 0, 1)
    Cmodel$phi0         <- phi0.init.chain
    Cmodel$mu.p         <- rnorm(1, 0, 1)
    Cmodel$beta.canopy  <- rnorm(1, 0, 1)
    Cmodel$beta.dwd.count     <- rnorm(1, 0, 1)
    Cmodel$beta.decay   <- rnorm(1, 0, 1)
    Cmodel$beta.char    <- rnorm(1, 0, 1)
    Cmodel$beta.fwd     <- rnorm(1, 0, 1)
    Cmodel$beta.veg     <- rnorm(1, 0, 1)
    Cmodel$beta.vol     <- rnorm(1, 0, 1)
    Cmodel$beta.temp    <- rnorm(1, 0, 1)
    Cmodel$beta.temp2   <- rnorm(1, 0, 1)
    Cmodel$beta.soil    <- rnorm(1, 0, 1)
    Cmodel$beta.days    <- rnorm(1, 0, 1)
    Cmodel$beta.jul     <- rnorm(1, 0, 1)
    Cmodel$sigma.stand        <- runif(1, 0.1, 1)
    Cmodel$alpha.stand        <- rnorm(nstand, 0, 0.5)
    Cmodel$beta.year[2:nyear] <- rnorm(nyear - 1, 0, 0.5)   # year 1 fixed at 0
    Cmodel$sigma.obs          <- runif(1, 0.1, 1)
    Cmodel$eps.obs            <- rnorm(nobs, 0, 0.5)         # all 6 free
    Cmodel$gam.canopy[2:3]    <- rnorm(2, 0, 1)
    Cmodel$gam.dwdcov[2:3]    <- rnorm(2, 0, 1)
    Cmodel$gam.soil[2:3]      <- rnorm(2, 0, 1)
    Cmodel$gam.fwd[2:3]       <- rnorm(2, 0, 1)

    Cmcmc$run(n.iter, reset = TRUE)
    chain_samples[[chain]] <- as.matrix(Cmcmc$mvSamples)
  }

  post.burn <- (n.burn/2 + 1):nrow(chain_samples[[1]])
  a <- as.mcmc.list(lapply(chain_samples, function(cs) mcmc(cs[post.burn, ])))

  
  out_dir <- "/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/abundance-ch2/data"
  saveRDS(chain_samples, file = file.path(out_dir, paste0("chain_samples_", dataset, ".rds")))
  saveRDS(a,             file = file.path(out_dir, paste0("mcmc_list_",     dataset, ".rds")))
  cat("Saved results for dataset:", dataset, "\n")


## ------------------------------------------------------------
## 10. CHECK
## ------------------------------------------------------------

  gelman <- gelman.diag(a, multivariate = FALSE)
  print(gelman)

  plot(a[, "Ntotal"])
  plot(a[, "p_age_baseline[3]"])

  samples <- do.call(rbind, lapply(a, as.matrix))

  cat("\nposterior mean Ntotal:", mean(samples[, "Ntotal"]), "\n")
  cat("posterior 95% CI Ntotal:", quantile(samples[, "Ntotal"], c(0.025, 0.975)), "\n\n")

  for (t in 1:ntrt) {
    cat(treatment_levels[t], "- posterior mean N:",
        round(mean(samples[, paste0("Ntrt[", t, "]")]), 1), "\n")
  }


## ------------------------------------------------------------
## 11. OUTPUT
## ------------------------------------------------------------

  cat("\n--- Baseline age composition by treatment (at avg covariates) ---\n")
  for (t in 1:ntrt) {
    cat("\n--", treatment_levels[t], "--\n")
    for (a in 1:3) {
      col <- paste0("pi_age_baseline[", t, ", ", a, "]")
      est <- mean(samples[, col])
      ci  <- quantile(samples[, col], probs = c(0.025, 0.975))
      cat(sprintf("  %s: est = %.2f, 95%% CI = [%.2f, %.2f]\n",
                   age_levels[a], est, ci[1], ci[2]))
    }
  }

  cat("\n--- Baseline detection by age class (avg covariates, observer = you) ---\n")
  for (a in 1:3) {
    col <- paste0("p_age_baseline[", a, "]")
    est <- mean(samples[, col])
    ci  <- quantile(samples[, col], probs = c(0.025, 0.975))
    cat(sprintf("  %s: est = %.2f, 95%% CI = [%.2f, %.2f]\n",
                 age_levels[a], est, ci[1], ci[2]))
  }

  cat("\n--- Detection covariates ---\n")
  for (nm in c("beta.temp", "beta.temp2", "beta.soil", "beta.days", "beta.jul")) {
    est <- mean(samples[, nm])
    ci  <- quantile(samples[, nm], probs = c(0.025, 0.975))
    cat(sprintf("  %s: est = %.2f, 95%% CI = [%.2f, %.2f]\n", nm, est, ci[1], ci[2]))
  }

  cat("\n--- Observer random effect (sigma + individual deviations) ---\n")
  est <- mean(samples[, "sigma.obs"])
  ci  <- quantile(samples[, "sigma.obs"], probs = c(0.025, 0.975))
  cat(sprintf("  sigma.obs: est = %.2f, 95%% CI = [%.2f, %.2f]\n", est, ci[1], ci[2]))
  for (o in 1:nobs) {
    col <- paste0("eps.obs[", o, "]")
    est <- mean(samples[, col])
    ci  <- quantile(samples[, col], probs = c(0.025, 0.975))
    cat(sprintf("  eps.obs[%d] (%s): est = %.2f, 95%% CI = [%.2f, %.2f]\n",
                 o, obs_labels[o], est, ci[1], ci[2]))
  }

  cat("\n--- Stand random effect variance ---\n")
  est <- mean(samples[, "sigma.stand"])
  ci  <- quantile(samples[, "sigma.stand"], probs = c(0.025, 0.975))
  cat(sprintf("  sigma.stand: est = %.2f, 95%% CI = [%.2f, %.2f]\n", est, ci[1], ci[2]))

  cat("\n--- Year effects (relative to year 1 =", year_levels[1], ") ---\n")
  for (yr in 2:nyear) {
    col <- paste0("beta.year[", yr, "]")
    est <- mean(samples[, col])
    ci  <- quantile(samples[, col], probs = c(0.025, 0.975))
    cat(sprintf("  year %d (%d): est = %.2f, 95%% CI = [%.2f, %.2f]\n",
                 yr, year_levels[yr], est, ci[1], ci[2]))
  }

  cat("\n--- Abundance covariates ---\n")
  for (nm in c("beta.canopy", "beta.dwd.count", "beta.decay", "beta.char", "beta.fwd", "beta.veg", "beta.vol")) {
    est <- mean(samples[, nm])
    ci  <- quantile(samples[, nm], probs = c(0.025, 0.975))
    cat(sprintf("  %s: est = %.2f, 95%% CI = [%.2f, %.2f]\n", nm, est, ci[1], ci[2]))
  }

  cat("\n--- Age composition covariates ---\n")
  for (nm in c("gam.canopy", "gam.dwdcov", "gam.soil", "gam.fwd")) {
    for (a in 2:3) {
      col <- paste0(nm, "[", a, "]")
      est <- mean(samples[, col])
      ci  <- quantile(samples[, col], probs = c(0.025, 0.975))
      cat(sprintf("  %s: est = %.2f, 95%% CI = [%.2f, %.2f]\n", col, est, ci[1], ci[2]))
    }
  }
