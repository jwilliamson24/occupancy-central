## =================================================
##
## Title: age_class_detection_simulator.R
##
## Author: Jasmine Williamson
## Date Created: 7/15/2026
##
## Description: Age-structured multinomial N-mixture (removal sampling)
## simulator, built to answer the question: how much does
## differential detection across age classes (small ones harder to see)
## distort estimated age composition (pi_age), and can the model still
## recover TRUE treatment differences in age structure once that
## detection difference is accounted for?
##
##  1. Simulates true site abundance (N), split into 3 age classes
##     (J = juvenile, SA = subadult, A = adult) via a treatment-specific
##     age-composition vector (pi_age)
##  2. Simulates removal-sampling observations per age class, using
##     age-specific detection probabilities (p_age)
##  3. Fits a Bayesian hierarchical model in NIMBLE to recover pi_age
##     (per treatment) and p_age (per age class) from the simulated data
##  4. Compares recovered estimates to the true simulated values
##  5. Repeats the whole thing across 3 detection-gap SCENARIOS
##     (no gap / moderate gap / large gap) at real-data sample sizes,
##     to see how much detection differences bias the age-composition
##     estimates, and whether treatment differences are still recoverable
##
##  =================================================


## settings ---------------------------------------------------

  rm(list=ls())
  library(nimble)
  library(coda)


## ------------------------------------------------------------
## 1. SIMULATE WORLD
## ------------------------------------------------------------

  # I = sites (site x year modeling units, matching the 889 used in your
  # treatment-only multinom_unmarked_simulator.R pipeline), K = removal
  # passes per site visit, 3 age classes (J, SA, A)
  I <- 889
  K <- 3

  # treatment levels - UU (uncut/unburned) listed first so it is the
  # reference level (beta.trt[1] <- 0 in the model below)
  treatment_levels <- c("UU", "BS", "BU", "HB", "HU")
  ntrt <- length(treatment_levels)

  # assign sites to treatments, roughly evenly split (~50 sites/trt)
  treatment <- sample(rep(treatment_levels, length.out = I))
  trt <- match(treatment, treatment_levels)   # numeric index 1:ntrt for NIMBLE
  table(treatment)


  ## true abundance model (log-scale, treatment effect on lambda):

  # these offsets are set (roughly) so that expected total captures per
  # treatment land near your real OSS numbers (UU=71, BU=74, HB=44, HU=36, BS=27)
  # spread across ~178 sites/treatment (889 total / 5 treatments) - lambda per
  # site has to be much lower than a 50-sites/treatment version to hit the
  # same real capture totals over more sites
  # ADJUST beta0.lambda / beta.trt below and re-run section 1-2 until
  # sum(y) by treatment matches your real per-treatment totals
  # (same trial-and-error approach as multinom_unmarked_simulator.R)

  beta0.lambda      <- log(0.55)         # mean site abundance, UU (reference)
  beta.trt          <- c(UU = 0,
                          BS = log(27/71),
                          BU = log(74/71),
                          HB = log(44/71),
                          HU = log(30/71))

  lambda <- matrix(ncol = I)
  for (i in 1:I) {
    lambda[i] <- exp(beta0.lambda + beta.trt[trt[i]])
  }

  N <- matrix(ncol = I)
  for (i in 1:I) {
    N[i] <- rpois(1, lambda[i])
  }

  sum(N)                     # true total abundance
  tapply(N, treatment, sum)  # true abundance by treatment


  ## true age composition per treatment (pi_age):

  # rows = treatment (same order as treatment_levels), columns = J, SA, A
  #
  # HYPOTHESIS: more disturbed treatments skew toward adults
  # adults may survive disturbance but be too impacted (body condition, habitat
  # quality) to successfully reproduce
  # UU has the most "normal"/balanced age structure; BS treated as most severe
  
  pi_age_true <- rbind(
    UU = c(J = 0.35, SA = 0.30, A = 0.35),   # least disturbed - balanced structure
    BU = c(J = 0.25, SA = 0.30, A = 0.45),   # burn only
    HU = c(J = 0.20, SA = 0.30, A = 0.50),   # harvest only
    HB = c(J = 0.15, SA = 0.25, A = 0.60),   # harvest + burn 
    BS = c(J = 0.10, SA = 0.25, A = 0.65)    # burn + salvage
  )
  pi_age_true <- pi_age_true[treatment_levels, ]  # keep order aligned w/ trt index
  pi_age_true


## ------------------------------------------------------------
## 2. SIMULATE OBSERVATIONS
## ------------------------------------------------------------

  # detection probability by age class - THIS IS THE THING WE ARE TESTING
  # start with a moderate gap as the main worked example below;
  # section 5 loops over three scenarios (none/moderate/large)
  p_age_true <- c(J = 0.25, SA = 0.35, A = 0.50)

  # split total N at each site into age classes using the site's
  # treatment-specific true age composition
  N_age <- matrix(0, nrow = I, ncol = 3, dimnames = list(NULL, c("J","SA","A")))
  for (i in 1:I) {
    N_age[i, ] <- rmultinom(1, N[i], pi_age_true[trt[i], ])
  }

  colSums(N_age)                       # true total captures-eligible by age
  tapply(rowSums(N_age), treatment, sum)


  # simulate removal-sampling observations, separately per age class
  # (same sequential removal logic as multinom_unmarked_simulator.R,
  # just repeated for each of the 3 age classes)
  y <- array(0, dim = c(I, K, 3), dimnames = list(NULL, NULL, c("J","SA","A")))
  for (i in 1:I) {
    for (a in 1:3) {
      remaining <- N_age[i, a]
      for (k in 1:K) {
        caught       <- rbinom(1, remaining, p_age_true[a])
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
    # estimates abundance per trt, each trt independently
    beta0 ~ dnorm(0, sd = 0.5)            # tight prior
    beta.trt[1] <- 0                      # UU = reference treatment
    for (t in 2:ntrt) {                   # trt priors
      beta.trt[t] ~ dnorm(0, sd = 0.5)
    }
    
    
    ## AGE COMPOSITION MODEL (pi_age) 
    # what fraction of sals are in each age class for each trt, sums to 1
    for (t in 1:ntrt) {
      phi[t,1] <- 0                       # J = reference age class
      phi[t,2] ~ dnorm(0, sd = 2)
      phi[t,3] ~ dnorm(0, sd = 2)
      for (a in 1:3) {
        exp.phi[t,a] <- exp(phi[t,a])
      }
      pi_age[t,1:3] <- exp.phi[t,1:3] / sum(exp.phi[t,1:3])
    }
    
    
    ## DETECTION MODEL (p_age) 
    # estimates detection for each age class, each pass, shared across trts
    mu.p ~ dnorm(-0.6, sd = 0.5)          # logit(0.35) =~ -0.6
    eps.p[1] <- 0                         # J = reference for the offset
    eps.p[2] ~ dnorm(0, sd = 1)
    eps.p[3] ~ dnorm(0, sd = 1)
    for (a in 1:3) {
      logit(p_age[a]) <- mu.p + eps.p[a]
    }
    
    
    ## LIKELIHOOD 
    # combines above 3 model pieces to say how many of each age exist per site,
    # and how many are caught across 3 removal passes, using simulated data y
    for (i in 1:I) {
      log(lambda[i]) <- beta0 + beta.trt[trt[i]]
      
      for (a in 1:3) {
        N_age[i,a] ~ dpois(lambda[i] * pi_age[trt[i],a])
        
        # sequential removal, pass 1-3, same logic as the simulation
        avail[i,1,a] <- N_age[i,a]
        y[i,1,a] ~ dbin(p_age[a], avail[i,1,a])
        avail[i,2,a] <- avail[i,1,a] - y[i,1,a]
        y[i,2,a] ~ dbin(p_age[a], avail[i,2,a])
        avail[i,3,a] <- avail[i,2,a] - y[i,2,a]
        y[i,3,a] ~ dbin(p_age[a], avail[i,3,a])
      }
      N[i] <- sum(N_age[i,1:3])
    }
    
    
    ## DIAGNOSTICS (not part of the model - monitoring only) -
    Ntotal <- sum(N[1:I])
    for (t in 1:ntrt) {
      Ntrt[t] <- inprod(N[1:I], trtmat[1:I, t])
    }
  })
  
  
  ## SETUP: constants, data, inits
  # trtmat: site x treatment 0/1 matrix, feeds the diagnostic Ntrt[] block only
  trtmat <- matrix(0, nrow = I, ncol = ntrt)
  for (t in 1:ntrt) trtmat[trt == t, t] <- 1
  
  constants <- list(I = I, K = K, ntrt = ntrt, trt = trt, trtmat = trtmat)
  Nimdata   <- list(y = y)
  
  # N_age inits from observed captures + buffer; N has no inits (derived)
  N_age.init <- apply(y, c(1,3), sum) + 2
  
  Niminits <- list(
    beta0 = 0, phi = matrix(0, ntrt, 3), mu.p = 0,
    N_age = N_age.init
  )
  Niminits$beta.trt <- rep(NA, ntrt)      # beta.trt[1] fixed in model code
  Niminits$eps.p    <- rep(NA, 3)         # eps.p[1] fixed in model code
  
  # Ntotal/Ntrt monitored alongside the real parameters purely as a
  # convergence check - not needed for the actual science
  parameters <- c("pi_age", "p_age", "beta0", "beta.trt", "Ntotal", "Ntrt")
  
  
  ## BUILD + COMPILE
  Rmodel <- nimbleModel(code = NimModel, constants = constants,
                        data = Nimdata, inits = Niminits, check = FALSE)
  conf   <- configureMCMC(Rmodel, monitors = parameters, thin = 2, useConjugacy = FALSE)
  Rmcmc  <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc  <- compileNimble(Rmcmc, project = Rmodel)
  
  
  ## FIT: 3 chains from different starting values
  # chains agreeing despite different starts = convergence; 
  # chains landing in different places = not converged
  n.iter <- 20000
  n.burn <- 5000
  n.chains <- 3
  chain_samples <- vector("list", n.chains)
  set.seed(NULL)
  
  for (chain in 1:n.chains) {
    
    N_age.init.chain   <- N_age.init * sample(1:3, 1)   # 1x/2x/3x perturbation
    phi.init.chain     <- matrix(rnorm(ntrt*3, 0, 1), ntrt, 3)
    phi.init.chain[,1] <- 0
    
    Cmodel$N_age <- N_age.init.chain   # N derived from this automatically
    Cmodel$beta0 <- rnorm(1, 0, 1)
    Cmodel$phi   <- phi.init.chain
    Cmodel$mu.p  <- rnorm(1, 0, 1)
    
    Cmcmc$run(n.iter, reset = TRUE)
    chain_samples[[chain]] <- as.matrix(Cmcmc$mvSamples)
  }
  
  # combine into mcmc.list + Gelman-Rubin (psrf ~1.0 = converged, >1.1 = not)
  post.burn <- (n.burn/2 + 1):nrow(chain_samples[[1]])
  a <- mcmc.list(
    mcmc(chain_samples[[1]][post.burn, ]),
    mcmc(chain_samples[[2]][post.burn, ]),
    mcmc(chain_samples[[3]][post.burn, ])
  )
  

  ## CHECK
  gelman <- gelman.diag(a, multivariate = FALSE)
  print(gelman)   # look for psrf column ~1.0; if > 1.1 = red flag

  # traceplots for the diagnostic Ntotal node + one pi_age/p_age param -
  plot(a[, "Ntotal"])
  plot(a[, "p_age[3]"])     # adult detection - the one that collapsed to ~0 before

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

  # recovered pi_age vs. true pi_age, by treatment
  for (t in 1:ntrt) {
    cat("\n--", treatment_levels[t], "--\n")
    for (a in 1:3) {
      col <- paste0("pi_age[", t, ", ", a, "]")
      est <- mean(samples[, col])
      ci  <- quantile(samples[, col], probs = c(0.025, 0.975))
      cat(sprintf("  age %d: true = %.2f, est = %.2f, 95%% CI = [%.2f, %.2f]\n",
                   a, pi_age_true[t, a], est, ci[1], ci[2]))
    }
  }

  # recovered p_age vs. true p_age
  for (a in 1:3) {
    col <- paste0("p_age[", a, "]")
    est <- mean(samples[, col])
    ci  <- quantile(samples[, col], probs = c(0.025, 0.975))
    cat(sprintf("p_age[%d]: true = %.2f, est = %.2f, 95%% CI = [%.2f, %.2f]\n",
                 a, p_age_true[a], est, ci[1], ci[2]))
  }


## ------------------------------------------------------------
## 5. REPEAT SIM ACROSS DETECTION SCENARIOS
## ------------------------------------------------------------

  # Goal: re-simulate + re-fit many times under three 
  # detection scenarios (no gap / moderate gap / large gap),
  # holding true pi_age (trt differences) fixed, and see:
  #   (a) does the model recover the true age composition regardless of
  #       how big the detection gap is?
  #   (b) does bias get worse in the sparser treatments (like BS)?

  detection_scenarios <- list(
    none     = c(J = 0.40, SA = 0.40, A = 0.40),   # no detection difference at all
    moderate = c(J = 0.25, SA = 0.35, A = 0.50),   # juveniles moderately harder to see
    large    = c(J = 0.10, SA = 0.30, A = 0.55)    # juveniles much harder to see
  )

  n_reps <- 20   # start small (5-10) to check timing, then increase

  # storage: one results table per scenario, rows = reps,
  # columns = pi_age[trt,age] flattened + p_age
  pi_age_names <- paste0(rep(treatment_levels, each = 3), "_", rep(c("J","SA","A"), ntrt))
  p_age_names  <- c("J","SA","A")

  results_all <- vector("list", length(detection_scenarios))
  names(results_all) <- names(detection_scenarios)

  for (scen in names(detection_scenarios)) {

    p_age_scenario <- detection_scenarios[[scen]]

    pi_age_results <- matrix(NA, nrow = n_reps, ncol = length(pi_age_names),
                              dimnames = list(NULL, pi_age_names))
    p_age_results  <- matrix(NA, nrow = n_reps, ncol = 3,
                              dimnames = list(NULL, p_age_names))

    for (rep in 1:n_reps) {

      ## re-simulate world (treatment assignment held fixed - only N,
      ## N_age, and y are re-drawn each rep)
      N <- matrix(ncol = I)
      for (i in 1:I) {
        N[i] <- rpois(1, lambda[i])   # lambda[] from section 1, unchanged
      }

      N_age <- matrix(0, nrow = I, ncol = 3)
      for (i in 1:I) {
        N_age[i, ] <- rmultinom(1, N[i], pi_age_true[trt[i], ])
      }

      y <- array(0, dim = c(I, K, 3))
      for (i in 1:I) {
        for (a in 1:3) {
          remaining <- N_age[i, a]
          for (k in 1:K) {
            caught     <- rbinom(1, remaining, p_age_scenario[a])
            y[i, k, a] <- caught
            remaining  <- remaining - caught
          }
        }
      }

      ## refresh data + latent starting values on the ALREADY COMPILED model,
      ## then re-run the chain from scratch (reset = TRUE)
      N_age.init <- apply(y, c(1,3), sum) + 2
      N.init     <- rowSums(N_age.init)

      Cmodel$setData(list(y = y))
      #Cmodel$N     <- N.init
      Cmodel$N_age <- N_age.init

      fit_ok <- tryCatch({
        Cmcmc$run(n.iter, reset = TRUE)
        TRUE
      }, error = function(e) FALSE)

      if (!fit_ok) next

      samples_rep <- as.matrix(Cmcmc$mvSamples)
      samples_rep <- samples_rep[(n.burn/2 + 1):nrow(samples_rep), ]

      # store posterior means
      for (t in 1:ntrt) {
        for (a in 1:3) {
          col <- paste0("pi_age[", t, ", ", a, "]")
          pi_age_results[rep, paste0(treatment_levels[t], "_", c("J","SA","A")[a])] <-
            mean(samples_rep[, col])
        }
      }
      for (a in 1:3) {
        col <- paste0("p_age[", a, "]")
        p_age_results[rep, c("J","SA","A")[a]] <- mean(samples_rep[, col])
      }
    }

    results_all[[scen]] <- list(pi_age = pi_age_results, p_age = p_age_results)
    cat("finished scenario:", scen, "\n")
  }


## ------------------------------------------------------------
## 6. COMPARE SCENARIOS
## ------------------------------------------------------------

  # bias = mean(estimate) - truth, per treatment/age, per scenario
  # this tells you directly whether a bigger detection gap pushes the
  # age-composition estimate away from the true value, and whether
  # your sparsest treatment (BS) is hit hardest

  for (scen in names(detection_scenarios)) {

    cat("\n=====", scen, "detection gap =====\n")
    pi_res <- results_all[[scen]]$pi_age

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

  # boxplots: one figure per scenario, true value marked in red,
  # so you can see at a glance whether boxes are centered on the red dots
  # and whether that gets worse as the detection gap widens
  for (scen in names(detection_scenarios)) {
    pi_res <- results_all[[scen]]$pi_age
    true_vec <- as.vector(t(pi_age_true))   # match column order above
    boxplot(pi_res, main = paste("pi_age recovery -", scen, "detection gap"),
            las = 2, ylab = "estimated proportion")
    points(1:length(pi_age_names), true_vec, col = "red", pch = 19)
  }


#### Notes --------------------------------------------------------------------

  # (1) does bias in pi_age increase from "none" -> "moderate" -> "large" gap?
  #     yes - juvi's bias low and adults bias high under large gap. the other two
  #     scenarios are decent recovery.
  
  # (2) is BS (sparsest treatment) more biased / wider CIs than UU/BU at the
  #     same detection gap?
  #     no - juvi proportion was already low (0.10), so maybe there's just less 
  #     room for the detection gap to pull it down further
  
  # (3) does the model still correctly show BS and UU as DIFFERENT from each
  #     other (the actual treatment-effect question), even under a large gap,
  #     or does the gap wash out the true treatment difference?
  #     yes - direction is preserved but the magnitude is compressed. so we could 
  #     say that BS has proportionally fewer juvi than UU, but if there really 
  #     is a large gap in detection, then we cant trust the effect sizes
  
  
  
  
  
  
  
