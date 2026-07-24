##
## Multi scale occupancy simulation from Josh Twining
## Copied and edited by JW
## 04-09-2025
##

# First section simulates environment and observations
# Second section models occupancy using simulated observations as input data
# Third section assesses convergence and extracts outputs
# Goal: see if we can recreate the simulated parameter values by modelling the simulated data

# Load packages
  library(nimble)
  library(coda)
  library(hesim)
  library(mltools)


#### Simulate ----

# I = sites
I <- 127
# J = plots within sites
J <- 7
# K = occasions
K <- 3


# simulate some covariates
canopycover <- rnorm(I) # simulates canopy cover values at I sites from a normal distribution with mean=0 and sd =1 (scaled)
downedwood <-array(runif(n=I*J, -2, 2), dim =c(I,J)) # simulates dwd values at I sites at J plots from a uniform distribution from -2 to 2

# # # # skip this bit see next section
# this is a hacky way of doing the dummy variables for the categorical, shouldn't impact anything in terms of how recoverable parameters are
# burnt <- rbinom(I, 1, 0.2)  # simulates whether site i was burned, probability set to 0.2
# logged <- rbinom(I, 1, 0.2) # simulates whether site i was logged, probability set to 0.2
# burntandlogged <- rbinom(I, 1, 0.2) # simulates whether site i was logged and then burned, probability set to 0.2
# loggedandburnt <- rbinom(I, 1, 0.2)# simulates whether site i was burned and then logged , probability set to 0.2

# this is the non-hacky way of doing the treatment using a categorical distribution
# install and load the package 'hesim' and 'mltools' for the required functions
# create probabilities that the different treatments occur at
p <- c(0.2, 0.2, 0.2, 0.2, 0.2) #<- set these probabilities to be the proportions the treatments occur in your real dataset
# create matrix for each site to use in the categorical distribution
pmat <- matrix(rep(p, I), nrow = I, ncol = length(p), byrow=TRUE)

# use categorical distribution to assign treatments to each site
treatment <- rcat(I, pmat)

# make dataframe
treatment <- data.frame(treatment)

# create dummy variables for treatment categories
treatment$burnt[treatment$treatment==1] <- 1
treatment$burnt[is.na(treatment$burnt)] <- 0
treatment$logged[treatment$treatment==2] <- 1
treatment$logged[is.na(treatment$logged)] <- 0
treatment$burntandlogged[treatment$treatment==3] <- 1
treatment$burntandlogged[is.na(treatment$burntandlogged)] <- 0
treatment$loggedandburnt[treatment$treatment==4] <- 1
treatment$loggedandburnt[is.na(treatment$loggedandburnt)] <- 0
treatment$control[treatment$treatment==5] <- 1
treatment$control[is.na(treatment$control)] <- 0

# name for the model
burnt <- treatment$burnt
logged <- treatment$logged
burntandlogged <- treatment$burntandlogged
loggedandburnt <- treatment$loggedandburnt

jasmineeffect <- rbinom(n=I*J*K, 1, 0.5)# simulates whether the observer was a pro salamander finder at site i plot j occasion k

dim(jasmineeffect) <- c(I,J,K) #set dimensions to be 3D inline with data sites (rows) x plots (columns) x occasions (slices) 

# look at the different slices
jasmineeffect[,,1]
jasmineeffect[,,2]
jasmineeffect[,,3]


# set coefficents for parameters (their impact on the world we are simulating)
# (on log-odds scale; use plogis(x) to see value on probability-scale, e.g. for beta0.psi/mean occupancy)
beta0.psi <- 1 # mean occupancy (the intercept)
beta1.psi <- 1 # slope term for canopy cover on site occupancy (positive)
beta2.psi <- -1 # slope term for burnt on site occupancy (negative)
beta3.psi <- -1.5 # slope term for logging on site occupancy (stronger negative effect)
beta4.psi <- -1.5 # slope term for burnt and logged on site occupancy 
beta5.psi <- -1.5 # slope term for logged and then burnt on site occupancy 
beta0.theta <- 0 # mean plot use (the intercept for plot usage)
sd.theta <- runif(1, 0.5, 1.5) # standard deviation for a random effect fit on plot use
beta1.theta <- 1 # slope term for downed wood on plot use
alpha0 <- 0 # mean detection probability (the intercept for detection)
alpha1 <- 1.5 # the slope term for jasmineeffect on detection probability


# calculate psi (occupancy probability) of each site using a logit-link
psi <- matrix(ncol = I)  
for (i in 1:I){  
  psi[i] <- plogis(beta0.psi + beta2.psi * burnt[i] + beta3.psi * logged[i] + 
                     beta4.psi * burntandlogged[i] + beta5.psi * loggedandburnt[i])
                  #+ beta1.psi * canopycover[i] 
}

range(psi)

# calculate true z states 
# (this is whether a salamander occurs at a site or not, and is the result of a bernoulli/binomial trail of probability psi)
z <- matrix(ncol = I)
for (i in 1:I){
  z[i] <- rbinom(1, 1, psi[i])
}

sum(z)  # (my data = 56 ENES and 77 OSS)

# mean centered random effect for plot use
plot.theta <- matrix(ncol=I)
for (i in 1:I){
  plot.theta[i] <- rnorm(1, beta0.theta, sd=sd.theta)
}

# calculate use of plots (theta) within sites 
# theta is our probability of use of plot j in site i
theta <- matrix(ncol = J, nrow=I)
for (i in 1:I){
  for (j in 1:J){
    theta[i,j] <- plogis(plot.theta[i]) #+ beta1.theta * downedwood[i,j]) 
  }
}

# calculate true use state (bernoulli/binomial trial of probability omega [w], conditional on z = 1, or site is occupied)
w <- matrix(ncol = J, nrow = I)
for (i in 1:I){
  for (j in 1:J){
    w[i,j] <- (rbinom(1,1, theta[i,j])*z[i])
  }
}

sum(w) # = 251, which means 251 plots are being used?

# calculate detection probability (p)
#3D matrix (i=site, j=plot, k=occasion)
p = array(0, dim = c(I, J, K))
for (i in 1:I){
  for (j in 1:J){
    for (k in 1:K){
      p[i,j,k] <- plogis(alpha0)# + alpha1 * jasmineeffect[i,j,k]) 
    }
  }
}

# simulate observations (sampling) 
# at site i, plot j, and occasion k, as a bernoulli/binomial trial of probability p, conditional on use of the plot 
y = array(0, dim = c(I, J, K))
for (i in 1:I){
  for (j in 1:J){
    for (k in 1:K){
      y[i,j,k] <- (rbinom(1, 1, p[i,j,k])*w[i,j])
    }
  }
}

sum(y[,,1]) # sum = 160
sum(y[,,2]) # sum = 160
y[,,3] # sum = 152

# sum of oss plot dets = 163, sum of enes plot dets = 107

# check out structure of observation data y
str(y)




#### Model ----


# run the chains
n.chains = 3
chains = vector("list", n.chains)
for(chain in 1:n.chains){
  
  #fit model
  constants <- list(I=I,J=J,K=K, logged = logged, burnt = burnt,
                    burntandlogged = burntandlogged, loggedandburnt = loggedandburnt) 
                    # dwd=downedwood, canopycover=canopycover, jasmineeffect=jasmineeffect)
  Nimdata <- list(y=y)
  
  # set initial values
  Niminits <- list(beta0.psi = 0, beta2.psi = 0, beta3.psi = 0, beta4.psi = 0, beta5.psi = 0, 
                   sd.theta=0.5, beta0.theta = 0, alpha0=0)
                  #beta1.psi = 0, alpha1 = 0, beta1.theta=0)
  
  # data summaries
  w.data <- 1*(apply(y,c(1,2),sum)>0)
  z.data <- 1*(rowSums(w.data)>0)
  w.data[w.data==0] <- NA
  z.data[z.data==0] <- NA
  
  # provide data for model
  Nimdata <- list(y=y,z=z.data,w=w.data)
  
  #set parameters to monitor
  parameters <- c("beta0.psi", "beta2.psi", "beta3.psi", "beta4.psi", "beta5.psi",
                  "sd.theta", "beta0.theta", "alpha0", 'zsum')
                  #"beta1.psi", "alpha1", "beta1.theta")
  # (zsum = number of subunits occupied at a site, given that the species is present at that site)
  
  # multi-scale occupancy model for estimating the occupancy state of site i, the use of plot j
  NimModel <- nimbleCode({
    
    # priors
    beta0.psi ~ dlogis(0,1)
    beta0.theta ~ dlogis(0,1)
    sd.theta ~ dunif(0,5)
    alpha0 ~ dlogis(0,1)
    #beta1.psi ~ dnorm(0, sd = 5)
    beta2.psi ~ dnorm(0, sd = 5)
    beta3.psi ~ dnorm(0, sd = 5)
    beta4.psi ~ dnorm(0, sd = 5)
    beta5.psi ~ dnorm(0, sd = 5)
    #beta1.theta ~ dnorm(0, sd = 5)
    #alpha1 ~ dnorm(0, sd = 5)
    
    # likelihood for state model
    for(i in 1:I){
      
      # estimate psi as function  of covs
      logit(psi[i]) <- beta0.psi + beta2.psi * burnt[i] + beta3.psi * logged[i] + 
        beta4.psi * burntandlogged[i] + beta5.psi * loggedandburnt[i] #+ beta1.psi * canopycover[i] 
      z[i] ~ dbern(psi[i]) #is site occupied? z = 1 is yes, z = 0 is no
      
      #centered Random Effect
      plot.theta[i] ~ dnorm(beta0.theta,sd=sd.theta)
      
      # estimate theta (plot usage) as function of covs
      for(j in 1:J){
        logit(theta[i,j]) <- plot.theta[i] #+ beta1.theta * dwd[i,j]
        w[i,j] ~ dbern(theta[i,j]*z[i]) #is plot used given the site occupied?
        
        # likelihood for observation model
        for(k in 1:K){
          logit(p[i,j,k]) <- alpha0 #+ alpha1 * jasmineeffect[i,j,k]
          y[i,j,k] ~ dbern(p[i,j,k]*w[i,j]) # observations as a function of detection probability conditional on the plot is used (w = 1)
        }
      }
    }
    zsum <- sum(z[1:I])
  })# end model
  
  
# Build the model, configure the mcmc, and compileConfigure
  start.time <- Sys.time()
  Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
  conf <- configureMCMC(Rmodel,monitors=parameters, thin=5, useConjugacy=FALSE) # this sets thinning interval to every 5th value
  
  # Build and compile
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  # Run the model
  start.time2 <- Sys.time()
  Cmcmc$run(100000,reset=FALSE) #Can keep extending the run by rerunning this line # this sets n.interations
  end.time <- Sys.time()
  end.time - start.time  # total time for compilation, replacing samplers, and fitting
  end.time - start.time2 # post-compilation run time
  
  burnin <- 10000
  
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



#### Output ----

  # check out traceplots
  plot(a)
  summary(a)
  
  # gelman rubin diags
  gelman <- gelman.diag(a) 
  gelman
  
  # combine the chains
  a=runjags::combine.mcmc(a)
  
  # extract point estimates and %95 CIs
  # can also get these from the summary code above
  
  mean(a[,"beta0.psi"]) # mean occupancy (intercept)
  quantile(a[,"beta0.psi"], probs = c(2.5, 97.5)/100)
  
  mean(a[,"beta1.psi"]) # effect of canopy cover
  quantile(a[,"beta1.psi"], probs = c(2.5, 97.5)/100)
  
  mean(a[,"beta2.psi"]) # burnt
  quantile(a[,"beta2.psi"], probs = c(2.5, 97.5)/100)
  
  mean(a[,"beta3.psi"]) # logged
  quantile(a[,"beta3.psi"], probs = c(2.5, 97.5)/100)
  
  mean(a[,"beta4.psi"]) # burntandlogged
  quantile(a[,"beta4.psi"], probs = c(2.5, 97.5)/100)
  
  mean(a[,"beta5.psi"]) # loggedandburnt
  quantile(a[,"beta5.psi"], probs = c(2.5, 97.5)/100)
  
  mean(a[,"alpha0"]) # mean det prob
  quantile(a[,"alpha0"], probs = c(2.5, 97.5)/100)
  
  mean(a[,"alpha1"]) # slope term for jasmineeffect on detection prob
  quantile(a[,"alpha1"], probs = c(2.5, 97.5)/100)
  
  # bias is calcualted as estimate - truth / truth
  (mean(a[,"alpha0"])-alpha0)
  (mean(a[,"alpha1"])-alpha1)/alpha1
  (mean(a[,"beta0.psi"])-beta0.psi)/beta0.psi


#save.image('salamander_test_run.RData')


  # >   summary(a)
  # 
  # Iterations = 1:10001
  # Thinning interval = 1 
  # Number of chains = 3 
  # Sample size per chain = 10001 
  # 
  # 1. Empirical mean and standard deviation for each variable,
  # plus standard error of the mean:
  #   
  #   Mean      SD  Naive SE Time-series SE
  # alpha0      -0.15175  0.1148 0.0006629       0.001947
  # beta0.psi    1.49659  0.8577 0.0049520       0.040998
  # beta0.theta -0.36909  0.5502 0.0031762       0.036348
  # beta2.psi   -1.96047  1.1716 0.0067638       0.042003
  # beta3.psi    0.03712  2.1776 0.0125719       0.101009
  # beta4.psi   -2.01749  0.9649 0.0055704       0.031870
  # beta5.psi   -2.25376  0.9780 0.0056462       0.031893
  # sd.theta     1.65385  0.6359 0.0036711       0.043851
  # zsum        67.09822 11.4260 0.0659648       0.875092
  # 
  # 2. Quantiles for each variable:
  #   
  #   2.5%     25%     50%      75%    97.5%
  # alpha0      -0.3777 -0.2292 -0.1505 -0.07355  0.07016
  # beta0.psi    0.2518  0.9057  1.3437  1.92000  3.67262
  # beta0.theta -1.6941 -0.6626 -0.2516  0.01998  0.42004
  # beta2.psi   -4.0719 -2.5582 -1.9505 -1.40995 -0.16335
  # beta3.psi   -2.3424 -1.0521 -0.4621  0.23417  7.00756
  # beta4.psi   -4.0304 -2.5300 -1.9643 -1.45887 -0.48902
  # beta5.psi   -4.3465 -2.7823 -2.1950 -1.66762 -0.62608
  # sd.theta     0.6976  1.1820  1.5524  2.04370  3.14502
  # zsum        54.0000 58.0000 64.0000 73.00000 96.00000
  
  