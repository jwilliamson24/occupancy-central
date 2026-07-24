##
## Multi scale occupancy simulation from Josh Twining
## Copied and edited by JW
## 04-09-2025
##

# First section simulates environment and observations
# Second section models occupancy using simulated observations as input data
# Third section assesses convergence and extracts outputs
# Goal: see if we can recreate the simulated parameter values by modelling the simulated data


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

# this is a hacky way of doing the dummy variables for the categorical, shouldn't impact anything in terms of how recoverable parameters are
burnt <- rbinom(I, 1, 0.2)  # simulates whether site i was burned, probability set to 0.2
logged <- rbinom(I, 1, 0.2) # simulates whether site i was logged, probability set to 0.2
burntandlogged <- rbinom(I, 1, 0.2) # simulates whether site i was logged and then burned, probability set to 0.2
loggedandburnt <- rbinom(I, 1, 0.2)# simulates whether site i was burned and then logged , probability set to 0.2

jasmineeffect <- rbinom(n=I*J*K, 1, 0.5)# simulates whether the observer was a pro salamander finder at site i plot j occasion k

dim(jasmineeffect) <- c(I,J,K) #set dimensions to be 3D inline with data sites (rows) x plots (columns) x occasions (slices) 

# look at the different slices
jasmineeffect[,,1]
jasmineeffect[,,2]
jasmineeffect[,,3]


# set coefficents for parameters (their impact on the world we are simulating)
# (on logs-odd scale; use plogis(x) to see value on probability-scale, e.g. for beta0.psi/mean occupancy)
beta0.psi <- 1 # mean occupancy (the intercept)
beta1.psi <- 1 # this is the slope term for canopy cover on site occupancy (positive)
beta2.psi <- -0.5 # this is the slope term for burnt on site occupancy (negative)
beta3.psi <- -1 # this is the slope term for logging on site occupancy (stronger negative effect)
beta4.psi <- -1 # this is the slope term for burnt and logged on site occupancy 
beta5.psi <- -1 # this is the slope term for logged and then burnt on site occupancy 
beta0.theta <- 0 # this is mean plot use (the intercept for plot usage)
sd.theta <- runif(1, 0.5, 1.5) # this is the standard deviation for a random effect fit on plot use
beta1.theta <- 1 # this is the slope term for downed wood on plot use
alpha0 <- 0 # this is mean detection probability (the intercept for detection)
alpha1 <- 1 # this is the slope term for jasmineeffect on detection probability


# calculate psi (occupancy probability) of each site using a logit-link
psi <- matrix(ncol = I)  
for (i in 1:I){  
  psi[i] <- plogis(beta0.psi + beta2.psi * burnt[i] + beta3.psi * logged[i] + 
                     beta4.psi * burntandlogged[i] + beta5.psi * loggedandburnt[i]) #+ beta1.psi * canopycover[i] 
}
            
range(psi)

# calculate true z states 
# (this is whether a salamander occurs at a site or not, and is the result of a bernoulli/binomial trail of probability psi)
z <- matrix(ncol = I)
for (i in 1:I){
  z[i] <- rbinom(1, 1, psi[i])
}

sum(z)  # 74 sites occupied, similar to my data (56 ENES and 77 OSS)

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
    theta[i,j] <- plogis(plot.theta[i] + beta1.theta * downedwood[i,j]) 
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
      p[i,j,k] <- plogis(alpha0 + alpha1 * jasmineeffect[i,j,k]) 
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

y[,,1] # sum = 160
y[,,2] # sum = 160
y[,,3] # sum = 152

# sum of oss plot dets = 163, sum of enes plot dets = 107

# check out structure of observation data y
str(y)




#### Model ----

# load nimble
library(nimble)
library(coda)


# run the chains
n.chains = 3
chains = vector("list", n.chains)
for(chain in 1:n.chains){
  
  #fit model
  constants <- list(I=I,J=J,K=K, logged = logged, burnt = burnt, dwd=downedwood,
                    burntandlogged = burntandlogged, loggedandburnt = loggedandburnt, jasmineeffect=jasmineeffect) 
                    #canopycover=canopycover)
  Nimdata <- list(y=y)
  
  # set initial values
  Niminits <- list(beta0.psi = 0, beta2.psi = 0, beta3.psi = 0, beta4.psi = 0, beta5.psi = 0, 
                   sd.theta=0.5, beta0.theta = 0, alpha0=0, alpha1 = 0, beta1.theta=0)
                  #beta1.psi = 0)
  
  # data summaries
  w.data <- 1*(apply(y,c(1,2),sum)>0)
  z.data <- 1*(rowSums(w.data)>0)
  w.data[w.data==0] <- NA
  z.data[z.data==0] <- NA
  
  # provide data for model
  Nimdata <- list(y=y,z=z.data,w=w.data)
  
  #set parameters to monitor
  parameters <- c("beta0.psi", "beta2.psi", "beta3.psi", "beta4.psi", "beta5.psi",
                  "sd.theta", "beta0.theta", "alpha0" , "alpha1", 'zsum', "beta1.theta")
                  #"beta1.psi")
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
    beta1.theta ~ dnorm(0, sd = 5)
    alpha1 ~ dnorm(0, sd = 5)
    
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
        logit(theta[i,j]) <- plot.theta[i] + beta1.theta * dwd[i,j]
        w[i,j] ~ dbern(theta[i,j]*z[i]) #is plot used given the site occupied?
        
        # likelihood for observation model
        for(k in 1:K){
          logit(p[i,j,k]) <- alpha0 + alpha1 * jasmineeffect[i,j,k]
          y[i,j,k] ~ dbern(p[i,j,k]*w[i,j]) # observations as a function of detection probability conditional on the plot is used (w = 1)
        }
      }
    }
    zsum <- sum(z[1:I])
  })# end model
  
  
# Build the model, configure the mcmc, and compileConfigure
  start.time <- Sys.time()
  Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
  conf <- configureMCMC(Rmodel,monitors=parameters, thin=5, useConjugacy=FALSE)
  
  # Build and compile
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  # Run the model
  start.time2 <- Sys.time()
  Cmcmc$run(50000,reset=FALSE) #Can keep extending the run by rerunning this line
  end.time <- Sys.time()
  end.time - start.time  # total time for compilation, replacing samplers, and fitting
  end.time - start.time2 # post-compilation run time
  
  burnin <- 10000
  
  #get the chains
  mvSamples = as.matrix(Cmcmc$mvSamples)
  chains[[chain]]=mvSamples
}

n.iter = 10000
n.burn = 5000


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
  
  # combine the chains
  a=runjags::combine.mcmc(a)
  
  # extract point estimates and %95 CIs
  
  mean(a[,"beta0.psi"]) # mean occupancy (intercept)
  #quantile(a[,"beta0.psi"], probs = c(2.5, 97.5)/100)
  
  mean(a[,"beta1.psi"]) # effect of canopy cover
  #quantile(a[,"beta1.psi"], probs = c(2.5, 97.5)/100)
  
  mean(a[,"beta2.psi"]) # burnt
  #quantile(a[,"beta2.psi"], probs = c(2.5, 97.5)/100)
  
  mean(a[,"beta3.psi"]) # logged
  #quantile(a[,"beta3.psi"], probs = c(2.5, 97.5)/100)
  
  mean(a[,"beta4.psi"]) # burntandlogged
  #quantile(a[,"beta4.psi"], probs = c(2.5, 97.5)/100)
  
  mean(a[,"beta5.psi"]) # loggedandburnt
  #quantile(a[,"beta5.psi"], probs = c(2.5, 97.5)/100)
  
  mean(a[,"alpha0"]) # mean det prob
  #quantile(a[,"alpha0"], probs = c(2.5, 97.5)/100)
  
  mean(a[,"alpha1"]) # slope term for jasmineeffect on detection prob
  #quantile(a[,"alpha1"], probs = c(2.5, 97.5)/100)
  
  # bias is calcualted as estimate - truth / truth
  (mean(a[,"alpha0"])-alpha0)
  (mean(a[,"alpha1"])-alpha1)/alpha1
  (mean(a[,"beta0.psi"])-beta0.psi)/beta0.psi


#save.image('salamander_test_run.RData')


# estimates with 1270 sites...
  # >   mean(a[,"beta0.psi"]) # mean occupancy (intercept)
  # [1] 0.8363193
  #   >   mean(a[,"beta1.psi"]) # effect of canopy cover
  # [1] 0.899209
  #   >   mean(a[,"beta2.psi"]) # burnt
  # [1] -0.6836093
  #   >   mean(a[,"beta3.psi"]) # logged
  # [1] -1.372336
  #   >   mean(a[,"beta4.psi"]) # burntandlogged
  # [1] -1.35614
  #   >   mean(a[,"beta5.psi"]) # loggedandburnt
  # [1] -1.403495
  #   >   mean(a[,"alpha0"]) # mean det prob
  # [1] 0.009684279
  #   >   mean(a[,"alpha1"]) # slope term for jasmineeffect on detection prob
  # [1] 1.417615

  
# estimates with 127 sites  
  # >   mean(a[,"beta0.psi"]) # mean occupancy (intercept)
  # [1] 1.469344
  # >   mean(a[,"beta1.psi"]) # effect of canopy cover
  # [1] 1.190795
  # >   mean(a[,"beta2.psi"]) # burnt
  # [1] -0.4060675
  # >   mean(a[,"beta3.psi"]) # logged
  # [1] -1.25438
  # >   mean(a[,"beta4.psi"]) # burntandlogged
  # [1] -1.816969
  # >   mean(a[,"beta5.psi"]) # loggedandburnt
  # [1] -1.520048
  # >   mean(a[,"alpha0"]) # mean det prob
  # [1] -0.05409016
  # >   mean(a[,"alpha1"]) # slope term for jasmineeffect on detection prob
  # [1] 1.4754
  
  
# what they should be  
  # beta0.psi <- 1 # mean occupancy (the intercept)
  # beta1.psi <- 1 # this is the slope term for canopy cover on site occupancy (positive)
  # beta2.psi <- -1 # this is the slope term for burnt on site occupancy (negative)
  # beta3.psi <- -1.5 # this is the slope term for logging on site occupancy (stronger negative effect)
  # beta4.psi <- -1.5 # this is the slope term for burnt and logged on site occupancy 
  # beta5.psi <- -1.5 # this is the slope term for logged and then burnt on site occupancy 
  # beta0.theta <- 0 # this is mean plot use (the intercept for plot usage)
  # sd.theta <- runif(1, 0.5, 1.5) # this is the standard deviation for a random effect fit on plot use
  # beta1.theta <- 1 # this is the slope term for downed wood on plot use
  # alpha0 <- 0 # this is mean detection probability (the intercept for detection)
  # alpha1 <- 1.5 # this is the slope term for jasmineeffect on detection probability 
  
  
# with 127 sites and changed effect sizes of betas as coded in the current model
  # >   mean(a[,"beta0.psi"]) # mean occupancy (intercept)
  # [1] 1.014734
  # >   mean(a[,"beta1.psi"]) # effect of canopy cover
  # [1] 1.376441
  # >   mean(a[,"beta2.psi"]) # burnt
  # [1] -0.9795046
  # >   mean(a[,"beta3.psi"]) # logged
  # [1] -1.333538
  # >   mean(a[,"beta4.psi"]) # burntandlogged
  # [1] -0.6592855
  # >   mean(a[,"beta5.psi"]) # loggedandburnt
  # [1] -0.798641
  # >   mean(a[,"alpha0"]) # mean det prob
  # [1] 0.1154041
  # >   mean(a[,"alpha1"]) # slope term for jasmineeffect on detection prob
  # [1] 1.150361
   
# taking out canopy cover so the model is smaller
  # >   mean(a[,"beta0.psi"]) # mean occupancy (intercept)
  # [1] 1.088093
  # >   mean(a[,"beta2.psi"]) # burnt
  # [1] -0.302911
  # >   mean(a[,"beta3.psi"]) # logged
  # [1] -0.6437466
  # >   mean(a[,"beta4.psi"]) # burntandlogged
  # [1] -0.8953914
  # >   mean(a[,"beta5.psi"]) # loggedandburnt
  # [1] -1.615628
  # >   mean(a[,"alpha0"]) # mean det prob
  # [1] 0.05032677
  # >   mean(a[,"alpha1"]) # slope term for jasmineeffect on detection prob
  # [1] 0.8938349
  