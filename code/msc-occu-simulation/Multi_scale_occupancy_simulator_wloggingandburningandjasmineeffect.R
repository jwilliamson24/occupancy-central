# I = sites
I <- 127
# J = plots within sites
J <- 7
# K = occasions
K <- 3


# simulate some covariates
canopycover <- rnorm(I) # simulates canopy cover values at I sites from a normal distribution with mean=0 and sd =1 (scaled)

downedwood <-array(runif(n=I*J, -2, 2), dim =c(I,J)) # simulates downed wood values at I sites at J plots from a uniform distribution from -2 to 2

# this is a hacky way of doing the dummy variables for the categorical, it shouldn't impact anything in terms of how recoverable parameters are
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

# set coefficents for parameters (on logs-odd scale, use plogis(x) to see value on probability-scale, e.g. for beta0.psi/mean occupancy)
beta0.psi <- 1 #is mean occupancy (the intercept)
beta1.psi <- 1 # this is the slope term for canopy cover on site occupancy (positive)
beta2.psi <- -1 # this is the slope term for burnt on site occupancy (negative)
beta3.psi <- -1.5 # this is the slope term for logging on site occupancy (negative)
beta4.psi <- -1.5 # this is the slope term for burnt and logged on site occupancy (set it to what ever you like)
beta5.psi <- -1.5 # this is the slope term for logged and then burnt on site occupancy (set it to what ever you like)
beta0.theta <- 0 #this is mean plot use (the intercept for plot usage)
sd.theta <- runif(1, 0.5, 1.5) #this is the standard deviation for a random effect fit on plot use
beta1.theta <- 1 #this is the slope term for downed wood on plot use
alpha0 <- 0 # this is mean detection probability (the intercept for detection)
alpha1 <- 1.5 #this is the slope term for jasmineeffect on detection probability


# calculate psi (occupancy probability) of each site using a logit-link
psi <- matrix(ncol = I)  
for (i in 1:I){  
  psi[i] <- plogis(beta0.psi + beta1.psi * canopycover[i] + beta2.psi * burnt[i] + beta3.psi * logged[i] + beta4.psi * burntandlogged[i])
}

range(psi)

# calculate true z states (this is wheter a salamander occurs at a site or not, and is the result of a bernoulli/binomial trail of probability psi)
z <- matrix(ncol = I)
for (i in 1:I){
  z[i] <- rbinom(1, 1, psi[i])
}

sum(z)

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

# calculate true use state (bernoulli/binomial trail of probability omega [w], conditional on z = 1, or site is occupied)
w <- matrix(ncol = J, nrow = I)
for (i in 1:I){
  for (j in 1:J){
    w[i,j] <- (rbinom(1,1, theta[i,j])*z[i])
  }
}

sum(w)

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

# simulate observations (sampling) at site i, plot j, and occasion k, as a bernoulli/binomial trail of probability p, conditional on use of the plot 
y = array(0, dim = c(I, J, K))
for (i in 1:I){
  for (j in 1:J){
    for (k in 1:K){
      y[i,j,k] <- (rbinom(1, 1, p[i,j,k])*w[i,j])
    }
  }
}
y[,,1]
y[,,2]
y[,,3]

# check out structure of observation data y
str(y)


# load nimble
library(nimble)
library(coda)

# run the chains
n.chains = 3
chains = vector("list", n.chains)
for(chain in 1:n.chains){
  
  #fit model
  constants <- list(I=I,J=J,K=K,canopycover=canopycover, dwd=downedwood, logged = logged, burnt = burnt, burntandlogged = burntandlogged, jasmineeffect=jasmineeffect)
  Nimdata <- list(y=y)
  
  # set initial values
  Niminits <- list(beta0.psi = 0, beta0.theta=0,alpha0=0, alpha1 = 0, beta1.psi = 0, beta2.psi = 0, beta3.psi = 0, beta4.psi = 0, beta1.theta = 0, sd.theta=0.5)
  
  # data summaries
  w.data <- 1*(apply(y,c(1,2),sum)>0)
  z.data <- 1*(rowSums(w.data)>0)
  w.data[w.data==0] <- NA
  z.data[z.data==0] <- NA
  
  # provide data for model
  Nimdata <- list(y=y,z=z.data,w=w.data)
  
  #set parameters to monitor
  parameters <- c("beta0.psi","beta0.theta","alpha0","alpha1", "beta1.psi", "beta2.psi","beta3.psi", "beta4.psi", "beta1.theta","alpha0" ,"sd.theta", 'zsum')
  
  # multi-scale occupancy model for estimating the occupancy state of site i, the use of plot j
  NimModel <- nimbleCode({
    # priors
    beta0.psi ~ dlogis(0,1)
    beta0.theta ~ dlogis(0,1)
    sd.theta ~ dunif(0,5)
    alpha0 ~ dlogis(0,1)
    beta1.psi ~ dnorm(0, sd = 5)
    beta2.psi ~ dnorm(0, sd = 5)
    beta3.psi ~ dnorm(0, sd = 5)
    beta4.psi ~ dnorm(0, sd = 5)
    beta1.theta ~ dnorm(0, sd = 5)
    alpha1 ~ dnorm(0, sd = 5)
    
    # likelihood for state model
    for(i in 1:I){
      # estimate psi as function  of covs
      logit(psi[i]) <- beta0.psi + beta1.psi * canopycover[i] + beta2.psi * burnt[i] + beta3.psi * logged[i]  + beta4.psi * burntandlogged[i]
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

# check out traceplots
plot(a)
summary(a)

# gelman rubin diags
gelman <- gelman.diag(a)

# combine the chains
a=runjags::combine.mcmc(a)

# extract point estimates and %95 CIs
mean(a[,"alpha0"])
quantile(a[,"alpha0"], probs = c(2.5, 97.5)/100)

mean(a[,"alpha0"])
quantile(a[,"alpha0"], probs = c(2.5, 97.5)/100)

mean(a[,"beta0.psi"])
quantile(a[,"beta0.psi"], probs = c(2.5, 97.5)/100)

mean(a[,"beta2.psi"])
quantile(a[,"beta2.psi"], probs = c(2.5, 97.5)/100)

mean(a[,"beta3.psi"])
quantile(a[,"beta3.psi"], probs = c(2.5, 97.5)/100)

mean(a[,"beta.psi"])
quantile(a[,"beta4.psi"], probs = c(2.5, 97.5)/100)

# bias is calcualted as estimate - truth / truth
(mean(a[,"alpha0"])-alpha0)
(mean(a[,"alpha1"])-alpha1)/alpha1
(mean(a[,"beta0.psi"])-beta0.psi)/beta0.psi


save.image('salamander_test_run.RData')



