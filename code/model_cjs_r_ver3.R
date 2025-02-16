model{
  # Model description:
  #   This model assumes heterogeneity in survival among life-history groups and among subsidy treatments
  #
  # Parameters:
  #   alpha - mean delta during the survey period (unit is m/per day)
  #   xi - detactability with two-pass electrofishing
  #   pi - survival probability during a capture-recapture interval (duration varies by occasion)
  #   p - daily survival
  #   mu.p - mean monthly survival
  #   sigma.p - sd of monthly survival among capture-recapture intervals
  #   phi - cumulative survival probability   
  #   theta - rate parameter for the dispersal model (Laplace)
  # Latent variables:
  #   zs - latent variable indicating whether an individidual remains in the study area or not (zs = 1, stay; zs = 0, leave)
  #   z - latent variable indicating whether an individidual survives or not (z = 1, survive; z = 0, dead)
  # Data:
  #   Y - capture history
  #   X - capture location history (measured as distance from the downstream end to the midpoint of each subsection)
  #   Nday - capture-recapture interval (unit: day)
  #   M - capture-recapture interval (unit: month)
  
  
# priors ------------------------------------------------------------------
  
  ninfo <- 0.01
  
  alpha ~ dnorm(0, ninfo) # Normal distribution
  xi ~ dbeta(1,1) # Beta distribution
  
  # prior distribution of monthly survival rate
  for(t in 1:(Nt-1)){
    for(j in 1:Ng){
      logit.p[j,t] ~ dnorm(logit.mu.p[j], tau.p)
      logit(p[j,t]) <- logit.p[j,t] # Logit transformation (prediction -> 0 ~ 1)
      pi[j,t] <- exp(Nday[t]*log(p[j,t])) # transform from p to pi (cumulative survival probability)
    }
  }

  ## Hyper parameters
  for(j in 1:Ng){
    logit.mu.p[j] <- logit(mu.p[j]) # survival rate 0~1
    mu.p[j] ~ dbeta(1,1) # Beta distribution
  }
  tau.p ~ dscaled.gamma(2.5, 1) # The scaled gamma distribution
  sigma.p <- sqrt(1/tau.p) # square root
  
  
# variable transformation -------------------------------------------------
  # Set start an opportunity for cohorts
  for(j in 1:Ng){
    phi[j,Fo[j]] <- 1
    for(t in Fo[j]:(Nt-1)){
      phi[j,t+1] <- exp(sum(log(pi[j,Fo[j]:t])))
    }
  }


# spatial CJS -------------------------------------------------------------
  
  for(i in 1:Nind){ # Individual replicate
    zs[i,ObsF[i]] <- 1 # latent variable when for first captured.
    for(t in ObsF[i]:(Nt-1) ){# Temporal replicate
      ## Observation process
      loglik[i,t+1] <- logdensity.bern(Y[i,t+1], nu[i,t+1])
      
      Y[i,t+1] ~ dbern(nu[i,t+1]) # Bernoulli distribution
      nu[i,t+1] <- xi*zs[i,t+1]*z[i,t+1]
      
      ## Survival process
      z[i,t+1] ~ dbern(pi[G[i],t]*z[i,t]) # Bernoulli distribution
      
      ## Dispersal process
      X[i,t+1] ~ ddexp(X[i,t], theta[t])T(,1350) # Laplace distribution (double exponential distribution)
      zs[i,t+1] <- step(X[i,t+1]) # step-wise regression
    }
  }
  
  for(t in 1:(Nt-1)){
    theta[t] <- 1/delta[t] # rate parameter for the dispersal model
    log(delta[t]) <- alpha + log(Nday[t])
  }
  
}
