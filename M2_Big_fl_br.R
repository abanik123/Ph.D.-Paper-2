library(nimble)
library(MCMCvis)
library(tidyverse)
library(dplyr)


#' Read in the data. 
#------------------------------------------------------------------------#

cap_data <- read.csv("Bigfish_capture_full.csv")
dd_fl <- read.csv("Bigfish_forklength_full.csv")
b_fl <- read.csv("Bigfish_breeding_full_wuk_sex.csv")

y <- cap_data %>%
  select(Fall_2009:Fall_2020) %>%
  as.matrix()
#head(y)

dd_fl <- dd_fl %>%
  select(Fall_2009: Fall_2020) %>%
  as.matrix()
dd_fl[dd_fl == 0] <- NA
#head(dd_fl)

b_fl <- b_fl %>%
  select(Fall_2009:Fall_2020) %>%
  as.matrix()
#head(y)

#Scaling the forklength
#-----------------------------------------------------------------------#

# Define a function to scale each row separately
logscale_row <- function(row) {
  # Filter out NA values
  row_without_nas <- row[!is.na(row)]
  # Scale the row
  logscaled_row <- log(row_without_nas)
  # Replace NA positions with NA in scaled row
  logscaled_row_with_nas <- rep(NA, length(row))
  logscaled_row_with_nas[!is.na(row)] <- logscaled_row
  return(logscaled_row_with_nas)
}

# Apply the scaling function to each row of the forklength
d_fl <- t(apply(dd_fl, 1, logscale_row))

#d_fl <- dd_fl

#head(d_fl)


#' Get occasion of first capture.
#------------------------------------------------------------------------#
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)
#first

# Nimble Code
#------------------------------------------------------------------------#

code_m <- nimbleCode({
  # Priors
  
  alpha[1] ~ dnorm(0,100)
  alpha[2] ~ dnorm(0,100)
  
  beta[1] ~ dnorm(0,100)
  beta[2] ~ dnorm(0,100)
  beta[3] ~ dnorm(0,100)
  
  psi[1] ~ dbeta(1,1)
  psi[2] ~ dbeta(1,1)
  
  #mu_o <- 0.001
  mu_sp_o <- 1
  mu_sl_o <- 1
  #sigma_mu_o <- 100
  sigma_musp_o <- 0.1
  sigma_musl_o <- 0.1
  
  a_sigma <- 0.001
  b_sigma <- 1000
  
  sigma ~ dinvgamma(a_sigma, b_sigma)
  
  for(t in 1:nyears){
    
    mu_sp[t] ~ dgamma(shape = mu_sp_o, scale = sigma_musp_o)
    mu_sl[t] ~ dgamma(shape = mu_sl_o, scale = sigma_musl_o)
    
  }
  
  for (i in 1:N) {
    
    for (t in (first[i]+1):nyears) {
      
      b[i,t] ~ dbern((psi[1]*b[i,t-1])+(psi[2]*(1-b[i,t-1])))
      
      mu[i,t] <- mu_sp[t]*b[i,t] + mu_sl[t]*(1-b[i,t])
      
      fl[i, t] ~ T(dnorm(fl[i, t-1] + mu[i,t], sd = sigma), fl[i,t-1], )
      
    }
    
    # Modelling
    
    # Initial state
    z[i,first[i]] <- y[i,first[i]] - 1
    
    for (t in (first[i]+1):nyears) {
      # State transition
      z[i, t] ~ dcat(px[z[i, t - 1], 1:2, i, t])
      
      # Observation process
      y[i, t] ~ dcat(po[z[i, t], 1:2, i, t])
    }
  }
  
  # probabilities of state z(t+1) given z(t)
  for (i in 1:N) {
    
    for (t in first[i]:nyears) {
      
      logit(phi[i, t]) <- alpha[1] + (alpha[2]*fl[i, t])
      logit(p[i, t]) <- beta[1] + (beta[2]*fl[i,t]) + (beta[3]*fl[i,t]*fl[i,t])
      
      
      # Latent state process probabilities
      px[1, 1, i, t] <- phi[i,t]
      px[1, 2, i, t] <- 1 - phi[i,t]
      
      px[2, 1, i, t] <- 0
      px[2, 2, i, t] <- 1
      
      # Observation process probabilities
      po[1, 1, i, t] <- 1 - p[i,t]
      po[1, 2, i, t] <- p[i,t]
      
      po[2, 1, i, t] <- 1
      po[2, 2, i, t] <- 0
      
    }
  }
})

#Input Data
#-----------------------------------------------------------------------#
my.data = list(y=y + 1,b = b_fl, fl = d_fl)

#Constants in a list. 
#-----------------------------------------------------------------------#
N = nrow(y)
nyears = ncol(y)

my.constants <- list(N = N, nyears = nyears, first = first)

# Generate inits for the latent states
#-----------------------------------------------------------------------#

# Initializing fl

fl_inits_temp <- dd_fl

replace_nas_increasing_order <- function(row) {
  na_indices <- which(is.na(row))
  
  for (i in na_indices) {
    if (i == 1) {
      row[i] <- row[i + 1] + 1
    } else if (i == length(row)) {
      row[i] <- row[i - 1] + 1
    } else {
      # Adjust the range to handle consecutive NAs
      start_range <- ifelse(!is.finite(row[i - 1]), 0, row[i - 1] + 1)
      end_range <- ifelse(!is.finite(row[i + 1]), start_range + 1, row[i + 1] - 1)
      
      if (!is.finite(start_range) || !is.finite(end_range)) {
        row[i] <- row[i - 1] + 1
      } else {
        row[i] <- sample(seq(start_range, end_range), 1)
      }
    }
  }
  
  return(row)
}

fl_inits_1 <- apply(fl_inits_temp, 1, replace_nas_increasing_order)

fl_inits_2 <- t(fl_inits_1)

#head(fl_inits_1)
#head(fl_inits_2)

# Function to replace values in Matrix B with NAs based on NAs in Matrix A
replace_values_with_nas <- function(matrix_a, matrix_b) {
  # Get the indices where Matrix A has NAs
  na_indices <- which(!is.na(matrix_a), arr.ind = TRUE)
  
  # Replace values in Matrix B with NAs at the corresponding positions
  matrix_b[na_indices] <- NA
  
  return(matrix_b)
}

d_fl_inits <- replace_values_with_nas(fl_inits_temp, fl_inits_2)

d_fl_inits <- as.matrix(d_fl_inits)

#head(d_fl_inits)

for (i in 1:N){
  if (first[i] == 1) next
  if (first[i] > 1) d_fl_inits[i,1:(first[i]-1)] <- NA
}

#head(d_fl_inits)

# Apply the scaling function to each row of the initial forklength
fl_inits <- t(apply(d_fl_inits, 1, logscale_row))
#fl_inits <- d_fl_inits
#head(fl_inits)
#tail(fl_inits)
#------------------------------------------------------------------------#
# Initializing alive state

x.init <- matrix(1, N, nyears)  # Set to 1 or another appropriate value for initial state
for (i in 1:N) {
  if (first[i] > 1) x.init[i, 1:(first[i]-1)] <- NA
}

z_inits <- x.init

#------------------------------------------------------------------------#

# Initialize Breeding

b_d <- b_fl
# Find the indices of the first occurrence of 0 or 1 in each row
first_01_indices <- apply(b_d, 1, function(x) which(!is.na(x) & x %in% c(0, 1))[1])

# Replace NAs with 0s and 1s, and 0s and 1s with NAs, keeping NAs before the first 0 or 1 unchanged
b_inits_d <- b_d
for (i in 1:nrow(b_d)) {
  if (!is.na(first_01_indices[i])) {
    b_inits_d[i, 1:(first_01_indices[i] - 1)] <- b_d[i, 1:(first_01_indices[i] - 1)]
    b_inits_d[i, first_01_indices[i]:ncol(b_d)] <- ifelse(is.na(b_d[i, first_01_indices[i]:ncol(b_d)]), sample(c(0, 1), sum(is.na(b_d[i, first_01_indices[i]:ncol(b_d)])), replace = TRUE), NA)
  } else {
    b_inits_d[i, ] <- ifelse(is.na(b_d[i, ]), sample(c(0, 1), sum(is.na(b_d[i, ])), replace = TRUE), NA)
  }
}

b_inits <- b_inits_d

#------------------------------------------------------------------------#

# initializing Drift

mu_sp_inits <- runif(nyears)
mu_sl_inits <- runif(nyears)

#head(mu_inits)

#------------------------ Final initialization --------------------------#

initial.values <- list(alpha = runif(2, 0, 1), beta = runif(3,0,1),
                       psi = runif(2, 0, 1),
                       mu_sp_o = 1, sigma_musp_o = 0.1,
                       mu_sl_o = 1, sigma_musl_o = 0.1,
                       sigma = 1, a_sigma = 0.001, b_sigma = 1000,
                       fl = fl_inits, z = z_inits, b = b_inits,
                       mu_sp = mu_sp_inits, mu_sl = mu_sl_inits)


#' Parameters to monitor. 
#------------------------------------------------------------------------#
parameters.to.save <- c("alpha", "beta", "psi")

#' MCMC settings.
#------------------------------------------------------------------------#
n.iter <- 300000
n.burnin <- 50000
n.chains <- 3

#--- CHECKED NIMBLE MODEL ---#

mcmc.multistate <- nimbleMCMC(code = code_m,
                              constants = my.constants,
                              data = my.data,
                              inits = initial.values,
                              monitors = parameters.to.save,
                              niter = n.iter,
                              nburnin = n.burnin,
                              nchains = n.chains,
                              summary = TRUE, WAIC = TRUE)

waic <- mcmc.multistate$WAIC
waic
#Model Summary
samples<- mcmc.multistate$samples

pdf(file = "Big_fl_mp_m2.pdf")
MCMCplot(samples, HPD = T)
dev.off()

s <- MCMCsummary(samples, round = 5, params = c("alpha", "beta", "psi"))
MCMCtrace(samples,pdf = T,open_pdf = F,filename = "Big_fl_m2", ind = TRUE,
          Rhat = FALSE, n.eff = FALSE)
write.csv(s, file = "Big_fl_m2_sum.csv")


