### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SEIRS discrete time model, exponential distributions

#Clear R environment
remove(list = ls())

#libraries
require(testthat)

#parameters
sigma_recip <- 4.5        #Average duration of latent period (days)     [From He, X., et al. Temporal dynamics in viral shedding and transmissibility of COVID-19. Nat Med (2020). https://doi.org/10.1038/s41591-020-0869-5]
gamma_recip <- 3.07       #Average duration of infectiousness (days)    [From He, X., et al.]
R0 <- 2.8
beta <- R0*(1/gamma_recip)
omega_recip <- 90

#vectors
pars <- c(beta=beta, sigma=1/sigma_recip, gamma=1/gamma_recip, omega=1/omega_recip)
time <- seq(1, 300, 1)

# State variables
yini <- c(S = 0.999, E = 0, I = 0.001, R = 0)
S = E = I = R = N = numeric(length(time))

#Pars
beta = pars["beta"]
sigma = pars["sigma"]
gamma = pars["gamma"]
omega = pars["omega"]

S[1] <- yini["S"]
E[1] <- yini["E"]
I[1] <- yini["I"]
R[1] <- yini["R"]

N[1] = S[1] + E[1] + I[1] + R[1]

for(i in 1:(length(time)-1)){
   t = time[i]
  
  S[(i+1)] = S[i] + omega*R[i] - beta*S[i]*I[i]
  E[(i+1)] = E[i] + beta*S[i]*I[i] - E[i]*sigma
  I[(i+1)] = I[i] + E[i]*sigma - I[i]*gamma
  R[(i+1)] = R[i] + I[i]*gamma - omega*R[i]
  
  #Check population size (N) = 1
  N[(i+1)] = S[(i+1)] + E[(i+1)] + I[(i+1)] + R[(i+1)]
  test_that("Initial Pop size (N) = 1", expect_equal(N[(i+1)], 1))
    
}

#Store as data.frame
out <- data.frame(time=time, S=S, E=E, I=I, R=R)
#plot
with(out, plot(S~time))