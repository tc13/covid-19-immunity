### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SEIRS discrete time model, gamma (Erlang) distributed waiting times

#Clear R environment
remove(list = ls())

#libraries
require(testthat)

#parameters
sigma_recip <- 4.5          #Average duration of latent period (days)     [From He, X., et al. Temporal dynamics in viral shedding and transmissibility of COVID-19. Nat Med (2020). https://doi.org/10.1038/s41591-020-0869-5]
gamma_recip <- 3.07         #Average duration of infectiousness (days)    [From He, X., et al.]
R0 <- 2.8                   #Basic reproduction number (in the absense of interventions) [From Petra]
beta <- R0*(1/gamma_recip)  #Transmission parameter
omega_recip <- 90           #Average duration of immunity (days) [user specified]
m <- 4                      #Shape paremeter, latent period
n <- 2                      #Shape parameter, infectious period
o <- 2                      #Shape parameter, immune period
dt <- 0.25                  #time period

#vectors
pars <- c(beta=beta, sigma=1/sigma_recip, gamma=1/gamma_recip, omega=1/omega_recip, m=m, n=n, o=o)
time <- seq(1, 300, by=dt)

# State variables
yini <- c(S=0.5, E1=0, E2=0, E3=0, E4=0, I1=0.5, I2=0, R1=0, R2=0) #Initial population size
S = E1 = E2 = E3 = E4 = I1 = I2 = R1 = R2 = N = numeric(length(time))

#Pars
beta = pars["beta"]
sigma = pars["sigma"]
gamma = pars["gamma"]
omega = pars["omega"]
m = pars["m"]
n = pars["m"]
o = pars["o"]

S[1] = yini["S"]
E1[1] = yini["E1"]
E2[1] = yini["E2"]
E3[1] = yini["E3"]
E4[1] = yini["E4"]
I1[1] = yini["I1"]
I2[1] = yini["I2"]
R1[1] = yini["R1"]
R2[1] = yini["R2"]

N[1] = S[1] + E1[1] + E2[1] + E3[1] + E4[1] + I1[1] + I2[1] + R1[1] + R2[1]

for(i in 1:(length(time)-1)){
  t = time[i]
  
  S[(i+1)] = S[i] + omega*dt*R2[i]*o - beta*dt*S[i]*(I1[i]+I2[i])
  E1[(i+1)] = E1[i] + beta*dt*S[i]*(I1[i]+I2[i]) - E1[i]*sigma*dt*m
  E2[(i+1)] = E2[i] + sigma*dt*m*(E1[i] - E2[i])
  E3[(i+1)] = E3[i] + sigma*dt*m*(E2[i] - E3[i])
  E4[(i+1)] = E4[i] + sigma*dt*m*(E3[i] - E4[i])
  I1[(i+1)] = I1[i] + E4[i]*sigma*dt*m - I1[i]*gamma*dt*n
  I2[(i+1)] = I2[i] + gamma*dt*n*(I1[i]-I2[i]) 
  R1[(i+1)] = R1[i] + I2[i]*gamma*dt*n - omega*dt*R1[i]*o
  R2[(i+1)] = R2[i] + omega*dt*o*(R1[i] - R2[i])
  
  #Check population size (N) = 1
  N[(i+1)] = S[(i+1)] + E1[(i+1)] + E2[(i+1)] + E3[(i+1)] + E4[(i+1)] + I1[(i+1)] + I2[(i+1)] + R1[(i+1)] + R2[(i+1)]
  test_that("Pop size (N) = 1", expect_equal(N[(i+1)], 1))
}

#Store as data.frame
out <- data.frame(time=time, S=S, E=E1+E2+E3+E4, I=I1+I2, R=R1+R2)
#plot
with(out, plot(S~time, type="l"))
with(out, plot(E~time, type="l"))
with(out, plot(I~time, type="l"))
with(out, plot(R~time, type="l"))

     