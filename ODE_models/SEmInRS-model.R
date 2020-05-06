### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SE^mI^nR^oS: differential equation model
### Gamma (Erlang) distributed infectious, latent and immune periods

#Clear R environment
#remove(list = ls())

#libraries
require(deSolve)
require(testthat)

################################
## SE^mI^nR^oS model function ## 
################################
SEmInRS <- function(t, y, params){
  
  #Check population size (N) = 1
  test_that("Initial Pop size (N) = 1", expect_equal(sum(y), 1))
  #Check parameter values are positive
  test_that("Positive param values", expect_true(all(params>=0)))
  
  #state variables
  S = y[1]
  E1 = y[2]
  E2 = y[3]
  E3 = y[4]
  E4 = y[5]
  I1 = y[6]
  I2 = y[7]
  R1 = y[8]
  R2 = y[9]
  
  #parameters
  beta = params[1]
  sigma = params[2]
  gamma = params[3]
  omega = params[4]
  m = params[5]
  n = params[6]
  o = params[7]
  
  #Differential equations (susceptible, exposed, infected, recovered)
  dSdt = omega*o*R2 - beta*S*(I1+I2)
  dE1dt = beta*S*(I1+I2) - sigma*m*E1
  dE2dt = sigma*m*(E1-E2)
  dE3dt = sigma*m*(E2-E3)
  dE4dt = sigma*m*(E3-E4)
  dI1dt = sigma*m*E4-gamma*n*I1
  dI2dt = gamma*n*(I1-I2)
  dR1dt = gamma*n*I2-omega*o*R1
  dR2dt = omega*o*(R1-R2)
  
  #return result as a list
  dxdt = c(dSdt, dE1dt, dE2dt, dE3dt, dE4dt, dI1dt, dI2dt, dR1dt, dR2dt)
  
  #Check rates of change sum to zero (as N is unchanged)
  test_that("Net rates of change = 0", expect_equal(sum(dxdt), 0))
  
  list(dxdt)
}
