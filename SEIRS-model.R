### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SEIRS: differential equation model
### exponentially distributed infectious period

#Clear R environment
#remove(list = ls())

#libraries
require(deSolve)
require(testthat)

##########################
## SEIRS model function ## 
##########################
SEIRS <- function(t, y, params){
  
  #Check population size (N) = 1
  test_that("Initial Pop size (N) = 1", expect_equal(sum(y), 1))
  #Check parameter values are positive
  test_that("Positive param values", expect_true(all(params>=0)))
  
  #state variables
  S = y[1]
  E = y[2]
  I = y[3]
  R = y[4]
  
  #parameters
  beta = params[1]
  sigma = params[2]
  gamma = params[3]
  omega = params[4]
  
  #Differential equations (susceptible, exposed, infected, recovered)
  dSdt = omega*R - beta*S*I
  dEdt = beta*S*I - sigma*E
  dIdt = sigma*E - gamma*I
  dRdt = gamma*I - omega*R
  
  #return result as a list
  dxdt = c(dSdt, dEdt, dIdt, dRdt)
  
  #Check rates of change sum to zero (as N is unchanged)
  test_that("Net rates of change = 0", expect_equal(sum(dxdt), 0))
  
  list(dxdt)
}
