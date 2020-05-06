### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SEIRS discrete time model, age structed for UK population, gamma (Erlang) distributed waiting times
### Model functions

#libraries
require(testthat)

#Sum matrix function
sum.mat <- function(mat1, mat2, mat3, mat4, c){
  # linear combination of four matrices; called in make.intervention.matrix
  sum <- c[1]*mat1 + c[2]*mat2 + c[3]*mat3 + c[4]*mat4
  return(sum)
}

#Combine and scale contact matrixes
make.intervention.matrix <- function(matrixlistin, # list of contact matrices by setting
                                     type= "all", # type of contact to use  -> "all" or "physical" 
                                     R0 = 2.8, # pre-intervention R0  (default 2.8)
                                     intervention = c(0.8,0.7,0,0.2)){#, # home, work, school, other
  # function to make intervention matrix from the list of contact matrices provided based on the type of contact ("all", "physical") specified in type
  # intervention vector specifying the proportion of contacts (compared to no intervention) by setting (home, work, school, other, in that order) that occur under the intervention
  
  context <- c("home", "work", "school", "other")
  matrix_labels <- apply(expand.grid(type, context), 1, paste, collapse = "_")
  home_mat <- as.matrix(matrixlistin[[matrix_labels[1]]][,-1])
  work_mat <- as.matrix(matrixlistin[[matrix_labels[2]]][,-1])
  school_mat <- as.matrix(matrixlistin[[matrix_labels[3]]][,-1])
  other_mat <- as.matrix(matrixlistin[[matrix_labels[4]]][,-1])
 
   # pre-epidemic matrix
  overall_mat <-  sum.mat(home_mat, work_mat, school_mat, other_mat, c = c(1, 1, 1, 1))
  rownames(overall_mat) <- colnames(home_mat)
  int_mat <- sum.mat(home_mat, work_mat, school_mat, other_mat, c = intervention)
  
  # the matrix has participant ages on x-axis (columns), and contact ages on y-axis (rows)
  rownames(int_mat) <-colnames(home_mat)
  return(int_mat)
  }

#Calculate transmission parameter
get_beta <- function(C, p_age, gamma, R0){ 
  #arguments: matrix with average contacts by age group (C), vector with age proportions, gamma, R0
  
  n = length(p_age)
  M = matrix(data=0, nrow=nrow(C), ncol=ncol(C))
  
  #Convert to probabilities (not sure if necessary - doesn't make a difference to beta)
  for(i in 1:n){
    for(j in 1:n){
      M[i,j] = C[i,j]*p_age[i]/p_age[j]}
    }
  
  #Dominant eigenvalue of M
  eigens = eigen(M)
  dominant = max(abs(Re(eigens$values)))
  
  #Calculate beta from R0 and gamma
  beta = (R0 * gamma) / dominant
  return(beta)
}

#Main model function
SEIRS_age_structure <- function(R0=2.8, latent_mean=4.5, infectious_mean=3.07, immune_mean=90,
                                latent_shape=4, infectious_shape=2, immune_shape=2, dt=1, days=400,
                                C, total_population=66435550, p_age){
#time vector
time <- seq(1, days, dt)

#Transmission probability
beta <- get_beta(C=C, p_age=p_age, gamma=(1/infectious_mean), R0=R0)

#Update parameter values (incoporating gamma shape and dt interval)
sigma <- (1/latent_mean) * latent_shape * dt               #probability of becoming infectious
gamma <- (1/infectious_mean) * infectious_shape * dt       #probability of recovery
omega <- (1/immune_mean) * immune_shape * dt               #probability of immunity waning

#number of age classes
classes <- length(p_age)

#Absolute numbers by age group
N_age <- p_age * total_population

#State array
state_names <- c("S", "E1", "E2", "E3", "E4", "I1", "I2", "R1", "R2")
n_states <- length(state_names)
state <- array(data=0, dim=c(length(time), classes, n_states), #rows are times, columns are age classes, matrixes are disease states
               dimnames = list(as.character(time), colnames(C), state_names))

#State variables at t=0
state[1,,"S"] <- N_age #susceptible absolute numbers

#Initialise infections in 100 working adults in each age classes
state[1,,"S"] <- state[1,,"S"]- 100
state[1,,"I1"] <- 50
state[1,,"I2"] <- 50

#Check initial population size (N) = UK Total Pop
N = sum(state[1,,])
test_that("Initial pop size (N) = 1", expect_equal(N, total_population))

#For each timestep t
for(t in 1:(length(time)-1)){

  #new infections per timestep  
  I <- apply(state[t,,c("I1","I2")], 1, sum) #sum I1 and I2 classes
  lambda = beta * (as.matrix(C) %*% I/N_age) #force of infection
  new_infection = lambda * state[t,,"S"] * dt #new infections per timestep

  #Epidemic transitions at time t+1
  state[(t+1),,"S"] = state[t,,"S"] + omega*state[t,,"R2"] - new_infection
  state[(t+1),,"E1"] = state[t,,"E1"] + new_infection - state[t,,"E1"]*sigma
  state[(t+1),,"E2"] = state[t,,"E2"] + sigma*(state[t,,"E1"]-state[t,,"E2"])
  state[(t+1),,"E3"] = state[t,,"E3"] + sigma*(state[t,,"E2"]-state[t,,"E3"])
  state[(t+1),,"E4"] = state[t,,"E4"] + sigma*(state[t,,"E3"]-state[t,,"E4"])
  state[(t+1),,"I1"] = state[t,,"I1"] + sigma*state[t,,"E4"] - state[t,,"I1"]*gamma
  state[(t+1),,"I2"] = state[t,,"I2"] + gamma*(state[t,,"I1"]-state[t,,"I2"]) 
  state[(t+1),,"R1"] = state[t,,"R1"] + gamma*state[t,,"I2"] - state[t,,"R1"]*omega
  state[(t+1),,"R2"] = state[t,,"R2"] + omega*(state[t,,"R1"] - state[t,,"R2"])

  #Check population size (N) =  UK total pop
  N = sum(state[(t+1),,])
  test_that("Pop size (N) = 1", expect_equal(N, total_population))
  
  #Test that no state value is negative
  test_that("States non-negative", expect_true(all(state[(t+1),,]>=0)))
}

#Combine E, I, R substates into new array
state_names_SEIR <- state_names <- c("S", "E", "I", "R")

out <- array(data=0, dim = c(length(time), classes, length(state_names_SEIR)),
             dimnames = list(as.character(time), colnames(C), state_names_SEIR))
out[,,"S"] <- state[,,"S"]
out[,,"E"] <- state[,,"E1"] + state[,,"E2"] + state[,,"E3"] + state[,, "E4"]
out[,,"I"] <- state[,,"I1"] + state[,,"I2"]
out[,,"R"] <- state[,,"R1"] + state[,,"R2"]

#Return out array
return(out)
}
