### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SEIRRS discrete time model, age structed for UK population, gamma distributed waiting times, two immunity classes
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

#Get omega function
get_omega <- function(mean, shape, dt){
  if(mean<=0) omega <- 0
  else omega <- (1/mean)*shape*dt
  return(omega)
}

#Text with numbers function
textNum <- function(text, num){
  v = c()
  for(i in 1:num){
    num.string = as.character(i)
    v[i] = paste(text, num.string, sep="")
  }
  return(v)
}

#Main model function
SEIRRS_age_structure <- function(R0=2.8, latent_mean=4.5, infectious_mean=3.07, immune_mean_1=90,
                                immune_mean_2=180, latent_shape=4, infectious_shape=2, immune_shape=2, 
                                dt=1, days=400, C, total_population=66435550, p_age, p_immunity, I_init){
  #time vector
  time <- seq(1, days, dt)
  
  #Transmission probability
  beta <- get_beta(C=C, p_age=p_age, gamma=(1/infectious_mean), R0=R0)
  
  #Update parameter values (incoporating gamma shape and dt interval)
  sigma <- (1/latent_mean) * latent_shape * dt               #probability of becoming infectious
  gamma <- (1/infectious_mean) * infectious_shape * dt       #probability of recovery
  omega_1 <- get_omega(mean=immune_mean_1, shape=immune_shape, dt=dt)   #Duration of immunity in R class 1
  omega_2 <- get_omega(mean=immune_mean_2, shape=immune_shape, dt=dt)   #Duration of immunity in R class 2
  #number of age classes
  classes <- length(p_age)
  
  #Absolute numbers by age group
  N_age <- p_age * total_population
  
  #State array
  Es <- textNum("E", latent_shape)
  Is <- textNum("I", infectious_shape)
  RAs <- textNum("RA", immune_shape)
  RBs <- textNum("RB", immune_shape)
  
  state_names <- c("S", Es, Is, RAs, RBs)
  n_states <- length(state_names)
  state <- array(data=0, dim=c(length(time), classes, n_states), #rows are times, columns are age classes, matrixes are disease states
                 dimnames = list(as.character(time), colnames(C), state_names))
  
  #State variables at t=0
  state[1,,"S"] <- N_age #susceptible absolute numbers
  
  #Initialise infections by age group with I_init vector
  state[1,,"S"] <- state[1,,"S"]- I_init
  for(j in 1:infectious_shape){
  state[1,,Is[j]] <- I_init/infectious_shape}
  
  #Check initial population size (N) = UK Total Pop
  N = sum(state[1,,])
  test_that("Initial pop size (N) = UK total pop", expect_equal(N, total_population))
  
  #For each timestep t
  for(t in 1:(length(time)-1)){
    
    #new infections per timestep  
    I <- apply(state[t,,Is], 1, sum) #sum I sub-classes
    lambda = beta * (as.matrix(C) %*% I/N_age) #force of infection
    new_infection = lambda * state[t,,"S"] * dt #new infections per timestep
    
    #Epidemic transitions at time t+1
    #Susceptibles
    state[(t+1),,"S"] = state[t,,"S"] + omega_1*state[t,,RAs[immune_shape]] + omega_2*state[t,,RBs[immune_shape]] - new_infection
    #First exposed class
    state[(t+1),,"E1"] = state[t,,"E1"] + new_infection - state[t,,"E1"]*sigma
    #Other exposed classes
    if(latent_shape>1){
      for(k in 2:latent_shape){
    state[(t+1),,Es[k]] = state[t,,Es[k]] + sigma*(state[t,,Es[(k-1)]]-state[t,,Es[k]])
        }
    }
    #first infectious class
    state[(t+1),,"I1"] = state[t,,"I1"] + sigma*state[t,,Es[latent_shape]] - state[t,,"I1"]*gamma
    #Other infectious classes
    if(infectious_shape>1){
      for(k in 2:infectious_shape){
        state[(t+1),,Is[k]] = state[t,,Is[k]] + gamma*(state[t,,Is[(k-1)]]-state[t,,Is[k]])
        }
    }
    #Duration of immunity, first classes
    state[(t+1),,"RA1"] = state[t,,"RA1"] + p_immunity*gamma*state[t,,Is[infectious_shape]] - state[t,,"RA1"]*omega_1
    state[(t+1),,"RB1"] = state[t,,"RB1"] + (1-p_immunity)*gamma*state[t,,Is[infectious_shape]] - state[t,,"RB1"]*omega_2
    if(immune_shape>1){
      for(k in 2:immune_shape){
        state[(t+1),,RAs[k]] = state[t,,RAs[k]] + omega_1*(state[t,,RAs[(k-1)]] - state[t,,RAs[k]])
        state[(t+1),,RBs[k]] = state[t,,RBs[k]] + omega_2*(state[t,,RBs[(k-1)]] - state[t,,RBs[k]])
      }
    }
    
    #Check population size (N) =  UK total pop
    N = sum(state[(t+1),,])
    test_that("Pop size (N) = Pop size", expect_equal(N, total_population))
    
    #Test that no state value is negative
    test_that("States non-negative", expect_true(all(state[(t+1),,]>=0)))
  }
  
  #Combine E, I, R substates into new array
  state_names_SEIR <- state_names <- c("S", "E", "I", "R1", "R2")
  
  out <- array(data=0, dim = c(length(time), classes, length(state_names_SEIR)),
               dimnames = list(as.character(time), colnames(C), state_names_SEIR))
  out[,,"S"] <- state[,,"S"]
  for(i in 1:latent_shape){out[,,"E"] <- out[,,"E"] + state[,,Es[i]]}
  for(i in 1:infectious_shape){out[,,"I"] <- out[,,"I"] + state[,,Is[i]]}  
  for(i in 1:immune_shape){out[,,"R1"] <- out[,,"R1"] + state[,,RAs[i]]}  
  for(i in 1:immune_shape){out[,,"R2"] <- out[,,"R2"] + state[,,RBs[i]]}  
  
  #Check population size (N) =  UK total pop
  N = sum(out[(t+1),,])
  test_that("Final pop size (N) = UK Total", expect_equal(N, total_population))
  
  #Return out array
  return(out)
}
