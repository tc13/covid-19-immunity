### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SEIIRRS discrete time model, age structed for UK population, gamma distributed waiting times, 
### Two infected (symtomatic & asymptomatic) and immunity classes, children are less susceptible to infection than adults (though equally infectious)
### Interventions: Model functions

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
                                     intervention = c(0.8,0.3,0.1,0.2)){#, # home, work, school, other
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

#check values are not less than 1 and integer
checkInteger <- function(x){
  if(x < 1){y <- 1}
  else if(x%%1!=0){y <- round(x, digits = 0)} 
  else y <- x
  return(y)
}

#Main model function
SEIIRRS_intervention <- function(R0=2.8, latent_mean=4.5, infectious_mean=3.07, immune_mean_1=90,
                                immune_mean_2=180, latent_shape=4, infectious_shape=2, immune_shape=2,
                                asymtomatic_relative_infectiousness=0.5, children_relative_susceptibility=0.4,
                                dt=1, days=400, BBC_contact_matrix, total_population=66435550, p_age, p_hospitalised, I_init,
                                Rt_partial=1.2, Rt_full=0.8, intervention_partial=c(1,0.8,0.85,0.75), 
                                intervention_full=c(0.8,0.3,0.1,0.2), threshold=10000, t_partial_intervention=14,
                                t_full_intervention=60, rho=c(rep(0.25, 4),rep(0.5,11))){
  
  #check shape parameters are integer values and >=1
  latent_shape = checkInteger(latent_shape)
  infectious_shape = checkInteger(infectious_shape)
  immune_shape = checkInteger(immune_shape)
  
  #susceptibility vector
  susceptibility <- c(rep(children_relative_susceptibility, 4), rep(1, 11))

  #time vector
  time <- seq(1, days, dt)
  
  #Transmission probability
  C <- make.intervention.matrix(BBC_contact_matrix, intervention=c(1,1,1,1), R0=R0)
  beta <- get_beta(C=C, p_age=p_age, gamma=(1/infectious_mean), R0=R0)
  
  #Transmission probability with intervention
  C_partial <- make.intervention.matrix(BBC_contact_matrix, intervention=intervention_partial, R0=Rt_partial)
  beta_partial <- get_beta(C=C_partial, p_age=p_age, gamma=(1/infectious_mean), R0=Rt_partial)
  C_full <- make.intervention.matrix(BBC_contact_matrix, intervention=intervention_full, R0=Rt_full)
  beta_full <- get_beta(C=C_full, p_age=p_age, gamma=(1/infectious_mean), R0=Rt_full)
  
  #Intervention time steps
  t_partial_intervention_step <- t_partial_intervention*(1/dt)
  t_full_intervention_step <- t_full_intervention*(1/dt)

  #Update parameter values (incoporating gamma shape and dt interval)
  sigma <- (1/latent_mean) * latent_shape * dt               #probability of becoming infectious
  gamma <- (1/infectious_mean) * infectious_shape * dt       #probability of recovery
  omega_1 <- get_omega(mean=immune_mean_1, shape=immune_shape, dt=dt)   #Duration of immunity in R class 1
  omega_2 <- get_omega(mean=immune_mean_2, shape=immune_shape, dt=dt)   #Duration of immunity in R class 2
  classes <- length(p_age)                #number of age classes
  N_age <- p_age * total_population       #Absolute numbers by age group
  
  #State array
  E.lab <- textNum("E", latent_shape)
  Ia.lab <- textNum("Ia", infectious_shape)
  Is.lab <- textNum("Is", infectious_shape)
  Rnh.lab <- textNum("Rnh", immune_shape)
  Rh.lab <- textNum("Rh", immune_shape)
  
  state_names <- c("S", E.lab, Ia.lab, Is.lab, Rnh.lab, Rh.lab)
  n_states <- length(state_names)
  state <- array(data=0, dim=c(length(time), classes, n_states), #rows are times, columns are age classes, matrixes are disease states
                 dimnames = list(as.character(time), colnames(C), state_names))
  
  #State variables at t=0
  state[1,,"S"] <- N_age #susceptible absolute numbers
  
  #Initialise infections by age group with I_init vector
  state[1,,"S"] <- state[1,,"S"]- I_init
  for(j in 1:infectious_shape){
    state[1,,Is.lab[j]] <- I_init/infectious_shape}
  
  #Check initial population size (N) = UK Total Pop
  N = sum(state[1,,])
  test_that("Initial pop size (N) = UK total pop", expect_equal(N, total_population))
  
  intervention_phase = rep(0, length(time))
  
  #For each timestep t
  for(t in 1:(length(time)-1)){
    
    #new symtomatic & asymtomatic infections per timestep
    if(infectious_shape>1){
    I_symtomatic = apply(state[t,,Is.lab], 1, sum) #sum Is sub-classes
    I_asymtomatic = apply(state[t,,Ia.lab], 1, sum) #sum Is sub-classes
    }else{
      I_symtomatic = state[t,,Is.lab]
      I_asymtomatic = state[t,,Ia.lab]
      }
    
    #Check if infection threshold has been reached, set intervention periods
    if(sum(I_symtomatic)>=threshold & intervention_phase[t]==0){
      #Partial Intervention period - triggered by threshold
      intervention_phase[(t+1):(t+t_partial_intervention_step)] <- 1
      #Full intervention period 
      intervention_phase[(t+t_partial_intervention_step+1):(t+t_partial_intervention_step+t_full_intervention_step)] <- 2
      #After full intervention; return to partial intervention Rt indefinately
      intervention_phase[(t+t_partial_intervention_step+t_full_intervention_step+1):length(time)] <- 1
    }
    
    #Transmission in absense of intervention
    if(intervention_phase[(t+1)]==0){
      lambda = beta * (as.matrix(C) %*% ((I_symtomatic/N_age) + ((I_asymtomatic*asymtomatic_relative_infectiousness)/N_age))) #force of infection (R0)
      new_infection = lambda * state[t,,"S"]* susceptibility * dt #new infections per timestep
    }else if(intervention_phase[(t+1)]==1){
      #Transmission during partial intervention period 
      lambda = beta_partial * (as.matrix(C_partial) %*% ((I_symtomatic/N_age) + ((I_asymtomatic*asymtomatic_relative_infectiousness)/N_age))) #force of infection (Rt)
      new_infection = lambda * state[t,,"S"]* susceptibility * dt #new infections per timestep
    }else{
      #Transmission during full intervention period 
      lambda = beta_full * (as.matrix(C_full) %*% ((I_symtomatic/N_age) + ((I_asymtomatic*asymtomatic_relative_infectiousness)/N_age))) #force of infection (Rt)
      new_infection = lambda * state[t,,"S"]* susceptibility * dt #new infections per timestep
    }
    
    #Epidemic transitions at time t+1
    #Susceptibles
    state[(t+1),,"S"] = state[t,,"S"] + omega_1*state[t,,Rh.lab[immune_shape]] + omega_2*state[t,,Rnh.lab[immune_shape]] - new_infection
    #First exposed class
    state[(t+1),,"E1"] = state[t,,"E1"] + new_infection - state[t,,"E1"]*sigma
    #Other exposed classes
    if(latent_shape>1){
      for(k in 2:latent_shape){
        state[(t+1),,E.lab[k]] = state[t,,E.lab[k]] + sigma*(state[t,,E.lab[(k-1)]]-state[t,,E.lab[k]])
      }
    }
    
    #first symtomatic infectious class
    state[(t+1),,"Is1"] = state[t,,"Is1"] + sigma*state[t,,E.lab[latent_shape]]*rho - state[t,,"Is1"]*gamma
    #Other symtomatic infectious classes
    if(infectious_shape>1){
      for(k in 2:infectious_shape){
        state[(t+1),,Is.lab[k]] = state[t,,Is.lab[k]] + gamma*(state[t,,Is.lab[(k-1)]]-state[t,,Is.lab[k]])
      }
    }
    
    #first asymtomatic infectious class
    state[(t+1),,"Ia1"] = state[t,,"Ia1"] + sigma*state[t,,E.lab[latent_shape]]*(1-rho) - state[t,,"Ia1"]*gamma
    #Other ainfectious classes
    if(infectious_shape>1){
      for(k in 2:infectious_shape){
        state[(t+1),,Ia.lab[k]] = state[t,,Ia.lab[k]] + gamma*(state[t,,Ia.lab[(k-1)]]-state[t,,Ia.lab[k]])
      }
    }
    
    #Duration of immunity, first classes
    state[(t+1),,"Rh1"] = state[t,,"Rh1"] + (p_hospitalised/rho)*gamma*state[t,,Is.lab[infectious_shape]] - state[t,,"Rh1"]*omega_1
    state[(t+1),,"Rnh1"] = state[t,,"Rnh1"] + gamma*state[t,,Ia.lab[infectious_shape]]+ (1-(p_hospitalised/rho))*gamma*state[t,,Is.lab[infectious_shape]] - state[t,,"Rnh1"]*omega_2
    if(immune_shape>1){
      for(k in 2:immune_shape){
        state[(t+1),,Rh.lab[k]] = state[t,,Rh.lab[k]] + omega_1*(state[t,,Rh.lab[(k-1)]] - state[t,,Rh.lab[k]])
        state[(t+1),,Rnh.lab[k]] = state[t,,Rnh.lab[k]] + omega_2*(state[t,,Rnh.lab[(k-1)]] - state[t,,Rnh.lab[k]])
      }
    }
    
    #Check population size (N) =  UK total pop
    N = sum(state[(t+1),,])
    test_that("Pop size (N) = Pop size", expect_equal(N, total_population))
    
    #Test that no state value is negative
    test_that("States non-negative", expect_true(all(state[(t+1),,]>=0)))
  }
  
  #Combine E, I, R substates into new array
  state_names_SEIR <- state_names <- c("S", "E", "Is", "Ia", "I", "Rh", "Rnh", "R")
  
  out <- array(data=0, dim = c(length(time), classes, length(state_names_SEIR)),
               dimnames = list(as.character(time), colnames(C), state_names_SEIR))
  out[,,"S"] <- state[,,"S"]
  for(i in 1:latent_shape){out[,,"E"] <- out[,,"E"] + state[,,E.lab[i]]}
  for(i in 1:infectious_shape){out[,,"Is"] <- out[,,"Is"] + state[,,Is.lab[i]]}  
  for(i in 1:infectious_shape){out[,,"Ia"] <- out[,,"Ia"] + state[,,Ia.lab[i]]}  
  for(i in 1:immune_shape){out[,,"Rh"] <- out[,,"Rh"] + state[,,Rh.lab[i]]}  
  for(i in 1:immune_shape){out[,,"Rnh"] <- out[,,"Rnh"] + state[,,Rnh.lab[i]]} 
  out[,,"I"] <- out[,,"Is"] + out[,,"Ia"]
  out[,,"R"] <- out[,,"Rh"] + out[,,"Rnh"]
  
  #Check population size (N) =  UK total pop
  N = sum(out[(t+1),,c("S", "E", "Is", "Ia", "Rh", "Rnh")])
  test_that("Final pop size (N) = UK Total", expect_equal(N, total_population))
  
  #Return out array
  return(out)
}
