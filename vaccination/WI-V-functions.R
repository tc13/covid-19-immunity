### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen, University of Oxford, thomas.crellen@bdi.ox.ac.uk, April 2020
### Adapted from code by Drs Kiesha Prem and Petra Klepac (LSHTM)
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
  #overall_mat <-  sum.mat(home_mat, work_mat, school_mat, other_mat, c = c(1, 1, 1, 1))
  #rownames(overall_mat) <- colnames(home_mat)
  int_mat <- sum.mat(home_mat, work_mat, school_mat, other_mat, c = intervention)
  
  # the matrix has participant ages on x-axis (columns), and contact ages on y-axis (rows)
  rownames(int_mat) <-colnames(home_mat)
  return(int_mat)
}

#Get dominant eigenvalue from matrix
get_eigen <- function(m){
  eigs <- eigen(m)
  dominant <- max(abs(Re(eigs$values)))
  return(dominant)
}

#Calculate K matrix (given beta)
get_K <- function(C, R0, beta=1,mean_infectious, 
                  susceptibility, prop_asymtomatic, 
                  asymtomatic_infectiousness){
  
  #Calculate next generation matrix K
  K = matrix(data=0, nrow=nrow(C), ncol=ncol(C))
  #Kij = infections in group i produced by individuals in group j
  for(i in 1:nrow(K)){  
    for(j in 1:nrow(K)){
      K[i,j] = mean_infectious*beta*(prop_asymtomatic[j]*asymtomatic_infectiousness*C[i,j]+C[i,j]*(1-prop_asymtomatic[j]))*susceptibility[i]
    }
  }
  return(K)
}

#Calculate transmission parameter
get_beta <- function(C, R0, mean_infectious, 
                     susceptibility, prop_asymtomatic, 
                     asymtomatic_infectiousness){ 
  
  #Calculate next generation matrix K
  K = get_K(C=C, R0=R0, mean_infectious = mean_infectious,
            susceptibility = susceptibility, 
            prop_asymtomatic=prop_asymtomatic,
            asymtomatic_infectiousness = asymtomatic_infectiousness)
  #Dominant eigenvalue of K
  dominant = get_eigen(K)
  
  #Calculate beta from R0 and K
  beta = R0/dominant
  
  #check new eigenvalue = R0
  NGM <- K*beta
  new_eigen <- get_eigen(NGM)
  test_that("NGM=R0", expect_equal(new_eigen, R0))
  return(beta)
}

#Get dominant eigenvector from NGM
get_eigenvector <- function(C, R0, mean_infectious, 
                     susceptibility, prop_asymtomatic, 
                     asymtomatic_infectiousness){ 
  
  #Calculate next generation matrix K
  K = get_K(C=C, R0=R0, mean_infectious = mean_infectious,
            susceptibility = susceptibility, 
            prop_asymtomatic=prop_asymtomatic,
            asymtomatic_infectiousness = asymtomatic_infectiousness)
  
  #Dominant eigenvalue of K
  dominant = get_eigen(K)
  
  #Calculate beta from R0 and K
  beta = R0/dominant
  
  #check new eigenvalue = R0
  NGM <- K*beta
  new_eigen <- get_eigen(NGM)
  
  #Get dominant eigenvector
  vectors <- abs(Re(eigen(NGM)$vectors))
  dom_eigenvec <- vectors[,1] #first column
  test_that("NGM=R0", expect_equal(new_eigen, R0))
  return(dom_eigenvec)
}

#Get omega function (rate of immunity decay)
get_omega <- function(mean, shape, dt, name){
  if(mean<=0){ 
    print(paste0("Immunity does not wane in class ", name))
    omega <- 0}
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
  else if(x%%1!=0){
    print("Shape parameter rounded to integer")
    y <- round(x, digits = 0)} 
  else y <- x
  return(y)
}

#Reporting function - compare model output against data
report <- function(mod, dates){
  R = apply(mod[,,"R"],1,sum)
  new = apply(mod[,,"new_infections"],1,sum)
  #New cases on 23rd March between 54-155k
  indx1 = which(dates=="2020-03-23")
  test1 = round(new[indx1])
  print(paste("New cases on 23rd March 2020= ", test1))
  #Proportion immune 13th July
  indx2 = which(dates=="2020-07-13")
  test2.abs = round(sum(mod[indx2,4:15,"R"]))
  adult_pop = sum(england_pop[14:15])
  test2.prop = signif(R[indx2]/total_pop)
  #print(paste("Adult pop with antibodies on 13th July 2020 = ", test2.abs))
  print(paste("Prop with antibodies on 13th July 2020 = ", test2.prop))
  #Number of cases at 2nd peak
  indx3 = which(dates=="2020-12-30")
  test3 = round(new[indx3])
  print(paste("New cases on 30th Dec 2020 = ", test3))
  #Number of cases April 4th 2021
  indx4 = which(dates=="2021-04-04")
  test4 = round(new[indx4])
  print(paste("New cases on 4th April 2021 = ", test4))
}

#Reporting function - compare model output against data
future <- function(mod, dates, current_date=as.Date("2021-05-01"), R, case_threshold, total_pop=56286961){
  #Get index of current date
  indx1 = which(dates==current_date)
  #New cases after current date
  new = apply(mod[,,"new_infections"],1,sum)
  future_infections = new[indx1:length(new)]
  #Check if threshold of new cases reached
  resurge = which(future_infections > case_threshold)
  if(length(resurge)>0){
    date_resurge = dates[min(resurge)+(indx1-1)]
    print(paste("New cases exceed", case_threshold, "on date", date_resurge))
  }else{
    print(paste("New cases do not exceed", case_threshold, "before end of simulation", max(dates)))}
  #Calculate Rt
  #R0 in presence of interventions after current date
  S = apply(mod[,,"S"],1,sum)
  Rt = R*(S[indx1:length(S)]/total_pop)
  Rt_threshold = which(Rt>1.00)
  if(length(Rt_threshold>0)){
    date_threshold = dates[min(Rt_threshold)+(indx1-1)]
    #print(paste("Rt exceeds 1 on date", date_threshold))
    #print(paste("Proportion susceptible at this time =", signif(S[min(Rt_threshold)+(indx1-1)]/total_pop, digits = 3)))
  }else{
   # print(paste("Rt never exceeds 1 before end of simulation", max(dates)))
  }
}

#Catch erros
catch <- function(x) if(any(!grepl('passed', x))) cat(x, sep='\n')

#Main model function
SEIIRRS_intervention <- function(
        #Mean durations
        latent_mean=4.5, infectious_mean=3.1, 
        immune_mean_1=90, immune_mean_2=180,
        #Shape parameters
        latent_shape=4, infectious_shape=2, immune_shape=2,
        #Vaccination
        vaccine_immune=365, vaccine_shape=2, 
        vaccine_start_days,  #start days by age group
        vaccine_rates_daily, #daily rates per age group
        vaccination_length,  #duration of vaccination period in days
        #Time steps
        dt=0.01, days=400, 
        #Population
        p_age, total_population=66435550, BBC_contact_matrix, 
        #Probabilities
        children_relative_susceptibility=0.4, 
        asymtomatic_relative_infectiousness=0.5,
        p_hospitalised, phi=c(rep(0.75, 4),rep(0.5,11)),
        #Intervention params
        R0=4.0, start_days=NULL, durations=NULL, R_values=NULL, contact_interventions=NULL,
        #Initial values
        I_init, E_init=rep(0,15)){
  
  #check shape parameters are integer values and >=1
  latent_shape = checkInteger(latent_shape)
  infectious_shape = checkInteger(infectious_shape)
  immune_shape = checkInteger(immune_shape)
  
  #Fix time step during development
  if(dt > 1){
    print("Time step set to 1 day")
    dt <- 1} #dt cannot be greater than 1
  
  #time vector
  time <- seq(0, days, dt)
  
  #Vaccine rates
  ### SET START TIMES ###
  vaccination_schedule <- list()
  for(i in 1:length(p_age)){
    vaccination_days <- numeric(length(time)) 
    onset = match(vaccine_start_days[i],time)
    term = match(vaccine_start_days[i]+vaccination_length, time)
    vaccination_days[onset:term] <-  vaccine_rates_daily[i]*dt
    vaccination_schedule[[i]] = vaccination_days
  }
  theta <- matrix(unlist(vaccination_schedule), byrow=F,
                  ncol=length(p_age), nrow=length(time))
  
  #Update gamma (Erlang) distribution parameter values (incoporating gamma shape and dt interval) - see Krylova & Earn 2013 - https://royalsocietypublishing.org/doi/pdf/10.1098/rsif.2013.0098
  sigma <- (1/latent_mean) * latent_shape * dt               #probability of becoming infectious
  gamma <- (1/infectious_mean) * infectious_shape * dt        #probability of recovery
  omega_1 <- get_omega(mean=immune_mean_1, shape=immune_shape, dt=dt, name="R1")   #Duration of immunity in R class 1
  omega_2 <- get_omega(mean=immune_mean_2, shape=immune_shape, dt=dt, name="R2")   #Duration of immunity in R class 2
  omega_v <- get_omega(mean=vaccine_immune, shape=vaccine_shape, dt=dt, name="vaccinated") #Duration of immunity from vaccination
  classes <- length(p_age)                #number of age classes
  N_age <- p_age * total_population       #Absolute numbers by age group
  
  #susceptibility vector
  susceptibility <- c(rep(children_relative_susceptibility, 4), rep(1, 11))
  
  #Contact matrix & transmission parameter (beta = probability of infection, given contact)
  C = make.intervention.matrix(BBC_contact_matrix, intervention=c(1,1,1,1))
  beta = get_beta(C=C, R0=R0, mean_infectious=infectious_mean,
                  susceptibility=susceptibility, prop_asymtomatic=phi,
                  asymtomatic_infectiousness=asymtomatic_relative_infectiousness)
  
  intervention_phase <- numeric(length(time))
    if(!is.null(start_days)){
      #Better intervention implementation
      contact_matrxs <- list()
      betas <- c()
  
      intervention_start <- match(start_days, time)
      intervention_end <- match(start_days+durations, time)
  
      for(i in 1:length(start_days)){
        #Update intervention vector
        intervention_phase[intervention_start[i]:intervention_end[i]] <- i
        #Generate contact matrix
        contact_matrxs[[i]] <- make.intervention.matrix(BBC_contact_matrix, intervention=contact_interventions[[i]])
        #Calculate transmission parameter
        betas[i] <- get_beta(C=contact_matrxs[[i]], R0=R_values[i], mean_infectious=infectious_mean,
                         susceptibility=susceptibility, prop_asymtomatic=phi,
                         asymtomatic_infectiousness=asymtomatic_relative_infectiousness)
        }
    }
  
  #State array
  E.lab <- textNum("E", latent_shape)
  Ia.lab <- textNum("Ia", infectious_shape)
  Is.lab <- textNum("Is", infectious_shape)
  Rnh.lab <- textNum("Rnh", immune_shape)
  Rh.lab <- textNum("Rh", immune_shape)
  V.lab <- textNum("V", vaccine_shape)
  
  state_names <- c("S", E.lab, Ia.lab, Is.lab, Rnh.lab, Rh.lab, V.lab)
  n_states <- length(state_names)
  state <- array(data=0, dim=c(length(time), classes, n_states), #rows are times, columns are age classes, matrices are disease states
                 dimnames = list(as.character(time), colnames(C), state_names))
  
  #Array for daily new infections and new hospitalisations
  daily_counts <- array(data=0, dim=c(length(time), classes, 2),
                dimnames = list(as.character(time), colnames(C), 
                            c("new_infections", "new_hospitalisations")))
  
  #State variables at t=0
  state[1,,"S"] <- N_age #susceptible absolute numbers
  
  #Initialise infections by age group with E_init and I_init (initial number of infections and exposed) vector
  state[1,,"S"] <- state[1,,"S"]- (E_init + I_init)
  for(j in 1:infectious_shape){
    state[1,,Is.lab[j]] <- I_init*(1-phi)/infectious_shape
    state[1,,Ia.lab[j]] <- I_init*phi/infectious_shape
    }
  
  for(j in 1:latent_shape){
    state[1,,E.lab[j]] <- E_init/latent_shape
  }
  
  #Check initial population size (N) = UK Total Pop
  N = sum(state[1,,])
  test_that("Initial pop size (N) = UK total pop", expect_equal(N, total_population))
  
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
         
    #Check for intervention period
    if(intervention_phase[(t+1)]==0){
      #Transmission at onset of epidemic
      lambda = beta * (C %*% (I_symtomatic + I_asymtomatic*asymtomatic_relative_infectiousness) * (1/N_age)) #force of infection (R0)
      new_infection = lambda * state[t,,"S"]* susceptibility * dt #new infections per timestep
    }else{
      indx = intervention_phase[(t+1)]
      #Transmission during intervention period i 
      lambda = betas[indx] * (contact_matrxs[[indx]] %*% (I_symtomatic + I_asymtomatic*asymtomatic_relative_infectiousness) * (1/N_age)) #force of infection (Rt)
      new_infection = lambda * state[t,,"S"] * susceptibility * dt #new infections per timestep
      }
    #update daily counts of new infections
    daily_counts[(t+1), ,"new_infections"] <- new_infection
    
    #Epidemic transitions at time t+1
    #Susceptibles
    state[(t+1),,"S"] = state[t,,"S"] + omega_1*state[t,,Rh.lab[immune_shape]] + omega_2*state[t,,Rnh.lab[immune_shape]] - new_infection - theta[t,]*(state[t,,"S"]-new_infection) + omega_v*state[t,,V.lab[vaccine_shape]]
    #Vaccinated class (first)
    state[(t+1),,"V1"] = state[t,,"V1"] + theta[t,]*(state[t,,"S"]-new_infection) - omega_v*state[t,,"V1"]
    #Other vaccinated classes
    if(vaccine_shape>1){
      for(k in 2:vaccine_shape){
        state[(t+1),,V.lab[k]] = state[t,,V.lab[k]] + omega_v*(state[t,,V.lab[(k-1)]]-state[t,,V.lab[k]])
      }
    }
    #First exposed class
    state[(t+1),,"E1"] = state[t,,"E1"] + new_infection - state[t,,"E1"]*sigma
    #Other exposed classes
    if(latent_shape>1){
      for(k in 2:latent_shape){
        state[(t+1),,E.lab[k]] = state[t,,E.lab[k]] + sigma*(state[t,,E.lab[(k-1)]]-state[t,,E.lab[k]])
      }
    }
    
    #first symtomatic infectious class
    state[(t+1),,"Is1"] = state[t,,"Is1"] + sigma*state[t,,E.lab[latent_shape]]*(1-phi) - state[t,,"Is1"]*gamma
    #Other symtomatic infectious classes
    if(infectious_shape>1){
      for(k in 2:infectious_shape){
        state[(t+1),,Is.lab[k]] = state[t,,Is.lab[k]] + gamma*(state[t,,Is.lab[(k-1)]]-state[t,,Is.lab[k]])
      }
    }
    
    #first asymtomatic infectious class
    state[(t+1),,"Ia1"] = state[t,,"Ia1"] + sigma*state[t,,E.lab[latent_shape]]*phi - state[t,,"Ia1"]*gamma
    #Other asymtomatic infectious classes
    if(infectious_shape>1){
      for(k in 2:infectious_shape){
        state[(t+1),,Ia.lab[k]] = state[t,,Ia.lab[k]] + gamma*(state[t,,Ia.lab[(k-1)]]-state[t,,Ia.lab[k]])
      }
    }
    
    #Duration of immunity, first classes
    state[(t+1),,"Rh1"] = state[t,,"Rh1"] + (p_hospitalised/(1-phi))*gamma*state[t,,Is.lab[infectious_shape]] - state[t,,"Rh1"]*omega_1
    state[(t+1),,"Rnh1"] = state[t,,"Rnh1"] + gamma*state[t,,Ia.lab[infectious_shape]]+ (1-(p_hospitalised/(1-phi)))*gamma*state[t,,Is.lab[infectious_shape]] - state[t,,"Rnh1"]*omega_2
    if(immune_shape>1){
      for(k in 2:immune_shape){
        state[(t+1),,Rh.lab[k]] = state[t,,Rh.lab[k]] + omega_1*(state[t,,Rh.lab[(k-1)]] - state[t,,Rh.lab[k]])
        state[(t+1),,Rnh.lab[k]] = state[t,,Rnh.lab[k]] + omega_2*(state[t,,Rnh.lab[(k-1)]] - state[t,,Rnh.lab[k]])
      }
    }
    
    #update daily counts of hospitalisation
    daily_counts[(t+1),,"new_hospitalisations"] <- (p_hospitalised/(1-phi))*gamma*state[t,,Is.lab[infectious_shape]]
    
    #Check population size (N) =  UK total pop
    N = sum(state[(t+1),,])
    pop_size_test = capture.output(test_that("Pop size (N) = Pop size", expect_equal(N, total_population)))
    catch(pop_size_test)
    
    #Test that state values are positive
    non_negative = capture.output(test_that("States non-negative", expect_true(all(state[(t+1),,]>=0))))
    catch(non_negative)
  }
  
  #Combine V, E, I, R substates into new array
  state_names_SEIR <- state_names <- c("S", "V","E", "Is", "Ia", "I", "Rh", "Rnh", "R", "new_infections", "new_hospitalisations")
  
  #Daily time points
  daily_indx <- time%%1==0
  
  out <- array(data=0, dim = c(length(time[daily_indx]), classes, length(state_names_SEIR)),
               dimnames = list(as.character(time[daily_indx]), colnames(C), state_names_SEIR))
  out[,,"S"] <- state[daily_indx,,"S"]
  for(i in 1:vaccine_shape){out[,,"V"] <- out[,,"V"] + state[daily_indx,,V.lab[i]]}
  for(i in 1:latent_shape){out[,,"E"] <- out[,,"E"] + state[daily_indx,,E.lab[i]]}
  for(i in 1:infectious_shape){out[,,"Is"] <- out[,,"Is"] + state[daily_indx,,Is.lab[i]]}  
  for(i in 1:infectious_shape){out[,,"Ia"] <- out[,,"Ia"] + state[daily_indx,,Ia.lab[i]]}  
  for(i in 1:immune_shape){out[,,"Rh"] <- out[,,"Rh"] + state[daily_indx,,Rh.lab[i]]}  
  for(i in 1:immune_shape){out[,,"Rnh"] <- out[,,"Rnh"] + state[daily_indx,,Rnh.lab[i]]} 
  out[,,"I"] <- out[,,"Is"] + out[,,"Ia"]
  out[,,"R"] <- out[,,"Rh"] + out[,,"Rnh"]
  
  out[,,"new_infections"] <- apply(daily_counts[,,"new_infections"],2,FUN=function(x){unname(tapply(x, (seq_along(x)-1) %/% (1/dt), sum))})
  out[,,"new_hospitalisations"] <- apply(daily_counts[,,"new_hospitalisations"],2,FUN=function(x){unname(tapply(x, (seq_along(x)-1) %/% (1/dt), sum))})
  
  #Check population size (N) =  UK total pop
  N = sum(out[(max(time)),,c("S", "V", "E", "Is", "Ia", "Rh", "Rnh")])
  test_that("Final pop size (N) = UK Total", expect_equal(N, total_population))
  
  #Return out array
  return(out)
}
