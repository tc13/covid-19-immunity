#Equilibirum analysis SARS-CoV-2

### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SEIIRRS discrete time model, age structed for UK population, gamma (Erlang) distributed waiting times, two immunity classes
### Model Scenarios with interventions

#Clear R environment
remove(list = ls())

#set path to parent folder
dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dirname(dir))

####################################################
## Source Model Functions, Parameters & Scenarios ##
####################################################
source("SEIIRRS-functions.R")
source("SEIIRRS-parameters.R")

######################################
## UK Population & Contact Matrixes ##
######################################
load("data/BBC_contact_matrices_population_vector.RData")
load("data/uk_population_size_by_age_2018.RData")
total_pop <- sum(uk.pop.2018.count$total)   #Get total UK population
p_age <- uk.pop.2018.count$total/total_pop  #Proportion of population in each age group

#####################
## Time  variables ##
#####################
dt <- 0.01                   #time period (days)
days <- 365*5             #number of days 
time <- seq(1, days, dt)    #time vector

##########################
## Run 5 Year Scenarios ##
##########################

#Rt post lockdown = 0.9
#Scenario 1: hospitalised & non-hospitalised patients have permanent immunity
S1.R09.eq <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=0, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=0.9, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 2: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 365 days
S2.R09.eq <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=365, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=0.9, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 3: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 180 days
S3.R09.eq <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=180, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=0.9, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 4: hospitalised patients have immunity for mean 365 days, non-hospitalised patients have immunity for mean 90 days
S4.R09.eq <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=365, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=0.9, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Rt post lockdown = 1.0
#Scenario 1: hospitalised & non-hospitalised patients have permanent immunity
S1.R10.eq <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=0, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.0, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 2: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 365 days
S2.R10.eq <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=365, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.0, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 3: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 180 days
S3.R10.eq <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=180, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.0, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 4: hospitalised patients have immunity for mean 365 days, non-hospitalised patients have immunity for mean 90 days
S4.R10.eq <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=365, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.0, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Rt post lockdown = 1.1
#Scenario 1: hospitalised & non-hospitalised patients have permanent immunity
S1.R11.eq <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=0, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.1, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 2: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 365 days
S2.R11.eq <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=365, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.1, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 3: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 180 days
S3.R11.eq <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=180, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.1, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 4: hospitalised patients have immunity for mean 365 days, non-hospitalised patients have immunity for mean 90 days
S4.R11.eq <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=365, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.1, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Rt post lockdown = 1.2
#Scenario 1: hospitalised & non-hospitalised patients have permanent immunity
S1.R12.eq <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=0, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.2, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 2: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 365 days
S2.R12.eq <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=365, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.2, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 3: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 180 days
S3.R12.eq <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=180, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.2, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 4: hospitalised patients have immunity for mean 365 days, non-hospitalised patients have immunity for mean 90 days
S4.R12.eq <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=365, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.2, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Date variables
date_UK_lockdown <- "2020-03-23"
start_date <- as.Date(date_UK_lockdown)-lockdown_day
hour_period <- 24*dt
dates <- seq(ymd(start_date),ymd(start_date+(days-1)), by = as.difftime(hours(hour_period)))

#Figures for Table 2 -function to extract relevant values
equilibrium <- function(model, dates, end_intervention=tpi+tfi+lockdown_day, population=total_pop){
  new_infections = apply(model[,,"new_infections"], 1, sum)
  new_infections = new_infections[-length(new_infections)]
  #dates= dates[-length(dates)]
  plot(new_infections~dates)
  max_days = length(new_infections)
  immune = apply(model[,,"R"], 1, sum)
  immune = immune[-length(immune)]
  #Check if new cases fall below zero
  if(!all(new_infections[end_intervention:max_days]>1)){
    day_zero_indx = which(new_infections<1)
    day_zero = min(day_zero_indx[day_zero_indx>end_intervention])
    date_zero = dates[day_zero]
    print(paste(as.character(date_zero), " is when new infections drop below zero"))
    immune_prop = round(immune[day_zero]/population, digits = 3)
    print(paste(as.character(immune_prop*100), "% is the percent immune", sep = ""))
  } else{
    #Else check if steady state reached
    delta = abs(diff(new_infections))
    delta2 = abs(diff(delta))
    eq1_indx = which(delta<0.1)
    eq2_indx = which(delta2<0.004)
    eq_both_indx = intersect(eq1_indx, eq2_indx)
    if(length(eq_both_indx)==0){day_eq=max_days
    
    }else{
    day_eq =  min(eq_both_indx[eq_both_indx>end_intervention])}
    date_eq = dates[day_eq]
    print(paste(as.character(date_eq), " is when endemic equilibirum reached"))
    new_inf_eq = round(new_infections[day_eq],digits = 0)
    print(paste(as.character(new_inf_eq), " daily new infections at equilibirum"))
    new_hospitalisations = apply(model[,,"new_hospitalisations"], 1, sum)
    new_hosp_eq = round(new_hospitalisations[day_eq],digits = 0)
    print(paste(as.character(new_hosp_eq), " daily new hospitalisations at equilibirum"))
    new_ICU_eq = round(new_hosp_eq*0.17)
    print(paste(as.character(new_ICU_eq), " daily new ICU admissions at equilibirum"))
    immune_prop = round(immune[day_eq]/population, digits = 3)
    print(paste(as.character(immune_prop*100), "% is the percent immune", sep = ""))
    }
}
##########################################
## Apply equilibirum function to models ##
##########################################

#Rt=0.9
equilibrium(model=S1.R09.eq, dates=dates)
equilibrium(model=S2.R09.eq, dates=dates)
equilibrium(model=S3.R09.eq, dates=dates)
equilibrium(model=S4.R09.eq, dates=dates)

#Rt=1.0
equilibrium(model=S1.R10.eq, dates=dates)
equilibrium(model=S2.R10.eq, dates=dates)
equilibrium(model=S3.R10.eq, dates=dates)
equilibrium(model=S4.R10.eq, dates=dates)

#Rt=1.1
equilibrium(model=S1.R11.eq, dates=dates)
equilibrium(model=S2.R11.eq, dates=dates)
equilibrium(model=S3.R11.eq, dates=dates)
equilibrium(model=S4.R11.eq, dates=dates)

#Rt=1.2
equilibrium(model=S1.R12.eq, dates=dates)
equilibrium(model=S2.R12.eq, dates=dates)
equilibrium(model=S3.R12.eq, dates=dates)
equilibrium(model=S4.R12.eq, dates=dates)
