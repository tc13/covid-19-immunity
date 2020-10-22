### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, June 2020
### SEIIRRS discrete time model, age structed for UK population, gamma (Erlang) distributed waiting times, two immunity classes
### Model Scenarios with interventions - Rt following lockdown 1.3-1.5

#Clear R environment
#remove(list = ls())

#packages
require(lubridate)

#set path to current folder
dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir)

#Load UK contact matrix and population data
load("data/BBC_contact_matrices_population_vector.RData")
load("data/uk_population_size_by_age_2018.RData")

#Source model functions
source("SEIIRRS-discrete-age-interventions-functions.R")

###########################################
## SARS-CoV-2 Natural history parameters ##
###########################################
R0 <- 2.8                 #Basic reproduction number (in the absense of interventions) [From Petra]
sigma_recip <- 4.5        #Average duration of latent period (days)     [From He, X., et al. Temporal dynamics in viral shedding and transmissibility of COVID-19. Nat Med (2020). https://doi.org/10.1038/s41591-020-0869-5]
gamma_recip <- 3.07       #Average duration of infectiousness (days)    [From He, X., et al.]
m <- 4                    #Shape paremeter, latent period
n <- 2                    #Shape parameter, infectious period
#p_hosp <- c(0, 0, 0.000408, 0.000408, 0.0104, 0.0104, 0.0343, 0.0343, 0.0425, 0.0425, 0.0816, 0.0816, 0.118, 0.118, 0.175) #Proportion of each age group that becomes hospitalised (Rh) [Verity et al. 2020 Lancet ID]
p_hosp <- c(0.001, 0.001, 0.003, 0.003, 0.012, 0.012, 0.032, 0.032, 0.049, 0.049, 0.102, 0.102, 0.166, 0.166, 0.26) #From Ferguson et al. Report 9
phi <- c(rep(0.75, 4),rep(0.5,11))  #proportion asymtomatic | infection by age group
crs <- 0.4                          #children relative susceptibility
ari <- 0.5                          #asymtomatic relative infectiousness

############################
## Intervention variables ##
############################
Rt_lockdown <- 0.8              #Effective reproduction number during lockdown (full intervention)
#Rt_half <- (R0+Rt_lockdown)/2   #Effective reproduction number at start of lockdown in tpi period (mid-point between R0 and Rt lockdown)
#Rt_post_lockdown <- 1.1         #Effective reproduction number after lockdown ends
Rt_partial <- 1.1
intervention_post_lockdown <- c(home=1,work=0.8,school=0.85, other=0.75)    #Scaling of contact matrix during partial intervention period
intervention_full <- c(home=0.8,work=0.3,school=0.1,other=0.2)        #Scaling of contact matrix during full intervention period
threshold <- 85000      #Number of symtomatic cases required to trigger partial intervention
tpi <- 21               #Length of partial intevention in days
tfi <- 60               #Length of full intervention in days
lockdown_day <- 49     #If triggering intervention by number of days, days after simulation starts

#####################
## Time  variables ##
#####################
dt <- 1                   #time period (days)
days <- 365*2             #number of days 
time <- seq(1, days, dt)    #time vector

######################################
## UK Population & Contact Matrixes ##
######################################
total_pop <- sum(uk.pop.2018.count$total)   #Get total UK population
p_age <- uk.pop.2018.count$total/total_pop  #Proportion of population in each age group
I_init <- c(0,0,0,0,25,25,25,25,25,25,25,25,0,0,0)  #Initial numbers infected: 25 individuals per working age group (20-59); total = 200

###################
## Run Scenarios ##
###################

#Rt post lockdown = 1.3
#Scenario 1: hospitalised & non-hospitalised patients have permanent immunity
S1.R13 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=0, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.3, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 2: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 365 days
S2.R13 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=365, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.3, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 3: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 180 days
S3.R13 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=180, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.3, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 4: hospitalised patients have immunity for mean 365 days, non-hospitalised patients have immunity for mean 90 days
S4.R13 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=365, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.3, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

S1.I.R13 <- apply(S1.R13[,,"I"], 1, sum)
S2.I.R13 <- apply(S2.R13[,,"I"], 1, sum)
S3.I.R13 <- apply(S3.R13[,,"I"], 1, sum)
S4.I.R13 <- apply(S4.R13[,,"I"], 1, sum)

plot(S1.I.R13~time, type="l", ylim=c(0,7.5E5))
lines(S2.I.R13~time)
lines(S3.I.R13~time)
lines(S4.I.R13~time)

#Rt post lockdown = 1.4
#Scenario 1: hospitalised & non-hospitalised patients have permanent immunity
S1.R14 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=0, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.4, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 2: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 365 days
S2.R14 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=365, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.4, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 3: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 180 days
S3.R14 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=180, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.4, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 4: hospitalised patients have immunity for mean 365 days, non-hospitalised patients have immunity for mean 90 days
S4.R14 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=365, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.4, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

S1.I.R14 <- apply(S1.R14[,,"I"], 1, sum)
S2.I.R14 <- apply(S2.R14[,,"I"], 1, sum)
S3.I.R14 <- apply(S3.R14[,,"I"], 1, sum)
S4.I.R14 <- apply(S4.R14[,,"I"], 1, sum)

plot(S1.I.R14~time, type="l", ylim=c(0,12E5))
lines(S2.I.R14~time)
lines(S3.I.R14~time)
lines(S4.I.R14~time)


#Rt post lockdown = 2
#Scenario 1: hospitalised & non-hospitalised patients have permanent immunity
S1.R20 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=0, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=2.0, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 2: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 365 days
S2.R20 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=365, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=2.0, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 3: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 180 days
S3.R20 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=180, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=2.0, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 4: hospitalised patients have immunity for mean 365 days, non-hospitalised patients have immunity for mean 90 days
S4.R20 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=365, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=2.0, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

S1.I.R20 <- apply(S1.R14[,,"I"], 1, sum)
S2.I.R20 <- apply(S2.R14[,,"I"], 1, sum)
S3.I.R20 <- apply(S3.R14[,,"I"], 1, sum)
S4.I.R20 <- apply(S4.R14[,,"I"], 1, sum)

plot(S1.I.R20~time, type="l", ylim=c(0,11E5))
lines(S2.I.R20~time)
lines(S3.I.R20~time)
lines(S4.I.R20~time)

