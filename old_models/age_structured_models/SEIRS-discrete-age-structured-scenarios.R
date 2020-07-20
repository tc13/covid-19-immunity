### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SEIRS discrete time model, age structed for UK population, gamma (Erlang) distributed waiting times
### Model Scenarios

#Clear R environment
remove(list = ls())

#packages
require(rstudioapi)

#set path
current_path <- getActiveDocumentContext()$path 
head_directory <- gsub(x=current_path, pattern="covid-19-immunity.*", replacement = "covid-19-immunity/")
setwd(head_directory)

#Load UK contact matrix and population data
load("data/BBC_contact_matrices_population_vector.RData")
load("data/uk_population_size_by_age_2018.RData")

#Source model functions
source("age_structured_models/SEIRS-discrete-age-structured-functions.R")

#parameters
R0 <- 2.8                   #Basic reproduction number (in the absense of interventions) [From Petra]
sigma_recip <- 4.5          #Average duration of latent period (days)     [From He, X., et al. Temporal dynamics in viral shedding and transmissibility of COVID-19. Nat Med (2020). https://doi.org/10.1038/s41591-020-0869-5]
gamma_recip <- 3.07         #Average duration of infectiousness (days)    [From He, X., et al.]
omega_recip <- 90           #Average duration of immunity (days) [user specified]
m <- 4                      #Shape paremeter, latent period
n <- 2                      #Shape parameter, infectious period
o <- 2                      #Shape parameter, immune period
dt <- 1                     #time period (days)
days <- 400                 #number of days 
time <- seq(1, days, dt)    #time vector

#Obtain overall contact matrix before and after interventions
C <- make.intervention.matrix(BBC_contact_matrix, intervention = c(1,1,1,1))

#Get total UK population
total_pop <- sum(uk.pop.2018.count$total)
p_age <- uk.pop.2018.count$total/total_pop

#initial infected vector; 50 individuals per working age group (20-59)
I_init <- c(0,0,0,0,50,50,50,50,50,50,50,50,0,0,0)

#SEIRS age structure model
out <- SEIRS_age_structure(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, immune_mean=90,
                    latent_shape=m, infectious_shape=n, immune_shape=o, dt=dt, days=days, C=C, 
                    total_population=total_pop, p_age=p_age, I_init=I_init )

#plot total S, I, R over time (sum over age classes)
S <- apply(out[,,"S"],1,sum)
E <- apply(out[,,"E"],1,sum)
I <- apply(out[,,"I"],1,sum)
R <- apply(out[,,"R"],1,sum)

plot(S~time, type="l", ylim=c(0,total_pop))
lines(I~time)
lines(R~time)

#plot I for seperate age classes
plot(out[,1,"I"]~time, type="l", ylim=c(0, 1000000))
for(k in 2:15){lines(out[,k,"I"]~time)}
