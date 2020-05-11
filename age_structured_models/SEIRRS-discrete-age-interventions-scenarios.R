### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SEIRRS discrete time model, age structed for UK population, gamma (Erlang) distributed waiting times, two immunity classes
### Model Scenarios with interventions

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
source("age_structured_models/SEIRRS-discrete-age-interventions-functions.R")

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
Rt <- 0.8 #Effective reproduction number during intervention
intervention <- c(1,0.3,0.1,0.2)

#Obtain overall contact matrix before and after interventions
C <- make.intervention.matrix(BBC_contact_matrix, intervention = c(1,1,1,1))

#Check Rt & beta after interventions
C_i <- make.intervention.matrix(BBC_contact_matrix, intervention = c(1,0.3,0.1,0.2), R0=Rt)
beta_Rt <- get_beta(C=C_i, p_age=p_age, gamma=1/gamma_recip, R0=Rt)
eigen_int <- eigen(C_i)
eig_int <- max(eigen_int$values)
beta_Rt*eig_int/(1/gamma_recip)

#Get total UK population
total_pop <- sum(uk.pop.2018.count$total)
p_age <- uk.pop.2018.count$total/total_pop

#initial infected vector; 50 individuals per working age group (20-59)
I_init <- c(0,0,0,0,50,50,50,50,50,50,50,50,0,0,0)

#Proportion of each age group in immunity class 1
p_immunity <- c(1, 1, 1, 1, 0.95, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.70, 0.65, 0.60, 0.55)

#SEIRRS age structure model
out <- SEIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, 
                            immune_mean_1=0, immune_mean_2=0, p_immunity=p_immunity,
                            latent_shape=4, infectious_shape=n, immune_shape=o, dt=dt, days=days, C=C, 
                            total_population=total_pop, p_age=p_age, I_init=I_init, Rt = Rt, 
                            intervention = intervention, threshold = 500000, t_intervention = 60 )

#plot total S, I, R over time (sum over age classes)
S <- apply(out[,,"S"],1,sum)
E <- apply(out[,,"E"],1,sum)
I <- apply(out[,,"I"],1,sum)
R1 <- apply(out[,,"R1"],1,sum)
R2 <- apply(out[,,"R2"],1,sum)
R <- R1 + R2

plot(S~time, type="l", ylim=c(0,total_pop))
lines(I~time)
lines(R~time)

#plot I for seperate age classes
plot(out[,1,"I"]~time, type="l", ylim=c(0,700000))
for(k in 2:15){lines(out[,k,"I"]~time)}

