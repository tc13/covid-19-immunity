### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SEIRRS discrete time model, age structed for UK population, gamma (Erlang) distributed waiting times, two immunity classes
### Model Scenarios with interventions

#Clear R environment
remove(list = ls())

#packages
require(rstudioapi)
require(lubridate)

#set path
current_path <- getActiveDocumentContext()$path 
head_directory <- gsub(x=current_path, pattern="covid-19-immunity.*", replacement = "covid-19-immunity/")
setwd(head_directory)

#Load UK contact matrix and population data
load("data/BBC_contact_matrices_population_vector.RData")
load("data/uk_population_size_by_age_2018.RData")

#Source model functions
source("age_structured_models/SEIIRRS-discrete-age-interventions-functions.R")

#parameters
R0 <- 2.8                   #Basic reproduction number (in the absense of interventions) [From Petra]
sigma_recip <- 4.5          #Average duration of latent period (days)     [From He, X., et al. Temporal dynamics in viral shedding and transmissibility of COVID-19. Nat Med (2020). https://doi.org/10.1038/s41591-020-0869-5]
gamma_recip <- 3.07         #Average duration of infectiousness (days)    [From He, X., et al.]
omega_recip <- 90           #Average duration of immunity (days) [user specified]
m <- 4                      #Shape paremeter, latent period
n <- 2                      #Shape parameter, infectious period
o <- 2                      #Shape parameter, immune period
dt <- 1                     #time period (days)
days <- 365                 #number of days 
time <- seq(1, days, dt)    #time vector
Rt <- 0.8 #Effective reproduction number during intervention
intervention <- c(0.8,0.3,0.1,0.2)

#start date
start_date <- as.Date("2020-02-28")
end_date <- start_date+(days-1)
dates <- seq.Date(from=start_date, to=end_date, by="1 day")

#Get total UK population
total_pop <- sum(uk.pop.2018.count$total)
p_age <- uk.pop.2018.count$total/total_pop

#Obtain overall contact matrix before and after interventions
C <- make.intervention.matrix(BBC_contact_matrix, intervention = c(1,1,1,1))

#Check Rt & beta after interventions
C_i <- make.intervention.matrix(BBC_contact_matrix, intervention = c(0.8,0.3,0.1,0.2), R0=Rt)
beta_Rt <- get_beta(C=C_i, p_age=p_age, gamma=1/gamma_recip, R0=Rt)

#initial infected vector; 10 individuals per working age group (20-59); total = 80
I_init <- c(0,0,0,0,10,10,10,10,10,10,10,10,0,0,0)

#Proportion of each age group that becomes hospitalised (Rh)
p_hospitalised <- c(0, 0, 0.000408, 0.000408, 0.0104, 0.0104, 0.0343, 0.0343, 0.0425, 0.0425, 0.0816, 0.0816, 0.118, 0.118, 0.175)

#Rho - proportion of infected symtomatic by age group
rho <- c(rep(0.25, 4),rep(0.5,11))

#SEIRRS age structure model
out <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, 
                            immune_mean_1=0, immune_mean_2=0,
                            latent_shape=4, infectious_shape=n, immune_shape=o, dt=1, days=days, C=C, 
                            total_population=total_pop, p_age=p_age, I_init=I_init, Rt = Rt, rho=rho,
                            p_hospitalised=p_hospitalised,
                            asymtomatic_relative_infectiousness=0.5, children_relative_susceptibility=0.4,
                            intervention = intervention, threshold = 10000, t_intervention = 60)

#plot total S, I, R over time (sum over age classes)
S <- apply(out[,,"S"],1,sum)
E <- apply(out[,,"E"],1,sum)
Is <- apply(out[,,"Is"],1,sum)
Ia <- apply(out[,,"Is"],1,sum)
Rh <- apply(out[,,"Rh"],1,sum)
Rnh <- apply(out[,,"Rnh"],1,sum)
R <- Rh + Rnh

plot(S~dates, type="l", ylim=c(0,total_pop), xlab="Time (days)",
     ylab="Population", cex.lab=1.4, cex.axis=1.4, lty=1, lwd=2, col="darkblue")
lines(Is~dates, lty=1, lwd=2, col="darkred")
lines(Ia~dates, lty=1, lwd=2, col="darkred")
lines(R~dates, lty=1, lwd=2, col="darkgreen")

#plot I for seperate age classes
plot(out[,1,"Is"]~dates, type="l", ylim=c(0,250000))
for(k in 2:15){lines(out[,k,"Is"]~dates)}

plot(out[,1,"Ia"]~dates, type="l", ylim=c(0,250000))
for(k in 2:15){lines(out[,k,"Ia"]~dates)}
