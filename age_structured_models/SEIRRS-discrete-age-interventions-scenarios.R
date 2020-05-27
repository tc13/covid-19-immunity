### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SEIIRRS discrete time model, age structed for UK population, gamma (Erlang) distributed waiting times, two immunity classes
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
rho <- c(rep(0.25, 4),rep(0.5,11))  #proportion of infected symtomatic by age group
crs <- 0.4                          #children relative susceptibility
ari <- 0.5                          #asymtomatic relative infectiousness

############################
## Intervention variables ##
############################
Rt_partial <- 1.6       #Effective reproduction number before and after main intervention
Rt_full <- 0.9          #Effective reproduction number during main intervention
intervention_partial <- c(home=1,work=0.8,school=0.85, other=0.75)    #Scaling of contact matrix during partial intervention period
intervention_full <- c(home=0.8,work=0.3,school=0.1,other=0.2)        #Scaling of contact matrix during full intervention period
threshold <- 50000      #Number of symtomatic cases required to trigger partial intervention
tpi <- 10               #Length of partial intevention in days
tfi <- 60               #Length of full intervention in days

#####################
## Time  variables ##
#####################
dt <- 1                    #time period (days)
days <- 365                 #number of days 
time <- seq(1, days, dt)    #time vector

######################################
## UK Population & Contact Matrixes ##
######################################
total_pop <- sum(uk.pop.2018.count$total)   #Get total UK population
p_age <- uk.pop.2018.count$total/total_pop  #Proportion of population in each age group
I_init <- c(0,0,0,0,20,20,20,20,20,20,20,20,0,0,0)  #Initial numbers infected: 20 individuals per working age group (20-59); total = 160

###################
## Run Scenarios ##
###################

#Scenario 1: hospitalised & non-hospitalised patients have permanent immunity
S1 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=1, infectious_shape=1, #Natural history parameters
                          p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, rho=rho, #Natural history parameters
                          immune_mean_1=0, immune_mean_2=0, immune_shape=1, #scenario specific immune duration
                          dt=dt, days=days, #times
                          BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                          Rt_partial=Rt_partial, intervention_partial=intervention_partial, t_partial_intervention=tpi, #Partial intervention params
                          Rt_full=Rt_full, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                          threshold = threshold)

#Scenario 2: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 365 days
S2 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                           p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, rho=rho, #Natural history parameters
                           immune_mean_1=0, immune_mean_2=365, immune_shape=2, #scenario specific immune duration
                           dt=dt, days=days, #times
                           BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                           Rt_partial=Rt_partial, intervention_partial=intervention_partial, t_partial_intervention=tpi, #Partial intervention params
                           Rt_full=Rt_full, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                           threshold = threshold)

#Scenario 3: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 180 days
S3 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                           p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, rho=rho, #Natural history parameters
                           immune_mean_1=0, immune_mean_2=180, immune_shape=2, #scenario specific immune duration
                           dt=dt, days=days, #times
                           BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                           Rt_partial=Rt_partial, intervention_partial=intervention_partial, t_partial_intervention=tpi, #Partial intervention params
                           Rt_full=Rt_full, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                           threshold = threshold)

#Scenario 4: hospitalised patients have immunity for mean 365 days, non-hospitalised patients have immunity for mean 90 days
S4 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                           p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, rho=rho, #Natural history parameters
                           immune_mean_1=365, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                           dt=dt, days=days, #times
                           BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                           Rt_partial=Rt_partial, intervention_partial=intervention_partial, t_partial_intervention=tpi, #Partial intervention params
                           Rt_full=Rt_full, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                           threshold = threshold)

#Plot and analyse scenario output (this could be a seperate script)
require(ggplot2)
require(reshape2)

#ggplot theme set
theme_set(theme_bw(base_size = 28, base_family = "Times"))
          
#Transform arrays in data.frames
S1.df <- melt(S1)
S2.df <- melt(S2)
S3.df <- melt(S3)
S4.df <- melt(S4)

#change column names
column.names <- c("time", "age_group", "state", "count")
colnames(S1.df) <- column.names
colnames(S2.df) <- column.names
colnames(S3.df) <- column.names
colnames(S4.df) <- column.names

#Date variables
start_date <- as.Date("2020-01-04")   #Simulation start date
end_date <- start_date+(days-1)       #Simulation end date
dates <- seq.Date(from=start_date, to=end_date, by="1 day")   #Date sequence
dates <- ymd(dates) #convert to lubridate format
S1.df$date <- dates
S2.df$date <- dates
S3.df$date <- dates
S4.df$date <- dates

#test
S <- apply(S1[,,"S"],1,sum)
I <- apply(S1[,,"I"],1,sum)
Is <- apply(S1[,,"Is"],1,sum)
Rh <- apply(S1[,,"Rh"],1,sum)
plot(S~dates)
plot(I~dates)
plot(Is~dates)
plot(Rh~dates)
R <- apply(S1[,,"R"],1,sum)

dates[min(which(Is>threshold))] #date initital threshold breached
begin_lockdown <- dates[min(which(Is>threshold))] + tpi #date main lockdown begins
end_lockdown <- dates[min(which(Is>threshold))] + tpi + tfi #date end of lockdown
which(dates==begin_lockdown)
which(dates==end_lockdown)
R[80]
R[150]/total_pop #fraction immune (but could be higher in older age groups)

#Plot S1
S1.plt <- ggplot(data=S1.df, aes(x=date, y=count, shape=age_group, colour=state))
S1.plt + geom_line() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())


######
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
