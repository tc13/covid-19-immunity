### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, July 2020
### SEIIRRS discrete time model, age structed for UK population, gamma (Erlang) distributed waiting times, two immunity classes
### Model parameters

###########################################
## SARS-CoV-2 Natural history parameters ##
###########################################
R0 <- 3.7                    #Basic reproduction number (in the absense of interventions) Flaxman et al. 2020 Nature. https://doi.org/10.1038/s41586-020-2405-7 
sig_mean <- 4.5              #Average duration of latent period (days)     [From He, X., et al. Temporal dynamics in viral shedding and transmissibility of COVID-19. Nat Med (2020). https://doi.org/10.1038/s41591-020-0869-5]
gam_mean <- 3.1              #Average duration of infectiousness (days)    [From He, X., et al.]
m <- 4                       #Shape paremeter, latent period
n <- 2                       #Shape parameter, infectious period
p_hosp <- c(0.001, 0.001, 0.003, 0.003, 0.012, 0.012, 0.032, 0.032, 0.049, 0.049, 0.102, 0.102, 0.166, 0.166, 0.26) #From Ferguson et al. Report 9
phi <- c(rep(0.25, 4), rep(0.15, 11)) #proportion asymtomatic | infection by age group
crs <- 0.5                          #children relative susceptibility
ari <- 0.2                          #asymtomatic relative infectiousness
susceptibility <- c(rep(crs, 4), rep(1, 11))

######################################
## UK Population & Contact Matrixes ##
######################################
load("Documents/covid-19-immunity/data/BBC_contact_matrices_population_vector.RData")
load("Documents/covid-19-immunity/data/uk_population_size_by_age_2018.RData")
total_pop <- 56286961   #Get total England population
england_pop <- c(3299637, 3538206, 3354246, 3090232, 3487863, 3801409, 
                 3807954, 3733642, 3414297, 3715812, 3907461, 3670651, 
                 3111835, 2796740, (2779326+ 1940686+ 1439913+ 879778+ 517273))
p_age <- england_pop/total_pop  #Proportion of population in each age group

#####################
## Time  variables ##
#####################
sim_start_date <- "2020-03-03"
days <- 365*3   #number of days
sim_dates <- seq.Date(from=as.Date(sim_start_date), to=as.Date(sim_start_date)+days,by="day")
sim_end_date <- max(sim_dates)
#lockdown dates
lockdown_1_start <- as.Date("2020-03-23")
lockdown_1_mid <- as.Date("2020-04-18")
lockdown_1_end <- as.Date("2020-07-15")
lockdown_2_start <- as.Date("2020-10-31")
lockdown_2_end <- as.Date("2020-12-02")
lockdown_3_start <- as.Date("2021-01-06")
lockdown_3_end <- as.Date("2021-04-12")
full_reopening <- as.Date("2021-06-21")

###################
## Interventions ##
##################
start_days <- match(c(lockdown_1_start, lockdown_1_mid, lockdown_1_end, lockdown_2_start, lockdown_2_end, lockdown_3_start, lockdown_3_end, full_reopening), sim_dates)
durations <- c((lockdown_1_mid-lockdown_1_start), c(lockdown_1_end-lockdown_1_mid), c(lockdown_2_start-lockdown_1_end), (lockdown_2_end-lockdown_2_start), (lockdown_3_start-lockdown_2_end),(lockdown_3_end-lockdown_3_start),(full_reopening-lockdown_3_end), (sim_end_date-full_reopening))-1
R_values=c(1.05, 0.85, 1.25, 1.15, 1.25, 1.12, 1.25, 1.5)
contact_interventions=list(
                   c(home=0.8,work=0.4,school=0.1,other=0.3),  #lockdown1 start
                   c(home=0.8,work=0.3,school=0.1,other=0.2),  #lockdown 1 later
                   c(home=1,work=0.6,school=0.85, other=0.75), #Summer 2020
                   c(home=0.8,work=0.4,school=0.1,other=0.2),  #lockdown 2
                   c(home=1,work=0.6,school=0.85, other=0.75), #Christmas 2020
                   c(home=0.8,work=0.4,school=0.2, other=0.2),  #lockdown 3
                   c(home=0.8,work=0.5,school=0.85, other=0.4),  #April-June 2021
                   c(home=0.9,work=0.6,school=1, other=0.8))  #After June 2021

#################
## vaccination ##
#################
vaccine_shape <- 2
vaccine_rates_weekly <-  c(rep(0,3),rep(0.133, 7), 0.396, 0.386, 
                           0.456, 0.457, 0.299) #rates of vaccination by age group and time
vaccine_rates_daily <- vaccine_rates_weekly/7
vaccine_start_dates <- c(rep("2021-03-21",3),rep("2021-03-21",7),"2021-03-14","2021-03-07","2021-02-28","2021-02-14","2020-12-06")
vaccine_start_days <- match(as.Date(vaccine_start_dates), sim_dates)
vaccination_length <- 7*16

##################################
##  Dominant eigenvector of NGM ##            
##################################
C <- make.intervention.matrix(BBC_contact_matrix, intervention=c(1,1,1,1))
eigenvec = get_eigenvector(C=C, R0=R0, mean_infectious=gam_mean,
                           susceptibility=susceptibility, prop_asymtomatic=phi,
                           asymtomatic_infectiousness=ari)
eigenvec_norm <- eigenvec/sum(eigenvec) #normalise

####################
## Initial values ##
####################
imported_lineages <- 2968  #from de Plessis et al. 2020. https://doi.org/10.1101/2020.10.23.20218446 
I_init <- imported_lineages*eigenvec_norm #distribute initial cases by age group according to dominant eigenvector

#Index intervention by time
intervention.indx <- c(rep(1, (start_days[1])), 
                    rep(c(2, 3, 4, 5, 6, 7, 8, 9), durations+1))
