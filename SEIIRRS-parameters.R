### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, July 2020
### SEIIRRS discrete time model, age structed for UK population, gamma (Erlang) distributed waiting times, two immunity classes
### Model parameters

###########################################
## SARS-CoV-2 Natural history parameters ##
###########################################
R0 <- 4.0                 #Basic reproduction number (in the absense of interventions) Flaxman et al. 2020 Nature. https://doi.org/10.1038/s41586-020-2405-7 
sigma_recip <- 4.5        #Average duration of latent period (days)     [From He, X., et al. Temporal dynamics in viral shedding and transmissibility of COVID-19. Nat Med (2020). https://doi.org/10.1038/s41591-020-0869-5]
gamma_recip <- 3.1       #Average duration of infectiousness (days)    [From He, X., et al.]
m <- 4                    #Shape paremeter, latent period
n <- 2                    #Shape parameter, infectious period
p_hosp <- c(0.001, 0.001, 0.003, 0.003, 0.012, 0.012, 0.032, 0.032, 0.049, 0.049, 0.102, 0.102, 0.166, 0.166, 0.26) #From Ferguson et al. Report 9
phi <- c(rep(0.75, 4),rep(0.5,11))  #proportion asymtomatic | infection by age group
crs <- 0.4                          #children relative susceptibility
ari <- 0.5                          #asymtomatic relative infectiousness
susceptibility <- c(rep(crs, 4), rep(1, 11))

############################
## Intervention variables ##
############################
Rt_lockdown <- 0.8              #Effective reproduction number during lockdown (full intervention)
Rt_partial <- 1.2               #Effective reproduction number at start of lockdown in tpi period 
intervention_post_lockdown <- c(home=1,work=0.8,school=0.85, other=0.75)    #Scaling of contact matrix during partial intervention period
intervention_full <- c(home=0.8,work=0.3,school=0.1,other=0.2)        #Scaling of contact matrix during full intervention period
tpi <- 21               #Length of partial intevention in days
tfi <- 60               #Length of full intervention in days
lockdown_day <- 20      #If triggering intervention by number of days, days after simulation starts

##################################
##  Dominant eigenvector of NGM ##            
##################################
C <- make.intervention.matrix(BBC_contact_matrix, intervention=c(1,1,1,1))
eigenvec = get_eigenvector(C=C, R0=R0, mean_infectious=gamma_recip,
                susceptibility=susceptibility, prop_asymtomatic=phi,
                asymtomatic_infectiousness=ari)
eigenvec_norm <- eigenvec/sum(eigenvec) #normalise

#####################
## Time  variables ##
#####################
dt <- 0.01                 #time step = 1/100 days
days <- 365*2              #number of days 
imported_lineages <- 2968  #from de Plessis et al. 2020. https://doi.org/10.1101/2020.10.23.20218446 
I_init <- imported_lineages*eigenvec_norm #distribute initial cases by age group according to dominant eigenvector
