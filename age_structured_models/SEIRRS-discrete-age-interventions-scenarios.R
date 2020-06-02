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
R0 <- 1                 #Basic reproduction number (in the absense of interventions) [From Petra]
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
Rt_partial <- 1       #Effective reproduction number before and after main intervention
Rt_full <- 1          #Effective reproduction number during main intervention
intervention_partial <- c(home=1,work=0.8,school=0.85, other=0.75)    #Scaling of contact matrix during partial intervention period
intervention_full <- c(home=0.8,work=0.3,school=0.1,other=0.2)        #Scaling of contact matrix during full intervention period
threshold <- 85000      #Number of symtomatic cases required to trigger partial intervention
tpi <- 14               #Length of partial intevention in days
tfi <- 30               #Length of full intervention in days
lockdown_day <- 50     #If triggering intervention by number of days, days after simulation starts

#####################
## Time  variables ##
#####################
dt <- 1                   #time period (days)
days <- 365#*2             #number of days 
time <- seq(1, days, dt)    #time vector

######################################
## UK Population & Contact Matrixes ##
######################################
total_pop <- sum(uk.pop.2018.count$total)   #Get total UK population
p_age <- uk.pop.2018.count$total/total_pop  #Proportion of population in each age group
I_init <- c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)  #Initial numbers infected: 25 individuals per working age group (20-59); total = 200
#I_init <- c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0)
###################
## Run Scenarios ##
###################
C = make.intervention.matrix(matrixlistin = BBC_contact_matrix, intervention=c(1,1,1,1))
get_beta(C=C, R0=R0, mean_infectious = gamma_recip, susceptibility = c(rep(crs, 4), rep(1, 11)),
         prop_asymtomatic = phi, asymtomatic_infectiousness = ari, p_age=p_age)

#Scenario 1: hospitalised & non-hospitalised patients have permanent immunity
S1 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                          p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                          immune_mean_1=0, immune_mean_2=0, immune_shape=2, #scenario specific immune duration
                          dt=1, days=days, #times
                          BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                          Rt_partial=Rt_partial, intervention_partial=intervention_partial, t_partial_intervention=tpi, #Partial intervention params
                          Rt_full=Rt_full, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                          trigger="none", lockdown_day=lockdown_day+100)

#Scenario 2: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 365 days
S2 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                           p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                           immune_mean_1=0, immune_mean_2=365, immune_shape=2, #scenario specific immune duration
                           dt=dt, days=days, #times
                           BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                           Rt_partial=Rt_partial, intervention_partial=intervention_partial, t_partial_intervention=tpi, #Partial intervention params
                           Rt_full=Rt_full, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                           trigger="days", lockdown_day=lockdown_day)

#Scenario 3: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 180 days
S3 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                           p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                           immune_mean_1=0, immune_mean_2=180, immune_shape=2, #scenario specific immune duration
                           dt=dt, days=days, #times
                           BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                           Rt_partial=Rt_partial, intervention_partial=intervention_partial, t_partial_intervention=tpi, #Partial intervention params
                           Rt_full=Rt_full, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                           trigger="days", lockdown_day=lockdown_day)

#Scenario 4: hospitalised patients have immunity for mean 365 days, non-hospitalised patients have immunity for mean 90 days
S4 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                           p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                           immune_mean_1=365, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                           dt=dt, days=days, #times
                           BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                           Rt_partial=Rt_partial, intervention_partial=intervention_partial, t_partial_intervention=tpi, #Partial intervention params
                           Rt_full=Rt_full, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                           trigger="days", lockdown_day=lockdown_day)

#Date variables
date_UK_lockdown <- "2020-03-23"
start_date <- as.Date(date_UK_lockdown)-lockdown_day
#dates <- seq.Date(from=start_date, to=start_date+(days-1), by="1/4 day")
hour_period <- 24*dt
dates <- seq(ymd(start_date),ymd(start_date+(days-1)), by = as.difftime(hours(hour_period)))

#New infection on lockdown day
new_infections <- apply(S1[,,"new_infections"],1,sum)
new_infections_adults <- apply(S1[,c(5:15),"new_infections"],1,sum)
new_infections[lockdown_day]

#cumulative infections
start_cum_date <- "2020-02-16"
which(dates==start_cum_date)
sum(new_infections[which(dates==start_cum_date):which(dates==date_UK_lockdown)])
sum(new_infections_adults[which(dates==start_cum_date):which(dates==date_UK_lockdown)])

#Recovered (ages 19+)
recovered_19plus <- apply(S1[,c(5:15),"R"],1,sum)
end_lockdown <- as.Date(date_UK_lockdown)+(tpi+tfi) 
which(dates==end_lockdown)
recovered_19plus[which(dates==end_lockdown)]/sum(uk.pop.2018.count$total[5:15])
sum(new_infections[1:which(dates==end_lockdown)])

recovered_all <- apply(S1[,,"R"],1,sum)
recovered_all[which(dates==end_lockdown)]/total_pop

#Base R plotting
I_S1 <- apply(S1[,,"I"],1,sum)
I_S2 <- apply(S2[,,"I"],1,sum)
I_S3 <- apply(S3[,,"I"],1,sum)
I_S4 <- apply(S4[,,"I"],1,sum)

plot(I_S1~dates, type="l", col="darkgrey",lwd=2)
lines(I_S2~dates, col="darkblue", lwd=2)
lines(I_S3~dates, col="darkgreen", lwd=2)
lines(I_S4~dates, col="darkred", lwd=2)

R_S1 <- apply(S1[,,"R"],1,sum)
R_S2 <- apply(S2[,,"R"],1,sum)
R_S3 <- apply(S3[,,"R"],1,sum)
R_S4 <- apply(S4[,,"R"],1,sum)

plot(R_S1~dates, type="l", col="darkgrey",lwd=2)
lines(R_S2~dates, col="darkblue", lwd=2)
lines(R_S3~dates, col="darkgreen", lwd=2)
lines(R_S4~dates, col="darkred", lwd=2)

##ggplotting
#Plot and analyse scenario output (this could be a seperate script)
require(ggplot2)
require(reshape2)
require(scales)

#ggplot theme set
theme_set(theme_bw(base_size = 28, base_family = "Times"))

#Function to make axis numbers in scientific notation
scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

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

dates <- ymd(dates) #convert to lubridate format
S1.df$date <- dates
S2.df$date <- dates
S3.df$date <- dates
S4.df$date <- dates

S1.df$scenario <- "1"
S2.df$scenario <- "2"
S3.df$scenario <- "3"
S4.df$scenario <- "4"

scenarios <- rbind(S1.df, S2.df, S3.df, S4.df)

##Rt over time
Rt <- ifelse(dates<=dates[lockdown_day], R0, ifelse(dates>dates[(lockdown_day+tpi)]&dates<=dates[(lockdown_day+tpi+tfi)], Rt_full, Rt_partial))
Rt.df <- data.frame(date=dates, Rt=Rt)
#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Figure3-Rt.pdf", height=7, width=11)
Rt.plt <- ggplot(data=Rt.df, aes(x=date, y=Rt))
Rt.plt + geom_line(size=0.9) + theme(panel.grid.minor=element_blank(),
                                     panel.grid.major=element_blank())+
  scale_x_date(name="Date", date_labels ="%b-%y", 
               breaks=seq.Date(from=as.Date("2020-01-01"), to=as.Date("2022-01-01"), by="4 months"))+
  scale_y_continuous(name="Effective Reproduction Number", limits=c(0,3),
                     breaks=c(0,1,2,3), labels=c(0,1,2,3), expand=c(0,0))+
  geom_vline(xintercept = c(as.Date(date_UK_lockdown), as.Date(date_UK_lockdown)+tpi+tfi), 
             size=0.9, linetype="dashed")
#dev.off()

#Plot infecteds (symtomatic & asymtomatic) over time
pal <- c("#1b9e77", "#d95f02", "#8856a7", "#2166ac")
#After summing over age groups
I_scenarios.wide <- data.frame(date=dates, S1=I_S1, S2=I_S2, S3=I_S3, S4=I_S4)
I_scenarios <- melt(I_scenarios.wide, id.vars = "date", variable.name="scenario", value.name="count")
y.vals <- c(0, 1E5, 2E5, 3E5)
y.labs <- scientific(y.vals)

#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Figure3-Infected.pdf", height=7, width=11)
I.plt <- ggplot(data=I_scenarios, aes(x=date, y=count, colour=scenario))
I.plt + geom_line(size=0.9) + theme(panel.grid.minor=element_blank(),
                            panel.grid.major=element_blank(),
                            legend.position = c(0.8, 0.8))+
  scale_x_date(name="Date", date_labels ="%b-%y", 
               breaks=seq.Date(from=as.Date("2020-01-01"), to=as.Date("2022-01-01"), by="4 months")) + 
  scale_y_continuous(name = "Number Infected")+ #,labels=y.labs, limits=c(0, 3.2E5), expand = c(0,0)) + 
  scale_color_manual(name="Model Scenario", values=pal, labels=c("1","2","3","4"))+
  geom_vline(xintercept = c(as.Date(date_UK_lockdown), as.Date(date_UK_lockdown)+tpi+tfi), 
             size=0.9, linetype="dashed")
#dev.off()

#Plot Recovered (hospitalised & Non-hospitalised) over time
#After summing over age groups
R_scenarios.wide <- data.frame(date=dates, S1=R_S1, S2=R_S2, S3=R_S3, S4=R_S4)
R_scenarios <- melt(R_scenarios.wide, id.vars = "date", variable.name="scenario", value.name="count")
y.vals <- c(0, 1E6, 2E6, 3E6, 4E6, 5E6)
y.labs <- scientific(y.vals)

#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Figure3-Recovered.pdf", height=7, width=11)
R.plt <- ggplot(data=R_scenarios, aes(x=date, y=count, colour=scenario))
R.plt + geom_line(size=0.9) + theme(panel.grid.minor=element_blank(),
                                    panel.grid.major=element_blank(),
                                    legend.position = "none")+
  scale_x_date(name="Date", date_labels ="%b-%y", 
               breaks=seq.Date(from=as.Date("2020-01-01"), to=as.Date("2022-01-01"), by="4 months")) + 
  scale_y_continuous(name = "Number Immune", labels=y.labs, limits=c(0, 5.2E6), expand = c(0,0)) + 
  scale_color_manual(values=pal, labels=c("1","2","3","4"))+
  geom_vline(xintercept = c(as.Date(date_UK_lockdown), as.Date(date_UK_lockdown)+tpi+tfi), 
             size=0.9, linetype="dashed")
#dev.off()

#Daily new hospitalisations
#Base R plotting
H_S1 <- apply(S1[,,"new_hospitalisations"],1,sum)
H_S2 <- apply(S2[,,"new_hospitalisations"],1,sum)
H_S3 <- apply(S3[,,"new_hospitalisations"],1,sum)
H_S4 <- apply(S4[,,"new_hospitalisations"],1,sum)

H_scenarios.wide <- data.frame(date=dates, S1=H_S1, S2=H_S2, S3=H_S3, S4=H_S4)
H_scenarios <- melt(H_scenarios.wide, id.vars = "date", variable.name="scenario", value.name="count")
#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Figure3-hospitalisations.pdf", height=7, width=11)
H.plt <- ggplot(data=H_scenarios, aes(x=date, y=count, colour=scenario))
H.plt + geom_line(size=0.9) + theme(panel.grid.minor=element_blank(),
                                    panel.grid.major=element_blank(),
                                    legend.position = "none")+
  scale_x_date(name="Date", date_labels ="%b-%y", 
               breaks=seq.Date(from=as.Date("2020-01-01"), to=as.Date("2022-01-01"), by="4 months")) + 
  scale_y_continuous(name = "Daily Hospitalisations", limits=c(0, 8700), expand = c(0,0)) + 
  scale_color_manual(values=pal, labels=c("1","2","3","4"))+
  geom_vline(xintercept = c(as.Date(date_UK_lockdown), as.Date(date_UK_lockdown)+tpi+tfi), 
             size=0.9, linetype="dashed")
#dev.off()
