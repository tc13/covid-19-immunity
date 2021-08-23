## Vaccine model - multiple scenarios for SPI-M

###############
## Scenarios ##
###############
#R0 =1.5 
V9.R015 <- SEIIRRS_intervention( #vaccine immunity=6 months
  #Mean durations
  latent_mean=4.5, infectious_mean=3.1, immune_mean_1=270, immune_mean_2=365,
  #Shape parameters
  latent_shape=4, infectious_shape=3, immune_shape=3,
  #Vaccination
  vaccine_immune=270, vaccine_shape=10, 
  vaccine_start_days=vaccine_start_days,  #start days by age group
  vaccine_rates_daily=vaccine_rates_daily, #daily rates per age group
  vaccination_length=vaccination_length,  #duration of vaccination period in days
  #Time steps
  dt=dt, days=365*3, 
  #Population
  p_age=p_age, total_population=total_pop, BBC_contact_matrix=BBC_contact_matrix, 
  #Probabilities
  children_relative_susceptibility=0.4, 
  asymtomatic_relative_infectiousness=0.2,
  p_hospitalised=p_hosp, phi=phi,
  #Intervention params
  R0=R0, start_days=start_days, durations=as.numeric(durations), 
  R_values=R_values, 
  contact_interventions=contact_interventions,
  #Initial values
  I_init=I_init, E_init=rep(0,15))

V12.R015 <- SEIIRRS_intervention( #vaccine immunity=12 months
  #Mean durations
  latent_mean=4.5, infectious_mean=3.1, immune_mean_1=270, immune_mean_2=365,
  #Shape parameters
  latent_shape=4, infectious_shape=3, immune_shape=3,
  #Vaccination
  vaccine_immune=365, vaccine_shape=10, 
  vaccine_start_days=vaccine_start_days,  #start days by age group
  vaccine_rates_daily=vaccine_rates_daily, #daily rates per age group
  vaccination_length=vaccination_length,  #duration of vaccination period in days
  #Time steps
  dt=dt, days=365*3, 
  #Population
  p_age=p_age, total_population=total_pop, BBC_contact_matrix=BBC_contact_matrix, 
  #Probabilities
  children_relative_susceptibility=0.4, 
  asymtomatic_relative_infectiousness=0.2,
  p_hospitalised=p_hosp, phi=phi,
  #Intervention params
  R0=R0, start_days=start_days, durations=as.numeric(durations), 
  R_values=c(R_values[1:7],  1.5), 
  contact_interventions=contact_interventions,
  #Initial values
  I_init=I_init, E_init=rep(0,15))

V18.R015 <- SEIIRRS_intervention( #vaccine immunity=18 months
  #Mean durations
  latent_mean=4.5, infectious_mean=3.1, immune_mean_1=270, immune_mean_2=365,
  #Shape parameters
  latent_shape=4, infectious_shape=3, immune_shape=3,
  #Vaccination
  vaccine_immune=540, vaccine_shape=10, 
  vaccine_start_days=vaccine_start_days,  #start days by age group
  vaccine_rates_daily=vaccine_rates_daily, #daily rates per age group
  vaccination_length=vaccination_length,  #duration of vaccination period in days
  #Time steps
  dt=dt, days=365*3, 
  #Population
  p_age=p_age, total_population=total_pop, BBC_contact_matrix=BBC_contact_matrix, 
  #Probabilities
  children_relative_susceptibility=0.4, 
  asymtomatic_relative_infectiousness=0.2,
  p_hospitalised=p_hosp, phi=phi,
  #Intervention params
  R0=R0, start_days=start_days, durations=as.numeric(durations), 
  R_values=c(R_values[1:7],  1.5), 
  contact_interventions=contact_interventions,
  #Initial values
  I_init=I_init, E_init=rep(0,15))

#R0 =2.0 
V9.R02 <- SEIIRRS_intervention( #vaccine immunity=6 months
  #Mean durations
  latent_mean=4.5, infectious_mean=3.1, immune_mean_1=270, immune_mean_2=365,
  #Shape parameters
  latent_shape=4, infectious_shape=3, immune_shape=3,
  #Vaccination
  vaccine_immune=270, vaccine_shape=10, 
  vaccine_start_days=vaccine_start_days,  #start days by age group
  vaccine_rates_daily=vaccine_rates_daily, #daily rates per age group
  vaccination_length=vaccination_length,  #duration of vaccination period in days
  #Time steps
  dt=dt, days=365*3, 
  #Population
  p_age=p_age, total_population=total_pop, BBC_contact_matrix=BBC_contact_matrix, 
  #Probabilities
  children_relative_susceptibility=0.4, 
  asymtomatic_relative_infectiousness=0.2,
  p_hospitalised=p_hosp, phi=phi,
  #Intervention params
  R0=R0, start_days=start_days, durations=as.numeric(durations), 
  R_values=c(R_values[1:7],  2.0), 
  contact_interventions=contact_interventions,
  #Initial values
  I_init=I_init, E_init=rep(0,15))

V12.R02 <- SEIIRRS_intervention( #vaccine immunity=12 months
  #Mean durations
  latent_mean=4.5, infectious_mean=3.1, immune_mean_1=270, immune_mean_2=365,
  #Shape parameters
  latent_shape=4, infectious_shape=3, immune_shape=3,
  #Vaccination
  vaccine_immune=365, vaccine_shape=10, 
  vaccine_start_days=vaccine_start_days,  #start days by age group
  vaccine_rates_daily=vaccine_rates_daily, #daily rates per age group
  vaccination_length=vaccination_length,  #duration of vaccination period in days
  #Time steps
  dt=dt, days=365*3, 
  #Population
  p_age=p_age, total_population=total_pop, BBC_contact_matrix=BBC_contact_matrix, 
  #Probabilities
  children_relative_susceptibility=0.4, 
  asymtomatic_relative_infectiousness=0.2,
  p_hospitalised=p_hosp, phi=phi,
  #Intervention params
  R0=R0, start_days=start_days, durations=as.numeric(durations), 
  R_values=c(R_values[1:7],  2), 
  contact_interventions=contact_interventions,
  #Initial values
  I_init=I_init, E_init=rep(0,15))

V18.R02 <- SEIIRRS_intervention( #vaccine immunity=18 months
  #Mean durations
  latent_mean=4.5, infectious_mean=3.1, immune_mean_1=270, immune_mean_2=365,
  #Shape parameters
  latent_shape=4, infectious_shape=3, immune_shape=3,
  #Vaccination
  vaccine_immune=540, vaccine_shape=10, 
  vaccine_start_days=vaccine_start_days,  #start days by age group
  vaccine_rates_daily=vaccine_rates_daily, #daily rates per age group
  vaccination_length=vaccination_length,  #duration of vaccination period in days
  #Time steps
  dt=dt, days=365*3, 
  #Population
  p_age=p_age, total_population=total_pop, BBC_contact_matrix=BBC_contact_matrix, 
  #Probabilities
  children_relative_susceptibility=0.4, 
  asymtomatic_relative_infectiousness=0.2,
  p_hospitalised=p_hosp, phi=phi,
  #Intervention params
  R0=R0, start_days=start_days, durations=as.numeric(durations), 
  R_values=c(R_values[1:7],  2), 
  contact_interventions=contact_interventions,
  #Initial values
  I_init=I_init, E_init=rep(0,15))

#R0 =2.5 
V9.R025 <- SEIIRRS_intervention( #vaccine immunity=6 months
  #Mean durations
  latent_mean=4.5, infectious_mean=3.1, immune_mean_1=270, immune_mean_2=365,
  #Shape parameters
  latent_shape=4, infectious_shape=3, immune_shape=3,
  #Vaccination
  vaccine_immune=270, vaccine_shape=10, 
  vaccine_start_days=vaccine_start_days,  #start days by age group
  vaccine_rates_daily=vaccine_rates_daily, #daily rates per age group
  vaccination_length=vaccination_length,  #duration of vaccination period in days
  #Time steps
  dt=dt, days=365*3, 
  #Population
  p_age=p_age, total_population=total_pop, BBC_contact_matrix=BBC_contact_matrix, 
  #Probabilities
  children_relative_susceptibility=0.4, 
  asymtomatic_relative_infectiousness=0.2,
  p_hospitalised=p_hosp, phi=phi,
  #Intervention params
  R0=R0, start_days=start_days, durations=as.numeric(durations), 
  R_values=c(R_values[1:7],  2.5), 
  contact_interventions=contact_interventions,
  #Initial values
  I_init=I_init, E_init=rep(0,15))

V12.R025 <- SEIIRRS_intervention( #vaccine immunity=12 months
  #Mean durations
  latent_mean=4.5, infectious_mean=3.1, immune_mean_1=270, immune_mean_2=365,
  #Shape parameters
  latent_shape=4, infectious_shape=3, immune_shape=3,
  #Vaccination
  vaccine_immune=365, vaccine_shape=10, 
  vaccine_start_days=vaccine_start_days,   #start days by age group
  vaccine_rates_daily=vaccine_rates_daily, #daily rates per age group
  vaccination_length=vaccination_length,   #duration of vaccination period in days
  #Time steps
  dt=dt, days=365*3, 
  #Population
  p_age=p_age, total_population=total_pop, BBC_contact_matrix=BBC_contact_matrix, 
  #Probabilities
  children_relative_susceptibility=0.4, 
  asymtomatic_relative_infectiousness=0.2,
  p_hospitalised=p_hosp, phi=phi,
  #Intervention params
  R0=R0, start_days=start_days, durations=as.numeric(durations), 
  R_values=c(R_values[1:7],  2.5), 
  contact_interventions=contact_interventions,
  #Initial values
  I_init=I_init, E_init=rep(0,15))

V18.R025 <- SEIIRRS_intervention( #vaccine immunity=18 months
  #Mean durations
  latent_mean=4.5, infectious_mean=3.1, immune_mean_1=270, immune_mean_2=365,
  #Shape parameters
  latent_shape=4, infectious_shape=3, immune_shape=3,
  #Vaccination
  vaccine_immune=540, vaccine_shape=10, 
  vaccine_start_days=vaccine_start_days,  #start days by age group
  vaccine_rates_daily=vaccine_rates_daily, #daily rates per age group
  vaccination_length=vaccination_length,  #duration of vaccination period in days
  #Time steps
  dt=dt, days=365*3, 
  #Population
  p_age=p_age, total_population=total_pop, BBC_contact_matrix=BBC_contact_matrix, 
  #Probabilities
  children_relative_susceptibility=0.4, 
  asymtomatic_relative_infectiousness=0.2,
  p_hospitalised=p_hosp, phi=phi,
  #Intervention params
  R0=R0, start_days=start_days, durations=as.numeric(durations), 
  R_values=c(R_values[1:7],  2.5), 
  contact_interventions=contact_interventions,
  #Initial values
  I_init=I_init, E_init=rep(0,15))

#R0 =3 
V9.R03 <- SEIIRRS_intervention( #vaccine immunity=6 months
  #Mean durations
  latent_mean=4.5, infectious_mean=3.1, immune_mean_1=270, immune_mean_2=365,
  #Shape parameters
  latent_shape=4, infectious_shape=3, immune_shape=3,
  #Vaccination
  vaccine_immune=270, vaccine_shape=10, 
  vaccine_start_days=vaccine_start_days,  #start days by age group
  vaccine_rates_daily=vaccine_rates_daily, #daily rates per age group
  vaccination_length=vaccination_length,  #duration of vaccination period in days
  #Time steps
  dt=dt, days=365*3, 
  #Population
  p_age=p_age, total_population=total_pop, BBC_contact_matrix=BBC_contact_matrix, 
  #Probabilities
  children_relative_susceptibility=0.4, 
  asymtomatic_relative_infectiousness=0.2,
  p_hospitalised=p_hosp, phi=phi,
  #Intervention params
  R0=R0, start_days=start_days, durations=as.numeric(durations), 
  R_values=c(R_values[1:7],  3), 
  contact_interventions=contact_interventions,
  #Initial values
  I_init=I_init, E_init=rep(0,15))

V12.R03 <- SEIIRRS_intervention( #vaccine immunity=12 months
  #Mean durations
  latent_mean=4.5, infectious_mean=3.1, immune_mean_1=270, immune_mean_2=365,
  #Shape parameters
  latent_shape=4, infectious_shape=3, immune_shape=3,
  #Vaccination
  vaccine_immune=365, vaccine_shape=10, 
  vaccine_start_days=vaccine_start_days,  #start days by age group
  vaccine_rates_daily=vaccine_rates_daily, #daily rates per age group
  vaccination_length=vaccination_length,  #duration of vaccination period in days
  #Time steps
  dt=dt, days=365*3, 
  #Population
  p_age=p_age, total_population=total_pop, BBC_contact_matrix=BBC_contact_matrix, 
  #Probabilities
  children_relative_susceptibility=0.4, 
  asymtomatic_relative_infectiousness=0.2,
  p_hospitalised=p_hosp, phi=phi,
  #Intervention params
  R0=R0, start_days=start_days, durations=as.numeric(durations), 
  R_values=c(R_values[1:7],  3), 
  contact_interventions=contact_interventions,
  #Initial values
  I_init=I_init, E_init=rep(0,15))

V18.R03 <- SEIIRRS_intervention( #vaccine immunity=18 months
  #Mean durations
  latent_mean=4.5, infectious_mean=3.1, immune_mean_1=270, immune_mean_2=365,
  #Shape parameters
  latent_shape=4, infectious_shape=3, immune_shape=3,
  #Vaccination
  vaccine_immune=540, vaccine_shape=10, 
  vaccine_start_days=vaccine_start_days,  #start days by age group
  vaccine_rates_daily=vaccine_rates_daily, #daily rates per age group
  vaccination_length=vaccination_length,  #duration of vaccination period in days
  #Time steps
  dt=dt, days=365*3, 
  #Population
  p_age=p_age, total_population=total_pop, BBC_contact_matrix=BBC_contact_matrix, 
  #Probabilities
  children_relative_susceptibility=0.4, 
  asymtomatic_relative_infectiousness=0.2,
  p_hospitalised=p_hosp, phi=phi,
  #Intervention params
  R0=R0, start_days=start_days, durations=as.numeric(durations), 
  R_values=c(R_values[1:7],  3), 
  contact_interventions=contact_interventions,
  #Initial values
  I_init=I_init, E_init=rep(0,15))
