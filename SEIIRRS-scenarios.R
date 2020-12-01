### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SEIIRRS discrete time model, age structed for UK population, gamma (Erlang) distributed waiting times, two immunity classes
### Model Scenarios with interventions

###################
## Run Scenarios ##
###################
#Rt post lockdown = 0.9
#Scenario 1: hospitalised & non-hospitalised patients have permanent immunity
S1.R09 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                          p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                          immune_mean_1=0, immune_mean_2=0, immune_shape=2, #scenario specific immune duration
                          days=days, #times
                          BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                          Rt_post_lockdown=0.9, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                          Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                          trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 2: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 365 days
S2.R09 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                           p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                           immune_mean_1=0, immune_mean_2=365, immune_shape=2, #scenario specific immune duration
                           days=days, #times
                           BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                           Rt_post_lockdown=0.9, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                           Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                           trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 3: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 180 days
S3.R09 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                           p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                           immune_mean_1=0, immune_mean_2=180, immune_shape=2, #scenario specific immune duration
                           days=days, #times
                           BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                           Rt_post_lockdown=0.9, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                           Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                           trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 4: hospitalised patients have immunity for mean 365 days, non-hospitalised patients have immunity for mean 90 days
S4.R09 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                           p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                           immune_mean_1=365, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                           days=days, #times
                           BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                           Rt_post_lockdown=0.9, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                           Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                           trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Rt post lockdown = 1.0
#Scenario 1: hospitalised & non-hospitalised patients have permanent immunity
S1.R10 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=0, immune_shape=2, #scenario specific immune duration
                               days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.0, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 2: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 365 days
S2.R10 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=365, immune_shape=2, #scenario specific immune duration
                               days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.0, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 3: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 180 days
S3.R10 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=180, immune_shape=2, #scenario specific immune duration
                               days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.0, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 4: hospitalised patients have immunity for mean 365 days, non-hospitalised patients have immunity for mean 90 days
S4.R10 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=365, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                               days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.0, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Rt post lockdown = 1.1
#Scenario 1: hospitalised & non-hospitalised patients have permanent immunity
S1.R11 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=0, immune_shape=2, #scenario specific immune duration
                               days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.1, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 2: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 365 days
S2.R11 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=365, immune_shape=2, #scenario specific immune duration
                               days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.1, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 3: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 180 days
S3.R11 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=180, immune_shape=2, #scenario specific immune duration
                               days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.1, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 4: hospitalised patients have immunity for mean 365 days, non-hospitalised patients have immunity for mean 90 days
S4.R11 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=365, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                               days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.1, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Rt post lockdown = 1.2
#Scenario 1: hospitalised & non-hospitalised patients have permanent immunity
S1.R12 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=0, immune_shape=2, #scenario specific immune duration
                               days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.2, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 2: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 365 days
S2.R12 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=365, immune_shape=2, #scenario specific immune duration
                               days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.2, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 3: hospitalised patients have permanent immunity, non-hospitalised patients have immunity for mean 180 days
S3.R12 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=0, immune_mean_2=180, immune_shape=2, #scenario specific immune duration
                               days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.2, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

#Scenario 4: hospitalised patients have immunity for mean 365 days, non-hospitalised patients have immunity for mean 90 days
S4.R12 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=365, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                               days=days, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.2, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)
