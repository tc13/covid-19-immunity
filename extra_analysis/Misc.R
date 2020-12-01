
##############################
## Height of secondary peak ##
##############################

S1.R12.ni <- apply(S1.R12[,,"new_infections"],1,sum)
plot(S1.R12.ni~time)
max(S1.R12.ni[200:400])

S1.R11.ni <- apply(S1.R11[,,"new_infections"],1,sum)
plot(S1.R11.ni~time)

S4.R11.ni <- apply(S4.R11[,,"new_infections"],1,sum)
plot(S4.R11.ni~time)
max(S4.R11.ni[200:600])
which(S4.R11.ni[200:600]==max(S4.R11.ni[200:600]))
dates[440]
S1.R11.ni[440]

#Height of total infected
R12.infected.S1 <- apply(S1.R12[,,"I"],1,sum)
R12.infected.S4 <- apply(S4.R12[,,"I"],1,sum)

plot(R12.infected.S1~time)
plot(R12.infected.S4~time)

max(R12.infected.S1[200:600])
max(R12.infected.S4[200:600])

R11.infected.S1 <- apply(S1.R11[,,"I"],1,sum)
R11.infected.S4 <- apply(S4.R11[,,"I"],1,sum)

plot(R11.infected.S1~time)
plot(R11.infected.S4~time)

max(R11.infected.S1[200:600])
max(R11.infected.S4[200:600])

################################################
#Susceptible and exposed vectors (when Rt=1.2) #
################################################

S_S1 <- apply(S1.R12[,,"S"],1,sum)/total_pop
S_S2 <- apply(S2.R12[,,"S"],1,sum)/total_pop
S_S3 <- apply(S3.R12[,,"S"],1,sum)/total_pop
S_S4 <- apply(S4.R12[,,"S"],1,sum)/total_pop

#find values closest to 1/Rt
equil <- 1/1.2
S1_eq_date = dates[min(which(S_S1<equil))]
S2_eq_date = dates[min(which(S_S2<equil))]
S3_eq_date = dates[min(which(S_S3<equil))]
S4_eq_date = dates[min(which(S_S4<equil))]

par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0))
plot(S_S1~dates, type="l", lwd=3, ylim=c(0.8,1), col="grey40", xaxt="n",cex.lab=2, 
     cex.axis=2, cex.main=2, ylab="Proportion susceptible", xlab="Date",
     main="Proportion susceptible; post-lockdown Rt = 1.2")
lines(S_S2~dates, lwd=3, col="blue4")
lines(S_S3~dates, lwd=3, col="green4")
lines(S_S4~dates, lwd=3, col="red3")
abline(h = 1/1.2, lty=2, lwd=3)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2) #x-axis dates
abline(v=c(S1_eq_date, S2_eq_date, S3_eq_date, S4_eq_date), lwd=2, lty=2, 
       col=c("grey40", "blue4", "green4", "red3"))
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2)

#Infection vectors
I_S1.p <- apply(S1.R12[,,"I"],1,sum)/total_pop
I_S2.p <- apply(S2.R12[,,"I"],1,sum)/total_pop
I_S3.p <- apply(S3.R12[,,"I"],1,sum)/total_pop
I_S4.p <- apply(S4.R12[,,"I"],1,sum)/total_pop

#Rt = 1.2
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(I_S1~dates, type="l", col="grey40",lwd=3, xaxt = "n", yaxt="n", cex.main=2,
     xlab="Date", ylab="Population Infected", cex.lab=2, ylim=c(0,45E4), yaxs="i",
     main="Population infectious; post-lockdown Rt = 1.2")
lines(I_S2~dates, col="blue4", lwd=3)
lines(I_S3~dates, col="green4", lwd=3)
lines(I_S4~dates, col="red3", lwd=3)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2) #x-axis dates
axis(2, y.vals, labels=y.labs, cex.axis=2) #y-axis scientific notation
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2)
abline(v=c(S1_eq_date, S2_eq_date, S3_eq_date, S4_eq_date), lwd=2, lty=2, 
       col=c("grey40", "blue4", "green4", "red3"))

#Plot of I vs S
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(S_S1,I_S1.p, type="l", lwd=3, xlim=c(0.75,1), ylim=c(0,0.0065), col="grey40", 
     cex.lab=2, cex.axis=2, cex.main=2,
     ylab="Proportion Infectious", xlab="Proportion Susceptible", main="Phase-portrait (10 years); post-lockdown Rt = 1.2")
lines(S_S2,I_S2.p, lwd=3, col="blue4")
lines(S_S3,I_S3.p, lwd=3, col="green4")
lines(S_S4,I_S4.p, lwd=3, col="red3")
abline(v = 1/1.2, lty=2, lwd=2)
legend(x=c(0.743,0.829), y=c(0.004,0.006), title="Immunity Scenario",
       legend=c(expression(paste("S1  ", omega^N, " = ", infinity, ",       ", omega^H," = ", infinity)), 
                expression(paste("S2  ", omega^N, " = ", 365^-1, ", ", omega^H," = ", infinity)), 
                expression(paste("S3  ", omega^N, " = ", 180^-1, ", ", omega^H," = ", infinity)), 
                expression(paste("S4  ", omega^N, " = ", 90^-1, ",   ", omega^H," = ", 365^-1))), 
       lty=1, lwd=c(3.2,3,3,3.2), cex=1.24,
       col=c("grey40","blue4","green4", "red3"))

#incidence of new cases
New_S1 <- apply(S1.R12[,,"new_infections"],1,sum)
New_S2 <- apply(S2.R12[,,"new_infections"],1,sum)
New_S3 <- apply(S3.R12[,,"new_infections"],1,sum)
New_S4 <- apply(S4.R12[,,"new_infections"],1,sum)

plot(New_S1~dates, type="l", lwd=2, ylim=c(0, 140000), col="grey40", xaxt="n",cex.lab=1.5, 
     cex.axis=1.5, cex.main=1.5, ylab="Daily new infections", xlab="Date",
     main="Daily new infections; post-lockdown Rt = 1.2")
lines(New_S2~dates, lwd=2, col="blue4")
lines(New_S3~dates, lwd=2, col="green4")
lines(New_S4~dates, lwd=2, col="red3")
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=1.5) #x-axis dates
abline(v=c(S1_eq_date, S2_eq_date, S3_eq_date, S4_eq_date), lwd=2, lty=2, 
       col=c("grey40", "blue4", "green4", "red3"))
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)

#Exposed
E_S1 <- apply(S1.R12[,,"E"],1,sum)/total_pop
E_S2 <- apply(S2.R12[,,"E"],1,sum)/total_pop
E_S3 <- apply(S3.R12[,,"E"],1,sum)/total_pop
E_S4 <- apply(S4.R12[,,"E"],1,sum)/total_pop

par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(E_S1~dates, type="l", lwd=3, ylim=c(0,0.01), col="grey40", 
     cex.lab=2, cex.axis=2, cex.main=2, xaxt = "n",
     ylab="Proportion Exposed", xlab="Dates", main="Proportion exposed; post-lockdown Rt = 1.2")
lines(E_S2~dates, lwd=3, col="blue4")
lines(E_S3~dates, lwd=3, col="green4")
lines(E_S4~dates, lwd=3, col="red3")
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2) #x-axis dates
abline(v=c(S1_eq_date, S2_eq_date, S3_eq_date, S4_eq_date), lwd=2, lty=2, 
       col=c("grey40", "blue4", "green4", "red3"))
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)

#Susceptible when Rt=1.1
S_S1 <- apply(S1.R11[,,"S"],1,sum)/total_pop
S_S2 <- apply(S2.R11[,,"S"],1,sum)/total_pop
S_S3 <- apply(S3.R11[,,"S"],1,sum)/total_pop
S_S4 <- apply(S4.R11[,,"S"],1,sum)/total_pop

plot(S_S1~dates, type="l", lwd=2, ylim=c(0.9,1), col="grey40")
lines(S_S2~dates, lwd=2, col="blue4")
lines(S_S3~dates, lwd=2, col="green4")
lines(S_S4~dates, lwd=2, col="red3")
abline(h = 1/1.1, lty=2)

########################
## Smaller time steps ##
########################

S4.R12.dt1 <- SEIIRRS_intervention(R0=1.2, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=4, infectious_shape=3, #Natural history parameters
                                   p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                                   immune_mean_1=0, immune_mean_2=0, immune_shape=1, #scenario specific immune duration
                                   dt=1, days=days, #times
                                   BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=rep(1,15), #population data
                                   Rt_post_lockdown=1.0, intervention_post_lockdown=c(1,1,1,1), t_partial_intervention=tpi, #Partial intervention params
                                   Rt_full=1.0, intervention_full=c(1,1,1,1), t_full_intervention = tfi, #Full intervention params
                                   trigger="NA", lockdown_day=lockdown_day, Rt_partial=1.0)

S4.R12.dt05 <- SEIIRRS_intervention(R0=1.2, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=4, infectious_shape=3, #Natural history parameters
                                    p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                                    immune_mean_1=0, immune_mean_2=0, immune_shape=1, #scenario specific immune duration
                                    dt=0.5, days=days, #times
                                    BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=rep(1,15), #population data
                                    Rt_post_lockdown=1.0, intervention_post_lockdown=c(1,1,1,1), t_partial_intervention=tpi, #Partial intervention params
                                    Rt_full=1.0, intervention_full=c(1,1,1,1), t_full_intervention = tfi, #Full intervention params
                                    trigger="NA", lockdown_day=lockdown_day, Rt_partial=1.0)


S4.R12.dt01 <- SEIIRRS_intervention(R0=1.2, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=4, infectious_shape=3, #Natural history parameters
                                    p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                                    immune_mean_1=0, immune_mean_2=0, immune_shape=1, #scenario specific immune duration
                                    dt=0.1, days=days, #times
                                    BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=rep(1,15), #population data
                                    Rt_post_lockdown=1.0, intervention_post_lockdown=c(1,1,1,1), t_partial_intervention=tpi, #Partial intervention params
                                    Rt_full=1.0, intervention_full=c(1,1,1,1), t_full_intervention = tfi, #Full intervention params
                                    trigger="NA", lockdown_day=lockdown_day, Rt_partial=1.0)


par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0), mfrow=c(3,1))
plot( apply(S4.R12.dt1[,,"I"],1,sum), ylab="Infected", main="Time step = 1 day", cex.lab=1.4, cex.axis=1.4, cex.main=1.6)
plot( apply(S4.R12.dt05[,,"I"],1,sum), ylab="Infected", main="Time step = 0.5 day", cex.lab=1.4, cex.axis=1.4, cex.main=1.6)
plot( apply(S4.R12.dt01[,,"I"],1,sum), ylab="Infected", main="Time step = 0.1 day", cex.lab=1.4, cex.axis=1.4, cex.main=1.6)

S4.R12.dt1 <- SEIIRRS_intervention(R0=1.6, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=4, infectious_shape=3, #Natural history parameters
                                   p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                                   immune_mean_1=0, immune_mean_2=0, immune_shape=1, #scenario specific immune duration
                                   dt=1, days=days, #times
                                   BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=rep(5,15), #population data
                                   Rt_post_lockdown=1.0, intervention_post_lockdown=c(1,1,1,1), t_partial_intervention=tpi, #Partial intervention params
                                   Rt_full=1.0, intervention_full=c(1,1,1,1), t_full_intervention = tfi, #Full intervention params
                                   trigger="days", lockdown_day=lockdown_day, Rt_partial=1.0)

S4.R12.dt05 <- SEIIRRS_intervention(R0=1.6, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=4, infectious_shape=3, #Natural history parameters
                                    p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                                    immune_mean_1=0, immune_mean_2=0, immune_shape=1, #scenario specific immune duration
                                    dt=0.5, days=days, #times
                                    BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=rep(5,15), #population data
                                    Rt_post_lockdown=1.0, intervention_post_lockdown=c(1,1,1,1), t_partial_intervention=tpi, #Partial intervention params
                                    Rt_full=1.0, intervention_full=c(1,1,1,1), t_full_intervention = tfi, #Full intervention params
                                    trigger="days", lockdown_day=lockdown_day, Rt_partial=1.0)


S4.R12.dt01 <- SEIIRRS_intervention(R0=1.6, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=4, infectious_shape=3, #Natural history parameters
                                    p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                                    immune_mean_1=0, immune_mean_2=0, immune_shape=1, #scenario specific immune duration
                                    dt=0.1, days=days, #times
                                    BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=rep(5,15), #population data
                                    Rt_post_lockdown=1.0, intervention_post_lockdown=c(1,1,1,1), t_partial_intervention=tpi, #Partial intervention params
                                    Rt_full=1.0, intervention_full=c(1,1,1,1), t_full_intervention = tfi, #Full intervention params
                                    trigger="days", lockdown_day=lockdown_day, Rt_partial=1.0)


plot( apply(S4.R12.dt1[,,"I"],1,sum), ylab="Infected", main="Time step = 1 day", cex.lab=1.4, cex.axis=1.4, cex.main=1.6)
plot( apply(S4.R12.dt05[,,"I"],1,sum), ylab="Infected", main="Time step = 0.5 day", cex.lab=1.4, cex.axis=1.4, cex.main=1.6)
plot( apply(S4.R12.dt01[,,"I"],1,sum), ylab="Infected", main="Time step = 0.1 day", cex.lab=1.4, cex.axis=1.4, cex.main=1.6)

S4.R12.dt1 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                                   p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                                   immune_mean_1=365, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                                   dt=1, days=days, #times
                                   BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                                   Rt_post_lockdown=1.2, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                                   Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                                   trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

S4.R12.dt05 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                                    p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                                    immune_mean_1=365, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                                    days=days, dt=0.5,#times
                                    BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                                    Rt_post_lockdown=1.2, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                                    Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                                    trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)

S4.R12.dt001 <- SEIIRRS_intervention(R0=R0, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                                     p_hospitalised=p_hosp, asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                                     immune_mean_1=365, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                                     days=days, dt=0.1, #times
                                     BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                                     Rt_post_lockdown=1.2, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                                     Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                                     trigger="days", lockdown_day=lockdown_day, Rt_partial=Rt_partial)


plot( apply(S4.R12.dt1[,,"I"],1,sum), ylab="Infected", main="Time step = 1 day", cex.lab=1.4, cex.axis=1.4, cex.main=1.6)
plot( apply(S4.R12.dt05[,,"I"],1,sum), ylab="Infected", main="Time step = 0.5 day", cex.lab=1.4, cex.axis=1.4, cex.main=1.6)
plot( apply(S4.R12.dt01[,,"I"],1,sum), ylab="Infected", main="Time step = 0.1 day", cex.lab=1.4, cex.axis=1.4, cex.main=1.6)
