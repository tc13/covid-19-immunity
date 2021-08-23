## Vaccine model - run scenarios

#Clear R environment
remove(list = ls())

####################################################
## Source Model Functions, Parameters & Scenarios ##
####################################################
source("Documents/covid-19-immunity/vaccination/WI-V-functions.R")
source("Documents/covid-19-immunity/vaccination/WI-V-parameters.R")

##################
## Extra params ##
##################
dt=0.25 #time step (days)
source("Documents/covid-19-immunity/vaccination/WI-V-scenarios.R")

###########
## checks ##
###########
report(V9.R015, dates = sim_dates)

################
## Beta & NGM ##
################
NGM.list <- list()
contact.mtrxs <- list()
betas <- c()

#When R0=4 at start
contact.mtrxs[[1]] = make.intervention.matrix(BBC_contact_matrix, intervention=c(1,1,1,1))
betas[1] <- get_beta(C=contact.mtrxs[[1]], R0=R0, mean_infectious = gam_mean, susceptibility = susceptibility, 
                     prop_asymtomatic = phi, asymtomatic_infectiousness = ari)
NGM.list[[1]] <- get_K(C=contact.mtrxs[[1]], R0=R0, mean_infectious = gam_mean, susceptibility = susceptibility, 
                       prop_asymtomatic = phi, asymtomatic_infectiousness = ari, beta = betas[1])

#Next generation matrixes for different R values
for(i in 2:(length(R_values)+1)){
  contact.mtrxs[[i]] = make.intervention.matrix(BBC_contact_matrix, intervention=contact_interventions[[(i-1)]])
  betas[i] <- get_beta(C=contact.mtrxs[[i]], R0=R_values[(i-1)],
                       mean_infectious = gam_mean, susceptibility = susceptibility, 
                       prop_asymtomatic = phi, asymtomatic_infectiousness = ari)
  NGM.list[[i]] <- get_K(C=contact.mtrxs[[i]], R0=R_values[(i-1)], mean_infectious = gam_mean, susceptibility = susceptibility, 
                         prop_asymtomatic = phi, asymtomatic_infectiousness = ari, beta = betas[i])
}

#############################
## R0 after June '21 = 1.5 ##
#############################
pdf("Google Drive/coronavirus/vaccination/V.R15.pdf", width=7.75, height=5.25)
par(mfrow=c(2,3), bty="n", cex.lab=1.2, cex.axis=1.1,
    mar=c(5,5,3,2))

#R0 in presence of NPIs
R.intervention <- c(rep(R0, (start_days[1])), 
                    rep(R_values, durations+1))
plot(R.intervention~sim_dates, type="l", xlab="Date", ylab="Rc (with NPIs)", lwd=1.2,
     ylim=c(0,4), main="Rc = 1.5 after 18th June '21")
abline(v=as.Date(full_reopening), lty=3, lwd=1.2)
abline(h=1, lty=3, lwd=1.2)

#Infected
plot(apply(V9.R015[,,"I"],1,sum)/total_pop~ sim_dates, type="l",
     lty=1, xlab="Date", ylab="Prop. infectious", ylim=c(0,0.005))
lines(apply(V12.R015[,,"I"],1,sum)/total_pop~ sim_dates, lty=2)
lines(apply(V18.R015[,,"I"],1,sum)/total_pop~ sim_dates, lty=3)

#Susceptible
plot(apply(V9.R015[,,"S"],1,sum)/total_pop~ sim_dates, type="l",
     lty=1, xlab="Date", ylab="Prop. susceptible", ylim=c(0,1))
lines(apply(V12.R015[,,"S"],1,sum)/total_pop~ sim_dates, lty=2)
lines(apply(V18.R015[,,"S"],1,sum)/total_pop~ sim_dates, lty=4)

#Rt - need next generation matrix
Rt.V9 = Rt.V12 = Rt.V18 = c()
for(t in 1:length(sim_dates)){
  St.V9 <- V9.R015[t,,"S"]/england_pop
  St.V12 <- V12.R015[t,,"S"]/england_pop
  St.V18 <- V18.R015[t,,"S"]/england_pop
  ngm.t <- NGM.list[[intervention.indx[t]]]
  Kij.V9 = Kij.V12 = Kij.V18 = matrix(data=0,nrow=nrow(ngm.t), ncol=ncol(ngm.t))
  for(i in 1:nrow(ngm.t)){
    for(j in 1:ncol(ngm.t)){
      Kij.V9[i,j] <- ngm.t[i,j]*St.V9[i]
      Kij.V12[i,j] <- ngm.t[i,j]*St.V12[i]
      Kij.V18[i,j] <- ngm.t[i,j]*St.V18[i]
    }
  }
  Rt.V9[t] <- get_eigen(Kij.V9)
  Rt.V12[t] <- get_eigen(Kij.V12)
  Rt.V18[t] <- get_eigen(Kij.V18)
}

plot(Rt.V9~sim_dates, type="l", xlab="Date", ylab="Rt", lwd=1,
     ylim=c(0,4))
lines(Rt.V12~sim_dates, lty=2)
lines(Rt.V18~sim_dates, lty=4)
abline(h=1, lty=3, lwd=1.2)
legend(x=c(as.Date("2020-05-01"),as.Date("2022-01-01")),
       y.intersp=0.85, x.intersp=0.8, bty="n", title.adj=0.3,
        y=c(2.2,4.2), title="Vaccine immunity",cex=1,
       lwd=c(1,1,1), lty=c(1,2,4), 
       legend=c("9 months","12 months", "18 months"))

#Vaccinated
plot(apply(V9.R015[,,"V"],1,sum)/total_pop~ sim_dates, type="l",
     lty=1, xlab="Date", ylab="Prop. vaccinated", ylim=c(0,0.8))
lines(apply(V12.R015[,,"V"],1,sum)/total_pop~ sim_dates, lty=2)
lines(apply(V18.R015[,,"V"],1,sum)/total_pop~ sim_dates, lty=4)

#Recovered
plot(apply(V9.R015[,,"R"],1,sum)/total_pop~ sim_dates, type="l",
     lty=1, xlab="Date", ylab="Prop. recovered",
     ylim=c(0,0.8))
lines(apply(V12.R015[,,"R"],1,sum)/total_pop~ sim_dates, lty=2)
lines(apply(V18.R015[,,"R"],1,sum)/total_pop~ sim_dates, lty=3)
dev.off()

############
## Output ##
############
date_indx <- which(sim_dates == full_reopening) #opening up
#Date / index Rt > 1
Rt.V9.indx = min(which(Rt.V9[date_indx:length(Rt.V9)]>1))
Rt.V12.indx = min(which(Rt.V12[date_indx:length(Rt.V12)]>1))
Rt.V18.indx = min(which(Rt.V18[date_indx:length(Rt.V18)]>1))

#Date Rt > 1
sim_dates[Rt.V9.indx+date_indx]
sim_dates[Rt.V12.indx+date_indx]
sim_dates[Rt.V18.indx+date_indx]

#Proportion susceptible
sum(V9.R015[Rt.V9.indx+date_indx,,"S"])/total_pop
sum(V12.R015[Rt.V12.indx+date_indx,,"S"])/total_pop
sum(V18.R015[Rt.V18.indx+date_indx,,"S"])/total_pop

#Resurgence
future(V9.R015, dates=sim_dates, R=1.5, case_threshold=500, current_date=full_reopening)
future(V12.R015, dates=sim_dates, R=1.5, case_threshold=500, current_date=full_reopening)
future(V18.R015, dates=sim_dates, R=1.5, case_threshold=500, current_date=full_reopening)

#############################
## R0 after June '21 = 2.0 ##
#############################
pdf("Google Drive/coronavirus/vaccination/V.R20.pdf", width=7.75, height=5.25)
par(mfrow=c(2,3), bty="n", cex.lab=1.2, cex.axis=1.1,
    mar=c(5,5,3,2))
#R0 in presence of NPIs
R.intervention <- c(rep(R0, (start_days[1])), 
                    rep(c(R_values[1:7],  2), durations+1))
plot(R.intervention~sim_dates, type="l", xlab="Date", ylab="Rc (with NPIs)", lwd=1.2,
     ylim=c(0,4), main="Rc = 2.0 after 18th June '21")
abline(v=as.Date(full_reopening), lty=3, lwd=1.2)

#Infected
plot(apply(V9.R02[,,"I"],1,sum)/total_pop~ sim_dates, type="l",
     lty=1, xlab="Date", ylab="Prop. infectious", ylim=c(0,0.005))
lines(apply(V12.R02[,,"I"],1,sum)/total_pop~ sim_dates, lty=2)
lines(apply(V18.R02[,,"I"],1,sum)/total_pop~ sim_dates, lty=4)

#Susceptible
plot(apply(V9.R02[,,"S"],1,sum)/total_pop~ sim_dates, type="l",
     lty=1, xlab="Date", ylab="Prop. susceptible", ylim=c(0,1))
lines(apply(V12.R02[,,"S"],1,sum)/total_pop~ sim_dates, lty=2)
lines(apply(V18.R02[,,"S"],1,sum)/total_pop~ sim_dates, lty=4)

#Rt - need next generation matrix
betas[9] <- get_beta(C=contact.mtrxs[[9]], R0=2,
                     mean_infectious = gam_mean, susceptibility = susceptibility, 
                     prop_asymtomatic = phi, asymtomatic_infectiousness = ari)
NGM.list[[9]] <- get_K(C=contact.mtrxs[[9]], R0=2, mean_infectious = gam_mean, susceptibility = susceptibility, 
                       prop_asymtomatic = phi, asymtomatic_infectiousness = ari, beta = betas[9])

Rt.V9 = Rt.V12 = Rt.V18 = c()
for(t in 1:length(sim_dates)){
  St.V9 <- V9.R02[t,,"S"]/england_pop
  St.V12 <- V12.R02[t,,"S"]/england_pop
  St.V18 <- V18.R02[t,,"S"]/england_pop
  ngm.t <- NGM.list[[intervention.indx[t]]]
  Kij.V9 = Kij.V12 = Kij.V18 = matrix(data=0,nrow=nrow(ngm.t), ncol=ncol(ngm.t))
  for(i in 1:nrow(ngm.t)){
    for(j in 1:ncol(ngm.t)){
      Kij.V9[i,j] <- ngm.t[i,j]*St.V9[i]
      Kij.V12[i,j] <- ngm.t[i,j]*St.V12[i]
      Kij.V18[i,j] <- ngm.t[i,j]*St.V18[i]
    }
  }
  Rt.V9[t] <- get_eigen(Kij.V9)
  Rt.V12[t] <- get_eigen(Kij.V12)
  Rt.V18[t] <- get_eigen(Kij.V18)
}

plot(Rt.V9~sim_dates, type="l", xlab="Date", ylab="Rt", lwd=1,
     ylim=c(0,4))
lines(Rt.V12~sim_dates, lty=2)
lines(Rt.V18~sim_dates, lty=4)
abline(h=1, lty=3, lwd=1.2)
legend(x=c(as.Date("2020-05-01"),as.Date("2022-01-01")),
       y.intersp=0.85, x.intersp=0.8, bty="n",title.adj=0.3,
       y=c(2.2,4.2), title="Vaccine immunity",cex=1,
       lwd=c(1,1,1), lty=c(1,2,4), 
       legend=c("9 months","12 months", "18 months"))

#Vaccinated
plot(apply(V9.R02[,,"V"],1,sum)/total_pop~ sim_dates, type="l",
     lty=1, xlab="Date", ylab="Prop. vaccinated", ylim=c(0,0.8))
lines(apply(V12.R02[,,"V"],1,sum)/total_pop~ sim_dates, lty=2)
lines(apply(V18.R02[,,"V"],1,sum)/total_pop~ sim_dates, lty=4)

#Recovered
plot(apply(V9.R02[,,"R"],1,sum)/total_pop~ sim_dates, type="l",
     lty=1, xlab="Date", ylab="Prop. recovered",
     ylim=c(0,0.8))
lines(apply(V12.R02[,,"R"],1,sum)/total_pop~ sim_dates, lty=2)
lines(apply(V18.R02[,,"R"],1,sum)/total_pop~ sim_dates, lty=4)

dev.off()

############
## Output ##
############
#Date / index Rt > 1
Rt.V9.indx = min(which(Rt.V9[date_indx:length(Rt.V9)]>1))
Rt.V12.indx = min(which(Rt.V12[date_indx:length(Rt.V12)]>1))
Rt.V18.indx = min(which(Rt.V18[date_indx:length(Rt.V18)]>1))

#Date Rt > 1
sim_dates[Rt.V9.indx+date_indx]
sim_dates[Rt.V12.indx+date_indx]
sim_dates[Rt.V18.indx+date_indx]

#Proportion susceptible
sum(V9.R02[Rt.V9.indx+date_indx,,"S"])/total_pop
sum(V12.R02[Rt.V12.indx+date_indx,,"S"])/total_pop
sum(V18.R02[Rt.V18.indx+date_indx,,"S"])/total_pop

#Resurgence
future(V9.R02, dates=sim_dates, R=2, case_threshold=500, current_date=full_reopening)
future(V12.R02, dates=sim_dates, R=2, case_threshold=500, current_date=full_reopening)
future(V18.R02, dates=sim_dates, R=2, case_threshold=500, current_date=full_reopening)

#############################
## R0 after June '21 = 2.5 ##
#############################
pdf("Google Drive/coronavirus/vaccination/V.R25.pdf", width=7.75, height=5.25)
par(mfrow=c(2,3), bty="n", cex.lab=1.2, cex.axis=1.1,
    mar=c(5,5,3,2))
#R0 in presence of NPIs
R.intervention <- c(rep(R0, (start_days[1])), 
                    rep(c(R_values[1:7],  2.5), durations+1))
plot(R.intervention~sim_dates, type="l", xlab="Date", ylab="Rc (with NPIs)", lwd=1.2,
     ylim=c(0,4), main="Rc = 2.5 after 18th June '21")
abline(v=as.Date(full_reopening), lty=3, lwd=1.2)

#Infected
plot(apply(V9.R025[,,"I"],1,sum)/total_pop~ sim_dates, type="l",
     lty=1, xlab="Date", ylab="Prop. infectious", ylim=c(0,0.005))
lines(apply(V12.R025[,,"I"],1,sum)/total_pop~ sim_dates, lty=2)
lines(apply(V18.R025[,,"I"],1,sum)/total_pop~ sim_dates, lty=4)

#Susceptible
plot(apply(V9.R025[,,"S"],1,sum)/total_pop~ sim_dates, type="l",
     lty=1, xlab="Date", ylab="Prop. susceptible", ylim=c(0,1))
lines(apply(V12.R025[,,"S"],1,sum)/total_pop~ sim_dates, lty=2)
lines(apply(V18.R025[,,"S"],1,sum)/total_pop~ sim_dates, lty=4)

#Rt - need next generation matrix
betas[9] <- get_beta(C=contact.mtrxs[[9]], R0=2.5,
                     mean_infectious = gam_mean, susceptibility = susceptibility, 
                     prop_asymtomatic = phi, asymtomatic_infectiousness = ari)
NGM.list[[9]] <- get_K(C=contact.mtrxs[[9]], R0=2.5, mean_infectious = gam_mean, susceptibility = susceptibility, 
                       prop_asymtomatic = phi, asymtomatic_infectiousness = ari, beta = betas[9])

Rt.V9 = Rt.V12 = Rt.V18 = c()
for(t in 1:length(sim_dates)){
  St.V9 <- V9.R025[t,,"S"]/england_pop
  St.V12 <- V12.R025[t,,"S"]/england_pop
  St.V18 <- V18.R025[t,,"S"]/england_pop
  ngm.t <- NGM.list[[intervention.indx[t]]]
  Kij.V9 = Kij.V12 = Kij.V18 = matrix(data=0,nrow=nrow(ngm.t), ncol=ncol(ngm.t))
  for(i in 1:nrow(ngm.t)){
    for(j in 1:ncol(ngm.t)){
      Kij.V9[i,j] <- ngm.t[i,j]*St.V9[i]
      Kij.V12[i,j] <- ngm.t[i,j]*St.V12[i]
      Kij.V18[i,j] <- ngm.t[i,j]*St.V18[i]
    }
  }
  Rt.V9[t] <- get_eigen(Kij.V9)
  Rt.V12[t] <- get_eigen(Kij.V12)
  Rt.V18[t] <- get_eigen(Kij.V18)
}

plot(Rt.V9~sim_dates, type="l", xlab="Date", ylab="Rt", lwd=1,
     ylim=c(0,4))
lines(Rt.V12~sim_dates, lty=2)
lines(Rt.V18~sim_dates, lty=4)
abline(h=1, lty=3, lwd=1.2)
legend(x=c(as.Date("2020-05-01"),as.Date("2022-01-01")),
       y.intersp=0.85, x.intersp=0.8, bty="n",title.adj=0.3,
       y=c(2.2,4.2), title="Vaccine immunity",cex=1,
       lwd=c(1,1,1), lty=c(1,2,4), 
       legend=c("9 months","12 months", "18 months"))

#Vaccinated
plot(apply(V9.R025[,,"V"],1,sum)/total_pop~ sim_dates, type="l",
     lty=1, xlab="Date", ylab="Prop. vaccinated", ylim=c(0,0.8))
lines(apply(V12.R025[,,"V"],1,sum)/total_pop~ sim_dates, lty=2)
lines(apply(V18.R025[,,"V"],1,sum)/total_pop~ sim_dates, lty=4)

#Recovered
plot(apply(V9.R025[,,"R"],1,sum)/total_pop~ sim_dates, type="l",
     lty=1, xlab="Date", ylab="Prop. recovered",
     ylim=c(0,0.8))
lines(apply(V12.R025[,,"R"],1,sum)/total_pop~ sim_dates, lty=2)
lines(apply(V18.R025[,,"R"],1,sum)/total_pop~ sim_dates, lty=4)

dev.off()

############
## Output ##
############
#Date / index Rt > 1
Rt.V9.indx = min(which(Rt.V9[date_indx:length(Rt.V9)]>1))
Rt.V12.indx = min(which(Rt.V12[date_indx:length(Rt.V12)]>1))
Rt.V18.indx = min(which(Rt.V18[date_indx:length(Rt.V18)]>1))

#Date Rt > 1
sim_dates[Rt.V9.indx+date_indx]
sim_dates[Rt.V12.indx+date_indx]
sim_dates[Rt.V18.indx+date_indx]

#Proportion susceptible
sum(V9.R025[Rt.V9.indx+date_indx,,"S"])/total_pop
sum(V12.R025[Rt.V12.indx+date_indx,,"S"])/total_pop
sum(V18.R025[Rt.V18.indx+date_indx,,"S"])/total_pop

#Resurgence
future(V9.R025, dates=sim_dates, R=2.5, case_threshold=500, current_date=full_reopening)
future(V12.R025, dates=sim_dates, R=2.5, case_threshold=500, current_date=full_reopening)
future(V18.R025, dates=sim_dates, R=2.5, case_threshold=500, current_date=full_reopening)

#############################
## R0 after June '21 = 3 ##
#############################
pdf("Google Drive/coronavirus/vaccination/V.R3.pdf", width=7.75, height=5.25)
par(mfrow=c(2,3), bty="n", cex.lab=1.2, cex.axis=1.1,
    mar=c(5,5,3,2))
#R0 in presence of NPIs
R.intervention <- c(rep(R0, (start_days[1])), 
                    rep(c(R_values[1:7],  3), durations+1))
plot(R.intervention~sim_dates, type="l", xlab="Date", ylab="Rc (with NPIs)", lwd=1.2,
     ylim=c(0,4), main="Rc = 3.0 after 18th June '21")
abline(v=as.Date(full_reopening), lty=3, lwd=1.2)

#Infected
plot(apply(V9.R03[,,"I"],1,sum)/total_pop~ sim_dates, type="l",
     lty=1, xlab="Date", ylab="Prop. infectious", ylim=c(0,0.005))
lines(apply(V12.R03[,,"I"],1,sum)/total_pop~ sim_dates, lty=2)
lines(apply(V18.R03[,,"I"],1,sum)/total_pop~ sim_dates, lty=4)

#Susceptible
plot(apply(V9.R03[,,"S"],1,sum)/total_pop~ sim_dates, type="l",
     lty=1, xlab="Date", ylab="Prop. susceptible", ylim=c(0,1))
lines(apply(V12.R03[,,"S"],1,sum)/total_pop~ sim_dates, lty=2)
lines(apply(V18.R03[,,"S"],1,sum)/total_pop~ sim_dates, lty=4)

#Rt - need next generation matrix
betas[9] <- get_beta(C=contact.mtrxs[[9]], R0=3,
                     mean_infectious = gam_mean, susceptibility = susceptibility, 
                     prop_asymtomatic = phi, asymtomatic_infectiousness = ari)
NGM.list[[9]] <- get_K(C=contact.mtrxs[[9]], R0=3, mean_infectious = gam_mean, susceptibility = susceptibility, 
                       prop_asymtomatic = phi, asymtomatic_infectiousness = ari, beta = betas[9])

Rt.V9 = Rt.V12 = Rt.V18 = c()
for(t in 1:length(sim_dates)){
  St.V9 <- V9.R03[t,,"S"]/england_pop
  St.V12 <- V12.R03[t,,"S"]/england_pop
  St.V18 <- V18.R03[t,,"S"]/england_pop
  ngm.t <- NGM.list[[intervention.indx[t]]]
  Kij.V9 = Kij.V12 = Kij.V18 = matrix(data=0,nrow=nrow(ngm.t), ncol=ncol(ngm.t))
  for(i in 1:nrow(ngm.t)){
    for(j in 1:ncol(ngm.t)){
      Kij.V9[i,j] <- ngm.t[i,j]*St.V9[i]
      Kij.V12[i,j] <- ngm.t[i,j]*St.V12[i]
      Kij.V18[i,j] <- ngm.t[i,j]*St.V18[i]
    }
  }
  Rt.V9[t] <- get_eigen(Kij.V9)
  Rt.V12[t] <- get_eigen(Kij.V12)
  Rt.V18[t] <- get_eigen(Kij.V18)
}

plot(Rt.V9~sim_dates, type="l", xlab="Date", ylab="Rt", lwd=1,
     ylim=c(0,4))
lines(Rt.V12~sim_dates, lty=2)
lines(Rt.V18~sim_dates, lty=4)
abline(h=1, lty=3, lwd=1.2)
legend(x=c(as.Date("2020-05-01"),as.Date("2022-01-01")),
       y.intersp=0.85, x.intersp=0.8, bty="n",title.adj=0.3,
       y=c(2.2,4.2), title="Vaccine immunity",cex=1,
       lwd=c(1,1,1), lty=c(1,2,4), 
       legend=c("9 months","12 months", "18 months"))

#Vaccinated
plot(apply(V9.R03[,,"V"],1,sum)/total_pop~ sim_dates, type="l",
     lty=1, xlab="Date", ylab="Prop. vaccinated", ylim=c(0,0.8))
lines(apply(V12.R03[,,"V"],1,sum)/total_pop~ sim_dates, lty=2)
lines(apply(V18.R03[,,"V"],1,sum)/total_pop~ sim_dates, lty=4)

#Recovered
plot(apply(V9.R03[,,"R"],1,sum)/total_pop~ sim_dates, type="l",
     lty=1, xlab="Date", ylab="Prop. recovered",
     ylim=c(0,0.8))
lines(apply(V12.R03[,,"R"],1,sum)/total_pop~ sim_dates, lty=2)
lines(apply(V18.R03[,,"R"],1,sum)/total_pop~ sim_dates, lty=4)

dev.off()

############
## Output ##
############
#Date / index Rt > 1
Rt.V9.indx = min(which(Rt.V9[date_indx:length(Rt.V9)]>1))
Rt.V12.indx = min(which(Rt.V12[date_indx:length(Rt.V12)]>1))
Rt.V18.indx = min(which(Rt.V18[date_indx:length(Rt.V18)]>1))

#Date Rt > 1
sim_dates[Rt.V9.indx+date_indx]
sim_dates[Rt.V12.indx+date_indx]
sim_dates[Rt.V18.indx+date_indx]

#Proportion susceptible
sum(V9.R03[Rt.V9.indx+date_indx,,"S"])/total_pop
sum(V12.R03[Rt.V12.indx+date_indx,,"S"])/total_pop
sum(V18.R03[Rt.V18.indx+date_indx,,"S"])/total_pop

#Resurgence
future(V9.R03, dates=sim_dates, R=3, case_threshold=500, current_date=full_reopening)
future(V12.R03, dates=sim_dates, R=3, case_threshold=500, current_date=full_reopening)
future(V18.R03, dates=sim_dates, R=3, case_threshold=500, current_date=full_reopening)

########################
#Vaccination
pop_size = p_age*total_pop
prop_vac <- matrix(data=0, nrow=dim(mod)[1], ncol=dim(mod)[2])
for(i in 1:length(p_age)){
  prop_vac[,i] = mod[,i,"V"]/pop_size[i]
}
#Palette
col.vec <- c("red","darkblue","darkgreen","purple","goldenrod1" ,rep("black",7))
plot(prop_vac[,15]~sim_dates, type="l", ylim=c(0,1), 
     axes=F,lwd=1.4, yaxs="i", xaxs="i",
     xlab="Date", ylab="Prop. of age group vaccinated", col="red",
     xlim=c(as.Date("2020-11-01"),as.Date("2023-03-01")))
axis(1,at=seq.Date(from=as.Date("2020-12-01"), to=as.Date("2023-03-01"),
                   by="3 months"), 
     labels = c("Dec 20", "March 21", "June 21", "Sept 21", "Dec 21",
     "March 22", "June 22", "Sept 22", "Dec 22", "March 23"))
axis(2, las=2)
for(i in 4:14){
lines(prop_vac[,i]~sim_dates, col=rev(col.vec)[(i-3)], lwd=1.4)
}
legend(x=c(as.Date("2022-05-01"),as.Date("2023-01-01")),
       y=c(0.55,1),
       title="Age group",
       legend=c("70+","65-69","60-64","55-59","50-54","under 50"), 
       col=c("red","darkblue","darkgreen","purple","goldenrod1","black"),
       cex=0.9, y.intersp=0.8, lty=1, #title.adj=0.3,
       xjust=0, lwd=1.4)

#Vaccination
prop_vac <- matrix(data=0, nrow=dim(mod)[1], ncol=dim(mod)[2])
for(i in 1:length(p_age)){
  prop_vac[,i] = mod.vs1[,i,"V"]/pop_size[i]
}
#Palette
col.vec <- c("red","darkblue","darkgreen","purple","goldenrod1" ,rep("black",7))
plot(prop_vac[,15]~sim_dates, type="l", ylim=c(0,1), 
     axes=F,lwd=1.4, yaxs="i", xaxs="i",
     xlab="Date", ylab="Prop. of age group vaccinated", col="red",
     xlim=c(as.Date("2020-11-01"),as.Date("2023-03-01")))
axis(1,at=seq.Date(from=as.Date("2020-12-01"), to=as.Date("2023-03-01"),
                   by="3 months"), 
     labels = c("Dec 20", "March 21", "June 21", "Sept 21", "Dec 21",
                "March 22", "June 22", "Sept 22", "Dec 22", "March 23"))
axis(2, las=2)
for(i in 4:14){
  lines(prop_vac[,i]~sim_dates, col=rev(col.vec)[(i-3)], lwd=1.4)
}
legend(x=c(as.Date("2022-05-01"),as.Date("2023-01-01")),
       y=c(0.55,1),
       title="Age group",
       legend=c("70+","65-69","60-64","55-59","50-54","under 50"), 
       col=c("red","darkblue","darkgreen","purple","goldenrod1","black"),
       cex=0.9, y.intersp=0.8, lty=1, title.adj=0.3,
       xjust=0, lwd=1.4)

#Vaccine shape parameter = 10
mod.vs10 <- SEIIRRS_intervention(
  #Mean durations
  latent_mean=4.5, infectious_mean=3.1, 
  immune_mean_1=150, immune_mean_2=180,
  #Shape parameters
  latent_shape=4, infectious_shape=3, immune_shape=3,
  #Vaccination
  vaccine_immune=365, vaccine_shape=10, 
  vaccine_start_days=vaccine_start_days,  #start days by age group
  vaccine_rates_daily=vaccine_rates_daily, #daily rates per age group
  vaccination_length=vaccination_length,  #duration of vaccination period in days
  #Time steps
  dt=1, days=365*3, 
  #Population
  p_age=p_age, total_population=66435550, BBC_contact_matrix=BBC_contact_matrix, 
  #Probabilities
  children_relative_susceptibility=0.4, 
  asymtomatic_relative_infectiousness=0.2,
  p_hospitalised=p_hosp, phi=phi,
  #Intervention params
  R0=4.0, start_days=start_days, durations=as.numeric(durations), R_values=R_values, contact_interventions=contact_interventions,
  #Initial values
  I_init=I_init, E_init=rep(0,15))

#Vaccination
prop_vac <- matrix(data=0, nrow=dim(mod)[1], ncol=dim(mod)[2])
for(i in 1:length(p_age)){
  prop_vac[,i] = mod.vs10[,i,"V"]/pop_size[i]
}
#Palette
col.vec <- c("red","darkblue","darkgreen","purple","goldenrod1" ,rep("black",7))
plot(prop_vac[,15]~sim_dates, type="l", ylim=c(0,1), 
     axes=F,lwd=1.4, yaxs="i", xaxs="i",
     xlab="Date", ylab="Prop. of age group vaccinated", col="red",
     xlim=c(as.Date("2020-11-01"),as.Date("2023-03-01")))
axis(1,at=seq.Date(from=as.Date("2020-12-01"), to=as.Date("2023-03-01"),
                   by="3 months"), 
     labels = c("Dec 20", "March 21", "June 21", "Sept 21", "Dec 21",
                "March 22", "June 22", "Sept 22", "Dec 22", "March 23"))
axis(2, las=2)
for(i in 4:14){
  lines(prop_vac[,i]~sim_dates, col=rev(col.vec)[(i-3)], lwd=1.4)
}
legend(x=c(as.Date("2022-05-01"),as.Date("2023-01-01")),
       y=c(0.55,1),
       title="Age group",
       legend=c("70+","65-69","60-64","55-59","50-54","under 50"), 
       col=c("red","darkblue","darkgreen","purple","goldenrod1","black"),
       cex=0.9, y.intersp=0.8, lty=1, title.adj=0.3,
       xjust=0, lwd=1.4)
