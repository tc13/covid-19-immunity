### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SE^mI^nR^oS: scenarios
### Epidemic dynamics in the context of (non-pharmaceutical) interventions

#Clear R environment
remove(list = ls())

#Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Get SEmInRoS model function
source("SEmInRS-model.R")
#Get intervention model function
source("NPI-model.R")

##########################
## Key model parameters ##
##########################

sigma_recip <- 4.5        #Average duration of latent period (days)     [From He, X., et al. Temporal dynamics in viral shedding and transmissibility of COVID-19. Nat Med (2020). https://doi.org/10.1038/s41591-020-0869-5]
gamma_recip <- 3.07       #Average duration of infectiousness (days)    [From He, X., et al.]
omega <- 0                #Immunity does not wane
m <- 4                    #Shape paremeter, latent period
n <- 2                    #Shape parameter, infectious period
o <- 2                    #Shape parameter, immune period

##############################
## User variable parameters ##
##############################

R0 <- 2.8               #Basic reproduction number (in the absense of interventions) [From Petra]
Rc <- 0.9               #Effective repoduction number during interventions 
C_i <- 0.015            #Critical proportion infected that triggers intervention
tau <- 60               #Duration of intervention (days)

#############
## Vectors ##
#############

time <- seq(from=1,to=300, by=1)       #time steps (days)
yini <- c(S=0.999, E1=0, E2=0, E3=0, E4=0, I1=0.0005, I2=0.0005, R1=0, R2=0) #Initial population size

##############################################
## Get intervention timepoint from SEmInRoS ##
##############################################

pars <- c(beta=R0*(1/gamma_recip), sigma=1/sigma_recip, gamma=1/gamma_recip, omega=1/90, m=m, n=n, o=o) #Parameters
baseline_temp <- as.data.frame(ode(func = SEmInRS, y=yini, parms=pars, times=time))
baseline <-  with(baseline_temp, data.frame(time=time, S=S, E=E1+E2+E3+E4, I=I1+I2, R=R1+R2)) #Sum E, I, and R compartments 
time_NPI <- min(which(baseline$I>=C_i))
rm(baseline_temp)

#################################
## Time dependent transmission ##
#################################

transmission <- data.frame(time=time, beta=rep(0, length(time)))
transmission$beta <- ifelse(transmission$time>=time_NPI & transmission$time <= (time_NPI+tau), Rc*(1/gamma_recip), R0*(1/gamma_recip))
beta_fn <- approxfun(transmission, rule=2) #interpolate function

########################
## Intervention model ##
########################

pars <- c(sigma=1/sigma_recip, gamma=1/gamma_recip, omega=1/90, m=m, n=n, o=o) #Parameters
temp <- as.data.frame(ode(y=yini,func=NPI,times=time,parms=pars))
out.NPI <- with(temp, data.frame(time=time, S=S, E=E1+E2+E3+E4, I=I1+I2, R=R1+R2)) #Sum E, I, and R compartments 
rm(temp)

########################
## Equilibrium values ##
########################

S.eq <- 1/R0
I.eq <- (pars["omega"]/(R0*pars["gamma"]))*(R0-1)
E.eq <- (pars["omega"]*(pars["omega"]+pars["gamma"])/(R0*pars["gamma"])*pars["sigma"])*(R0-1)
R.eq <- 1-S.eq-I.eq-E.eq

##############
## Plotting ##
##############

#pdf("/Users/tomc/Google Drive/coronavirus/plots/SEmInRoS-NPI1.pdf", height=7, width=9)
with(out.NPI, plot(S~time, xlab="Days", ylab="Proportion of Population", cex.lab=1.4, cex.axis=1.4, 
                      type="l", ylim=c(0,1), col="darkblue", lwd=2))
with(out.NPI, lines(I~time, col="darkred", lwd=2))
with(out.NPI, lines(R~time, col="darkgreen", lwd=2))
legend(x = c(220,300), y=c(0.75,0.95), legend=c("Susceptible", "Infected", "Recovered"), lty=c(1,1,1), 
       lwd=c(2,2,2), col=c("darkblue", "darkred", "darkgreen"), cex=1.2)
#abline(h=S.eq, lty="dashed", lwd=2, col="darkblue") #equilbrium proportion susceptible
#abline(h=I.eq, lty="dashed", lwd=2, col="darkred") #equilbrium proportion infected
#abline(h=R.eq, lty="dashed", lwd=2, col="darkgreen") #equilibrium proportion recovered
abline(v=c(time_NPI, time_NPI+tau), lty=c(2,2), lwd=c(2,2)) #intervention period
#dev.off()

#pdf("/Users/tomc/Google Drive/coronavirus/plots/I-interventions.pdf", height=7, width=9)
with(out.NPI, plot(I~time, xlab="Days", ylab="Proportion of Population", cex.lab=1.4, cex.axis=1.4, 
                     type="l", col="darkred", lwd=2, ylim=c(0,0.14), lty=4))
abline(v=c(time_NPI, time_NPI+tau), lty=c(2,2), lwd=c(2,2)) #intervention period
with(baseline, lines(I~time, col="darkred", lwd=2, lty=1))
legend(x = c(200,300), y=c(0.12,0.14), legend=c("Infected (baseline)", "Infected (with NPI)"), lty=c(1,4), 
       lwd=c(2,2), col=c("darkred", "darkred"), cex=1.2)
abline(h=I.eq, lty="dashed", lwd=2, col="darkred") #equilbrium proportion infected
#dev.off()
