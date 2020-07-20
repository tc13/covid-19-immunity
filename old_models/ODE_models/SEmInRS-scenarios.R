### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SE^mI^nR^oS: Scenarios
### Gamma (Erlang) distributed infectious, latent and immune periods

#Clear R environment
remove(list = ls())

#Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Get model function
source("SEmInRS-model.R")

##################################################
## Scenario 1: Recovery with permanent immunity ##
##################################################

#Parameters
R0 <- 2.8                 #Basic reproduction number                    [From Petra]
sigma_recip <- 4.5        #Average duration of latent period (days)     [From He, X., et al. Temporal dynamics in viral shedding and transmissibility of COVID-19. Nat Med (2020). https://doi.org/10.1038/s41591-020-0869-5]
gamma_recip <- 3.07       #Average duration of infectiousness (days)    [From He, X., et al.]
omega <- 0                #Immunity does not wane
m <- 4                    #Shape paremeter, latent period
n <- 2                    #Shape parameter, infectious period
o <- 2                    #Shape parameter, immune period

#Vectors
time <- seq(from=1,to=100, by=1)       #time steps (days)
pars <- c(beta=R0*(1/gamma_recip), sigma=1/sigma_recip, gamma=1/gamma_recip, omega=omega, m=m, n=n, o=o) #Parameters
Y <- c(S=0.999, E1=0, E2=0, E3=0, E4=0, I1=0.0005, I2=0.0005, R1=0, R2=0) #Initial population size

#Run ODE solver function
out.gam.t.0 <- as.data.frame(
  ode(y=Y,
      func=SEmInRS,
      times=time,
      parms=pars)
)

#Sum E, I, and R compartments 
out.gam.0 <- with(out.gam.t.0, 
                  data.frame(
                    time=time, 
                    S=S, 
                    E=E1+E2+E3+E4, 
                    I=I1+I2, 
                    R=R1+R2))

#plot SE^mI^nR^oS
with(out.gam.0, plot(S~time, xlab="Days", ylab="Proportion of Population", cex.lab=1.4, cex.axis=1.4, 
                     type="l", ylim=c(0,1), col="darkblue", lwd=2))
with(out.gam.0, lines(I~time, col="darkred", lwd=2))
with(out.gam.0, lines(R~time, col="darkgreen", lwd=2))
legend(x = c(70,95), y=c(0.6,0.8), legend=c("Susceptible", "Infected", "Recovered"), lty=c(1,1,1), 
       lwd=c(2,2,2), col=c("darkblue", "darkred", "darkgreen"), cex=1.2)

#Assuming SEIR model, scenario 1 is in R environment
#Curve of I, compared between SEIR and SE^mI^nR
pdf("/Users/tomc/Google Drive/coronavirus/plots/I-curve-gam-exp.pdf", height=7, width=9)
plot(out.0$I~out.0$time, type="l", ylim=c(0,0.14), xlab="Time (days)", 
     ylab="Proportion of Population", cex.lab=1.4, cex.axis=1.4, lty=4, lwd=2, col="darkred")
lines(out.gam.0$I~out.gam.0$time, lty=1, lwd=2, col="darkred")
legend(x=c(70,100), y=c(0.1, 0.12), legend=c("Gamma", "Exponential"), lty=c(1, 4), lwd=c(2,2), cex=1.2,
       col=c("darkred", "darkred"))
dev.off()

###################################################
## Scenario 2: Immunity lasts average of 90 days ##
###################################################

#Parameters
omega_recip <- 90       #Average duration of immunity (days) [scenario assumption]

#Vectors
time <- seq(from=1,to=300, by=1)       #time steps (days)
pars <- c(beta=R0*(1/gamma_recip), sigma=1/sigma_recip, gamma=1/gamma_recip, omega=1/omega_recip, m=m, n=n, o=o) #Parameters
Y <- c(S=0.999, E1=0, E2=0, E3=0, E4=0, I1=0.0005, I2=0.0005, R1=0, R2=0) #Initial population size

#Run ODE solver function
out.gam.t.90 <- as.data.frame(
  ode(y=Y,
      func=SEmInRS,
      times=time,
      parms=pars)
)

#Sum E, I, and R compartments 
out.gam.90 <- with(out.gam.t.90, 
                   data.frame(
                     time=time, 
                     S=S, 
                     E=E1+E2+E3+E4, 
                     I=I1+I2, 
                     R=R1+R2))

#Equilibrium values
S.eq <- 1/R0
I.eq <- (pars["omega"]/pars["beta"])*(R0-1)
E.eq <- (pars["omega"]*(pars["omega"]+pars["gamma"])/pars["beta"]*pars["sigma"])*(R0-1)
R.eq <- 1-S.eq-I.eq-E.eq

#plot SE^mI^nR^oS
#pdf("/Users/tomc/Google Drive/coronavirus/plots/SEmInRoS-scenario2.pdf", height=7, width=9)
with(out.gam.90, plot(S~time, xlab="Days", ylab="Proportion of Population", cex.lab=1.4, cex.axis=1.4, 
                      type="l", ylim=c(0,1), col="darkblue", lwd=2))
with(out.gam.90, lines(I~time, col="darkred", lwd=2))
with(out.gam.90, lines(R~time, col="darkgreen", lwd=2))
legend(x = c(220,300), y=c(0.75,0.95), legend=c("Susceptible", "Infected", "Recovered"), lty=c(1,1,1), 
       lwd=c(2,2,2), col=c("darkblue", "darkred", "darkgreen"), cex=1.2)
abline(h=S.eq, lty="dashed", lwd=2, col="darkblue") #equilbrium proportion susceptible
abline(h=I.eq, lty="dashed", lwd=2, col="darkred") #equilbrium proportion infected
abline(h=R.eq, lty="dashed", lwd=2, col="darkgreen") #equilibrium proportion recovered
#dev.off()

#Assuming SEIR model, scenario 2 is in R environment
#pdf("/Users/tomc/Google Drive/coronavirus/plots/S-curve-gam-exp.pdf", height=7, width=9)
plot(out.90$S~out.90$time, type="l", xlab="Days", ylim=c(0,1), col="darkblue",
     ylab="Proportion of Population", cex.lab=1.4, cex.axis=1.4, lty=4, lwd=2)
lines(out.gam.90$S~out.gam.90$time, lty=1, lwd=2, col="darkblue")
legend(x = c(220,300), y=c(0.75,0.9), legend=c("Gamma", "Exponential"), lty=c(1, 4), 
       lwd=c(2,2), cex=1.2, col=c("darkblue", "darkblue"))
#dev.off()

####################################################
## Scenario 3: Immunity lasts average of 180 days ##
####################################################

#Parameters
omega_recip <- 180       #Average duration of immunity (days) [scenario assumption]

#Vectors
time <- seq(from=1,to=500, by=1)       #time steps (days)
pars <- c(beta=R0*(1/gamma_recip), sigma=1/sigma_recip, gamma=1/gamma_recip, omega=1/omega_recip, m=m, n=n, o=o) #Parameters
Y <- c(S=0.999, E1=0, E2=0, E3=0, E4=0, I1=0.0005, I2=0.0005, R1=0, R2=0) #Initial population size

#Run ODE solver function
out.gam.t.180 <- as.data.frame(
  ode(y=Y,
      func=SEmInRS,
      times=time,
      parms=pars)
)

#Sum E, I, and R compartments 
out.gam.180 <- with(out.gam.t.180, 
                    data.frame(
                      time=time, 
                      S=S, 
                      E=E1+E2+E3+E4, 
                      I=I1+I2, 
                      R=R1+R2))

#Equilibrium values
S.eq <- 1/R0
I.eq <- (pars["omega"]/pars["beta"])*(R0-1)
E.eq <- (pars["omega"]*(pars["omega"]+pars["gamma"])/pars["beta"]*pars["sigma"])*(R0-1)
R.eq <- 1-S.eq-I.eq-E.eq

#plot SE^mI^nR^oS
#pdf("/Users/tomc/Google Drive/coronavirus/plots/SEmInRoS-scenario3.pdf", height=7, width=9)
with(out.gam.180, plot(S~time, xlab="Days", ylab="Proportion of Population", cex.lab=1.4, cex.axis=1.4, 
                       type="l", ylim=c(0,1), col="darkblue", lwd=2))
with(out.gam.180, lines(I~time, col="darkred", lwd=2))
with(out.gam.180, lines(R~time, col="darkgreen", lwd=2))
legend(x = c(380,510), y=c(0.75,0.95), legend=c("Susceptible", "Infected", "Recovered"), lty=c(1,1,1), 
       lwd=c(2,2,2), col=c("darkblue", "darkred", "darkgreen"), cex=1.2)
abline(h=S.eq, lty="dashed", lwd=2, col="darkblue") #equilbrium proportion susceptible
abline(h=I.eq, lty="dashed", lwd=2, col="darkred") #equilbrium proportion infected
abline(h=R.eq, lty="dashed", lwd=2, col="darkgreen") #equilibrium proportion recovered
#dev.off()

#Assuming SEIR model, scenario 3 is in R environment
#Secondary peak in I compared between SEIR and SE^mI^nR^oS
#pdf("/Users/tomc/Google Drive/coronavirus/plots/I-secondary-curve.pdf", height=7, width=9)
plot(out.gam.180$I~out.gam.180$time, type="l", xlab="Days", ylim=c(0,0.08),
     ylab="Proportion of Population", cex.lab=1.4, cex.axis=1.4, lty=1, lwd=2, 
     xlim=c(200,350), col="darkred")
lines(out.180$I~out.180$time, lty=4, lwd=2, col="darkred")
legend(x = c(310,350), y=c(0.06,0.0725), legend=c("Gamma", "Exponential"), lty=c(1, 4), lwd=c(2,2), 
       cex=1.2, col=c("darkred", "darkred"))
#dev.off()
