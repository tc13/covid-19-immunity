### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SEIRS: Scenarios
### exponentially distributed infectious period

#Clear R environment
#remove(list = ls())

#Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Get SEIRS model function
source("SEIRS-model.R")

##################################################
## Scenario 1: Recovery with permanent immunity ##
##################################################

#Parameters
R0 <- 2.8               #Basic reproduction number                    [From Petra]
sigma_recip <- 4.5      #Average duration of latent period (days)     [From He, X., et al. Temporal dynamics in viral shedding and transmissibility of COVID-19. Nat Med (2020). https://doi.org/10.1038/s41591-020-0869-5]
gamma_recip <- 3.1      #Average duration of infectiousness (days)    [From He, X., et al.]
omega <- 0              #Average duration of immunity (days) [scenario assumption]

#Vectors
time <- seq(from=1,to=100, by=1)       #time steps (days)
pars <- c(beta=R0*(1/gamma_recip), sigma=1/sigma_recip, gamma=1/gamma_recip, omega=omega) #parameters
Y <- c(S=0.999, E=0, I=0.001, R=0)      #Initial population size

#Run ODE solver function
out.0 <- as.data.frame(
  ode(y=Y,
      func=SEIRS,
      times=time,
      parms=pars))

#plot SEIRS
with(out.0, plot(S~time, xlab="Days", ylab="Proportion of Population", cex.lab=1.4, cex.axis=1.4, 
                 type="l", ylim=c(0,1), col="darkblue", lwd=2))
with(out.0, lines(I~time, col="darkred", lwd=2))
with(out.0, lines(R~time, col="darkgreen", lwd=2))
legend(x = c(70,95), y=c(0.6,0.8), legend=c("Susceptible", "Infected", "Recovered"), lty=c(1,1,1), 
       lwd=c(2,2,2), col=c("darkblue", "darkred", "darkgreen"), cex=1.2)

###################################################
## Scenario 2: Immunity lasts average of 90 days ##
###################################################

#Parameters
omega_recip <- 90       #Average duration of immunity (days) [scenario assumption]

#Vectors
time <- seq(from=1,to=300, by=1)       #time steps (days)
pars <- c(beta=R0*(1/gamma_recip), sigma=1/sigma_recip, gamma=1/gamma_recip, omega=1/omega_recip) #parameters
Y <- c(S=0.999, E=0, I=0.001, R=0)      #Initial population size

#Run ODE solver function
out.90 <- as.data.frame(
  ode(y=Y,
      func=SEIRS,
      times=time,
      parms=pars))

#plot SEIRS
with(out.90, plot(S~time, xlab="Days", ylab="Proportion of Population", cex.lab=1.4, cex.axis=1.4, 
                  type="l", ylim=c(0,1), col="darkblue", lwd=2))
with(out.90, lines(I~time, col="darkred", lwd=2))
with(out.90, lines(R~time, col="darkgreen", lwd=2))
legend(x = c(220,300), y=c(0.75,0.95), legend=c("Susceptible", "Infected", "Recovered"), lty=c(1,1,1), 
       lwd=c(2,2,2), col=c("darkblue", "darkred", "darkgreen"), cex=1.2)
abline(h=1/R0, lty="dashed") #equilbrium proportion of susceptibles

####################################################
## Scenario 3: Immunity lasts average of 180 days ##
####################################################

#Parameters
omega_recip <- 180       #Average duration of immunity (days) [scenario assumption]

#Vectors
time <- seq(from=1,to=500, by=1)       #time steps (days)
pars <- c(beta=R0*(1/gamma_recip), sigma=1/sigma_recip, gamma=1/gamma_recip, omega=1/omega_recip) #parameters
Y <- c(S=0.999, E=0, I=0.001, R=0)      #Initial population size

#Run ODE solver function
out.180 <- as.data.frame(
  ode(y=Y,
      func=SEIRS,
      times=time,
      parms=pars))

#plot SEIRS
with(out.180, plot(S~time, xlab="Days", ylab="Proportion of Population", cex.lab=1.4, cex.axis=1.4, 
                   type="l", ylim=c(0,1), col="darkblue", lwd=2))
with(out.180, lines(I~time, col="darkred", lwd=2))
with(out.180, lines(R~time, col="darkgreen", lwd=2))
legend(x = c(380,510), y=c(0.75,0.95), legend=c("Susceptible", "Infected", "Recovered"), lty=c(1,1,1), 
       lwd=c(2,2,2), col=c("darkblue", "darkred", "darkgreen"), cex=1.2)
abline(h=1/R0, lty="dashed") #equilbrium proportion of susceptibles

#############################################
## Scenario 4: No immunity after infection ##
#############################################

#Parameters
omega_recip <- 1        #Average duration of immunity (days) [scenario assumption]

#Vectors
time <- seq(from=1,to=100, by=1)       #time steps (days)
pars <- c(beta=R0*(1/gamma_recip), sigma=1/sigma_recip, gamma=1/gamma_recip, omega=1/omega_recip) #parameters
Y <- c(S=0.999, E=0, I=0.001, R=0)      #Initial population size

#Run ODE solver function
out.1 <- as.data.frame(
  ode(y=Y,
      func=SEIRS,
      times=time,
      parms=pars))

#plot
#pdf("Google Drive/coronavirus/SEIRS-scenario4.pdf", height=7, width=9)
with(out.1, plot(S~time, xlab="Days", ylab="Proportion of Population", cex.lab=1.4, cex.axis=1.4, 
                 type="l", ylim=c(0,1), col="darkblue", lwd=2))
with(out.1, lines(I~time, col="darkred", lwd=2))
with(out.1, lines(R~time, col="darkgreen", lwd=2))
legend(x = c(70,100), y=c(0.75,0.95), legend=c("Susceptible", "Infected", "Recovered"), lty=c(1,1,1), 
       lwd=c(2,2,2), col=c("darkblue", "darkred", "darkgreen"), cex=1.2)
abline(h=1/R0, lty="dashed") #equilbrium proportion of susceptibles
#dev.off()
