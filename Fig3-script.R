#Plot Figures 3 (model output) in base R
#coronavirus waning immunity

#Clear R environment
remove(list = ls())

#packages
require(scales)
require(lubridate)
library(dplyr)

#set path to current folder
dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir)

#Model scenerios
source("SEIIRRS-discrete-age-interventions-scenarios.R")

#Function to make axis numbers in scientific notation
scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

#Date variables
date_UK_lockdown <- "2020-03-23"
start_date <- as.Date(date_UK_lockdown)-lockdown_day
hour_period <- 24*dt
dates <- seq(ymd(start_date),ymd(start_date+(days-1)), by = as.difftime(hours(hour_period)))

#New infection on lockdown day
new_infections <- apply(S1.R10[,,"new_infections"],1,sum)
new_infections_adults <- apply(S1.R10[,c(5:15),"new_infections"],1,sum)
new_infections[lockdown_day]

#cumulative infections
start_cum_date <- "2020-02-16"
which(dates==start_cum_date)
sum(new_infections[which(dates==start_cum_date):which(dates==date_UK_lockdown)])
sum(new_infections_adults[which(dates==start_cum_date):which(dates==date_UK_lockdown)])

#Recovered (ages 19+)
recovered_19plus <- apply(S1.R10[,c(5:15),"R"],1,sum)
end_lockdown <- as.Date(date_UK_lockdown)+(tpi+tfi) 
which(dates==end_lockdown)
recovered_19plus[which(dates==end_lockdown)]/sum(uk.pop.2018.count$total[5:15])
sum(new_infections[1:which(dates==end_lockdown)])
recovered_all <- apply(S1.R10[,,"R"],1,sum)
recovered_all[which(dates==end_lockdown)]/total_pop

#percentage of pop infected - Imperial Report 13
imp_date <- "2020-03-28"
date_num <- which(dates==imp_date)
I_S1[date_num]/total_pop*100

##########################
## Infected by Rt value ##
##########################
#axis values
y.vals <- c(0, 5E4, 15E4, 25E4)
y.labs <- scientific(y.vals)
date_vec <- seq.Date(from=as.Date("2020-02-01"), to=as.Date("2022-02-01"), by="4 months")

#Infection vectors
I_S1 <- apply(S1.R09[,,"I"],1,sum)
I_S2 <- apply(S2.R09[,,"I"],1,sum)
I_S3 <- apply(S3.R09[,,"I"],1,sum)
I_S4 <- apply(S4.R09[,,"I"],1,sum)

#Rt = 0.9
#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-I-Rt-09.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(I_S1~dates, type="l", col="grey25",lwd=3, xaxt = "n", yaxt="n", cex.main=3,
     xlab="Date", ylab="Population Infected", cex.lab=2.5, ylim=c(0,26E4), yaxs="i",
     main=expression(paste("Post-intervention R"[t]," = 0.9")))
lines(I_S2~dates, col="darkblue", lwd=3)
lines(I_S3~dates, col="darkgreen", lwd=3)
lines(I_S4~dates, col="darkred", lwd=3)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
axis(2, y.vals, labels=y.labs, cex.axis=2.0) #y-axis scientific notation
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
legend(x=c(as.Date("2020-11-15"), as.Date("2022-02-01")), y=c(1.4E5, 2.5E5), title="Immunity Scenario",
       legend=c(expression(paste("S1  ", omega^N, " = ", infinity, ",       ", omega^H," = ", infinity)), 
                expression(paste("S2  ", omega^N, " = ", 365^-1, ", ", omega^H," = ", infinity)), 
                expression(paste("S3  ", omega^N, " = ", 180^-1, ", ", omega^H," = ", infinity)), 
                expression(paste("S4  ", omega^N, " = ", 90^-1, ",   ", omega^H," = ", 365^-1))), 
       lty=1, lwd=c(3.2,3,3,3.2), cex=2.2,
       col=c("grey25","darkblue","darkgreen", "darkred"))
#dev.off()


#Infection vectors
I_S1 <- apply(S1.R10[,,"I"],1,sum)
I_S2 <- apply(S2.R10[,,"I"],1,sum)
I_S3 <- apply(S3.R10[,,"I"],1,sum)
I_S4 <- apply(S4.R10[,,"I"],1,sum)

#Rt = 1.0
#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-I-Rt-10.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(I_S1~dates, type="l", col="grey25",lwd=3.0, xaxt = "n", yaxt="n", cex.main=3,
     xlab="Date", ylab="Population Infected", cex.lab=2.5, ylim=c(0,26E4), yaxs="i",
     main=expression(paste("Post-intervention R"[t]," = 1.0")))
lines(I_S2~dates, col="darkblue", lwd=3)
lines(I_S3~dates, col="darkgreen", lwd=3)
lines(I_S4~dates, col="darkred", lwd=3)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
axis(2, y.vals, labels=y.labs, cex.axis=2.0) #y-axis scientific notation
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
#dev.off()


#Infection vectors
I_S1 <- apply(S1.R11[,,"I"],1,sum)
I_S2 <- apply(S2.R11[,,"I"],1,sum)
I_S3 <- apply(S3.R11[,,"I"],1,sum)
I_S4 <- apply(S4.R11[,,"I"],1,sum)

#Rt = 1.1
#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-I-Rt-11.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(I_S1~dates, type="l", col="grey25",lwd=3.0, xaxt = "n", yaxt="n", cex.main=3,
     xlab="Date", ylab="Population Infected", cex.lab=2.5, ylim=c(0,26E4), yaxs="i",
     main=expression(paste("Post-intervention R"[t]," = 1.1")))
lines(I_S2~dates, col="darkblue", lwd=3)
lines(I_S3~dates, col="darkgreen", lwd=3)
lines(I_S4~dates, col="darkred", lwd=3)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
axis(2, y.vals, labels=y.labs, cex.axis=2.0) #y-axis scientific notation
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
#dev.off()

#Infection vectors
I_S1 <- apply(S1.R12[,,"I"],1,sum)
I_S2 <- apply(S2.R12[,,"I"],1,sum)
I_S3 <- apply(S3.R12[,,"I"],1,sum)
I_S4 <- apply(S4.R12[,,"I"],1,sum)

#Rt = 1.2
y.vals <- c(0, 5E4, 15E4, 25E4, 35E4, 45E4)
y.labs <- scientific(y.vals)
#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-I-Rt-12.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(I_S1~dates, type="l", col="grey25",lwd=3.0, xaxt = "n", yaxt="n", cex.main=3,
     xlab="Date", ylab="Population Infected", cex.lab=2.5, ylim=c(0,45E4), yaxs="i",
     main=expression(paste("Post-intervention R"[t]," = 1.2")))
lines(I_S2~dates, col="darkblue", lwd=3)
lines(I_S3~dates, col="darkgreen", lwd=3)
lines(I_S4~dates, col="darkred", lwd=3)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
axis(2, y.vals, labels=y.labs, cex.axis=2.0) #y-axis scientific notation
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
#dev.off()

########################
## Recovered (immune) ##
########################

#Recovered vector
R_S1 <- apply(S1.R09[,,"R"],1,sum)/total_pop
R_S2 <- apply(S2.R09[,,"R"],1,sum)/total_pop
R_S3 <- apply(S3.R09[,,"R"],1,sum)/total_pop
R_S4 <- apply(S4.R09[,,"R"],1,sum)/total_pop

#Rt = 0.9
#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-R-Rt-09.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(R_S1~dates, type="l", col="grey25",lwd=3.0, xaxt = "n", cex.main=3, cex.axis=2.0,
     xlab="Date", ylab="Proportion Immune", cex.lab=2.5, ylim=c(0,0.21), yaxs="i",
     main=expression(paste("Post-intervention R"[t]," = 0.9")))
lines(R_S2~dates, col="darkblue", lwd=3)
lines(R_S3~dates, col="darkgreen", lwd=3)
lines(R_S4~dates, col="darkred", lwd=3)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
legend(x=c(as.Date("2020-11-15"), as.Date("2022-02-01")), y=c(0.113, 0.202), title="Immunity Scenario",
       legend=c(expression(paste("S1  ", omega^N, " = ", infinity, ",       ", omega^H," = ", infinity)), 
                expression(paste("S2  ", omega^N, " = ", 365^-1, ", ", omega^H," = ", infinity)), 
                expression(paste("S3  ", omega^N, " = ", 180^-1, ", ", omega^H," = ", infinity)), 
                expression(paste("S4  ", omega^N, " = ", 90^-1, ",   ", omega^H," = ", 365^-1))), 
       lty=1, lwd=c(3.2,3,3,3.2), cex=2.2,
       col=c("grey25","darkblue","darkgreen", "darkred"))
#dev.off()

#Recovered vector
R_S1 <- apply(S1.R10[,,"R"],1,sum)/total_pop
R_S2 <- apply(S2.R10[,,"R"],1,sum)/total_pop
R_S3 <- apply(S3.R10[,,"R"],1,sum)/total_pop
R_S4 <- apply(S4.R10[,,"R"],1,sum)/total_pop

#Rt = 1.0
#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-R-Rt-10.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(R_S1~dates, type="l", col="grey25",lwd=3.0, xaxt = "n", cex.main=3, cex.axis=2.0,
     xlab="Date", ylab="Proportion Immune", cex.lab=2.5, ylim=c(0,0.21), yaxs="i",
     main=expression(paste("Post-intervention R"[t]," = 1.0")))
lines(R_S2~dates, col="darkblue", lwd=3)
lines(R_S3~dates, col="darkgreen", lwd=3)
lines(R_S4~dates, col="darkred", lwd=3)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
#dev.off()

#Recovered vector
R_S1 <- apply(S1.R11[,,"R"],1,sum)/total_pop
R_S2 <- apply(S2.R11[,,"R"],1,sum)/total_pop
R_S3 <- apply(S3.R11[,,"R"],1,sum)/total_pop
R_S4 <- apply(S4.R11[,,"R"],1,sum)/total_pop

#Rt = 1.1
#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-R-Rt-11.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(R_S1~dates, type="l", col="grey25",lwd=3.0, xaxt = "n", cex.main=3, cex.axis=2.0,
     xlab="Date", ylab="Proportion Immune", cex.lab=2.5, ylim=c(0,0.21), yaxs="i",
     main=expression(paste("Post-intervention R"[t]," = 1.1")))
lines(R_S2~dates, col="darkblue", lwd=3)
lines(R_S3~dates, col="darkgreen", lwd=3)
lines(R_S4~dates, col="darkred", lwd=3)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
#dev.off()

#Recovered vector
R_S1 <- apply(S1.R12[,,"R"],1,sum)/total_pop
R_S2 <- apply(S2.R12[,,"R"],1,sum)/total_pop
R_S3 <- apply(S3.R12[,,"R"],1,sum)/total_pop
R_S4 <- apply(S4.R12[,,"R"],1,sum)/total_pop

#Rt = 1.2
#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-R-Rt-12.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(R_S1~dates, type="l", col="grey25",lwd=3.0, xaxt = "n", cex.main=3, cex.axis=2.0,
     xlab="Date", ylab="Proportion Immune", cex.lab=2.5, ylim=c(0,0.21), yaxs="i",
     main=expression(paste("Post-intervention R"[t]," = 1.2")))
lines(R_S2~dates, col="darkblue", lwd=3)
lines(R_S3~dates, col="darkgreen", lwd=3)
lines(R_S4~dates, col="darkred", lwd=3)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
#dev.off()

col.vec <- c(rep("black",5),rep("darkred", 3), rep("black",7))

#Immunity by age group when Rt = 1.2
#scenario 1
plot(S1.R12[,1,"R"]/uk.pop.2018.count$total[1]~time, type="l", ylim=c(0,0.35), col=col.vec[1])
for(a in 2:15){
  lines(S1.R12[,a,"R"]/uk.pop.2018.count$total[a]~time, col=col.vec[a])
}

#scenario 2
plot(S2.R12[,1,"R"]/uk.pop.2018.count$total[1]~time, type="l", ylim=c(0,0.35))
for(a in 2:15){
  lines(S2.R12[,a,"R"]/uk.pop.2018.count$total[a]~time)
}

#scenario 3
plot(S3.R12[,1,"R"]/uk.pop.2018.count$total[1]~time, type="l", ylim=c(0,0.35))
for(a in 2:15){
  lines(S3.R12[,a,"R"]/uk.pop.2018.count$total[a]~time)
}

#scenario 3
plot(S3.R12[,1,"R"]/uk.pop.2018.count$total[1]~time, type="l", ylim=c(0,0.35))
for(a in 2:15){
  lines(S3.R12[,a,"R"]/uk.pop.2018.count$total[a]~time)
}

#scenario 4
plot(S4.R12[,1,"R"]/uk.pop.2018.count$total[1]~time, type="l", ylim=c(0,0.35), col=col.vec[1])
for(a in 2:15){
  lines(S4.R12[,a,"R"]/uk.pop.2018.count$total[a]~time, col=col.vec[a])
}

#equilibrium analysis
Eq.R12 <- SEIIRRS_intervention(R0=2.8, latent_mean=sigma_recip, infectious_mean=gamma_recip, latent_shape=m, infectious_shape=n, #Natural history parameters
                               p_hospitalised=rep(0,15), asymtomatic_relative_infectiousness=ari, children_relative_susceptibility=crs, phi=phi, #Natural history parameters
                               immune_mean_1=90, immune_mean_2=90, immune_shape=2, #scenario specific immune duration
                               dt=dt, days=1500, #times
                               BBC_contact_matrix=BBC_contact_matrix, total_population=total_pop, p_age=p_age, I_init=I_init, #population data
                               Rt_post_lockdown=1.2, intervention_post_lockdown=intervention_post_lockdown, t_partial_intervention=tpi, #Partial intervention params
                               Rt_full=Rt_lockdown, intervention_full=intervention_full, t_full_intervention = tfi, #Full intervention params
                               trigger="days", lockdown_day=lockdown_day, Rt_partial=1.1)

times <- seq(1,1500,dt)
eq.age.S <- apply(Eq.R12[,,"S"], 1, sum)
plot(eq.age.S~times, type="l")
EE.S <- (1/1.2)*total_pop 
abline(h=EE.S, lty=2) #S equilibrium doesn't hold

eq.age.I <- apply(Eq.R12[,,"I"], 1, sum)
plot(eq.age.I~times, type="l")
EE.I <- (((1/90)/(1.2*(1/gamma_recip)))*(1.2-1))*total_pop
abline(h=EE.I, lty=2) #I equilibrium doesn't hold
