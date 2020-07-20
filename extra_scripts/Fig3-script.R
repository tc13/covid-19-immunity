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
y.vals <- c(0, 5E4, 15E4, 25E4, 35E4, 45E4)
y.labs <- scientific(y.vals)
date_vec <- seq.Date(from=as.Date("2020-02-01"), to=as.Date("2022-02-01"), by="4 months")

#Infection vectors
I_S1 <- apply(S1.R09[,,"I"],1,sum)
I_S2 <- apply(S2.R09[,,"I"],1,sum)
I_S3 <- apply(S3.R09[,,"I"],1,sum)
I_S4 <- apply(S4.R09[,,"I"],1,sum)

#Rt = 0.9
pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-I-Rt-09.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(I_S1~dates, type="l", col="grey25",lwd=3.5, xaxt = "n", yaxt="n", cex.main=3,
     xlab="Date", ylab="Population Infected", cex.lab=2.5, ylim=c(0,45E4), yaxs="i",
     main=expression(paste("Post-lockdown R"[t]," = 0.9")))
lines(I_S2~dates, col="darkblue", lwd=3.5)
lines(I_S3~dates, col="darkgreen", lwd=3.5)
lines(I_S4~dates, col="darkred", lwd=3.5)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
axis(2, y.vals, labels=y.labs, cex.axis=2.0) #y-axis scientific notation
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
#legend(x=c(as.Date("2020-11-15"), as.Date("2022-02-01")), y=c(242308, 432692), title="Immunity Scenario",
#       legend=c(expression(paste("S1  ", omega^N, " = ", infinity, ",       ", omega^H," = ", infinity)), 
#                expression(paste("S2  ", omega^N, " = ", 365^-1, ", ", omega^H," = ", infinity)), 
#                expression(paste("S3  ", omega^N, " = ", 180^-1, ", ", omega^H," = ", infinity)), 
#                expression(paste("S4  ", omega^N, " = ", 90^-1, ",   ", omega^H," = ", 365^-1))), 
#       lty=1, lwd=c(3.2,3,3,3.2), cex=2.2,
#       col=c("grey25","darkblue","darkgreen", "darkred"))
legend(x=c(as.Date("2020-11-15"), as.Date("2022-02-01")), y=c(242308, 432692), title="Immunity Scenario",
       legend=c("S1: Complete", "S2: Waning (12 months)", "S3: Waning (6 months)", "S4: Short-lived"), 
       lty=1, lwd=c(3.7,3.5,3.5,3.7), cex=2.2,
       col=c("grey25","darkblue","darkgreen", "darkred"))
dev.off()

#Infection vectors
I_S1 <- apply(S1.R10[,,"I"],1,sum)
I_S2 <- apply(S2.R10[,,"I"],1,sum)
I_S3 <- apply(S3.R10[,,"I"],1,sum)
I_S4 <- apply(S4.R10[,,"I"],1,sum)

#Rt = 1.0
pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-I-Rt-10.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(I_S1~dates, type="l", col="grey25",lwd=3.5, xaxt = "n", yaxt="n", cex.main=3,
     xlab="Date", ylab="Population Infected", cex.lab=2.5, ylim=c(0,45E4), yaxs="i",
     main=expression(paste("Post-lockdown R"[t]," = 1.0")))
lines(I_S2~dates, col="darkblue", lwd=3.5)
lines(I_S3~dates, col="darkgreen", lwd=3.5)
lines(I_S4~dates, col="darkred", lwd=3.5)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
axis(2, y.vals, labels=y.labs, cex.axis=2.0) #y-axis scientific notation
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
dev.off()

#Infection vectors
I_S1 <- apply(S1.R11[,,"I"],1,sum)
I_S2 <- apply(S2.R11[,,"I"],1,sum)
I_S3 <- apply(S3.R11[,,"I"],1,sum)
I_S4 <- apply(S4.R11[,,"I"],1,sum)

#Rt = 1.1
pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-I-Rt-11.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(I_S1~dates, type="l", col="grey25",lwd=3.5, xaxt = "n", yaxt="n", cex.main=3,
     xlab="Date", ylab="Population Infected", cex.lab=2.5, ylim=c(0,45E4), yaxs="i",
     main=expression(paste("Post-lockdown R"[t]," = 1.1")))
lines(I_S2~dates, col="darkblue", lwd=3.5)
lines(I_S3~dates, col="darkgreen", lwd=3.5)
lines(I_S4~dates, col="darkred", lwd=3.5)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
axis(2, y.vals, labels=y.labs, cex.axis=2.0) #y-axis scientific notation
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
dev.off()

#Infection vectors
I_S1 <- apply(S1.R12[,,"I"],1,sum)
I_S2 <- apply(S2.R12[,,"I"],1,sum)
I_S3 <- apply(S3.R12[,,"I"],1,sum)
I_S4 <- apply(S4.R12[,,"I"],1,sum)

#Rt = 1.2
pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-I-Rt-12.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(I_S1~dates, type="l", col="grey25",lwd=3.5, xaxt = "n", yaxt="n", cex.main=3,
     xlab="Date", ylab="Population Infected", cex.lab=2.5, ylim=c(0,45E4), yaxs="i",
     main=expression(paste("Post-lockdown R"[t]," = 1.2")))
lines(I_S2~dates, col="darkblue", lwd=3.5)
lines(I_S3~dates, col="darkgreen", lwd=3.5)
lines(I_S4~dates, col="darkred", lwd=3.5)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
axis(2, y.vals, labels=y.labs, cex.axis=2.0) #y-axis scientific notation
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
dev.off()

########################
## Recovered (immune) ##
########################

#Recovered vector
R_S1 <- apply(S1.R09[,,"R"],1,sum)/total_pop
R_S2 <- apply(S2.R09[,,"R"],1,sum)/total_pop
R_S3 <- apply(S3.R09[,,"R"],1,sum)/total_pop
R_S4 <- apply(S4.R09[,,"R"],1,sum)/total_pop

#Rt = 0.9
pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-R-Rt-09.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(R_S1~dates, type="l", col="grey25",lwd=3.5, xaxt = "n", cex.main=3, cex.axis=2.0,
     xlab="Date", ylab="Proportion Immune", cex.lab=2.5, ylim=c(0,0.21), yaxs="i",
     main=expression(paste("Post-lockdown R"[t]," = 0.9")))
lines(R_S2~dates, col="darkblue", lwd=3.5)
lines(R_S3~dates, col="darkgreen", lwd=3.5)
lines(R_S4~dates, col="darkred", lwd=3.5)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
#legend(x=c(as.Date("2020-11-15"), as.Date("2022-02-01")), y=c(0.113, 0.202), title="Immunity Scenario",
#       legend=c(expression(paste("S1  ", omega^N, " = ", infinity, ",       ", omega^H," = ", infinity)), 
#                expression(paste("S2  ", omega^N, " = ", 365^-1, ", ", omega^H," = ", infinity)), 
#                expression(paste("S3  ", omega^N, " = ", 180^-1, ", ", omega^H," = ", infinity)), 
#                expression(paste("S4  ", omega^N, " = ", 90^-1, ",   ", omega^H," = ", 365^-1))), 
#       lty=1, lwd=c(3.2,3,3,3.2), cex=2.2,
#       col=c("grey25","darkblue","darkgreen", "darkred"))
legend(x=c(as.Date("2020-11-15"), as.Date("2022-02-01")), y=c(0.113, 0.202), title="Immunity Scenario",
       legend=c("S1: Complete", "S2: Waning (12 months)", "S3: Waning (6 months)", "S4: Short-lived"), 
       lty=1, lwd=c(3.7,3.5,3.5,3.7), cex=2.2,
       col=c("grey25","darkblue","darkgreen", "darkred"))
dev.off()

#Recovered vector
R_S1 <- apply(S1.R10[,,"R"],1,sum)/total_pop
R_S2 <- apply(S2.R10[,,"R"],1,sum)/total_pop
R_S3 <- apply(S3.R10[,,"R"],1,sum)/total_pop
R_S4 <- apply(S4.R10[,,"R"],1,sum)/total_pop

#Rt = 1.0
pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-R-Rt-10.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(R_S1~dates, type="l", col="grey25",lwd=3.5, xaxt = "n", cex.main=3, cex.axis=2.0,
     xlab="Date", ylab="Proportion Immune", cex.lab=2.5, ylim=c(0,0.21), yaxs="i",
     main=expression(paste("Post-lockdown R"[t]," = 1.0")))
lines(R_S2~dates, col="darkblue", lwd=3.5)
lines(R_S3~dates, col="darkgreen", lwd=3.5)
lines(R_S4~dates, col="darkred", lwd=3.5)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
dev.off()

#Recovered vector
R_S1 <- apply(S1.R11[,,"R"],1,sum)/total_pop
R_S2 <- apply(S2.R11[,,"R"],1,sum)/total_pop
R_S3 <- apply(S3.R11[,,"R"],1,sum)/total_pop
R_S4 <- apply(S4.R11[,,"R"],1,sum)/total_pop

#Rt = 1.1
pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-R-Rt-11.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(R_S1~dates, type="l", col="grey25",lwd=3.5, xaxt = "n", cex.main=3, cex.axis=2.0,
     xlab="Date", ylab="Proportion Immune", cex.lab=2.5, ylim=c(0,0.21), yaxs="i",
     main=expression(paste("Post-lockdown R"[t]," = 1.1")))
lines(R_S2~dates, col="darkblue", lwd=3.5)
lines(R_S3~dates, col="darkgreen", lwd=3.5)
lines(R_S4~dates, col="darkred", lwd=3.5)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
dev.off()

#Recovered vector
R_S1 <- apply(S1.R12[,,"R"],1,sum)/total_pop
R_S2 <- apply(S2.R12[,,"R"],1,sum)/total_pop
R_S3 <- apply(S3.R12[,,"R"],1,sum)/total_pop
R_S4 <- apply(S4.R12[,,"R"],1,sum)/total_pop

#Rt = 1.2
pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-R-Rt-12.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(R_S1~dates, type="l", col="grey25",lwd=3.5, xaxt = "n", cex.main=3, cex.axis=2.0,
     xlab="Date", ylab="Proportion Immune", cex.lab=2.5, ylim=c(0,0.21), yaxs="i",
     main=expression(paste("Post-lockdown R"[t]," = 1.2")))
lines(R_S2~dates, col="darkblue", lwd=3.5)
lines(R_S3~dates, col="darkgreen", lwd=3.5)
lines(R_S4~dates, col="darkred", lwd=3.5)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
dev.off()

#########################################
## Immunity by age group when Rt = 1.2 ##
#########################################
col.vec <- c(rep("grey50",5),"darkorange","deepskyblue3","forestgreen", rep("grey50",6), "darkviolet")
lwd.vec <- c(rep(1.5,5), 3, 3, 3, rep(1.5,6), 2.5)
#Scenario 1, recovered
pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-S1-R12-R-age-groups.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0))
plot(S1.R12[,1,"R"]/uk.pop.2018.count$total[1]~dates, type="l", ylim=c(0,0.35), 
     col=col.vec[1], lwd=lwd.vec[1], xaxt = "n", cex.main=3, cex.axis=2.0,
     xlab="Date", ylab="Proportion of age group immune", cex.lab=2.5, yaxs="i",
     main=expression(paste("Scenario 1, Post-lockdown R"[t]," = 1.2")))
for(a in 2:15){
  lines(S1.R12[,a,"R"]/uk.pop.2018.count$total[a]~dates, col=col.vec[a], lwd=lwd.vec[a])
}
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
legend(x=c(as.Date("2020-07-01"), as.Date("2021-01-01")),
       y=c(0.225, 0.34), title = "Age groups", cex=1.5,
  legend = c("25-29 years", "30-34 years", "35-39 years", "70+ years", "Others"),
       lwd=c(3.2,3.2,3.2,3,2), col=c("darkorange","deepskyblue3","forestgreen", "darkviolet", "grey50"))
dev.off()

#Infected
pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-S1-R12-I-age-groups.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0))
plot(S1.R12[,1,"I"]/uk.pop.2018.count$total[1]~dates, type="l", ylim=c(0, 0.01), 
     col=col.vec[1], lwd=lwd.vec[1], xaxt = "n", cex.main=3, cex.axis=2.0, yaxs="i",
     xlab="Date", ylab="Proportion of age group infected", cex.lab=2.5,
     main=expression(paste("Scenario 1, Post-lockdown R"[t]," = 1.2")))
for(a in 2:15){
  lines(S1.R12[,a,"I"]/uk.pop.2018.count$total[a]~dates, col=col.vec[a], lwd=lwd.vec[a])
}
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
legend(x=c(as.Date("2020-07-01"), as.Date("2021-01-01")),
       y=c(0.006428571, 0.009714286), title = "Age groups", cex=1.5,
       legend = c("25-29 years", "30-34 years", "35-39 years", "70+ years", "Others"),
       lwd=c(3.2,3.2,3.2,3,2), col=c("darkorange","deepskyblue3","forestgreen", "darkviolet", "grey50"))
dev.off()

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

#scenario 4
pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-S4-R12-R-age-groups.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0))
plot(S4.R12[,1,"R"]/uk.pop.2018.count$total[1]~dates, type="l", ylim=c(0,0.35), 
     col=col.vec[1], lwd=lwd.vec[1], xaxt = "n", cex.main=3, cex.axis=2.0,
     xlab="Date", ylab="Proportion of age group immune", cex.lab=2.5, yaxs="i",
     main=expression(paste("Scenario 4, Post-lockdown R"[t]," = 1.2")))
for(a in 2:15){
  lines(S4.R12[,a,"R"]/uk.pop.2018.count$total[a]~dates, col=col.vec[a], lwd=lwd.vec[a])
}
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
dev.off()

#infected
pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-S4-R12-I-age-groups.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0))
plot(S4.R12[,1,"I"]/uk.pop.2018.count$total[1]~dates, type="l", ylim=c(0, 0.01), 
     col=col.vec[1], lwd=lwd.vec[1], xaxt = "n", cex.main=3, cex.axis=2.0, yaxs="i",
     xlab="Date", ylab="Proportion of age group infected", cex.lab=2.5,
     main=expression(paste("Scenario 4, Post-lockdown R"[t]," = 1.2")))
for(a in 2:15){
  lines(S4.R12[,a,"I"]/uk.pop.2018.count$total[a]~dates, col=col.vec[a], lwd=lwd.vec[a])
}
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
dev.off()

#height of secondary peak
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

pdf("/Users/tomc/Google Drive/coronavirus/Figures/Rt12-S-turnover.pdf", height=7, width=10)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0))
plot(S_S1~dates, type="l", lwd=3, ylim=c(0.8,1), col="grey25", xaxt="n",cex.lab=2, 
     cex.axis=2, cex.main=2, ylab="Proportion susceptible", xlab="Date",
     main="Proportion susceptible; post-lockdown Rt = 1.2")
lines(S_S2~dates, lwd=3, col="darkblue")
lines(S_S3~dates, lwd=3, col="darkgreen")
lines(S_S4~dates, lwd=3, col="darkred")
abline(h = 1/1.2, lty=2, lwd=3)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2) #x-axis dates
abline(v=c(S1_eq_date, S2_eq_date, S3_eq_date, S4_eq_date), lwd=2, lty=2, 
       col=c("grey25", "darkblue", "darkgreen", "darkred"))
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2)
dev.off()

#Infection vectors
I_S1.p <- apply(S1.R12[,,"I"],1,sum)/total_pop
I_S2.p <- apply(S2.R12[,,"I"],1,sum)/total_pop
I_S3.p <- apply(S3.R12[,,"I"],1,sum)/total_pop
I_S4.p <- apply(S4.R12[,,"I"],1,sum)/total_pop

#Rt = 1.2
pdf("/Users/tomc/Google Drive/coronavirus/Figures/Rt12-I-turnover.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(I_S1~dates, type="l", col="grey25",lwd=3, xaxt = "n", yaxt="n", cex.main=2,
     xlab="Date", ylab="Population Infected", cex.lab=2, ylim=c(0,45E4), yaxs="i",
     main="Population infectious; post-lockdown Rt = 1.2")
lines(I_S2~dates, col="darkblue", lwd=3)
lines(I_S3~dates, col="darkgreen", lwd=3)
lines(I_S4~dates, col="darkred", lwd=3)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2) #x-axis dates
axis(2, y.vals, labels=y.labs, cex.axis=2) #y-axis scientific notation
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2)
abline(v=c(S1_eq_date, S2_eq_date, S3_eq_date, S4_eq_date), lwd=2, lty=2, 
       col=c("grey25", "darkblue", "darkgreen", "darkred"))
dev.off()

#Plot of I vs S
pdf("/Users/tomc/Google Drive/coronavirus/Figures/phase-plot-Rt12.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(S_S1,I_S1.p, type="l", lwd=3, xlim=c(0.75,1), ylim=c(0,0.0065), col="grey25", 
     cex.lab=2, cex.axis=2, cex.main=2,
     ylab="Proportion Infectious", xlab="Proportion Susceptible", main="Phase-portrait (10 years); post-lockdown Rt = 1.2")
lines(S_S2,I_S2.p, lwd=3, col="darkblue")
lines(S_S3,I_S3.p, lwd=3, col="darkgreen")
lines(S_S4,I_S4.p, lwd=3, col="darkred")
abline(v = 1/1.2, lty=2, lwd=2)
legend(x=c(0.743,0.829), y=c(0.004,0.006), title="Immunity Scenario",
       legend=c(expression(paste("S1  ", omega^N, " = ", infinity, ",       ", omega^H," = ", infinity)), 
                expression(paste("S2  ", omega^N, " = ", 365^-1, ", ", omega^H," = ", infinity)), 
                expression(paste("S3  ", omega^N, " = ", 180^-1, ", ", omega^H," = ", infinity)), 
                expression(paste("S4  ", omega^N, " = ", 90^-1, ",   ", omega^H," = ", 365^-1))), 
       lty=1, lwd=c(3.2,3,3,3.2), cex=1.24,
       col=c("grey25","darkblue","darkgreen", "darkred"))
dev.off()

#incidence of new cases
New_S1 <- apply(S1.R12[,,"new_infections"],1,sum)
New_S2 <- apply(S2.R12[,,"new_infections"],1,sum)
New_S3 <- apply(S3.R12[,,"new_infections"],1,sum)
New_S4 <- apply(S4.R12[,,"new_infections"],1,sum)

plot(New_S1~dates, type="l", lwd=2, ylim=c(0, 140000), col="grey25", xaxt="n",cex.lab=1.5, 
     cex.axis=1.5, cex.main=1.5, ylab="Daily new infections", xlab="Date",
     main="Daily new infections; post-lockdown Rt = 1.2")
lines(New_S2~dates, lwd=2, col="darkblue")
lines(New_S3~dates, lwd=2, col="darkgreen")
lines(New_S4~dates, lwd=2, col="darkred")
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=1.5) #x-axis dates
abline(v=c(S1_eq_date, S2_eq_date, S3_eq_date, S4_eq_date), lwd=2, lty=2, 
       col=c("grey25", "darkblue", "darkgreen", "darkred"))
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)

#Exposed
E_S1 <- apply(S1.R12[,,"E"],1,sum)/total_pop
E_S2 <- apply(S2.R12[,,"E"],1,sum)/total_pop
E_S3 <- apply(S3.R12[,,"E"],1,sum)/total_pop
E_S4 <- apply(S4.R12[,,"E"],1,sum)/total_pop

pdf("/Users/tomc/Google Drive/coronavirus/Figures/Rt12-E-turnover.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(E_S1~dates, type="l", lwd=3, ylim=c(0,0.01), col="grey25", 
     cex.lab=2, cex.axis=2, cex.main=2, xaxt = "n",
     ylab="Proportion Exposed", xlab="Dates", main="Proportion exposed; post-lockdown Rt = 1.2")
lines(E_S2~dates, lwd=3, col="darkblue")
lines(E_S3~dates, lwd=3, col="darkgreen")
lines(E_S4~dates, lwd=3, col="darkred")
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2) #x-axis dates
abline(v=c(S1_eq_date, S2_eq_date, S3_eq_date, S4_eq_date), lwd=2, lty=2, 
       col=c("grey25", "darkblue", "darkgreen", "darkred"))
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
dev.off()

#Susceptible when Rt=1.1
S_S1 <- apply(S1.R11[,,"S"],1,sum)/total_pop
S_S2 <- apply(S2.R11[,,"S"],1,sum)/total_pop
S_S3 <- apply(S3.R11[,,"S"],1,sum)/total_pop
S_S4 <- apply(S4.R11[,,"S"],1,sum)/total_pop

plot(S_S1~dates, type="l", lwd=2, ylim=c(0.9,1), col="grey25")
lines(S_S2~dates, lwd=2, col="darkblue")
lines(S_S3~dates, lwd=2, col="darkgreen")
lines(S_S4~dates, lwd=2, col="darkred")
abline(h = 1/1.1, lty=2)

#for higher values of Rt
source("/Users/tomc/Documents/covid-19-immunity/SEIIRRS-discrete-age-interventions-scenarios-higher-Rt.R")

#Susceptible when Rt=1.3
S_S1 <- apply(S1.R13[,,"S"],1,sum)/total_pop
S_S2 <- apply(S2.R13[,,"S"],1,sum)/total_pop
S_S3 <- apply(S3.R13[,,"S"],1,sum)/total_pop
S_S4 <- apply(S4.R13[,,"S"],1,sum)/total_pop

plot(S_S1~dates, type="l", lwd=2, ylim=c(0.7,1), col="grey25")
lines(S_S2~dates, lwd=2, col="darkblue")
lines(S_S3~dates, lwd=2, col="darkgreen")
lines(S_S4~dates, lwd=2, col="darkred")
abline(h = 1/1.3, lty=2)

