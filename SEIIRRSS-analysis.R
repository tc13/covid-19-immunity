### Dynamics of SARS-CoV-2 with waning immunity
### Thomas Crellen thomas.crellen@bdi.ox.ac.uk, April 2020
### SEIIRRS discrete time model, age structed for UK population, gamma (Erlang) distributed waiting times, two immunity classes
### Analysis of simulations - results and figures 

#Clear R environment
remove(list = ls())

#packages
require(scales)
require(lubridate)
library(dplyr)

#set path to current folder
dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir)

######################################
## UK Population & Contact Matrixes ##
######################################
load("data/BBC_contact_matrices_population_vector.RData")
load("data/uk_population_size_by_age_2018.RData")
total_pop <- sum(uk.pop.2018.count$total)   #Get total UK population
p_age <- uk.pop.2018.count$total/total_pop  #Proportion of population in each age group

####################################################
## Source Model Functions, Parameters & Scenarios ##
####################################################
source("SEIIRRS-functions.R")
source("SEIIRRS-parameters.R")
source("SEIIRRS-scenarios.R") #Takes around 50 mins to run

#Function to make axis numbers in scientific notation
scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

#Date variables
time <- seq(0, days, 1)
date_UK_lockdown <- "2020-03-23"
start_date <- as.Date(date_UK_lockdown)-lockdown_day
hour_period <- 24
dates <- seq(ymd(start_date),ymd(start_date+days), by = as.difftime(hours(hour_period)))

#New infection on lockdown day
new_infections <- apply(S1.R10[,,"new_infections"],1,sum)
new_infections_adults <- apply(S1.R10[,c(5:15),"new_infections"],1,sum)
new_infections[lockdown_day] #107,669 (108,000)

#Total infectious individuals
I <- apply(S1.R10[,,"I"],1,sum)
I[lockdown_day] #75,000

#cumulative infections
#start_cum_date <- "2020-02-16"
#which(dates==start_cum_date)
sum(new_infections[1:which(dates==date_UK_lockdown)]) #559,000
sum(new_infections_adults[1:which(dates==date_UK_lockdown)]) #524,000

#Recovered (ages 19+)
recovered_19plus <- apply(S1.R10[,c(5:15),"R"],1,sum)
end_lockdown <- as.Date(date_UK_lockdown)+(tpi+tfi) 
which(dates==end_lockdown)
recovered_19plus[which(dates==end_lockdown)]/sum(uk.pop.2018.count$total[5:15]) #6.8% of adults with immunity
sum(new_infections[1:which(dates==end_lockdown)])
recovered_all <- apply(S1.R10[,,"R"],1,sum)
recovered_all[which(dates==end_lockdown)]/total_pop #5.7% of the total population with immunity

####################
## Secondary Peak ##
####################

#If Rt=1.2
#In immunity scenario 4
max(apply(S4.R12[,,"I"],1,sum)) #387,000
max(apply(S4.R12[,,"new_infections"],1,sum)) #125,000
#In immunity scenario 1, restrict time to after lockdown- days 130-730
max(apply(S1.R12[,,"I"],1,sum)[130:730]) #126,000
max(apply(S1.R12[,,"new_infections"],1,sum)[130:730]) #41,000

#If Rt=1.1, examine figures in April 2021
april_2021_indx <- which(dates>="2021-04-01"&dates<="2021-04-30")
#In immunity scenario 4
max(apply(S4.R11[,,"I"],1,sum)[april_2021_indx]) #154,000
max(apply(S4.R11[,,"new_infections"],1,sum)[april_2021_indx]) #50,000
#In immunity scenario 1
max(apply(S1.R11[,,"I"],1,sum)[april_2021_indx]) #15,000
max(apply(S1.R11[,,"new_infections"],1,sum)[april_2021_indx]) #5,000

153742.5/14995.98 #nine-fold difference in number infected
49602.97/4661.453 #nine-fold difference in new cases

####################################
## Proportion of Age group immune ##
####################################

#1st October 2020
oct1 <- which(dates=="2020-10-01")
round(S1.R12[oct1,,"R"]/uk.pop.2018.count$total*100, digits = 1)[5:8] #percentage of age group immune (S1)
round(S4.R12[oct1,,"R"]/uk.pop.2018.count$total*100, digits = 1)[5:8] #percentage of age group immune (S4)

##########################
## Infected by Rt value ##
##########################
#axis values
y.vals <- c(0, 5E4, 15E4, 25E4, 35E4, 45E4)
y.labs <- scientific(y.vals)
date_vec <- seq.Date(from=as.Date("2020-03-01"), to=as.Date("2022-03-01"), by="4 months")

#Infection vectors
I_S1 <- apply(S1.R09[,,"I"],1,sum)
I_S2 <- apply(S2.R09[,,"I"],1,sum)
I_S3 <- apply(S3.R09[,,"I"],1,sum)
I_S4 <- apply(S4.R09[,,"I"],1,sum)

#Rt = 0.9
#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-I-Rt-09.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(I_S1~dates, type="l", col="grey40",lwd=4, xaxt = "n", yaxt="n", cex.main=3,
     xlab="Date", ylab="Population Infected", cex.lab=2.5, ylim=c(0,45E4), yaxs="i",
     main=expression(paste("Post-lockdown R"[t]," = 0.9")))
lines(I_S2~dates, col="blue4", lwd=4)
lines(I_S3~dates, col="green4", lwd=4)
lines(I_S4~dates, col="red3", lwd=4)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
axis(2, y.vals, labels=y.labs, cex.axis=2.0) #y-axis scientific notation
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
legend(x=c(as.Date("2020-12-15"), as.Date("2022-03-01")), y=c(242308, 432692), title="Immunity Scenario",
       legend=c("S1: Permanent", "S2: Waning (12 months)", "S3: Waning (6 months)", "S4: Short-lived"), 
       lty=1, lwd=c(3.7,3.5,3.5,3.7), cex=2.2,
       col=c("grey40","blue4","green4", "red3"))
#dev.off()

#Infection vectors
I_S1 <- apply(S1.R10[,,"I"],1,sum)
I_S2 <- apply(S2.R10[,,"I"],1,sum)
I_S3 <- apply(S3.R10[,,"I"],1,sum)
I_S4 <- apply(S4.R10[,,"I"],1,sum)

#Rt = 1.0
#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-I-Rt-10.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(I_S1~dates, type="l", col="grey40",lwd=4, xaxt = "n", yaxt="n", cex.main=3,
     xlab="Date", ylab="Population Infected", cex.lab=2.5, ylim=c(0,45E4), yaxs="i",
     main=expression(paste("Post-lockdown R"[t]," = 1.0")))
lines(I_S2~dates, col="blue4", lwd=4)
lines(I_S3~dates, col="green4", lwd=4)
lines(I_S4~dates, col="red3", lwd=4)
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
plot(I_S1~dates, type="l", col="grey40",lwd=4, xaxt = "n", yaxt="n", cex.main=3,
     xlab="Date", ylab="Population Infected", cex.lab=2.5, ylim=c(0,45E4), yaxs="i",
     main=expression(paste("Post-lockdown R"[t]," = 1.1")))
lines(I_S2~dates, col="blue4", lwd=4)
lines(I_S3~dates, col="green4", lwd=4)
lines(I_S4~dates, col="red3", lwd=4)
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
#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-I-Rt-12.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(I_S1~dates, type="l", col="grey40",lwd=4, xaxt = "n", yaxt="n", cex.main=3,
     xlab="Date", ylab="Population Infected", cex.lab=2.5, ylim=c(0,45E4), yaxs="i",
     main=expression(paste("Post-lockdown R"[t]," = 1.2")))
lines(I_S2~dates, col="blue4", lwd=4)
lines(I_S3~dates, col="green4", lwd=4)
lines(I_S4~dates, col="red3", lwd=4)
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
plot(R_S1~dates, type="l", col="grey40",lwd=4, xaxt = "n", cex.main=3, cex.axis=2.0,
     xlab="Date", ylab="Proportion Immune", cex.lab=2.5, ylim=c(0,0.21), yaxs="i",
     main=expression(paste("Post-lockdown R"[t]," = 0.9")))
lines(R_S2~dates, col="blue4", lwd=4)
lines(R_S3~dates, col="green4", lwd=4)
lines(R_S4~dates, col="red3", lwd=4)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
legend(x=c(as.Date("2020-12-15"), as.Date("2022-03-01")), y=c(0.113, 0.202), title="Immunity Scenario",
       legend=c("S1: Permanent", "S2: Waning (12 months)", "S3: Waning (6 months)", "S4: Short-lived"), 
       lty=1, lwd=c(3.7,3.5,3.5,3.7), cex=2.2,
       col=c("grey40","blue4","green4", "red3"))
#dev.off()

#Recovered vector
R_S1 <- apply(S1.R10[,,"R"],1,sum)/total_pop
R_S2 <- apply(S2.R10[,,"R"],1,sum)/total_pop
R_S3 <- apply(S3.R10[,,"R"],1,sum)/total_pop
R_S4 <- apply(S4.R10[,,"R"],1,sum)/total_pop

#Rt = 1.0
#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-3-R-Rt-10.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(R_S1~dates, type="l", col="grey40",lwd=4, xaxt = "n", cex.main=3, cex.axis=2.0,
     xlab="Date", ylab="Proportion Immune", cex.lab=2.5, ylim=c(0,0.21), yaxs="i",
     main=expression(paste("Post-lockdown R"[t]," = 1.0")))
lines(R_S2~dates, col="blue4", lwd=4)
lines(R_S3~dates, col="green4", lwd=4)
lines(R_S4~dates, col="red3", lwd=4)
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
plot(R_S1~dates, type="l", col="grey40",lwd=4, xaxt = "n", cex.main=3, cex.axis=2.0,
     xlab="Date", ylab="Proportion Immune", cex.lab=2.5, ylim=c(0,0.21), yaxs="i",
     main=expression(paste("Post-lockdown R"[t]," = 1.1")))
lines(R_S2~dates, col="blue4", lwd=4)
lines(R_S3~dates, col="green4", lwd=4)
lines(R_S4~dates, col="red3", lwd=4)
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
plot(R_S1~dates, type="l", col="grey40",lwd=4, xaxt = "n", cex.main=3, cex.axis=2.0,
     xlab="Date", ylab="Proportion Immune", cex.lab=2.5, ylim=c(0,0.21), yaxs="i",
     main=expression(paste("Post-lockdown R"[t]," = 1.2")))
lines(R_S2~dates, col="blue4", lwd=4)
lines(R_S3~dates, col="green4", lwd=4)
lines(R_S4~dates, col="red3", lwd=4)
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
#dev.off()

#########################################
## Immunity by age group when Rt = 1.2 ##
#########################################
#Binning age groups into 10 years- immune
S1.R12.AG <- list()
S1.R12.AG[[1]] <- (S1.R12[,1,"R"]+S1.R12[,2,"R"])/(uk.pop.2018.count$total[1]+uk.pop.2018.count$total[2])
S1.R12.AG[[2]] <- (S1.R12[,3,"R"]+S1.R12[,4,"R"])/(uk.pop.2018.count$total[3]+uk.pop.2018.count$total[4])
S1.R12.AG[[3]] <- (S1.R12[,5,"R"]+S1.R12[,6,"R"])/(uk.pop.2018.count$total[5]+uk.pop.2018.count$total[6])
S1.R12.AG[[4]]<- (S1.R12[,7,"R"]+S1.R12[,8,"R"])/(uk.pop.2018.count$total[7]+uk.pop.2018.count$total[8])
S1.R12.AG[[5]]<- (S1.R12[,9,"R"]+S1.R12[,10,"R"])/(uk.pop.2018.count$total[9]+uk.pop.2018.count$total[10])
S1.R12.AG[[6]] <- (S1.R12[,11,"R"]+S1.R12[,12,"R"])/(uk.pop.2018.count$total[11]+uk.pop.2018.count$total[12])
S1.R12.AG[[7]] <- (S1.R12[,13,"R"]+S1.R12[,14,"R"])/(uk.pop.2018.count$total[13]+uk.pop.2018.count$total[14])
S1.R12.AG[[8]] <- S1.R12[,15,"R"]/uk.pop.2018.count$total[15]

#colour vectors
col.vec.8 <- c(rep("grey50",2),"darkorange","forestgreen", rep("grey50",3), "darkviolet")
lwd.vec.8 <- c(rep(1.5,2), 3.5, 3.5, rep(1.5,3), 3.5)

#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-S1-R12-R-age-groups.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0))
plot(S1.R12.AG[[1]]~dates, type="l", ylim=c(0,0.35), 
     col=col.vec.8[1], lwd=lwd.vec.8[1], xaxt = "n", cex.main=3, cex.axis=2.0,
     xlab="Date", ylab="Proportion of age group immune", cex.lab=2.5, yaxs="i",
     main=expression(paste("Immunity Scenario 1, Post-lockdown R"[t]," = 1.2")))
for(a in 2:8){
  lines(S1.R12.AG[[a]]~dates, col=col.vec.8[a], lwd=lwd.vec.8[a])
}
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
legend(x=c(as.Date("2020-07-01"), as.Date("2021-01-01")),
       y=c(0.24, 0.34), title = "Age groups", cex=1.5,
       legend = c("20-29 years", "30-39 years", "70+ years", "Others"),
       lwd=c(3.2,3.2,3.2,3,2), col=c("darkorange","forestgreen", "darkviolet", "grey50"))
#dev.off()

#Binning age groups into 10 years- immune
S4.R12.AG <- list()
S4.R12.AG[[1]] <- (S4.R12[,1,"R"]+S4.R12[,2,"R"])/(uk.pop.2018.count$total[1]+uk.pop.2018.count$total[2])
S4.R12.AG[[2]] <- (S4.R12[,3,"R"]+S4.R12[,4,"R"])/(uk.pop.2018.count$total[3]+uk.pop.2018.count$total[4])
S4.R12.AG[[3]] <- (S4.R12[,5,"R"]+S4.R12[,6,"R"])/(uk.pop.2018.count$total[5]+uk.pop.2018.count$total[6])
S4.R12.AG[[4]]<- (S4.R12[,7,"R"]+S4.R12[,8,"R"])/(uk.pop.2018.count$total[7]+uk.pop.2018.count$total[8])
S4.R12.AG[[5]]<- (S4.R12[,9,"R"]+S4.R12[,10,"R"])/(uk.pop.2018.count$total[9]+uk.pop.2018.count$total[10])
S4.R12.AG[[6]] <- (S4.R12[,11,"R"]+S4.R12[,12,"R"])/(uk.pop.2018.count$total[11]+uk.pop.2018.count$total[12])
S4.R12.AG[[7]] <- (S4.R12[,13,"R"]+S4.R12[,14,"R"])/(uk.pop.2018.count$total[13]+uk.pop.2018.count$total[14])
S4.R12.AG[[8]] <- S4.R12[,15,"R"]/uk.pop.2018.count$total[15]

#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-S4-R12-R-age-groups.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0))
plot(S4.R12.AG[[1]]~dates, type="l", ylim=c(0,0.35), 
     col=col.vec.8[1], lwd=lwd.vec.8[1], xaxt = "n", cex.main=3, cex.axis=2.0,
     xlab="Date", ylab="Proportion of age group immune", cex.lab=2.5, yaxs="i",
     main=expression(paste("Immunity Scenario 4, Post-lockdown R"[t]," = 1.2")))
for(a in 2:8){
  lines(S4.R12.AG[[a]]~dates, col=col.vec.8[a], lwd=lwd.vec.8[a])
}
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
#dev.off()

#Binning age groups into 10 years- infected
S1.R12.AG[[9]] <- (S1.R12[,1,"I"]+S1.R12[,2,"I"])/(uk.pop.2018.count$total[1]+uk.pop.2018.count$total[2])
S1.R12.AG[[10]] <- (S1.R12[,3,"I"]+S1.R12[,4,"I"])/(uk.pop.2018.count$total[3]+uk.pop.2018.count$total[4])
S1.R12.AG[[11]] <- (S1.R12[,5,"I"]+S1.R12[,6,"I"])/(uk.pop.2018.count$total[5]+uk.pop.2018.count$total[6])
S1.R12.AG[[12]]<- (S1.R12[,7,"I"]+S1.R12[,8,"I"])/(uk.pop.2018.count$total[7]+uk.pop.2018.count$total[8])
S1.R12.AG[[13]]<- (S1.R12[,9,"I"]+S1.R12[,10,"I"])/(uk.pop.2018.count$total[9]+uk.pop.2018.count$total[10])
S1.R12.AG[[14]] <- (S1.R12[,11,"I"]+S1.R12[,12,"I"])/(uk.pop.2018.count$total[11]+uk.pop.2018.count$total[12])
S1.R12.AG[[15]] <- (S1.R12[,13,"I"]+S1.R12[,14,"I"])/(uk.pop.2018.count$total[13]+uk.pop.2018.count$total[14])
S1.R12.AG[[16]] <- S1.R12[,15,"I"]/uk.pop.2018.count$total[15]

#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-S1-R12-I-age-groups.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0))
plot(S1.R12.AG[[9]]~dates, type="l", ylim=c(0, 0.01), 
     col=col.vec.8[1], lwd=lwd.vec.8[1], xaxt = "n", cex.main=3, cex.axis=2.0, yaxs="i",
     xlab="Date", ylab="Proportion of age group infected", cex.lab=2.5,
     main=expression(paste("Immunity Scenario 1, Post-lockdown R"[t]," = 1.2")))
for(a in 2:8){
  lines(S1.R12.AG[[(a+8)]]~dates, col=col.vec.8[a], lwd=lwd.vec.8[a])
}
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
legend(x=c(as.Date("2020-07-01"), as.Date("2021-01-01")),
       y=c(0.006857143, 0.009714286), title = "Age groups", cex=1.5,
       legend = c("20-29 years", "30-39 years", "70+ years", "Others"),
       lwd=c(3.2,3.2,3.2,3,2), col=c("darkorange","forestgreen", "darkviolet", "grey50"))
#dev.off()

#Binning age groups into 10 years- infected
S4.R12.AG[[9]] <- (S4.R12[,1,"I"]+S4.R12[,2,"I"])/(uk.pop.2018.count$total[1]+uk.pop.2018.count$total[2])
S4.R12.AG[[10]] <- (S4.R12[,3,"I"]+S4.R12[,4,"I"])/(uk.pop.2018.count$total[3]+uk.pop.2018.count$total[4])
S4.R12.AG[[11]] <- (S4.R12[,5,"I"]+S4.R12[,6,"I"])/(uk.pop.2018.count$total[5]+uk.pop.2018.count$total[6])
S4.R12.AG[[12]]<- (S4.R12[,7,"I"]+S4.R12[,8,"I"])/(uk.pop.2018.count$total[7]+uk.pop.2018.count$total[8])
S4.R12.AG[[13]]<- (S4.R12[,9,"I"]+S4.R12[,10,"I"])/(uk.pop.2018.count$total[9]+uk.pop.2018.count$total[10])
S4.R12.AG[[14]] <- (S4.R12[,11,"I"]+S4.R12[,12,"I"])/(uk.pop.2018.count$total[11]+uk.pop.2018.count$total[12])
S4.R12.AG[[15]] <- (S4.R12[,13,"I"]+S4.R12[,14,"I"])/(uk.pop.2018.count$total[13]+uk.pop.2018.count$total[14])
S4.R12.AG[[16]] <- S4.R12[,15,"I"]/uk.pop.2018.count$total[15]

#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-S4-R12-I-age-groups.pdf", height=8, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0))
plot(S4.R12.AG[[9]]~dates, type="l", ylim=c(0, 0.01), 
     col=col.vec.8[1], lwd=lwd.vec.8[1], xaxt = "n", cex.main=3, cex.axis=2.0, yaxs="i",
     xlab="Date", ylab="Proportion of age group infected", cex.lab=2.5,
     main=expression(paste("Immunity Scenario 4, Post-lockdown R"[t]," = 1.2")))
for(a in 2:8){
  lines(S4.R12.AG[[(a+8)]]~dates, col=col.vec.8[a], lwd=lwd.vec.8[a])
}
axis(1, date_vec, format(date_vec, "%b-%y"), cex.axis=2.0) #x-axis dates
abline(v=c(as.Date(date_UK_lockdown), end_lockdown), lty=2, lwd=2.5)
#dev.off()
