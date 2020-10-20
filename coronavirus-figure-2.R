#Figures for coronavirus modelling immunity paper

#Clear R environment
remove(list = ls())

#libraries
require(reshape2)
require(ggplot2)
require(lubridate)
require(scales)

#set path to parent folder
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
source("SEIIRRS-scenarios.R")

#ggplot theme set
theme_set(theme_bw(base_size = 28))

#Function to make axis numbers in scientific notation
scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

#Next generation NGM
NGM <- function(C, beta, R0, mean_infectious, 
                            susceptibility, prop_asymtomatic, 
                            asymtomatic_infectiousness, p_age){ 
  
  #Calculate next generation matrix K
  K = matrix(data=0, nrow=nrow(C), ncol=ncol(C))
  #Kij = infections in group i produced by individuals in group j
  for(i in 1:nrow(K)){  
    for(j in 1:nrow(K)){
      K[i,j] = beta*mean_infectious*(prop_asymtomatic[j]*asymtomatic_infectiousness*C[i,j]+C[i,j]*(1-prop_asymtomatic[j]))*susceptibility[i]
    }
  }
  #Dominant eigenvalue of K
  dominant = get_eigen(K)
  
  test_that("NGM=R0", expect_equal(dominant, R0))
  return(K)
}

#Extra params
o=2
immune_mean_1=360
immune_mean_2=180
immune_mean_3=90

#density plots
curve(dgamma(x, shape=m, rate= m/sigma_recip), from=0,to=20)
curve(dgamma(x, shape=n, rate=n/gamma_recip), from=0,to=20)
curve(dgamma(x, shape=o, rate=o/immune_mean_1), from=0,to=365*2)
curve(dgamma(x, shape=o, rate=o/immune_mean_2), from=0,to=365*2)
curve(dgamma(x, shape=o, rate=o/immune_mean_3), from=0,to=365*2)

#simulation
N <- 100000
set.seed(1234)
latent_periods <- rgamma(n=N, shape=m, rate=m/sigma_recip)
set.seed(1234)
infectious_periods <- rgamma(n=N, shape=n, rate=m/gamma_recip)
max_day = max(latent_periods+infectious_periods) #23
n_latent = n_infectious = c()
n_latent[1] = N
n_infectious[1] = 0
for(i in 1:max_day){
  #still latent
  n_latent[(i+1)] <- length(which(latent_periods>=i))
  #no longer latent
  inft_indx <- which(latent_periods<i)
  #still infectious
  n_infectious[(i+1)] <- length(which((infectious_periods[inft_indx]+latent_periods[inft_indx])>=i))
}
time <- seq(1,max_day+1,1)
period.df <- data.frame(time=time, latent=n_latent, infectious=n_infectious)
period.m <- melt(data = period.df, id.vars = "time")

#Base R
#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-E-I-state-plots.pdf", height=9, width=11)
par(family = "serif", mar=c(5.0, 5.0, 3.0, 2.0))
plot(n_latent/N~time, type="l", lwd=3, xlab="Time since infection (days)", ylim=c(0.02,0.98), 
     ylab="Proportion", xlim=c(1,20), cex.main=2, cex.lab=2.5, cex.axis=2)
lines(n_infectious/N~time, lwd=3, lty=4)
legend(x=c(14,20), y=c(0.75, 1), legend = c("Exposed", "Infectious"), 
       lty=c(1,4), lwd=c(3,3), title="Disease State", cex=2)
#dev.off()

#waning immunity
x <- seq(0,365*4,by=0.5)
im.90 <- sapply(x, FUN=function(x) 1-pgamma(q = x, shape = 2, rate=2/90))
im.180 <- sapply(x, FUN=function(x) 1-pgamma(q = x, shape = 2, rate=2/180))
im.365 <- sapply(x, FUN=function(x) 1-pgamma(q = x, shape = 2, rate=2/365))

#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Fig-immune-wane.pdf", height=9, width=11)
par(family = "serif", mar=c(5.0, 4.6, 4.0, 2.0)) #set font and plotting margin
plot(im.365~x, type="l", lwd=3, xaxt="n", cex.lab=2.5, lty=2, ylim=c(0.02,0.98),
     xlab="Time since recovery (months)", ylab="Proportion immune", cex.axis=2.0, xlim=c(30,1430))
lines(im.180 ~ x, lwd=3, lty=3)
lines(im.90 ~ x, lwd=3, lty=1)
axis(1, at = seq(from=0,to=1500, by=182), 
     labels=as.character(seq(0,48,6)), cex.axis=2.0) #x-axis dates
legend(x=c(25*30,47*30),y=c(0.70,1), 
       legend=c("90 days", "180 days", "365 days"), 
       lty=c(1,3,2), lwd=c(3,3,3), title="Mean duration of immunity", cex=2)
#dev.off()

##Contact matrixes
C0 <- make.intervention.matrix(BBC_contact_matrix, intervention=c(1,1,1,1)) #baseline
C1 <- make.intervention.matrix(BBC_contact_matrix, intervention=c(home=0.8,work=0.3,school=0.1,other=0.2)) #full lockdown
C2 <- make.intervention.matrix(BBC_contact_matrix, intervention=c(home=1,work=0.8,school=0.85, other=0.75)) #post lockdown

#Get beta
beta_baseline <- get_beta(C=C0, R0=R0, mean_infectious=gamma_recip, susceptibility=susceptibility, 
                 prop_asymtomatic=phi, asymtomatic_infectiousness=ari, p_age=p_age)
NGM_baseline <- NGM(C=C0, beta=beta_baseline, R0=R0, mean_infectious = gamma_recip, susceptibility = susceptibility,
           prop_asymtomatic=phi, asymtomatic_infectiousness = ari, p_age=p_age)

beta_intervention <- get_beta(C=C1, R0=0.8, mean_infectious=gamma_recip, susceptibility=susceptibility, 
                              prop_asymtomatic=phi, asymtomatic_infectiousness=ari, p_age=p_age)
NGM_intervention <- NGM(C=C1, beta=beta_intervention, R0=0.8, mean_infectious = gamma_recip, susceptibility = susceptibility,
                        prop_asymtomatic=phi, asymtomatic_infectiousness = ari, p_age=p_age)

beta_post_09 <- get_beta(C=C2, R0=0.9, mean_infectious=gamma_recip, susceptibility=susceptibility, 
                         prop_asymtomatic=phi, asymtomatic_infectiousness=ari, p_age=p_age)

beta_post_10 <- get_beta(C=C2, R0=1.0, mean_infectious=gamma_recip, susceptibility=susceptibility, 
                         prop_asymtomatic=phi, asymtomatic_infectiousness=ari, p_age=p_age)

beta_post_11 <- get_beta(C=C2, R0=1.1, mean_infectious=gamma_recip, susceptibility=susceptibility, 
                         prop_asymtomatic=phi, asymtomatic_infectiousness=ari, p_age=p_age)

beta_post_12 <- get_beta(C=C2, R0=1.2, mean_infectious=gamma_recip, susceptibility=susceptibility, 
                         prop_asymtomatic=phi, asymtomatic_infectiousness=ari, p_age=p_age)

NGM_post_09 <- NGM(C=C2, beta=beta_post_09, R0=0.9, mean_infectious = gamma_recip, susceptibility = susceptibility,
                   prop_asymtomatic=phi, asymtomatic_infectiousness = ari, p_age=p_age)

NGM_post_10 <- NGM(C=C2, beta=beta_post_10, R0=1.0, mean_infectious = gamma_recip, susceptibility = susceptibility,
                   prop_asymtomatic=phi, asymtomatic_infectiousness = ari, p_age=p_age)

NGM_post_11 <- NGM(C=C2, beta=beta_post_11, R0=1.1, mean_infectious = gamma_recip, susceptibility = susceptibility,
                   prop_asymtomatic=phi, asymtomatic_infectiousness = ari, p_age=p_age)

NGM_post_12 <- NGM(C=C2, beta=beta_post_12, R0=1.2, mean_infectious = gamma_recip, susceptibility = susceptibility,
                   prop_asymtomatic=phi, asymtomatic_infectiousness = ari, p_age=p_age)

##Figure 2A - Population of UK by age group
x.lb <- c(0, 2000000, 4000000, 6000000, 8000000, 10000000)
xlabs_bar <- scientific(x.lb)

#pdf("/Users/tomc/Google Drive/coronavirus/Figures/Figure2A.pdf", height=7, width=10)
barplot(uk.pop.2018.count$total, xlab="Age Group", names.arg=uk.pop.2018.count$age_group, 
        ylim=c(0,10000000), yaxt = "n",
        ylab="Population Size", cex.axis = 1.4, cex.lab=1.4, cex.names = 1.4)
axis(side = 2, at = x.lb, labels = xlabs_bar, cex.axis=1.4)
#dev.off()

#Heatmap from contacts
#work
work <- melt(BBC_contact_matrix$all_work)
colnames(work) <- c("Var1", "Var2", "value")
work$location <- "work"
work$time <- "baseline"

work_intervention <- work
work_intervention$value <- work$value*0.3
work_intervention$location <- "work"
work_intervention$time <- "intervention"

work_post_intervention <- work
work_post_intervention$value <- work$value*0.8
work_post_intervention$location <- "work"
work_post_intervention$time <- "post_intervention"

#school
school <- melt(BBC_contact_matrix$all_school)
colnames(school) <- c("Var1", "Var2", "value")
school$location <- "school"
school$time <- "baseline"

school_intervention <- school
school_intervention$value <- school$value*0.1
school_intervention$location <- "school"
school_intervention$time <- "intervention"

school_post_intervention <- school
school_post_intervention$value <- school$value*0.85
school_post_intervention$location <- "school"
school_post_intervention$time <- "post_intervention"

#home
home <- melt(BBC_contact_matrix$all_home)
colnames(home) <- c("Var1", "Var2", "value")
home$location <- "home"
home$time <- "baseline"

home_intervention <- home
home_intervention$value <- home$value*0.8
home_intervention$location <- "home"
home_intervention$time <- "intervention"

home_post_intervention <- home
home_post_intervention$value <- home$value*1
home_post_intervention$location <- "home"
home_post_intervention$time <- "post_intervention"

#other
other <- melt(BBC_contact_matrix$all_other)
colnames(other) <- c("Var1", "Var2", "value")
other$location <- "other"
other$time <- "baseline"

other_intervention <- other
other_intervention$value <- other$value*0.2
other_intervention$location <- "other"
other_intervention$time <- "intervention"

other_post_intervention <- other
other_post_intervention$value <- other$value*0.75
other_post_intervention$location <- "other"
other_post_intervention$time <- "post_intervention"

#all
all <- melt(C0)
all$location <- "all"
all$time <- "baseline"

all_intervention <- melt(C1)
all_intervention$location <- "all"
all_intervention$time <- "intervention"

all_post_intervention <- melt(C2)
all_post_intervention$location <- "all"
all_post_intervention$time <- "post_intervention"

#contact matrixes combined
contacts <- rbind(home, home_intervention, home_post_intervention,
                  work, work_intervention, work_post_intervention,
                  school, school_intervention, school_post_intervention,
                  other, other_intervention, other_post_intervention,
                  all, all_intervention, all_post_intervention)

ggplot(data = contacts, aes(x=Var1, y=Var2, fill=value)) + 
  geom_raster(hjust=1, vjust=1)+
  scale_fill_distiller(type="seq", palette=1, direction = 1, values=c(0.075,1), na.value ="#ffffff")+
  theme(legend.position = "none", text=element_text(family = "Times"))+
  xlab("Age of Individual (years)")+ylab("Age of Individual (years)")+
  scale_x_discrete(breaks=c("0-4", "10-14", "20-24", "30-34", "40-44", "50-54", "60-64", "70+"), 
                   labels=c("0", "10", "20", "30", "40", "50", "60", "70"), expand=c(0,0))+
  scale_y_discrete(breaks=c("0-4", "10-14", "20-24", "30-34", "40-44", "50-54", "60-64", "70+"), 
                   labels=c("0", "10", "20", "30", "40", "50", "60", "70"), expand=c(0,0))+
  facet_grid(time~location)

#Melt NGMS
NGM_baseline_m <- melt(NGM_baseline)
NGM_baseline_m$time <- "Baseline~(R[0]==2.8)"

NGM_intervention_m <- melt(NGM_intervention)
NGM_intervention_m$time <- "Lockdown~(R[t]==0.8)"
  
NGM_post_09_m <- melt(NGM_post_09)
NGM_post_09_m$time <- "Post-lockdown~(R[t]==0.9)"

NGM_post_10_m <- melt(NGM_post_10)
NGM_post_10_m$time <- "Post-lockdown~(R[t]==1.0)"

NGM_post_11_m <- melt(NGM_post_11)
NGM_post_11_m$time <- "Post-lockdown~(R[t]==1.1)"

NGM_post_12_m <- melt(NGM_post_12)
NGM_post_12_m$time <- "Post-lockdown~(R[t]==1.2)"

#Next generation matrixes combined
NGMs <- rbind(NGM_baseline_m, NGM_intervention_m, NGM_post_09_m, 
              NGM_post_10_m, NGM_post_11_m, NGM_post_12_m)

#png("/Users/tomc/Google Drive/coronavirus/Figures/Fig-NGM.png", width=14, height=11.5, units="in", res=450)
ggplot(data = NGMs, aes(x=Var1, y=Var2, fill=value))+ 
  geom_raster(hjust=1, vjust=1)+
  scale_fill_distiller(name="Secondary cases from index case", breaks=c(0, 0.25, 0.5, 0.75, 1), 
                       labels=c("0", "0.25", "0.5", "0.75", "1"),
  limit=c(0,1),type="seq", palette=1, direction = 1, values=c(0.1,1), na.value ="#ffffff",
    guide=guide_colourbar(title.position="top", barwidth = 12, frame.colour="black", 
                          label.theme=element_text(size=16, family="Times"),
                          ticks.colour = "black", ticks.linewidth = 2, label.vjust=2.5, title.vjust = 1))+ 
  theme(legend.position = "top", text=element_text(family = "Times"), legend.title = element_text(size = 18), 
        legend.box.spacing =grid::unit(-0.5, "cm"), legend.justification='left', 
        strip.background = element_rect(fill="white", color="white"))+
  ylab("Age of index case (years)")+xlab("Age of secondary cases (years)")+
  scale_x_continuous(breaks=c(1, 3, 5, 7, 9, 11, 13, 15), 
                   labels=c("0", "10", "20", "30", "40", "50", "60", "70"), expand=c(0,0))+
  scale_y_continuous(breaks=c(1, 3, 5, 7, 9, 11, 13, 15), 
                   labels=c("0", "10", "20", "30", "40", "50", "60", "70"), expand=c(0,0))+
  facet_wrap(time~., labeller = label_parsed)
#dev.off()
