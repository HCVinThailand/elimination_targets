##############################
# HEP C AGE STRUCTURED MODEL #
##############################

# Cleanup - CAUTION clears R environment ####
rm(list=ls())

# Setup ####
library(pacman)
#library(MIMSunit)
p_load(deSolve, tidyverse, doParallel, manipulate, readxl, gridExtra, grid, scales, tictoc, Hmisc, viridis)

#setwd("~/MGH/thailand/dissertation/")

scaleFUN <- function(x) sprintf("%.2f", x)
scaleFUN2 <- function(x) sprintf("%.4f", x)

age_group_vector<-c('0 - 4', '5 - 9', '10 - 14', '15 - 19', '20 - 24', '25 - 29', '30 - 34', '35 -39',
                    '40 - 44', '45 - 49', '50 - 54', '55 - 59', '60 - 64', '65 - 69', '70 - 74', '75 - 79',
                    '80 - 84', '85 - 89', '90 - 94', '95 - 99', '100 and over')

get_only_legend <- function(plot) {
  # get tabular interpretation of plot
  plot_table <- ggplot_gtable(ggplot_build(plot)) 
  #  Mark only legend in plot
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box") 
  # extract legend
  legend <- plot_table$grobs[[legend_plot]]
  # return legend
  return(legend) 
}

`%nin%` = Negate(`%in%`)

# Read Excel data ####

# HCV deaths data
hcv_deaths_data <- as.data.frame(read_excel("data/pop_structure_hcv.xlsx", sheet="HCV_deaths", range="A1:D9", col_names=TRUE)) #%>% 

# Read mortality and population structure data

scenario_base <- as_tibble(as.data.frame(read_excel("data/pop_structure_hcv.xlsx", sheet="mortality_scenarios", range="C3:W3", col_names=FALSE)))
scenario_1 <- as_tibble(as.data.frame(read_excel("data/pop_structure_hcv.xlsx", sheet="mortality_scenarios", range="C4:W4", col_names=FALSE)))
scenario_2 <- as_tibble(as.data.frame(read_excel("data/pop_structure_hcv.xlsx", sheet="mortality_scenarios", range="C5:W5", col_names=FALSE)))
scenario_3 <- as_tibble(as.data.frame(read_excel("data/pop_structure_hcv.xlsx", sheet="mortality_scenarios", range="C6:W6", col_names=FALSE)))

# Contact matrix data
contact0 <- as.data.frame(read_excel("data/pop_structure_hcv.xlsx", sheet="sexual_contact_matrix", range="D3:X23", col_names=FALSE))

# Initial conditions values
init_condits <- as_tibble(as.data.frame(read_excel("data/pop_structure_hcv.xlsx", sheet="initial_conditions_new", range="B5:V31", col_names=FALSE)))

# Birth rate values
birth.approx <- as_tibble(as.data.frame(read_excel("data/pop_structure_hcv.xlsx", sheet="birthrate", range="A1:B38", col_names=TRUE)))

# Reading mortality baseline (empty from 2023 onwards)
mortality_base_empty <- as_tibble(as.data.frame(read_excel("data/pop_structure_hcv.xlsx", sheet="mortality_rates", range="B3:W40", col_names=TRUE)))

# Population structure data
age_struc <- as_tibble(as.data.frame(read_excel("data/pop_structure_hcv.xlsx", sheet="Population_data_proportion", range="A2:V20", col_names=TRUE)))

# Population and birth rate data and projection - UN
pop_data <- cbind(c(2004:2021),as_tibble(as.data.frame(read_excel("data/pop_structure_hcv.xlsx", sheet="Population_data_absolute", range="W3:W20", col_names=FALSE))))
pop_proj <- cbind(c(2022:2040),as_tibble(as.data.frame(read_excel("data/pop_structure_hcv.xlsx", sheet="Population_data_absolute", range="W21:W39", col_names=FALSE))))

pop_data1 <- cbind(c(2004:2040),as_tibble(as.data.frame(read_excel("data/pop_structure_hcv.xlsx", sheet="Population_data_absolute", range="W3:X39", col_names=FALSE))))

# Population structure data and projection
age_struc_proportion_data <- cbind(2004:2021,as_tibble(as.data.frame(read_excel("data/pop_structure_hcv.xlsx", sheet="Population_data_absolute", range="B2:V20", col_names=TRUE))))
age_struc_proportion_proj <- cbind(2004:2040,as_tibble(as.data.frame(read_excel("data/pop_structure_hcv.xlsx", sheet="Population_data_absolute", range="B2:V39", col_names=TRUE))))

# Prevalence data
data_prev <- as_tibble(as.data.frame(read_excel("data/pop_structure_hcv.xlsx", sheet="prevalence_data_by_age", range="B1:G16", col_names=TRUE)))

# Aging process ####

age_groups <-  c(seq(1,21,length.out=21))
groups <- length(age_groups)
da <- rep(5,groups) # difference between consecutive age groups in Years

aging <- diag(-1/da)
aging[row(aging)-col(aging)==1] <- 1/head(da,-1) # time in years to enter and leave each age group 


# Contact matrix ####

contact0[contact0==0] <- min(contact0[contact0!=0])/50

beta <- matrix(NA,nrow=groups,ncol=groups)

beta_multiplier <- 0.025

for(j in 1:groups){
  for(i in 1:groups){
    beta[i,j] <- beta_multiplier*contact0[i,j]
  }
}

# Screening by age group ####

scr_start <- 19 # screening programme start Year 19=2023
scr_dur <- 7 # screening programme duration

# Baseline screening coverage

cov_mean <- 0.07
cov_sd <- 0.03
cov_lower <- cov_mean + 1.96*cov_sd
cov_upper <- cov_mean - 1.96*cov_sd

groupBASE <- rep(cov_mean,groups)
groupBASE_lower <- rep(cov_lower,groups)
groupBASE_upper <- rep(cov_upper,groups)

# Targeted screening strategies 

coverageBASE <- cov_mean
coverage1 <- 0.5
coverage2 <- 0.9

groupA <- c(rep(coverageBASE, 6), rep(coverage1, 2), rep(coverageBASE,13)) # 30 - 39 at 50% per year
groupB <- c(rep(coverageBASE, 6), rep(coverage2, 2), rep(coverageBASE,13)) # 30 - 39 at 90% per year
groupC <- c(rep(coverageBASE, 8), rep(coverage1, 2), rep(coverageBASE,11)) # 40 - 49 at 50% per year
groupD <- c(rep(coverageBASE, 8), rep(coverage2, 2), rep(coverageBASE,11)) # 40 - 49 at 90% per year
groupE <- c(rep(coverageBASE, 10), rep(coverage1, 2), rep(coverageBASE,9)) # 50 - 59 at 50% per year
groupF <- c(rep(coverageBASE, 6), rep(coverage1, 6), rep(coverageBASE,9)) # 30 - 59 at 50% per year

# Set up empty indices for storing results ####

Sindex <- 1:groups
F0index <- (groups+1):(2*groups)
F1index <- (2*groups+1):(3*groups)
F2index <- (3*groups+1):(4*groups)
F3index <- (4*groups+1):(5*groups)
C1index <- (5*groups+1):(6*groups)
C2index <- (6*groups+1):(7*groups)
C3index <- (7*groups+1):(8*groups)
C4index <- (8*groups+1):(9*groups)

HCCAindex <- (9*groups+1):(10*groups)
HCCBindex <- (10*groups+1):(11*groups)
HCCCindex <- (11*groups+1):(12*groups)
HCCDindex <- (12*groups+1):(13*groups)

F0cureindex <- (13*groups+1):(14*groups)
F1cureindex <- (14*groups+1):(15*groups)
F2cureindex <- (15*groups+1):(16*groups)
F3cureindex <- (16*groups+1):(17*groups)
C1cureindex <- (17*groups+1):(18*groups)
C2cureindex <- (18*groups+1):(19*groups)
C3cureindex <- (19*groups+1):(20*groups)
C4cureindex <- (20*groups+1):(21*groups)

Dindex <- (21*groups+1):(22*groups)
dthC14index <- (22*groups+1):(23*groups)
dthHCCindex <- (23*groups+1):(24*groups)

translivindex <- (24*groups+1):(25*groups)

aliveindex <- (1:21*groups)
infectindex <- (groups+1):(13*groups)
totalHCCindex <- (9*groups+1):(13*groups)
totalHCVindex <- (groups+1):(9*groups)

newdeathindex <- (25*groups+1):(26*groups)
incHCCindex <- (26*groups+1):(27*groups)

CIncindex <- (27*groups+1):(28*groups)
CScrindex <- (28*groups+1):(29*groups)

infectindex <- (groups+1):(13*groups)

totalHCCindex <- (9*groups+1):(13*groups)
totalHCVindex <- (groups+1):(9*groups)


# Parameters ####

start_year <- 2004
# Compartment flow rates
f0f1 <- 0.117
f1f2 <- 0.085
f2f3 <- 0.12
f3c1 <- 0.116
c1c2 <- 0.044
c2c3 <- 0.044
c3c4 <- 0.076
c1bA <- 0.0068
c1bB <- 0.0099
c1bC <- 0.0029
c1bD <- 0.0068
c2bA <- 0.0068
c2bB <- 0.0099
c2bC <- 0.0029
c2bD <- 0.0068
c3bD <- 0.0664
c4bD <- 0.0664
# Compartment mortality rates
dthc1 <- 0.01 
dthc2 <- 0.01
dthc3 <- 0.2
dthc4 <- 0.57
dthbA <- 1/(36/12) # 3 years to die from HCCA
dthbB <- 1/(16/12) # 1.3 years to die from HCCB
dthbC <- 1/(6/12) # 6 months to die from HCCC
dthbD <- 1/(3/12) # 3 months to die from HCCD
dthtrn <- 1/(240/12)
tranc4 <- 0.0015 # liver transplant rates
tranbA <- 0.0015
tranbB <- 0.0015
trt_start <- 15 # DAA treatment starting in 2019
sens <- 0.985 # sensitivity - combined

pF0scr<-1 # proportion of compartment screened
pF1scr<-1
pF2scr<-1
pF3scr<-1
pC1scr<-1
pC2scr<-1
pC3scr<-1
pC4scr<-1

std_cureF0<-0.7 # efficacy of standard treatment
std_cureF1<-0.7
std_cureF2<-0.7
std_cureF3<-0.7
std_cureC1<-0.7
new_cureF0<-0.985
new_cureF1<-0.985 # efficacy of new treatment - proportion cured
new_cureF2<-0.985
new_cureF3<-0.985
new_cureC1<-0.985
new_cureC2<-0.985
new_cureC3<-0.985
new_cureC4<-0.985

# set up parameters
parameters<- c(
  start_year,
  f0f1,
  f1f2,
  f2f3,
  f3c1,
  c1c2,
  c2c3,
  c3c4,
  c1bA,
  c1bB,
  c1bC,
  c1bD,
  c2bA,
  c2bB,
  c2bC,
  c2bD,
  c3bD,
  c4bD,
  dthc1,
  dthc2,
  dthc3,
  dthc4,
  dthbA,
  dthbB,
  dthbC,
  dthbD,
  dthtrn,
  tranc4,
  tranbA,
  tranbB,
  trt_start,
  sens,
  pF0scr,
  pF1scr,
  pF2scr,
  pF3scr,
  pC1scr,
  pC2scr,
  pC3scr,
  pC4scr,
  std_cureF0,
  std_cureF1,
  std_cureF2,
  std_cureF3,
  std_cureC1,
  new_cureF0,
  new_cureF1,
  new_cureF2,
  new_cureF3,
  new_cureC1,
  new_cureC2,
  new_cureC3,
  new_cureC4)

# Age dependent transition vectors. Currently uniform across age groups
f0f1_vec <- rep(f0f1,groups)
f3c1_vec <- rep(f3c1,groups)

# Set up time ####

simu.time <- seq(0, 36, by=1) 

# Initial conditions ####

initS <- as.numeric(init_condits[1,]) # S row

initF0 <- as.numeric(init_condits[2,]) # F0 row
initF1 <- as.numeric(init_condits[3,])
initF2 <- as.numeric(init_condits[4,])
initF3 <- as.numeric(init_condits[5,])

initC1 <- as.numeric(init_condits[6,])
initC2 <- as.numeric(init_condits[7,])
initC3 <- as.numeric(init_condits[8,])
initC4 <- as.numeric(init_condits[9,])

initHCCA <- as.numeric(init_condits[10,])
initHCCB <- as.numeric(init_condits[11,])
initHCCC <- as.numeric(init_condits[12,])
initHCCD <- as.numeric(init_condits[13,])

initF0cure <- as.numeric(init_condits[14,])
initF1cure <- as.numeric(init_condits[15,])
initF2cure <- as.numeric(init_condits[16,])
initF3cure <- as.numeric(init_condits[17,])

initC1cure <- as.numeric(init_condits[18,])
initC2cure <- as.numeric(init_condits[19,])
initC3cure <- as.numeric(init_condits[20,])
initC4cure <- as.numeric(init_condits[21,])

#set up initial death
initD <- as.numeric(init_condits[22,])
initdthC14 <- as.numeric(init_condits[23,])
initdthHCC <- as.numeric(init_condits[24,])

initDHCC <- as.numeric(init_condits[25,])
initDC14 <- as.numeric(init_condits[26,])
inittransliv <- as.numeric(init_condits[27,])

initCInc <- rep(0,groups)
initCScr <- rep(0,groups)

# Set up initial values

init <- c(S=initS,F0=initF0,F1=initF1,F2=initF2,F3=initF3,C1=initC1,C2=initC2,C3=initC3,C4=initC4,
          HCCA=initHCCA,HCCB=initHCCB,HCCC=initHCCC,HCCD=initHCCD,
          F0cure=initF0cure,F1cure=initF1cure,F2cure=initF2cure,F3cure=initF3cure,
          C1cure=initC1cure,C2cure=initC2cure,C3cure=initC3cure,C4cure=initC4cure,
          D=initD,dthC14=initdthC14,dthHCC=initdthHCC,transliv=inittransliv, CInc=initCInc, CScr=initCScr)



# Birth rate function ####

# Birth rate multiplier (brm) with uncertainty
brm_mean <- 1.1
brm_sd <- 0.025

brm_upper <- brm_mean + 1.96*brm_sd
brm_lower <- brm_mean - 1.96*brm_sd

# Define model ####
hepC.mod<- function(t, y, parameters)
{
  with(as.list(c(y, parameters)),
       {
         
         S <- y[Sindex]
         F0 <- y[F0index]
         F1 <- y[F1index]
         F2 <- y[F2index]
         F3 <- y[F3index]
         C1 <- y[C1index]
         C2 <- y[C2index]
         C3 <- y[C3index]
         C4 <- y[C4index]
         
         HCCA <- y[HCCAindex]
         HCCB <- y[HCCBindex]
         HCCC <- y[HCCCindex]
         HCCD <- y[HCCDindex]
         
         F0cure <- y[F0cureindex]
         F1cure <- y[F1cureindex]
         F2cure <- y[F2cureindex]
         F3cure <- y[F3cureindex]
         C1cure <- y[C1cureindex]
         C2cure <- y[C2cureindex]
         C3cure <- y[C3cureindex]
         C4cure <- y[C4cureindex]
         
         D <- y[Dindex]
         dthC14 <- y[dthC14index]
         dthHCC <- y[dthHCCindex]
         
         transliv <- y[translivindex]
         
         alive <- y[aliveindex]
         infect <- y[infectindex]
         totalHCC <- y[totalHCCindex]
         totalHCV <- y[totalHCVindex]
         
         newdeath <- y[newdeathindex]
         incHCC <- y[incHCCindex]
         
         dCInc <- y[CIncindex]
         
         dCScr <- y[CScrindex]

         
         # Screening

         if((t>scr_start)&(t<scr_start+scr_dur)){scr<-scr_group}else{scr<-groupBASE}

         alive <- S+F0+F1+F2+F3+C1+C2+C3+C4+HCCA+HCCB+HCCC+HCCD+F0cure+F1cure+F2cure+F3cure+C1cure+C2cure+C3cure+C4cure
         pop <- sum(alive)
         
         infect <- (F0+F1+F2+F3+C1+C2+C3+C4+HCCA+HCCB+HCCC+HCCD)
         lambda <- beta%*%infect/pop
         totalHCC <- HCCA+HCCB+HCCC+HCCD
         totalHCV  <- F0+F1+F2+F3+C1+C2+C3+C4
         prevalence <- 100*(infect/pop)
         
         
         flowin <- birth.func(t)*pop
         birth <- c(flowin, rep(0,groups-1))
         
         cureF0 <- ifelse(t<=trt_start, std_cureF0, new_cureF0) # efficacy of old and new treatment
         cureF1 <- ifelse(t<=trt_start, std_cureF1, new_cureF1)
         cureF2 <- ifelse(t<=trt_start, std_cureF2, new_cureF2)
         cureF3 <- ifelse(t<=trt_start, std_cureF3, new_cureF3)
         cureC1 <- ifelse(t<=trt_start, std_cureC1, new_cureC1)
         cureC2 <- ifelse(t<=trt_start, 0, new_cureC2)
         cureC3 <- ifelse(t<=trt_start, 0, new_cureC3)
         cureC4 <- ifelse(t<=trt_start, 0, new_cureC4)
         
         natdeath <- c(mortality_func1(t),mortality_func2(t),mortality_func3(t),mortality_func4(t),mortality_func5(t),mortality_func6(t),mortality_func7(t),mortality_func8(t),mortality_func9(t),mortality_func10(t),
                       mortality_func11(t),mortality_func12(t),mortality_func13(t),mortality_func14(t),mortality_func15(t),mortality_func16(t),mortality_func17(t),mortality_func18(t),mortality_func19(t),mortality_func20(t),mortality_func21(t))

         # ODEs
         dS <- birth -lambda*S -natdeath*S + aging %*% S
         dF0 <- ifelse(F0>=0,-f0f1_vec*F0 + lambda*S - scr*pF0scr*sens*cureF0*F0 -natdeath*F0 + aging %*% F0,lambda*S + aging %*% F0)
         dF1 <- ifelse(F1>=0,f0f1_vec*F0 -f1f2*F1 - scr*pF1scr*sens*cureF1*F1 -natdeath*F1 + aging %*% F1,rep(0,groups))
         dF2 <- ifelse(F2>=0,f1f2*F1 -f2f3*F2 - scr*pF2scr*sens*cureF2*F2  -natdeath*F2 + aging %*% F2,rep(0,groups))
         dF3 <- ifelse(F3>=0,f2f3*F2 -f3c1_vec*F3 - scr*pF3scr*sens*cureF3*F3 -natdeath*F3 + aging %*% F3,rep(0,groups))
         dC1 <- ifelse(C1>=0,f3c1_vec*F3 -dthc1*C1 -c1c2*C1 -scr*pC1scr*sens*cureC1*C1 -(c1bA + c1bB + c1bC + c1bD)*C1 -natdeath*C1 + aging %*% C1,rep(0,groups))
         dC2 <- ifelse(C2>=0,c1c2*C1 -dthc2*C2 -c2c3*C2 - scr*pC2scr*sens*cureC2*C2 -(c2bA + c2bB + c2bC + c2bD)*C2 -natdeath*C2 + aging %*% C2,rep(0,groups))
         dC3 <- ifelse(C3>=0,c2c3*C2 -dthc3*C3 -c3c4*C3 - scr*pC3scr*sens*cureC3*C3 -c3bD*C3 -natdeath*C3 + aging %*% C3,rep(0,groups))
         dC4 <- ifelse(C4>=0,c3c4*C3 -dthc4*C4 - scr*pC4scr*sens*cureC4*C4 - c4bD*C4 - tranc4*C4 -natdeath*C4 + aging %*% C4,rep(0,groups))
         
         dHCCA <- c1bA*(C1+C1cure) + c2bA*(C2+C2cure) -dthbA*HCCA -tranbA*HCCA -natdeath*HCCA + aging %*% HCCA
         dHCCB <- c1bB*(C1+C1cure) + c2bB*(C2+C2cure) -dthbB*HCCB -tranbB*HCCB -natdeath*HCCB + aging %*% HCCB
         dHCCC <- c1bC*(C1+C1cure) + c2bC*(C2+C2cure) -dthbC*HCCC -natdeath*HCCC + aging %*% HCCC
         dHCCD <- c1bD*(C1+C1cure) + c2bD*(C2+C2cure) + c3bD*(C3+C3cure) + c4bD*(C4+C4cure) -dthbD*HCCD -natdeath*HCCD + aging %*% HCCD
         
         dF0cure <- ifelse(F0cure>=0, scr*pF0scr*sens*cureF0*F0 - natdeath*F0cure + aging %*% F0cure, rep(0,groups))
         dF1cure <- ifelse(F1cure>=0, scr*pF1scr*sens*cureF1*F1 - natdeath*F1cure + aging %*% F1cure, rep(0,groups))
         dF2cure <- ifelse(F2cure>=0, scr*pF2scr*sens*cureF2*F2 - natdeath*F2cure + aging %*% F2cure, rep(0,groups))
         dF3cure <- ifelse(F3cure>=0, scr*pF3scr*sens*cureF3*F3 - natdeath*F3cure + aging %*% F3cure, rep(0,groups))
         dC1cure <- ifelse(C1cure>=0, scr*pC1scr*sens*cureC1*C1 - natdeath*C1cure -(c1bA + c1bB + c1bC + c1bD)*C1cure + aging %*% C1cure, rep(0,groups))
         dC2cure <- ifelse(C2cure>=0, scr*pC2scr*sens*cureC2*C2 - natdeath*C2cure -(c2bA + c2bB + c2bC + c2bD)*C2cure + aging %*% C2cure, rep(0,groups))
         dC3cure <- ifelse(C3cure>=0, scr*pC3scr*sens*cureC3*C3 - natdeath*C3cure - c3bD*C3cure + aging %*% C3cure, rep(0,groups))
         dC4cure <- ifelse(C4cure>=0, scr*pC4scr*sens*cureC4*C4 - natdeath*C4cure - c4bD*C4cure + aging %*% C4cure, rep(0,groups))
         
         dD <- dthc1*C1 + dthc2*C2 + dthc3*C3 + dthc4*C4 + dthbA*HCCA + dthbB*HCCB + dthbC*HCCC + dthbD*HCCD
         ddthC14 <- dthc1*C1 + dthc2*C2 + dthc3*C3 + dthc4*C4
         ddthHCC <- dthbA*HCCA + dthbB*HCCB + dthbC*HCCC + dthbD*HCCD
         
         dtransliv <- tranbA*HCCA + tranbB*HCCB + tranc4*C4
         dCInc <- lambda * S
         
         # Cumulative number of individuals screened/treated
         dCScr <- scr*sens*(pF0scr*cureF0*F0 + 
                              pF1scr*cureF1*F1 + 
                              pF2scr*cureF2*F2 + 
                              pF3scr*cureF3*F3 + 
                              pC1scr*cureC1*C1 +
                              pC2scr*cureC2*C2 + 
                              pC3scr*cureC3*C3 + 
                              pC4scr*cureC4*C4)
         
         newdeath <- ifelse(C1+C2+C3+totalHCC>0, dthc1*C1 + dthc2*C2 + dthc3*C3 + dthc4*C4 + dthbA*HCCA + dthbB*HCCB + dthbC*HCCC + dthbD*HCCD,0)
         incHCC <- c1bA*(C1+C1cure)+c2bA*(C2+C2cure)+c1bB*(C1+C1cure)+
           c2bB*(C2+C2cure)+c1bC*(C1+C1cure)+c2bC*(C2+C2cure)+c1bD*(C1+C1cure)+
           c2bD*(C2+C2cure)+c3bD*(C3+C3cure)+c4bD*(C4+C4cure)

         list(c(dS,dF0,dF1,dF2,dF3,dC1,dC2,dC3,dC4,
                dHCCA,dHCCB,dHCCC,dHCCD,
                dF0cure,dF1cure,dF2cure,dF3cure,dC1cure,dC2cure,dC3cure,dC4cure,
                dD,ddthC14,ddthHCC,dtransliv, dCInc, dCScr),
              incidenceHCC=incHCC,newdeath=newdeath,prevalence=prevalence,population=pop,total.infection=infect,
              totalHCC=totalHCC,totalHCV=totalHCV)
       }
  )
}

# Mutate data function ####

mutate_data <- function(ode_output){
  
  df1<-as_tibble(as.data.frame(ode_output)) %>% 
    mutate(
      
      # Totals across all age groups
      S=(S1+S2+S3+S4+S5+S6+S7+S8+S9+S10+S11+S12+S3+S14+S15+S16+S17+S18+S19+S20+S21),
      
      # Total Screened
      CScr=(CScr1+CScr2+CScr3+CScr4+CScr5+CScr6+CScr7+CScr8+CScr9+CScr10+CScr11+CScr12+CScr3+CScr14+CScr15+CScr16+CScr17+CScr18+CScr19+CScr20+CScr21),
      
      # Total Cases
      CInc=(CInc1+CInc2+CInc3+CInc4+CInc5+CInc6+CInc7+CInc8+CInc9+CInc10+CInc11+CInc12+CInc3+CInc14+CInc15+CInc16+CInc17+CInc18+CInc19+CInc20+CInc21),
      
      # Cured
      F0cure=(F0cure1+F0cure2+F0cure3+F0cure4+F0cure5+F0cure6+F0cure7+F0cure8+F0cure9+F0cure10+F0cure11+
                F0cure12+F0cure13+F0cure14+F0cure15+F0cure16+F0cure17+F0cure18+F0cure19+F0cure20+F0cure21),
      F1cure=(F1cure1+F1cure2+F1cure3+F1cure4+F1cure5+F1cure6+F1cure7+F1cure8+F1cure9+F1cure10+F1cure11+
                F1cure12+F1cure13+F1cure14+F1cure15+F1cure16+F1cure17+F1cure18+F1cure19+F1cure20+F1cure21),
      F2cure=(F2cure1+F2cure2+F2cure3+F2cure4+F2cure3+F2cure6+F2cure7+F2cure8+F2cure9+F2cure10+F2cure11+
                F2cure12+F2cure13+F2cure14+F2cure15+F2cure16+F2cure17+F2cure18+F2cure19+F2cure20+F2cure21),
      F3cure=(F3cure1+F3cure2+F3cure3+F3cure4+F3cure5+F3cure6+F3cure7+F3cure8+F3cure9+F3cure10+F3cure11+
                F3cure12+F3cure13+F3cure14+F3cure15+F3cure16+F3cure17+F3cure18+F3cure19+F3cure20+F3cure21),
      C1cure=(C1cure1+C1cure2+C1cure3+C1cure4+C1cure5+C1cure6+C1cure7+C1cure8+C1cure9+C1cure10+C1cure11+
                C1cure12+C1cure13+C1cure14+C1cure15+C1cure16+C1cure17+C1cure18+C1cure19+C1cure20+C1cure21),
      C2cure=(C2cure1+C2cure2+C2cure3+C2cure4+C2cure5+C2cure6+C2cure7+C2cure8+C2cure9+C2cure10+C2cure11+
                C2cure12+C2cure13+C2cure14+C2cure15+C2cure16+C2cure17+C2cure18+C2cure19+C2cure20+C2cure21),
      C3cure=(C3cure1+C3cure2+C3cure3+C3cure4+C3cure5+C3cure6+C3cure7+C3cure8+C3cure9+C3cure10+C3cure11+
                C3cure12+C3cure13+C3cure14+C3cure15+C3cure16+C3cure17+C3cure18+C3cure19+C3cure20+C3cure21),
      C4cure=(C4cure1+C4cure2+C4cure3+C4cure4+C4cure5+C4cure6+C4cure7+C4cure8+C4cure9+C4cure10+C4cure11+
                C4cure12+C4cure13+C4cure14+C4cure15+C4cure16+C4cure17+C4cure18+C4cure19+C4cure20+C4cure21),
      
      cured = F0cure+F1cure+F2cure+F3cure+C1cure+C2cure+C3cure+C4cure,

      # Fibrosis, cirrhosis and HCC totals (sum of all age groups)
      F0=(F01+F02+F03+F04+F05+F06+F07+F08+F09+F010+F011+F012+F013+F014+F015+F016+F017+F018+F019+F020+F021),
      F1=(F11+F12+F13+F14+F15+F16+F17+F18+F19+F110+F111+F112+F113+F114+F115+F116+F117+F118+F119+F120+F121),
      F2=(F21+F22+F23+F24+F23+F26+F27+F28+F29+F210+F211+F212+F213+F214+F215+F216+F217+F218+F219+F220+F221),
      F3=(F31+F32+F33+F34+F35+F36+F37+F38+F39+F310+F311+F312+F313+F314+F315+F316+F317+F318+F319+F320+F321),
      C1=(C11+C12+C13+C14+C15+C16+C17+C18+C19+C110+C111+C112+C113+C114+C115+C116+C117+C118+C119+C120+C121),
      C2=(C21+C22+C23+C24+C25+C26+C27+C28+C29+C210+C211+C212+C213+C214+C215+C216+C217+C218+C219+C220+C221),
      C3=(C31+C32+C33+C34+C35+C36+C37+C38+C39+C310+C311+C312+C313+C314+C315+C316+C317+C318+C319+C320+C321),
      C4=(C41+C42+C43+C44+C45+C46+C47+C48+C49+C410+C411+C412+C413+C414+C415+C416+C417+C418+C419+C420+C421),
      HCCA=(HCCA1+HCCA2+HCCA3+HCCA4+HCCA5+HCCA6+HCCA7+HCCA8+HCCA9+HCCA10+HCCA11+HCCA12+HCCA13+HCCA14+HCCA15+HCCA16+HCCA17+HCCA18+HCCA19+HCCA20+HCCA21),
      HCCB=(HCCB1+HCCB2+HCCB3+HCCB4+HCCB5+HCCB6+HCCB7+HCCB8+HCCB9+HCCB10+HCCB11+HCCB12+HCCB13+HCCB14+HCCB15+HCCB16+HCCB17+HCCB18+HCCB19+HCCB20+HCCB21),
      HCCC=(HCCC1+HCCC2+HCCC3+HCCC4+HCCC5+HCCC6+HCCC7+HCCC8+HCCC9+HCCC10+HCCC11+HCCC12+HCCC13+HCCC14+HCCC15+HCCC16+HCCC17+HCCC18+HCCC19+HCCC20+HCCC21),
      HCCD=(HCCD1+HCCD2+HCCD3+HCCD4+HCCD5+HCCD6+HCCD7+HCCD8+HCCD9+HCCD10+HCCD11+HCCD12+HCCD13+HCCD14+HCCD15+HCCD16+HCCD17+HCCD18+HCCD19+HCCD20+HCCD21),
      
      # Incidence and mortality
      CInc=(CInc1+CInc2+CInc3+CInc4+CInc5+CInc6+CInc7+CInc8+CInc9+CInc10+CInc11+CInc12+CInc13+CInc14+CInc15+CInc16+CInc17+CInc18+CInc19+CInc20+CInc21),
      D=D1+D2+D3+D4+D5+D6+D7+D8+D9+D10+D11+D12+D13+D14+D15+D16+D17+D18+D19+D20+D21,
      # Yearly mortality
      Deaths = c(0, diff(D)),
      # Yearly incidence
      Inc = c(0, diff(CInc)),
      
      Inc5 = c(0, diff(CInc5)),
      Inc6 = c(0, diff(CInc6)),
      Inc7 = c(0, diff(CInc7)),
      Inc8 = c(0, diff(CInc8)),
      Inc9 = c(0, diff(CInc9)),
      Inc10 = c(0, diff(CInc10)),
      Inc11 = c(0, diff(CInc11)),
      Inc12 = c(0, diff(CInc12)),
      Inc13 = c(0, diff(CInc13)),
      Inc14 = c(0, diff(CInc14)),
      
      # Total population of entire system (should fit Thai population data)
      total = population,
      
      # Total infections, HCC and HCV (sum across all age groups)
      infect = (F0+F1+F2+F3+C1+C2+C3+C4+HCCA+HCCB+HCCC+HCCD),
      totalHCC = (HCCA+HCCB+HCCC+HCCD),
      totalHCV  = (F0+F1+F2+F3+C1+C2+C3+C4),
      
      # Prevalance of the above (%)
      infectprev = 100*(infect/total),
      HCCprev = 100*(totalHCC/total),
      HCVprev = 100*(totalHCV/total),
      
      
      # Total infections for each individual age group
      infect1 = (F01+F11+F21+F31+C11+C21+C31+C41+HCCA1+HCCB1+HCCC1+HCCD1),
      infect2 = (F02+F12+F22+F32+C12+C22+C32+C42+HCCA2+HCCB2+HCCC2+HCCD2),
      infect3 = (F03+F13+F23+F33+C13+C23+C33+C43+HCCA3+HCCB3+HCCC3+HCCD3),
      infect4 = (F04+F14+F24+F34+C14+C24+C34+C44+HCCA4+HCCB4+HCCC4+HCCD4),
      infect5 = (F05+F15+F25+F35+C15+C25+C35+C45+HCCA5+HCCB5+HCCC5+HCCD5),
      infect6 = (F06+F16+F26+F36+C16+C26+C36+C46+HCCA6+HCCB6+HCCC6+HCCD6),
      infect7 = (F07+F17+F27+F37+C17+C27+C37+C47+HCCA7+HCCB7+HCCC7+HCCD7),
      infect8 = (F08+F18+F28+F38+C18+C28+C38+C48+HCCA8+HCCB8+HCCC8+HCCD8),
      infect9 = (F09+F19+F29+F39+C19+C29+C39+C49+HCCA9+HCCB9+HCCC9+HCCD9),
      infect10 = (F010+F110+F210+F310+C110+C210+C310+C410+HCCA10+HCCB10+HCCC10+HCCD10),
      infect11 = (F011+F111+F211+F311+C111+C211+C311+C411+HCCA11+HCCB11+HCCC11+HCCD11),
      infect12 = (F012+F112+F212+F312+C112+C212+C312+C412+HCCA12+HCCB12+HCCC12+HCCD12),
      infect13 = (F013+F113+F213+F313+C113+C213+C313+C413+HCCA13+HCCB13+HCCC13+HCCD13),
      infect14 = (F014+F114+F214+F314+C114+C214+C314+C414+HCCA14+HCCB14+HCCC14+HCCD14),
      infect15 = (F015+F115+F215+F315+C115+C215+C315+C415+HCCA15+HCCB15+HCCC15+HCCD15),
      infect16 = (F016+F116+F216+F316+C116+C216+C316+C416+HCCA16+HCCB16+HCCC16+HCCD16),
      infect17 = (F017+F117+F217+F317+C117+C217+C317+C417+HCCA17+HCCB17+HCCC17+HCCD17),
      infect18 = (F018+F118+F218+F318+C118+C218+C318+C418+HCCA18+HCCB18+HCCC18+HCCD18),
      infect19 = (F019+F119+F219+F319+C119+C219+C319+C419+HCCA19+HCCB19+HCCC19+HCCD19),
      infect20 = (F020+F120+F220+F320+C120+C220+C320+C420+HCCA20+HCCB20+HCCC20+HCCD20),
      infect21 = (F021+F121+F221+F321+C121+C221+C321+C421+HCCA21+HCCB21+HCCC21+HCCD21),
      
      # Total proportion in each age group (sum of all compartments per age group)
      group1 = (S1+F01+F11+F21+F31+C11+C21+C31+C41+HCCA1+HCCB1+HCCC1+HCCD1+F0cure1+F1cure1+F2cure1+F3cure1+C1cure1+C2cure1+C3cure1+C4cure1), # which compartments contribute? not the cumulative ones but do I need std and new cured or will this result in overcounting?
      group2 = (S2+F02+F12+F22+F32+C12+C22+C32+C42+HCCA2+HCCB2+HCCC2+HCCD2+F0cure2+F1cure2+F2cure2+F3cure2+C1cure2+C2cure2+C3cure2+C4cure2),
      group3 = (S3+F03+F13+F23+F33+C13+C23+C33+C43+HCCA3+HCCB3+HCCC3+HCCD3+F0cure3+F1cure3+F2cure3+F3cure3+C1cure3+C2cure3+C3cure3+C4cure3),
      group4 = (S4+F04+F14+F24+F34+C14+C24+C34+C44+HCCA4+HCCB4+HCCC4+HCCD4+F0cure4+F1cure4+F2cure4+F3cure4+C1cure4+C2cure4+C3cure4+C4cure4),
      group5 = (S5+F05+F15+F25+F35+C15+C25+C35+C45+HCCA5+HCCB5+HCCC5+HCCD5+F0cure5+F1cure5+F2cure5+F3cure5+C1cure5+C2cure5+C3cure5+C4cure5),
      group6 = (S6+F06+F16+F26+F36+C16+C26+C36+C46+HCCA6+HCCB6+HCCC6+HCCD6+F0cure6+F1cure6+F2cure6+F3cure6+C1cure6+C2cure6+C3cure6+C4cure6),
      group7 = (S7+F07+F17+F27+F37+C17+C27+C37+C47+HCCA7+HCCB7+HCCC7+HCCD7+F0cure7+F1cure7+F2cure7+F3cure7+C1cure7+C2cure7+C3cure7+C4cure7),
      group8 = (S8+F08+F18+F28+F38+C18+C28+C38+C48+HCCA8+HCCB8+HCCC8+HCCD8+F0cure8+F1cure8+F2cure8+F3cure8+C1cure8+C2cure8+C3cure8+C4cure8),
      group9 = (S9+F09+F19+F29+F39+C19+C29+C39+C49+HCCA9+HCCB9+HCCC9+HCCD9+F0cure9+F1cure9+F2cure9+F3cure9+C1cure9+C2cure9+C3cure9+C4cure9),
      group10 = (S10+F010+F110+F210+F310+C110+C210+C310+C410+HCCA10+HCCB10+HCCC10+HCCD10+F0cure10+F1cure10+F2cure10+F3cure10+C1cure10+C2cure10+C3cure10+C4cure10), 
      group11 = (S11+F011+F111+F211+F311+C111+C211+C311+C411+HCCA11+HCCB11+HCCC11+HCCD11+F0cure11+F1cure11+F2cure11+F3cure11+C1cure11+C2cure11+C3cure11+C4cure11),
      group12 = (S12+F012+F112+F212+F312+C112+C212+C312+C412+HCCA12+HCCB12+HCCC12+HCCD12+F0cure12+F1cure12+F2cure12+F3cure12+C1cure12+C2cure12+C3cure12+C4cure12),
      group13 = (S13+F013+F113+F213+F313+C113+C213+C313+C413+HCCA13+HCCB13+HCCC13+HCCD13+F0cure13+F1cure13+F2cure13+F3cure13+C1cure13+C2cure13+C3cure13+C4cure13),
      group14 = (S14+F014+F114+F214+F314+C114+C214+C314+C414+HCCA14+HCCB14+HCCC14+HCCD14+F0cure14+F1cure14+F2cure14+F3cure14+C1cure14+C2cure14+C3cure14+C4cure14),
      group15 = (S15+F015+F115+F215+F315+C115+C215+C315+C415+HCCA15+HCCB15+HCCC15+HCCD15+F0cure15+F1cure15+F2cure15+F3cure15+C1cure15+C2cure15+C3cure15+C4cure15),
      group16 = (S16+F016+F116+F216+F316+C116+C216+C316+C416+HCCA16+HCCB16+HCCC16+HCCD16+F0cure16+F1cure16+F2cure16+F3cure16+C1cure16+C2cure16+C3cure16+C4cure16),
      group17 = (S17+F017+F117+F217+F317+C117+C217+C317+C417+HCCA17+HCCB17+HCCC17+HCCD17+F0cure17+F1cure17+F2cure17+F3cure17+C1cure17+C2cure17+C3cure17+C4cure17),
      group18 = (S18+F018+F118+F218+F318+C118+C218+C318+C418+HCCA18+HCCB18+HCCC18+HCCD18+F0cure18+F1cure18+F2cure18+F3cure18+C1cure18+C2cure18+C3cure18+C4cure18),
      group19 = (S19+F019+F119+F219+F319+C119+C219+C319+C419+HCCA19+HCCB19+HCCC19+HCCD19+F0cure19+F1cure19+F2cure19+F3cure19+C1cure19+C2cure19+C3cure19+C4cure19),
      group20 = (S20+F020+F120+F220+F320+C120+C220+C320+C420+HCCA20+HCCB20+HCCC20+HCCD20+F0cure20+F1cure20+F2cure20+F3cure20+C1cure20+C2cure20+C3cure20+C4cure20),
      group21 = (S21+F021+F121+F221+F321+C121+C221+C321+C421+HCCA21+HCCB21+HCCC21+HCCD21+F0cure21+F1cure21+F2cure21+F3cure21+C1cure21+C2cure21+C3cure21+C4cure21),
      
      # Prevalence by age group (%)
      prev1 = (infect1 * 100) / (group1),
      prev2 = (infect2 * 100) / (group2),
      prev3 = (infect3 * 100) / (group3),
      prev4 = (infect4 * 100) / (group4),
      prev5 = (infect5 * 100) / (group5),
      prev6 = (infect6 * 100) / (group6),
      prev7 = (infect7 * 100) / (group7),
      prev8 = (infect8 * 100) / (group8),
      prev9 = (infect9 * 100) / (group9),
      prev10 = (infect10 * 100) / (group10),
      prev11 = (infect11 * 100) / (group11),
      prev12 = (infect12 * 100) / (group12),
      prev13 = (infect13 * 100) / (group13),
      prev14 = (infect14 * 100) / (group14),
      prev15 = (infect15 * 100) / (group15),
      prev16 = (infect16 * 100) / (group16),
      prev17 = (infect17 * 100) / (group17),
      prev18 = (infect18 * 100) / (group18),
      prev19 = (infect19 * 100) / (group19),
      prev20 = (infect20 * 100) / (group20),
      prev21 = (infect21 * 100) / (group21),
      
      # Aggregated prevalence by age group from puclished data
      prev1_2 = 100*(infect1+infect2)/(group1+group2), # 0-9
      prev3_4 = 100*(infect3+infect4)/(group3+group4), # 10-19
      prev5_6 = 100*(infect5+infect6)/(group5+group6), # 20-29
      prev7_8 = 100*(infect7+infect8)/(group7+group8), # 30-39
      prev9_10 = 100*(infect9+infect10)/(group9+group10), #40-49
      prev11_21 = 100*(infect11+infect12+infect13+infect14+infect15+infect16+infect17+infect18+infect19+infect20+infect21)/(group11+group12+group13+group14+group15+group16+group17+group18+group19+group20+group21) # Over 50
      
    ) %>% 
    pivot_longer(names_to = "variable", cols = !1) 
  return(df1)
}


# Collect results function ####

Collect_Results <- function(dataframe){
  
  inc_2030 <- (dataframe %>% filter(Year %in% c("2030")) %>% filter(variable %in% c("Inc")))$value
  inc_2015 <- as.numeric((dataframe %>% filter(variable %in% c("Inc")) %>% filter(Year %in% c("2015")))[,3])
  inc_target <- 0.1*inc_2015
  inc_target_diff <- inc_2030 - inc_target
  
  mort_2030 <- (dataframe %>% filter(Year %in% c("2030")) %>% filter(variable %in% c("Deaths")))$value
  mort_2015 <- as.numeric((dataframe %>% filter(variable %in% c("Deaths")) %>% filter(Year %in% c("2015")))[,3])
  mort_target <- 0.35*mort_2015
  mort_target_diff <- mort_2030 - mort_target
  
  CScr_2023 <- (dataframe %>% filter(Year %in% c("2023")) %>% filter(variable %in% c("CScr")))$value
  CScr_2030 <- (dataframe %>% filter(Year %in% c("2030")) %>% filter(variable %in% c("CScr")))$value
  total_screened <- CScr_2030 - CScr_2023
  
  CInc_2023 <- (dataframe %>% filter(Year %in% c("2023")) %>% filter(variable %in% c("CInc")))$value
  CInc_2030 <- (dataframe %>% filter(Year %in% c("2030")) %>% filter(variable %in% c("CInc")))$value
  total_cases <- CInc_2030 - CInc_2023
  
  CD_2023 <- (dataframe %>% filter(Year %in% c("2023")) %>% filter(variable %in% c("D")))$value
  CD_2030 <- (dataframe %>% filter(Year %in% c("2030")) %>% filter(variable %in% c("D")))$value
  total_deaths <- CD_2030 - CD_2023
  
  inc_elim_year <- 0
  inc_only <- dataframe %>% filter(variable %in% c("Inc"))
  for(i in 2:nrow(inc_only)){
    if(inc_only[i,]$value <= inc_target){
      inc_elim_year <- inc_only[i,]$Year
      break
    }
  }
  
  if(inc_elim_year==2040){
    inc_elim_year <- "Beyond simulation" }
  
  mort_elim_year <- 0
  mort_only <- dataframe %>% filter(variable %in% c("Deaths"))
  for(i in 2:nrow(mort_only)){
    if(mort_only[i,]$value <= mort_target){
      mort_elim_year <- mort_only[i,]$Year
      break
    }
  }
  
  if(mort_elim_year==0){
    mort_elim_year <- "Beyond simulation"
  }
  
  
  
  if(inc_target_diff<=0){
    Inc_Elim = "Yes"}
  else{Inc_Elim = "No"}
  
  if(mort_target_diff<=0){
    Mort_Elim = "Yes"}
  else{Mort_Elim = "No"}
  
  
  return(c("Incidence 2030" = as.numeric(inc_2030), "Incidence Difference to Target" = as.numeric(inc_target_diff), "Deaths 2030" = as.numeric(mort_2030), "Deaths Difference to Target" = as.numeric(mort_target_diff), "Total Cases" = as.numeric(total_cases), "Cases Averted" = NA, "Total Deaths" = total_deaths, "Deaths Averted" = NA, "Total Screened" =  as.numeric(total_screened), "Extra Screened" = NA, "Year Incidence Elimination Reached" = inc_elim_year, "Year Deaths Elimination Reached" = mort_elim_year))
  
}
# Running baseline mortality scenario ####

mort_scenario <- scenario_base
mortality_base <- mortality_base_empty

for(i in 19:37){
  mortality_base[i,2:22] <- mortality_base[i-1,2:22]*mort_scenario
}

mortality_approx <- mortality_base

mortality_func1 <- approxfun(mortality_approx$t,mortality_approx$rate1,method="linear")
mortality_func2 <- approxfun(mortality_approx$t,mortality_approx$rate2,method="linear")
mortality_func3 <- approxfun(mortality_approx$t,mortality_approx$rate3,method="linear")
mortality_func4 <- approxfun(mortality_approx$t,mortality_approx$rate4,method="linear")
mortality_func5 <- approxfun(mortality_approx$t,mortality_approx$rate5,method="linear")
mortality_func6 <- approxfun(mortality_approx$t,mortality_approx$rate6,method="linear")
mortality_func7 <- approxfun(mortality_approx$t,mortality_approx$rate7,method="linear")
mortality_func8 <- approxfun(mortality_approx$t,mortality_approx$rate8,method="linear")
mortality_func9 <- approxfun(mortality_approx$t,mortality_approx$rate9,method="linear")
mortality_func10 <- approxfun(mortality_approx$t,mortality_approx$rate10,method="linear")
mortality_func11 <- approxfun(mortality_approx$t,mortality_approx$rate11,method="linear")
mortality_func12 <- approxfun(mortality_approx$t,mortality_approx$rate12,method="linear")
mortality_func13 <- approxfun(mortality_approx$t,mortality_approx$rate13,method="linear")
mortality_func14 <- approxfun(mortality_approx$t,mortality_approx$rate14,method="linear")
mortality_func15 <- approxfun(mortality_approx$t,mortality_approx$rate15,method="linear")
mortality_func16 <- approxfun(mortality_approx$t,mortality_approx$rate16,method="linear")
mortality_func17 <- approxfun(mortality_approx$t,mortality_approx$rate17,method="linear")
mortality_func18 <- approxfun(mortality_approx$t,mortality_approx$rate18,method="linear")
mortality_func19 <- approxfun(mortality_approx$t,mortality_approx$rate19,method="linear")
mortality_func20 <- approxfun(mortality_approx$t,mortality_approx$rate20,method="linear")
mortality_func21 <- approxfun(mortality_approx$t,mortality_approx$rate21,method="linear")


# Baseline screening strategy - whole population at fitted coverage
scr_group <- groupBASE
birth.approx_mean <- birth.approx
birthrate_multiplier <- brm_mean
birth.approx_mean$birth <- birthrate_multiplier * birth.approx_mean$birth
birth.func <- approxfun(birth.approx_mean$t,birth.approx_mean$birth,method="linear")

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_base_scr_base <- mutate_data(out)
results_base_base <- Collect_Results(df1_base_scr_base)

inc_2015 <- as.numeric((df1_base_scr_base %>% filter(variable %in% c("Inc")) %>% filter(Year %in% c("2015")))[,3])
inc_target <- 0.1*inc_2015

mort_2015 <- as.numeric((df1_base_scr_base %>% filter(variable %in% c("Deaths")) %>% filter(Year %in% c("2015")))[,3])
mort_target <- 0.35*mort_2015

# Lower bound for birth rate multiplier
birth.approx_lower <- birth.approx
birthrate_multiplier <- brm_lower
birth.approx_lower$birth <- birthrate_multiplier * birth.approx_lower$birth
birth.func <- approxfun(birth.approx_lower$t,birth.approx_lower$birth,method="linear")

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_base_scr_base_brm_lower <- mutate_data(out)

# Upper bound for birth rate multiplier
birth.approx_upper <- birth.approx
birthrate_multiplier <- brm_upper
birth.approx_upper$birth <- birthrate_multiplier * birth.approx_upper$birth
birth.func <- approxfun(birth.approx_upper$t,birth.approx_upper$birth,method="linear")

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_base_scr_base_brm_upper <- mutate_data(out)

# Lower bound for screening
# Return birth to mean for baseline
birth.approx_mean <- birth.approx
birthrate_multiplier <- brm_mean
birth.approx_mean$birth <- birthrate_multiplier * birth.approx_mean$birth
birth.func <- approxfun(birth.approx_mean$t,birth.approx_mean$birth,method="linear")

groupBASE <- groupBASE_lower # Lower level for baseline
scr_group <- groupBASE

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_base_scr_base_lower <- mutate_data(out)

# Upper bound for screening
groupBASE <- groupBASE_upper # Upper level for baseline
scr_group <- groupBASE

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_base_scr_base_upper <- mutate_data(out)


# Screening strategy A 
# Back to mean baseline
groupBASE <- rep(cov_mean,groups)
scr_group <- groupA

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_base_scr_A <- mutate_data(out)
results_base_A <- Collect_Results(df1_base_scr_A)

# Screening strategy B 

scr_group <- groupB

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_base_scr_B <- mutate_data(out)
results_base_B <- Collect_Results(df1_base_scr_B)

# Screening strategy C 

scr_group <- groupC

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_base_scr_C <- mutate_data(out)
results_base_C <- Collect_Results(df1_base_scr_C)

# Screening strategy D 

scr_group <- groupD

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_base_scr_D <- mutate_data(out)
results_base_D <- Collect_Results(df1_base_scr_D)

# Screening strategy E 

scr_group <- groupE

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_base_scr_E <- mutate_data(out)
results_base_E <- Collect_Results(df1_base_scr_E)

# Screening strategy F 

scr_group <- groupF

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_base_scr_F <- mutate_data(out)
results_base_F <- Collect_Results(df1_base_scr_F)

# Running mortality scenario 1 ####

mort_scenario <- scenario_1
mortality_1 <- mortality_base_empty

for(i in 19:37){
  mortality_1[i,2:22] <- mortality_1[i-1,2:22]*mort_scenario
}

mortality_approx <- mortality_1

mortality_func1 <- approxfun(mortality_approx$t,mortality_approx$rate1,method="linear")
mortality_func2 <- approxfun(mortality_approx$t,mortality_approx$rate2,method="linear")
mortality_func3 <- approxfun(mortality_approx$t,mortality_approx$rate3,method="linear")
mortality_func4 <- approxfun(mortality_approx$t,mortality_approx$rate4,method="linear")
mortality_func5 <- approxfun(mortality_approx$t,mortality_approx$rate5,method="linear")
mortality_func6 <- approxfun(mortality_approx$t,mortality_approx$rate6,method="linear")
mortality_func7 <- approxfun(mortality_approx$t,mortality_approx$rate7,method="linear")
mortality_func8 <- approxfun(mortality_approx$t,mortality_approx$rate8,method="linear")
mortality_func9 <- approxfun(mortality_approx$t,mortality_approx$rate9,method="linear")
mortality_func10 <- approxfun(mortality_approx$t,mortality_approx$rate10,method="linear")
mortality_func11 <- approxfun(mortality_approx$t,mortality_approx$rate11,method="linear")
mortality_func12 <- approxfun(mortality_approx$t,mortality_approx$rate12,method="linear")
mortality_func13 <- approxfun(mortality_approx$t,mortality_approx$rate13,method="linear")
mortality_func14 <- approxfun(mortality_approx$t,mortality_approx$rate14,method="linear")
mortality_func15 <- approxfun(mortality_approx$t,mortality_approx$rate15,method="linear")
mortality_func16 <- approxfun(mortality_approx$t,mortality_approx$rate16,method="linear")
mortality_func17 <- approxfun(mortality_approx$t,mortality_approx$rate17,method="linear")
mortality_func18 <- approxfun(mortality_approx$t,mortality_approx$rate18,method="linear")
mortality_func19 <- approxfun(mortality_approx$t,mortality_approx$rate19,method="linear")
mortality_func20 <- approxfun(mortality_approx$t,mortality_approx$rate20,method="linear")
mortality_func21 <- approxfun(mortality_approx$t,mortality_approx$rate21,method="linear")

# Baseline screening strategy
scr_group <- groupBASE

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_1_scr_base <- mutate_data(out)
results_1_base <- Collect_Results(df1_mort_1_scr_base)

# Screening strategy A 

scr_group <- groupA

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_1_scr_A <- mutate_data(out)
results_1_A <- Collect_Results(df1_mort_1_scr_A)

# Screening strategy B 

scr_group <- groupB

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_1_scr_B <- mutate_data(out)
results_1_B <- Collect_Results(df1_mort_1_scr_B)

# Screening strategy C 

scr_group <- groupC

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_1_scr_C <- mutate_data(out)
results_1_C <- Collect_Results(df1_mort_1_scr_C)

# Screening strategy D 

scr_group <- groupD

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_1_scr_D <- mutate_data(out)
results_1_D <- Collect_Results(df1_mort_1_scr_D)

# Screening strategy E 

scr_group <- groupE

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_1_scr_E <- mutate_data(out)
results_1_E <- Collect_Results(df1_mort_1_scr_E)

# Screening strategy F 

scr_group <- groupF

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_1_scr_F <- mutate_data(out)
results_1_F <- Collect_Results(df1_mort_1_scr_F)

# Running mortality scenario 2 ####

mort_scenario <- scenario_2
mortality_2 <- mortality_base_empty

for(i in 19:37){
  mortality_2[i,2:22] <- mortality_2[i-1,2:22]*mort_scenario
}

mortality_approx <- mortality_2

mortality_func1 <- approxfun(mortality_approx$t,mortality_approx$rate1,method="linear")
mortality_func2 <- approxfun(mortality_approx$t,mortality_approx$rate2,method="linear")
mortality_func3 <- approxfun(mortality_approx$t,mortality_approx$rate3,method="linear")
mortality_func4 <- approxfun(mortality_approx$t,mortality_approx$rate4,method="linear")
mortality_func5 <- approxfun(mortality_approx$t,mortality_approx$rate5,method="linear")
mortality_func6 <- approxfun(mortality_approx$t,mortality_approx$rate6,method="linear")
mortality_func7 <- approxfun(mortality_approx$t,mortality_approx$rate7,method="linear")
mortality_func8 <- approxfun(mortality_approx$t,mortality_approx$rate8,method="linear")
mortality_func9 <- approxfun(mortality_approx$t,mortality_approx$rate9,method="linear")
mortality_func10 <- approxfun(mortality_approx$t,mortality_approx$rate10,method="linear")
mortality_func11 <- approxfun(mortality_approx$t,mortality_approx$rate11,method="linear")
mortality_func12 <- approxfun(mortality_approx$t,mortality_approx$rate12,method="linear")
mortality_func13 <- approxfun(mortality_approx$t,mortality_approx$rate13,method="linear")
mortality_func14 <- approxfun(mortality_approx$t,mortality_approx$rate14,method="linear")
mortality_func15 <- approxfun(mortality_approx$t,mortality_approx$rate15,method="linear")
mortality_func16 <- approxfun(mortality_approx$t,mortality_approx$rate16,method="linear")
mortality_func17 <- approxfun(mortality_approx$t,mortality_approx$rate17,method="linear")
mortality_func18 <- approxfun(mortality_approx$t,mortality_approx$rate18,method="linear")
mortality_func19 <- approxfun(mortality_approx$t,mortality_approx$rate19,method="linear")
mortality_func20 <- approxfun(mortality_approx$t,mortality_approx$rate20,method="linear")
mortality_func21 <- approxfun(mortality_approx$t,mortality_approx$rate21,method="linear")

# Baseline screening strategy
scr_group <- groupBASE

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_2_scr_base <- mutate_data(out)
results_2_base <- Collect_Results(df1_mort_2_scr_base)

# Screening strategy A 

scr_group <- groupA

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_2_scr_A <- mutate_data(out)
results_2_A <- Collect_Results(df1_mort_2_scr_A)

# Screening strategy B 

scr_group <- groupB

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_2_scr_B <- mutate_data(out)
results_2_B <- Collect_Results(df1_mort_2_scr_B)

# Screening strategy C 

scr_group <- groupC

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_2_scr_C <- mutate_data(out)
results_2_C <- Collect_Results(df1_mort_2_scr_C)

# Screening strategy D 

scr_group <- groupD

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_2_scr_D <- mutate_data(out)
results_2_D <- Collect_Results(df1_mort_2_scr_D)

# Screening strategy E 

scr_group <- groupE

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_2_scr_E <- mutate_data(out)
results_2_E <- Collect_Results(df1_mort_2_scr_E)

# Screening strategy F 

scr_group <- groupF

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_2_scr_F <- mutate_data(out)
results_2_F <- Collect_Results(df1_mort_2_scr_F)

# Running mortality scenario 3 ####

mort_scenario <- scenario_3
mortality_3 <- mortality_base_empty

for(i in 19:37){
  mortality_3[i,2:22] <- mortality_3[i-1,2:22]*mort_scenario
}

mortality_approx <- mortality_3

mortality_func1 <- approxfun(mortality_approx$t,mortality_approx$rate1,method="linear")
mortality_func2 <- approxfun(mortality_approx$t,mortality_approx$rate2,method="linear")
mortality_func3 <- approxfun(mortality_approx$t,mortality_approx$rate3,method="linear")
mortality_func4 <- approxfun(mortality_approx$t,mortality_approx$rate4,method="linear")
mortality_func5 <- approxfun(mortality_approx$t,mortality_approx$rate5,method="linear")
mortality_func6 <- approxfun(mortality_approx$t,mortality_approx$rate6,method="linear")
mortality_func7 <- approxfun(mortality_approx$t,mortality_approx$rate7,method="linear")
mortality_func8 <- approxfun(mortality_approx$t,mortality_approx$rate8,method="linear")
mortality_func9 <- approxfun(mortality_approx$t,mortality_approx$rate9,method="linear")
mortality_func10 <- approxfun(mortality_approx$t,mortality_approx$rate10,method="linear")
mortality_func11 <- approxfun(mortality_approx$t,mortality_approx$rate11,method="linear")
mortality_func12 <- approxfun(mortality_approx$t,mortality_approx$rate12,method="linear")
mortality_func13 <- approxfun(mortality_approx$t,mortality_approx$rate13,method="linear")
mortality_func14 <- approxfun(mortality_approx$t,mortality_approx$rate14,method="linear")
mortality_func15 <- approxfun(mortality_approx$t,mortality_approx$rate15,method="linear")
mortality_func16 <- approxfun(mortality_approx$t,mortality_approx$rate16,method="linear")
mortality_func17 <- approxfun(mortality_approx$t,mortality_approx$rate17,method="linear")
mortality_func18 <- approxfun(mortality_approx$t,mortality_approx$rate18,method="linear")
mortality_func19 <- approxfun(mortality_approx$t,mortality_approx$rate19,method="linear")
mortality_func20 <- approxfun(mortality_approx$t,mortality_approx$rate20,method="linear")
mortality_func21 <- approxfun(mortality_approx$t,mortality_approx$rate21,method="linear")

# Baseline screening strategy
scr_group <- groupBASE

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_3_scr_base <- mutate_data(out)
results_3_base <- Collect_Results(df1_mort_3_scr_base)

# Screening strategy A 

scr_group <- groupA

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_3_scr_A <- mutate_data(out)
results_3_A <- Collect_Results(df1_mort_3_scr_A)

# Screening strategy B 

scr_group <- groupB

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_3_scr_B <- mutate_data(out)
results_3_B <- Collect_Results(df1_mort_3_scr_B)

# Screening strategy C 

scr_group <- groupC

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_3_scr_C <- mutate_data(out)
results_3_C <- Collect_Results(df1_mort_3_scr_C)

# Screening strategy D 

scr_group <- groupD

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_3_scr_D <- mutate_data(out)
results_3_D <- Collect_Results(df1_mort_3_scr_D)

# Screening strategy E 

scr_group <- groupE

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_3_scr_E <- mutate_data(out)
results_3_E <- Collect_Results(df1_mort_3_scr_E)

# Screening strategy F 

scr_group <- groupF

tic("modelrun")
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
toc()

out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]

df1_mort_3_scr_F <- mutate_data(out)
results_3_F <- Collect_Results(df1_mort_3_scr_F)

# Collect and view results ####

Results_Table <- rbind(results_base_base, results_base_A, results_base_B, results_base_C, results_base_D, results_base_E, results_base_F,
                       results_1_base, results_1_A, results_1_B, results_1_C, results_1_D, results_1_E, results_1_F,
                       results_2_base, results_2_A, results_2_B, results_2_C, results_2_D, results_2_E, results_2_F,
                       results_3_base, results_3_A, results_3_B, results_3_C, results_3_D, results_3_E, results_3_F)

for(j in 1:nrow(Results_Table)){
  for(i in 1:5){
    Results_Table[j,i] <- round(as.numeric(Results_Table[j,i]),0)
    Results_Table[j,7] <- round(as.numeric(Results_Table[j,7]),0)
    Results_Table[j,9] <- round(as.numeric(Results_Table[j,9]),0)
  }
}

for(i in 1:nrow(Results_Table)){
  Results_Table[i,6] <- as.numeric(Results_Table[1,5]) - as.numeric(Results_Table[i,5])
  Results_Table[i,8] <- as.numeric(Results_Table[1,7]) - as.numeric(Results_Table[i,7])
  Results_Table[i,10] <- as.numeric(Results_Table[i,9]) - as.numeric(Results_Table[1,9])
}


view(Results_Table)

# Plot mortality scenarios ####

mortality_base$t <- (mortality_base$t+2004)
mortality_base <- mortality_base %>%  pivot_longer(names_to = "group", cols = !1)
mortality_base$scenario <- rep("baseline",nrow(mortality_base))

mortality_1$t <- (mortality_1$t+2004)
mortality_1 <- mortality_1 %>%  pivot_longer(names_to = "group", cols = !1)
mortality_1$scenario <- rep("1",nrow(mortality_1))

mortality_2$t <- (mortality_2$t+2004)
mortality_2 <- mortality_2%>%  pivot_longer(names_to = "group", cols = !1)
mortality_2$scenario <- rep("2",nrow(mortality_2))

mortality_3$t <- (mortality_3$t+2004)
mortality_3 <- mortality_3 %>%  pivot_longer(names_to = "group", cols = !1)
mortality_3$scenario <- rep("3",nrow(mortality_3))

mortality_compare <- rbind(mortality_base, mortality_1, mortality_2, mortality_3) 

mortality_compare[mortality_compare == "rate1"] <- "0 - 4"
mortality_compare[mortality_compare == "rate2"] <- "5 - 9"
mortality_compare[mortality_compare == "rate3"] <- "10 - 14"
mortality_compare[mortality_compare == "rate4"] <- "15 - 19"
mortality_compare[mortality_compare == "rate5"] <- "20 - 24"
mortality_compare[mortality_compare == "rate6"] <- "25 - 29"
mortality_compare[mortality_compare == "rate7"] <- "30 - 34"
mortality_compare[mortality_compare == "rate8"] <- "35 - 39"
mortality_compare[mortality_compare == "rate9"] <- "40 - 44"
mortality_compare[mortality_compare == "rate10"] <- "45 - 49"
mortality_compare[mortality_compare == "rate11"] <- "50 - 54"
mortality_compare[mortality_compare == "rate12"] <- "55 - 59"
mortality_compare[mortality_compare == "rate13"] <- "60 - 64"
mortality_compare[mortality_compare == "rate14"] <- "65 - 69"
mortality_compare[mortality_compare == "rate15"] <- "70 - 74"
mortality_compare[mortality_compare == "rate16"] <- "75 - 79"
mortality_compare[mortality_compare == "rate17"] <- "80 - 84"
mortality_compare[mortality_compare == "rate18"] <- "85 - 89"
mortality_compare[mortality_compare == "rate19"] <- "90 - 94"
mortality_compare[mortality_compare == "rate20"] <- "95 - 99"
mortality_compare[mortality_compare == "rate21"] <- "Over 100"

mortality_compare$group <- factor(mortality_compare$group, levels = c('0 - 4', '5 - 9', '10 - 14', '15 - 19', '20 - 24', '25 - 29', '30 - 34', '35 - 39', '40 - 44', '45 - 49', '50 - 54', '55 - 59', '60 - 64', '65 - 69', '70 - 74', '75 - 79', '80 - 84', '85 - 89', '90 - 94', '95 - 99', 'Over 100'))


mortality_compare %>%
  ggplot()+
  geom_line(size = 1, aes(x = t, y=value, colour = scenario))+
  theme_minimal(base_size = 12) +
  facet_wrap(~group, scales="free_y") +
  theme(legend.key.width=unit(0.4,"cm"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  scale_colour_discrete(labels=c('Decline', 'Plateau', 'Growth', 'Baseline')) +
  scale_y_continuous(labels=scaleFUN2) +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "italic"
    ),
    legend.text=element_text(size=14)
  )+
  labs(title = "Population Scenarios by Age Group", x="Year", y =("Mortality Rate (Deaths per Person in Group per Year"), colour="Population Scenario")




# Plot population structure over time ####

age_struc_long <- age_struc %>% pivot_longer(names_to = "age_group", cols = !1)

age_struc_long[age_struc_long == "group1"] <- "0 - 4"
age_struc_long[age_struc_long == "group2"] <- "5 - 9"
age_struc_long[age_struc_long == "group3"] <- "10 - 14"
age_struc_long[age_struc_long == "group4"] <- "15 - 19"
age_struc_long[age_struc_long == "group5"] <- "20 - 24"
age_struc_long[age_struc_long == "group6"] <- "25 - 29"
age_struc_long[age_struc_long == "group7"] <- "30 - 34"
age_struc_long[age_struc_long == "group8"] <- "35 - 39"
age_struc_long[age_struc_long == "group9"] <- "40 - 44"
age_struc_long[age_struc_long == "group10"] <- "45 - 49"
age_struc_long[age_struc_long == "group11"] <- "50 - 54"
age_struc_long[age_struc_long == "group12"] <- "55 - 59"
age_struc_long[age_struc_long == "group13"] <- "60 - 64"
age_struc_long[age_struc_long == "group14"] <- "65 - 69"
age_struc_long[age_struc_long == "group15"] <- "70 - 74"
age_struc_long[age_struc_long == "group16"] <- "75 - 79"
age_struc_long[age_struc_long == "group17"] <- "80 - 84"
age_struc_long[age_struc_long == "group18"] <- "85 - 89"
age_struc_long[age_struc_long == "group19"] <- "90 - 94"
age_struc_long[age_struc_long == "group20"] <- "95 - 99"
age_struc_long[age_struc_long == "group21"] <- "Over 100"

age_struc_long$age_group <- factor(age_struc_long$age_group, levels = c('0 - 4', '5 - 9', '10 - 14', '15 - 19', '20 - 24', '25 - 29', '30 - 34', '35 - 39', '40 - 44', '45 - 49', '50 - 54', '55 - 59', '60 - 64', '65 - 69', '70 - 74', '75 - 79', '80 - 84', '85 - 89', '90 - 94', '95 - 99', 'Over 100'))

age_struc_long %>% 
  ggplot(aes(x=age_group,y=value,fill=as.factor(Year))) +
  geom_bar(stat="identity",position="dodge2") +
  labs(title="Population Structure Over Time", x="Age Group",y="Proportion of Population")+
  theme_minimal(base_size = 13)+
  guides(fill=guide_legend(nrow=10,byrow=TRUE,title="Year"))+
  scale_fill_viridis_d(direction=-1,option="D")


# Plot birthrate over time with population projection ####

pop_data <- cbind(pop_data,rep("data",nrow(pop_data)))
pop_proj <- cbind(pop_proj,rep("proj",nrow(pop_proj)))
names(pop_data)[names(pop_data) == colnames(pop_data)[1]] <- "Year"
names(pop_data)[names(pop_data) == colnames(pop_data)[2]] <- "value"
names(pop_data)[names(pop_data) == colnames(pop_data)[3]] <- "type"
names(pop_proj)[names(pop_proj) == colnames(pop_proj)[1]] <- "Year"
names(pop_proj)[names(pop_proj) == colnames(pop_proj)[2]] <- "value"
names(pop_proj)[names(pop_proj) == colnames(pop_proj)[3]] <- "type"
pop_data$type2 <- rep("pop", nrow(pop_data))
pop_proj$type2 <- rep("pop", nrow(pop_proj))

birth.func <- approxfun(birth.approx$t,birth.approx$birth,method="linear")
birthrate_data <- birth.func(0:18)
birthrate_proj <- birth.func(19:36)
birthrate_data <- cbind(c(2004:2022),birthrate_data,rep("data",length(birthrate_data)),rep("birth",length(birthrate_data)))
birthrate_proj <- cbind(c(2023:2040),birthrate_proj,rep("proj",length(birthrate_proj)),rep("birth",length(birthrate_proj)))
colnames(birthrate_data) <- c("Year", "value", "type", "type2")
colnames(birthrate_proj) <- c("Year"," value", "type", "type2")


birthrate_compare <- as.data.frame(rbind(birthrate_data,birthrate_proj))
birthrate_compare[,1] <- as.numeric(birthrate_compare[,1])
birthrate_compare[,2] <- as.numeric(birthrate_compare[,2])

pop_and_birth <- rbind(pop_data,pop_proj,birthrate_compare)

ggplot(pop_and_birth, aes(x = Year, y = value, colour = type)) +
  geom_point(size=2) +
  theme_minimal(base_size=14) +
  theme(legend.title=element_blank())+
  scale_y_continuous(labels=scales::comma) +
  facet_wrap(~type2, scales = "free_y", nrow = 2, 
             strip.position = "left", 
             labeller = as_labeller(c(pop = "Population", birth = "Birth Rate (Births per Person per Year)") ) )  +
  ylab(NULL) +
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  theme(legend.text=element_text(size=15))+
  scale_shape_discrete(labels=c('Birth Rate','Total Population'))+
  scale_colour_discrete(labels=c('Data', 'Projection')) +
  labs(title = "Birth Rate and Population", x="Year", y =(""), colour="type")

# Beta matrix heat map plot ####

beta_matrix <- as.data.frame(cbind(c(1:21),beta))

beta_matrix  <- beta_matrix  %>% pivot_longer(names_to = "group", cols = !1)
beta_matrix  <- cbind(beta_matrix [,1],rep(1:21,21),beta_matrix [,3])
names(beta_matrix )[names(beta_matrix ) == colnames(beta_matrix )[1]] <- "X"
names(beta_matrix )[names(beta_matrix ) == colnames(beta_matrix )[2]] <- "Y"
names(beta_matrix )[names(beta_matrix ) == colnames(beta_matrix )[3]] <- "Z"
beta_matrix  <- beta_matrix  %>% mutate(X=5*X-1, Y=5*Y-1)

ggplot(beta_matrix , aes(Y, X, fill= Z)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE, direction=1) +
  labs(title = "Beta Matrix: Derived from Sexual Contact and HPV in Laos", x="Age", y =("Age"), fill = "Transmission Coefficient") +
  theme_minimal(base_size = 18)

# Plot total population of model to compare with data (4 mortality scenarios) ####

names(pop_data1)[names(pop_data1) == colnames(pop_data1)[1]] <- "Year"
names(pop_data1)[names(pop_data1) == colnames(pop_data1)[2]] <- "value"
names(pop_data1)[names(pop_data1) == colnames(pop_data1)[3]] <- "type"
pop_data1$mort_scenario <- rep("none",nrow(pop_data1))

pop_model_base <-  df1_base_scr_base %>% filter(variable %in% c("population"))
pop_model_base <- cbind(pop_model_base,rep("model",nrow(pop_model_base)))
pop_model_base <- cbind(pop_model_base[,1],pop_model_base[,3:4])
names(pop_model_base)[names(pop_model_base) == colnames(pop_model_base)[1]] <- "Year"
names(pop_model_base)[names(pop_model_base) == colnames(pop_model_base)[3]] <- "type"
pop_model_base$mort_scenario <- rep("baseline",nrow(pop_model_base))

pop_model_1 <-  df1_mort_1_scr_base %>% filter(variable %in% c("population"))
pop_model_1 <- cbind(pop_model_1,rep("model",nrow(pop_model_1)))
pop_model_1 <- cbind(pop_model_1[,1],pop_model_1[,3:4])
names(pop_model_1)[names(pop_model_1) == colnames(pop_model_1)[1]] <- "Year"
names(pop_model_1)[names(pop_model_1) == colnames(pop_model_1)[3]] <- "type"
pop_model_1$mort_scenario <- rep("1",nrow(pop_model_1))

pop_model_2 <-  df1_mort_2_scr_base %>% filter(variable %in% c("population"))
pop_model_2 <- cbind(pop_model_2,rep("model",nrow(pop_model_2)))
pop_model_2 <- cbind(pop_model_2[,1],pop_model_2[,3:4])
names(pop_model_2)[names(pop_model_2) == colnames(pop_model_2)[1]] <- "Year"
names(pop_model_2)[names(pop_model_2) == colnames(pop_model_2)[3]] <- "type"
pop_model_2$mort_scenario <- rep("2",nrow(pop_model_2))

pop_model_3 <-  df1_mort_3_scr_base %>% filter(variable %in% c("population"))
pop_model_3 <- cbind(pop_model_3,rep("model",nrow(pop_model_3)))
pop_model_3 <- cbind(pop_model_3[,1],pop_model_3[,3:4])
names(pop_model_3)[names(pop_model_3) == colnames(pop_model_3)[1]] <- "Year"
names(pop_model_3)[names(pop_model_3) == colnames(pop_model_3)[3]] <- "type"
pop_model_3$mort_scenario <- rep("3",nrow(pop_model_3))

pop_compare <- rbind(pop_data1, pop_model_base,pop_model_1, pop_model_2, pop_model_3)

pop_brm_lower <- df1_base_scr_base_brm_lower %>% filter(variable %in% c("population")) %>% select(-2)
pop_brm_upper <- (df1_base_scr_base_brm_upper %>% filter(variable %in% c("population")))$value
pop_brm_limits <- cbind(pop_brm_lower,pop_brm_upper)
names(pop_brm_limits)[names(pop_brm_limits) == colnames(pop_brm_limits)[2]] <- "lower"
names(pop_brm_limits)[names(pop_brm_limits) == colnames(pop_brm_limits)[3]] <- "upper"


pop_compare %>% 
  ggplot() +
  geom_line(size=1.3, data = subset(pop_compare, type == "model"), aes(x = Year, y = value, colour = mort_scenario)) +
  geom_point(size=3, data = subset(pop_compare, type != "model"), aes(x = Year, y = value, colour = mort_scenario, shape=type)) +
  theme_minimal(base_size = 18) +
  scale_y_continuous(labels=scales::comma) +
  geom_ribbon(pop_brm_limits, mapping=aes(x=Year,
                                          ymin = lower,
                                          ymax = upper),
              fill = "azure4",
              alpha = 0.3
  ) +
  scale_shape_manual(values=c(1,2),
                     limits=c('data','projection'),
                     labels=c('Total Population Data','Total Population Projection'))+
  scale_colour_manual(values=c("#F8766D","#C49A00","#53B400","#00C094"),
                      limits=c('1','2','3','baseline'),
                      labels=c('Decline','Plateau','Growth','Baseline'))+
  labs(title = "Population of Thailand", x="Year", y =("Population"), colour="Model Output", shape="Data Point")


# Plot population structure of model to compare with data for each age group (4 mortality scenarios) ####

names(age_struc_proportion_data)[names(age_struc_proportion_data) == colnames(age_struc_proportion_data)[1]] <- "Year"
UN_data_groups <- age_struc_proportion_data %>% pivot_longer(names_to = "variable", cols = !1)
UN_data_groups$type <- rep("data",nrow(UN_data_groups))

names(age_struc_proportion_proj)[names(age_struc_proportion_proj) == colnames(age_struc_proportion_proj)[1]] <- "Year"
age_struc_proportion_proj <- age_struc_proportion_proj %>% filter(Year %in% c("2022","2023","2024","2025","2026","2027","2028","2029","2030","2031","2032","2033","2034","2035","2036","2037","2038","2039","2040"))
proj_groups <- age_struc_proportion_proj %>% pivot_longer(names_to = "variable", cols = !1)
proj_groups$type <- rep("projection",nrow(proj_groups))

data_groups <- rbind(UN_data_groups,proj_groups)

model_groups_base <- df1_base_scr_base %>% filter(variable %in% c("group1","group2","group3","group4","group5","group6","group7","group8","group9","group10","group11","group12","group13","group14","group15","group16","group17","group18","group19","group20","group21"))
model_groups_base$type <- rep("baseline",nrow(model_groups_base))

groups_compare <- rbind(data_groups, model_groups_base)

groups_compare[groups_compare == "group1"] <- "0 - 4"
groups_compare[groups_compare == "group2"] <- "5 - 9"
groups_compare[groups_compare == "group3"] <- "10 - 14"
groups_compare[groups_compare == "group4"] <- "15 - 19"
groups_compare[groups_compare == "group5"] <- "20 - 24"
groups_compare[groups_compare == "group6"] <- "25 - 29"
groups_compare[groups_compare == "group7"] <- "30 - 34"
groups_compare[groups_compare == "group8"] <- "35 - 39"
groups_compare[groups_compare == "group9"] <- "40 - 44"
groups_compare[groups_compare == "group10"] <- "45 - 49"
groups_compare[groups_compare == "group11"] <- "50 - 54"
groups_compare[groups_compare == "group12"] <- "55 - 59"
groups_compare[groups_compare == "group13"] <- "60 - 64"
groups_compare[groups_compare == "group14"] <- "65 - 69"
groups_compare[groups_compare == "group15"] <- "70 - 74"
groups_compare[groups_compare == "group16"] <- "75 - 79"
groups_compare[groups_compare == "group17"] <- "80 - 84"
groups_compare[groups_compare == "group18"] <- "85 - 89"
groups_compare[groups_compare == "group19"] <- "90 - 94"
groups_compare[groups_compare == "group20"] <- "95 - 99"
groups_compare[groups_compare == "group21"] <- "Over 100"

groups_compare$variable <- factor(groups_compare$variable, levels = c('0 - 4', '5 - 9', '10 - 14', '15 - 19', '20 - 24', '25 - 29', '30 - 34', '35 - 39', '40 - 44', '45 - 49', '50 - 54', '55 - 59', '60 - 64', '65 - 69', '70 - 74', '75 - 79', '80 - 84', '85 - 89', '90 - 94', '95 - 99', 'Over 100'))

pop_groups_brm_lower <- df1_base_scr_base_brm_lower %>% filter(variable %in% c("group1","group2","group3","group4","group5","group6","group7","group8","group9","group10","group11","group12","group13","group14","group15","group16","group17","group18","group19","group20","group21"))
pop_groups_brm_upper <- (df1_base_scr_base_brm_upper %>% filter(variable %in% c("group1","group2","group3","group4","group5","group6","group7","group8","group9","group10","group11","group12","group13","group14","group15","group16","group17","group18","group19","group20","group21")))$value
pop_groups_brm_limits <- cbind(pop_groups_brm_lower,pop_groups_brm_upper)
names(pop_groups_brm_limits)[names(pop_groups_brm_limits) == colnames(pop_groups_brm_limits)[3]] <- "lower"
names(pop_groups_brm_limits)[names(pop_groups_brm_limits) == colnames(pop_groups_brm_limits)[4]] <- "upper"

pop_groups_brm_limits[pop_groups_brm_limits == "group1"] <- "0 - 4"
pop_groups_brm_limits[pop_groups_brm_limits == "group2"] <- "5 - 9"
pop_groups_brm_limits[pop_groups_brm_limits == "group3"] <- "10 - 14"
pop_groups_brm_limits[pop_groups_brm_limits == "group4"] <- "15 - 19"
pop_groups_brm_limits[pop_groups_brm_limits == "group5"] <- "20 - 24"
pop_groups_brm_limits[pop_groups_brm_limits == "group6"] <- "25 - 29"
pop_groups_brm_limits[pop_groups_brm_limits == "group7"] <- "30 - 34"
pop_groups_brm_limits[pop_groups_brm_limits == "group8"] <- "35 - 39"
pop_groups_brm_limits[pop_groups_brm_limits == "group9"] <- "40 - 44"
pop_groups_brm_limits[pop_groups_brm_limits == "group10"] <- "45 - 49"
pop_groups_brm_limits[pop_groups_brm_limits == "group11"] <- "50 - 54"
pop_groups_brm_limits[pop_groups_brm_limits == "group12"] <- "55 - 59"
pop_groups_brm_limits[pop_groups_brm_limits == "group13"] <- "60 - 64"
pop_groups_brm_limits[pop_groups_brm_limits == "group14"] <- "65 - 69"
pop_groups_brm_limits[pop_groups_brm_limits == "group15"] <- "70 - 74"
pop_groups_brm_limits[pop_groups_brm_limits == "group16"] <- "75 - 79"
pop_groups_brm_limits[pop_groups_brm_limits == "group17"] <- "80 - 84"
pop_groups_brm_limits[pop_groups_brm_limits == "group18"] <- "85 - 89"
pop_groups_brm_limits[pop_groups_brm_limits == "group19"] <- "90 - 94"
pop_groups_brm_limits[pop_groups_brm_limits == "group20"] <- "95 - 99"
pop_groups_brm_limits[pop_groups_brm_limits == "group21"] <- "Over 100"

pop_groups_brm_limits$variable <- factor(pop_groups_brm_limits$variable, levels = c('0 - 4', '5 - 9', '10 - 14', '15 - 19', '20 - 24', '25 - 29', '30 - 34', '35 - 39', '40 - 44', '45 - 49', '50 - 54', '55 - 59', '60 - 64', '65 - 69', '70 - 74', '75 - 79', '80 - 84', '85 - 89', '90 - 94', '95 - 99', 'Over 100'))

groups_compare %>% 
  ggplot()+
  geom_line(size=1.8, data = groups_compare %>% filter(type %nin% c("data","projection")),aes(x = Year, y = value, colour = type)) +
  geom_point(size=2, shape = 20, data = groups_compare %>% filter(type %in% c("data","projection")),aes(x = Year, y = value, colour = type)) +
  guides(color = guide_legend(override.aes = list(linetype = c(1, NA, NA),
                                                  size = 3,
                                                  shape    = c(NA, 20, 20)))) +
  theme_minimal(base_size = 16) +
  facet_wrap(~variable) +
  scale_colour_discrete(labels=c('Baseline', 'Data', 'Projection')) +
  scale_y_continuous(labels=scales::comma) +
  theme(
    strip.text.x = element_text(
      size = 16, color = "black", face = "italic"
    ),
    legend.text=element_text(size=16)
  )+
  labs(title = "Model vs. Data by Age Group", x="Year", y =("Size of Population"), colour="Population Scenario")


# Plot prevalence of model to compare with prevalence data from 2004 and 2014 ####

data_prev$lower_95 <- data_prev$value-data_prev$CI95
data_prev$upper_95 <- data_prev$value+data_prev$CI95
data_prev <- select(data_prev, -4,-5,-6)
data_prev$type <- rep("data", nrow(data_prev))

model_prev <- df1_base_scr_base %>% filter(variable %in% c("prev1_2","prev3_4","prev5_6","prev7_8","prev9_10","prev11_21", "infectprev"))
model_prev_upper <- df1_base_scr_base_upper %>% filter(variable %in% c("prev1_2","prev3_4","prev5_6","prev7_8","prev9_10","prev11_21", "infectprev"))
model_prev_lower <- df1_base_scr_base_lower %>% filter(variable %in% c("prev1_2","prev3_4","prev5_6","prev7_8","prev9_10","prev11_21", "infectprev"))
model_prev$lower_95 <- model_prev_lower$value
model_prev$upper_95 <- model_prev_upper$value

model_prev$type <- rep("model",nrow(model_prev))

prev_compare_2004 <- rbind(data_prev %>% filter(Year %in% c("2004")), model_prev %>% filter(Year %in% c("2004")))

prev_compare_2004[prev_compare_2004 == "infectprev"] <- "Total"
prev_compare_2004[prev_compare_2004 == "prev1_2"] <- "0 - 9"
prev_compare_2004[prev_compare_2004 == "prev3_4"] <- "10 - 19"
prev_compare_2004[prev_compare_2004 == "prev5_6"] <- "20 - 29"
prev_compare_2004[prev_compare_2004 == "prev7_8"] <- "30 - 39"
prev_compare_2004[prev_compare_2004 == "prev9_10"] <- "40 - 49"
prev_compare_2004[prev_compare_2004 == "prev11_21"] <- "Over 50"

prev_compare_2004$variable <- factor(prev_compare_2004$variable, levels = c('0 - 9','10 - 19','20 - 29','30 - 39','40 - 49','Over 50','Total'))


plot_2004_legend <- ggplot(prev_compare_2004, aes(x=variable, y=value, fill=type)) + 
  geom_bar(stat="identity",position="dodge2") +
  geom_errorbar(aes(ymin=upper_95, ymax=lower_95), width=.2,
                position=position_dodge(.9))+
  scale_y_continuous(labels=scaleFUN, limits=c(0,4)) +
  scale_fill_discrete(labels=c("Data","Model"))+
  labs(title="2004", x="Age Group",y="Thailand Prevalence (%)",fill="")+
  theme_minimal(base_size = 18)

plot_2004 <- ggplot(prev_compare_2004, aes(x=variable, y=value, fill=type)) + 
  geom_bar(stat="identity",position="dodge2") +
  geom_errorbar(aes(ymin=upper_95, ymax=lower_95), width=.2,
                position=position_dodge(.9))+
  scale_y_continuous(labels=scaleFUN, limits=c(0,4)) +
  scale_fill_discrete(labels=c("Data","Model"))+
  labs(title="2004", x="Age Group",y="Thailand Prevalence (%)",fill="")+
  theme_minimal(base_size = 13)+
  theme(legend.position = "none")

prev_compare_2014 <- rbind(data_prev %>% filter(Year %in% c("2014")), model_prev %>% filter(Year %in% c("2014")))

prev_compare_2014[prev_compare_2014 == "infectprev"] <- "Total"
prev_compare_2014[prev_compare_2014 == "prev1_2"] <- "0 - 9"
prev_compare_2014[prev_compare_2014 == "prev3_4"] <- "10 - 19"
prev_compare_2014[prev_compare_2014 == "prev5_6"] <- "20 - 29"
prev_compare_2014[prev_compare_2014 == "prev7_8"] <- "30 - 39"
prev_compare_2014[prev_compare_2014 == "prev9_10"] <- "40 - 49"
prev_compare_2014[prev_compare_2014 == "prev11_21"] <- "Over 50"

prev_compare_2014$variable <- factor(prev_compare_2014$variable, levels = c('0 - 9','10 - 19','20 - 29','30 - 39','40 - 49','Over 50','Total'))

prev_2022_model <- (df1_base_scr_base %>% filter(Year %in% c("2022")) %>% filter(variable%in%c("infectprev")))$value
prev_2022_model_upper <- signif((df1_base_scr_base_upper %>% filter(Year %in% c("2022")) %>% filter(variable%in%c("infectprev")))$value,3)
prev_2022_model_lower <- signif((df1_base_scr_base_lower %>% filter(Year %in% c("2022")) %>% filter(variable%in%c("infectprev")))$value,3)
prev_2022_data <- (data_prev %>% filter(Year %in% c("2022")))$value

labels <- as.character(c(signif(prev_2022_data,3),signif(prev_2022_model,3),prev_2022_model_lower,prev_2022_model_upper))
labels <- c(labels[1], paste(labels[2], " (",labels[3], " - ", labels[4], ")", sep=""))

plot_2014 <- ggplot(prev_compare_2014, aes(x=variable, y=value, fill=type)) + 
  geom_bar(stat="identity",position="dodge2") +
  annotate("text", x=2.6, y=3.4, label="2022 total prevalence = ", colour = "#00BFC4", fontface=2, size=5)+
  annotate("text", x=5.7, y=3.4, label=labels[2], colour = "#00BFC4", fontface=2, size=5)+
  annotate("text", x=2.6, y=3.7, label="2022 total prevalence = ", colour = "#F8766D", fontface=2, size=5)+
  annotate("text", x=4.6, y=3.7, label=labels[1], colour = "#F8766D", fontface=2, size=5)+
  geom_errorbar(aes(ymin=upper_95, ymax=lower_95), width=.2,
                position=position_dodge(.9))+
  scale_y_continuous(labels=scaleFUN, limits=c(0,4)) +
  labs(title="2014", x="Age Group",y="Thailand Prevalence (%)",fill="")+
  theme_minimal(base_size = 13)+
  theme(legend.position = "none")


legend_2004 <- get_only_legend(plot_2004_legend)  

lay <- rbind(c(1,1,1,1,2,2,2,2,3),
             c(1,1,1,1,2,2,2,2,3))

grid.arrange(plot_2004,plot_2014, legend_2004, layout_matrix = lay, top=textGrob("Prevalence: Model Baseline vs. Data", gp=gpar(fontsize=18)))


# Plot yearly incidence  of model to compare with 2030 target ####

# Baseline population scenario, all screening scenarios
model_inc_base_base <- df1_base_scr_base %>% filter(variable %in% c("Inc"))
model_inc_base_base$mort_scenario <- rep("mort_baseline",nrow(model_inc_base_base))
model_inc_base_base$scr_scenario <- rep("scr_baseline",nrow(model_inc_base_base))

model_inc_base_A <- df1_base_scr_A %>% filter(variable %in% c("Inc"))
model_inc_base_A$mort_scenario <- rep("mort_baseline",nrow(model_inc_base_A))
model_inc_base_A$scr_scenario <- rep("A",nrow(model_inc_base_A))

model_inc_base_B <- df1_base_scr_B %>% filter(variable %in% c("Inc"))
model_inc_base_B$mort_scenario <- rep("mort_baseline",nrow(model_inc_base_B))
model_inc_base_B$scr_scenario <- rep("B",nrow(model_inc_base_B))

model_inc_base_C <- df1_base_scr_C %>% filter(variable %in% c("Inc"))
model_inc_base_C$mort_scenario <- rep("mort_baseline",nrow(model_inc_base_C))
model_inc_base_C$scr_scenario <- rep("C",nrow(model_inc_base_C))

model_inc_base_D <- df1_base_scr_D %>% filter(variable %in% c("Inc"))
model_inc_base_D$mort_scenario <- rep("mort_baseline",nrow(model_inc_base_D))
model_inc_base_D$scr_scenario <- rep("D",nrow(model_inc_base_D))

model_inc_base_E <- df1_base_scr_E %>% filter(variable %in% c("Inc"))
model_inc_base_E$mort_scenario <- rep("mort_baseline",nrow(model_inc_base_E))
model_inc_base_E$scr_scenario <- rep("E",nrow(model_inc_base_E))

model_inc_base_F <- df1_base_scr_F %>% filter(variable %in% c("Inc"))
model_inc_base_F$mort_scenario <- rep("mort_baseline",nrow(model_inc_base_F))
model_inc_base_F$scr_scenario <- rep("F",nrow(model_inc_base_F))

# Population scenario 1, all screening scenarios
model_inc_mort_1_base <- df1_mort_1_scr_base %>% filter(variable %in% c("Inc"))
model_inc_mort_1_base$mort_scenario <- rep("1",nrow(model_inc_mort_1_base))
model_inc_mort_1_base$scr_scenario <- rep("scr_baseline",nrow(model_inc_mort_1_base))

model_inc_mort_1_A <- df1_mort_1_scr_A %>% filter(variable %in% c("Inc"))
model_inc_mort_1_A$mort_scenario <- rep("1",nrow(model_inc_mort_1_A))
model_inc_mort_1_A$scr_scenario <- rep("A",nrow(model_inc_mort_1_A))

model_inc_mort_1_B <- df1_mort_1_scr_B %>% filter(variable %in% c("Inc"))
model_inc_mort_1_B$mort_scenario <- rep("1",nrow(model_inc_mort_1_B))
model_inc_mort_1_B$scr_scenario <- rep("B",nrow(model_inc_mort_1_B))

model_inc_mort_1_C <- df1_mort_1_scr_C %>% filter(variable %in% c("Inc"))
model_inc_mort_1_C$mort_scenario <- rep("1",nrow(model_inc_mort_1_C))
model_inc_mort_1_C$scr_scenario <- rep("C",nrow(model_inc_mort_1_C))

model_inc_mort_1_D <- df1_mort_1_scr_D %>% filter(variable %in% c("Inc"))
model_inc_mort_1_D$mort_scenario <- rep("1",nrow(model_inc_mort_1_D))
model_inc_mort_1_D$scr_scenario <- rep("D",nrow(model_inc_mort_1_D))

model_inc_mort_1_E <- df1_mort_1_scr_E %>% filter(variable %in% c("Inc"))
model_inc_mort_1_E$mort_scenario <- rep("1",nrow(model_inc_mort_1_E))
model_inc_mort_1_E$scr_scenario <- rep("E",nrow(model_inc_mort_1_E))

model_inc_mort_1_F <- df1_mort_1_scr_F %>% filter(variable %in% c("Inc"))
model_inc_mort_1_F$mort_scenario <- rep("1",nrow(model_inc_mort_1_F))
model_inc_mort_1_F$scr_scenario <- rep("F",nrow(model_inc_mort_1_F))

# Population scenario 2, all screening scenarios
model_inc_mort_2_base <- df1_mort_2_scr_base %>% filter(variable %in% c("Inc"))
model_inc_mort_2_base$mort_scenario <- rep("2",nrow(model_inc_mort_2_base))
model_inc_mort_2_base$scr_scenario <- rep("scr_baseline",nrow(model_inc_mort_2_base))

model_inc_mort_2_A <- df1_mort_2_scr_A %>% filter(variable %in% c("Inc"))
model_inc_mort_2_A$mort_scenario <- rep("2",nrow(model_inc_mort_2_A))
model_inc_mort_2_A$scr_scenario <- rep("A",nrow(model_inc_mort_2_A))

model_inc_mort_2_B <- df1_mort_2_scr_B %>% filter(variable %in% c("Inc"))
model_inc_mort_2_B$mort_scenario <- rep("2",nrow(model_inc_mort_2_B))
model_inc_mort_2_B$scr_scenario <- rep("B",nrow(model_inc_mort_2_B))

model_inc_mort_2_C <- df1_mort_2_scr_C %>% filter(variable %in% c("Inc"))
model_inc_mort_2_C$mort_scenario <- rep("2",nrow(model_inc_mort_2_C))
model_inc_mort_2_C$scr_scenario <- rep("C",nrow(model_inc_mort_2_C))

model_inc_mort_2_D <- df1_mort_2_scr_D %>% filter(variable %in% c("Inc"))
model_inc_mort_2_D$mort_scenario <- rep("2",nrow(model_inc_mort_2_D))
model_inc_mort_2_D$scr_scenario <- rep("D",nrow(model_inc_mort_2_D))

model_inc_mort_2_E <- df1_mort_2_scr_E %>% filter(variable %in% c("Inc"))
model_inc_mort_2_E$mort_scenario <- rep("2",nrow(model_inc_mort_2_E))
model_inc_mort_2_E$scr_scenario <- rep("E",nrow(model_inc_mort_2_E))

model_inc_mort_2_F <- df1_mort_2_scr_F %>% filter(variable %in% c("Inc"))
model_inc_mort_2_F$mort_scenario <- rep("2",nrow(model_inc_mort_2_F))
model_inc_mort_2_F$scr_scenario <- rep("F",nrow(model_inc_mort_2_F))

# Population scenario 3, all screening scenarios
model_inc_mort_3_base <- df1_mort_3_scr_base %>% filter(variable %in% c("Inc"))
model_inc_mort_3_base$mort_scenario <- rep("3",nrow(model_inc_mort_3_base))
model_inc_mort_3_base$scr_scenario <- rep("scr_baseline",nrow(model_inc_mort_3_base))

model_inc_mort_3_A <- df1_mort_3_scr_A %>% filter(variable %in% c("Inc"))
model_inc_mort_3_A$mort_scenario <- rep("3",nrow(model_inc_mort_3_A))
model_inc_mort_3_A$scr_scenario <- rep("A",nrow(model_inc_mort_3_A))

model_inc_mort_3_B <- df1_mort_3_scr_B %>% filter(variable %in% c("Inc"))
model_inc_mort_3_B$mort_scenario <- rep("3",nrow(model_inc_mort_3_B))
model_inc_mort_3_B$scr_scenario <- rep("B",nrow(model_inc_mort_3_B))

model_inc_mort_3_C <- df1_mort_3_scr_C %>% filter(variable %in% c("Inc"))
model_inc_mort_3_C$mort_scenario <- rep("3",nrow(model_inc_mort_3_C))
model_inc_mort_3_C$scr_scenario <- rep("C",nrow(model_inc_mort_3_C))

model_inc_mort_3_D <- df1_mort_3_scr_D %>% filter(variable %in% c("Inc"))
model_inc_mort_3_D$mort_scenario <- rep("3",nrow(model_inc_mort_3_D))
model_inc_mort_3_D$scr_scenario <- rep("D",nrow(model_inc_mort_3_D))

model_inc_mort_3_E <- df1_mort_3_scr_E %>% filter(variable %in% c("Inc"))
model_inc_mort_3_E$mort_scenario <- rep("3",nrow(model_inc_mort_3_E))
model_inc_mort_3_E$scr_scenario <- rep("E",nrow(model_inc_mort_3_E))

model_inc_mort_3_F <- df1_mort_3_scr_F %>% filter(variable %in% c("Inc"))
model_inc_mort_3_F$mort_scenario <- rep("3",nrow(model_inc_mort_3_F))
model_inc_mort_3_F$scr_scenario <- rep("F",nrow(model_inc_mort_3_F))

model_inc_compare <- rbind(model_inc_base_base, model_inc_base_A, model_inc_base_B, model_inc_base_C, model_inc_base_D, model_inc_base_E, model_inc_base_F,
                           model_inc_mort_1_base, model_inc_mort_1_A, model_inc_mort_1_B, model_inc_mort_1_C, model_inc_mort_1_D, model_inc_mort_1_E, model_inc_mort_1_F,
                           model_inc_mort_2_base, model_inc_mort_2_A, model_inc_mort_2_B, model_inc_mort_2_C, model_inc_mort_2_D, model_inc_mort_2_E, model_inc_mort_2_F,
                           model_inc_mort_3_base, model_inc_mort_3_A, model_inc_mort_3_B, model_inc_mort_3_C, model_inc_mort_3_D, model_inc_mort_3_E, model_inc_mort_3_F)
model_inc_compare <- select(model_inc_compare,-2)

model_inc_compare$mort_scenario[model_inc_compare$mort_scenario == "mort_baseline"] <- "Baseline Population"
model_inc_compare$mort_scenario[model_inc_compare$mort_scenario == "1"] <- "Population Decline"
model_inc_compare$mort_scenario[model_inc_compare$mort_scenario == "2"] <- "Population Plateau"
model_inc_compare$mort_scenario[model_inc_compare$mort_scenario == "3"] <- "Population Growth"

model_inc_lower <- df1_base_scr_base_lower %>% filter(variable %in% c("Inc")) %>% select(-2)
model_inc_upper <- df1_base_scr_base_upper %>% filter(variable %in% c("Inc")) %>% select(-2)
model_inc_limits <- model_inc_lower
model_inc_limits$upper <- model_inc_upper$value
names(model_inc_limits)[names(model_inc_limits) == colnames(model_inc_limits)[2]] <- "lower"
model_inc_limits$upper[model_inc_limits$upper >= 6000] <- 6000
model_inc_limits <- model_inc_limits %>% filter(Year >= 2018)

model_inc_compare %>% filter(Year >= 2010) %>% 
  ggplot()+
  geom_line(size=1.1, aes(x = Year, y=value, colour = scr_scenario))+
  facet_wrap(~mort_scenario)+
  coord_cartesian(xlim = c(2010,2040), ylim = c(0, 6000))+
  geom_ribbon(model_inc_limits %>% filter(Year >= 2010),
              mapping=aes(x=Year,
                          ymin = lower,
                          ymax = upper),
              fill = "azure4",
              alpha = 0.2
  ) +
  geom_hline(yintercept=inc_target, linetype="dashed", color = "gray30") +
  geom_vline(xintercept=2030, linetype="dashed", color = "gray30") +
  geom_text(aes(2018, inc_target, label = "2030 Elimination Target", vjust = -1, hjust = 0.5), size = 4.5, color="gray30") +
  theme_minimal(base_size = 18) +
  scale_colour_discrete(labels=c('30 - 39 at 50%', '30 - 39 at 90%', '40 - 49 at 50%', '40 - 49 at 90%', '50 - 59 at 50%', '30 - 59 at 50%', 'Baseline')) +
  theme(
    strip.text.x = element_text(
      size = 18, color = "black", face = "italic"
    ),
    legend.text=element_text(size=16),
    legend.position = "bottom"
  )+
  labs(title = "Yearly Incidence", x="Year", y =("Number of Individuals"), colour="Screening Strategy")


# Plot Yearly deaths over time with 2030 target ####

# Baseline population scenario, all screening strategies
model_mort_base_base <- df1_base_scr_base %>% filter(variable %in% c("Deaths"))
model_mort_base_base$mort_scenario <- rep("mort_baseline",nrow(model_mort_base_base))
model_mort_base_base$scr_scenario <- rep("scr_baseline",nrow(model_mort_base_base))

model_mort_base_A <- df1_base_scr_A %>% filter(variable %in% c("Deaths"))
model_mort_base_A$mort_scenario <- rep("mort_baseline",nrow(model_mort_base_A))
model_mort_base_A$scr_scenario <- rep("A",nrow(model_mort_base_A))

model_mort_base_B <- df1_base_scr_B %>% filter(variable %in% c("Deaths"))
model_mort_base_B$mort_scenario <- rep("mort_baseline",nrow(model_mort_base_B))
model_mort_base_B$scr_scenario <- rep("B",nrow(model_mort_base_B))

model_mort_base_C <- df1_base_scr_C %>% filter(variable %in% c("Deaths"))
model_mort_base_C$mort_scenario <- rep("mort_baseline",nrow(model_mort_base_C))
model_mort_base_C$scr_scenario <- rep("C",nrow(model_mort_base_C))

model_mort_base_D <- df1_base_scr_D %>% filter(variable %in% c("Deaths"))
model_mort_base_D$mort_scenario <- rep("mort_baseline",nrow(model_mort_base_D))
model_mort_base_D$scr_scenario <- rep("D",nrow(model_mort_base_D))

model_mort_base_E <- df1_base_scr_E %>% filter(variable %in% c("Deaths"))
model_mort_base_E$mort_scenario <- rep("mort_baseline",nrow(model_mort_base_E))
model_mort_base_E$scr_scenario <- rep("E",nrow(model_mort_base_E))

model_mort_base_F <- df1_base_scr_F %>% filter(variable %in% c("Deaths"))
model_mort_base_F$mort_scenario <- rep("mort_baseline",nrow(model_mort_base_F))
model_mort_base_F$scr_scenario <- rep("F",nrow(model_mort_base_F))

# Population scenario 1, all screening scenarios
model_mort_mort_1_base <- df1_mort_1_scr_base %>% filter(variable %in% c("Deaths"))
model_mort_mort_1_base$mort_scenario <- rep("1",nrow(model_mort_mort_1_base))
model_mort_mort_1_base$scr_scenario <- rep("scr_baseline",nrow(model_mort_mort_1_base))

model_mort_mort_1_A <- df1_mort_1_scr_A %>% filter(variable %in% c("Deaths"))
model_mort_mort_1_A$mort_scenario <- rep("1",nrow(model_mort_mort_1_A))
model_mort_mort_1_A$scr_scenario <- rep("A",nrow(model_mort_mort_1_A))

model_mort_mort_1_B <- df1_mort_1_scr_B %>% filter(variable %in% c("Deaths"))
model_mort_mort_1_B$mort_scenario <- rep("1",nrow(model_mort_mort_1_B))
model_mort_mort_1_B$scr_scenario <- rep("B",nrow(model_mort_mort_1_B))

model_mort_mort_1_C <- df1_mort_1_scr_C %>% filter(variable %in% c("Deaths"))
model_mort_mort_1_C$mort_scenario <- rep("1",nrow(model_mort_mort_1_C))
model_mort_mort_1_C$scr_scenario <- rep("C",nrow(model_mort_mort_1_C))

model_mort_mort_1_D <- df1_mort_1_scr_D %>% filter(variable %in% c("Deaths"))
model_mort_mort_1_D$mort_scenario <- rep("1",nrow(model_mort_mort_1_D))
model_mort_mort_1_D$scr_scenario <- rep("D",nrow(model_mort_mort_1_D))

model_mort_mort_1_E <- df1_mort_1_scr_E %>% filter(variable %in% c("Deaths"))
model_mort_mort_1_E$mort_scenario <- rep("1",nrow(model_mort_mort_1_E))
model_mort_mort_1_E$scr_scenario <- rep("E",nrow(model_mort_mort_1_E))

model_mort_mort_1_F <- df1_mort_1_scr_F %>% filter(variable %in% c("Deaths"))
model_mort_mort_1_F$mort_scenario <- rep("1",nrow(model_mort_mort_1_F))
model_mort_mort_1_F$scr_scenario <- rep("F",nrow(model_mort_mort_1_F))

# Population scenario 2, all screening scenarios
model_mort_mort_2_base <- df1_mort_2_scr_base %>% filter(variable %in% c("Deaths"))
model_mort_mort_2_base$mort_scenario <- rep("2",nrow(model_mort_mort_2_base))
model_mort_mort_2_base$scr_scenario <- rep("scr_baseline",nrow(model_mort_mort_2_base))

model_mort_mort_2_A <- df1_mort_2_scr_A %>% filter(variable %in% c("Deaths"))
model_mort_mort_2_A$mort_scenario <- rep("2",nrow(model_mort_mort_2_A))
model_mort_mort_2_A$scr_scenario <- rep("A",nrow(model_mort_mort_2_A))

model_mort_mort_2_B <- df1_mort_2_scr_B %>% filter(variable %in% c("Deaths"))
model_mort_mort_2_B$mort_scenario <- rep("2",nrow(model_mort_mort_2_B))
model_mort_mort_2_B$scr_scenario <- rep("B",nrow(model_mort_mort_2_B))

model_mort_mort_2_C <- df1_mort_2_scr_C %>% filter(variable %in% c("Deaths"))
model_mort_mort_2_C$mort_scenario <- rep("2",nrow(model_mort_mort_2_C))
model_mort_mort_2_C$scr_scenario <- rep("C",nrow(model_mort_mort_2_C))

model_mort_mort_2_D <- df1_mort_2_scr_D %>% filter(variable %in% c("Deaths"))
model_mort_mort_2_D$mort_scenario <- rep("2",nrow(model_mort_mort_2_D))
model_mort_mort_2_D$scr_scenario <- rep("D",nrow(model_mort_mort_2_D))

model_mort_mort_2_E <- df1_mort_2_scr_E %>% filter(variable %in% c("Deaths"))
model_mort_mort_2_E$mort_scenario <- rep("2",nrow(model_mort_mort_2_E))
model_mort_mort_2_E$scr_scenario <- rep("E",nrow(model_mort_mort_2_E))

model_mort_mort_2_F <- df1_mort_2_scr_F %>% filter(variable %in% c("Deaths"))
model_mort_mort_2_F$mort_scenario <- rep("2",nrow(model_mort_mort_2_F))
model_mort_mort_2_F$scr_scenario <- rep("F",nrow(model_mort_mort_2_F))

# Population scenario 3, all screening scenarios
model_mort_mort_3_base <- df1_mort_3_scr_base %>% filter(variable %in% c("Deaths"))
model_mort_mort_3_base$mort_scenario <- rep("3",nrow(model_mort_mort_3_base))
model_mort_mort_3_base$scr_scenario <- rep("scr_baseline",nrow(model_mort_mort_3_base))

model_mort_mort_3_A <- df1_mort_3_scr_A %>% filter(variable %in% c("Deaths"))
model_mort_mort_3_A$mort_scenario <- rep("3",nrow(model_mort_mort_3_A))
model_mort_mort_3_A$scr_scenario <- rep("A",nrow(model_mort_mort_3_A))

model_mort_mort_3_B <- df1_mort_3_scr_B %>% filter(variable %in% c("Deaths"))
model_mort_mort_3_B$mort_scenario <- rep("3",nrow(model_mort_mort_3_B))
model_mort_mort_3_B$scr_scenario <- rep("B",nrow(model_mort_mort_3_B))

model_mort_mort_3_C <- df1_mort_3_scr_C %>% filter(variable %in% c("Deaths"))
model_mort_mort_3_C$mort_scenario <- rep("3",nrow(model_mort_mort_3_C))
model_mort_mort_3_C$scr_scenario <- rep("C",nrow(model_mort_mort_3_C))

model_mort_mort_3_D <- df1_mort_3_scr_D %>% filter(variable %in% c("Deaths"))
model_mort_mort_3_D$mort_scenario <- rep("3",nrow(model_mort_mort_3_D))
model_mort_mort_3_D$scr_scenario <- rep("D",nrow(model_mort_mort_3_D))

model_mort_mort_3_E <- df1_mort_3_scr_E %>% filter(variable %in% c("Deaths"))
model_mort_mort_3_E$mort_scenario <- rep("3",nrow(model_mort_mort_3_E))
model_mort_mort_3_E$scr_scenario <- rep("E",nrow(model_mort_mort_3_E))

model_mort_mort_3_F <- df1_mort_3_scr_F %>% filter(variable %in% c("Deaths"))
model_mort_mort_3_F$mort_scenario <- rep("3",nrow(model_mort_mort_3_F))
model_mort_mort_3_F$scr_scenario <- rep("F",nrow(model_mort_mort_3_F))

model_mort_compare <- rbind(model_mort_base_base, model_mort_base_A, model_mort_base_B, model_mort_base_C, model_mort_base_D, model_mort_base_E, model_mort_base_F,
                            model_mort_mort_1_base, model_mort_mort_1_A, model_mort_mort_1_B, model_mort_mort_1_C, model_mort_mort_1_D, model_mort_mort_1_E, model_mort_mort_1_F,
                            model_mort_mort_2_base, model_mort_mort_2_A, model_mort_mort_2_B, model_mort_mort_2_C, model_mort_mort_2_D, model_mort_mort_2_E, model_mort_mort_2_F,
                            model_mort_mort_3_base, model_mort_mort_3_A, model_mort_mort_3_B, model_mort_mort_3_C, model_mort_mort_3_D, model_mort_mort_3_E, model_mort_mort_3_F)
model_mort_compare <- select(model_mort_compare,-2)

model_mort_compare$mort_scenario[model_mort_compare$mort_scenario == "mort_baseline"] <- "Baseline Population"
model_mort_compare$mort_scenario[model_mort_compare$mort_scenario == "1"] <- "Population Decline"
model_mort_compare$mort_scenario[model_mort_compare$mort_scenario == "2"] <- "Population Plateau"
model_mort_compare$mort_scenario[model_mort_compare$mort_scenario == "3"] <- "Population Growth"

model_mort_lower <- df1_base_scr_base_lower %>% filter(variable %in% c("Deaths")) %>% select(-2)
model_mort_upper <- df1_base_scr_base_upper %>% filter(variable %in% c("Deaths")) %>% select(-2)
model_mort_limits <- model_mort_lower
model_mort_limits$upper <- model_mort_upper$value
names(model_mort_limits)[names(model_mort_limits) == colnames(model_mort_limits)[2]] <- "lower"

model_mort_compare %>% 
  ggplot()+
  geom_line(size=1.1, aes(x = Year, y=value, colour = scr_scenario))+
  geom_point(hcv_deaths_data%>% filter(Year %in% c("2012","2013","2014","2015","2016","2017","2018","2019")), mapping = aes(x = Year, y = value))+
  geom_errorbar(hcv_deaths_data %>% filter(Year %in% c("2012","2013","2014","2015","2016","2017","2018","2019")), mapping = aes(x = Year, ymin = lower_95, ymax= upper_95)) +
  geom_ribbon(model_mort_limits, mapping=aes(x=Year,
                                             ymin = lower,
                                             ymax = upper),
              fill = "azure4",
              alpha = 0.3
  ) +
  facet_wrap(~mort_scenario)+
  ylim(0,15000)+
  coord_cartesian( xlim = c(2010,2040))+
  geom_hline(yintercept=mort_target, linetype="dashed", color = "gray30") +
  geom_vline(xintercept=2030, linetype="dashed", color = "gray30") +
  geom_text(aes(2018, mort_target, label = "2030 Elimination Target", vjust = -1, hjust = 0.5), size = 5, color="gray30") +
  theme_minimal(base_size = 18) +
  scale_colour_discrete(labels=c('30 - 39 at 50%', '30 - 39 at 90%', '40 - 49 at 50%', '40 - 49 at 90%', '50 - 59 at 50%', '30 - 59 at 50%', 'Baseline')) +
  theme(
    strip.text.x = element_text(
      size = 18, color = "black", face = "italic"
    ),
    legend.text=element_text(size=16),
    legend.position = "bottom"
  )+
  labs(title = "Yearly Deaths from HCV", x="Year", y =("Number of Deaths per Year"), colour="Screening Strategy")


