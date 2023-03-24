#================================================================#
# Stochastic modeling of swine influenza transmission
# Using GillespieSSA package
# https://cran.r-project.org/web/packages/GillespieSSA/GillespieSSA.pdf
# Created by EK on 2/2/2023 
# builds from Vaccination r scripts

# Using this file to test the combination of all CMs
# Adding: 1. Maternal Antibody protection for the first 3 weeks when piglets enter model (we cite 6 weeks of MDA protection in the introduction)
#         2. Death rate for infected pigs - setting this equal to the natural death rate
#================================================================#
#### LIBRARIES ####
# library(tidyverse)
# library(GillespieSSA)


#### MODEL CONSTANTS ####
#Hog constants
# beta_Hog <- 0.285   #high transmision
# beta_Hog <- 0.001   #low transmission
beta_Hog <- (10/(1000*5)) #R0 = 10
beta_Hog_ind <- beta_Hog / 178 #(following Etbaigha et al 2018 paper)
beta_Hog_MDA <- (10/(1000*5))/10
beta_Hog_ind_MDA <- beta_Hog_MDA / 178
# beta_Hog_ind <- beta_Hog / 500
beta_QHog_ind <- 0
Hog_Vacc_Eff <- 0.6
u_Hog <- 0.00028     #natural death rate (following Etbaigha et al 2018 paper)
u_inf_Hog <- u_Hog#0.1     #infected death rate (0 for now to problem solve)
w_Hog <- 0 #1/180         #(following Etbaigha et al 2018 paper)  White et al shows a range from 56 - 112 days
sigma_Hog <- 1/2   #Latent period h1n1 and h3n2
delta_Hog <- 1/5   #Infectious period h1n1 and h3n2
qrate <- 1/3       #Quarantine Rate (either 1, 1/2, or 1/3)
q_cap <- 10   #Quarantine Capcity per room

#Interspecies transmission
##applied the beta = Ro / N*D from the modeling book
# Ro from the following website for swine to human infection 
# https://bmcmedicine.biomedcentral.com/articles/10.1186/1741-7015-7-30#Sec6   (from super-strain figure) possible the high end of transmisilibty for this parameter
# beta_S2WF <- (2.3 / (4002*5))  #  #Hogs + #Workforce * Duration of Hogs
#Alt beta_S2WF from fair outbreak
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3557885/
#assuming one hog got the first three confirmed cases sick and infectious period of 5 days Ro = 3/5 = 0.6
beta_S2WF <- (0.6 / (4002*5))
beta_QS2WF <- 0

#when Ro WF2S = Ro S2WF
# beta_WF2S <- (2.3 / (4002*3))  #  #Hogs + #Workforce * Duration of Hogs
# WF2S calculated from Swine Outbreak of Pandemic Influenza A Virus on a Canadian Research Farm Supports Human-to-Swine Transmission paper
#r0 = 10 (hogs infected) / 4(the four sicke workers) / 3(infectious period of workers)
beta_WF2S <- (0.83 / (4002*3))     #Hogs + #Workforce * Duration of WF

#Human constants
#applied the beta = Ro / N*D from the modeling book
beta_WF2WF <- (1.92 / (2*3)) # http://m-hikari.com/ams/ams-2013/ams-41-44-2013/joseAMS41-44-2013.pdf 
Hum_Vacc_Eff <- 0.0 
u_WF <- 0     #Natural death rate of human (0 for now to problem solve)
u_inf_WF <- 0   #Infected death rate of human (0 for now to problem solve)
w_WF <- 0      #Recovery rate for human (0 for now to problem solve)
sigma_WF <- 1/2 #Latent period for humans
delta_WF <- 1/3 #Infectious Period for humans

#### MODEL SETUP ####
# Model Parameters
TransmissionModel.par <- c(beta_Hog = beta_Hog, beta_Hog_ind = beta_Hog_ind, beta_Hog_MDA = beta_Hog_MDA, beta_Hog_ind_MDA = beta_Hog_ind_MDA,
                           beta_S2WF = beta_S2WF, beta_WF2S = beta_WF2S,
                           qrate = qrate, q_cap = q_cap, beta_QS2WF = beta_QS2WF, beta_QHog_ind = beta_QHog_ind,
                           w_Hog = w_Hog, sigma_Hog = sigma_Hog, delta_Hog = delta_Hog, u_Hog = u_Hog, u_inf_Hog = u_inf_Hog,
                           beta_WF2WF = beta_WF2WF,
                           w_WF = w_WF, sigma_WF = sigma_WF, delta_WF = delta_WF, u_WF = u_WF, u_inf_WF = u_inf_WF)
# The Vacc Eff constants are not added to the parameters because they are only used to set initial values which occurs before the ssa function called

#Differential Equations
TransmissionModel_Room1MDA.ODE <- c(
  #Room 1
  #Susceptible
  "w_Hog*R1", "(beta_Hog_MDA*I1*S1)", "(beta_Hog_ind_MDA*(I2+I3+I4)*S1)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S1)", "(beta_WF2S*HI*S1)", "u_Hog*S1",
  #Exposed
  "sigma_Hog*E1", "u_Hog*E1",
  #Infected
  "ifelse(Q1<q_cap, (qrate*I1), 0)", "(delta_Hog*I1)", "u_inf_Hog*I1",
  #Quarantine
  "(delta_Hog*Q1)", "u_inf_Hog*Q1",
  #Recovered
  "u_Hog*R1", #"w_Hog*R1" "(delta_Hog*I1),
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E1"
  #True Recovered
  # "(delta_Hog*I1)"
  #Total Quarantined
  # "qrate*E1"
  #Natural Death Rate
  # "u_Hog*S1", "u_Hog*E1","u_Hog*R1",
  #Infected Death Rate
  #"u_inf_Hog*I1",
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S1)"
  
  #Room 2
  #Susceptible
  "w_Hog*R2", "(beta_Hog*I2*S2)", "(beta_Hog_ind*(I3+I4)*S2)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S2)", "(beta_WF2S*HI*S2)", "u_Hog*S2",
  #Exposed
  "sigma_Hog*E2", "u_Hog*E2",
  #Infected
  "ifelse(Q2<q_cap, (qrate*I2), 0)", "(delta_Hog*I2)", "u_inf_Hog*I2",
  #Quarantine
  "(delta_Hog*Q2)", "u_inf_Hog*Q2",
  #Recovered
  "u_Hog*R2",# "w_Hog*R2" "(delta_Hog*I2),  "(delta_Hog*Q2)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E1"
  #True Recovered
  # "(delta_Hog*I1)"
  #Total Quarantined
  # "qrate*E2"
  #Natural Death Rate
  # "u_Hog*S2", "u_Hog*E2","u_Hog*R2",
  #Infected Death Rate
  # "u_inf_Hog*I2", "u_inf_Hog*Q1"
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S2)"
  
  #Room 3
  #Susceptible
  "w_Hog*R3", "(beta_Hog*I3*S3)", "(beta_Hog_ind*(I4)*S3)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S3)", "(beta_WF2S*HI*S3)", "u_Hog*S3",
  #Exposed
  "sigma_Hog*E3", "u_Hog*E3",
  #Infected
  "ifelse(Q3<q_cap, (qrate*I3), 0)", "(delta_Hog*I3)", "u_inf_Hog*I3",
  #Quarantine
  "(delta_Hog*Q3)", "u_inf_Hog*Q3",
  #Recovered
  "u_Hog*R3",#"w_Hog*R3" "(delta_Hog*I3),"(delta_Hog*Q3)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E3"
  #True Recovered
  # "(delta_Hog*I3)"
  #Total Quarantined
  # "qrate*E3"
  #Natural Death Rate
  # "u_Hog*S3", "u_Hog*E3","u_Hog*R1",
  #Infected Death Rate
  # "u_inf_Hog*I3", "u_inf_Hog*Q3",
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S3)"
  
  #Room 4
  #Susceptible
  "w_Hog*R4", "(beta_Hog*I4*S4)", "(beta_Hog_ind*(0)*S4)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S4)", "(beta_WF2S*HI*S4)", "u_Hog*S4",
  #Exposed
  "sigma_Hog*E4", "u_Hog*E4",
  #Infected
  "ifelse(Q4<q_cap, (qrate*I4), 0)", "(delta_Hog*I4)", "u_inf_Hog*I4",
  #Quarantine
  "(delta_Hog*Q4)", "u_inf_Hog*Q4",
  #Recovered
  "u_Hog*R4", #"(delta_Hog*I4), "(delta_Hog*Q4)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E4"
  #True Recovered
  # "(delta_Hog*I4)"
  #Total Quarantined
  # "qrate*E4"
  #Natural Death Rate
  # "u_Hog*S4", "u_Hog*E4","u_Hog*R4",
  #Infected Death Rate
  # "u_inf_Hog*I4", "u_inf_Hog*Q4",
  #Count WF-to-Hog infections
  #"(beta_WF2S*HI*S4)"
  
  #Workforce
  #Susceptible
  "w_WF*HR", "(beta_S2WF*(I1+I2+I3+I4)*HS)", "(beta_QS2WF*(Q1+Q2+Q3+Q4)*HS)", "(beta_WF2WF*HI*HS)", "u_WF*HS",
  #Exposed
  "sigma_WF*HE", "u_WF*HE",
  #Infected
  "(delta_WF*HI)", "u_inf_WF*HI",
  #Recovered
  "u_WF*HR" #"w_Hog*R4" "(delta_hum*HI)"
  # #Natural Death Rate
  # "u_WF*HS", "u_WF*HE", "u_WF*HR",
  # #Infected Death Rate
  # "u_inf_WF*HI"
  #Count Hog-to-WF infections
  # "(beta_S2WF*(I1+I2+I3+I4)*HS)"
)

TransmissionModel_Room2MDA.ODE <- c(
  #Room 1
  #Susceptible
  "w_Hog*R1", "(beta_Hog*I1*S1)", "(beta_Hog_ind*(I2+I3+I4)*S1)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S1)", "(beta_WF2S*HI*S1)", "u_Hog*S1",
  #Exposed
  "sigma_Hog*E1", "u_Hog*E1",
  #Infected
  "ifelse(Q1<q_cap, (qrate*I1), 0)", "(delta_Hog*I1)", "u_inf_Hog*I1",
  #Quarantine
  "(delta_Hog*Q1)", "u_inf_Hog*Q1",
  #Recovered
  "u_Hog*R1", #"w_Hog*R1" "(delta_Hog*I1), "(delta_Hog*Q1)", 
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E1"
  #True Recovered
  # "(delta_Hog*I1)"
  #Total Quarantined
  # "qrate*E1"
  #Natural Death Rate
  # "u_Hog*S1", "u_Hog*E1","u_Hog*R1",
  #Infected Death Rate
  #"u_inf_Hog*I1", "u_inf_Hog*Q1"
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S1)"
  
  #Room 2
  #Susceptible
  "w_Hog*R2", "(beta_Hog_MDA*I2*S2)", "(beta_Hog_ind_MDA*(I3+I4)*S2)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S2)", "(beta_WF2S*HI*S2)", "u_Hog*S2",
  #Exposed
  "sigma_Hog*E2", "u_Hog*E2",
  #Infected
  "ifelse(Q2<q_cap, (qrate*I2), 0)", "(delta_Hog*I2)", "u_inf_Hog*I2",
  #Quarantine
  "(delta_Hog*Q2)", "u_inf_Hog*Q2",
  #Recovered
  "u_Hog*R2",# "w_Hog*R2" "(delta_Hog*I2),  "(delta_Hog*Q2)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E1"
  #True Recovered
  # "(delta_Hog*I1)"
  #Total Quarantined
  # "qrate*E2"
  #Natural Death Rate
  # "u_Hog*S2", "u_Hog*E2","u_Hog*R2",
  #Infected Death Rate
  # "u_inf_Hog*I2", "u_inf_Hog*Q1"
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S2)"
  
  #Room 3
  #Susceptible
  "w_Hog*R3", "(beta_Hog*I3*S3)", "(beta_Hog_ind*(I4)*S3)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S3)", "(beta_WF2S*HI*S3)", "u_Hog*S3",
  #Exposed
  "sigma_Hog*E3", "u_Hog*E3",
  #Infected
  "ifelse(Q3<q_cap, (qrate*I3), 0)", "(delta_Hog*I3)", "u_inf_Hog*I3",
  #Quarantine
  "(delta_Hog*Q3)", "u_inf_Hog*Q3",
  #Recovered
  "u_Hog*R3",#"w_Hog*R3" "(delta_Hog*I3),"(delta_Hog*Q3)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E3"
  #True Recovered
  # "(delta_Hog*I3)"
  #Total Quarantined
  # "qrate*E3"
  #Natural Death Rate
  # "u_Hog*S3", "u_Hog*E3","u_Hog*R1",
  #Infected Death Rate
  # "u_inf_Hog*I3", "u_inf_Hog*Q3",
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S3)"
  
  #Room 4
  #Susceptible
  "w_Hog*R4", "(beta_Hog*I4*S4)", "(beta_Hog_ind*(0)*S4)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S4)", "(beta_WF2S*HI*S4)", "u_Hog*S4",
  #Exposed
  "sigma_Hog*E4", "u_Hog*E4",
  #Infected
  "ifelse(Q4<q_cap, (qrate*I4), 0)", "(delta_Hog*I4)", "u_inf_Hog*I4",
  #Quarantine
  "(delta_Hog*Q4)", "u_inf_Hog*Q4",
  #Recovered
  "u_Hog*R4", #"(delta_Hog*I4), "(delta_Hog*Q4)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E4"
  #True Recovered
  # "(delta_Hog*I4)"
  #Total Quarantined
  # "qrate*E4"
  #Natural Death Rate
  # "u_Hog*S4", "u_Hog*E4","u_Hog*R4",
  #Infected Death Rate
  # "u_inf_Hog*I4", "u_inf_Hog*Q4",
  #Count WF-to-Hog infections
  #"(beta_WF2S*HI*S4)"
  
  #Workforce
  #Susceptible
  "w_WF*HR", "(beta_S2WF*(I1+I2+I3+I4)*HS)", "(beta_QS2WF*(Q1+Q2+Q3+Q4)*HS)", "(beta_WF2WF*HI*HS)", "u_WF*HS",
  #Exposed
  "sigma_WF*HE", "u_WF*HE",
  #Infected
  "(delta_WF*HI)", "u_inf_WF*HI",
  #Recovered
  "u_WF*HR" #"w_Hog*R4" "(delta_hum*HI)"
  # #Natural Death Rate
  # "u_WF*HS", "u_WF*HE", "u_WF*HR",
  # #Infected Death Rate
  # "u_inf_WF*HI"
  #Count Hog-to-WF infections
  # "(beta_S2WF*(I1+I2+I3+I4)*HS)"
)
TransmissionModel_Room3MDA.ODE <- c(
  #Room 1
  #Susceptible
  "w_Hog*R1", "(beta_Hog*I1*S1)", "(beta_Hog_ind*(I2+I3+I4)*S1)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S1)", "(beta_WF2S*HI*S1)", "u_Hog*S1",
  #Exposed
  "sigma_Hog*E1", "u_Hog*E1",
  #Infected
  "ifelse(Q1<q_cap, (qrate*I1), 0)", "(delta_Hog*I1)", "u_inf_Hog*I1",
  #Quarantine
  "(delta_Hog*Q1)", "u_inf_Hog*Q1",
  #Recovered
  "u_Hog*R1", #"w_Hog*R1" "(delta_Hog*I1), "(delta_Hog*Q1)", 
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E1"
  #True Recovered
  # "(delta_Hog*I1)"
  #Total Quarantined
  # "qrate*E1"
  #Natural Death Rate
  # "u_Hog*S1", "u_Hog*E1","u_Hog*R1",
  #Infected Death Rate
  #"u_inf_Hog*I1", "u_inf_Hog*Q1"
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S1)"
  
  #Room 2
  #Susceptible
  "w_Hog*R2", "(beta_Hog*I2*S2)", "(beta_Hog_ind*(I3+I4)*S2)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S2)", "(beta_WF2S*HI*S2)", "u_Hog*S2",
  #Exposed
  "sigma_Hog*E2", "u_Hog*E2",
  #Infected
  "ifelse(Q2<q_cap, (qrate*I2), 0)", "(delta_Hog*I2)", "u_inf_Hog*I2",
  #Quarantine
  "(delta_Hog*Q2)", "u_inf_Hog*Q2",
  #Recovered
  "u_Hog*R2",# "w_Hog*R2" "(delta_Hog*I2),  "(delta_Hog*Q2)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E1"
  #True Recovered
  # "(delta_Hog*I1)"
  #Total Quarantined
  # "qrate*E2"
  #Natural Death Rate
  # "u_Hog*S2", "u_Hog*E2","u_Hog*R2",
  #Infected Death Rate
  # "u_inf_Hog*I2", "u_inf_Hog*Q1"
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S2)"
  
  #Room 3
  #Susceptible
  "w_Hog*R3", "(beta_Hog_MDA*I3*S3)", "(beta_Hog_ind_MDA*(I4)*S3)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S3)", "(beta_WF2S*HI*S3)", "u_Hog*S3",
  #Exposed
  "sigma_Hog*E3", "u_Hog*E3",
  #Infected
  "ifelse(Q3<q_cap, (qrate*I3), 0)", "(delta_Hog*I3)", "u_inf_Hog*I3",
  #Quarantine
  "(delta_Hog*Q3)", "u_inf_Hog*Q3",
  #Recovered
  "u_Hog*R3",#"w_Hog*R3" "(delta_Hog*I3),"(delta_Hog*Q3)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E3"
  #True Recovered
  # "(delta_Hog*I3)"
  #Total Quarantined
  # "qrate*E3"
  #Natural Death Rate
  # "u_Hog*S3", "u_Hog*E3","u_Hog*R1",
  #Infected Death Rate
  # "u_inf_Hog*I3", "u_inf_Hog*Q3",
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S3)"
  
  #Room 4
  #Susceptible
  "w_Hog*R4", "(beta_Hog*I4*S4)", "(beta_Hog_ind*(0)*S4)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S4)", "(beta_WF2S*HI*S4)", "u_Hog*S4",
  #Exposed
  "sigma_Hog*E4", "u_Hog*E4",
  #Infected
  "ifelse(Q4<q_cap, (qrate*I4), 0)", "(delta_Hog*I4)", "u_inf_Hog*I4",
  #Quarantine
  "(delta_Hog*Q4)", "u_inf_Hog*Q4",
  #Recovered
  "u_Hog*R4", #"(delta_Hog*I4), "(delta_Hog*Q4)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E4"
  #True Recovered
  # "(delta_Hog*I4)"
  #Total Quarantined
  # "qrate*E4"
  #Natural Death Rate
  # "u_Hog*S4", "u_Hog*E4","u_Hog*R4",
  #Infected Death Rate
  # "u_inf_Hog*I4", "u_inf_Hog*Q4",
  #Count WF-to-Hog infections
  #"(beta_WF2S*HI*S4)"
  
  #Workforce
  #Susceptible
  "w_WF*HR", "(beta_S2WF*(I1+I2+I3+I4)*HS)", "(beta_QS2WF*(Q1+Q2+Q3+Q4)*HS)", "(beta_WF2WF*HI*HS)", "u_WF*HS",
  #Exposed
  "sigma_WF*HE", "u_WF*HE",
  #Infected
  "(delta_WF*HI)", "u_inf_WF*HI",
  #Recovered
  "u_WF*HR" #"w_Hog*R4" "(delta_hum*HI)"
  # #Natural Death Rate
  # "u_WF*HS", "u_WF*HE", "u_WF*HR",
  # #Infected Death Rate
  # "u_inf_WF*HI"
  #Count Hog-to-WF infections
  # "(beta_S2WF*(I1+I2+I3+I4)*HS)"
)
TransmissionModel_Room4MDA.ODE <- c(
  #Room 1
  #Susceptible
  "w_Hog*R1", "(beta_Hog*I1*S1)", "(beta_Hog_ind*(I2+I3+I4)*S1)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S1)", "(beta_WF2S*HI*S1)", "u_Hog*S1",
  #Exposed
  "sigma_Hog*E1", "u_Hog*E1",
  #Infected
  "ifelse(Q1<q_cap, (qrate*I1), 0)", "(delta_Hog*I1)", "u_inf_Hog*I1",
  #Quarantine
  "(delta_Hog*Q1)", "u_inf_Hog*Q1",
  #Recovered
  "u_Hog*R1", #"w_Hog*R1" "(delta_Hog*I1), "(delta_Hog*Q1)", 
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E1"
  #True Recovered
  # "(delta_Hog*I1)"
  #Total Quarantined
  # "qrate*E1"
  #Natural Death Rate
  # "u_Hog*S1", "u_Hog*E1","u_Hog*R1",
  #Infected Death Rate
  #"u_inf_Hog*I1", "u_inf_Hog*Q1"
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S1)"
  
  #Room 2
  #Susceptible
  "w_Hog*R2", "(beta_Hog*I2*S2)", "(beta_Hog_ind*(I3+I4)*S2)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S2)", "(beta_WF2S*HI*S2)", "u_Hog*S2",
  #Exposed
  "sigma_Hog*E2", "u_Hog*E2",
  #Infected
  "ifelse(Q2<q_cap, (qrate*I2), 0)", "(delta_Hog*I2)", "u_inf_Hog*I2",
  #Quarantine
  "(delta_Hog*Q2)", "u_inf_Hog*Q2",
  #Recovered
  "u_Hog*R2",# "w_Hog*R2" "(delta_Hog*I2),  "(delta_Hog*Q2)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E1"
  #True Recovered
  # "(delta_Hog*I1)"
  #Total Quarantined
  # "qrate*E2"
  #Natural Death Rate
  # "u_Hog*S2", "u_Hog*E2","u_Hog*R2",
  #Infected Death Rate
  # "u_inf_Hog*I2", "u_inf_Hog*Q1"
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S2)"
  
  #Room 3
  #Susceptible
  "w_Hog*R3", "(beta_Hog*I3*S3)", "(beta_Hog_ind*(I4)*S3)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S3)", "(beta_WF2S*HI*S3)", "u_Hog*S3",
  #Exposed
  "sigma_Hog*E3", "u_Hog*E3",
  #Infected
  "ifelse(Q3<q_cap, (qrate*I3), 0)", "(delta_Hog*I3)", "u_inf_Hog*I3",
  #Quarantine
  "(delta_Hog*Q3)", "u_inf_Hog*Q3",
  #Recovered
  "u_Hog*R3",#"w_Hog*R3" "(delta_Hog*I3),"(delta_Hog*Q3)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E3"
  #True Recovered
  # "(delta_Hog*I3)"
  #Total Quarantined
  # "qrate*E3"
  #Natural Death Rate
  # "u_Hog*S3", "u_Hog*E3","u_Hog*R1",
  #Infected Death Rate
  # "u_inf_Hog*I3", "u_inf_Hog*Q3",
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S3)"
  
  #Room 4
  #Susceptible
  "w_Hog*R4", "(beta_Hog_MDA*I4*S4)", "(beta_Hog_ind_MDA*(0)*S4)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S4)", "(beta_WF2S*HI*S4)", "u_Hog*S4",
  #Exposed
  "sigma_Hog*E4", "u_Hog*E4",
  #Infected
  "ifelse(Q4<q_cap, (qrate*I4), 0)", "(delta_Hog*I4)", "u_inf_Hog*I4",
  #Quarantine
  "(delta_Hog*Q4)", "u_inf_Hog*Q4",
  #Recovered
  "u_Hog*R4", #"(delta_Hog*I4), "(delta_Hog*Q4)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E4"
  #True Recovered
  # "(delta_Hog*I4)"
  #Total Quarantined
  # "qrate*E4"
  #Natural Death Rate
  # "u_Hog*S4", "u_Hog*E4","u_Hog*R4",
  #Infected Death Rate
  # "u_inf_Hog*I4", "u_inf_Hog*Q4",
  #Count WF-to-Hog infections
  #"(beta_WF2S*HI*S4)"
  
  #Workforce
  #Susceptible
  "w_WF*HR", "(beta_S2WF*(I1+I2+I3+I4)*HS)", "(beta_QS2WF*(Q1+Q2+Q3+Q4)*HS)", "(beta_WF2WF*HI*HS)", "u_WF*HS",
  #Exposed
  "sigma_WF*HE", "u_WF*HE",
  #Infected
  "(delta_WF*HI)", "u_inf_WF*HI",
  #Recovered
  "u_WF*HR" #"w_Hog*R4" "(delta_hum*HI)"
  # #Natural Death Rate
  # "u_WF*HS", "u_WF*HE", "u_WF*HR",
  # #Infected Death Rate
  # "u_inf_WF*HI"
  #Count Hog-to-WF infections
  # "(beta_S2WF*(I1+I2+I3+I4)*HS)"
)

TransmissionModel_Room1and2MDA.ODE <- c(
  #Room 1
  #Susceptible
  "w_Hog*R1", "(beta_Hog_MDA*I1*S1)", "(beta_Hog_ind_MDA*(I2+I3+I4)*S1)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S1)", "(beta_WF2S*HI*S1)", "u_Hog*S1",
  #Exposed
  "sigma_Hog*E1", "u_Hog*E1",
  #Infected
  "ifelse(Q1<q_cap, (qrate*I1), 0)", "(delta_Hog*I1)", "u_inf_Hog*I1",
  #Quarantine
  "(delta_Hog*Q1)", "u_inf_Hog*Q1",
  #Recovered
  "u_Hog*R1", #"w_Hog*R1" "(delta_Hog*I1),
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E1"
  #True Recovered
  # "(delta_Hog*I1)"
  #Total Quarantined
  # "qrate*E1"
  #Natural Death Rate
  # "u_Hog*S1", "u_Hog*E1","u_Hog*R1",
  #Infected Death Rate
  #"u_inf_Hog*I1",
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S1)"
  
  #Room 2
  #Susceptible
  "w_Hog*R2", "(beta_Hog_MDA*I2*S2)", "(beta_Hog_ind_MDA*(I3+I4)*S2)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S2)", "(beta_WF2S*HI*S2)", "u_Hog*S2",
  #Exposed
  "sigma_Hog*E2", "u_Hog*E2",
  #Infected
  "ifelse(Q2<q_cap, (qrate*I2), 0)", "(delta_Hog*I2)", "u_inf_Hog*I2",
  #Quarantine
  "(delta_Hog*Q2)", "u_inf_Hog*Q2",
  #Recovered
  "u_Hog*R2",# "w_Hog*R2" "(delta_Hog*I2),  "(delta_Hog*Q2)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E1"
  #True Recovered
  # "(delta_Hog*I1)"
  #Total Quarantined
  # "qrate*E2"
  #Natural Death Rate
  # "u_Hog*S2", "u_Hog*E2","u_Hog*R2",
  #Infected Death Rate
  # "u_inf_Hog*I2", "u_inf_Hog*Q1"
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S2)"
  
  #Room 3
  #Susceptible
  "w_Hog*R3", "(beta_Hog*I3*S3)", "(beta_Hog_ind*(I4)*S3)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S3)", "(beta_WF2S*HI*S3)", "u_Hog*S3",
  #Exposed
  "sigma_Hog*E3", "u_Hog*E3",
  #Infected
  "ifelse(Q3<q_cap, (qrate*I3), 0)", "(delta_Hog*I3)", "u_inf_Hog*I3",
  #Quarantine
  "(delta_Hog*Q3)", "u_inf_Hog*Q3",
  #Recovered
  "u_Hog*R3",#"w_Hog*R3" "(delta_Hog*I3),"(delta_Hog*Q3)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E3"
  #True Recovered
  # "(delta_Hog*I3)"
  #Total Quarantined
  # "qrate*E3"
  #Natural Death Rate
  # "u_Hog*S3", "u_Hog*E3","u_Hog*R1",
  #Infected Death Rate
  # "u_inf_Hog*I3", "u_inf_Hog*Q3",
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S3)"
  
  #Room 4
  #Susceptible
  "w_Hog*R4", "(beta_Hog*I4*S4)", "(beta_Hog_ind*(0)*S4)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S4)", "(beta_WF2S*HI*S4)", "u_Hog*S4",
  #Exposed
  "sigma_Hog*E4", "u_Hog*E4",
  #Infected
  "ifelse(Q4<q_cap, (qrate*I4), 0)", "(delta_Hog*I4)", "u_inf_Hog*I4",
  #Quarantine
  "(delta_Hog*Q4)", "u_inf_Hog*Q4",
  #Recovered
  "u_Hog*R4", #"(delta_Hog*I4), "(delta_Hog*Q4)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E4"
  #True Recovered
  # "(delta_Hog*I4)"
  #Total Quarantined
  # "qrate*E4"
  #Natural Death Rate
  # "u_Hog*S4", "u_Hog*E4","u_Hog*R4",
  #Infected Death Rate
  # "u_inf_Hog*I4", "u_inf_Hog*Q4",
  #Count WF-to-Hog infections
  #"(beta_WF2S*HI*S4)"
  
  #Workforce
  #Susceptible
  "w_WF*HR", "(beta_S2WF*(I1+I2+I3+I4)*HS)", "(beta_QS2WF*(Q1+Q2+Q3+Q4)*HS)", "(beta_WF2WF*HI*HS)", "u_WF*HS",
  #Exposed
  "sigma_WF*HE", "u_WF*HE",
  #Infected
  "(delta_WF*HI)", "u_inf_WF*HI",
  #Recovered
  "u_WF*HR" #"w_Hog*R4" "(delta_hum*HI)"
  # #Natural Death Rate
  # "u_WF*HS", "u_WF*HE", "u_WF*HR",
  # #Infected Death Rate
  # "u_inf_WF*HI"
  #Count Hog-to-WF infections
  # "(beta_S2WF*(I1+I2+I3+I4)*HS)"
)

TransmissionModel_Room2and3MDA.ODE <- c(
  #Room 1
  #Susceptible
  "w_Hog*R1", "(beta_Hog*I1*S1)", "(beta_Hog_ind*(I2+I3+I4)*S1)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S1)", "(beta_WF2S*HI*S1)", "u_Hog*S1",
  #Exposed
  "sigma_Hog*E1", "u_Hog*E1",
  #Infected
  "ifelse(Q1<q_cap, (qrate*I1), 0)", "(delta_Hog*I1)", "u_inf_Hog*I1",
  #Quarantine
  "(delta_Hog*Q1)", "u_inf_Hog*Q1",
  #Recovered
  "u_Hog*R1", #"w_Hog*R1" "(delta_Hog*I1), "(delta_Hog*Q1)", 
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E1"
  #True Recovered
  # "(delta_Hog*I1)"
  #Total Quarantined
  # "qrate*E1"
  #Natural Death Rate
  # "u_Hog*S1", "u_Hog*E1","u_Hog*R1",
  #Infected Death Rate
  #"u_inf_Hog*I1", "u_inf_Hog*Q1"
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S1)"
  
  #Room 2
  #Susceptible
  "w_Hog*R2", "(beta_Hog_MDA*I2*S2)", "(beta_Hog_ind_MDA*(I3+I4)*S2)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S2)", "(beta_WF2S*HI*S2)", "u_Hog*S2",
  #Exposed
  "sigma_Hog*E2", "u_Hog*E2",
  #Infected
  "ifelse(Q2<q_cap, (qrate*I2), 0)", "(delta_Hog*I2)", "u_inf_Hog*I2",
  #Quarantine
  "(delta_Hog*Q2)", "u_inf_Hog*Q2",
  #Recovered
  "u_Hog*R2",# "w_Hog*R2" "(delta_Hog*I2),  "(delta_Hog*Q2)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E1"
  #True Recovered
  # "(delta_Hog*I1)"
  #Total Quarantined
  # "qrate*E2"
  #Natural Death Rate
  # "u_Hog*S2", "u_Hog*E2","u_Hog*R2",
  #Infected Death Rate
  # "u_inf_Hog*I2", "u_inf_Hog*Q1"
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S2)"
  
  #Room 3
  #Susceptible
  "w_Hog*R3", "(beta_Hog_MDA*I3*S3)", "(beta_Hog_ind_MDA*(I4)*S3)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S3)", "(beta_WF2S*HI*S3)", "u_Hog*S3",
  #Exposed
  "sigma_Hog*E3", "u_Hog*E3",
  #Infected
  "ifelse(Q3<q_cap, (qrate*I3), 0)", "(delta_Hog*I3)", "u_inf_Hog*I3",
  #Quarantine
  "(delta_Hog*Q3)", "u_inf_Hog*Q3",
  #Recovered
  "u_Hog*R3",#"w_Hog*R3" "(delta_Hog*I3),"(delta_Hog*Q3)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E3"
  #True Recovered
  # "(delta_Hog*I3)"
  #Total Quarantined
  # "qrate*E3"
  #Natural Death Rate
  # "u_Hog*S3", "u_Hog*E3","u_Hog*R1",
  #Infected Death Rate
  # "u_inf_Hog*I3", "u_inf_Hog*Q3",
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S3)"
  
  #Room 4
  #Susceptible
  "w_Hog*R4", "(beta_Hog*I4*S4)", "(beta_Hog_ind*(0)*S4)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S4)", "(beta_WF2S*HI*S4)", "u_Hog*S4",
  #Exposed
  "sigma_Hog*E4", "u_Hog*E4",
  #Infected
  "ifelse(Q4<q_cap, (qrate*I4), 0)", "(delta_Hog*I4)", "u_inf_Hog*I4",
  #Quarantine
  "(delta_Hog*Q4)", "u_inf_Hog*Q4",
  #Recovered
  "u_Hog*R4", #"(delta_Hog*I4), "(delta_Hog*Q4)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E4"
  #True Recovered
  # "(delta_Hog*I4)"
  #Total Quarantined
  # "qrate*E4"
  #Natural Death Rate
  # "u_Hog*S4", "u_Hog*E4","u_Hog*R4",
  #Infected Death Rate
  # "u_inf_Hog*I4", "u_inf_Hog*Q4",
  #Count WF-to-Hog infections
  #"(beta_WF2S*HI*S4)"
  
  #Workforce
  #Susceptible
  "w_WF*HR", "(beta_S2WF*(I1+I2+I3+I4)*HS)", "(beta_QS2WF*(Q1+Q2+Q3+Q4)*HS)", "(beta_WF2WF*HI*HS)", "u_WF*HS",
  #Exposed
  "sigma_WF*HE", "u_WF*HE",
  #Infected
  "(delta_WF*HI)", "u_inf_WF*HI",
  #Recovered
  "u_WF*HR" #"w_Hog*R4" "(delta_hum*HI)"
  # #Natural Death Rate
  # "u_WF*HS", "u_WF*HE", "u_WF*HR",
  # #Infected Death Rate
  # "u_inf_WF*HI"
  #Count Hog-to-WF infections
  # "(beta_S2WF*(I1+I2+I3+I4)*HS)"
)

TransmissionModel_Room3and4MDA.ODE <- c(
  #Room 1
  #Susceptible
  "w_Hog*R1", "(beta_Hog*I1*S1)", "(beta_Hog_ind*(I2+I3+I4)*S1)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S1)", "(beta_WF2S*HI*S1)", "u_Hog*S1",
  #Exposed
  "sigma_Hog*E1", "u_Hog*E1",
  #Infected
  "ifelse(Q1<q_cap, (qrate*I1), 0)", "(delta_Hog*I1)", "u_inf_Hog*I1",
  #Quarantine
  "(delta_Hog*Q1)", "u_inf_Hog*Q1",
  #Recovered
  "u_Hog*R1", #"w_Hog*R1" "(delta_Hog*I1), "(delta_Hog*Q1)", 
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E1"
  #True Recovered
  # "(delta_Hog*I1)"
  #Total Quarantined
  # "qrate*E1"
  #Natural Death Rate
  # "u_Hog*S1", "u_Hog*E1","u_Hog*R1",
  #Infected Death Rate
  #"u_inf_Hog*I1", "u_inf_Hog*Q1"
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S1)"
  
  #Room 2
  #Susceptible
  "w_Hog*R2", "(beta_Hog*I2*S2)", "(beta_Hog_ind*(I3+I4)*S2)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S2)", "(beta_WF2S*HI*S2)", "u_Hog*S2",
  #Exposed
  "sigma_Hog*E2", "u_Hog*E2",
  #Infected
  "ifelse(Q2<q_cap, (qrate*I2), 0)", "(delta_Hog*I2)", "u_inf_Hog*I2",
  #Quarantine
  "(delta_Hog*Q2)", "u_inf_Hog*Q2",
  #Recovered
  "u_Hog*R2",# "w_Hog*R2" "(delta_Hog*I2),  "(delta_Hog*Q2)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E1"
  #True Recovered
  # "(delta_Hog*I1)"
  #Total Quarantined
  # "qrate*E2"
  #Natural Death Rate
  # "u_Hog*S2", "u_Hog*E2","u_Hog*R2",
  #Infected Death Rate
  # "u_inf_Hog*I2", "u_inf_Hog*Q1"
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S2)"
  
  #Room 3
  #Susceptible
  "w_Hog*R3", "(beta_Hog_MDA*I3*S3)", "(beta_Hog_ind_MDA*(I4)*S3)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S3)", "(beta_WF2S*HI*S3)", "u_Hog*S3",
  #Exposed
  "sigma_Hog*E3", "u_Hog*E3",
  #Infected
  "ifelse(Q3<q_cap, (qrate*I3), 0)", "(delta_Hog*I3)", "u_inf_Hog*I3",
  #Quarantine
  "(delta_Hog*Q3)", "u_inf_Hog*Q3",
  #Recovered
  "u_Hog*R3",#"w_Hog*R3" "(delta_Hog*I3),"(delta_Hog*Q3)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E3"
  #True Recovered
  # "(delta_Hog*I3)"
  #Total Quarantined
  # "qrate*E3"
  #Natural Death Rate
  # "u_Hog*S3", "u_Hog*E3","u_Hog*R1",
  #Infected Death Rate
  # "u_inf_Hog*I3", "u_inf_Hog*Q3",
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S3)"
  
  #Room 4
  #Susceptible
  "w_Hog*R4", "(beta_Hog_MDA*I4*S4)", "(beta_Hog_ind_MDA*(0)*S4)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S4)", "(beta_WF2S*HI*S4)", "u_Hog*S4",
  #Exposed
  "sigma_Hog*E4", "u_Hog*E4",
  #Infected
  "ifelse(Q4<q_cap, (qrate*I4), 0)", "(delta_Hog*I4)", "u_inf_Hog*I4",
  #Quarantine
  "(delta_Hog*Q4)", "u_inf_Hog*Q4",
  #Recovered
  "u_Hog*R4", #"(delta_Hog*I4), "(delta_Hog*Q4)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E4"
  #True Recovered
  # "(delta_Hog*I4)"
  #Total Quarantined
  # "qrate*E4"
  #Natural Death Rate
  # "u_Hog*S4", "u_Hog*E4","u_Hog*R4",
  #Infected Death Rate
  # "u_inf_Hog*I4", "u_inf_Hog*Q4",
  #Count WF-to-Hog infections
  #"(beta_WF2S*HI*S4)"
  
  #Workforce
  #Susceptible
  "w_WF*HR", "(beta_S2WF*(I1+I2+I3+I4)*HS)", "(beta_QS2WF*(Q1+Q2+Q3+Q4)*HS)", "(beta_WF2WF*HI*HS)", "u_WF*HS",
  #Exposed
  "sigma_WF*HE", "u_WF*HE",
  #Infected
  "(delta_WF*HI)", "u_inf_WF*HI",
  #Recovered
  "u_WF*HR" #"w_Hog*R4" "(delta_hum*HI)"
  # #Natural Death Rate
  # "u_WF*HS", "u_WF*HE", "u_WF*HR",
  # #Infected Death Rate
  # "u_inf_WF*HI"
  #Count Hog-to-WF infections
  # "(beta_S2WF*(I1+I2+I3+I4)*HS)"
)
TransmissionModel_NoMDA.ODE <- c(
  #Room 1
  #Susceptible
  "w_Hog*R1", "(beta_Hog*I1*S1)", "(beta_Hog_ind*(I2+I3+I4)*S1)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S1)", "(beta_WF2S*HI*S1)", "u_Hog*S1",
  #Exposed
  "sigma_Hog*E1", "u_Hog*E1",
  #Infected
  "ifelse(Q1<q_cap, (qrate*I1), 0)", "(delta_Hog*I1)", "u_inf_Hog*I1",
  #Quarantine
  "(delta_Hog*Q1)", "u_inf_Hog*Q1",
  #Recovered
  "u_Hog*R1", #"w_Hog*R1" "(delta_Hog*I1), "(delta_Hog*Q1)", 
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E1"
  #True Recovered
  # "(delta_Hog*I1)"
  #Total Quarantined
  # "qrate*E1"
  #Natural Death Rate
  # "u_Hog*S1", "u_Hog*E1","u_Hog*R1",
  #Infected Death Rate
  #"u_inf_Hog*I1", "u_inf_Hog*Q1"
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S1)"
  
  #Room 2
  #Susceptible
  "w_Hog*R2", "(beta_Hog*I2*S2)", "(beta_Hog_ind*(I3+I4)*S2)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S2)", "(beta_WF2S*HI*S2)", "u_Hog*S2",
  #Exposed
  "sigma_Hog*E2", "u_Hog*E2",
  #Infected
  "ifelse(Q2<q_cap, (qrate*I2), 0)", "(delta_Hog*I2)", "u_inf_Hog*I2",
  #Quarantine
  "(delta_Hog*Q2)", "u_inf_Hog*Q2",
  #Recovered
  "u_Hog*R2",# "w_Hog*R2" "(delta_Hog*I2),  "(delta_Hog*Q2)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E1"
  #True Recovered
  # "(delta_Hog*I1)"
  #Total Quarantined
  # "qrate*E2"
  #Natural Death Rate
  # "u_Hog*S2", "u_Hog*E2","u_Hog*R2",
  #Infected Death Rate
  # "u_inf_Hog*I2", "u_inf_Hog*Q1"
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S2)"
  
  #Room 3
  #Susceptible
  "w_Hog*R3", "(beta_Hog*I3*S3)", "(beta_Hog_ind*(I4)*S3)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S3)", "(beta_WF2S*HI*S3)", "u_Hog*S3",
  #Exposed
  "sigma_Hog*E3", "u_Hog*E3",
  #Infected
  "ifelse(Q3<q_cap, (qrate*I3), 0)", "(delta_Hog*I3)", "u_inf_Hog*I3",
  #Quarantine
  "(delta_Hog*Q3)", "u_inf_Hog*Q3",
  #Recovered
  "u_Hog*R3",#"w_Hog*R3" "(delta_Hog*I3),"(delta_Hog*Q3)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E3"
  #True Recovered
  # "(delta_Hog*I3)"
  #Total Quarantined
  # "qrate*E3"
  #Natural Death Rate
  # "u_Hog*S3", "u_Hog*E3","u_Hog*R1",
  #Infected Death Rate
  # "u_inf_Hog*I3", "u_inf_Hog*Q3",
  #Count WF-to-Hog  infections
  #"(beta_WF2S*HI*S3)"
  
  #Room 4
  #Susceptible
  "w_Hog*R4", "(beta_Hog*I4*S4)", "(beta_Hog_ind*(0)*S4)", "(beta_QHog_ind*(Q1+Q2+Q3+Q4)*S4)", "(beta_WF2S*HI*S4)", "u_Hog*S4",
  #Exposed
  "sigma_Hog*E4", "u_Hog*E4",
  #Infected
  "ifelse(Q4<q_cap, (qrate*I4), 0)", "(delta_Hog*I4)", "u_inf_Hog*I4",
  #Quarantine
  "(delta_Hog*Q4)", "u_inf_Hog*Q4",
  #Recovered
  "u_Hog*R4", #"(delta_Hog*I4), "(delta_Hog*Q4)",
  #Total Infected (not Vaccinated hogs)
  # "sigma_Hog*E4"
  #True Recovered
  # "(delta_Hog*I4)"
  #Total Quarantined
  # "qrate*E4"
  #Natural Death Rate
  # "u_Hog*S4", "u_Hog*E4","u_Hog*R4",
  #Infected Death Rate
  # "u_inf_Hog*I4", "u_inf_Hog*Q4",
  #Count WF-to-Hog infections
  #"(beta_WF2S*HI*S4)"
  
  #Workforce
  #Susceptible
  "w_WF*HR", "(beta_S2WF*(I1+I2+I3+I4)*HS)", "(beta_QS2WF*(Q1+Q2+Q3+Q4)*HS)", "(beta_WF2WF*HI*HS)", "u_WF*HS",
  #Exposed
  "sigma_WF*HE", "u_WF*HE",
  #Infected
  "(delta_WF*HI)", "u_inf_WF*HI",
  #Recovered
  "u_WF*HR" #"w_Hog*R4" "(delta_hum*HI)"
  # #Natural Death Rate
  # "u_WF*HS", "u_WF*HE", "u_WF*HR",
  # #Infected Death Rate
  # "u_inf_WF*HI"
  #Count Hog-to-WF infections
  # "(beta_S2WF*(I1+I2+I3+I4)*HS)"
)

TransmissionModel.Matrix<- matrix(c(
  #Room 1
  #Susceptible
  1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 10),
  #Exposed
  0, 1, 1, 1, 1, 0, -1, -1, 0, 0, 0, 0, 0, 0, rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 10),
  #Infected
  0, 0, 0, 0, 0, 0, 1, 0, -1, -1, -1, 0, 0, 0, rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 10),
  #Quarantine
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, -1, 0, rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 10),
  #Recovered
  -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -1, rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 10),
  #Total Infected (not Vaccinated hogs)
  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 10),
  #True Recovered
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 10),
  #Total Quarantined
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 10), 
  #Natural Death Rate
  0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 10),
  #Infected Death Rate
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 10),
  #Count WF-to-Hog infections
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 10),
  #Room 2
  #Susceptible
  rep(0, 14), 1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, rep(0, 14), rep(0, 14), rep(0, 10),
  #Exposed
  rep(0, 14), 0, 1, 1, 1, 1, 0, -1, -1, 0, 0, 0, 0, 0, 0, rep(0, 14), rep(0, 14), rep(0, 10),
  #Infected
  rep(0, 14), 0, 0, 0, 0, 0, 0, 1, 0, -1, -1, -1, 0, 0, 0, rep(0, 14), rep(0, 14), rep(0, 10),
  #Quarantine
  rep(0, 14), 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, -1, 0, rep(0, 14), rep(0, 14), rep(0, 10),
  #Recovered
  rep(0, 14), -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -1, rep(0, 14), rep(0, 14), rep(0, 10),
  #Total Infected (not Vaccinated hogs)
  rep(0, 14), 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, rep(0, 14), rep(0, 14), rep(0, 10),
  #True Recovered
  rep(0, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, rep(0, 14), rep(0, 14), rep(0, 10),
  #Total Quarantined
  rep(0, 14), 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, rep(0, 14), rep(0, 14), rep(0, 10), 
  #Natural Death Rate
  rep(0, 14), 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, rep(0, 14), rep(0, 14), rep(0, 10),
  #Infected Death Rate
  rep(0, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, rep(0, 14), rep(0, 14), rep(0, 10),
  #Count WF-to-Hog infections
  rep(0, 14), 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, rep(0, 14), rep(0, 14), rep(0, 10),
  #Room 3
  #Susceptible
  rep(0, 14), rep(0, 14), 1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, rep(0, 14), rep(0, 10),
  #Exposed
  rep(0, 14), rep(0, 14), 0, 1, 1, 1, 1, 0, -1, -1, 0, 0, 0, 0, 0, 0, rep(0, 14), rep(0, 10),
  #Infected
  rep(0, 14), rep(0, 14), 0, 0, 0, 0, 0, 0, 1, 0, -1, -1, -1, 0, 0, 0, rep(0, 14), rep(0, 10),
  #Quarantine
  rep(0, 14), rep(0, 14), 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, -1, 0, rep(0, 14), rep(0, 10),
  #Recovered
  rep(0, 14), rep(0, 14), -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -1, rep(0, 14), rep(0, 10),
  #Total Infected (not Vaccinated hogs)
  rep(0, 14), rep(0, 14), 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, rep(0, 14), rep(0, 10),
  #True Recovered
  rep(0, 14), rep(0, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, rep(0, 14), rep(0, 10),
  #Total Quarantined
  rep(0, 14), rep(0, 14), 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, rep(0, 14), rep(0, 10), 
  #Natural Death Rate
  rep(0, 14), rep(0, 14), 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, rep(0, 14), rep(0, 10),
  #Infected Death Rate
  rep(0, 14), rep(0, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, rep(0, 14), rep(0, 10),
  #Count WF-to-Hog infections
  rep(0, 14), rep(0, 14), 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, rep(0, 14), rep(0, 10),
  #Room 4
  #Susceptible
  rep(0, 14), rep(0, 14), rep(0, 14), 1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, rep(0, 10),
  #Exposed
  rep(0, 14), rep(0, 14), rep(0, 14), 0, 1, 1, 1, 1, 0, -1, -1, 0, 0, 0, 0, 0, 0, rep(0, 10),
  #Infected
  rep(0, 14), rep(0, 14), rep(0, 14), 0, 0, 0, 0, 0, 0, 1, 0, -1, -1, -1, 0, 0, 0, rep(0, 10),
  #Quarantine
  rep(0, 14), rep(0, 14), rep(0, 14), 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, -1, 0, rep(0, 10),
  #Recovered
  rep(0, 14), rep(0, 14), rep(0, 14), -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -1, rep(0, 10),
  #Total Infected (not Vaccinated hogs)
  rep(0, 14), rep(0, 14), rep(0, 14), 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, rep(0, 10),
  #True Recovered
  rep(0, 14), rep(0, 14), rep(0, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, rep(0, 10),
  #Total Quarantined
  rep(0, 14), rep(0, 14), rep(0, 14), 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, rep(0, 10), 
  #Natural Death Rate
  rep(0, 14), rep(0, 14), rep(0, 14), 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, rep(0, 10),
  #Infected Death Rate
  rep(0, 14), rep(0, 14), rep(0, 14), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, rep(0, 10),
  #Count WF-to-Hog infections
  rep(0, 14), rep(0, 14), rep(0, 14), 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, rep(0, 10),
  #Workforce 
  #Susceptible
  rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 14), 1, -1, -1, -1, -1, 0, 0, 0, 0, 0,
  #Exposed
  rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 14), 0, 1, 1, 1, 0, -1, -1, 0, 0, 0,
  #Infected
  rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 14), 0, 0, 0, 0, 0, 1, 0, -1, -1, 0,
  #Recovered
  rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 14), -1, 0, 0, 0, 0, 0, 0, 1, 0, -1,
  #Natural Death Rate
  rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 14), 0, 0, 0, 0, 1, 0, 1, 0, 0, 1,
  #Infected Death Rate
  rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 14), 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
  #Count Hog-to-WF infections
  rep(0, 14), rep(0, 14), rep(0, 14), rep(0, 14), 0, 1, 1, 0, 0, 0, 0, 0, 0, 0
),
nrow = 51, byrow = TRUE)
# write.csv(TransmissionModel.Matrix, file = "v2_TransmissionModelTransistionMatrix.csv", row.names = F)
#### MODEL SETUP ####
#Time to model end (in days)
OneWeek.tf <- (1 * 7)  #Run for 23 weeks (till hogs reach market weight)
TwoWeek.tf <- (2 * 7)
FourteenWeek.tf <- (14 * 7)

#Model Run
results <- NULL
n.itr <- 5000
# seed.num <- 616
seed.nums <- c(100:(100+n.itr))
for (i in c(1:n.itr)) {
  print(i)
  tmp1 <- NA
  tmp2 <- NA
  tmp3 <- NA
  tmp4 <- NA
  tmp5 <- NA
  tmp6 <- NA
  tmp7 <- NA
  tmp8 <- NA
  tmp9 <- NA
  tmp10 <- NA
  tmp11 <- NA
  tmp12 <- NA
  ### WEEKS 0-2 Room 1 fill | Room 1 with MDA
  #Initials
  TransmissionModel.Int.0to2 <- c(
    #Room1
    S1 = round((1-Hog_Vacc_Eff )*1000, digits = 0),
    E1 = 0,
    I1 = 0,
    Q1 = 0,
    R1 = round(Hog_Vacc_Eff *1000, digits = 0),
    Total_inf1 = 0,
    True_Rec1 = 0,
    Total_Q1 = 0,
    D1_Nat = 0,
    D1_Inf = 0,
    WF2Hog_Count1 = 0,
    
    #Room2
    S2 = 0,
    E2 = 0,
    I2 = 0,
    Q2 = 0,
    R2 = 0,
    Total_inf2 = 0,
    True_Rec2 = 0,
    Total_Q2 = 0,
    D2_Nat = 0,
    D2_Inf = 0,
    WF2Hog_Count2 = 0,
    
    #Room3
    S3 = 0,
    E3 = 0,
    I3 = 0,
    Q3 = 0,
    R3 = 0,
    Total_inf3 = 0,
    True_Rec3 = 0,
    Total_Q3 = 0,
    D3_Nat = 0,
    D3_Inf = 0,
    WF2Hog_Count3 = 0,
    
    #Room4
    S4 = 0,
    E4 = 0,
    I4 = 0,
    Q4 = 0,
    R4 = 0,
    Total_inf4 = 0,
    True_Rec4 = 0,
    Total_Q4 = 0,
    D4_Nat = 0,
    D4_Inf = 0,
    WF2Hog_Count4 = 0,
    
    #Workforce
    HS = round((1-Hum_Vacc_Eff )*2, digits = 0), #need to add a +1 for 5 workers b/c round function drops one 
    HE = 0, 
    HI = 0, 
    HR = round(Hum_Vacc_Eff *2, digits = 0), # +1 for 5 workers to either susceptible or recovered
    H_DN = 0, 
    H_DI = 0,
    Hog2WF_Count = 0
  )
  set.seed(seed.nums[i])
  FullRooms.0to2 <- ssa(
    x0 = TransmissionModel.Int.0to2,
    a = TransmissionModel_Room1MDA.ODE,
    nu = TransmissionModel.Matrix,
    parms = TransmissionModel.par,
    tf = TwoWeek.tf,
    method = ssa.d(),
    verbose = FALSE,
    censusInterval = 0.1,
    consoleInterval = 1
  ) 
  tmp1 <- as.data.frame(FullRooms.0to2$data)
  tmp1 <- tmp1[tmp1$t < 14,]
  
  ### WEEKS 2-3 Room 2 fill | Room 1 and Room 2 with MDAs
  #Initials
  TransmissionModel.Int.2to3 <- c(
    #Room1
    S1 = tmp1$S1[nrow(tmp1)],
    E1 = tmp1$E1[nrow(tmp1)],
    I1 = tmp1$I1[nrow(tmp1)],
    Q1 = tmp1$Q1[nrow(tmp1)],
    R1 = tmp1$R1[nrow(tmp1)],
    Total_inf1 = tmp1$Total_inf1[nrow(tmp1)],
    True_Rec1 = tmp1$True_Rec1[nrow(tmp1)],
    Total_Q1 = tmp1$Total_Q1[nrow(tmp1)],
    D1_Nat = tmp1$D1_Nat[nrow(tmp1)],
    D1_Inf = tmp1$D1_Inf[nrow(tmp1)],
    WF2Hog_Count1 = tmp1$WF2Hog_Count1[nrow(tmp1)],
    
    #Room2
    S2 = round((1-Hog_Vacc_Eff )*1000, digits = 0),
    E2 = 0,
    I2 = 0,
    Q2 = 0,
    R2 = round(Hog_Vacc_Eff *1000, digits = 0),
    Total_inf2 = 0,
    True_Rec2 = 0,
    Total_Q2 = 0,
    D2_Nat = 0,
    D2_Inf = 0,
    WF2Hog_Count2 = 0,
    
    #Room3
    S3 = 0,
    E3 = 0,
    I3 = 0,
    Q3 = 0,
    R3 = 0,
    Total_inf3 = 0,
    True_Rec3 = 0,
    Total_Q3 = 0,
    D3_Nat = 0,
    D3_Inf = 0,
    WF2Hog_Count3 = 0,
    
    #Room4
    S4 = 0,
    E4 = 0,
    I4 = 0,
    Q4 = 0,
    R4 = 0,
    Total_inf4 = 0,
    True_Rec4 = 0,
    Total_Q4 = 0,
    D4_Nat = 0,
    D4_Inf = 0,
    WF2Hog_Count4 = 0,
    
    #Workforce
    HS = tmp1$HS[nrow(tmp1)], 
    HE = tmp1$HE[nrow(tmp1)], 
    HI = tmp1$HI[nrow(tmp1)], 
    HR = tmp1$HR[nrow(tmp1)], 
    H_DN = tmp1$H_DN[nrow(tmp1)], 
    H_DI = tmp1$H_DI[nrow(tmp1)],
    Hog2WF_Count = tmp1$Hog2WF_Count[nrow(tmp1)]
  )
  set.seed(seed.nums[i])
  FullRooms.2to3 <- ssa(
    x0 = TransmissionModel.Int.2to3,
    a = TransmissionModel_Room1and2MDA.ODE,
    nu = TransmissionModel.Matrix,
    parms = TransmissionModel.par,
    tf = OneWeek.tf,
    method = ssa.d(),
    verbose = FALSE,
    censusInterval = 0.1,
    consoleInterval = 1
  ) 
  tmp2 <- as.data.frame(FullRooms.2to3$data)
  tmp2$t <- tmp2$t +14
  tmp2 <- tmp2[tmp2$t < 21,]
  
  ### WEEKS 3-4 Room 1 No MDA | Room 2 with MDAs
  #Initials
  TransmissionModel.Int.3to4 <- c(
    #Room1
    S1 = tmp2$S1[nrow(tmp2)],
    E1 = tmp2$E1[nrow(tmp2)],
    I1 = tmp2$I1[nrow(tmp2)],
    Q1 = tmp2$Q1[nrow(tmp2)],
    R1 = tmp2$R1[nrow(tmp2)],
    Total_inf1 = tmp2$Total_inf1[nrow(tmp2)],
    True_Rec1 = tmp2$True_Rec1[nrow(tmp2)],
    Total_Q1 = tmp2$Total_Q1[nrow(tmp2)],
    D1_Nat = tmp2$D1_Nat[nrow(tmp2)],
    D1_Inf = tmp2$D1_Inf[nrow(tmp2)],
    WF2Hog_Count1 = tmp2$WF2Hog_Count1[nrow(tmp2)],
    
    #Room2
    S2 = tmp2$S2[nrow(tmp2)],
    E2 = tmp2$E2[nrow(tmp2)],
    I2 = tmp2$I2[nrow(tmp2)],
    Q2 = tmp2$Q2[nrow(tmp2)],
    R2 = tmp2$R2[nrow(tmp2)],
    Total_inf2 = tmp2$Total_inf2[nrow(tmp2)],
    True_Rec2 = tmp2$True_Rec2[nrow(tmp2)],
    Total_Q2 = tmp2$Total_Q2[nrow(tmp2)],
    D2_Nat = tmp2$D2_Nat[nrow(tmp2)],
    D2_Inf = tmp2$D2_Inf[nrow(tmp2)],
    WF2Hog_Count2 = tmp2$WF2Hog_Count2[nrow(tmp2)],
    
    #Room3
    S3 = 0,
    E3 = 0,
    I3 = 0,
    Q3 = 0,
    R3 = 0,
    Total_inf3 = 0,
    True_Rec3 = 0,
    Total_Q3 = 0,
    D3_Nat = 0,
    D3_Inf = 0,
    WF2Hog_Count3 = 0,
    
    #Room4
    S4 = 0,
    E4 = 0,
    I4 = 0,
    Q4 = 0,
    R4 = 0,
    Total_inf4 = 0,
    True_Rec4 = 0,
    Total_Q4 = 0,
    D4_Nat = 0,
    D4_Inf = 0,
    WF2Hog_Count4 = 0,
    
    #Workforce
    #Workforce
    HS = tmp2$HS[nrow(tmp2)], 
    HE = tmp2$HE[nrow(tmp2)], 
    HI = tmp2$HI[nrow(tmp2)], 
    HR = tmp2$HR[nrow(tmp2)], 
    H_DN = tmp2$H_DN[nrow(tmp2)], 
    H_DI = tmp2$H_DI[nrow(tmp2)],
    Hog2WF_Count = tmp2$Hog2WF_Count[nrow(tmp2)]
  )
  set.seed(seed.nums[i])
  FullRooms.3to4 <- ssa(
    x0 = TransmissionModel.Int.3to4,
    a = TransmissionModel_Room2MDA.ODE,
    nu = TransmissionModel.Matrix,
    parms = TransmissionModel.par,
    tf = OneWeek.tf,
    method = ssa.d(),
    verbose = FALSE,
    censusInterval = 0.1,
    consoleInterval = 1
  ) 
  tmp3 <- as.data.frame(FullRooms.3to4$data)
  tmp3$t <- tmp3$t +21
  tmp3 <- tmp3[tmp3$t < 28,]
  
  ### WEEKS 4-5 Room 3 Fill | Room 2 and Room 3 with MDA | Room 1 no MDA
  #Initials
  TransmissionModel.Int.4to5 <- c(
    #Room1
    S1 = tmp3$S1[nrow(tmp3)],
    E1 = tmp3$E1[nrow(tmp3)],
    I1 = tmp3$I1[nrow(tmp3)],
    Q1 = tmp3$Q1[nrow(tmp3)],
    R1 = tmp3$R1[nrow(tmp3)],
    Total_inf1 = tmp3$Total_inf1[nrow(tmp3)],
    True_Rec1 = tmp3$True_Rec1[nrow(tmp3)],
    Total_Q1 = tmp3$Total_Q1[nrow(tmp3)],
    D1_Nat = tmp3$D1_Nat[nrow(tmp3)],
    D1_Inf = tmp3$D1_Inf[nrow(tmp3)],
    WF2Hog_Count1 = tmp3$WF2Hog_Count1[nrow(tmp3)],
    
    #Room2
    S2 = tmp3$S2[nrow(tmp3)],
    E2 = tmp3$E2[nrow(tmp3)],
    I2 = tmp3$I2[nrow(tmp3)],
    Q2 = tmp3$Q2[nrow(tmp3)],
    R2 = tmp3$R2[nrow(tmp3)],
    Total_inf2 = tmp3$Total_inf2[nrow(tmp3)],
    True_Rec2 = tmp3$True_Rec2[nrow(tmp3)],
    Total_Q2 = tmp3$Total_Q2[nrow(tmp3)],
    D2_Nat = tmp3$D2_Nat[nrow(tmp3)],
    D2_Inf = tmp3$D2_Inf[nrow(tmp3)],
    WF2Hog_Count2 = tmp3$WF2Hog_Count2[nrow(tmp3)],
    
    #Room3
    S3 = round((1-Hog_Vacc_Eff )*1000, digits = 0),
    E3 = 0,
    I3 = 0,
    Q3 = 0,
    R3 = round(Hog_Vacc_Eff *1000, digits = 0),
    Total_inf3 = 0,
    True_Rec3 = 0,
    Total_Q3 = 0,
    D3_Nat = 0,
    D3_Inf = 0,
    WF2Hog_Count3 = 0,
    
    #Room4
    S4 = 0,
    E4 = 0,
    I4 = 0,
    Q4 = 0,
    R4 = 0,
    Total_inf4 = 0,
    True_Rec4 = 0,
    Total_Q4 = 0,
    D4_Nat = 0,
    D4_Inf = 0,
    WF2Hog_Count4 = 0,
    
    #Workforce
    HS = tmp3$HS[nrow(tmp3)], 
    HE = tmp3$HE[nrow(tmp3)], 
    HI = tmp3$HI[nrow(tmp3)], 
    HR = tmp3$HR[nrow(tmp3)], 
    H_DN = tmp3$H_DN[nrow(tmp3)], 
    H_DI = tmp3$H_DI[nrow(tmp3)],
    Hog2WF_Count = tmp3$Hog2WF_Count[nrow(tmp3)]
  )
  set.seed(seed.nums[i])
  FullRooms.4to5 <- ssa(
    x0 = TransmissionModel.Int.4to5,
    a = TransmissionModel_Room2and3MDA.ODE,
    nu = TransmissionModel.Matrix,
    parms = TransmissionModel.par,
    tf = OneWeek.tf,
    method = ssa.d(),
    verbose = FALSE,
    censusInterval = 0.1,
    consoleInterval = 1
  ) 
  tmp4 <- as.data.frame(FullRooms.4to5$data)
  tmp4$t <- tmp4$t +28
  tmp4 <- tmp4[tmp4$t < 35,] 
  
  ### WEEKS 5-6 Room 3 with MDA | Room 1 and Room 2 no MDA
  #Initials
  TransmissionModel.Int.5to6 <- c(
    #Room1
    S1 = tmp4$S1[nrow(tmp4)],
    E1 = tmp4$E1[nrow(tmp4)],
    I1 = tmp4$I1[nrow(tmp4)],
    Q1 = tmp4$Q1[nrow(tmp4)],
    R1 = tmp4$R1[nrow(tmp4)],
    Total_inf1 = tmp4$Total_inf1[nrow(tmp4)],
    True_Rec1 = tmp4$True_Rec1[nrow(tmp4)],
    Total_Q1 = tmp4$Total_Q1[nrow(tmp4)],
    D1_Nat = tmp4$D1_Nat[nrow(tmp4)],
    D1_Inf = tmp4$D1_Inf[nrow(tmp4)],
    WF2Hog_Count1 = tmp4$WF2Hog_Count1[nrow(tmp4)],
    
    #Room2
    S2 = tmp4$S2[nrow(tmp4)],
    E2 = tmp4$E2[nrow(tmp4)],
    I2 = tmp4$I2[nrow(tmp4)],
    Q2 = tmp4$Q2[nrow(tmp4)],
    R2 = tmp4$R2[nrow(tmp4)],
    Total_inf2 = tmp4$Total_inf2[nrow(tmp4)],
    True_Rec2 = tmp4$True_Rec2[nrow(tmp4)],
    Total_Q2 = tmp4$Total_Q2[nrow(tmp4)],
    D2_Nat = tmp4$D2_Nat[nrow(tmp4)],
    D2_Inf = tmp4$D2_Inf[nrow(tmp4)],
    WF2Hog_Count2 = tmp4$WF2Hog_Count2[nrow(tmp4)],
    
    #Room3
    S3 = tmp4$S3[nrow(tmp4)],
    E3 = tmp4$E3[nrow(tmp4)],
    I3 = tmp4$I3[nrow(tmp4)],
    Q3 = tmp4$Q3[nrow(tmp4)],
    R3 = tmp4$R3[nrow(tmp4)],
    Total_inf3 = tmp4$Total_inf3[nrow(tmp4)],
    True_Rec3 = tmp4$True_Rec3[nrow(tmp4)],
    Total_Q3 = tmp4$Total_Q3[nrow(tmp4)],
    D3_Nat = tmp4$D3_Nat[nrow(tmp4)],
    D3_Inf = tmp4$D3_Inf[nrow(tmp4)],
    WF2Hog_Count3 = tmp4$WF2Hog_Count3[nrow(tmp4)],
    
    #Room4
    S4 = 0,
    E4 = 0,
    I4 = 0,
    Q4 = 0,
    R4 = 0,
    Total_inf4 = 0,
    True_Rec4 = 0,
    Total_Q4 = 0,
    D4_Nat = 0,
    D4_Inf = 0,
    WF2Hog_Count4 = 0,
    
    #Workforce
    HS = tmp4$HS[nrow(tmp4)], 
    HE = tmp4$HE[nrow(tmp4)], 
    HI = tmp4$HI[nrow(tmp4)], 
    HR = tmp4$HR[nrow(tmp4)], 
    H_DN = tmp4$H_DN[nrow(tmp4)], 
    H_DI = tmp4$H_DI[nrow(tmp4)],
    Hog2WF_Count = tmp4$Hog2WF_Count[nrow(tmp4)]
  )
  set.seed(seed.nums[i])
  FullRooms.5to6 <- ssa(
    x0 = TransmissionModel.Int.5to6,
    a = TransmissionModel_Room3MDA.ODE,
    nu = TransmissionModel.Matrix,
    parms = TransmissionModel.par,
    tf = OneWeek.tf,
    method = ssa.d(),
    verbose = FALSE,
    censusInterval = 0.1,
    consoleInterval = 1
  ) 
  tmp5 <- as.data.frame(FullRooms.5to6$data)
  tmp5$t <- tmp5$t +35
  tmp5 <- tmp5[tmp5$t < 42,]  
  
  ### WEEKS 6-7 Room 4 Fill | Room 3 and Room 4 with MDA | Room 1 and Room 2 no MDA
  #Initials
  TransmissionModel.Int.6to7 <- c(
    #Room1
    S1 = tmp5$S1[nrow(tmp5)],
    E1 = tmp5$E1[nrow(tmp5)],
    I1 = tmp5$I1[nrow(tmp5)],
    Q1 = tmp5$Q1[nrow(tmp5)],
    R1 = tmp5$R1[nrow(tmp5)],
    Total_inf1 = tmp5$Total_inf1[nrow(tmp5)],
    True_Rec1 = tmp5$True_Rec1[nrow(tmp5)],
    Total_Q1 = tmp5$Total_Q1[nrow(tmp5)],
    D1_Nat = tmp5$D1_Nat[nrow(tmp5)],
    D1_Inf = tmp5$D1_Inf[nrow(tmp5)],
    WF2Hog_Count1 = tmp5$WF2Hog_Count1[nrow(tmp5)],
    
    #Room2
    S2 = tmp5$S2[nrow(tmp5)],
    E2 = tmp5$E2[nrow(tmp5)],
    I2 = tmp5$I2[nrow(tmp5)],
    Q2 = tmp5$Q2[nrow(tmp5)],
    R2 = tmp5$R2[nrow(tmp5)],
    Total_inf2 = tmp5$Total_inf2[nrow(tmp5)],
    True_Rec2 = tmp5$True_Rec2[nrow(tmp5)],
    Total_Q2 = tmp5$Total_Q2[nrow(tmp5)],
    D2_Nat = tmp5$D2_Nat[nrow(tmp5)],
    D2_Inf = tmp5$D2_Inf[nrow(tmp5)],
    WF2Hog_Count2 = tmp5$WF2Hog_Count2[nrow(tmp5)],
    
    #Room3
    S3 = tmp5$S3[nrow(tmp5)],
    E3 = tmp5$E3[nrow(tmp5)],
    I3 = tmp5$I3[nrow(tmp5)],
    Q3 = tmp5$Q3[nrow(tmp5)],
    R3 = tmp5$R3[nrow(tmp5)],
    Total_inf3 = tmp5$Total_inf3[nrow(tmp5)],
    True_Rec3 = tmp5$True_Rec3[nrow(tmp5)],
    Total_Q3 = tmp5$Total_Q3[nrow(tmp5)],
    D3_Nat = tmp5$D3_Nat[nrow(tmp5)],
    D3_Inf = tmp5$D3_Inf[nrow(tmp5)],
    WF2Hog_Count3 = tmp5$WF2Hog_Count3[nrow(tmp5)],
    
    #Room4
    S4 = round((1-Hog_Vacc_Eff )*999, digits = 0),
    E4 = 0,
    I4 = 1,
    Q4 = 0,
    R4 = round(Hog_Vacc_Eff *999, digits = 0),
    Total_inf4 = 0,
    True_Rec4 = 0,
    Total_Q4 = 0,
    D4_Nat = 0,
    D4_Inf = 0,
    WF2Hog_Count4 = 0,
    
    #Workforce
    HS = tmp5$HS[nrow(tmp5)], 
    HE = tmp5$HE[nrow(tmp5)], 
    HI = tmp5$HI[nrow(tmp5)], 
    HR = tmp5$HR[nrow(tmp5)], 
    H_DN = tmp5$H_DN[nrow(tmp5)], 
    H_DI = tmp5$H_DI[nrow(tmp5)],
    Hog2WF_Count = tmp5$Hog2WF_Count[nrow(tmp5)]
  )
  set.seed(seed.nums[i])
  FullRooms.6to7 <- ssa(
    x0 = TransmissionModel.Int.6to7,
    a = TransmissionModel_Room3and4MDA.ODE,
    nu = TransmissionModel.Matrix,
    parms = TransmissionModel.par,
    tf = OneWeek.tf,
    method = ssa.d(),
    verbose = FALSE,
    censusInterval = 0.1,
    consoleInterval = 1
  ) 
  tmp6 <- as.data.frame(FullRooms.6to7$data)
  tmp6$t <- tmp6$t +42
  tmp6 <- tmp6[tmp6$t < 49,]   
  
  ### WEEKS 7-89 Room 4 with MDA | Room 1, Room 2, and Room 3 no MDA
  #Initials
  TransmissionModel.Int.7to9 <- c(
    #Room1
    S1 = tmp6$S1[nrow(tmp6)],
    E1 = tmp6$E1[nrow(tmp6)],
    I1 = tmp6$I1[nrow(tmp6)],
    Q1 = tmp6$Q1[nrow(tmp6)],
    R1 = tmp6$R1[nrow(tmp6)],
    Total_inf1 = tmp6$Total_inf1[nrow(tmp6)],
    True_Rec1 = tmp6$True_Rec1[nrow(tmp6)],
    Total_Q1 = tmp6$Total_Q1[nrow(tmp6)],
    D1_Nat = tmp6$D1_Nat[nrow(tmp6)],
    D1_Inf = tmp6$D1_Inf[nrow(tmp6)],
    WF2Hog_Count1 = tmp6$WF2Hog_Count1[nrow(tmp6)],
    
    #Room2
    S2 = tmp6$S2[nrow(tmp6)],
    E2 = tmp6$E2[nrow(tmp6)],
    I2 = tmp6$I2[nrow(tmp6)],
    Q2 = tmp6$Q2[nrow(tmp6)],
    R2 = tmp6$R2[nrow(tmp6)],
    Total_inf2 = tmp6$Total_inf2[nrow(tmp6)],
    True_Rec2 = tmp6$True_Rec2[nrow(tmp6)],
    Total_Q2 = tmp6$Total_Q2[nrow(tmp6)],
    D2_Nat = tmp6$D2_Nat[nrow(tmp6)],
    D2_Inf = tmp6$D2_Inf[nrow(tmp6)],
    WF2Hog_Count2 = tmp6$WF2Hog_Count2[nrow(tmp6)],
    
    #Room3
    S3 = tmp6$S3[nrow(tmp6)],
    E3 = tmp6$E3[nrow(tmp6)],
    I3 = tmp6$I3[nrow(tmp6)],
    Q3 = tmp6$Q3[nrow(tmp6)],
    R3 = tmp6$R3[nrow(tmp6)],
    Total_inf3 = tmp6$Total_inf3[nrow(tmp6)],
    True_Rec3 = tmp6$True_Rec3[nrow(tmp6)],
    Total_Q3 = tmp6$Total_Q3[nrow(tmp6)],
    D3_Nat = tmp6$D3_Nat[nrow(tmp6)],
    D3_Inf = tmp6$D3_Inf[nrow(tmp6)],
    WF2Hog_Count3 = tmp6$WF2Hog_Count3[nrow(tmp6)],
    
    #Room4
    S4 = tmp6$S4[nrow(tmp6)],
    E4 = tmp6$E4[nrow(tmp6)],
    I4 = tmp6$I4[nrow(tmp6)],
    Q4 = tmp6$Q4[nrow(tmp6)],
    R4 = tmp6$R4[nrow(tmp6)],
    Total_inf4 = tmp6$Total_inf4[nrow(tmp6)],
    True_Rec4 = tmp6$True_Rec4[nrow(tmp6)],
    Total_Q4 = tmp6$Total_Q4[nrow(tmp6)],
    D4_Nat = tmp6$D4_Nat[nrow(tmp6)],
    D4_Inf = tmp6$D4_Inf[nrow(tmp6)],
    WF2Hog_Count4 = tmp6$WF2Hog_Count4[nrow(tmp6)],
    
    #Workforce
    HS = tmp6$HS[nrow(tmp6)], 
    HE = tmp6$HE[nrow(tmp6)], 
    HI = tmp6$HI[nrow(tmp6)], 
    HR = tmp6$HR[nrow(tmp6)], 
    H_DN = tmp6$H_DN[nrow(tmp6)], 
    H_DI = tmp6$H_DI[nrow(tmp6)],
    Hog2WF_Count = tmp6$Hog2WF_Count[nrow(tmp6)]
  )
  set.seed(seed.nums[i])
  FullRooms.7to9 <- ssa(
    x0 = TransmissionModel.Int.7to9,
    a = TransmissionModel_Room4MDA.ODE,
    nu = TransmissionModel.Matrix,
    parms = TransmissionModel.par,
    tf = TwoWeek.tf,
    method = ssa.d(),
    verbose = FALSE,
    censusInterval = 0.1,
    consoleInterval = 1
  ) 
  tmp7 <- as.data.frame(FullRooms.7to9$data)
  tmp7$t <- tmp7$t +49
  tmp7 <- tmp7[tmp7$t < 63,]    
  
  ### WEEKS 9-23 no MDA in any rooms
  #Initials
  TransmissionModel.Int.9to23 <- c(
    #Room1
    S1 = tmp7$S1[nrow(tmp7)],
    E1 = tmp7$E1[nrow(tmp7)],
    I1 = tmp7$I1[nrow(tmp7)],
    Q1 = tmp7$Q1[nrow(tmp7)],
    R1 = tmp7$R1[nrow(tmp7)],
    Total_inf1 = tmp7$Total_inf1[nrow(tmp7)],
    True_Rec1 = tmp7$True_Rec1[nrow(tmp7)],
    Total_Q1 = tmp7$Total_Q1[nrow(tmp7)],
    D1_Nat = tmp7$D1_Nat[nrow(tmp7)],
    D1_Inf = tmp7$D1_Inf[nrow(tmp7)],
    WF2Hog_Count1 = tmp7$WF2Hog_Count1[nrow(tmp7)],
    
    #Room2
    S2 = tmp7$S2[nrow(tmp7)],
    E2 = tmp7$E2[nrow(tmp7)],
    I2 = tmp7$I2[nrow(tmp7)],
    Q2 = tmp7$Q2[nrow(tmp7)],
    R2 = tmp7$R2[nrow(tmp7)],
    Total_inf2 = tmp7$Total_inf2[nrow(tmp7)],
    True_Rec2 = tmp7$True_Rec2[nrow(tmp7)],
    Total_Q2 = tmp7$Total_Q2[nrow(tmp7)],
    D2_Nat = tmp7$D2_Nat[nrow(tmp7)],
    D2_Inf = tmp7$D2_Inf[nrow(tmp7)],
    WF2Hog_Count2 = tmp7$WF2Hog_Count2[nrow(tmp7)],
    
    #Room3
    S3 = tmp7$S3[nrow(tmp7)],
    E3 = tmp7$E3[nrow(tmp7)],
    I3 = tmp7$I3[nrow(tmp7)],
    Q3 = tmp7$Q3[nrow(tmp7)],
    R3 = tmp7$R3[nrow(tmp7)],
    Total_inf3 = tmp7$Total_inf3[nrow(tmp7)],
    True_Rec3 = tmp7$True_Rec3[nrow(tmp7)],
    Total_Q3 = tmp7$Total_Q3[nrow(tmp7)],
    D3_Nat = tmp7$D3_Nat[nrow(tmp7)],
    D3_Inf = tmp7$D3_Inf[nrow(tmp7)],
    WF2Hog_Count3 = tmp7$WF2Hog_Count3[nrow(tmp7)],
    
    #Room4
    S4 = tmp7$S4[nrow(tmp7)],
    E4 = tmp7$E4[nrow(tmp7)],
    I4 = tmp7$I4[nrow(tmp7)],
    Q4 = tmp7$Q4[nrow(tmp7)],
    R4 = tmp7$R4[nrow(tmp7)],
    Total_inf4 = tmp7$Total_inf4[nrow(tmp7)],
    True_Rec4 = tmp7$True_Rec4[nrow(tmp7)],
    Total_Q4 = tmp7$Total_Q4[nrow(tmp7)],
    D4_Nat = tmp7$D4_Nat[nrow(tmp7)],
    D4_Inf = tmp7$D4_Inf[nrow(tmp7)],
    WF2Hog_Count4 = tmp7$WF2Hog_Count4[nrow(tmp7)],
    
    #Workforce
    HS = tmp7$HS[nrow(tmp7)], 
    HE = tmp7$HE[nrow(tmp7)], 
    HI = tmp7$HI[nrow(tmp7)], 
    HR = tmp7$HR[nrow(tmp7)], 
    H_DN = tmp7$H_DN[nrow(tmp7)], 
    H_DI = tmp7$H_DI[nrow(tmp7)],
    Hog2WF_Count = tmp7$Hog2WF_Count[nrow(tmp7)]
  )
  set.seed(seed.nums[i])
  FullRooms.9to23 <- ssa(
    x0 = TransmissionModel.Int.9to23,
    a = TransmissionModel_NoMDA.ODE,
    nu = TransmissionModel.Matrix,
    parms = TransmissionModel.par,
    tf = FourteenWeek.tf,
    method = ssa.d(),
    verbose = FALSE,
    censusInterval = 0.1,
    consoleInterval = 1
  ) 
  tmp8 <- as.data.frame(FullRooms.9to23$data)
  tmp8$t <- tmp8$t +63
  tmp8 <- tmp8[tmp8$t < 161,] 
  
  ### Room 1 Exit | No MDAs 
  TransmissionModel.Int.23to25 <- c(
    #Room1
    S1 = 0,
    E1 = 0,
    I1 = 0,
    Q1 = 0,
    R1 = 0,
    Total_inf1 = tmp8$Total_inf1[nrow(tmp8)],
    True_Rec1 = tmp8$True_Rec1[nrow(tmp8)],
    Total_Q1 = tmp8$Total_Q1[nrow(tmp8)],
    D1_Nat = tmp8$D1_Nat[nrow(tmp8)],
    D1_Inf = tmp8$D1_Inf[nrow(tmp8)],
    WF2Hog_Count1 = tmp8$WF2Hog_Count1[nrow(tmp8)],
    
    #Room2
    S2 = tmp8$S2[nrow(tmp8)],
    E2 = tmp8$E2[nrow(tmp8)],
    I2 = tmp8$I2[nrow(tmp8)],
    Q2 = tmp8$Q2[nrow(tmp8)],
    R2 = tmp8$R2[nrow(tmp8)],
    Total_inf2 = tmp8$Total_inf2[nrow(tmp8)],
    True_Rec2 = tmp8$True_Rec2[nrow(tmp8)],
    Total_Q2 = tmp8$Total_Q2[nrow(tmp8)],
    D2_Nat = tmp8$D2_Nat[nrow(tmp8)],
    D2_Inf = tmp8$D2_Inf[nrow(tmp8)],
    WF2Hog_Count2 = tmp8$WF2Hog_Count2[nrow(tmp8)],
    
    #Room3
    S3 = tmp8$S3[nrow(tmp8)],
    E3 = tmp8$E3[nrow(tmp8)],
    I3 = tmp8$I3[nrow(tmp8)],
    Q3 = tmp8$Q3[nrow(tmp8)],
    R3 = tmp8$R3[nrow(tmp8)],
    Total_inf3 = tmp8$Total_inf3[nrow(tmp8)],
    True_Rec3 = tmp8$True_Rec3[nrow(tmp8)],
    Total_Q3 = tmp8$Total_Q3[nrow(tmp8)],
    D3_Nat = tmp8$D3_Nat[nrow(tmp8)],
    D3_Inf = tmp8$D3_Inf[nrow(tmp8)],
    WF2Hog_Count3 = tmp8$WF2Hog_Count3[nrow(tmp8)],
    
    #Room4
    S4 = tmp8$S4[nrow(tmp8)],
    E4 = tmp8$E4[nrow(tmp8)],
    I4 = tmp8$I4[nrow(tmp8)],
    Q4 = tmp8$Q4[nrow(tmp8)],
    R4 = tmp8$R4[nrow(tmp8)],
    Total_inf4 = tmp8$Total_inf4[nrow(tmp8)],
    True_Rec4 = tmp8$True_Rec4[nrow(tmp8)],
    Total_Q4 = tmp8$Total_Q4[nrow(tmp8)],
    D4_Nat = tmp8$D4_Nat[nrow(tmp8)],
    D4_Inf = tmp8$D4_Inf[nrow(tmp8)],
    WF2Hog_Count4 = tmp8$WF2Hog_Count4[nrow(tmp8)],
    
    #Workforce
    HS = tmp8$HS[nrow(tmp8)], 
    HE = tmp8$HE[nrow(tmp8)], 
    HI = tmp8$HI[nrow(tmp8)], 
    HR = tmp8$HR[nrow(tmp8)], 
    H_DN = tmp8$H_DN[nrow(tmp8)], 
    H_DI = tmp8$H_DI[nrow(tmp8)],
    Hog2WF_Count = tmp8$Hog2WF_Count[nrow(tmp8)]
  )
  set.seed(seed.nums[i])
  FullRooms.23to25 <- ssa(
    x0 = TransmissionModel.Int.23to25,
    a = TransmissionModel_NoMDA.ODE,
    nu = TransmissionModel.Matrix,
    parms = TransmissionModel.par,
    tf = TwoWeek.tf,
    method = ssa.d(),
    verbose = FALSE,
    censusInterval = 0.1,
    consoleInterval = 1
  ) 
  tmp9 <- as.data.frame(FullRooms.23to25$data)
  tmp9$t <- tmp9$t + 161
  tmp9 <- tmp9[tmp9$t < 175,]
  
  ### Room 2 Exit | No MDAs
  TransmissionModel.Int.25to27 <- c(
    #Room1
    S1 = 0,
    E1 = 0,
    I1 = 0,
    Q1 = 0,
    R1 = 0,
    Total_inf1 = tmp9$Total_inf1[nrow(tmp9)],
    True_Rec1 = tmp9$True_Rec1[nrow(tmp9)],
    Total_Q1 = tmp9$Total_Q1[nrow(tmp9)],
    D1_Nat = tmp9$D1_Nat[nrow(tmp9)],
    D1_Inf = tmp9$D1_Inf[nrow(tmp9)],
    WF2Hog_Count1 = tmp9$WF2Hog_Count1[nrow(tmp9)],
    
    #Room2
    S2 = 0,
    E2 = 0,
    I2 = 0,
    Q2 = 0,
    R2 = 0,
    Total_inf2 = tmp9$Total_inf2[nrow(tmp9)],
    True_Rec2 = tmp9$True_Rec2[nrow(tmp9)],
    Total_Q2 = tmp9$Total_Q2[nrow(tmp9)],
    D2_Nat = tmp9$D2_Nat[nrow(tmp9)],
    D2_Inf = tmp9$D2_Inf[nrow(tmp9)],
    WF2Hog_Count2 = tmp9$WF2Hog_Count2[nrow(tmp9)],
    
    #Room3
    S3 = tmp9$S3[nrow(tmp9)],
    E3 = tmp9$E3[nrow(tmp9)],
    I3 = tmp9$I3[nrow(tmp9)],
    Q3 = tmp9$Q3[nrow(tmp9)],
    R3 = tmp9$R3[nrow(tmp9)],
    Total_inf3 = tmp9$Total_inf3[nrow(tmp9)],
    True_Rec3 = tmp9$True_Rec3[nrow(tmp9)],
    Total_Q3 = tmp9$Total_Q3[nrow(tmp9)],
    D3_Nat = tmp9$D3_Nat[nrow(tmp9)],
    D3_Inf = tmp9$D3_Inf[nrow(tmp9)],
    WF2Hog_Count3 = tmp9$WF2Hog_Count3[nrow(tmp9)],
    
    #Room4
    S4 = tmp9$S4[nrow(tmp9)],
    E4 = tmp9$E4[nrow(tmp9)],
    I4 = tmp9$I4[nrow(tmp9)],
    Q4 = tmp9$Q4[nrow(tmp9)],
    R4 = tmp9$R4[nrow(tmp9)],
    Total_inf4 = tmp9$Total_inf4[nrow(tmp9)],
    True_Rec4 = tmp9$True_Rec4[nrow(tmp9)],
    Total_Q4 = tmp9$Total_Q4[nrow(tmp9)],
    D4_Nat = tmp9$D4_Nat[nrow(tmp9)],
    D4_Inf = tmp9$D4_Inf[nrow(tmp9)],
    WF2Hog_Count4 = tmp9$WF2Hog_Count4[nrow(tmp9)],
    
    #Workforce
    HS = tmp9$HS[nrow(tmp9)], 
    HE = tmp9$HE[nrow(tmp9)], 
    HI = tmp9$HI[nrow(tmp9)], 
    HR = tmp9$HR[nrow(tmp9)], 
    H_DN = tmp9$H_DN[nrow(tmp9)], 
    H_DI = tmp9$H_DI[nrow(tmp9)],
    Hog2WF_Count = tmp9$Hog2WF_Count[nrow(tmp9)]
  )
  set.seed(seed.nums[i])
  FullRooms.25to27 <- ssa(
    x0 = TransmissionModel.Int.25to27,
    a = TransmissionModel_NoMDA.ODE,
    nu = TransmissionModel.Matrix,
    parms = TransmissionModel.par,
    tf = TwoWeek.tf,
    method = ssa.d(),
    verbose = FALSE,
    censusInterval = 0.1,
    consoleInterval = 1
  ) 
  tmp10 <- as.data.frame(FullRooms.25to27 $data)
  tmp10$t <- tmp10$t + 175
  tmp10 <- tmp10[tmp10$t < 189,]
  
  ### Room 3 Exit | No MDAs
  TransmissionModel.Int.27to29 <- c(
    S1 = 0,
    E1 = 0,
    I1 = 0,
    Q1 = 0,
    R1 = 0,
    Total_inf1 = tmp10$Total_inf1[nrow(tmp10)],
    True_Rec1 = tmp10$True_Rec1[nrow(tmp10)],
    Total_Q1 = tmp10$Total_Q1[nrow(tmp10)],
    D1_Nat = tmp10$D1_Nat[nrow(tmp10)],
    D1_Inf = tmp10$D1_Inf[nrow(tmp10)],
    WF2Hog_Count1 = tmp10$WF2Hog_Count1[nrow(tmp10)],
    
    #Room2
    S2 = 0,
    E2 = 0,
    I2 = 0,
    Q2 = 0,
    R2 = 0,
    Total_inf2 = tmp10$Total_inf2[nrow(tmp10)],
    True_Rec2 = tmp10$True_Rec2[nrow(tmp10)],
    Total_Q2 = tmp10$Total_Q2[nrow(tmp10)],
    D2_Nat = tmp10$D2_Nat[nrow(tmp10)],
    D2_Inf = tmp10$D2_Inf[nrow(tmp10)],
    WF2Hog_Count2 = tmp10$WF2Hog_Count2[nrow(tmp10)],
    
    #Room3
    S3 = 0,
    E3 = 0,
    I3 = 0,
    Q3 = 0,
    R3 = 0,
    Total_inf3 = tmp10$Total_inf3[nrow(tmp10)],
    True_Rec3 = tmp10$True_Rec3[nrow(tmp10)],
    Total_Q3 = tmp10$Total_Q3[nrow(tmp10)],
    D3_Nat = tmp10$D3_Nat[nrow(tmp10)],
    D3_Inf = tmp10$D3_Inf[nrow(tmp10)],
    WF2Hog_Count3 = tmp10$WF2Hog_Count3[nrow(tmp10)],
    
    #Room4
    S4 = tmp10$S4[nrow(tmp10)],
    E4 = tmp10$E4[nrow(tmp10)],
    I4 = tmp10$I4[nrow(tmp10)],
    Q4 = tmp10$Q4[nrow(tmp10)],
    R4 = tmp10$R4[nrow(tmp10)],
    Total_inf4 = tmp10$Total_inf4[nrow(tmp10)],
    True_Rec4 = tmp10$True_Rec4[nrow(tmp10)],
    Total_Q4 = tmp10$Total_Q4[nrow(tmp10)],
    D4_Nat = tmp10$D4_Nat[nrow(tmp10)],
    D4_Inf = tmp10$D4_Inf[nrow(tmp10)],
    WF2Hog_Count4 = tmp10$WF2Hog_Count4[nrow(tmp10)],
    
    #Workforce
    HS = tmp10$HS[nrow(tmp10)], 
    HE = tmp10$HE[nrow(tmp10)], 
    HI = tmp10$HI[nrow(tmp10)], 
    HR = tmp10$HR[nrow(tmp10)], 
    H_DN = tmp10$H_DN[nrow(tmp10)], 
    H_DI = tmp10$H_DI[nrow(tmp10)],
    Hog2WF_Count = tmp10$Hog2WF_Count[nrow(tmp10)]
  )
  set.seed(seed.nums[i])
  FullRooms.27to29 <- ssa(
    x0 = TransmissionModel.Int.27to29,
    a = TransmissionModel_NoMDA.ODE,
    nu = TransmissionModel.Matrix,
    parms = TransmissionModel.par,
    tf = TwoWeek.tf,
    method = ssa.d(),
    verbose = FALSE,
    censusInterval = 0.1,
    consoleInterval = 1
  ) 
  tmp11 <- as.data.frame(FullRooms.27to29 $data)
  tmp11$t <- tmp11$t + 189
  tmp11 <- tmp11[tmp11$t < 203,]
  
  ### Room 4 Exit | No MDAs
  TransmissionModel.Int.29to31 <- c(
    #Room1
    S1 = 0,
    E1 = 0,
    I1 = 0,
    Q1 = 0,
    R1 = 0,
    Total_inf1 = tmp11$Total_inf1[nrow(tmp11)],
    True_Rec1 = tmp11$True_Rec1[nrow(tmp11)],
    Total_Q1 = tmp11$Total_Q1[nrow(tmp11)],
    D1_Nat = tmp11$D1_Nat[nrow(tmp11)],
    D1_Inf = tmp11$D1_Inf[nrow(tmp11)],
    WF2Hog_Count1 = tmp11$WF2Hog_Count1[nrow(tmp11)],
    
    #Room2
    S2 = 0,
    E2 = 0,
    I2 = 0,
    Q2 = 0,
    R2 = 0,
    Total_inf2 = tmp11$Total_inf2[nrow(tmp11)],
    True_Rec2 = tmp11$True_Rec2[nrow(tmp11)],
    Total_Q2 = tmp11$Total_Q2[nrow(tmp11)],
    D2_Nat = tmp11$D2_Nat[nrow(tmp11)],
    D2_Inf = tmp11$D2_Inf[nrow(tmp11)],
    WF2Hog_Count2 = tmp11$WF2Hog_Count2[nrow(tmp11)],
    
    #Room3
    S3 = 0,
    E3 = 0,
    I3 = 0,
    Q3 = 0,
    R3 = 0,
    Total_inf3 = tmp11$Total_inf3[nrow(tmp11)],
    True_Rec3 = tmp11$True_Rec3[nrow(tmp11)],
    Total_Q3 = tmp11$Total_Q3[nrow(tmp11)],
    D3_Nat = tmp11$D3_Nat[nrow(tmp11)],
    D3_Inf = tmp11$D3_Inf[nrow(tmp11)],
    WF2Hog_Count3 = tmp11$WF2Hog_Count3[nrow(tmp11)],
    
    #Room4
    S4 = 0,
    E4 = 0,
    I4 = 0,
    Q4 = 0,
    R4 = 0,
    Total_inf4 = tmp11$Total_inf4[nrow(tmp11)],
    True_Rec4 = tmp11$True_Rec4[nrow(tmp11)],
    Total_Q4 = tmp11$Total_Q4[nrow(tmp11)],
    D4_Nat = tmp11$D4_Nat[nrow(tmp11)],
    D4_Inf = tmp11$D4_Inf[nrow(tmp11)],
    WF2Hog_Count4 = tmp11$WF2Hog_Count4[nrow(tmp11)],
    
    #Workforce
    HS = tmp11$HS[nrow(tmp11)], 
    HE = tmp11$HE[nrow(tmp11)], 
    HI = tmp11$HI[nrow(tmp11)], 
    HR = tmp11$HR[nrow(tmp11)], 
    H_DN = tmp11$H_DN[nrow(tmp11)], 
    H_DI = tmp11$H_DI[nrow(tmp11)],
    Hog2WF_Count = tmp11$Hog2WF_Count[nrow(tmp11)]
  )
  set.seed(seed.nums[i])
  FullRooms.29to31 <- ssa(
    x0 = TransmissionModel.Int.29to31,
    a = TransmissionModel_NoMDA.ODE,
    nu = TransmissionModel.Matrix,
    parms = TransmissionModel.par,
    tf = TwoWeek.tf,
    method = ssa.d(),
    verbose = FALSE,
    censusInterval = 0.1,
    consoleInterval = 1
  ) 
  tmp12 <- as.data.frame(FullRooms.29to31 $data)
  tmp12$t <- tmp12$t + 203
  
  tmp.final <- rbind(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12)
  tmp.final$iter <- as.factor(i)
  results <- rbind(results, tmp.final)
}
saveRDS(results, file = "Current code/Reviewer Update Code and Models/AllCMs/Indirect by 178/R10_Room4_.33Qrate_60VaccEff.rds")

### Reporting Values ####
results <- readRDS()
## % iterations without WF infection
PerIterInf <- results %>%
  group_by(iter) %>%
  summarise(total_I = sum(HI)) %>%
  filter(total_I != 0) %>%
  select(iter)
## % witout infection calculation
(1 - (nrow(PerIterInf)/n.itr))*100

## Time to First WF infection
results %>%
  # HT_0PerVacc %>%
  arrange(iter, t) %>%
  filter(HI == 1) %>%
  group_by(iter) %>%
  filter(row_number() == 1) %>%
  mutate(TimeSineInfIntro = t - 42) %>%
  ungroup() %>%
  summarise(#mean_time = mean(t),
    #mean_realtime = mean(TimeSineInfIntro),
    #med_time = median(t),
    med_realtime = median(TimeSineInfIntro),
    LB_realtime = quantile(TimeSineInfIntro, probs = .05),
    UB_realtime = quantile(TimeSineInfIntro, probs = .95))

## Time to Second Room Infection
results %>% 
  group_by(iter) %>% 
  summarise(t_cross = t[which(E1 >= 1 | 
                                E2 >= 1 |
                                E3 >= 1)]) %>% 
  summarise(min_t_cross = min(t_cross)) %>% 
  ungroup() %>% 
  summarise(med_t_cross = median(min_t_cross) - 42,
            LB_t_cross = quantile(min_t_cross, probs = 0.05) - 42,
            UB_t_cross = quantile(min_t_cross, probs = 0.95) - 42)
#Total infected and Recovered
results %>% 
  #select(t, iter, Total_inf1, Total_inf2, Total_inf3, Total_inf4) %>% 
  mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
         Tot_rec = True_Rec1 + True_Rec2 + True_Rec3 + True_Rec4) %>% 
  group_by(iter) %>% 
  summarise(Tot_inf2 = max(Tot_inf),
            Tot_rec2 = max(Tot_rec)) %>% 
  ungroup() %>% 
  summarise(#med_inf = median(Tot_inf2),
    # LB_med_inf = quantile(Tot_inf2, probs = 0.05),
    # UB_med_inf = quantile(Tot_inf2, probs = 0.95),
    med_rec = median(Tot_rec2),
    LB_med_rec = quantile(Tot_rec2, probs = 0.05),
    UB_med_rec = quantile(Tot_rec2, probs = 0.95))

## Minimum Time to Peak Infection
results %>% 
  mutate(S_tot = S1+S2+S3+S4,                             
         E_tot = E1+E2+E3+E4,
         I_tot = I1+I2+I3+I4,
         R_tot = R1+R2+R3+R4) %>%    
  group_by(iter) %>%                                      
  summarise(max_I = max(I_tot),
            t_peak = t[which(I_tot == max(I_tot))]) %>% 
  # filter(max_I != 0) %>%  #had to add this to remove when no hogs are infected. occurs when WF is the only source of initial infection
  summarise(min_t_peak = min(t_peak)) %>% # doing this because there are some iterations that have a couple time points that have max infectd hogs
  ungroup() %>% 
  summarise(med_t_peak = median(min_t_peak) - 42,
            LB_t_peak = quantile(min_t_peak, probs = 0.05) - 42,
            UB_t_peak = quantile(min_t_peak, probs = 0.95) - 42)

#### % With WF Infection Bar Graph ####
#Want to make columns for DT type (high or low), Vaccine Eff, and Indirect transmission (Baseline vs. Improved)
HT_0PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/HT_0VaccEff.rds") %>% 
  select(iter, t, HI)
HT_40PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/HT_40VaccEff.rds")%>% 
  select(iter, t, HI)
HT_80PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/HT_80VaccEff.rds")%>% 
  select(iter, t, HI)
#High Direct Transmission, 0 Vacc Eff, Baseline Indirect
HT_0PerVacc_IDby178$Trnsm<- as.factor("High")
HT_0PerVacc_IDby178$Vacc_Eff <- as.factor(0)
HT_0PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")
#High Direct Transmission, 40 Vacc Eff, Baseline Indirect
HT_40PerVacc_IDby178$Trnsm<- as.factor("High")
HT_40PerVacc_IDby178$Vacc_Eff <- as.factor(40)
HT_40PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")
#High Direct Transmission, 80 Vacc Eff, Baseline Indirect
HT_80PerVacc_IDby178$Trnsm<- as.factor("High")
HT_80PerVacc_IDby178$Vacc_Eff <- as.factor(80)
HT_80PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")

HT_Baseline <- rbind(HT_0PerVacc_IDby178, HT_40PerVacc_IDby178, HT_80PerVacc_IDby178)
rm(HT_0PerVacc_IDby178, HT_40PerVacc_IDby178, HT_80PerVacc_IDby178)

LT_0PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/LT_0VaccEff.rds")%>% 
  select(iter, t, HI)
LT_40PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/LT_40VaccEff.rds")%>% 
  select(iter, t, HI)
LT_80PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/LT_80VaccEff.rds")%>% 
  select(iter, t, HI)
#Low Direct Transmission, 0 Vacc Eff, Baseline Indirect
LT_0PerVacc_IDby178$Trnsm<- as.factor("Low")
LT_0PerVacc_IDby178$Vacc_Eff <- as.factor(0)
LT_0PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")
#Low Direct Transmission, 40 Vacc Eff, Baseline Indirect
LT_40PerVacc_IDby178$Trnsm<- as.factor("Low")
LT_40PerVacc_IDby178$Vacc_Eff <- as.factor(40)
LT_40PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")
#Low Direct Transmission, 80 Vacc Eff, Baseline Indirect
LT_80PerVacc_IDby178$Trnsm<- as.factor("Low")
LT_80PerVacc_IDby178$Vacc_Eff <- as.factor(80)
LT_80PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")

LT_Baseline <- rbind(LT_0PerVacc_IDby178, LT_40PerVacc_IDby178, LT_80PerVacc_IDby178)
rm(LT_0PerVacc_IDby178, LT_40PerVacc_IDby178, LT_80PerVacc_IDby178)

HT_0PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/HT_0VaccEff.rds")%>% 
  select(iter, t, HI)
HT_40PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/HT_40VaccEff.rds")%>% 
  select(iter, t, HI)
HT_80PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/HT_80VaccEff.rds")%>% 
  select(iter, t, HI)
#High Direct Transmission, 0 Vacc Eff, Improved Indirect
HT_0PerVacc_IDby500$Trnsm<- as.factor("High")
HT_0PerVacc_IDby500$Vacc_Eff <- as.factor(0)
HT_0PerVacc_IDby500$BtwnRoom <- as.factor("Improved")
#High Direct Transmission, 40 Vacc Eff, Improved Indirect
HT_40PerVacc_IDby500$Trnsm<- as.factor("High")
HT_40PerVacc_IDby500$Vacc_Eff <- as.factor(40)
HT_40PerVacc_IDby500$BtwnRoom <- as.factor("Improved")
#High Direct Transmission, 80 Vacc Eff, Improved Indirect
HT_80PerVacc_IDby500$Trnsm<- as.factor("High")
HT_80PerVacc_IDby500$Vacc_Eff <- as.factor(80)
HT_80PerVacc_IDby500$BtwnRoom <- as.factor("Improved")

HT_improved <- rbind(HT_0PerVacc_IDby500, HT_40PerVacc_IDby500, HT_80PerVacc_IDby500)
rm(HT_0PerVacc_IDby500, HT_40PerVacc_IDby500, HT_80PerVacc_IDby500)

LT_0PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/LT_0VaccEff.rds")%>% 
  select(iter, t, HI)
LT_40PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/LT_40VaccEff.rds")%>% 
  select(iter, t, HI)
LT_80PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/LT_80VaccEff.rds")%>% 
  select(iter, t, HI)
#Low Direct Transmission, 0 Vacc Eff, Improved Indirect
LT_0PerVacc_IDby500$Trnsm<- as.factor("Low")
LT_0PerVacc_IDby500$Vacc_Eff <- as.factor(0)
LT_0PerVacc_IDby500$BtwnRoom <- as.factor("Improved")
#Low Direct Transmission, 40 Vacc Eff, Improved Indirect
LT_40PerVacc_IDby500$Trnsm<- as.factor("Low")
LT_40PerVacc_IDby500$Vacc_Eff <- as.factor(40)
LT_40PerVacc_IDby500$BtwnRoom <- as.factor("Improved")
#Low Direct Transmission, 80 Vacc Eff, Improved Indirect
LT_80PerVacc_IDby500$Trnsm<- as.factor("Low")
LT_80PerVacc_IDby500$Vacc_Eff <- as.factor(80)
LT_80PerVacc_IDby500$BtwnRoom <- as.factor("Improved")

LT_improved <- rbind(LT_0PerVacc_IDby500, LT_40PerVacc_IDby500, LT_80PerVacc_IDby500)
rm(LT_0PerVacc_IDby500, LT_40PerVacc_IDby500, LT_80PerVacc_IDby500)

All <- rbind(HT_Baseline, LT_Baseline, HT_improved, LT_improved)

All %>% 
  group_by(Trnsm, BtwnRoom, Vacc_Eff, iter) %>% 
  summarise(total_I = sum(HI)) %>%
  filter(total_I != 0) %>% 
  ungroup() %>% 
  group_by(Trnsm, BtwnRoom, Vacc_Eff) %>% 
  summarise(countwithWFInf = n(),
            PerwithWFInf = (n()/5000)*100,
            PerWOWFinf = (1-(n()/5000))*100) %>% 
  mutate(Vacc_Eff = factor(Vacc_Eff, levels = c("0","40", "80"), labels = c("0%", "40%", "80%")),
         BtwnRoom = factor(BtwnRoom, levels = c("Baseline", "Improved")),
         Trnsm = factor(Trnsm, levels = c("High", "Low"), labels = c("High Transmission\nStrain", "Low Transmission\nStrain"))) %>%
  ungroup() %>% 
  ggplot(aes(x = Vacc_Eff, y = PerwithWFInf, fill = BtwnRoom))+
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  scale_fill_manual(name = "Between Room\nTransmission",
                    values = c("#7D0900", "#FFFFBF"),
                    labels = c("Baseline", "Improved"))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0" = "0%", "25" = "25%", "50" = "50%", "75" = "75%", "100" = "100%"))+
  facet_wrap(~Trnsm)+
  labs(x = "Vaccine Efficacy",
       y = "Iterations with Workforce Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        legend.title = element_text(face = "bold", size = 28),
        legend.text = element_text(face = "bold", size = 26),
        legend.title.align = 0.5,
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 24),
        axis.title = element_text(face = "bold", size = 28),
        axis.text = element_text(face = "bold", size = 26, color = "black"))
#1356 x 863

#############################
#### Time to 1st WF inf ####
#############################
HT_0PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/HT_0VaccEff.rds") %>% 
  select(iter, t, HI)
HT_40PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/HT_40VaccEff.rds")%>% 
  select(iter, t, HI)
HT_80PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/HT_80VaccEff.rds")%>% 
  select(iter, t, HI)
#High Direct Transmission, 0 Vacc Eff, Baseline Indirect
HT_0PerVacc_IDby178$Trnsm<- as.factor("High")
HT_0PerVacc_IDby178$Vacc_Eff <- as.factor(0)
HT_0PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")
#High Direct Transmission, 40 Vacc Eff, Baseline Indirect
HT_40PerVacc_IDby178$Trnsm<- as.factor("High")
HT_40PerVacc_IDby178$Vacc_Eff <- as.factor(40)
HT_40PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")
#High Direct Transmission, 80 Vacc Eff, Baseline Indirect
HT_80PerVacc_IDby178$Trnsm<- as.factor("High")
HT_80PerVacc_IDby178$Vacc_Eff <- as.factor(80)
HT_80PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")

HT_Baseline <- rbind(HT_0PerVacc_IDby178, HT_40PerVacc_IDby178, HT_80PerVacc_IDby178)
rm(HT_0PerVacc_IDby178, HT_40PerVacc_IDby178, HT_80PerVacc_IDby178)

LT_0PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/LT_0VaccEff.rds")%>% 
  select(iter, t, HI)
LT_40PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/LT_40VaccEff.rds")%>% 
  select(iter, t, HI)
LT_80PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/LT_80VaccEff.rds")%>% 
  select(iter, t, HI)
#Low Direct Transmission, 0 Vacc Eff, Baseline Indirect
LT_0PerVacc_IDby178$Trnsm<- as.factor("Low")
LT_0PerVacc_IDby178$Vacc_Eff <- as.factor(0)
LT_0PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")
#Low Direct Transmission, 40 Vacc Eff, Baseline Indirect
LT_40PerVacc_IDby178$Trnsm<- as.factor("Low")
LT_40PerVacc_IDby178$Vacc_Eff <- as.factor(40)
LT_40PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")
#Low Direct Transmission, 80 Vacc Eff, Baseline Indirect
LT_80PerVacc_IDby178$Trnsm<- as.factor("Low")
LT_80PerVacc_IDby178$Vacc_Eff <- as.factor(80)
LT_80PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")

LT_Baseline <- rbind(LT_0PerVacc_IDby178, LT_40PerVacc_IDby178, LT_80PerVacc_IDby178)
rm(LT_0PerVacc_IDby178, LT_40PerVacc_IDby178, LT_80PerVacc_IDby178)

HT_0PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/HT_0VaccEff.rds")%>% 
  select(iter, t, HI)
HT_40PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/HT_40VaccEff.rds")%>% 
  select(iter, t, HI)
HT_80PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/HT_80VaccEff.rds")%>% 
  select(iter, t, HI)
#High Direct Transmission, 0 Vacc Eff, Improved Indirect
HT_0PerVacc_IDby500$Trnsm<- as.factor("High")
HT_0PerVacc_IDby500$Vacc_Eff <- as.factor(0)
HT_0PerVacc_IDby500$BtwnRoom <- as.factor("Improved")
#High Direct Transmission, 40 Vacc Eff, Improved Indirect
HT_40PerVacc_IDby500$Trnsm<- as.factor("High")
HT_40PerVacc_IDby500$Vacc_Eff <- as.factor(40)
HT_40PerVacc_IDby500$BtwnRoom <- as.factor("Improved")
#High Direct Transmission, 80 Vacc Eff, Improved Indirect
HT_80PerVacc_IDby500$Trnsm<- as.factor("High")
HT_80PerVacc_IDby500$Vacc_Eff <- as.factor(80)
HT_80PerVacc_IDby500$BtwnRoom <- as.factor("Improved")

HT_improved <- rbind(HT_0PerVacc_IDby500, HT_40PerVacc_IDby500, HT_80PerVacc_IDby500)
rm(HT_0PerVacc_IDby500, HT_40PerVacc_IDby500, HT_80PerVacc_IDby500)

LT_0PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/LT_0VaccEff.rds")%>% 
  select(iter, t, HI)
LT_40PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/LT_40VaccEff.rds")%>% 
  select(iter, t, HI)
LT_80PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/LT_80VaccEff.rds")%>% 
  select(iter, t, HI)
#Low Direct Transmission, 0 Vacc Eff, Improved Indirect
LT_0PerVacc_IDby500$Trnsm<- as.factor("Low")
LT_0PerVacc_IDby500$Vacc_Eff <- as.factor(0)
LT_0PerVacc_IDby500$BtwnRoom <- as.factor("Improved")
#Low Direct Transmission, 40 Vacc Eff, Improved Indirect
LT_40PerVacc_IDby500$Trnsm<- as.factor("Low")
LT_40PerVacc_IDby500$Vacc_Eff <- as.factor(40)
LT_40PerVacc_IDby500$BtwnRoom <- as.factor("Improved")
#Low Direct Transmission, 80 Vacc Eff, Improved Indirect
LT_80PerVacc_IDby500$Trnsm<- as.factor("Low")
LT_80PerVacc_IDby500$Vacc_Eff <- as.factor(80)
LT_80PerVacc_IDby500$BtwnRoom <- as.factor("Improved")

LT_improved <- rbind(LT_0PerVacc_IDby500, LT_40PerVacc_IDby500, LT_80PerVacc_IDby500)
rm(LT_0PerVacc_IDby500, LT_40PerVacc_IDby500, LT_80PerVacc_IDby500)

All <- rbind(HT_Baseline, LT_Baseline, HT_improved, LT_improved)

All %>% 
  arrange(Trnsm, BtwnRoom, Vacc_Eff, iter, t) %>%
  filter(HI == 1) %>%
  group_by(Trnsm, BtwnRoom, Vacc_Eff, iter) %>%
  filter(row_number() == 1) %>%
  mutate(TimeSineInfIntro = t - 42) %>% 
  mutate(Vacc_Eff = factor(Vacc_Eff, levels = c("0","40", "80"), labels = c("0%", "40%", "80%")),
         BtwnRoom = factor(BtwnRoom, levels = c("Baseline", "Improved")),
         Trnsm = factor(Trnsm, levels = c("High", "Low"), labels = c("High Transmission\nStrain", "Low Transmission\nStrain"))) %>%
  ggplot(aes(x = Vacc_Eff, y = TimeSineInfIntro, fill = BtwnRoom))+
  geom_boxplot(outlier.shape = 21, outlier.fill =  NULL, outlier.size = 2.5, notch = F, color = "black")+
  scale_fill_manual(name = "Between Room\nTransmission",
                    values = c("#7D0900", "#FFFFBF"),
                    labels = c("Baseline", "Improved"))+ 
  facet_wrap(~Trnsm,
             scales = "free")+
  # facet_wrap(~Trnsm)+
  labs(x = "Vaccine Efficacy",
       y = "Time to 1st Workforce Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        legend.title = element_text(face = "bold", size = 28),
        legend.text = element_text(face = "bold", size = 26),
        legend.title.align = 0.5,
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 24),
        axis.title = element_text(face = "bold", size = 28),
        axis.text = element_text(face = "bold", size = 26, color = "black"))
#1356 x 863

#############################
#### Total Hogs Infected ####
#############################
HT_0PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/HT_0VaccEff.rds") %>% 
  select(iter, t, True_Rec1, True_Rec2, True_Rec3, True_Rec4)
HT_40PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/HT_40VaccEff.rds")%>% 
  select(iter, t, True_Rec1, True_Rec2, True_Rec3, True_Rec4)
HT_80PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/HT_80VaccEff.rds")%>% 
  select(iter, t, True_Rec1, True_Rec2, True_Rec3, True_Rec4)
#High Direct Transmission, 0 Vacc Eff, Baseline Indirect
HT_0PerVacc_IDby178$Trnsm<- as.factor("High")
HT_0PerVacc_IDby178$Vacc_Eff <- as.factor(0)
HT_0PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")
#High Direct Transmission, 40 Vacc Eff, Baseline Indirect
HT_40PerVacc_IDby178$Trnsm<- as.factor("High")
HT_40PerVacc_IDby178$Vacc_Eff <- as.factor(40)
HT_40PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")
#High Direct Transmission, 80 Vacc Eff, Baseline Indirect
HT_80PerVacc_IDby178$Trnsm<- as.factor("High")
HT_80PerVacc_IDby178$Vacc_Eff <- as.factor(80)
HT_80PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")

HT_Baseline <- rbind(HT_0PerVacc_IDby178, HT_40PerVacc_IDby178, HT_80PerVacc_IDby178)
rm(HT_0PerVacc_IDby178, HT_40PerVacc_IDby178, HT_80PerVacc_IDby178)

LT_0PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/LT_0VaccEff.rds")%>% 
  select(iter, t, True_Rec1, True_Rec2, True_Rec3, True_Rec4)
LT_40PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/LT_40VaccEff.rds")%>% 
  select(iter, t, True_Rec1, True_Rec2, True_Rec3, True_Rec4)
LT_80PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/LT_80VaccEff.rds")%>% 
  select(iter, t, True_Rec1, True_Rec2, True_Rec3, True_Rec4)
#Low Direct Transmission, 0 Vacc Eff, Baseline Indirect
LT_0PerVacc_IDby178$Trnsm<- as.factor("Low")
LT_0PerVacc_IDby178$Vacc_Eff <- as.factor(0)
LT_0PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")
#Low Direct Transmission, 40 Vacc Eff, Baseline Indirect
LT_40PerVacc_IDby178$Trnsm<- as.factor("Low")
LT_40PerVacc_IDby178$Vacc_Eff <- as.factor(40)
LT_40PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")
#Low Direct Transmission, 80 Vacc Eff, Baseline Indirect
LT_80PerVacc_IDby178$Trnsm<- as.factor("Low")
LT_80PerVacc_IDby178$Vacc_Eff <- as.factor(80)
LT_80PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")

LT_Baseline <- rbind(LT_0PerVacc_IDby178, LT_40PerVacc_IDby178, LT_80PerVacc_IDby178)
rm(LT_0PerVacc_IDby178, LT_40PerVacc_IDby178, LT_80PerVacc_IDby178)

HT_0PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/HT_0VaccEff.rds")%>% 
  select(iter, t, True_Rec1, True_Rec2, True_Rec3, True_Rec4)
HT_40PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/HT_40VaccEff.rds")%>% 
  select(iter, t, True_Rec1, True_Rec2, True_Rec3, True_Rec4)
HT_80PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/HT_80VaccEff.rds")%>% 
  select(iter, t, True_Rec1, True_Rec2, True_Rec3, True_Rec4)
#High Direct Transmission, 0 Vacc Eff, Improved Indirect
HT_0PerVacc_IDby500$Trnsm<- as.factor("High")
HT_0PerVacc_IDby500$Vacc_Eff <- as.factor(0)
HT_0PerVacc_IDby500$BtwnRoom <- as.factor("Improved")
#High Direct Transmission, 40 Vacc Eff, Improved Indirect
HT_40PerVacc_IDby500$Trnsm<- as.factor("High")
HT_40PerVacc_IDby500$Vacc_Eff <- as.factor(40)
HT_40PerVacc_IDby500$BtwnRoom <- as.factor("Improved")
#High Direct Transmission, 80 Vacc Eff, Improved Indirect
HT_80PerVacc_IDby500$Trnsm<- as.factor("High")
HT_80PerVacc_IDby500$Vacc_Eff <- as.factor(80)
HT_80PerVacc_IDby500$BtwnRoom <- as.factor("Improved")

HT_improved <- rbind(HT_0PerVacc_IDby500, HT_40PerVacc_IDby500, HT_80PerVacc_IDby500)
rm(HT_0PerVacc_IDby500, HT_40PerVacc_IDby500, HT_80PerVacc_IDby500)

LT_0PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/LT_0VaccEff.rds")%>% 
  select(iter, t, True_Rec1, True_Rec2, True_Rec3, True_Rec4)
LT_40PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/LT_40VaccEff.rds")%>% 
  select(iter, t, True_Rec1, True_Rec2, True_Rec3, True_Rec4)
LT_80PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/LT_80VaccEff.rds")%>% 
  select(iter, t, True_Rec1, True_Rec2, True_Rec3, True_Rec4)
#Low Direct Transmission, 0 Vacc Eff, Improved Indirect
LT_0PerVacc_IDby500$Trnsm<- as.factor("Low")
LT_0PerVacc_IDby500$Vacc_Eff <- as.factor(0)
LT_0PerVacc_IDby500$BtwnRoom <- as.factor("Improved")
#Low Direct Transmission, 40 Vacc Eff, Improved Indirect
LT_40PerVacc_IDby500$Trnsm<- as.factor("Low")
LT_40PerVacc_IDby500$Vacc_Eff <- as.factor(40)
LT_40PerVacc_IDby500$BtwnRoom <- as.factor("Improved")
#Low Direct Transmission, 80 Vacc Eff, Improved Indirect
LT_80PerVacc_IDby500$Trnsm<- as.factor("Low")
LT_80PerVacc_IDby500$Vacc_Eff <- as.factor(80)
LT_80PerVacc_IDby500$BtwnRoom <- as.factor("Improved")

LT_improved <- rbind(LT_0PerVacc_IDby500, LT_40PerVacc_IDby500, LT_80PerVacc_IDby500)
rm(LT_0PerVacc_IDby500, LT_40PerVacc_IDby500, LT_80PerVacc_IDby500)

All <- rbind(HT_Baseline, LT_Baseline, HT_improved, LT_improved)

All %>% 
  mutate(Tot_rec = True_Rec1 + True_Rec2 + True_Rec3 + True_Rec4) %>% 
  group_by(Trnsm, BtwnRoom, Vacc_Eff, iter) %>% 
  summarise(Tot_rec2 = max(Tot_rec)) %>% 
  mutate(Vacc_Eff = factor(Vacc_Eff, levels = c("0","40", "80"), labels = c("0%", "40%", "80%")),
         BtwnRoom = factor(BtwnRoom, levels = c("Baseline", "Improved")),
         Trnsm = factor(Trnsm, levels = c("High", "Low"), labels = c("High Transmission\nStrain", "Low Transmission\nStrain"))) %>%
  ungroup() %>% 
  ggplot(aes(x = Vacc_Eff, y = Tot_rec2, fill = BtwnRoom))+
  geom_boxplot(outlier.shape = 21, outlier.fill =  NULL, outlier.size = 2.5, notch = F, color = "black")+
  scale_fill_manual(name = "Between Room\nTransmission",
                    values = c("#7D0900", "#FFFFBF"),
                    labels = c("Baseline", "Improved"))+ 
  scale_y_continuous(expand = c(0.02, 0), limits = c(0, 4000))+
  facet_wrap(~Trnsm)+
  labs(x = "Vaccine Efficacy",
       y = "Total Infected Hogs")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        legend.title = element_text(face = "bold", size = 28),
        legend.text = element_text(face = "bold", size = 26),
        legend.title.align = 0.5,
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 24),
        axis.title = element_text(face = "bold", size = 28),
        axis.text = element_text(face = "bold", size = 26, color = "black"))
#1356 x 863

#############################
#### Time to 2nd room Inf ####
#############################
HT_0PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/HT_0VaccEff.rds") %>% 
  select(iter, t, E1, E2, E3)
HT_40PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/HT_40VaccEff.rds")%>% 
  select(iter, t, E1, E2, E3)
HT_80PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/HT_80VaccEff.rds")%>% 
  select(iter, t, E1, E2, E3)
#High Direct Transmission, 0 Vacc Eff, Baseline Indirect
HT_0PerVacc_IDby178$Trnsm<- as.factor("High")
HT_0PerVacc_IDby178$Vacc_Eff <- as.factor(0)
HT_0PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")
#High Direct Transmission, 40 Vacc Eff, Baseline Indirect
HT_40PerVacc_IDby178$Trnsm<- as.factor("High")
HT_40PerVacc_IDby178$Vacc_Eff <- as.factor(40)
HT_40PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")
#High Direct Transmission, 80 Vacc Eff, Baseline Indirect
HT_80PerVacc_IDby178$Trnsm<- as.factor("High")
HT_80PerVacc_IDby178$Vacc_Eff <- as.factor(80)
HT_80PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")

HT_Baseline <- rbind(HT_0PerVacc_IDby178, HT_40PerVacc_IDby178, HT_80PerVacc_IDby178)
rm(HT_0PerVacc_IDby178, HT_40PerVacc_IDby178, HT_80PerVacc_IDby178)

LT_0PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/LT_0VaccEff.rds")%>% 
  select(iter, t, E1, E2, E3)
LT_40PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/LT_40VaccEff.rds")%>% 
  select(iter, t, E1, E2, E3)
LT_80PerVacc_IDby178 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 178/LT_80VaccEff.rds")%>% 
  select(iter, t, E1, E2, E3)
#Low Direct Transmission, 0 Vacc Eff, Baseline Indirect
LT_0PerVacc_IDby178$Trnsm<- as.factor("Low")
LT_0PerVacc_IDby178$Vacc_Eff <- as.factor(0)
LT_0PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")
#Low Direct Transmission, 40 Vacc Eff, Baseline Indirect
LT_40PerVacc_IDby178$Trnsm<- as.factor("Low")
LT_40PerVacc_IDby178$Vacc_Eff <- as.factor(40)
LT_40PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")
#Low Direct Transmission, 80 Vacc Eff, Baseline Indirect
LT_80PerVacc_IDby178$Trnsm<- as.factor("Low")
LT_80PerVacc_IDby178$Vacc_Eff <- as.factor(80)
LT_80PerVacc_IDby178$BtwnRoom <- as.factor("Baseline")

LT_Baseline <- rbind(LT_0PerVacc_IDby178, LT_40PerVacc_IDby178, LT_80PerVacc_IDby178)
rm(LT_0PerVacc_IDby178, LT_40PerVacc_IDby178, LT_80PerVacc_IDby178)

HT_0PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/HT_0VaccEff.rds")%>% 
  select(iter, t, E1, E2, E3)
HT_40PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/HT_40VaccEff.rds")%>% 
  select(iter, t, E1, E2, E3)
HT_80PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/HT_80VaccEff.rds")%>% 
  select(iter, t, E1, E2, E3)
#High Direct Transmission, 0 Vacc Eff, Improved Indirect
HT_0PerVacc_IDby500$Trnsm<- as.factor("High")
HT_0PerVacc_IDby500$Vacc_Eff <- as.factor(0)
HT_0PerVacc_IDby500$BtwnRoom <- as.factor("Improved")
#High Direct Transmission, 40 Vacc Eff, Improved Indirect
HT_40PerVacc_IDby500$Trnsm<- as.factor("High")
HT_40PerVacc_IDby500$Vacc_Eff <- as.factor(40)
HT_40PerVacc_IDby500$BtwnRoom <- as.factor("Improved")
#High Direct Transmission, 80 Vacc Eff, Improved Indirect
HT_80PerVacc_IDby500$Trnsm<- as.factor("High")
HT_80PerVacc_IDby500$Vacc_Eff <- as.factor(80)
HT_80PerVacc_IDby500$BtwnRoom <- as.factor("Improved")

HT_improved <- rbind(HT_0PerVacc_IDby500, HT_40PerVacc_IDby500, HT_80PerVacc_IDby500)
rm(HT_0PerVacc_IDby500, HT_40PerVacc_IDby500, HT_80PerVacc_IDby500)

LT_0PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/LT_0VaccEff.rds")%>% 
  select(iter, t, E1, E2, E3)
LT_40PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/LT_40VaccEff.rds")%>% 
  select(iter, t, E1, E2, E3)
LT_80PerVacc_IDby500 <- readRDS(file = "Current code/Control Meausres/Vaccination/Indirect is DT by 500/LT_80VaccEff.rds")%>% 
  select(iter, t, E1, E2, E3)
#Low Direct Transmission, 0 Vacc Eff, Improved Indirect
LT_0PerVacc_IDby500$Trnsm<- as.factor("Low")
LT_0PerVacc_IDby500$Vacc_Eff <- as.factor(0)
LT_0PerVacc_IDby500$BtwnRoom <- as.factor("Improved")
#Low Direct Transmission, 40 Vacc Eff, Improved Indirect
LT_40PerVacc_IDby500$Trnsm<- as.factor("Low")
LT_40PerVacc_IDby500$Vacc_Eff <- as.factor(40)
LT_40PerVacc_IDby500$BtwnRoom <- as.factor("Improved")
#Low Direct Transmission, 80 Vacc Eff, Improved Indirect
LT_80PerVacc_IDby500$Trnsm<- as.factor("Low")
LT_80PerVacc_IDby500$Vacc_Eff <- as.factor(80)
LT_80PerVacc_IDby500$BtwnRoom <- as.factor("Improved")

LT_improved <- rbind(LT_0PerVacc_IDby500, LT_40PerVacc_IDby500, LT_80PerVacc_IDby500)
rm(LT_0PerVacc_IDby500, LT_40PerVacc_IDby500, LT_80PerVacc_IDby500)

All <- rbind(HT_Baseline, LT_Baseline, HT_improved, LT_improved)

All %>% 
  group_by(Trnsm, BtwnRoom, Vacc_Eff, iter) %>% 
  summarise(t_cross = t[which(E1 >= 1 | 
                                E2 >= 1 |
                                E3 >= 1)]) %>% 
  summarise(min_t_cross = min(t_cross)) %>%
  ungroup() %>% 
  mutate(min_t_cross2 = min_t_cross - 42) %>% 
  mutate(Vacc_Eff = factor(Vacc_Eff, levels = c("0","40", "80"), labels = c("0%", "40%", "80%")),
         BtwnRoom = factor(BtwnRoom, levels = c("Baseline", "Improved")),
         Trnsm = factor(Trnsm, levels = c("High", "Low"), labels = c("High Transmission\nStrain", "Low Transmission\nStrain"))) %>%
  ggplot(aes(x = Vacc_Eff, y = min_t_cross2, fill = BtwnRoom))+
  geom_boxplot(outlier.shape = 21, outlier.fill =  NULL, outlier.size = 2.5, notch = F, color = "black")+
  scale_fill_manual(name = "Between Room\nTransmission",
                    values = c("#7D0900", "#FFFFBF"),
                    labels = c("Baseline", "Improved"))+ 
  facet_wrap(~Trnsm,
             scales = "free")+
  # facet_wrap(~Trnsm)+
  labs(x = "Vaccine Efficacy",
       y = "Time to 2nd Room Infection")+
  theme_classic()+
  theme(panel.spacing.x = unit(1.25, "lines"),
        legend.title = element_text(face = "bold", size = 28),
        legend.text = element_text(face = "bold", size = 26),
        legend.title.align = 0.5,
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 24),
        axis.title = element_text(face = "bold", size = 28),
        axis.text = element_text(face = "bold", size = 26, color = "black"))
#1356 x 863




