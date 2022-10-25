#====================================================================================#
# Use for Sensativity Analysis
# Will be using a deterministic model sampling from distrubtions for model parameters
# Created by EK on 6/30/2022
#====================================================================================#
#### LIBRARIES ####
# library(tidyverse)
# library(deSolve)
# library(EnvStats)

#### Vacc Only ####
results <- NULL
n.itr <- 5000
# seed.num <- 616
seed.nums <- c(100:(100+n.itr))
for (i in c(1:n.itr)) {
  print(i)

  #### MODEL CONSTANTS ####
  #Hog constants
  beta_Hog <- rtri(n = 1, min = 0.001, max = 0.1, mode = (10/(1000*5)))
  beta_Hog_ind <- beta_Hog / 178 #(following Etbaigha et al 2018 paper)
  # beta_Hog_ind <- beta_Hog / 500
  Hog_Vacc_Eff <- 0.8
  u_Hog <- 0.00028     #natural death rate (following Etbaigha et al 2018 paper)
  u_inf_Hog <- 0 #0.1     #infected death rate (0 for now to problem solve)
  w_Hog <- 0 #1/180         #(following Etbaigha et al 2018 paper)  White et al shows a range from 56 - 112 days
  sigma_Hog <- abs(1 / rnorm(n = 1, mean = 2, sd = 1))   #Latent period h1n1 and h3n2
  delta_Hog <- abs(1 / rnorm(n = 1, mean = 5, sd = 1))   #Infectious period h1n1 and h3n2
  
  #Interspecies transmission
  ##applied the beta = Ro / N*D from the modeling book
  # Ro from the following website for swine to human infection 
  # https://bmcmedicine.biomedcentral.com/articles/10.1186/1741-7015-7-30#Sec6   (from super-strain figure) possible the high end of transmisilibty for this parameter
  # beta_S2WF <- (2.3 / (4002*5))  #  #Hogs + #Workforce * Duration of Hogs
  #Alt beta_S2WF from fair outbreak
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3557885/
  #assuming one hog got the first three confirmed cases sick and infectious period of 5 days Ro = 3/5 = 0.6
  beta_S2WF <- runif(1, min = 0.000001, max = 0.0001)
  #when Ro WF2S = Ro S2WF
  # beta_WF2S <- (2.3 / (4002*3))  #  #Hogs + #Workforce * Duration of Hogs
  # WF2S calculated from Swine Outbreak of Pandemic Influenza A Virus on a Canadian Research Farm Supports Human-to-Swine Transmission paper
  beta_WF2S <- runif(1, min = 0.000001, max = 0.0001)    #Hogs + #Workforce * Duration of WF
  #Human constants
  #applied the beta = Ro / N*D from the modeling book
  beta_WF2WF <-   rtri(n = 1, min = 0.0000001, max = 0.64, mode = 0.32) # http://m-hikari.com/ams/ams-2013/ams-41-44-2013/joseAMS41-44-2013.pdf 
  Hum_Vacc_Eff <- 0.0 
  u_WF <- 0     #Natural death rate of human (0 for now to problem solve)
  u_inf_WF <- 0   #Infected death rate of human (0 for now to problem solve)
  w_WF <- 0      #Recovery rate for human (0 for now to problem solve)
  sigma_WF <- abs(1 / rnorm(n = 1, mean = 2, sd = 1))   #Latent period for humans
  delta_WF <- abs(1 / rnorm(n = 1, mean = 3, sd = 1))   #Infectious Period for humans
  
  #### MODEL SETUP ####
  # Model Parameters
  TransmissionModel.par <- c(beta_Hog = beta_Hog, beta_Hog_ind = beta_Hog_ind, beta_S2WF = beta_S2WF, beta_WF2S = beta_WF2S,
                             w_Hog = w_Hog, sigma_Hog = sigma_Hog, delta_Hog = delta_Hog, u_Hog = u_Hog, u_inf_Hog = u_inf_Hog,
                             beta_WF2WF = beta_WF2WF,
                             w_WF = w_WF, sigma_WF = sigma_WF, delta_WF = delta_WF, u_WF = u_WF, u_inf_WF = u_inf_WF)
  #=====================================================================================================
  #Model Equations 
  #=====================================================================================================
  TransmissionModel.dyn <- function(t,var,par) {
    # Rename the parameters 
    beta_Hog <- par[1]
    beta_Hog_ind <- par[2]
    beta_S2WF <- par[3]
    beta_WF2S <- par[4]
    w_Hog <- par[5]
    sigma_Hog <- par[6]
    delta_Hog <- par[7]
    u_Hog <- par[8]
    u_inf_Hog<- par[9]
    beta_WF2WF <- par[10]
    w_WF <- par[11]
    sigma_WF <- par[12]
    delta_WF <- par[13]
    u_WF <- par[14]
    u_inf_WF <- par[15]
    
    S1 <- var[1]
    E1 <- var[2]
    I1 <- var[3]
    R1 <- var[4]
    Total_inf1 <- var[5] 
    True_Rec1 <- var[6]
    D1_Nat <- var[7]
    D1_Inf <- var[8]
    WF2Hog_Count1 <- var[9] 
    
    S2 <- var[10]
    E2 <- var[11]
    I2 <- var[12]
    R2 <- var[13]
    Total_inf2 <- var[14]
    True_Rec2 <- var[15]
    D2_Nat <- var[16]
    D2_Inf <- var[17]
    WF2Hog_Count2 <- var[18]
    
    S3 <- var[19]
    E3 <- var[20]
    I3 <- var[21]
    R3 <- var[22]
    Total_inf3 <- var[23]
    True_Rec3 <- var[24]
    D3_Nat <- var[25]
    D3_Inf <- var[26]
    WF2Hog_Count3 <- var[27]   
    
    S4 <- var[28]
    E4 <- var[29]
    I4 <- var[30]
    R4 <- var[31]
    Total_inf4 <- var[32]
    True_Rec4 <- var[33]
    D4_Nat <- var[34]
    D4_Inf <- var[35]
    WF2Hog_Count4 <- var[36]
    
    HS <- var[37]
    HE <- var[38]
    HI <- var[39] 
    HR <- var[40]
    H_DN <- var[41]
    H_DI <- var[42]
    Hog2WF_Count <- var[43] 
    
    # Calculate the derivatives
    #Room 1
    dS1 <- (w_Hog*R1) - (beta_Hog*I1*S1) - (beta_Hog_ind*(I2+I3+I4)*S1) - (beta_WF2S*HI*S1) - (u_Hog*S1)
    dE1 <-(beta_Hog*I1*S1) + (beta_Hog_ind*(I2+I3+I4)*S1) + (beta_WF2S*HI*S1) - sigma_Hog*E1 - (u_Hog*E1)
    dI1 <- (sigma_Hog*E1) - (delta_Hog*I1) - (u_inf_Hog*I1)
    dR1 <- (delta_Hog*I1) - (w_Hog*R1) - (u_Hog*R1)
    dTotal_inf1 <- (sigma_Hog*E1)
    dTotal_Rec1 <- (delta_Hog*I1)
    dD1_Nat <- (u_Hog*S1) + (u_Hog*E1) + (u_Hog*R1)
    dD1_Inf <- (u_inf_Hog*I1)
    dWF2Hog_Count1 <- (beta_WF2S*HI*S1)
    
    #Room 2
    dS2 <- (w_Hog*R2) - (beta_Hog*I2*S2) - (beta_Hog_ind*(I1+I3+I4)*S2) - (beta_WF2S*HI*S2) - (u_Hog*S2)
    dE2 <-(beta_Hog*I2*S2) + (beta_Hog_ind*(I1+I3+I4)*S2) + (beta_WF2S*HI*S2) - (sigma_Hog*E2) - (u_Hog*E2) 
    dI2 <- (sigma_Hog*E2) - (delta_Hog*I2) - (u_inf_Hog*I2)
    dR2 <- (delta_Hog*I2) - (w_Hog*R2) - (u_Hog*R2)
    dTotal_inf2 <- (sigma_Hog*E2)
    dTotal_Rec2 <- (delta_Hog*I2)
    dD2_Nat <- (u_Hog*S2) + (u_Hog*E2) + (u_Hog*R2)
    dD2_Inf <- (u_inf_Hog*I2)
    dWF2Hog_Count2 <- (beta_WF2S*HI*S2)
    
    #Room 3
    dS3 <- (w_Hog*R3) - (beta_Hog*I3*S3) - (beta_Hog_ind*(I1+I2+I4)*S3) - (beta_WF2S*HI*S3) - (u_Hog*S3)
    dE3 <-(beta_Hog*I3*S3) + (beta_Hog_ind*(I1+I2+I4)*S3) + (beta_WF2S*HI*S3) - (sigma_Hog*E3) - (u_Hog*E3)
    dI3 <- (sigma_Hog*E3) - (delta_Hog*I3) - (u_inf_Hog*I3)
    dR3 <- (delta_Hog*I3) - (w_Hog*R3) - (u_Hog*R3)
    dTotal_inf3 <- (sigma_Hog*E3)
    dTotal_Rec3 <- (delta_Hog*I3)
    dD3_Nat <- (u_Hog*S3) + (u_Hog*E3) + (u_Hog*R3)
    dD3_Inf <- (u_inf_Hog*I3)
    dWF2Hog_Count3 <- (beta_WF2S*HI*S3)
    
    #Room 4
    dS4 <- w_Hog*R4 - (beta_Hog*I4*S4) - (beta_Hog_ind*(I1+I2+I3)*S4) - (beta_WF2S*HI*S4) - (u_Hog*S4)
    dE4 <-(beta_Hog*I4*S4) + (beta_Hog_ind*(I1+I2+I3)*S4) + (beta_WF2S*HI*S4) - (sigma_Hog*E4) - (u_Hog*E4) 
    dI4 <- (sigma_Hog*E4) - (delta_Hog*I4) - (u_inf_Hog*I4)
    dR4 <- (delta_Hog*I4) - (w_Hog*R4) - (u_Hog*R4)
    dTotal_inf4 <- (sigma_Hog*E4)
    dTotal_Rec4 <- (delta_Hog*I4)
    dD4_Nat <- (u_Hog*S4) + (u_Hog*E4) + (u_Hog*R4)
    dD4_Inf <- (u_inf_Hog*I4)
    dWF2Hog_Count4 <- (beta_WF2S*HI*S4)
    
    #Humans
    dHS <- (w_WF*HR)- (beta_WF2WF*HI*HS) - (beta_S2WF*(I1+I2+I3+I4)*HS) - (u_WF*HS)
    dHE <- (beta_WF2WF*HI*HS) + (beta_S2WF*(I1+I2+I3+I4)*HS) - (sigma_WF*HE) - (u_WF*HE)
    dHI <- (sigma_WF*HE) - (delta_WF*HI) - (u_inf_WF*HI)
    dHR <- (delta_WF*HI) - (w_WF*HR) - u_WF*HR
    dH_DN <- (u_WF*HS) + (u_WF*HE) + (u_WF*HR)
    dH_DI <- (u_inf_WF*HI)
    dHog2WF_Count <- (beta_S2WF*(I1+I2+I3+I4)*HS)
    
    
    # Last instruction: return a list 
    
    return(list(c(dS1, dE1, dI1, dR1, dTotal_inf1, dTotal_Rec1, dD1_Nat, dD1_Inf, dWF2Hog_Count1,
                  dS2, dE2, dI2, dR2, dTotal_inf2, dTotal_Rec2, dD2_Nat, dD2_Inf, dWF2Hog_Count2,
                  dS3, dE3, dI3, dR3, dTotal_inf3, dTotal_Rec3, dD3_Nat, dD3_Inf, dWF2Hog_Count3,
                  dS4, dE4, dI4, dR4, dTotal_inf4, dTotal_Rec4, dD4_Nat, dD4_Inf, dWF2Hog_Count4,
                  dHS, dHE, dHI, dHR, dH_DN, dH_DI, dHog2WF_Count)))
    
  }
  
  #### MODEL ITERATIONS ####
  #### WEEK 0 TO 2 (Room 1 fill) ####
  Week0to2 <- 14
  #====================================================================================================
  # Model Initial Values  Weeks 0 to 2
  #===========================================================================
  S1_int <- (1-Hog_Vacc_Eff)*1000
  E1_int <- 0
  I1_int <- 0 #1
  R1_int <- Hog_Vacc_Eff*1000
  Total_inf1_int <- 0
  Total_Rec1_int <- 0
  D1_Nat_int <- 0
  D1_Inf_int <- 0
  WF2Hog_Count1_int <- 0
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- 0
  Total_Rec2_int <- 0
  D2_Nat_int <- 0
  D2_Inf_int <- 0
  WF2Hog_Count2_int <- 0
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <- 2
  HE_int <- 0
  HI_int <- 0
  HR_int <- 0
  H_DN_int <- 0
  H_DI_int <- 0
  Hog2WF_Count_int <- 0
  
  #====================================================================================================
  # Model Run  Weeks 0 to 2
  #====================================================================================================
  Week0to2.t <- seq(0,Week0to2,0.1)   
  Week0to2.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week0to2.sol <- lsoda(Week0to2.init, Week0to2.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 2 TO 4 (Room 2 Fill)####
  Week2to4 <- 28
  #====================================================================================================
  # Model Initial Values  Weeks 2 to 4
  #===========================================================================
  S1_int <- Week0to2.sol[nrow(Week0to2.sol),2]
  E1_int <- Week0to2.sol[nrow(Week0to2.sol),3]
  I1_int <- Week0to2.sol[nrow(Week0to2.sol),4]
  R1_int <- Week0to2.sol[nrow(Week0to2.sol),5]
  Total_inf1_int <-  Week0to2.sol[nrow(Week0to2.sol),6]
  Total_Rec1_int <-  Week0to2.sol[nrow(Week0to2.sol),7]
  D1_Nat_int <-  Week0to2.sol[nrow(Week0to2.sol),8]
  D1_Inf_int <-  Week0to2.sol[nrow(Week0to2.sol),9]
  WF2Hog_Count1_int <-  Week0to2.sol[nrow(Week0to2.sol),10]
  
  S2_int <- (1-Hog_Vacc_Eff)*1000
  E2_int <- 0
  I2_int <- 0 #1
  R2_int <- Hog_Vacc_Eff*1000
  Total_inf2_int <- 0
  Total_Rec2_int <- 0
  D2_Nat_int <- 0
  D2_Inf_int <- 0
  WF2Hog_Count2_int <- 0
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week0to2.sol[nrow(Week0to2.sol),38]
  HE_int <-  Week0to2.sol[nrow(Week0to2.sol),39]
  HI_int <-  Week0to2.sol[nrow(Week0to2.sol),40]
  HR_int <-  Week0to2.sol[nrow(Week0to2.sol),41]
  H_DN_int <-  Week0to2.sol[nrow(Week0to2.sol),42]
  H_DI_int <-  Week0to2.sol[nrow(Week0to2.sol),43]
  Hog2WF_Count_int <-  Week0to2.sol[nrow(Week0to2.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 2 to 4
  #====================================================================================================
  Week2to4.t <- seq(14.1,Week2to4,0.1)          
  Week2to4.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week2to4.sol <- lsoda(Week2to4.init, Week2to4.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 4 TO 6 (Room 3 Fill)####
  Week4to6 <- 42
  #====================================================================================================
  # Model Initial Values  Weeks 4 to 6
  #===========================================================================
  S1_int <- Week2to4.sol[nrow(Week2to4.sol),2]
  E1_int <- Week2to4.sol[nrow(Week2to4.sol),3]
  I1_int <- Week2to4.sol[nrow(Week2to4.sol),4]
  R1_int <- Week2to4.sol[nrow(Week2to4.sol),5]
  Total_inf1_int <- Week2to4.sol[nrow(Week2to4.sol),6]
  Total_Rec1_int <- Week2to4.sol[nrow(Week2to4.sol),7]
  D1_Nat_int <-  Week2to4.sol[nrow(Week2to4.sol),8]
  D1_Inf_int <-  Week2to4.sol[nrow(Week2to4.sol),9]
  WF2Hog_Count1_int <-  Week2to4.sol[nrow(Week2to4.sol),10]
  
  S2_int <- Week2to4.sol[nrow(Week2to4.sol),11]
  E2_int <- Week2to4.sol[nrow(Week2to4.sol),12]
  I2_int <- Week2to4.sol[nrow(Week2to4.sol),13]
  R2_int <- Week2to4.sol[nrow(Week2to4.sol),14]
  Total_inf2_int <- Week2to4.sol[nrow(Week2to4.sol),15]
  Total_Rec2_int <- Week2to4.sol[nrow(Week2to4.sol),16]
  D2_Nat_int <- Week2to4.sol[nrow(Week2to4.sol),17]
  D2_Inf_int <- Week2to4.sol[nrow(Week2to4.sol),18]
  WF2Hog_Count2_int <- Week2to4.sol[nrow(Week2to4.sol),19]
  
  S3_int <- (1-Hog_Vacc_Eff)*1000
  E3_int <- 0
  I3_int <- 0
  R3_int <-  Hog_Vacc_Eff*1000
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week2to4.sol[nrow(Week2to4.sol),38]
  HE_int <-  Week2to4.sol[nrow(Week2to4.sol),39]
  HI_int <-  Week2to4.sol[nrow(Week2to4.sol),40]
  HR_int <-  Week2to4.sol[nrow(Week2to4.sol),41]
  H_DN_int <-  Week2to4.sol[nrow(Week2to4.sol),42]
  H_DI_int <-  Week2to4.sol[nrow(Week2to4.sol),43]
  Hog2WF_Count_int <-  Week2to4.sol[nrow(Week2to4.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 4 to 6
  #====================================================================================================
  Week4to6.t <- seq(28.1,Week4to6,0.1)         
  Week4to6.init <-c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                    S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                    S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                    S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                    HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week4to6.sol <- lsoda(Week4to6.init, Week4to6.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 6 TO 8 (Room 4 Fill)####
  Week6to8 <- 56
  #====================================================================================================
  # Model Initial Values  Weeks 6 to 8
  #===========================================================================
  S1_int <- Week4to6.sol[nrow(Week4to6.sol),2]
  E1_int <- Week4to6.sol[nrow(Week4to6.sol),3]
  I1_int <- Week4to6.sol[nrow(Week4to6.sol),4]
  R1_int <- Week4to6.sol[nrow(Week4to6.sol),5]
  Total_inf1_int <- Week4to6.sol[nrow(Week4to6.sol),6]
  Total_Rec1_int <- Week4to6.sol[nrow(Week4to6.sol),7]
  D1_Nat_int <-  Week4to6.sol[nrow(Week4to6.sol),8]
  D1_Inf_int <-  Week4to6.sol[nrow(Week4to6.sol),9]
  WF2Hog_Count1_int <-  Week4to6.sol[nrow(Week4to6.sol),10]
  
  S2_int <- Week4to6.sol[nrow(Week4to6.sol),11]
  E2_int <- Week4to6.sol[nrow(Week4to6.sol),12]
  I2_int <- Week4to6.sol[nrow(Week4to6.sol),13]
  R2_int <- Week4to6.sol[nrow(Week4to6.sol),14]
  Total_inf2_int <- Week4to6.sol[nrow(Week4to6.sol),15]
  Total_Rec2_int <- Week4to6.sol[nrow(Week4to6.sol),16]
  D2_Nat_int <- Week4to6.sol[nrow(Week4to6.sol),17]
  D2_Inf_int <- Week4to6.sol[nrow(Week4to6.sol),18]
  WF2Hog_Count2_int <- Week4to6.sol[nrow(Week4to6.sol),19]
  
  S3_int <- Week4to6.sol[nrow(Week4to6.sol),20]
  E3_int <- Week4to6.sol[nrow(Week4to6.sol),21]
  I3_int <- Week4to6.sol[nrow(Week4to6.sol),22]
  R3_int <-  Week4to6.sol[nrow(Week4to6.sol),23]
  Total_inf3_int <- Week4to6.sol[nrow(Week4to6.sol),24]
  Total_Rec3_int <- Week4to6.sol[nrow(Week4to6.sol),25]
  D3_Nat_int <- Week4to6.sol[nrow(Week4to6.sol),26]
  D3_Inf_int <- Week4to6.sol[nrow(Week4to6.sol),27]
  WF2Hog_Count3_int <- Week4to6.sol[nrow(Week4to6.sol),28]
  
  S4_int <- (1-Hog_Vacc_Eff)*999
  E4_int <- 0
  I4_int <- 1
  R4_int <- Hog_Vacc_Eff*999
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week4to6.sol[nrow(Week4to6.sol),38]
  HE_int <-  Week4to6.sol[nrow(Week4to6.sol),39]
  HI_int <-  Week4to6.sol[nrow(Week4to6.sol),40]
  HR_int <-  Week4to6.sol[nrow(Week4to6.sol),41]
  H_DN_int <-  Week4to6.sol[nrow(Week4to6.sol),42]
  H_DI_int <-  Week4to6.sol[nrow(Week4to6.sol),43]
  Hog2WF_Count_int <-  Week4to6.sol[nrow(Week4to6.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 6 to 8
  #====================================================================================================
  Week6to8.t <- seq(42.1,Week6to8,0.1)           
  Week6to8.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week6to8.sol <- lsoda(Week6to8.init, Week6to8.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 8 TO 23 (Room 1 Exit) ####
  Week8to23 <- 161
  #====================================================================================================
  # Model Initial Values  Weeks 8 to 23
  #===========================================================================
  S1_int <- Week6to8.sol[nrow(Week6to8.sol),2]
  E1_int <- Week6to8.sol[nrow(Week6to8.sol),3]
  I1_int <- Week6to8.sol[nrow(Week6to8.sol),4]
  R1_int <- Week6to8.sol[nrow(Week6to8.sol),5]
  Total_inf1_int <- Week6to8.sol[nrow(Week6to8.sol),6]
  Total_Rec1_int <- Week6to8.sol[nrow(Week6to8.sol),7]
  D1_Nat_int <-  Week6to8.sol[nrow(Week6to8.sol),8]
  D1_Inf_int <-  Week6to8.sol[nrow(Week6to8.sol),9]
  WF2Hog_Count1_int <-  Week6to8.sol[nrow(Week6to8.sol),10]
  
  
  S2_int <- Week6to8.sol[nrow(Week6to8.sol),11]
  E2_int <- Week6to8.sol[nrow(Week6to8.sol),12]
  I2_int <- Week6to8.sol[nrow(Week6to8.sol),13]
  R2_int <- Week6to8.sol[nrow(Week6to8.sol),14]
  Total_inf2_int <- Week6to8.sol[nrow(Week6to8.sol),15]
  Total_Rec2_int <- Week6to8.sol[nrow(Week6to8.sol),16]
  D2_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),17]
  D2_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),18]
  WF2Hog_Count2_int <- Week6to8.sol[nrow(Week6to8.sol),19]
  
  S3_int <- Week6to8.sol[nrow(Week6to8.sol),20]
  E3_int <- Week6to8.sol[nrow(Week6to8.sol),21]
  I3_int <- Week6to8.sol[nrow(Week6to8.sol),22]
  R3_int <-  Week6to8.sol[nrow(Week6to8.sol),23]
  Total_inf3_int <- Week6to8.sol[nrow(Week6to8.sol),24]
  Total_Rec3_int <- Week6to8.sol[nrow(Week6to8.sol),25]
  D3_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),26]
  D3_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),27]
  WF2Hog_Count3_int <- Week6to8.sol[nrow(Week6to8.sol),28]
  
  S4_int <- Week6to8.sol[nrow(Week6to8.sol),29]
  E4_int <- Week6to8.sol[nrow(Week6to8.sol),30]
  I4_int <- Week6to8.sol[nrow(Week6to8.sol),31]
  R4_int <- Week6to8.sol[nrow(Week6to8.sol),32]
  Total_inf4_int <- Week6to8.sol[nrow(Week6to8.sol),33]
  Total_Rec4_int <- Week6to8.sol[nrow(Week6to8.sol),34]
  D4_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),35]
  D4_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),36]
  WF2Hog_Count4_int <- Week6to8.sol[nrow(Week6to8.sol),37]
  
  HS_int <-  Week6to8.sol[nrow(Week6to8.sol),38]
  HE_int <-  Week6to8.sol[nrow(Week6to8.sol),39]
  HI_int <-  Week6to8.sol[nrow(Week6to8.sol),40]
  HR_int <-  Week6to8.sol[nrow(Week6to8.sol),41]
  H_DN_int <-  Week6to8.sol[nrow(Week6to8.sol),42]
  H_DI_int <-  Week6to8.sol[nrow(Week6to8.sol),43]
  Hog2WF_Count_int <-  Week6to8.sol[nrow(Week6to8.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 8 to 23
  #====================================================================================================
  Week8to23.t <- seq(56.1,Week8to23,0.1)          
  Week8to23.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                      S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                      S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                      S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                      HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week8to23.sol <- lsoda(Week8to23.init, Week8to23.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 23 TO 25 (Room 2 Exit) ####
  Week23to25 <- 175
  #====================================================================================================
  # Model Initial Values  Weeks 23 to 25
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week8to23.sol[nrow(Week8to23.sol),6]
  Total_Rec1_int <- Week8to23.sol[nrow(Week8to23.sol),7]
  D1_Nat_int <-  Week8to23.sol[nrow(Week8to23.sol),8]
  D1_Inf_int <-  Week8to23.sol[nrow(Week8to23.sol),9]
  WF2Hog_Count1_int <-  Week8to23.sol[nrow(Week8to23.sol),10]
  
  
  S2_int <- Week8to23.sol[nrow(Week8to23.sol),11]
  E2_int <- Week8to23.sol[nrow(Week8to23.sol),12]
  I2_int <- Week8to23.sol[nrow(Week8to23.sol),13]
  R2_int <- Week8to23.sol[nrow(Week8to23.sol),14]
  Total_inf2_int <- Week8to23.sol[nrow(Week8to23.sol),15]
  Total_Rec2_int <- Week8to23.sol[nrow(Week8to23.sol),16]
  D2_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),17]
  D2_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),18]
  WF2Hog_Count2_int <- Week8to23.sol[nrow(Week8to23.sol),19]
  
  S3_int <- Week8to23.sol[nrow(Week8to23.sol),20]
  E3_int <- Week8to23.sol[nrow(Week8to23.sol),21]
  I3_int <- Week8to23.sol[nrow(Week8to23.sol),22]
  R3_int <-  Week8to23.sol[nrow(Week8to23.sol),23]
  Total_inf3_int <- Week8to23.sol[nrow(Week8to23.sol),24]
  Total_Rec3_int <- Week8to23.sol[nrow(Week8to23.sol),25]
  D3_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),26]
  D3_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),27]
  WF2Hog_Count3_int <- Week8to23.sol[nrow(Week8to23.sol),28]
  
  S4_int <- Week8to23.sol[nrow(Week8to23.sol),29]
  E4_int <- Week8to23.sol[nrow(Week8to23.sol),30]
  I4_int <- Week8to23.sol[nrow(Week8to23.sol),31]
  R4_int <- Week8to23.sol[nrow(Week8to23.sol),32]
  Total_inf4_int <- Week8to23.sol[nrow(Week8to23.sol),33]
  Total_Rec4_int <- Week8to23.sol[nrow(Week8to23.sol),34]
  D4_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),35]
  D4_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),36]
  WF2Hog_Count4_int <- Week8to23.sol[nrow(Week8to23.sol),37]
  
  HS_int <-  Week8to23.sol[nrow(Week8to23.sol),38]
  HE_int <-  Week8to23.sol[nrow(Week8to23.sol),39]
  HI_int <-  Week8to23.sol[nrow(Week8to23.sol),40]
  HR_int <-  Week8to23.sol[nrow(Week8to23.sol),41]
  H_DN_int <-  Week8to23.sol[nrow(Week8to23.sol),42]
  H_DI_int <-  Week8to23.sol[nrow(Week8to23.sol),43]
  Hog2WF_Count_int <-  Week8to23.sol[nrow(Week8to23.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 23 to 25
  #====================================================================================================
  Week23to25.t <- seq(161.1,Week23to25,0.1)         
  Week23to25.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week23to25.sol <- lsoda(Week23to25.init, Week23to25.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 25 TO 27 (Room 3 Exit) ####
  Week25to27 <- 189
  #====================================================================================================
  # Model Initial Values  Weeks 25 to 27
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week23to25.sol[nrow(Week23to25.sol),6]
  Total_Rec1_int <- Week23to25.sol[nrow(Week23to25.sol),7]
  D1_Nat_int <-  Week23to25.sol[nrow(Week23to25.sol),8]
  D1_Inf_int <-  Week23to25.sol[nrow(Week23to25.sol),9]
  WF2Hog_Count1_int <-  Week23to25.sol[nrow(Week23to25.sol),10]
  
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- Week23to25.sol[nrow(Week23to25.sol),15]
  Total_Rec2_int <- Week23to25.sol[nrow(Week23to25.sol),16]
  D2_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),17]
  D2_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),18]
  WF2Hog_Count2_int <- Week23to25.sol[nrow(Week23to25.sol),19]
  
  S3_int <- Week23to25.sol[nrow(Week23to25.sol),20]
  E3_int <- Week23to25.sol[nrow(Week23to25.sol),21]
  I3_int <- Week23to25.sol[nrow(Week23to25.sol),22]
  R3_int <-  Week23to25.sol[nrow(Week23to25.sol),23]
  Total_inf3_int <- Week23to25.sol[nrow(Week23to25.sol),24]
  Total_Rec3_int <- Week23to25.sol[nrow(Week23to25.sol),25]
  D3_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),26]
  D3_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),27]
  WF2Hog_Count3_int <- Week23to25.sol[nrow(Week23to25.sol),28]
  
  S4_int <- Week23to25.sol[nrow(Week23to25.sol),29]
  E4_int <- Week23to25.sol[nrow(Week23to25.sol),30]
  I4_int <- Week23to25.sol[nrow(Week23to25.sol),31]
  R4_int <- Week23to25.sol[nrow(Week23to25.sol),32]
  Total_inf4_int <- Week23to25.sol[nrow(Week23to25.sol),33]
  Total_Rec4_int <- Week23to25.sol[nrow(Week23to25.sol),34]
  D4_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),35]
  D4_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),36]
  WF2Hog_Count4_int <- Week23to25.sol[nrow(Week23to25.sol),37]
  
  HS_int <-  Week23to25.sol[nrow(Week23to25.sol),38]
  HE_int <-  Week23to25.sol[nrow(Week23to25.sol),39]
  HI_int <-  Week23to25.sol[nrow(Week23to25.sol),40]
  HR_int <-  Week23to25.sol[nrow(Week23to25.sol),41]
  H_DN_int <-  Week23to25.sol[nrow(Week23to25.sol),42]
  H_DI_int <-  Week23to25.sol[nrow(Week23to25.sol),43]
  Hog2WF_Count_int <-  Week23to25.sol[nrow(Week23to25.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 25 to 27
  #====================================================================================================
  Week25to27.t <- seq(175.1,Week25to27,0.1)          
  Week25to27.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week25to27.sol <- lsoda(Week25to27.init, Week25to27.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 27 TO 29 (Room 4 Exit) ####
  Week27to29 <- 203
  #====================================================================================================
  # Model Initial Values  Weeks 27to 29
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week25to27.sol[nrow(Week25to27.sol),6]
  Total_Rec1_int <- Week25to27.sol[nrow(Week25to27.sol),7]
  D1_Nat_int <-  Week25to27.sol[nrow(Week25to27.sol),8]
  D1_Inf_int <-  Week25to27.sol[nrow(Week25to27.sol),9]
  WF2Hog_Count1_int <-  Week25to27.sol[nrow(Week25to27.sol),10]
  
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- Week25to27.sol[nrow(Week25to27.sol),15]
  Total_Rec2_int <- Week25to27.sol[nrow(Week25to27.sol),16]
  D2_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),17]
  D2_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),18]
  WF2Hog_Count2_int <- Week25to27.sol[nrow(Week25to27.sol),19]
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- Week25to27.sol[nrow(Week25to27.sol),24]
  Total_Rec3_int <- Week25to27.sol[nrow(Week25to27.sol),25]
  D3_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),26]
  D3_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),27]
  WF2Hog_Count3_int <- Week25to27.sol[nrow(Week25to27.sol),28]
  
  S4_int <- Week25to27.sol[nrow(Week25to27.sol),29]
  E4_int <- Week25to27.sol[nrow(Week25to27.sol),30]
  I4_int <- Week25to27.sol[nrow(Week25to27.sol),31]
  R4_int <- Week25to27.sol[nrow(Week25to27.sol),32]
  Total_inf4_int <- Week25to27.sol[nrow(Week25to27.sol),33]
  Total_Rec4_int <- Week25to27.sol[nrow(Week25to27.sol),34]
  D4_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),35]
  D4_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),36]
  WF2Hog_Count4_int <- Week25to27.sol[nrow(Week25to27.sol),37]
  
  HS_int <-  Week25to27.sol[nrow(Week25to27.sol),38]
  HE_int <-  Week25to27.sol[nrow(Week25to27.sol),39]
  HI_int <-  Week25to27.sol[nrow(Week25to27.sol),40]
  HR_int <-  Week25to27.sol[nrow(Week25to27.sol),41]
  H_DN_int <-  Week25to27.sol[nrow(Week25to27.sol),42]
  H_DI_int <-  Week25to27.sol[nrow(Week25to27.sol),43]
  Hog2WF_Count_int <-  Week25to27.sol[nrow(Week25to27.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 27 to 29
  #====================================================================================================
  Week27to29.t <- seq(189.1,Week27to29,0.1)            
  Week27to29.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week27to29.sol <- lsoda(Week27to29.init, Week27to29.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### COMBINING THE DATASETS ####
  Week0to2.sol.df <- as.data.frame(Week0to2.sol)
  colnames(Week0to2.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week2to4.sol.df <- as.data.frame(Week2to4.sol)
  colnames(Week2to4.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week4to6.sol.df <- as.data.frame(Week4to6.sol)
  colnames(Week4to6.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week6to8.sol.df <- as.data.frame(Week6to8.sol)
  colnames(Week6to8.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week8to23.sol.df <- as.data.frame(Week8to23.sol)
  colnames(Week8to23.sol.df) <- c("time",
                                  "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                  "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                  "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                  "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                  "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week23to25.sol.df <- as.data.frame(Week23to25.sol)
  colnames(Week23to25.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week25to27.sol.df <- as.data.frame(Week25to27.sol)
  colnames(Week25to27.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week27to29.sol.df <- as.data.frame(Week27to29.sol)
  colnames(Week27to29.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  
  fulldata.df <- rbind(Week0to2.sol.df, Week2to4.sol.df, Week4to6.sol.df, Week6to8.sol.df, 
                       Week8to23.sol.df, Week23to25.sol.df, Week25to27.sol.df, Week27to29.sol.df)
  
  ## Is there a WF infection
  HI_inf <- fulldata.df %>%
    summarise(total_I = sum(HI)) %>%
    mutate(HI_Inf = ifelse(total_I > 0, 1, 0)) %>%
    select(HI_Inf) %>% as.numeric()
  #Time to 1st WF infection
  Time2WFinf0 <- fulldata.df %>%
    arrange(time) %>%
    filter(HI > 0) %>%
    filter(row_number() == 1) %>%
    mutate(TimeSineInfIntro = time - 42) %>%
    select(TimeSineInfIntro) %>% as.numeric()
  
  Time2WFinf.5 <- fulldata.df %>%
    arrange(time) %>%
    filter(HI > 0.5) %>%
    filter(row_number() == 1) %>%
    mutate(TimeSineInfIntro = time - 42) %>%
    select(TimeSineInfIntro) %>% as.numeric()
  ## Time to Second Room Infection
  Timeto2room <- fulldata.df %>% 
    summarise(t_cross = time[which(E1 > 1 | 
                                     E2 > 1 |
                                     E3 > 1)]) %>% 
    summarise(min_t_cross = min(t_cross)) %>% as.numeric()
  ## Total Hogs infected
  TotHogsInf <- fulldata.df %>% 
    #select(t, iter, Total_inf1, Total_inf2, Total_inf3, Total_inf4) %>% 
    mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
           Tot_rec = Total_Rec1 + Total_Rec2 + Total_Rec3 + Total_Rec4) %>% 
    summarise(Tot_inf2 = max(Tot_inf),
              Tot_rec2 = max(Tot_rec))
  ## Minimum Time to Peak Infection
  Time2PeakInf <- fulldata.df %>% 
    mutate(S_tot = S1+S2+S3+S4,                             
           E_tot = E1+E2+E3+E4,
           I_tot = I1+I2+I3+I4,
           R_tot = R1+R2+R3+R4) %>%                                     
    summarise(max_I = max(I_tot),
              t_peak = time[which(I_tot == max(I_tot))]) %>% 
    summarise(min_t_peak = min(t_peak)) %>% as.numeric()
  ##Ro
  Ro <- fulldata.df %>% 
    arrange(time) %>% 
    mutate(Tot_rec = Total_Rec1 + Total_Rec2 + Total_Rec3 + Total_Rec4) %>%  
    filter(row_number() == n()) %>% 
    select(Tot_rec) %>% 
    mutate(Ro = ((log((4000 - Tot_rec)/4000) - log((4000 - (Hog_Vacc_Eff*4000))/4000)) /
                   ((((4000 - Tot_rec) / 4000) - (4000 - (Hog_Vacc_Eff*4000)) / 4000)))) %>% 
    select(Ro) %>% as.numeric()
  
  #Create Summary data frame
  fulldata.df.sum <- data.frame(iter = as.factor(i),
                                HI_inf = HI_inf,
                                Time2WFinf0 = Time2WFinf0,
                                Time2WFinf.5 = Time2WFinf.5,
                                Timeto2room = Timeto2room,
                                TotHogsInf = TotHogsInf$Tot_inf2,
                                TotHogsRec = TotHogsInf$Tot_rec2,
                                Time2PeakInf = Time2PeakInf,
                                Ro = Ro,
                                beta_Hog = beta_Hog,
                                beta_Hog_ind = beta_Hog_ind,
                                beta_S2WF = beta_S2WF,
                                beta_WF2S = beta_WF2S,
                                sigma_Hog = sigma_Hog,
                                delta_Hog = delta_Hog,
                                beta_WF2WF = beta_WF2WF,
                                sigma_WF = sigma_WF,
                                delta_WF = delta_WF)
  
  results <- rbind(results, fulldata.df.sum)
  
  rm(list=ls()[! ls() %in% c("results", "n.itr")])  #Deletes everything from environment expect growing results dataframe
}

saveRDS(results, file = "Current code/Control Meausres/Sensativity Analysis/80VaccEff_Sum_SensAnlaysis.rds")

#### Quarantine ####
results <- NULL
n.itr <- 5000
# seed.num <- 616
seed.nums <- c(100:(100+n.itr))
for (i in c(1:n.itr)) {
  print(i)
  
  #### MODEL CONSTANTS ####
  #Hog constants
  beta_Hog <- rtri(n = 1, min = 0.001, max = 0.1, mode = (10/(1000*5)))
  beta_Hog_ind <- beta_Hog / 178 #(following Etbaigha et al 2018 paper)
  # beta_Hog_ind <- beta_Hog / 500
  Hog_Vacc_Eff <- 0
  u_Hog <- 0.00028     #natural death rate (following Etbaigha et al 2018 paper)
  u_inf_Hog <- 0 #0.1     #infected death rate (0 for now to problem solve)
  w_Hog <- 0 #1/180         #(following Etbaigha et al 2018 paper)  White et al shows a range from 56 - 112 days
  sigma_Hog <- abs(1 / rnorm(n = 1, mean = 2, sd = 1))   #Latent period h1n1 and h3n2
  delta_Hog <- abs(1 / rnorm(n = 1, mean = 5, sd = 1))   #Infectious period h1n1 and h3n2
  qrate <- 0.5
  beta_qHog_ind <- 0
  q_cap <- 10
  
  
  #Interspecies transmission
  ##applied the beta = Ro / N*D from the modeling book
  # Ro from the following website for swine to human infection 
  # https://bmcmedicine.biomedcentral.com/articles/10.1186/1741-7015-7-30#Sec6   (from super-strain figure) possible the high end of transmisilibty for this parameter
  # beta_S2WF <- (2.3 / (4002*5))  #  #Hogs + #Workforce * Duration of Hogs
  #Alt beta_S2WF from fair outbreak
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3557885/
  #assuming one hog got the first three confirmed cases sick and infectious period of 5 days Ro = 3/5 = 0.6
  beta_S2WF <- runif(1, min = 0.000001, max = 0.0001)
  beta_qS2WF <- 0
  #when Ro WF2S = Ro S2WF
  # beta_WF2S <- (2.3 / (4002*3))  #  #Hogs + #Workforce * Duration of Hogs
  # WF2S calculated from Swine Outbreak of Pandemic Influenza A Virus on a Canadian Research Farm Supports Human-to-Swine Transmission paper
  beta_WF2S <- runif(1, min = 0.000001, max = 0.0001)    #Hogs + #Workforce * Duration of WF
  #Human constants
  #applied the beta = Ro / N*D from the modeling book
  beta_WF2WF <-   rtri(n = 1, min = 0.0000001, max = 0.64, mode = 0.32) # http://m-hikari.com/ams/ams-2013/ams-41-44-2013/joseAMS41-44-2013.pdf 
  Hum_Vacc_Eff <- 0.0 
  u_WF <- 0     #Natural death rate of human (0 for now to problem solve)
  u_inf_WF <- 0   #Infected death rate of human (0 for now to problem solve)
  w_WF <- 0      #Recovery rate for human (0 for now to problem solve)
  sigma_WF <- abs(1 / rnorm(n = 1, mean = 2, sd = 1))   #Latent period for humans
  delta_WF <- abs(1 / rnorm(n = 1, mean = 3, sd = 1))   #Infectious Period for humans
  
  #### MODEL SETUP ####
  # Model Parameters
  TransmissionModel.par <- c(beta_Hog = beta_Hog, beta_Hog_ind = beta_Hog_ind, beta_qHog_ind = beta_qHog_ind, 
                             beta_S2WF = beta_S2WF, beta_qS2WF = beta_qS2WF, beta_WF2S = beta_WF2S,
                             w_Hog = w_Hog, sigma_Hog = sigma_Hog, delta_Hog = delta_Hog, qrate = qrate, 
                             u_Hog = u_Hog, u_inf_Hog = u_inf_Hog, beta_WF2WF = beta_WF2WF,
                             w_WF = w_WF, sigma_WF = sigma_WF, delta_WF = delta_WF, u_WF = u_WF, u_inf_WF = u_inf_WF,
                             q_cap = q_cap)
  #=====================================================================================================
  #Model Equations 
  #=====================================================================================================
  TransmissionModel.dyn <- function(t,var,par) {
    # Rename the parameters 
    beta_Hog <- par[1]
    beta_Hog_ind <- par[2]
    beta_qHog_ind <- par[3]
    beta_S2WF <- par[4]
    beta_qS2WF <- par[5]
    beta_WF2S <- par[6]
    w_Hog <- par[7]
    sigma_Hog <- par[8]
    delta_Hog <- par[9]
    qrate <- par[10]
    u_Hog <- par[11]
    u_inf_Hog<- par[12]
    beta_WF2WF <- par[13]
    w_WF <- par[14]
    sigma_WF <- par[15]
    delta_WF <- par[16]
    u_WF <- par[17]
    u_inf_WF <- par[18]
    qcap <- par[19]
    
    S1 <- var[1]
    E1 <- var[2]
    Q1 <- var[3]
    I1 <- var[4]
    R1 <- var[5]
    Total_inf1 <- var[6] 
    True_Rec1 <- var[7]
    D1_Nat <- var[8]
    D1_Inf <- var[9]
    WF2Hog_Count1 <- var[10] 
    
    S2 <- var[11]
    E2 <- var[12]
    Q2 <- var[13]
    I2 <- var[14]
    R2 <- var[15]
    Total_inf2 <- var[16]
    True_Rec2 <- var[17]
    D2_Nat <- var[18]
    D2_Inf <- var[19]
    WF2Hog_Count2 <- var[20]
    
    S3 <- var[21]
    E3 <- var[22]
    Q3 <- var[23]
    I3 <- var[24]
    R3 <- var[25]
    Total_inf3 <- var[26]
    True_Rec3 <- var[27]
    D3_Nat <- var[28]
    D3_Inf <- var[29]
    WF2Hog_Count3 <- var[30]   
    
    S4 <- var[31]
    E4 <- var[32]
    Q4 <- var[33]
    I4 <- var[34]
    R4 <- var[35]
    Total_inf4 <- var[36]
    True_Rec4 <- var[37]
    D4_Nat <- var[38]
    D4_Inf <- var[39]
    WF2Hog_Count4 <- var[40]
    
    HS <- var[41]
    HE <- var[42]
    HI <- var[43] 
    HR <- var[44]
    H_DN <- var[45]
    H_DI <- var[46]
    Hog2WF_Count <- var[47] 
    
    # Calculate the derivatives
    #Room 1
    dS1 <- (w_Hog*R1) - (beta_Hog*I1*S1) - (beta_Hog_ind*(I2+I3+I4)*S1) - (beta_qHog_ind*(Q1+Q2+Q3+Q4)*S1) - (beta_WF2S*HI*S1) - (u_Hog*S1)
    dE1 <-(beta_Hog*I1*S1) + (beta_Hog_ind*(I2+I3+I4)*S1) + (beta_qHog_ind*(Q1+Q2+Q3+Q4)*S1) + (beta_WF2S*HI*S1) - sigma_Hog*E1 - (u_Hog*E1)
    dI1 <- (sigma_Hog*E1) - (delta_Hog*I1) - (ifelse(Q1<q_cap, (qrate*I1), 0)) - (u_inf_Hog*I1)
    dQ1 <- (ifelse(Q1<q_cap, (qrate*I1), 0)) - (delta_Hog*Q1) - (u_inf_Hog*Q1)
    dR1 <- (delta_Hog*I1) + (delta_Hog*Q1) - (w_Hog*R1) - (u_Hog*R1)
    dTotal_inf1 <- (sigma_Hog*E1)
    dTotal_Rec1 <- (delta_Hog*I1) + (delta_Hog*Q1)
    dD1_Nat <- (u_Hog*S1) + (u_Hog*E1) + (u_Hog*R1)
    dD1_Inf <- (u_inf_Hog*I1)
    dWF2Hog_Count1 <- (beta_WF2S*HI*S1)
    
    #Room 2
    dS2 <- (w_Hog*R2) - (beta_Hog*I2*S2) - (beta_Hog_ind*(I1+I3+I4)*S2) - (beta_qHog_ind*(Q1+Q2+Q3+Q4)*S2) - (beta_WF2S*HI*S2) - (u_Hog*S2)
    dE2 <-(beta_Hog*I2*S2) + (beta_Hog_ind*(I1+I3+I4)*S2) + (beta_qHog_ind*(Q1+Q2+Q3+Q4)*S2) + (beta_WF2S*HI*S2) - (sigma_Hog*E2) - (u_Hog*E2) 
    dI2 <- (sigma_Hog*E2) - (delta_Hog*I2) - (ifelse(Q2<q_cap, (qrate*I2), 0)) - (u_inf_Hog*I2)
    dQ2 <- (ifelse(Q2<q_cap, (qrate*I2), 0)) - (delta_Hog*Q2) - (u_inf_Hog*Q2)
    dR2 <- (delta_Hog*I2) + (delta_Hog*Q2) - (w_Hog*R2) - (u_Hog*R2)
    dTotal_inf2 <- (sigma_Hog*E2)
    dTotal_Rec2 <- (delta_Hog*I2) + (delta_Hog*Q2)
    dD2_Nat <- (u_Hog*S2) + (u_Hog*E2) + (u_Hog*R2)
    dD2_Inf <- (u_inf_Hog*I2)
    dWF2Hog_Count2 <- (beta_WF2S*HI*S2)
    
    #Room 3
    dS3 <- (w_Hog*R3) - (beta_Hog*I3*S3) - (beta_Hog_ind*(I1+I2+I4)*S3) - (beta_qHog_ind*(Q1+Q2+Q3+Q4)*S3) - (beta_WF2S*HI*S3) - (u_Hog*S3)
    dE3 <-(beta_Hog*I3*S3) + (beta_Hog_ind*(I1+I2+I4)*S3) + (beta_qHog_ind*(Q1+Q2+Q3+Q4)*S3) + (beta_WF2S*HI*S3) - (sigma_Hog*E3) - (u_Hog*E3)
    dI3 <- (sigma_Hog*E3) - (delta_Hog*I3) - (ifelse(Q3<q_cap, (qrate*I3), 0)) - (u_inf_Hog*I3)
    dQ3 <- (ifelse(Q3<q_cap, (qrate*I3), 0)) - (delta_Hog*Q3) - (u_inf_Hog*Q3)
    dR3 <- (delta_Hog*I3) + (delta_Hog*Q3) - (w_Hog*R3) - (u_Hog*R3)
    dTotal_inf3 <- (sigma_Hog*E3)
    dTotal_Rec3 <- (delta_Hog*I3) + (delta_Hog*Q3)
    dD3_Nat <- (u_Hog*S3) + (u_Hog*E3) + (u_Hog*R3)
    dD3_Inf <- (u_inf_Hog*I3)
    dWF2Hog_Count3 <- (beta_WF2S*HI*S3)
    
    #Room 4
    dS4 <- w_Hog*R4 - (beta_Hog*I4*S4) - (beta_Hog_ind*(I1+I2+I3)*S4) - (beta_qHog_ind*(Q1+Q2+Q3+Q4)*S4) - (beta_WF2S*HI*S4) - (u_Hog*S4)
    dE4 <-(beta_Hog*I4*S4) + (beta_Hog_ind*(I1+I2+I3)*S4) + (beta_qHog_ind*(Q1+Q2+Q3+Q4)*S4) + (beta_WF2S*HI*S4) - (sigma_Hog*E4) - (u_Hog*E4) 
    dI4 <- (sigma_Hog*E4) - (delta_Hog*I4) - (ifelse(Q4<q_cap, (qrate*I4), 0)) - (u_inf_Hog*I4)
    dQ4 <- (ifelse(Q4<q_cap, (qrate*I4), 0)) - (delta_Hog*Q4) - (u_inf_Hog*Q4)
    dR4 <- (delta_Hog*I4) + (delta_Hog*Q4) - (w_Hog*R4) - (u_Hog*R4)
    dTotal_inf4 <- (sigma_Hog*E4)
    dTotal_Rec4 <- (delta_Hog*I4) + (delta_Hog*Q4)
    dD4_Nat <- (u_Hog*S4) + (u_Hog*E4) + (u_Hog*R4)
    dD4_Inf <- (u_inf_Hog*I4)
    dWF2Hog_Count4 <- (beta_WF2S*HI*S4)
    
    #Humans
    dHS <- (w_WF*HR)- (beta_WF2WF*HI*HS) - (beta_S2WF*(I1+I2+I3+I4)*HS) - (beta_qS2WF*(Q1+Q2+Q3+Q4)*HS) - (u_WF*HS)
    dHE <- (beta_WF2WF*HI*HS) + (beta_S2WF*(I1+I2+I3+I4)*HS) + (beta_qS2WF*(Q1+Q2+Q3+Q4)*HS) - (sigma_WF*HE) - (u_WF*HE)
    dHI <- (sigma_WF*HE) - (delta_WF*HI) - (u_inf_WF*HI)
    dHR <- (delta_WF*HI) - (w_WF*HR) - u_WF*HR
    dH_DN <- (u_WF*HS) + (u_WF*HE) + (u_WF*HR)
    dH_DI <- (u_inf_WF*HI)
    dHog2WF_Count <- (beta_S2WF*(I1+I2+I3+I4)*HS)
    
    
    # Last instruction: return a list 
    
    return(list(c(dS1, dE1, dI1, dQ1, dR1, dTotal_inf1, dTotal_Rec1, dD1_Nat, dD1_Inf, dWF2Hog_Count1,
                  dS2, dE2, dI2, dQ2, dR2, dTotal_inf2, dTotal_Rec2, dD2_Nat, dD2_Inf, dWF2Hog_Count2,
                  dS3, dE3, dI3, dQ3, dR3, dTotal_inf3, dTotal_Rec3, dD3_Nat, dD3_Inf, dWF2Hog_Count3,
                  dS4, dE4, dI4, dQ4, dR4, dTotal_inf4, dTotal_Rec4, dD4_Nat, dD4_Inf, dWF2Hog_Count4,
                  dHS, dHE, dHI, dHR, dH_DN, dH_DI, dHog2WF_Count)))
    
  }
  
  #### MODEL ITERATIONS ####
  #### WEEK 0 TO 2 (Room 1 fill) ####
  Week0to2 <- 14
  #====================================================================================================
  # Model Initial Values  Weeks 0 to 2
  #===========================================================================
  S1_int <- (1-Hog_Vacc_Eff)*1000
  E1_int <- 0
  I1_int <- 0 
  Q1_int <- 0
  R1_int <- Hog_Vacc_Eff*1000
  Total_inf1_int <- 0
  Total_Rec1_int <- 0
  D1_Nat_int <- 0
  D1_Inf_int <- 0
  WF2Hog_Count1_int <- 0
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  Q2_int <- 0
  R2_int <- 0
  Total_inf2_int <- 0
  Total_Rec2_int <- 0
  D2_Nat_int <- 0
  D2_Inf_int <- 0
  WF2Hog_Count2_int <- 0
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  Q3_int <- 0
  R3_int <- 0
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  Q4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <- 2
  HE_int <- 0
  HI_int <- 0
  HR_int <- 0
  H_DN_int <- 0
  H_DI_int <- 0
  Hog2WF_Count_int <- 0
  
  #====================================================================================================
  # Model Run  Weeks 0 to 2
  #====================================================================================================
  Week0to2.t <- seq(0,Week0to2,0.1)   
  Week0to2.init <- c(S1_int, E1_int, I1_int, Q1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, Q2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, Q3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, Q4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week0to2.sol <- lsoda(Week0to2.init, Week0to2.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 2 TO 4 (Room 2 Fill)####
  Week2to4 <- 28
  #====================================================================================================
  # Model Initial Values  Weeks 2 to 4
  #===========================================================================
  S1_int <- Week0to2.sol[nrow(Week0to2.sol),2]
  E1_int <- Week0to2.sol[nrow(Week0to2.sol),3]
  I1_int <- Week0to2.sol[nrow(Week0to2.sol),4]
  Q1_int <- Week0to2.sol[nrow(Week0to2.sol),5]
  R1_int <- Week0to2.sol[nrow(Week0to2.sol),6]
  Total_inf1_int <-  Week0to2.sol[nrow(Week0to2.sol),7]
  Total_Rec1_int <-  Week0to2.sol[nrow(Week0to2.sol),8]
  D1_Nat_int <-  Week0to2.sol[nrow(Week0to2.sol),9]
  D1_Inf_int <-  Week0to2.sol[nrow(Week0to2.sol),10]
  WF2Hog_Count1_int <-  Week0to2.sol[nrow(Week0to2.sol),11]
  
  S2_int <- (1-Hog_Vacc_Eff)*1000
  E2_int <- 0
  I2_int <- 0 #1
  Q2_int <- 0
  R2_int <- Hog_Vacc_Eff*1000
  Total_inf2_int <- 0
  Total_Rec2_int <- 0
  D2_Nat_int <- 0
  D2_Inf_int <- 0
  WF2Hog_Count2_int <- 0
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  Q3_int <- 0
  R3_int <- 0
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  Q4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week0to2.sol[nrow(Week0to2.sol),42]
  HE_int <-  Week0to2.sol[nrow(Week0to2.sol),43]
  HI_int <-  Week0to2.sol[nrow(Week0to2.sol),44]
  HR_int <-  Week0to2.sol[nrow(Week0to2.sol),45]
  H_DN_int <-  Week0to2.sol[nrow(Week0to2.sol),46]
  H_DI_int <-  Week0to2.sol[nrow(Week0to2.sol),47]
  Hog2WF_Count_int <-  Week0to2.sol[nrow(Week0to2.sol),48]
  
  #====================================================================================================
  # Model Run  Weeks 2 to 4
  #====================================================================================================
  Week2to4.t <- seq(14.1,Week2to4,0.1)          
  Week2to4.init <- c(S1_int, E1_int, I1_int, Q1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, Q2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, Q3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, Q4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week2to4.sol <- lsoda(Week2to4.init, Week2to4.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 4 TO 6 (Room 3 Fill)####
  Week4to6 <- 42
  #====================================================================================================
  # Model Initial Values  Weeks 4 to 6
  #===========================================================================
  S1_int <- Week2to4.sol[nrow(Week2to4.sol),2]
  E1_int <- Week2to4.sol[nrow(Week2to4.sol),3]
  I1_int <- Week2to4.sol[nrow(Week2to4.sol),4]
  Q1_int <- Week2to4.sol[nrow(Week2to4.sol),5]
  R1_int <- Week2to4.sol[nrow(Week2to4.sol),6]
  Total_inf1_int <- Week2to4.sol[nrow(Week2to4.sol),7]
  Total_Rec1_int <- Week2to4.sol[nrow(Week2to4.sol),8]
  D1_Nat_int <-  Week2to4.sol[nrow(Week2to4.sol),9]
  D1_Inf_int <-  Week2to4.sol[nrow(Week2to4.sol),10]
  WF2Hog_Count1_int <-  Week2to4.sol[nrow(Week2to4.sol),11]
  
  S2_int <- Week2to4.sol[nrow(Week2to4.sol),12]
  E2_int <- Week2to4.sol[nrow(Week2to4.sol),13]
  I2_int <- Week2to4.sol[nrow(Week2to4.sol),14]
  Q2_int <- Week2to4.sol[nrow(Week2to4.sol),15]
  R2_int <- Week2to4.sol[nrow(Week2to4.sol),16]
  Total_inf2_int <- Week2to4.sol[nrow(Week2to4.sol),17]
  Total_Rec2_int <- Week2to4.sol[nrow(Week2to4.sol),18]
  D2_Nat_int <- Week2to4.sol[nrow(Week2to4.sol),19]
  D2_Inf_int <- Week2to4.sol[nrow(Week2to4.sol),20]
  WF2Hog_Count2_int <- Week2to4.sol[nrow(Week2to4.sol),21]
  
  S3_int <- (1-Hog_Vacc_Eff)*1000
  E3_int <- 0
  I3_int <- 0
  Q3_int <- 0
  R3_int <-  Hog_Vacc_Eff*1000
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  Q4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week2to4.sol[nrow(Week2to4.sol),42]
  HE_int <-  Week2to4.sol[nrow(Week2to4.sol),43]
  HI_int <-  Week2to4.sol[nrow(Week2to4.sol),44]
  HR_int <-  Week2to4.sol[nrow(Week2to4.sol),45]
  H_DN_int <-  Week2to4.sol[nrow(Week2to4.sol),46]
  H_DI_int <-  Week2to4.sol[nrow(Week2to4.sol),47]
  Hog2WF_Count_int <-  Week2to4.sol[nrow(Week2to4.sol),48]
  
  #====================================================================================================
  # Model Run  Weeks 4 to 6
  #====================================================================================================
  Week4to6.t <- seq(28.1,Week4to6,0.1)         
  Week4to6.init <-c(S1_int, E1_int, I1_int, Q1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                    S2_int, E2_int, I2_int, Q2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                    S3_int, E3_int, I3_int, Q3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                    S4_int, E4_int, I4_int, Q4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                    HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week4to6.sol <- lsoda(Week4to6.init, Week4to6.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 6 TO 8 (Room 4 Fill)####
  Week6to8 <- 56
  #====================================================================================================
  # Model Initial Values  Weeks 6 to 8
  #===========================================================================
  S1_int <- Week4to6.sol[nrow(Week4to6.sol),2]
  E1_int <- Week4to6.sol[nrow(Week4to6.sol),3]
  I1_int <- Week4to6.sol[nrow(Week4to6.sol),4]
  Q1_int <- Week4to6.sol[nrow(Week4to6.sol),5]
  R1_int <- Week4to6.sol[nrow(Week4to6.sol),6]
  Total_inf1_int <- Week4to6.sol[nrow(Week4to6.sol),7]
  Total_Rec1_int <- Week4to6.sol[nrow(Week4to6.sol),8]
  D1_Nat_int <-  Week4to6.sol[nrow(Week4to6.sol),9]
  D1_Inf_int <-  Week4to6.sol[nrow(Week4to6.sol),10]
  WF2Hog_Count1_int <-  Week4to6.sol[nrow(Week4to6.sol),11]
  
  S2_int <- Week4to6.sol[nrow(Week4to6.sol),12]
  E2_int <- Week4to6.sol[nrow(Week4to6.sol),13]
  I2_int <- Week4to6.sol[nrow(Week4to6.sol),14]
  Q2_int <- Week4to6.sol[nrow(Week4to6.sol),15]
  R2_int <- Week4to6.sol[nrow(Week4to6.sol),16]
  Total_inf2_int <- Week4to6.sol[nrow(Week4to6.sol),17]
  Total_Rec2_int <- Week4to6.sol[nrow(Week4to6.sol),18]
  D2_Nat_int <- Week4to6.sol[nrow(Week4to6.sol),19]
  D2_Inf_int <- Week4to6.sol[nrow(Week4to6.sol),20]
  WF2Hog_Count2_int <- Week4to6.sol[nrow(Week4to6.sol),21]
  
  S3_int <- Week4to6.sol[nrow(Week4to6.sol),22]
  E3_int <- Week4to6.sol[nrow(Week4to6.sol),23]
  I3_int <- Week4to6.sol[nrow(Week4to6.sol),24]
  Q3_int <- Week4to6.sol[nrow(Week4to6.sol),25]
  R3_int <-  Week4to6.sol[nrow(Week4to6.sol),26]
  Total_inf3_int <- Week4to6.sol[nrow(Week4to6.sol),27]
  Total_Rec3_int <- Week4to6.sol[nrow(Week4to6.sol),28]
  D3_Nat_int <- Week4to6.sol[nrow(Week4to6.sol),29]
  D3_Inf_int <- Week4to6.sol[nrow(Week4to6.sol),30]
  WF2Hog_Count3_int <- Week4to6.sol[nrow(Week4to6.sol),31]
  
  S4_int <- (1-Hog_Vacc_Eff)*999
  E4_int <- 0
  I4_int <- 1
  Q4_int <- 0
  R4_int <- Hog_Vacc_Eff*999
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week4to6.sol[nrow(Week4to6.sol),42]
  HE_int <-  Week4to6.sol[nrow(Week4to6.sol),43]
  HI_int <-  Week4to6.sol[nrow(Week4to6.sol),44]
  HR_int <-  Week4to6.sol[nrow(Week4to6.sol),45]
  H_DN_int <-  Week4to6.sol[nrow(Week4to6.sol),46]
  H_DI_int <-  Week4to6.sol[nrow(Week4to6.sol),47]
  Hog2WF_Count_int <-  Week4to6.sol[nrow(Week4to6.sol),48]
  
  #====================================================================================================
  # Model Run  Weeks 6 to 8
  #====================================================================================================
  Week6to8.t <- seq(42.1,Week6to8,0.1)           
  Week6to8.init <- c(S1_int, E1_int, I1_int, Q1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, Q2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, Q3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, Q4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week6to8.sol <- lsoda(Week6to8.init, Week6to8.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 8 TO 23 (Room 1 Exit) ####
  Week8to23 <- 161
  #====================================================================================================
  # Model Initial Values  Weeks 8 to 23
  #===========================================================================
  S1_int <- Week6to8.sol[nrow(Week6to8.sol),2]
  E1_int <- Week6to8.sol[nrow(Week6to8.sol),3]
  I1_int <- Week6to8.sol[nrow(Week6to8.sol),4]
  Q1_int <- Week6to8.sol[nrow(Week6to8.sol),5]
  R1_int <- Week6to8.sol[nrow(Week6to8.sol),6]
  Total_inf1_int <- Week6to8.sol[nrow(Week6to8.sol),7]
  Total_Rec1_int <- Week6to8.sol[nrow(Week6to8.sol),8]
  D1_Nat_int <-  Week6to8.sol[nrow(Week6to8.sol),9]
  D1_Inf_int <-  Week6to8.sol[nrow(Week6to8.sol),10]
  WF2Hog_Count1_int <-  Week6to8.sol[nrow(Week6to8.sol),11]
  
  
  S2_int <- Week6to8.sol[nrow(Week6to8.sol),12]
  E2_int <- Week6to8.sol[nrow(Week6to8.sol),13]
  I2_int <- Week6to8.sol[nrow(Week6to8.sol),14]
  Q2_int <- Week6to8.sol[nrow(Week6to8.sol),15]
  R2_int <- Week6to8.sol[nrow(Week6to8.sol),16]
  Total_inf2_int <- Week6to8.sol[nrow(Week6to8.sol),17]
  Total_Rec2_int <- Week6to8.sol[nrow(Week6to8.sol),18]
  D2_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),19]
  D2_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),20]
  WF2Hog_Count2_int <- Week6to8.sol[nrow(Week6to8.sol),21]
  
  S3_int <- Week6to8.sol[nrow(Week6to8.sol),22]
  E3_int <- Week6to8.sol[nrow(Week6to8.sol),23]
  I3_int <- Week6to8.sol[nrow(Week6to8.sol),24]
  Q3_int <- Week6to8.sol[nrow(Week6to8.sol),25]
  R3_int <-  Week6to8.sol[nrow(Week6to8.sol),26]
  Total_inf3_int <- Week6to8.sol[nrow(Week6to8.sol),27]
  Total_Rec3_int <- Week6to8.sol[nrow(Week6to8.sol),28]
  D3_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),29]
  D3_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),30]
  WF2Hog_Count3_int <- Week6to8.sol[nrow(Week6to8.sol),31]
  
  S4_int <- Week6to8.sol[nrow(Week6to8.sol),32]
  E4_int <- Week6to8.sol[nrow(Week6to8.sol),33]
  I4_int <- Week6to8.sol[nrow(Week6to8.sol),34]
  Q4_int <- Week6to8.sol[nrow(Week6to8.sol),35]
  R4_int <- Week6to8.sol[nrow(Week6to8.sol),36]
  Total_inf4_int <- Week6to8.sol[nrow(Week6to8.sol),37]
  Total_Rec4_int <- Week6to8.sol[nrow(Week6to8.sol),38]
  D4_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),39]
  D4_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),40]
  WF2Hog_Count4_int <- Week6to8.sol[nrow(Week6to8.sol),41]
  
  HS_int <-  Week6to8.sol[nrow(Week6to8.sol),42]
  HE_int <-  Week6to8.sol[nrow(Week6to8.sol),43]
  HI_int <-  Week6to8.sol[nrow(Week6to8.sol),44]
  HR_int <-  Week6to8.sol[nrow(Week6to8.sol),45]
  H_DN_int <-  Week6to8.sol[nrow(Week6to8.sol),46]
  H_DI_int <-  Week6to8.sol[nrow(Week6to8.sol),47]
  Hog2WF_Count_int <-  Week6to8.sol[nrow(Week6to8.sol),48]
  
  #====================================================================================================
  # Model Run  Weeks 8 to 23
  #====================================================================================================
  Week8to23.t <- seq(56.1,Week8to23,0.1)          
  Week8to23.init <- c(S1_int, E1_int, I1_int, Q1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                      S2_int, E2_int, I2_int, Q2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                      S3_int, E3_int, I3_int, Q3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                      S4_int, E4_int, I4_int, Q4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                      HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week8to23.sol <- lsoda(Week8to23.init, Week8to23.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 23 TO 25 (Room 2 Exit) ####
  Week23to25 <- 175
  #====================================================================================================
  # Model Initial Values  Weeks 23 to 25
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  Q1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week8to23.sol[nrow(Week8to23.sol),7]
  Total_Rec1_int <- Week8to23.sol[nrow(Week8to23.sol),8]
  D1_Nat_int <-  Week8to23.sol[nrow(Week8to23.sol),9]
  D1_Inf_int <-  Week8to23.sol[nrow(Week8to23.sol),10]
  WF2Hog_Count1_int <-  Week8to23.sol[nrow(Week8to23.sol),11]
  
  
  S2_int <- Week8to23.sol[nrow(Week8to23.sol),12]
  E2_int <- Week8to23.sol[nrow(Week8to23.sol),13]
  I2_int <- Week8to23.sol[nrow(Week8to23.sol),14]
  Q2_int <- Week8to23.sol[nrow(Week8to23.sol),15]
  R2_int <- Week8to23.sol[nrow(Week8to23.sol),16]
  Total_inf2_int <- Week8to23.sol[nrow(Week8to23.sol),17]
  Total_Rec2_int <- Week8to23.sol[nrow(Week8to23.sol),18]
  D2_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),19]
  D2_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),20]
  WF2Hog_Count2_int <- Week8to23.sol[nrow(Week8to23.sol),21]
  
  S3_int <- Week8to23.sol[nrow(Week8to23.sol),22]
  E3_int <- Week8to23.sol[nrow(Week8to23.sol),23]
  I3_int <- Week8to23.sol[nrow(Week8to23.sol),24]
  Q3_int <- Week8to23.sol[nrow(Week8to23.sol),25]
  R3_int <-  Week8to23.sol[nrow(Week8to23.sol),26]
  Total_inf3_int <- Week8to23.sol[nrow(Week8to23.sol),27]
  Total_Rec3_int <- Week8to23.sol[nrow(Week8to23.sol),28]
  D3_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),29]
  D3_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),30]
  WF2Hog_Count3_int <- Week8to23.sol[nrow(Week8to23.sol),31]
  
  S4_int <- Week8to23.sol[nrow(Week8to23.sol),32]
  E4_int <- Week8to23.sol[nrow(Week8to23.sol),33]
  I4_int <- Week8to23.sol[nrow(Week8to23.sol),34]
  Q4_int <- Week8to23.sol[nrow(Week8to23.sol),35]
  R4_int <- Week8to23.sol[nrow(Week8to23.sol),36]
  Total_inf4_int <- Week8to23.sol[nrow(Week8to23.sol),37]
  Total_Rec4_int <- Week8to23.sol[nrow(Week8to23.sol),38]
  D4_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),39]
  D4_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),40]
  WF2Hog_Count4_int <- Week8to23.sol[nrow(Week8to23.sol),41]
  
  HS_int <-  Week8to23.sol[nrow(Week8to23.sol),42]
  HE_int <-  Week8to23.sol[nrow(Week8to23.sol),43]
  HI_int <-  Week8to23.sol[nrow(Week8to23.sol),44]
  HR_int <-  Week8to23.sol[nrow(Week8to23.sol),45]
  H_DN_int <-  Week8to23.sol[nrow(Week8to23.sol),46]
  H_DI_int <-  Week8to23.sol[nrow(Week8to23.sol),47]
  Hog2WF_Count_int <-  Week8to23.sol[nrow(Week8to23.sol),48]
  
  #====================================================================================================
  # Model Run  Weeks 23 to 25
  #====================================================================================================
  Week23to25.t <- seq(161.1,Week23to25,0.1)         
  Week23to25.init <- c(S1_int, E1_int, I1_int, Q1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, Q2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, Q3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, Q4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week23to25.sol <- lsoda(Week23to25.init, Week23to25.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 25 TO 27 (Room 3 Exit) ####
  Week25to27 <- 189
  #====================================================================================================
  # Model Initial Values  Weeks 25 to 27
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  Q1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week23to25.sol[nrow(Week23to25.sol),7]
  Total_Rec1_int <- Week23to25.sol[nrow(Week23to25.sol),8]
  D1_Nat_int <-  Week23to25.sol[nrow(Week23to25.sol),9]
  D1_Inf_int <-  Week23to25.sol[nrow(Week23to25.sol),10]
  WF2Hog_Count1_int <-  Week23to25.sol[nrow(Week23to25.sol),11]
  
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  Q2_int <- 0
  R2_int <- 0
  Total_inf2_int <- Week23to25.sol[nrow(Week23to25.sol),17]
  Total_Rec2_int <- Week23to25.sol[nrow(Week23to25.sol),18]
  D2_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),19]
  D2_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),20]
  WF2Hog_Count2_int <- Week23to25.sol[nrow(Week23to25.sol),21]
  
  S3_int <- Week23to25.sol[nrow(Week23to25.sol),22]
  E3_int <- Week23to25.sol[nrow(Week23to25.sol),23]
  I3_int <- Week23to25.sol[nrow(Week23to25.sol),24]
  Q3_int <- Week23to25.sol[nrow(Week23to25.sol),25]
  R3_int <-  Week23to25.sol[nrow(Week23to25.sol),26]
  Total_inf3_int <- Week23to25.sol[nrow(Week23to25.sol),27]
  Total_Rec3_int <- Week23to25.sol[nrow(Week23to25.sol),28]
  D3_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),29]
  D3_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),30]
  WF2Hog_Count3_int <- Week23to25.sol[nrow(Week23to25.sol),31]
  
  S4_int <- Week23to25.sol[nrow(Week23to25.sol),32]
  E4_int <- Week23to25.sol[nrow(Week23to25.sol),33]
  I4_int <- Week23to25.sol[nrow(Week23to25.sol),34]
  Q4_int <- Week23to25.sol[nrow(Week23to25.sol),35]
  R4_int <- Week23to25.sol[nrow(Week23to25.sol),36]
  Total_inf4_int <- Week23to25.sol[nrow(Week23to25.sol),37]
  Total_Rec4_int <- Week23to25.sol[nrow(Week23to25.sol),38]
  D4_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),39]
  D4_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),40]
  WF2Hog_Count4_int <- Week23to25.sol[nrow(Week23to25.sol),41]
  
  HS_int <-  Week23to25.sol[nrow(Week23to25.sol),42]
  HE_int <-  Week23to25.sol[nrow(Week23to25.sol),43]
  HI_int <-  Week23to25.sol[nrow(Week23to25.sol),44]
  HR_int <-  Week23to25.sol[nrow(Week23to25.sol),45]
  H_DN_int <-  Week23to25.sol[nrow(Week23to25.sol),46]
  H_DI_int <-  Week23to25.sol[nrow(Week23to25.sol),47]
  Hog2WF_Count_int <-  Week23to25.sol[nrow(Week23to25.sol),48]
  
  #====================================================================================================
  # Model Run  Weeks 25 to 27
  #====================================================================================================
  Week25to27.t <- seq(175.1,Week25to27,0.1)          
  Week25to27.init <- c(S1_int, E1_int, I1_int, Q1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, Q2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, Q3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, Q4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week25to27.sol <- lsoda(Week25to27.init, Week25to27.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 27 TO 29 (Room 4 Exit) ####
  Week27to29 <- 203
  #====================================================================================================
  # Model Initial Values  Weeks 27to 29
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  Q1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week25to27.sol[nrow(Week25to27.sol),7]
  Total_Rec1_int <- Week25to27.sol[nrow(Week25to27.sol),8]
  D1_Nat_int <-  Week25to27.sol[nrow(Week25to27.sol),9]
  D1_Inf_int <-  Week25to27.sol[nrow(Week25to27.sol),10]
  WF2Hog_Count1_int <-  Week25to27.sol[nrow(Week25to27.sol),11]
  
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  Q2_int <- 0
  R2_int <- 0
  Total_inf2_int <- Week25to27.sol[nrow(Week25to27.sol),17]
  Total_Rec2_int <- Week25to27.sol[nrow(Week25to27.sol),18]
  D2_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),19]
  D2_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),20]
  WF2Hog_Count2_int <- Week25to27.sol[nrow(Week25to27.sol),21]
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  Q3_int <- 0
  R3_int <- 0
  Total_inf3_int <- Week25to27.sol[nrow(Week25to27.sol),27]
  Total_Rec3_int <- Week25to27.sol[nrow(Week25to27.sol),28]
  D3_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),29]
  D3_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),30]
  WF2Hog_Count3_int <- Week25to27.sol[nrow(Week25to27.sol),31]
  
  S4_int <- Week25to27.sol[nrow(Week25to27.sol),32]
  E4_int <- Week25to27.sol[nrow(Week25to27.sol),33]
  I4_int <- Week25to27.sol[nrow(Week25to27.sol),34]
  Q4_int <- Week25to27.sol[nrow(Week25to27.sol),35]
  R4_int <- Week25to27.sol[nrow(Week25to27.sol),36]
  Total_inf4_int <- Week25to27.sol[nrow(Week25to27.sol),37]
  Total_Rec4_int <- Week25to27.sol[nrow(Week25to27.sol),38]
  D4_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),39]
  D4_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),40]
  WF2Hog_Count4_int <- Week25to27.sol[nrow(Week25to27.sol),41]
  
  HS_int <-  Week25to27.sol[nrow(Week25to27.sol),42]
  HE_int <-  Week25to27.sol[nrow(Week25to27.sol),43]
  HI_int <-  Week25to27.sol[nrow(Week25to27.sol),44]
  HR_int <-  Week25to27.sol[nrow(Week25to27.sol),45]
  H_DN_int <-  Week25to27.sol[nrow(Week25to27.sol),46]
  H_DI_int <-  Week25to27.sol[nrow(Week25to27.sol),47]
  Hog2WF_Count_int <-  Week25to27.sol[nrow(Week25to27.sol),48]
  
  #====================================================================================================
  # Model Run  Weeks 27 to 29
  #====================================================================================================
  Week27to29.t <- seq(189.1,Week27to29,0.1)            
  Week27to29.init <- c(S1_int, E1_int, I1_int, Q1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, Q2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, Q3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, Q4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week27to29.sol <- lsoda(Week27to29.init, Week27to29.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### COMBINING THE DATASETS ####
  Week0to2.sol.df <- as.data.frame(Week0to2.sol)
  colnames(Week0to2.sol.df) <- c("time",
                                 "S1", "E1", "I1", "Q1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "Q2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "Q3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "Q4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week2to4.sol.df <- as.data.frame(Week2to4.sol)
  colnames(Week2to4.sol.df) <- c("time",
                                 "S1", "E1", "I1", "Q1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "Q2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "Q3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "Q4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week4to6.sol.df <- as.data.frame(Week4to6.sol)
  colnames(Week4to6.sol.df) <- c("time",
                                 "S1", "E1", "I1", "Q1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "Q2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "Q3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "Q4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week6to8.sol.df <- as.data.frame(Week6to8.sol)
  colnames(Week6to8.sol.df) <- c("time",
                                 "S1", "E1", "I1", "Q1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "Q2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "Q3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "Q4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week8to23.sol.df <- as.data.frame(Week8to23.sol)
  colnames(Week8to23.sol.df) <- c("time",
                                  "S1", "E1", "I1", "Q1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                  "S2", "E2", "I2", "Q2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                  "S3", "E3", "I3", "Q3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                  "S4", "E4", "I4", "Q4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                  "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week23to25.sol.df <- as.data.frame(Week23to25.sol)
  colnames(Week23to25.sol.df) <- c("time",
                                   "S1", "E1", "I1", "Q1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "Q2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "Q3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "Q4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week25to27.sol.df <- as.data.frame(Week25to27.sol)
  colnames(Week25to27.sol.df) <- c("time",
                                   "S1", "E1", "I1", "Q1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "Q2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "Q3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "Q4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week27to29.sol.df <- as.data.frame(Week27to29.sol)
  colnames(Week27to29.sol.df) <- c("time",
                                   "S1", "E1", "I1", "Q1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "Q2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "Q3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "Q4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  
  fulldata.df <- rbind(Week0to2.sol.df, Week2to4.sol.df, Week4to6.sol.df, Week6to8.sol.df, 
                       Week8to23.sol.df, Week23to25.sol.df, Week25to27.sol.df, Week27to29.sol.df)
  
  ## Is there a WF infection
  HI_inf <- fulldata.df %>%
    summarise(total_I = sum(HI)) %>%
    mutate(HI_Inf = ifelse(total_I>1, 1, 0)) %>%
    select(HI_Inf) %>% as.numeric()
  #Time to 1st WF infection
  Time2WFinf0 <- fulldata.df %>%
    arrange(time) %>%
    filter(HI > 0) %>%
    filter(row_number() == 1) %>%
    mutate(TimeSineInfIntro = time - 42) %>%
    select(TimeSineInfIntro) %>% as.numeric()
  
  Time2WFinf.5 <- fulldata.df %>%
    arrange(time) %>%
    filter(HI > 0.5) %>%
    filter(row_number() == 1) %>%
    mutate(TimeSineInfIntro = time - 42) %>%
    select(TimeSineInfIntro) %>% as.numeric()
  ## Time to Second Room Infection
  Timeto2room <- fulldata.df %>% 
    summarise(t_cross = time[which(E1 > 1 | 
                                     E2 > 1 |
                                     E3 > 1)]) %>% 
    summarise(min_t_cross = min(t_cross)) %>% as.numeric()
  ## Total Hogs infected
  TotHogsInf <- fulldata.df %>% 
    #select(t, iter, Total_inf1, Total_inf2, Total_inf3, Total_inf4) %>% 
    mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
           Tot_rec = Total_Rec1 + Total_Rec2 + Total_Rec3 + Total_Rec4) %>% 
    summarise(Tot_inf2 = max(Tot_inf),
              Tot_rec2 = max(Tot_rec))
  ## Minimum Time to Peak Infection
  Time2PeakInf <- fulldata.df %>% 
    mutate(S_tot = S1+S2+S3+S4,                             
           E_tot = E1+E2+E3+E4,
           I_tot = I1+I2+I3+I4,
           R_tot = R1+R2+R3+R4) %>%                                     
    summarise(max_I = max(I_tot),
              t_peak = time[which(I_tot == max(I_tot))]) %>% 
    summarise(min_t_peak = min(t_peak)) %>% as.numeric()
  ##Ro
  Ro <- fulldata.df %>% 
    arrange(time) %>% 
    mutate(Tot_rec = Total_Rec1 + Total_Rec2 + Total_Rec3 + Total_Rec4) %>%  
    filter(row_number() == n()) %>% 
    select(Tot_rec) %>% 
    mutate(Ro = ((log((4000 - Tot_rec)/4000) - log((4000 - (Hog_Vacc_Eff*4000))/4000)) /
                   ((((4000 - Tot_rec) / 4000) - (4000 - (Hog_Vacc_Eff*4000)) / 4000)))) %>% 
    select(Ro) %>% as.numeric()
  
  #Create Summary data frame
  fulldata.df.sum <- data.frame(iter = as.factor(i),
                                HI_inf = HI_inf,
                                Time2WFinf = Time2WFinf,
                                Timeto2room = Timeto2room,
                                TotHogsInf = TotHogsInf$Tot_inf2,
                                TotHogsRec = TotHogsInf$Tot_rec2,
                                Time2PeakInf = Time2PeakInf,
                                Ro = Ro,
                                beta_Hog = beta_Hog,
                                beta_Hog_ind = beta_Hog_ind,
                                beta_S2WF = beta_S2WF,
                                beta_WF2S = beta_WF2S,
                                sigma_Hog = sigma_Hog,
                                delta_Hog = delta_Hog,
                                beta_WF2WF = beta_WF2WF,
                                sigma_WF = sigma_WF,
                                delta_WF = delta_WF)
  
  results <- rbind(results, fulldata.df.sum)
  
  rm(list=ls()[! ls() %in% c("results", "n.itr")])  #Deletes everything from environment expect growing results dataframe
}


#### WF Flow Only (Room 4 Intro) ####
results <- NULL
n.itr <- 5000
# seed.num <- 616
seed.nums <- c(100:(100+n.itr))
for (i in c(1:n.itr)) {
  print(i)
  
  #### MODEL CONSTANTS ####
  #Hog constants
  beta_Hog <- rtri(n = 1, min = 0.001, max = 0.1, mode = (10/(1000*5)))
  beta_Hog_ind <- beta_Hog / 178 #(following Etbaigha et al 2018 paper)
  # beta_Hog_ind <- beta_Hog / 500
  Hog_Vacc_Eff <- 0
  u_Hog <- 0.00028     #natural death rate (following Etbaigha et al 2018 paper)
  u_inf_Hog <- 0 #0.1     #infected death rate (0 for now to problem solve)
  w_Hog <- 0 #1/180         #(following Etbaigha et al 2018 paper)  White et al shows a range from 56 - 112 days
  sigma_Hog <- abs(1 / rnorm(n = 1, mean = 2, sd = 1))   #Latent period h1n1 and h3n2
  delta_Hog <- abs(1 / rnorm(n = 1, mean = 5, sd = 1))   #Infectious period h1n1 and h3n2
  
  #Interspecies transmission
  ##applied the beta = Ro / N*D from the modeling book
  # Ro from the following website for swine to human infection 
  # https://bmcmedicine.biomedcentral.com/articles/10.1186/1741-7015-7-30#Sec6   (from super-strain figure) possible the high end of transmisilibty for this parameter
  # beta_S2WF <- (2.3 / (4002*5))  #  #Hogs + #Workforce * Duration of Hogs
  #Alt beta_S2WF from fair outbreak
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3557885/
  #assuming one hog got the first three confirmed cases sick and infectious period of 5 days Ro = 3/5 = 0.6
  beta_S2WF <- runif(1, min = 0.000001, max = 0.0001)
  #when Ro WF2S = Ro S2WF
  # beta_WF2S <- (2.3 / (4002*3))  #  #Hogs + #Workforce * Duration of Hogs
  # WF2S calculated from Swine Outbreak of Pandemic Influenza A Virus on a Canadian Research Farm Supports Human-to-Swine Transmission paper
  beta_WF2S <- runif(1, min = 0.000001, max = 0.0001)    #Hogs + #Workforce * Duration of WF
  #Human constants
  #applied the beta = Ro / N*D from the modeling book
  beta_WF2WF <-   rtri(n = 1, min = 0.0000001, max = 0.64, mode = 0.32) # http://m-hikari.com/ams/ams-2013/ams-41-44-2013/joseAMS41-44-2013.pdf 
  Hum_Vacc_Eff <- 0.0 
  u_WF <- 0     #Natural death rate of human (0 for now to problem solve)
  u_inf_WF <- 0   #Infected death rate of human (0 for now to problem solve)
  w_WF <- 0      #Recovery rate for human (0 for now to problem solve)
  sigma_WF <- abs(1 / rnorm(n = 1, mean = 2, sd = 1))   #Latent period for humans
  delta_WF <- abs(1 / rnorm(n = 1, mean = 3, sd = 1))   #Infectious Period for humans
  
  #### MODEL SETUP ####
  # Model Parameters
  TransmissionModel.par <- c(beta_Hog = beta_Hog, beta_Hog_ind = beta_Hog_ind, beta_S2WF = beta_S2WF, beta_WF2S = beta_WF2S,
                             w_Hog = w_Hog, sigma_Hog = sigma_Hog, delta_Hog = delta_Hog, u_Hog = u_Hog, u_inf_Hog = u_inf_Hog,
                             beta_WF2WF = beta_WF2WF,
                             w_WF = w_WF, sigma_WF = sigma_WF, delta_WF = delta_WF, u_WF = u_WF, u_inf_WF = u_inf_WF)
  #=====================================================================================================
  #Model Equations 
  #=====================================================================================================
  TransmissionModel.dyn <- function(t,var,par) {
    # Rename the parameters 
    beta_Hog <- par[1]
    beta_Hog_ind <- par[2]
    beta_S2WF <- par[3]
    beta_WF2S <- par[4]
    w_Hog <- par[5]
    sigma_Hog <- par[6]
    delta_Hog <- par[7]
    u_Hog <- par[8]
    u_inf_Hog<- par[9]
    beta_WF2WF <- par[10]
    w_WF <- par[11]
    sigma_WF <- par[12]
    delta_WF <- par[13]
    u_WF <- par[14]
    u_inf_WF <- par[15]
    
    S1 <- var[1]
    E1 <- var[2]
    I1 <- var[3]
    R1 <- var[4]
    Total_inf1 <- var[5] 
    True_Rec1 <- var[6]
    D1_Nat <- var[7]
    D1_Inf <- var[8]
    WF2Hog_Count1 <- var[9] 
    
    S2 <- var[10]
    E2 <- var[11]
    I2 <- var[12]
    R2 <- var[13]
    Total_inf2 <- var[14]
    True_Rec2 <- var[15]
    D2_Nat <- var[16]
    D2_Inf <- var[17]
    WF2Hog_Count2 <- var[18]
    
    S3 <- var[19]
    E3 <- var[20]
    I3 <- var[21]
    R3 <- var[22]
    Total_inf3 <- var[23]
    True_Rec3 <- var[24]
    D3_Nat <- var[25]
    D3_Inf <- var[26]
    WF2Hog_Count3 <- var[27]   
    
    S4 <- var[28]
    E4 <- var[29]
    I4 <- var[30]
    R4 <- var[31]
    Total_inf4 <- var[32]
    True_Rec4 <- var[33]
    D4_Nat <- var[34]
    D4_Inf <- var[35]
    WF2Hog_Count4 <- var[36]
    
    HS <- var[37]
    HE <- var[38]
    HI <- var[39] 
    HR <- var[40]
    H_DN <- var[41]
    H_DI <- var[42]
    Hog2WF_Count <- var[43] 
    
    # Calculate the derivatives
    #Room 1
    dS1 <- (w_Hog*R1) - (beta_Hog*I1*S1) - (beta_Hog_ind*(I2+I3+I4)*S1) - (beta_WF2S*HI*S1) - (u_Hog*S1)
    dE1 <-(beta_Hog*I1*S1) + (beta_Hog_ind*(I2+I3+I4)*S1) + (beta_WF2S*HI*S1) - sigma_Hog*E1 - (u_Hog*E1)
    dI1 <- (sigma_Hog*E1) - (delta_Hog*I1) - (u_inf_Hog*I1)
    dR1 <- (delta_Hog*I1) - (w_Hog*R1) - (u_Hog*R1)
    dTotal_inf1 <- (sigma_Hog*E1)
    dTotal_Rec1 <- (delta_Hog*I1)
    dD1_Nat <- (u_Hog*S1) + (u_Hog*E1) + (u_Hog*R1)
    dD1_Inf <- (u_inf_Hog*I1)
    dWF2Hog_Count1 <- (beta_WF2S*HI*S1)
    
    #Room 2
    dS2 <- (w_Hog*R2) - (beta_Hog*I2*S2) - (beta_Hog_ind*(I3+I4)*S2) - (beta_WF2S*HI*S2) - (u_Hog*S2)
    dE2 <-(beta_Hog*I2*S2) + (beta_Hog_ind*(I3+I4)*S2) + (beta_WF2S*HI*S2) - (sigma_Hog*E2) - (u_Hog*E2) 
    dI2 <- (sigma_Hog*E2) - (delta_Hog*I2) - (u_inf_Hog*I2)
    dR2 <- (delta_Hog*I2) - (w_Hog*R2) - (u_Hog*R2)
    dTotal_inf2 <- (sigma_Hog*E2)
    dTotal_Rec2 <- (delta_Hog*I2)
    dD2_Nat <- (u_Hog*S2) + (u_Hog*E2) + (u_Hog*R2)
    dD2_Inf <- (u_inf_Hog*I2)
    dWF2Hog_Count2 <- (beta_WF2S*HI*S2)
    
    #Room 3
    dS3 <- (w_Hog*R3) - (beta_Hog*I3*S3) - (beta_Hog_ind*(I4)*S3) - (beta_WF2S*HI*S3) - (u_Hog*S3)
    dE3 <-(beta_Hog*I3*S3) + (beta_Hog_ind*(I4)*S3) + (beta_WF2S*HI*S3) - (sigma_Hog*E3) - (u_Hog*E3)
    dI3 <- (sigma_Hog*E3) - (delta_Hog*I3) - (u_inf_Hog*I3)
    dR3 <- (delta_Hog*I3) - (w_Hog*R3) - (u_Hog*R3)
    dTotal_inf3 <- (sigma_Hog*E3)
    dTotal_Rec3 <- (delta_Hog*I3)
    dD3_Nat <- (u_Hog*S3) + (u_Hog*E3) + (u_Hog*R3)
    dD3_Inf <- (u_inf_Hog*I3)
    dWF2Hog_Count3 <- (beta_WF2S*HI*S3)
    
    #Room 4
    dS4 <- w_Hog*R4 - (beta_Hog*I4*S4) - (beta_Hog_ind*(0)*S4) - (beta_WF2S*HI*S4) - (u_Hog*S4)
    dE4 <-(beta_Hog*I4*S4) + (beta_Hog_ind*(0)*S4) + (beta_WF2S*HI*S4) - (sigma_Hog*E4) - (u_Hog*E4) 
    dI4 <- (sigma_Hog*E4) - (delta_Hog*I4) - (u_inf_Hog*I4)
    dR4 <- (delta_Hog*I4) - (w_Hog*R4) - (u_Hog*R4)
    dTotal_inf4 <- (sigma_Hog*E4)
    dTotal_Rec4 <- (delta_Hog*I4)
    dD4_Nat <- (u_Hog*S4) + (u_Hog*E4) + (u_Hog*R4)
    dD4_Inf <- (u_inf_Hog*I4)
    dWF2Hog_Count4 <- (beta_WF2S*HI*S4)
    
    #Humans
    dHS <- (w_WF*HR)- (beta_WF2WF*HI*HS) - (beta_S2WF*(I1+I2+I3+I4)*HS) - (u_WF*HS)
    dHE <- (beta_WF2WF*HI*HS) + (beta_S2WF*(I1+I2+I3+I4)*HS) - (sigma_WF*HE) - (u_WF*HE)
    dHI <- (sigma_WF*HE) - (delta_WF*HI) - (u_inf_WF*HI)
    dHR <- (delta_WF*HI) - (w_WF*HR) - u_WF*HR
    dH_DN <- (u_WF*HS) + (u_WF*HE) + (u_WF*HR)
    dH_DI <- (u_inf_WF*HI)
    dHog2WF_Count <- (beta_S2WF*(I1+I2+I3+I4)*HS)
    
    
    # Last instruction: return a list 
    
    return(list(c(dS1, dE1, dI1, dR1, dTotal_inf1, dTotal_Rec1, dD1_Nat, dD1_Inf, dWF2Hog_Count1,
                  dS2, dE2, dI2, dR2, dTotal_inf2, dTotal_Rec2, dD2_Nat, dD2_Inf, dWF2Hog_Count2,
                  dS3, dE3, dI3, dR3, dTotal_inf3, dTotal_Rec3, dD3_Nat, dD3_Inf, dWF2Hog_Count3,
                  dS4, dE4, dI4, dR4, dTotal_inf4, dTotal_Rec4, dD4_Nat, dD4_Inf, dWF2Hog_Count4,
                  dHS, dHE, dHI, dHR, dH_DN, dH_DI, dHog2WF_Count)))
    
  }
  
  #### MODEL ITERATIONS ####
  #### WEEK 0 TO 2 (Room 1 fill) ####
  Week0to2 <- 14
  #====================================================================================================
  # Model Initial Values  Weeks 0 to 2
  #===========================================================================
  S1_int <- (1-Hog_Vacc_Eff)*1000
  E1_int <- 0
  I1_int <- 0 #1
  R1_int <- Hog_Vacc_Eff*1000
  Total_inf1_int <- 0
  Total_Rec1_int <- 0
  D1_Nat_int <- 0
  D1_Inf_int <- 0
  WF2Hog_Count1_int <- 0
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- 0
  Total_Rec2_int <- 0
  D2_Nat_int <- 0
  D2_Inf_int <- 0
  WF2Hog_Count2_int <- 0
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <- 2
  HE_int <- 0
  HI_int <- 0
  HR_int <- 0
  H_DN_int <- 0
  H_DI_int <- 0
  Hog2WF_Count_int <- 0
  
  #====================================================================================================
  # Model Run  Weeks 0 to 2
  #====================================================================================================
  Week0to2.t <- seq(0,Week0to2,0.1)   
  Week0to2.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week0to2.sol <- lsoda(Week0to2.init, Week0to2.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 2 TO 4 (Room 2 Fill)####
  Week2to4 <- 28
  #====================================================================================================
  # Model Initial Values  Weeks 2 to 4
  #===========================================================================
  S1_int <- Week0to2.sol[nrow(Week0to2.sol),2]
  E1_int <- Week0to2.sol[nrow(Week0to2.sol),3]
  I1_int <- Week0to2.sol[nrow(Week0to2.sol),4]
  R1_int <- Week0to2.sol[nrow(Week0to2.sol),5]
  Total_inf1_int <-  Week0to2.sol[nrow(Week0to2.sol),6]
  Total_Rec1_int <-  Week0to2.sol[nrow(Week0to2.sol),7]
  D1_Nat_int <-  Week0to2.sol[nrow(Week0to2.sol),8]
  D1_Inf_int <-  Week0to2.sol[nrow(Week0to2.sol),9]
  WF2Hog_Count1_int <-  Week0to2.sol[nrow(Week0to2.sol),10]
  
  S2_int <- (1-Hog_Vacc_Eff)*1000
  E2_int <- 0
  I2_int <- 0 #1
  R2_int <- Hog_Vacc_Eff*1000
  Total_inf2_int <- 0
  Total_Rec2_int <- 0
  D2_Nat_int <- 0
  D2_Inf_int <- 0
  WF2Hog_Count2_int <- 0
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week0to2.sol[nrow(Week0to2.sol),38]
  HE_int <-  Week0to2.sol[nrow(Week0to2.sol),39]
  HI_int <-  Week0to2.sol[nrow(Week0to2.sol),40]
  HR_int <-  Week0to2.sol[nrow(Week0to2.sol),41]
  H_DN_int <-  Week0to2.sol[nrow(Week0to2.sol),42]
  H_DI_int <-  Week0to2.sol[nrow(Week0to2.sol),43]
  Hog2WF_Count_int <-  Week0to2.sol[nrow(Week0to2.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 2 to 4
  #====================================================================================================
  Week2to4.t <- seq(14.1,Week2to4,0.1)          
  Week2to4.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week2to4.sol <- lsoda(Week2to4.init, Week2to4.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 4 TO 6 (Room 3 Fill)####
  Week4to6 <- 42
  #====================================================================================================
  # Model Initial Values  Weeks 4 to 6
  #===========================================================================
  S1_int <- Week2to4.sol[nrow(Week2to4.sol),2]
  E1_int <- Week2to4.sol[nrow(Week2to4.sol),3]
  I1_int <- Week2to4.sol[nrow(Week2to4.sol),4]
  R1_int <- Week2to4.sol[nrow(Week2to4.sol),5]
  Total_inf1_int <- Week2to4.sol[nrow(Week2to4.sol),6]
  Total_Rec1_int <- Week2to4.sol[nrow(Week2to4.sol),7]
  D1_Nat_int <-  Week2to4.sol[nrow(Week2to4.sol),8]
  D1_Inf_int <-  Week2to4.sol[nrow(Week2to4.sol),9]
  WF2Hog_Count1_int <-  Week2to4.sol[nrow(Week2to4.sol),10]
  
  S2_int <- Week2to4.sol[nrow(Week2to4.sol),11]
  E2_int <- Week2to4.sol[nrow(Week2to4.sol),12]
  I2_int <- Week2to4.sol[nrow(Week2to4.sol),13]
  R2_int <- Week2to4.sol[nrow(Week2to4.sol),14]
  Total_inf2_int <- Week2to4.sol[nrow(Week2to4.sol),15]
  Total_Rec2_int <- Week2to4.sol[nrow(Week2to4.sol),16]
  D2_Nat_int <- Week2to4.sol[nrow(Week2to4.sol),17]
  D2_Inf_int <- Week2to4.sol[nrow(Week2to4.sol),18]
  WF2Hog_Count2_int <- Week2to4.sol[nrow(Week2to4.sol),19]
  
  S3_int <- (1-Hog_Vacc_Eff)*1000
  E3_int <- 0
  I3_int <- 0
  R3_int <-  Hog_Vacc_Eff*1000
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week2to4.sol[nrow(Week2to4.sol),38]
  HE_int <-  Week2to4.sol[nrow(Week2to4.sol),39]
  HI_int <-  Week2to4.sol[nrow(Week2to4.sol),40]
  HR_int <-  Week2to4.sol[nrow(Week2to4.sol),41]
  H_DN_int <-  Week2to4.sol[nrow(Week2to4.sol),42]
  H_DI_int <-  Week2to4.sol[nrow(Week2to4.sol),43]
  Hog2WF_Count_int <-  Week2to4.sol[nrow(Week2to4.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 4 to 6
  #====================================================================================================
  Week4to6.t <- seq(28.1,Week4to6,0.1)         
  Week4to6.init <-c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                    S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                    S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                    S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                    HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week4to6.sol <- lsoda(Week4to6.init, Week4to6.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 6 TO 8 (Room 4 Fill)####
  Week6to8 <- 56
  #====================================================================================================
  # Model Initial Values  Weeks 6 to 8
  #===========================================================================
  S1_int <- Week4to6.sol[nrow(Week4to6.sol),2]
  E1_int <- Week4to6.sol[nrow(Week4to6.sol),3]
  I1_int <- Week4to6.sol[nrow(Week4to6.sol),4]
  R1_int <- Week4to6.sol[nrow(Week4to6.sol),5]
  Total_inf1_int <- Week4to6.sol[nrow(Week4to6.sol),6]
  Total_Rec1_int <- Week4to6.sol[nrow(Week4to6.sol),7]
  D1_Nat_int <-  Week4to6.sol[nrow(Week4to6.sol),8]
  D1_Inf_int <-  Week4to6.sol[nrow(Week4to6.sol),9]
  WF2Hog_Count1_int <-  Week4to6.sol[nrow(Week4to6.sol),10]
  
  S2_int <- Week4to6.sol[nrow(Week4to6.sol),11]
  E2_int <- Week4to6.sol[nrow(Week4to6.sol),12]
  I2_int <- Week4to6.sol[nrow(Week4to6.sol),13]
  R2_int <- Week4to6.sol[nrow(Week4to6.sol),14]
  Total_inf2_int <- Week4to6.sol[nrow(Week4to6.sol),15]
  Total_Rec2_int <- Week4to6.sol[nrow(Week4to6.sol),16]
  D2_Nat_int <- Week4to6.sol[nrow(Week4to6.sol),17]
  D2_Inf_int <- Week4to6.sol[nrow(Week4to6.sol),18]
  WF2Hog_Count2_int <- Week4to6.sol[nrow(Week4to6.sol),19]
  
  S3_int <- Week4to6.sol[nrow(Week4to6.sol),20]
  E3_int <- Week4to6.sol[nrow(Week4to6.sol),21]
  I3_int <- Week4to6.sol[nrow(Week4to6.sol),22]
  R3_int <-  Week4to6.sol[nrow(Week4to6.sol),23]
  Total_inf3_int <- Week4to6.sol[nrow(Week4to6.sol),24]
  Total_Rec3_int <- Week4to6.sol[nrow(Week4to6.sol),25]
  D3_Nat_int <- Week4to6.sol[nrow(Week4to6.sol),26]
  D3_Inf_int <- Week4to6.sol[nrow(Week4to6.sol),27]
  WF2Hog_Count3_int <- Week4to6.sol[nrow(Week4to6.sol),28]
  
  S4_int <- (1-Hog_Vacc_Eff)*999
  E4_int <- 0
  I4_int <- 1
  R4_int <- Hog_Vacc_Eff*999
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week4to6.sol[nrow(Week4to6.sol),38]
  HE_int <-  Week4to6.sol[nrow(Week4to6.sol),39]
  HI_int <-  Week4to6.sol[nrow(Week4to6.sol),40]
  HR_int <-  Week4to6.sol[nrow(Week4to6.sol),41]
  H_DN_int <-  Week4to6.sol[nrow(Week4to6.sol),42]
  H_DI_int <-  Week4to6.sol[nrow(Week4to6.sol),43]
  Hog2WF_Count_int <-  Week4to6.sol[nrow(Week4to6.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 6 to 8
  #====================================================================================================
  Week6to8.t <- seq(42.1,Week6to8,0.1)           
  Week6to8.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week6to8.sol <- lsoda(Week6to8.init, Week6to8.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 8 TO 23 (Room 1 Exit) ####
  Week8to23 <- 161
  #====================================================================================================
  # Model Initial Values  Weeks 8 to 23
  #===========================================================================
  S1_int <- Week6to8.sol[nrow(Week6to8.sol),2]
  E1_int <- Week6to8.sol[nrow(Week6to8.sol),3]
  I1_int <- Week6to8.sol[nrow(Week6to8.sol),4]
  R1_int <- Week6to8.sol[nrow(Week6to8.sol),5]
  Total_inf1_int <- Week6to8.sol[nrow(Week6to8.sol),6]
  Total_Rec1_int <- Week6to8.sol[nrow(Week6to8.sol),7]
  D1_Nat_int <-  Week6to8.sol[nrow(Week6to8.sol),8]
  D1_Inf_int <-  Week6to8.sol[nrow(Week6to8.sol),9]
  WF2Hog_Count1_int <-  Week6to8.sol[nrow(Week6to8.sol),10]
  
  
  S2_int <- Week6to8.sol[nrow(Week6to8.sol),11]
  E2_int <- Week6to8.sol[nrow(Week6to8.sol),12]
  I2_int <- Week6to8.sol[nrow(Week6to8.sol),13]
  R2_int <- Week6to8.sol[nrow(Week6to8.sol),14]
  Total_inf2_int <- Week6to8.sol[nrow(Week6to8.sol),15]
  Total_Rec2_int <- Week6to8.sol[nrow(Week6to8.sol),16]
  D2_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),17]
  D2_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),18]
  WF2Hog_Count2_int <- Week6to8.sol[nrow(Week6to8.sol),19]
  
  S3_int <- Week6to8.sol[nrow(Week6to8.sol),20]
  E3_int <- Week6to8.sol[nrow(Week6to8.sol),21]
  I3_int <- Week6to8.sol[nrow(Week6to8.sol),22]
  R3_int <-  Week6to8.sol[nrow(Week6to8.sol),23]
  Total_inf3_int <- Week6to8.sol[nrow(Week6to8.sol),24]
  Total_Rec3_int <- Week6to8.sol[nrow(Week6to8.sol),25]
  D3_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),26]
  D3_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),27]
  WF2Hog_Count3_int <- Week6to8.sol[nrow(Week6to8.sol),28]
  
  S4_int <- Week6to8.sol[nrow(Week6to8.sol),29]
  E4_int <- Week6to8.sol[nrow(Week6to8.sol),30]
  I4_int <- Week6to8.sol[nrow(Week6to8.sol),31]
  R4_int <- Week6to8.sol[nrow(Week6to8.sol),32]
  Total_inf4_int <- Week6to8.sol[nrow(Week6to8.sol),33]
  Total_Rec4_int <- Week6to8.sol[nrow(Week6to8.sol),34]
  D4_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),35]
  D4_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),36]
  WF2Hog_Count4_int <- Week6to8.sol[nrow(Week6to8.sol),37]
  
  HS_int <-  Week6to8.sol[nrow(Week6to8.sol),38]
  HE_int <-  Week6to8.sol[nrow(Week6to8.sol),39]
  HI_int <-  Week6to8.sol[nrow(Week6to8.sol),40]
  HR_int <-  Week6to8.sol[nrow(Week6to8.sol),41]
  H_DN_int <-  Week6to8.sol[nrow(Week6to8.sol),42]
  H_DI_int <-  Week6to8.sol[nrow(Week6to8.sol),43]
  Hog2WF_Count_int <-  Week6to8.sol[nrow(Week6to8.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 8 to 23
  #====================================================================================================
  Week8to23.t <- seq(56.1,Week8to23,0.1)          
  Week8to23.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                      S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                      S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                      S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                      HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week8to23.sol <- lsoda(Week8to23.init, Week8to23.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 23 TO 25 (Room 2 Exit) ####
  Week23to25 <- 175
  #====================================================================================================
  # Model Initial Values  Weeks 23 to 25
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week8to23.sol[nrow(Week8to23.sol),6]
  Total_Rec1_int <- Week8to23.sol[nrow(Week8to23.sol),7]
  D1_Nat_int <-  Week8to23.sol[nrow(Week8to23.sol),8]
  D1_Inf_int <-  Week8to23.sol[nrow(Week8to23.sol),9]
  WF2Hog_Count1_int <-  Week8to23.sol[nrow(Week8to23.sol),10]
  
  
  S2_int <- Week8to23.sol[nrow(Week8to23.sol),11]
  E2_int <- Week8to23.sol[nrow(Week8to23.sol),12]
  I2_int <- Week8to23.sol[nrow(Week8to23.sol),13]
  R2_int <- Week8to23.sol[nrow(Week8to23.sol),14]
  Total_inf2_int <- Week8to23.sol[nrow(Week8to23.sol),15]
  Total_Rec2_int <- Week8to23.sol[nrow(Week8to23.sol),16]
  D2_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),17]
  D2_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),18]
  WF2Hog_Count2_int <- Week8to23.sol[nrow(Week8to23.sol),19]
  
  S3_int <- Week8to23.sol[nrow(Week8to23.sol),20]
  E3_int <- Week8to23.sol[nrow(Week8to23.sol),21]
  I3_int <- Week8to23.sol[nrow(Week8to23.sol),22]
  R3_int <-  Week8to23.sol[nrow(Week8to23.sol),23]
  Total_inf3_int <- Week8to23.sol[nrow(Week8to23.sol),24]
  Total_Rec3_int <- Week8to23.sol[nrow(Week8to23.sol),25]
  D3_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),26]
  D3_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),27]
  WF2Hog_Count3_int <- Week8to23.sol[nrow(Week8to23.sol),28]
  
  S4_int <- Week8to23.sol[nrow(Week8to23.sol),29]
  E4_int <- Week8to23.sol[nrow(Week8to23.sol),30]
  I4_int <- Week8to23.sol[nrow(Week8to23.sol),31]
  R4_int <- Week8to23.sol[nrow(Week8to23.sol),32]
  Total_inf4_int <- Week8to23.sol[nrow(Week8to23.sol),33]
  Total_Rec4_int <- Week8to23.sol[nrow(Week8to23.sol),34]
  D4_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),35]
  D4_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),36]
  WF2Hog_Count4_int <- Week8to23.sol[nrow(Week8to23.sol),37]
  
  HS_int <-  Week8to23.sol[nrow(Week8to23.sol),38]
  HE_int <-  Week8to23.sol[nrow(Week8to23.sol),39]
  HI_int <-  Week8to23.sol[nrow(Week8to23.sol),40]
  HR_int <-  Week8to23.sol[nrow(Week8to23.sol),41]
  H_DN_int <-  Week8to23.sol[nrow(Week8to23.sol),42]
  H_DI_int <-  Week8to23.sol[nrow(Week8to23.sol),43]
  Hog2WF_Count_int <-  Week8to23.sol[nrow(Week8to23.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 23 to 25
  #====================================================================================================
  Week23to25.t <- seq(161.1,Week23to25,0.1)         
  Week23to25.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week23to25.sol <- lsoda(Week23to25.init, Week23to25.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 25 TO 27 (Room 3 Exit) ####
  Week25to27 <- 189
  #====================================================================================================
  # Model Initial Values  Weeks 25 to 27
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week23to25.sol[nrow(Week23to25.sol),6]
  Total_Rec1_int <- Week23to25.sol[nrow(Week23to25.sol),7]
  D1_Nat_int <-  Week23to25.sol[nrow(Week23to25.sol),8]
  D1_Inf_int <-  Week23to25.sol[nrow(Week23to25.sol),9]
  WF2Hog_Count1_int <-  Week23to25.sol[nrow(Week23to25.sol),10]
  
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- Week23to25.sol[nrow(Week23to25.sol),15]
  Total_Rec2_int <- Week23to25.sol[nrow(Week23to25.sol),16]
  D2_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),17]
  D2_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),18]
  WF2Hog_Count2_int <- Week23to25.sol[nrow(Week23to25.sol),19]
  
  S3_int <- Week23to25.sol[nrow(Week23to25.sol),20]
  E3_int <- Week23to25.sol[nrow(Week23to25.sol),21]
  I3_int <- Week23to25.sol[nrow(Week23to25.sol),22]
  R3_int <-  Week23to25.sol[nrow(Week23to25.sol),23]
  Total_inf3_int <- Week23to25.sol[nrow(Week23to25.sol),24]
  Total_Rec3_int <- Week23to25.sol[nrow(Week23to25.sol),25]
  D3_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),26]
  D3_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),27]
  WF2Hog_Count3_int <- Week23to25.sol[nrow(Week23to25.sol),28]
  
  S4_int <- Week23to25.sol[nrow(Week23to25.sol),29]
  E4_int <- Week23to25.sol[nrow(Week23to25.sol),30]
  I4_int <- Week23to25.sol[nrow(Week23to25.sol),31]
  R4_int <- Week23to25.sol[nrow(Week23to25.sol),32]
  Total_inf4_int <- Week23to25.sol[nrow(Week23to25.sol),33]
  Total_Rec4_int <- Week23to25.sol[nrow(Week23to25.sol),34]
  D4_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),35]
  D4_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),36]
  WF2Hog_Count4_int <- Week23to25.sol[nrow(Week23to25.sol),37]
  
  HS_int <-  Week23to25.sol[nrow(Week23to25.sol),38]
  HE_int <-  Week23to25.sol[nrow(Week23to25.sol),39]
  HI_int <-  Week23to25.sol[nrow(Week23to25.sol),40]
  HR_int <-  Week23to25.sol[nrow(Week23to25.sol),41]
  H_DN_int <-  Week23to25.sol[nrow(Week23to25.sol),42]
  H_DI_int <-  Week23to25.sol[nrow(Week23to25.sol),43]
  Hog2WF_Count_int <-  Week23to25.sol[nrow(Week23to25.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 25 to 27
  #====================================================================================================
  Week25to27.t <- seq(175.1,Week25to27,0.1)          
  Week25to27.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week25to27.sol <- lsoda(Week25to27.init, Week25to27.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 27 TO 29 (Room 4 Exit) ####
  Week27to29 <- 203
  #====================================================================================================
  # Model Initial Values  Weeks 27to 29
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week25to27.sol[nrow(Week25to27.sol),6]
  Total_Rec1_int <- Week25to27.sol[nrow(Week25to27.sol),7]
  D1_Nat_int <-  Week25to27.sol[nrow(Week25to27.sol),8]
  D1_Inf_int <-  Week25to27.sol[nrow(Week25to27.sol),9]
  WF2Hog_Count1_int <-  Week25to27.sol[nrow(Week25to27.sol),10]
  
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- Week25to27.sol[nrow(Week25to27.sol),15]
  Total_Rec2_int <- Week25to27.sol[nrow(Week25to27.sol),16]
  D2_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),17]
  D2_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),18]
  WF2Hog_Count2_int <- Week25to27.sol[nrow(Week25to27.sol),19]
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- Week25to27.sol[nrow(Week25to27.sol),24]
  Total_Rec3_int <- Week25to27.sol[nrow(Week25to27.sol),25]
  D3_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),26]
  D3_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),27]
  WF2Hog_Count3_int <- Week25to27.sol[nrow(Week25to27.sol),28]
  
  S4_int <- Week25to27.sol[nrow(Week25to27.sol),29]
  E4_int <- Week25to27.sol[nrow(Week25to27.sol),30]
  I4_int <- Week25to27.sol[nrow(Week25to27.sol),31]
  R4_int <- Week25to27.sol[nrow(Week25to27.sol),32]
  Total_inf4_int <- Week25to27.sol[nrow(Week25to27.sol),33]
  Total_Rec4_int <- Week25to27.sol[nrow(Week25to27.sol),34]
  D4_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),35]
  D4_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),36]
  WF2Hog_Count4_int <- Week25to27.sol[nrow(Week25to27.sol),37]
  
  HS_int <-  Week25to27.sol[nrow(Week25to27.sol),38]
  HE_int <-  Week25to27.sol[nrow(Week25to27.sol),39]
  HI_int <-  Week25to27.sol[nrow(Week25to27.sol),40]
  HR_int <-  Week25to27.sol[nrow(Week25to27.sol),41]
  H_DN_int <-  Week25to27.sol[nrow(Week25to27.sol),42]
  H_DI_int <-  Week25to27.sol[nrow(Week25to27.sol),43]
  Hog2WF_Count_int <-  Week25to27.sol[nrow(Week25to27.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 27 to 29
  #====================================================================================================
  Week27to29.t <- seq(189.1,Week27to29,0.1)            
  Week27to29.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week27to29.sol <- lsoda(Week27to29.init, Week27to29.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### COMBINING THE DATASETS ####
  Week0to2.sol.df <- as.data.frame(Week0to2.sol)
  colnames(Week0to2.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week2to4.sol.df <- as.data.frame(Week2to4.sol)
  colnames(Week2to4.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week4to6.sol.df <- as.data.frame(Week4to6.sol)
  colnames(Week4to6.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week6to8.sol.df <- as.data.frame(Week6to8.sol)
  colnames(Week6to8.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week8to23.sol.df <- as.data.frame(Week8to23.sol)
  colnames(Week8to23.sol.df) <- c("time",
                                  "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                  "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                  "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                  "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                  "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week23to25.sol.df <- as.data.frame(Week23to25.sol)
  colnames(Week23to25.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week25to27.sol.df <- as.data.frame(Week25to27.sol)
  colnames(Week25to27.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week27to29.sol.df <- as.data.frame(Week27to29.sol)
  colnames(Week27to29.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  
  fulldata.df <- rbind(Week0to2.sol.df, Week2to4.sol.df, Week4to6.sol.df, Week6to8.sol.df, 
                       Week8to23.sol.df, Week23to25.sol.df, Week25to27.sol.df, Week27to29.sol.df)
  
  ## Is there a WF infection
  HI_inf <- fulldata.df %>%
    summarise(total_I = sum(HI)) %>%
    mutate(HI_Inf = ifelse(total_I>1, 1, 0)) %>%
    select(HI_Inf) %>% as.numeric()
  #Time to 1st WF infection
  Time2WFinf0 <- fulldata.df %>%
    arrange(time) %>%
    filter(HI > 0) %>%
    filter(row_number() == 1) %>%
    mutate(TimeSineInfIntro = time - 42) %>%
    select(TimeSineInfIntro) %>% as.numeric()
  
  Time2WFinf.5 <- fulldata.df %>%
    arrange(time) %>%
    filter(HI > 0.5) %>%
    filter(row_number() == 1) %>%
    mutate(TimeSineInfIntro = time - 42) %>%
    select(TimeSineInfIntro) %>% as.numeric()
  ## Time to Second Room Infection
  Timeto2room <- fulldata.df %>% 
    summarise(t_cross = time[which(E1 > 1 | 
                                     E2 > 1 |
                                     E3 > 1)]) %>% 
    summarise(min_t_cross = min(t_cross)) %>% as.numeric()
  ## Total Hogs infected
  TotHogsInf <- fulldata.df %>% 
    #select(t, iter, Total_inf1, Total_inf2, Total_inf3, Total_inf4) %>% 
    mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
           Tot_rec = Total_Rec1 + Total_Rec2 + Total_Rec3 + Total_Rec4) %>% 
    summarise(Tot_inf2 = max(Tot_inf),
              Tot_rec2 = max(Tot_rec))
  ## Minimum Time to Peak Infection
  Time2PeakInf <- fulldata.df %>% 
    mutate(S_tot = S1+S2+S3+S4,                             
           E_tot = E1+E2+E3+E4,
           I_tot = I1+I2+I3+I4,
           R_tot = R1+R2+R3+R4) %>%                                     
    summarise(max_I = max(I_tot),
              t_peak = time[which(I_tot == max(I_tot))]) %>% 
    summarise(min_t_peak = min(t_peak)) %>% as.numeric()
  ##Ro
  Ro <- fulldata.df %>% 
    arrange(time) %>% 
    mutate(Tot_rec = Total_Rec1 + Total_Rec2 + Total_Rec3 + Total_Rec4) %>%  
    filter(row_number() == n()) %>% 
    select(Tot_rec) %>% 
    mutate(Ro = ((log((4000 - Tot_rec)/4000) - log((4000 - (Hog_Vacc_Eff*4000))/4000)) /
                   ((((4000 - Tot_rec) / 4000) - (4000 - (Hog_Vacc_Eff*4000)) / 4000)))) %>% 
    select(Ro) %>% as.numeric()
  
  #Create Summary data frame
  fulldata.df.sum <- data.frame(iter = as.factor(i),
                                HI_inf = HI_inf,
                                Time2WFinf0 = Time2WFinf0,
                                Time2WFinf.5 = Time2WFinf.5,
                                Timeto2room = Timeto2room,
                                TotHogsInf = TotHogsInf$Tot_inf2,
                                TotHogsRec = TotHogsInf$Tot_rec2,
                                Time2PeakInf = Time2PeakInf,
                                Ro = Ro,
                                beta_Hog = beta_Hog,
                                beta_Hog_ind = beta_Hog_ind,
                                beta_S2WF = beta_S2WF,
                                beta_WF2S = beta_WF2S,
                                sigma_Hog = sigma_Hog,
                                delta_Hog = delta_Hog,
                                beta_WF2WF = beta_WF2WF,
                                sigma_WF = sigma_WF,
                                delta_WF = delta_WF)
  
  results <- rbind(results, fulldata.df.sum)
  
  rm(list=ls()[! ls() %in% c("results", "n.itr")])  #Deletes everything from environment expect growing results dataframe
}

saveRDS(results, file = "Current code/Control Meausres/Sensativity Analysis/WFFlowRoom4_Sum_SensAnlaysis.rds")

#### WF Flow Only (Room 3 Intro) ####
results <- NULL
n.itr <- 5000
# seed.num <- 616
seed.nums <- c(100:(100+n.itr))
for (i in c(1:n.itr)) {
  print(i)
  
  #### MODEL CONSTANTS ####
  #Hog constants
  beta_Hog <- rtri(n = 1, min = 0.001, max = 0.1, mode = (10/(1000*5)))
  beta_Hog_ind <- beta_Hog / 178 #(following Etbaigha et al 2018 paper)
  # beta_Hog_ind <- beta_Hog / 500
  Hog_Vacc_Eff <- 0
  u_Hog <- 0.00028     #natural death rate (following Etbaigha et al 2018 paper)
  u_inf_Hog <- 0 #0.1     #infected death rate (0 for now to problem solve)
  w_Hog <- 0 #1/180         #(following Etbaigha et al 2018 paper)  White et al shows a range from 56 - 112 days
  sigma_Hog <- abs(1 / rnorm(n = 1, mean = 2, sd = 1))   #Latent period h1n1 and h3n2
  delta_Hog <- abs(1 / rnorm(n = 1, mean = 5, sd = 1))   #Infectious period h1n1 and h3n2
  
  #Interspecies transmission
  ##applied the beta = Ro / N*D from the modeling book
  # Ro from the following website for swine to human infection 
  # https://bmcmedicine.biomedcentral.com/articles/10.1186/1741-7015-7-30#Sec6   (from super-strain figure) possible the high end of transmisilibty for this parameter
  # beta_S2WF <- (2.3 / (4002*5))  #  #Hogs + #Workforce * Duration of Hogs
  #Alt beta_S2WF from fair outbreak
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3557885/
  #assuming one hog got the first three confirmed cases sick and infectious period of 5 days Ro = 3/5 = 0.6
  beta_S2WF <- runif(1, min = 0.000001, max = 0.0001)
  #when Ro WF2S = Ro S2WF
  # beta_WF2S <- (2.3 / (4002*3))  #  #Hogs + #Workforce * Duration of Hogs
  # WF2S calculated from Swine Outbreak of Pandemic Influenza A Virus on a Canadian Research Farm Supports Human-to-Swine Transmission paper
  beta_WF2S <- runif(1, min = 0.000001, max = 0.0001)    #Hogs + #Workforce * Duration of WF
  #Human constants
  #applied the beta = Ro / N*D from the modeling book
  beta_WF2WF <-   rtri(n = 1, min = 0.0000001, max = 0.64, mode = 0.32) # http://m-hikari.com/ams/ams-2013/ams-41-44-2013/joseAMS41-44-2013.pdf 
  Hum_Vacc_Eff <- 0.0 
  u_WF <- 0     #Natural death rate of human (0 for now to problem solve)
  u_inf_WF <- 0   #Infected death rate of human (0 for now to problem solve)
  w_WF <- 0      #Recovery rate for human (0 for now to problem solve)
  sigma_WF <- abs(1 / rnorm(n = 1, mean = 2, sd = 1))   #Latent period for humans
  delta_WF <- abs(1 / rnorm(n = 1, mean = 3, sd = 1))   #Infectious Period for humans
  
  #### MODEL SETUP ####
  # Model Parameters
  TransmissionModel.par <- c(beta_Hog = beta_Hog, beta_Hog_ind = beta_Hog_ind, beta_S2WF = beta_S2WF, beta_WF2S = beta_WF2S,
                             w_Hog = w_Hog, sigma_Hog = sigma_Hog, delta_Hog = delta_Hog, u_Hog = u_Hog, u_inf_Hog = u_inf_Hog,
                             beta_WF2WF = beta_WF2WF,
                             w_WF = w_WF, sigma_WF = sigma_WF, delta_WF = delta_WF, u_WF = u_WF, u_inf_WF = u_inf_WF)
  #=====================================================================================================
  #Model Equations 
  #=====================================================================================================
  TransmissionModel.dyn <- function(t,var,par) {
    # Rename the parameters 
    beta_Hog <- par[1]
    beta_Hog_ind <- par[2]
    beta_S2WF <- par[3]
    beta_WF2S <- par[4]
    w_Hog <- par[5]
    sigma_Hog <- par[6]
    delta_Hog <- par[7]
    u_Hog <- par[8]
    u_inf_Hog<- par[9]
    beta_WF2WF <- par[10]
    w_WF <- par[11]
    sigma_WF <- par[12]
    delta_WF <- par[13]
    u_WF <- par[14]
    u_inf_WF <- par[15]
    
    S1 <- var[1]
    E1 <- var[2]
    I1 <- var[3]
    R1 <- var[4]
    Total_inf1 <- var[5] 
    True_Rec1 <- var[6]
    D1_Nat <- var[7]
    D1_Inf <- var[8]
    WF2Hog_Count1 <- var[9] 
    
    S2 <- var[10]
    E2 <- var[11]
    I2 <- var[12]
    R2 <- var[13]
    Total_inf2 <- var[14]
    True_Rec2 <- var[15]
    D2_Nat <- var[16]
    D2_Inf <- var[17]
    WF2Hog_Count2 <- var[18]
    
    S3 <- var[19]
    E3 <- var[20]
    I3 <- var[21]
    R3 <- var[22]
    Total_inf3 <- var[23]
    True_Rec3 <- var[24]
    D3_Nat <- var[25]
    D3_Inf <- var[26]
    WF2Hog_Count3 <- var[27]   
    
    S4 <- var[28]
    E4 <- var[29]
    I4 <- var[30]
    R4 <- var[31]
    Total_inf4 <- var[32]
    True_Rec4 <- var[33]
    D4_Nat <- var[34]
    D4_Inf <- var[35]
    WF2Hog_Count4 <- var[36]
    
    HS <- var[37]
    HE <- var[38]
    HI <- var[39] 
    HR <- var[40]
    H_DN <- var[41]
    H_DI <- var[42]
    Hog2WF_Count <- var[43] 
    
    # Calculate the derivatives
    #Room 1
    dS1 <- (w_Hog*R1) - (beta_Hog*I1*S1) - (beta_Hog_ind*(I2+I3+I4)*S1) - (beta_WF2S*HI*S1) - (u_Hog*S1)
    dE1 <-(beta_Hog*I1*S1) + (beta_Hog_ind*(I2+I3+I4)*S1) + (beta_WF2S*HI*S1) - sigma_Hog*E1 - (u_Hog*E1)
    dI1 <- (sigma_Hog*E1) - (delta_Hog*I1) - (u_inf_Hog*I1)
    dR1 <- (delta_Hog*I1) - (w_Hog*R1) - (u_Hog*R1)
    dTotal_inf1 <- (sigma_Hog*E1)
    dTotal_Rec1 <- (delta_Hog*I1)
    dD1_Nat <- (u_Hog*S1) + (u_Hog*E1) + (u_Hog*R1)
    dD1_Inf <- (u_inf_Hog*I1)
    dWF2Hog_Count1 <- (beta_WF2S*HI*S1)
    
    #Room 2
    dS2 <- (w_Hog*R2) - (beta_Hog*I2*S2) - (beta_Hog_ind*(I3+I4)*S2) - (beta_WF2S*HI*S2) - (u_Hog*S2)
    dE2 <-(beta_Hog*I2*S2) + (beta_Hog_ind*(I3+I4)*S2) + (beta_WF2S*HI*S2) - (sigma_Hog*E2) - (u_Hog*E2) 
    dI2 <- (sigma_Hog*E2) - (delta_Hog*I2) - (u_inf_Hog*I2)
    dR2 <- (delta_Hog*I2) - (w_Hog*R2) - (u_Hog*R2)
    dTotal_inf2 <- (sigma_Hog*E2)
    dTotal_Rec2 <- (delta_Hog*I2)
    dD2_Nat <- (u_Hog*S2) + (u_Hog*E2) + (u_Hog*R2)
    dD2_Inf <- (u_inf_Hog*I2)
    dWF2Hog_Count2 <- (beta_WF2S*HI*S2)
    
    #Room 3
    dS3 <- (w_Hog*R3) - (beta_Hog*I3*S3) - (beta_Hog_ind*(I4)*S3) - (beta_WF2S*HI*S3) - (u_Hog*S3)
    dE3 <-(beta_Hog*I3*S3) + (beta_Hog_ind*(I4)*S3) + (beta_WF2S*HI*S3) - (sigma_Hog*E3) - (u_Hog*E3)
    dI3 <- (sigma_Hog*E3) - (delta_Hog*I3) - (u_inf_Hog*I3)
    dR3 <- (delta_Hog*I3) - (w_Hog*R3) - (u_Hog*R3)
    dTotal_inf3 <- (sigma_Hog*E3)
    dTotal_Rec3 <- (delta_Hog*I3)
    dD3_Nat <- (u_Hog*S3) + (u_Hog*E3) + (u_Hog*R3)
    dD3_Inf <- (u_inf_Hog*I3)
    dWF2Hog_Count3 <- (beta_WF2S*HI*S3)
    
    #Room 4
    dS4 <- w_Hog*R4 - (beta_Hog*I4*S4) - (beta_Hog_ind*(0)*S4) - (beta_WF2S*HI*S4) - (u_Hog*S4)
    dE4 <-(beta_Hog*I4*S4) + (beta_Hog_ind*(0)*S4) + (beta_WF2S*HI*S4) - (sigma_Hog*E4) - (u_Hog*E4) 
    dI4 <- (sigma_Hog*E4) - (delta_Hog*I4) - (u_inf_Hog*I4)
    dR4 <- (delta_Hog*I4) - (w_Hog*R4) - (u_Hog*R4)
    dTotal_inf4 <- (sigma_Hog*E4)
    dTotal_Rec4 <- (delta_Hog*I4)
    dD4_Nat <- (u_Hog*S4) + (u_Hog*E4) + (u_Hog*R4)
    dD4_Inf <- (u_inf_Hog*I4)
    dWF2Hog_Count4 <- (beta_WF2S*HI*S4)
    
    #Humans
    dHS <- (w_WF*HR)- (beta_WF2WF*HI*HS) - (beta_S2WF*(I1+I2+I3+I4)*HS) - (u_WF*HS)
    dHE <- (beta_WF2WF*HI*HS) + (beta_S2WF*(I1+I2+I3+I4)*HS) - (sigma_WF*HE) - (u_WF*HE)
    dHI <- (sigma_WF*HE) - (delta_WF*HI) - (u_inf_WF*HI)
    dHR <- (delta_WF*HI) - (w_WF*HR) - u_WF*HR
    dH_DN <- (u_WF*HS) + (u_WF*HE) + (u_WF*HR)
    dH_DI <- (u_inf_WF*HI)
    dHog2WF_Count <- (beta_S2WF*(I1+I2+I3+I4)*HS)
    
    
    # Last instruction: return a list 
    
    return(list(c(dS1, dE1, dI1, dR1, dTotal_inf1, dTotal_Rec1, dD1_Nat, dD1_Inf, dWF2Hog_Count1,
                  dS2, dE2, dI2, dR2, dTotal_inf2, dTotal_Rec2, dD2_Nat, dD2_Inf, dWF2Hog_Count2,
                  dS3, dE3, dI3, dR3, dTotal_inf3, dTotal_Rec3, dD3_Nat, dD3_Inf, dWF2Hog_Count3,
                  dS4, dE4, dI4, dR4, dTotal_inf4, dTotal_Rec4, dD4_Nat, dD4_Inf, dWF2Hog_Count4,
                  dHS, dHE, dHI, dHR, dH_DN, dH_DI, dHog2WF_Count)))
    
  }
  
  #### MODEL ITERATIONS ####
  #### WEEK 0 TO 2 (Room 1 fill) ####
  Week0to2 <- 14
  #====================================================================================================
  # Model Initial Values  Weeks 0 to 2
  #===========================================================================
  S1_int <- (1-Hog_Vacc_Eff)*1000
  E1_int <- 0
  I1_int <- 0 #1
  R1_int <- Hog_Vacc_Eff*1000
  Total_inf1_int <- 0
  Total_Rec1_int <- 0
  D1_Nat_int <- 0
  D1_Inf_int <- 0
  WF2Hog_Count1_int <- 0
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- 0
  Total_Rec2_int <- 0
  D2_Nat_int <- 0
  D2_Inf_int <- 0
  WF2Hog_Count2_int <- 0
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <- 2
  HE_int <- 0
  HI_int <- 0
  HR_int <- 0
  H_DN_int <- 0
  H_DI_int <- 0
  Hog2WF_Count_int <- 0
  
  #====================================================================================================
  # Model Run  Weeks 0 to 2
  #====================================================================================================
  Week0to2.t <- seq(0,Week0to2,0.1)   
  Week0to2.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week0to2.sol <- lsoda(Week0to2.init, Week0to2.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 2 TO 4 (Room 2 Fill)####
  Week2to4 <- 28
  #====================================================================================================
  # Model Initial Values  Weeks 2 to 4
  #===========================================================================
  S1_int <- Week0to2.sol[nrow(Week0to2.sol),2]
  E1_int <- Week0to2.sol[nrow(Week0to2.sol),3]
  I1_int <- Week0to2.sol[nrow(Week0to2.sol),4]
  R1_int <- Week0to2.sol[nrow(Week0to2.sol),5]
  Total_inf1_int <-  Week0to2.sol[nrow(Week0to2.sol),6]
  Total_Rec1_int <-  Week0to2.sol[nrow(Week0to2.sol),7]
  D1_Nat_int <-  Week0to2.sol[nrow(Week0to2.sol),8]
  D1_Inf_int <-  Week0to2.sol[nrow(Week0to2.sol),9]
  WF2Hog_Count1_int <-  Week0to2.sol[nrow(Week0to2.sol),10]
  
  S2_int <- (1-Hog_Vacc_Eff)*1000
  E2_int <- 0
  I2_int <- 0 #1
  R2_int <- Hog_Vacc_Eff*1000
  Total_inf2_int <- 0
  Total_Rec2_int <- 0
  D2_Nat_int <- 0
  D2_Inf_int <- 0
  WF2Hog_Count2_int <- 0
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week0to2.sol[nrow(Week0to2.sol),38]
  HE_int <-  Week0to2.sol[nrow(Week0to2.sol),39]
  HI_int <-  Week0to2.sol[nrow(Week0to2.sol),40]
  HR_int <-  Week0to2.sol[nrow(Week0to2.sol),41]
  H_DN_int <-  Week0to2.sol[nrow(Week0to2.sol),42]
  H_DI_int <-  Week0to2.sol[nrow(Week0to2.sol),43]
  Hog2WF_Count_int <-  Week0to2.sol[nrow(Week0to2.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 2 to 4
  #====================================================================================================
  Week2to4.t <- seq(14.1,Week2to4,0.1)          
  Week2to4.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week2to4.sol <- lsoda(Week2to4.init, Week2to4.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 4 TO 6 (Room 3 Fill)####
  Week4to6 <- 42
  #====================================================================================================
  # Model Initial Values  Weeks 4 to 6
  #===========================================================================
  S1_int <- Week2to4.sol[nrow(Week2to4.sol),2]
  E1_int <- Week2to4.sol[nrow(Week2to4.sol),3]
  I1_int <- Week2to4.sol[nrow(Week2to4.sol),4]
  R1_int <- Week2to4.sol[nrow(Week2to4.sol),5]
  Total_inf1_int <- Week2to4.sol[nrow(Week2to4.sol),6]
  Total_Rec1_int <- Week2to4.sol[nrow(Week2to4.sol),7]
  D1_Nat_int <-  Week2to4.sol[nrow(Week2to4.sol),8]
  D1_Inf_int <-  Week2to4.sol[nrow(Week2to4.sol),9]
  WF2Hog_Count1_int <-  Week2to4.sol[nrow(Week2to4.sol),10]
  
  S2_int <- Week2to4.sol[nrow(Week2to4.sol),11]
  E2_int <- Week2to4.sol[nrow(Week2to4.sol),12]
  I2_int <- Week2to4.sol[nrow(Week2to4.sol),13]
  R2_int <- Week2to4.sol[nrow(Week2to4.sol),14]
  Total_inf2_int <- Week2to4.sol[nrow(Week2to4.sol),15]
  Total_Rec2_int <- Week2to4.sol[nrow(Week2to4.sol),16]
  D2_Nat_int <- Week2to4.sol[nrow(Week2to4.sol),17]
  D2_Inf_int <- Week2to4.sol[nrow(Week2to4.sol),18]
  WF2Hog_Count2_int <- Week2to4.sol[nrow(Week2to4.sol),19]
  
  S3_int <- (1-Hog_Vacc_Eff)*999
  E3_int <- 0
  I3_int <- 1
  R3_int <-  Hog_Vacc_Eff*999
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week2to4.sol[nrow(Week2to4.sol),38]
  HE_int <-  Week2to4.sol[nrow(Week2to4.sol),39]
  HI_int <-  Week2to4.sol[nrow(Week2to4.sol),40]
  HR_int <-  Week2to4.sol[nrow(Week2to4.sol),41]
  H_DN_int <-  Week2to4.sol[nrow(Week2to4.sol),42]
  H_DI_int <-  Week2to4.sol[nrow(Week2to4.sol),43]
  Hog2WF_Count_int <-  Week2to4.sol[nrow(Week2to4.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 4 to 6
  #====================================================================================================
  Week4to6.t <- seq(28.1,Week4to6,0.1)         
  Week4to6.init <-c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                    S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                    S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                    S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                    HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week4to6.sol <- lsoda(Week4to6.init, Week4to6.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 6 TO 8 (Room 4 Fill)####
  Week6to8 <- 56
  #====================================================================================================
  # Model Initial Values  Weeks 6 to 8
  #===========================================================================
  S1_int <- Week4to6.sol[nrow(Week4to6.sol),2]
  E1_int <- Week4to6.sol[nrow(Week4to6.sol),3]
  I1_int <- Week4to6.sol[nrow(Week4to6.sol),4]
  R1_int <- Week4to6.sol[nrow(Week4to6.sol),5]
  Total_inf1_int <- Week4to6.sol[nrow(Week4to6.sol),6]
  Total_Rec1_int <- Week4to6.sol[nrow(Week4to6.sol),7]
  D1_Nat_int <-  Week4to6.sol[nrow(Week4to6.sol),8]
  D1_Inf_int <-  Week4to6.sol[nrow(Week4to6.sol),9]
  WF2Hog_Count1_int <-  Week4to6.sol[nrow(Week4to6.sol),10]
  
  S2_int <- Week4to6.sol[nrow(Week4to6.sol),11]
  E2_int <- Week4to6.sol[nrow(Week4to6.sol),12]
  I2_int <- Week4to6.sol[nrow(Week4to6.sol),13]
  R2_int <- Week4to6.sol[nrow(Week4to6.sol),14]
  Total_inf2_int <- Week4to6.sol[nrow(Week4to6.sol),15]
  Total_Rec2_int <- Week4to6.sol[nrow(Week4to6.sol),16]
  D2_Nat_int <- Week4to6.sol[nrow(Week4to6.sol),17]
  D2_Inf_int <- Week4to6.sol[nrow(Week4to6.sol),18]
  WF2Hog_Count2_int <- Week4to6.sol[nrow(Week4to6.sol),19]
  
  S3_int <- Week4to6.sol[nrow(Week4to6.sol),20]
  E3_int <- Week4to6.sol[nrow(Week4to6.sol),21]
  I3_int <- Week4to6.sol[nrow(Week4to6.sol),22]
  R3_int <-  Week4to6.sol[nrow(Week4to6.sol),23]
  Total_inf3_int <- Week4to6.sol[nrow(Week4to6.sol),24]
  Total_Rec3_int <- Week4to6.sol[nrow(Week4to6.sol),25]
  D3_Nat_int <- Week4to6.sol[nrow(Week4to6.sol),26]
  D3_Inf_int <- Week4to6.sol[nrow(Week4to6.sol),27]
  WF2Hog_Count3_int <- Week4to6.sol[nrow(Week4to6.sol),28]
  
  S4_int <- (1-Hog_Vacc_Eff)*1000
  E4_int <- 0
  I4_int <- 0
  R4_int <- Hog_Vacc_Eff*1000
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week4to6.sol[nrow(Week4to6.sol),38]
  HE_int <-  Week4to6.sol[nrow(Week4to6.sol),39]
  HI_int <-  Week4to6.sol[nrow(Week4to6.sol),40]
  HR_int <-  Week4to6.sol[nrow(Week4to6.sol),41]
  H_DN_int <-  Week4to6.sol[nrow(Week4to6.sol),42]
  H_DI_int <-  Week4to6.sol[nrow(Week4to6.sol),43]
  Hog2WF_Count_int <-  Week4to6.sol[nrow(Week4to6.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 6 to 8
  #====================================================================================================
  Week6to8.t <- seq(42.1,Week6to8,0.1)           
  Week6to8.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week6to8.sol <- lsoda(Week6to8.init, Week6to8.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 8 TO 23 (Room 1 Exit) ####
  Week8to23 <- 161
  #====================================================================================================
  # Model Initial Values  Weeks 8 to 23
  #===========================================================================
  S1_int <- Week6to8.sol[nrow(Week6to8.sol),2]
  E1_int <- Week6to8.sol[nrow(Week6to8.sol),3]
  I1_int <- Week6to8.sol[nrow(Week6to8.sol),4]
  R1_int <- Week6to8.sol[nrow(Week6to8.sol),5]
  Total_inf1_int <- Week6to8.sol[nrow(Week6to8.sol),6]
  Total_Rec1_int <- Week6to8.sol[nrow(Week6to8.sol),7]
  D1_Nat_int <-  Week6to8.sol[nrow(Week6to8.sol),8]
  D1_Inf_int <-  Week6to8.sol[nrow(Week6to8.sol),9]
  WF2Hog_Count1_int <-  Week6to8.sol[nrow(Week6to8.sol),10]
  
  
  S2_int <- Week6to8.sol[nrow(Week6to8.sol),11]
  E2_int <- Week6to8.sol[nrow(Week6to8.sol),12]
  I2_int <- Week6to8.sol[nrow(Week6to8.sol),13]
  R2_int <- Week6to8.sol[nrow(Week6to8.sol),14]
  Total_inf2_int <- Week6to8.sol[nrow(Week6to8.sol),15]
  Total_Rec2_int <- Week6to8.sol[nrow(Week6to8.sol),16]
  D2_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),17]
  D2_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),18]
  WF2Hog_Count2_int <- Week6to8.sol[nrow(Week6to8.sol),19]
  
  S3_int <- Week6to8.sol[nrow(Week6to8.sol),20]
  E3_int <- Week6to8.sol[nrow(Week6to8.sol),21]
  I3_int <- Week6to8.sol[nrow(Week6to8.sol),22]
  R3_int <-  Week6to8.sol[nrow(Week6to8.sol),23]
  Total_inf3_int <- Week6to8.sol[nrow(Week6to8.sol),24]
  Total_Rec3_int <- Week6to8.sol[nrow(Week6to8.sol),25]
  D3_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),26]
  D3_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),27]
  WF2Hog_Count3_int <- Week6to8.sol[nrow(Week6to8.sol),28]
  
  S4_int <- Week6to8.sol[nrow(Week6to8.sol),29]
  E4_int <- Week6to8.sol[nrow(Week6to8.sol),30]
  I4_int <- Week6to8.sol[nrow(Week6to8.sol),31]
  R4_int <- Week6to8.sol[nrow(Week6to8.sol),32]
  Total_inf4_int <- Week6to8.sol[nrow(Week6to8.sol),33]
  Total_Rec4_int <- Week6to8.sol[nrow(Week6to8.sol),34]
  D4_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),35]
  D4_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),36]
  WF2Hog_Count4_int <- Week6to8.sol[nrow(Week6to8.sol),37]
  
  HS_int <-  Week6to8.sol[nrow(Week6to8.sol),38]
  HE_int <-  Week6to8.sol[nrow(Week6to8.sol),39]
  HI_int <-  Week6to8.sol[nrow(Week6to8.sol),40]
  HR_int <-  Week6to8.sol[nrow(Week6to8.sol),41]
  H_DN_int <-  Week6to8.sol[nrow(Week6to8.sol),42]
  H_DI_int <-  Week6to8.sol[nrow(Week6to8.sol),43]
  Hog2WF_Count_int <-  Week6to8.sol[nrow(Week6to8.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 8 to 23
  #====================================================================================================
  Week8to23.t <- seq(56.1,Week8to23,0.1)          
  Week8to23.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                      S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                      S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                      S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                      HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week8to23.sol <- lsoda(Week8to23.init, Week8to23.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 23 TO 25 (Room 2 Exit) ####
  Week23to25 <- 175
  #====================================================================================================
  # Model Initial Values  Weeks 23 to 25
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week8to23.sol[nrow(Week8to23.sol),6]
  Total_Rec1_int <- Week8to23.sol[nrow(Week8to23.sol),7]
  D1_Nat_int <-  Week8to23.sol[nrow(Week8to23.sol),8]
  D1_Inf_int <-  Week8to23.sol[nrow(Week8to23.sol),9]
  WF2Hog_Count1_int <-  Week8to23.sol[nrow(Week8to23.sol),10]
  
  
  S2_int <- Week8to23.sol[nrow(Week8to23.sol),11]
  E2_int <- Week8to23.sol[nrow(Week8to23.sol),12]
  I2_int <- Week8to23.sol[nrow(Week8to23.sol),13]
  R2_int <- Week8to23.sol[nrow(Week8to23.sol),14]
  Total_inf2_int <- Week8to23.sol[nrow(Week8to23.sol),15]
  Total_Rec2_int <- Week8to23.sol[nrow(Week8to23.sol),16]
  D2_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),17]
  D2_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),18]
  WF2Hog_Count2_int <- Week8to23.sol[nrow(Week8to23.sol),19]
  
  S3_int <- Week8to23.sol[nrow(Week8to23.sol),20]
  E3_int <- Week8to23.sol[nrow(Week8to23.sol),21]
  I3_int <- Week8to23.sol[nrow(Week8to23.sol),22]
  R3_int <-  Week8to23.sol[nrow(Week8to23.sol),23]
  Total_inf3_int <- Week8to23.sol[nrow(Week8to23.sol),24]
  Total_Rec3_int <- Week8to23.sol[nrow(Week8to23.sol),25]
  D3_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),26]
  D3_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),27]
  WF2Hog_Count3_int <- Week8to23.sol[nrow(Week8to23.sol),28]
  
  S4_int <- Week8to23.sol[nrow(Week8to23.sol),29]
  E4_int <- Week8to23.sol[nrow(Week8to23.sol),30]
  I4_int <- Week8to23.sol[nrow(Week8to23.sol),31]
  R4_int <- Week8to23.sol[nrow(Week8to23.sol),32]
  Total_inf4_int <- Week8to23.sol[nrow(Week8to23.sol),33]
  Total_Rec4_int <- Week8to23.sol[nrow(Week8to23.sol),34]
  D4_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),35]
  D4_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),36]
  WF2Hog_Count4_int <- Week8to23.sol[nrow(Week8to23.sol),37]
  
  HS_int <-  Week8to23.sol[nrow(Week8to23.sol),38]
  HE_int <-  Week8to23.sol[nrow(Week8to23.sol),39]
  HI_int <-  Week8to23.sol[nrow(Week8to23.sol),40]
  HR_int <-  Week8to23.sol[nrow(Week8to23.sol),41]
  H_DN_int <-  Week8to23.sol[nrow(Week8to23.sol),42]
  H_DI_int <-  Week8to23.sol[nrow(Week8to23.sol),43]
  Hog2WF_Count_int <-  Week8to23.sol[nrow(Week8to23.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 23 to 25
  #====================================================================================================
  Week23to25.t <- seq(161.1,Week23to25,0.1)         
  Week23to25.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week23to25.sol <- lsoda(Week23to25.init, Week23to25.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 25 TO 27 (Room 3 Exit) ####
  Week25to27 <- 189
  #====================================================================================================
  # Model Initial Values  Weeks 25 to 27
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week23to25.sol[nrow(Week23to25.sol),6]
  Total_Rec1_int <- Week23to25.sol[nrow(Week23to25.sol),7]
  D1_Nat_int <-  Week23to25.sol[nrow(Week23to25.sol),8]
  D1_Inf_int <-  Week23to25.sol[nrow(Week23to25.sol),9]
  WF2Hog_Count1_int <-  Week23to25.sol[nrow(Week23to25.sol),10]
  
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- Week23to25.sol[nrow(Week23to25.sol),15]
  Total_Rec2_int <- Week23to25.sol[nrow(Week23to25.sol),16]
  D2_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),17]
  D2_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),18]
  WF2Hog_Count2_int <- Week23to25.sol[nrow(Week23to25.sol),19]
  
  S3_int <- Week23to25.sol[nrow(Week23to25.sol),20]
  E3_int <- Week23to25.sol[nrow(Week23to25.sol),21]
  I3_int <- Week23to25.sol[nrow(Week23to25.sol),22]
  R3_int <-  Week23to25.sol[nrow(Week23to25.sol),23]
  Total_inf3_int <- Week23to25.sol[nrow(Week23to25.sol),24]
  Total_Rec3_int <- Week23to25.sol[nrow(Week23to25.sol),25]
  D3_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),26]
  D3_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),27]
  WF2Hog_Count3_int <- Week23to25.sol[nrow(Week23to25.sol),28]
  
  S4_int <- Week23to25.sol[nrow(Week23to25.sol),29]
  E4_int <- Week23to25.sol[nrow(Week23to25.sol),30]
  I4_int <- Week23to25.sol[nrow(Week23to25.sol),31]
  R4_int <- Week23to25.sol[nrow(Week23to25.sol),32]
  Total_inf4_int <- Week23to25.sol[nrow(Week23to25.sol),33]
  Total_Rec4_int <- Week23to25.sol[nrow(Week23to25.sol),34]
  D4_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),35]
  D4_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),36]
  WF2Hog_Count4_int <- Week23to25.sol[nrow(Week23to25.sol),37]
  
  HS_int <-  Week23to25.sol[nrow(Week23to25.sol),38]
  HE_int <-  Week23to25.sol[nrow(Week23to25.sol),39]
  HI_int <-  Week23to25.sol[nrow(Week23to25.sol),40]
  HR_int <-  Week23to25.sol[nrow(Week23to25.sol),41]
  H_DN_int <-  Week23to25.sol[nrow(Week23to25.sol),42]
  H_DI_int <-  Week23to25.sol[nrow(Week23to25.sol),43]
  Hog2WF_Count_int <-  Week23to25.sol[nrow(Week23to25.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 25 to 27
  #====================================================================================================
  Week25to27.t <- seq(175.1,Week25to27,0.1)          
  Week25to27.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week25to27.sol <- lsoda(Week25to27.init, Week25to27.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 27 TO 29 (Room 4 Exit) ####
  Week27to29 <- 203
  #====================================================================================================
  # Model Initial Values  Weeks 27to 29
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week25to27.sol[nrow(Week25to27.sol),6]
  Total_Rec1_int <- Week25to27.sol[nrow(Week25to27.sol),7]
  D1_Nat_int <-  Week25to27.sol[nrow(Week25to27.sol),8]
  D1_Inf_int <-  Week25to27.sol[nrow(Week25to27.sol),9]
  WF2Hog_Count1_int <-  Week25to27.sol[nrow(Week25to27.sol),10]
  
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- Week25to27.sol[nrow(Week25to27.sol),15]
  Total_Rec2_int <- Week25to27.sol[nrow(Week25to27.sol),16]
  D2_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),17]
  D2_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),18]
  WF2Hog_Count2_int <- Week25to27.sol[nrow(Week25to27.sol),19]
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- Week25to27.sol[nrow(Week25to27.sol),24]
  Total_Rec3_int <- Week25to27.sol[nrow(Week25to27.sol),25]
  D3_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),26]
  D3_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),27]
  WF2Hog_Count3_int <- Week25to27.sol[nrow(Week25to27.sol),28]
  
  S4_int <- Week25to27.sol[nrow(Week25to27.sol),29]
  E4_int <- Week25to27.sol[nrow(Week25to27.sol),30]
  I4_int <- Week25to27.sol[nrow(Week25to27.sol),31]
  R4_int <- Week25to27.sol[nrow(Week25to27.sol),32]
  Total_inf4_int <- Week25to27.sol[nrow(Week25to27.sol),33]
  Total_Rec4_int <- Week25to27.sol[nrow(Week25to27.sol),34]
  D4_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),35]
  D4_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),36]
  WF2Hog_Count4_int <- Week25to27.sol[nrow(Week25to27.sol),37]
  
  HS_int <-  Week25to27.sol[nrow(Week25to27.sol),38]
  HE_int <-  Week25to27.sol[nrow(Week25to27.sol),39]
  HI_int <-  Week25to27.sol[nrow(Week25to27.sol),40]
  HR_int <-  Week25to27.sol[nrow(Week25to27.sol),41]
  H_DN_int <-  Week25to27.sol[nrow(Week25to27.sol),42]
  H_DI_int <-  Week25to27.sol[nrow(Week25to27.sol),43]
  Hog2WF_Count_int <-  Week25to27.sol[nrow(Week25to27.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 27 to 29
  #====================================================================================================
  Week27to29.t <- seq(189.1,Week27to29,0.1)            
  Week27to29.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week27to29.sol <- lsoda(Week27to29.init, Week27to29.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### COMBINING THE DATASETS ####
  Week0to2.sol.df <- as.data.frame(Week0to2.sol)
  colnames(Week0to2.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week2to4.sol.df <- as.data.frame(Week2to4.sol)
  colnames(Week2to4.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week4to6.sol.df <- as.data.frame(Week4to6.sol)
  colnames(Week4to6.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week6to8.sol.df <- as.data.frame(Week6to8.sol)
  colnames(Week6to8.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week8to23.sol.df <- as.data.frame(Week8to23.sol)
  colnames(Week8to23.sol.df) <- c("time",
                                  "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                  "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                  "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                  "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                  "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week23to25.sol.df <- as.data.frame(Week23to25.sol)
  colnames(Week23to25.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week25to27.sol.df <- as.data.frame(Week25to27.sol)
  colnames(Week25to27.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week27to29.sol.df <- as.data.frame(Week27to29.sol)
  colnames(Week27to29.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  
  fulldata.df <- rbind(Week0to2.sol.df, Week2to4.sol.df, Week4to6.sol.df, Week6to8.sol.df, 
                       Week8to23.sol.df, Week23to25.sol.df, Week25to27.sol.df, Week27to29.sol.df)
  
  ## Is there a WF infection
  HI_inf <- fulldata.df %>%
    summarise(total_I = sum(HI)) %>%
    mutate(HI_Inf = ifelse(total_I>1, 1, 0)) %>%
    select(HI_Inf) %>% as.numeric()
  #Time to 1st WF infection
  Time2WFinf0 <- fulldata.df %>%
    arrange(time) %>%
    filter(HI > 0) %>%
    filter(row_number() == 1) %>%
    mutate(TimeSineInfIntro = time - 28) %>%
    select(TimeSineInfIntro) %>% as.numeric()
  
  Time2WFinf.5 <- fulldata.df %>%
    arrange(time) %>%
    filter(HI > 0.5) %>%
    filter(row_number() == 1) %>%
    mutate(TimeSineInfIntro = time - 28) %>%
    select(TimeSineInfIntro) %>% as.numeric()
  ## Time to Second Room Infection
  Timeto2room <- fulldata.df %>% 
    summarise(t_cross = time[which(E1 > 1 | 
                                     E2 > 1 |
                                     E4 > 1)]) %>% 
    summarise(min_t_cross = min(t_cross)) %>% as.numeric()
  ## Total Hogs infected
  TotHogsInf <- fulldata.df %>% 
    #select(t, iter, Total_inf1, Total_inf2, Total_inf3, Total_inf4) %>% 
    mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
           Tot_rec = Total_Rec1 + Total_Rec2 + Total_Rec3 + Total_Rec4) %>% 
    summarise(Tot_inf2 = max(Tot_inf),
              Tot_rec2 = max(Tot_rec))
  ## Minimum Time to Peak Infection
  Time2PeakInf <- fulldata.df %>% 
    mutate(S_tot = S1+S2+S3+S4,                             
           E_tot = E1+E2+E3+E4,
           I_tot = I1+I2+I3+I4,
           R_tot = R1+R2+R3+R4) %>%                                     
    summarise(max_I = max(I_tot),
              t_peak = time[which(I_tot == max(I_tot))]) %>% 
    summarise(min_t_peak = min(t_peak)) %>% as.numeric()
  ##Ro
  Ro <- fulldata.df %>% 
    arrange(time) %>% 
    mutate(Tot_rec = Total_Rec1 + Total_Rec2 + Total_Rec3 + Total_Rec4) %>%  
    filter(row_number() == n()) %>% 
    select(Tot_rec) %>% 
    mutate(Ro = ((log((4000 - Tot_rec)/4000) - log((4000 - (Hog_Vacc_Eff*4000))/4000)) /
                   ((((4000 - Tot_rec) / 4000) - (4000 - (Hog_Vacc_Eff*4000)) / 4000)))) %>% 
    select(Ro) %>% as.numeric()
  
  #Create Summary data frame
  fulldata.df.sum <- data.frame(iter = as.factor(i),
                                HI_inf = HI_inf,
                                Time2WFinf0 = Time2WFinf0,
                                Time2WFinf.5 = Time2WFinf.5,
                                Timeto2room = Timeto2room,
                                TotHogsInf = TotHogsInf$Tot_inf2,
                                TotHogsRec = TotHogsInf$Tot_rec2,
                                Time2PeakInf = Time2PeakInf,
                                Ro = Ro,
                                beta_Hog = beta_Hog,
                                beta_Hog_ind = beta_Hog_ind,
                                beta_S2WF = beta_S2WF,
                                beta_WF2S = beta_WF2S,
                                sigma_Hog = sigma_Hog,
                                delta_Hog = delta_Hog,
                                beta_WF2WF = beta_WF2WF,
                                sigma_WF = sigma_WF,
                                delta_WF = delta_WF)
  
  results <- rbind(results, fulldata.df.sum)
  
  rm(list=ls()[! ls() %in% c("results", "n.itr")])  #Deletes everything from environment expect growing results dataframe
}

saveRDS(results, file = "Current code/Control Meausres/Sensativity Analysis/WFFlowRoom3_Sum_SensAnlaysis.rds")

#### WF Flow Only (Room 2 Intro) ####
results <- NULL
n.itr <- 5000
# seed.num <- 616
seed.nums <- c(100:(100+n.itr))
for (i in c(1:n.itr)) {
  print(i)
  
  #### MODEL CONSTANTS ####
  #Hog constants
  beta_Hog <- rtri(n = 1, min = 0.001, max = 0.1, mode = (10/(1000*5)))
  beta_Hog_ind <- beta_Hog / 178 #(following Etbaigha et al 2018 paper)
  # beta_Hog_ind <- beta_Hog / 500
  Hog_Vacc_Eff <- 0
  u_Hog <- 0.00028     #natural death rate (following Etbaigha et al 2018 paper)
  u_inf_Hog <- 0 #0.1     #infected death rate (0 for now to problem solve)
  w_Hog <- 0 #1/180         #(following Etbaigha et al 2018 paper)  White et al shows a range from 56 - 112 days
  sigma_Hog <- abs(1 / rnorm(n = 1, mean = 2, sd = 1))   #Latent period h1n1 and h3n2
  delta_Hog <- abs(1 / rnorm(n = 1, mean = 5, sd = 1))   #Infectious period h1n1 and h3n2
  
  #Interspecies transmission
  ##applied the beta = Ro / N*D from the modeling book
  # Ro from the following website for swine to human infection 
  # https://bmcmedicine.biomedcentral.com/articles/10.1186/1741-7015-7-30#Sec6   (from super-strain figure) possible the high end of transmisilibty for this parameter
  # beta_S2WF <- (2.3 / (4002*5))  #  #Hogs + #Workforce * Duration of Hogs
  #Alt beta_S2WF from fair outbreak
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3557885/
  #assuming one hog got the first three confirmed cases sick and infectious period of 5 days Ro = 3/5 = 0.6
  beta_S2WF <- runif(1, min = 0.000001, max = 0.0001)
  #when Ro WF2S = Ro S2WF
  # beta_WF2S <- (2.3 / (4002*3))  #  #Hogs + #Workforce * Duration of Hogs
  # WF2S calculated from Swine Outbreak of Pandemic Influenza A Virus on a Canadian Research Farm Supports Human-to-Swine Transmission paper
  beta_WF2S <- runif(1, min = 0.000001, max = 0.0001)    #Hogs + #Workforce * Duration of WF
  #Human constants
  #applied the beta = Ro / N*D from the modeling book
  beta_WF2WF <-   rtri(n = 1, min = 0.0000001, max = 0.64, mode = 0.32) # http://m-hikari.com/ams/ams-2013/ams-41-44-2013/joseAMS41-44-2013.pdf 
  Hum_Vacc_Eff <- 0.0 
  u_WF <- 0     #Natural death rate of human (0 for now to problem solve)
  u_inf_WF <- 0   #Infected death rate of human (0 for now to problem solve)
  w_WF <- 0      #Recovery rate for human (0 for now to problem solve)
  sigma_WF <- abs(1 / rnorm(n = 1, mean = 2, sd = 1))   #Latent period for humans
  delta_WF <- abs(1 / rnorm(n = 1, mean = 3, sd = 1))   #Infectious Period for humans
  
  #### MODEL SETUP ####
  # Model Parameters
  TransmissionModel.par <- c(beta_Hog = beta_Hog, beta_Hog_ind = beta_Hog_ind, beta_S2WF = beta_S2WF, beta_WF2S = beta_WF2S,
                             w_Hog = w_Hog, sigma_Hog = sigma_Hog, delta_Hog = delta_Hog, u_Hog = u_Hog, u_inf_Hog = u_inf_Hog,
                             beta_WF2WF = beta_WF2WF,
                             w_WF = w_WF, sigma_WF = sigma_WF, delta_WF = delta_WF, u_WF = u_WF, u_inf_WF = u_inf_WF)
  #=====================================================================================================
  #Model Equations 
  #=====================================================================================================
  TransmissionModel.dyn <- function(t,var,par) {
    # Rename the parameters 
    beta_Hog <- par[1]
    beta_Hog_ind <- par[2]
    beta_S2WF <- par[3]
    beta_WF2S <- par[4]
    w_Hog <- par[5]
    sigma_Hog <- par[6]
    delta_Hog <- par[7]
    u_Hog <- par[8]
    u_inf_Hog<- par[9]
    beta_WF2WF <- par[10]
    w_WF <- par[11]
    sigma_WF <- par[12]
    delta_WF <- par[13]
    u_WF <- par[14]
    u_inf_WF <- par[15]
    
    S1 <- var[1]
    E1 <- var[2]
    I1 <- var[3]
    R1 <- var[4]
    Total_inf1 <- var[5] 
    True_Rec1 <- var[6]
    D1_Nat <- var[7]
    D1_Inf <- var[8]
    WF2Hog_Count1 <- var[9] 
    
    S2 <- var[10]
    E2 <- var[11]
    I2 <- var[12]
    R2 <- var[13]
    Total_inf2 <- var[14]
    True_Rec2 <- var[15]
    D2_Nat <- var[16]
    D2_Inf <- var[17]
    WF2Hog_Count2 <- var[18]
    
    S3 <- var[19]
    E3 <- var[20]
    I3 <- var[21]
    R3 <- var[22]
    Total_inf3 <- var[23]
    True_Rec3 <- var[24]
    D3_Nat <- var[25]
    D3_Inf <- var[26]
    WF2Hog_Count3 <- var[27]   
    
    S4 <- var[28]
    E4 <- var[29]
    I4 <- var[30]
    R4 <- var[31]
    Total_inf4 <- var[32]
    True_Rec4 <- var[33]
    D4_Nat <- var[34]
    D4_Inf <- var[35]
    WF2Hog_Count4 <- var[36]
    
    HS <- var[37]
    HE <- var[38]
    HI <- var[39] 
    HR <- var[40]
    H_DN <- var[41]
    H_DI <- var[42]
    Hog2WF_Count <- var[43] 
    
    # Calculate the derivatives
    #Room 1
    dS1 <- (w_Hog*R1) - (beta_Hog*I1*S1) - (beta_Hog_ind*(I2+I3+I4)*S1) - (beta_WF2S*HI*S1) - (u_Hog*S1)
    dE1 <-(beta_Hog*I1*S1) + (beta_Hog_ind*(I2+I3+I4)*S1) + (beta_WF2S*HI*S1) - sigma_Hog*E1 - (u_Hog*E1)
    dI1 <- (sigma_Hog*E1) - (delta_Hog*I1) - (u_inf_Hog*I1)
    dR1 <- (delta_Hog*I1) - (w_Hog*R1) - (u_Hog*R1)
    dTotal_inf1 <- (sigma_Hog*E1)
    dTotal_Rec1 <- (delta_Hog*I1)
    dD1_Nat <- (u_Hog*S1) + (u_Hog*E1) + (u_Hog*R1)
    dD1_Inf <- (u_inf_Hog*I1)
    dWF2Hog_Count1 <- (beta_WF2S*HI*S1)
    
    #Room 2
    dS2 <- (w_Hog*R2) - (beta_Hog*I2*S2) - (beta_Hog_ind*(I3+I4)*S2) - (beta_WF2S*HI*S2) - (u_Hog*S2)
    dE2 <-(beta_Hog*I2*S2) + (beta_Hog_ind*(I3+I4)*S2) + (beta_WF2S*HI*S2) - (sigma_Hog*E2) - (u_Hog*E2) 
    dI2 <- (sigma_Hog*E2) - (delta_Hog*I2) - (u_inf_Hog*I2)
    dR2 <- (delta_Hog*I2) - (w_Hog*R2) - (u_Hog*R2)
    dTotal_inf2 <- (sigma_Hog*E2)
    dTotal_Rec2 <- (delta_Hog*I2)
    dD2_Nat <- (u_Hog*S2) + (u_Hog*E2) + (u_Hog*R2)
    dD2_Inf <- (u_inf_Hog*I2)
    dWF2Hog_Count2 <- (beta_WF2S*HI*S2)
    
    #Room 3
    dS3 <- (w_Hog*R3) - (beta_Hog*I3*S3) - (beta_Hog_ind*(I4)*S3) - (beta_WF2S*HI*S3) - (u_Hog*S3)
    dE3 <-(beta_Hog*I3*S3) + (beta_Hog_ind*(I4)*S3) + (beta_WF2S*HI*S3) - (sigma_Hog*E3) - (u_Hog*E3)
    dI3 <- (sigma_Hog*E3) - (delta_Hog*I3) - (u_inf_Hog*I3)
    dR3 <- (delta_Hog*I3) - (w_Hog*R3) - (u_Hog*R3)
    dTotal_inf3 <- (sigma_Hog*E3)
    dTotal_Rec3 <- (delta_Hog*I3)
    dD3_Nat <- (u_Hog*S3) + (u_Hog*E3) + (u_Hog*R3)
    dD3_Inf <- (u_inf_Hog*I3)
    dWF2Hog_Count3 <- (beta_WF2S*HI*S3)
    
    #Room 4
    dS4 <- w_Hog*R4 - (beta_Hog*I4*S4) - (beta_Hog_ind*(0)*S4) - (beta_WF2S*HI*S4) - (u_Hog*S4)
    dE4 <-(beta_Hog*I4*S4) + (beta_Hog_ind*(0)*S4) + (beta_WF2S*HI*S4) - (sigma_Hog*E4) - (u_Hog*E4) 
    dI4 <- (sigma_Hog*E4) - (delta_Hog*I4) - (u_inf_Hog*I4)
    dR4 <- (delta_Hog*I4) - (w_Hog*R4) - (u_Hog*R4)
    dTotal_inf4 <- (sigma_Hog*E4)
    dTotal_Rec4 <- (delta_Hog*I4)
    dD4_Nat <- (u_Hog*S4) + (u_Hog*E4) + (u_Hog*R4)
    dD4_Inf <- (u_inf_Hog*I4)
    dWF2Hog_Count4 <- (beta_WF2S*HI*S4)
    
    #Humans
    dHS <- (w_WF*HR)- (beta_WF2WF*HI*HS) - (beta_S2WF*(I1+I2+I3+I4)*HS) - (u_WF*HS)
    dHE <- (beta_WF2WF*HI*HS) + (beta_S2WF*(I1+I2+I3+I4)*HS) - (sigma_WF*HE) - (u_WF*HE)
    dHI <- (sigma_WF*HE) - (delta_WF*HI) - (u_inf_WF*HI)
    dHR <- (delta_WF*HI) - (w_WF*HR) - u_WF*HR
    dH_DN <- (u_WF*HS) + (u_WF*HE) + (u_WF*HR)
    dH_DI <- (u_inf_WF*HI)
    dHog2WF_Count <- (beta_S2WF*(I1+I2+I3+I4)*HS)
    
    
    # Last instruction: return a list 
    
    return(list(c(dS1, dE1, dI1, dR1, dTotal_inf1, dTotal_Rec1, dD1_Nat, dD1_Inf, dWF2Hog_Count1,
                  dS2, dE2, dI2, dR2, dTotal_inf2, dTotal_Rec2, dD2_Nat, dD2_Inf, dWF2Hog_Count2,
                  dS3, dE3, dI3, dR3, dTotal_inf3, dTotal_Rec3, dD3_Nat, dD3_Inf, dWF2Hog_Count3,
                  dS4, dE4, dI4, dR4, dTotal_inf4, dTotal_Rec4, dD4_Nat, dD4_Inf, dWF2Hog_Count4,
                  dHS, dHE, dHI, dHR, dH_DN, dH_DI, dHog2WF_Count)))
    
  }
  
  #### MODEL ITERATIONS ####
  #### WEEK 0 TO 2 (Room 1 fill) ####
  Week0to2 <- 14
  #====================================================================================================
  # Model Initial Values  Weeks 0 to 2
  #===========================================================================
  S1_int <- (1-Hog_Vacc_Eff)*1000
  E1_int <- 0
  I1_int <- 0 #1
  R1_int <- Hog_Vacc_Eff*1000
  Total_inf1_int <- 0
  Total_Rec1_int <- 0
  D1_Nat_int <- 0
  D1_Inf_int <- 0
  WF2Hog_Count1_int <- 0
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- 0
  Total_Rec2_int <- 0
  D2_Nat_int <- 0
  D2_Inf_int <- 0
  WF2Hog_Count2_int <- 0
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <- 2
  HE_int <- 0
  HI_int <- 0
  HR_int <- 0
  H_DN_int <- 0
  H_DI_int <- 0
  Hog2WF_Count_int <- 0
  
  #====================================================================================================
  # Model Run  Weeks 0 to 2
  #====================================================================================================
  Week0to2.t <- seq(0,Week0to2,0.1)   
  Week0to2.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week0to2.sol <- lsoda(Week0to2.init, Week0to2.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 2 TO 4 (Room 2 Fill)####
  Week2to4 <- 28
  #====================================================================================================
  # Model Initial Values  Weeks 2 to 4
  #===========================================================================
  S1_int <- Week0to2.sol[nrow(Week0to2.sol),2]
  E1_int <- Week0to2.sol[nrow(Week0to2.sol),3]
  I1_int <- Week0to2.sol[nrow(Week0to2.sol),4]
  R1_int <- Week0to2.sol[nrow(Week0to2.sol),5]
  Total_inf1_int <-  Week0to2.sol[nrow(Week0to2.sol),6]
  Total_Rec1_int <-  Week0to2.sol[nrow(Week0to2.sol),7]
  D1_Nat_int <-  Week0to2.sol[nrow(Week0to2.sol),8]
  D1_Inf_int <-  Week0to2.sol[nrow(Week0to2.sol),9]
  WF2Hog_Count1_int <-  Week0to2.sol[nrow(Week0to2.sol),10]
  
  S2_int <- (1-Hog_Vacc_Eff)*999
  E2_int <- 0
  I2_int <- 1
  R2_int <- Hog_Vacc_Eff*999
  Total_inf2_int <- 0
  Total_Rec2_int <- 0
  D2_Nat_int <- 0
  D2_Inf_int <- 0
  WF2Hog_Count2_int <- 0
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week0to2.sol[nrow(Week0to2.sol),38]
  HE_int <-  Week0to2.sol[nrow(Week0to2.sol),39]
  HI_int <-  Week0to2.sol[nrow(Week0to2.sol),40]
  HR_int <-  Week0to2.sol[nrow(Week0to2.sol),41]
  H_DN_int <-  Week0to2.sol[nrow(Week0to2.sol),42]
  H_DI_int <-  Week0to2.sol[nrow(Week0to2.sol),43]
  Hog2WF_Count_int <-  Week0to2.sol[nrow(Week0to2.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 2 to 4
  #====================================================================================================
  Week2to4.t <- seq(14.1,Week2to4,0.1)          
  Week2to4.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week2to4.sol <- lsoda(Week2to4.init, Week2to4.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 4 TO 6 (Room 3 Fill)####
  Week4to6 <- 42
  #====================================================================================================
  # Model Initial Values  Weeks 4 to 6
  #===========================================================================
  S1_int <- Week2to4.sol[nrow(Week2to4.sol),2]
  E1_int <- Week2to4.sol[nrow(Week2to4.sol),3]
  I1_int <- Week2to4.sol[nrow(Week2to4.sol),4]
  R1_int <- Week2to4.sol[nrow(Week2to4.sol),5]
  Total_inf1_int <- Week2to4.sol[nrow(Week2to4.sol),6]
  Total_Rec1_int <- Week2to4.sol[nrow(Week2to4.sol),7]
  D1_Nat_int <-  Week2to4.sol[nrow(Week2to4.sol),8]
  D1_Inf_int <-  Week2to4.sol[nrow(Week2to4.sol),9]
  WF2Hog_Count1_int <-  Week2to4.sol[nrow(Week2to4.sol),10]
  
  S2_int <- Week2to4.sol[nrow(Week2to4.sol),11]
  E2_int <- Week2to4.sol[nrow(Week2to4.sol),12]
  I2_int <- Week2to4.sol[nrow(Week2to4.sol),13]
  R2_int <- Week2to4.sol[nrow(Week2to4.sol),14]
  Total_inf2_int <- Week2to4.sol[nrow(Week2to4.sol),15]
  Total_Rec2_int <- Week2to4.sol[nrow(Week2to4.sol),16]
  D2_Nat_int <- Week2to4.sol[nrow(Week2to4.sol),17]
  D2_Inf_int <- Week2to4.sol[nrow(Week2to4.sol),18]
  WF2Hog_Count2_int <- Week2to4.sol[nrow(Week2to4.sol),19]
  
  S3_int <- (1-Hog_Vacc_Eff)*1000
  E3_int <- 0
  I3_int <- 0
  R3_int <-  Hog_Vacc_Eff*1000
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week2to4.sol[nrow(Week2to4.sol),38]
  HE_int <-  Week2to4.sol[nrow(Week2to4.sol),39]
  HI_int <-  Week2to4.sol[nrow(Week2to4.sol),40]
  HR_int <-  Week2to4.sol[nrow(Week2to4.sol),41]
  H_DN_int <-  Week2to4.sol[nrow(Week2to4.sol),42]
  H_DI_int <-  Week2to4.sol[nrow(Week2to4.sol),43]
  Hog2WF_Count_int <-  Week2to4.sol[nrow(Week2to4.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 4 to 6
  #====================================================================================================
  Week4to6.t <- seq(28.1,Week4to6,0.1)         
  Week4to6.init <-c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                    S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                    S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                    S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                    HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week4to6.sol <- lsoda(Week4to6.init, Week4to6.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 6 TO 8 (Room 4 Fill)####
  Week6to8 <- 56
  #====================================================================================================
  # Model Initial Values  Weeks 6 to 8
  #===========================================================================
  S1_int <- Week4to6.sol[nrow(Week4to6.sol),2]
  E1_int <- Week4to6.sol[nrow(Week4to6.sol),3]
  I1_int <- Week4to6.sol[nrow(Week4to6.sol),4]
  R1_int <- Week4to6.sol[nrow(Week4to6.sol),5]
  Total_inf1_int <- Week4to6.sol[nrow(Week4to6.sol),6]
  Total_Rec1_int <- Week4to6.sol[nrow(Week4to6.sol),7]
  D1_Nat_int <-  Week4to6.sol[nrow(Week4to6.sol),8]
  D1_Inf_int <-  Week4to6.sol[nrow(Week4to6.sol),9]
  WF2Hog_Count1_int <-  Week4to6.sol[nrow(Week4to6.sol),10]
  
  S2_int <- Week4to6.sol[nrow(Week4to6.sol),11]
  E2_int <- Week4to6.sol[nrow(Week4to6.sol),12]
  I2_int <- Week4to6.sol[nrow(Week4to6.sol),13]
  R2_int <- Week4to6.sol[nrow(Week4to6.sol),14]
  Total_inf2_int <- Week4to6.sol[nrow(Week4to6.sol),15]
  Total_Rec2_int <- Week4to6.sol[nrow(Week4to6.sol),16]
  D2_Nat_int <- Week4to6.sol[nrow(Week4to6.sol),17]
  D2_Inf_int <- Week4to6.sol[nrow(Week4to6.sol),18]
  WF2Hog_Count2_int <- Week4to6.sol[nrow(Week4to6.sol),19]
  
  S3_int <- Week4to6.sol[nrow(Week4to6.sol),20]
  E3_int <- Week4to6.sol[nrow(Week4to6.sol),21]
  I3_int <- Week4to6.sol[nrow(Week4to6.sol),22]
  R3_int <-  Week4to6.sol[nrow(Week4to6.sol),23]
  Total_inf3_int <- Week4to6.sol[nrow(Week4to6.sol),24]
  Total_Rec3_int <- Week4to6.sol[nrow(Week4to6.sol),25]
  D3_Nat_int <- Week4to6.sol[nrow(Week4to6.sol),26]
  D3_Inf_int <- Week4to6.sol[nrow(Week4to6.sol),27]
  WF2Hog_Count3_int <- Week4to6.sol[nrow(Week4to6.sol),28]
  
  S4_int <- (1-Hog_Vacc_Eff)*1000
  E4_int <- 0
  I4_int <- 0
  R4_int <- Hog_Vacc_Eff*1000
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week4to6.sol[nrow(Week4to6.sol),38]
  HE_int <-  Week4to6.sol[nrow(Week4to6.sol),39]
  HI_int <-  Week4to6.sol[nrow(Week4to6.sol),40]
  HR_int <-  Week4to6.sol[nrow(Week4to6.sol),41]
  H_DN_int <-  Week4to6.sol[nrow(Week4to6.sol),42]
  H_DI_int <-  Week4to6.sol[nrow(Week4to6.sol),43]
  Hog2WF_Count_int <-  Week4to6.sol[nrow(Week4to6.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 6 to 8
  #====================================================================================================
  Week6to8.t <- seq(42.1,Week6to8,0.1)           
  Week6to8.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week6to8.sol <- lsoda(Week6to8.init, Week6to8.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 8 TO 23 (Room 1 Exit) ####
  Week8to23 <- 161
  #====================================================================================================
  # Model Initial Values  Weeks 8 to 23
  #===========================================================================
  S1_int <- Week6to8.sol[nrow(Week6to8.sol),2]
  E1_int <- Week6to8.sol[nrow(Week6to8.sol),3]
  I1_int <- Week6to8.sol[nrow(Week6to8.sol),4]
  R1_int <- Week6to8.sol[nrow(Week6to8.sol),5]
  Total_inf1_int <- Week6to8.sol[nrow(Week6to8.sol),6]
  Total_Rec1_int <- Week6to8.sol[nrow(Week6to8.sol),7]
  D1_Nat_int <-  Week6to8.sol[nrow(Week6to8.sol),8]
  D1_Inf_int <-  Week6to8.sol[nrow(Week6to8.sol),9]
  WF2Hog_Count1_int <-  Week6to8.sol[nrow(Week6to8.sol),10]
  
  
  S2_int <- Week6to8.sol[nrow(Week6to8.sol),11]
  E2_int <- Week6to8.sol[nrow(Week6to8.sol),12]
  I2_int <- Week6to8.sol[nrow(Week6to8.sol),13]
  R2_int <- Week6to8.sol[nrow(Week6to8.sol),14]
  Total_inf2_int <- Week6to8.sol[nrow(Week6to8.sol),15]
  Total_Rec2_int <- Week6to8.sol[nrow(Week6to8.sol),16]
  D2_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),17]
  D2_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),18]
  WF2Hog_Count2_int <- Week6to8.sol[nrow(Week6to8.sol),19]
  
  S3_int <- Week6to8.sol[nrow(Week6to8.sol),20]
  E3_int <- Week6to8.sol[nrow(Week6to8.sol),21]
  I3_int <- Week6to8.sol[nrow(Week6to8.sol),22]
  R3_int <-  Week6to8.sol[nrow(Week6to8.sol),23]
  Total_inf3_int <- Week6to8.sol[nrow(Week6to8.sol),24]
  Total_Rec3_int <- Week6to8.sol[nrow(Week6to8.sol),25]
  D3_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),26]
  D3_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),27]
  WF2Hog_Count3_int <- Week6to8.sol[nrow(Week6to8.sol),28]
  
  S4_int <- Week6to8.sol[nrow(Week6to8.sol),29]
  E4_int <- Week6to8.sol[nrow(Week6to8.sol),30]
  I4_int <- Week6to8.sol[nrow(Week6to8.sol),31]
  R4_int <- Week6to8.sol[nrow(Week6to8.sol),32]
  Total_inf4_int <- Week6to8.sol[nrow(Week6to8.sol),33]
  Total_Rec4_int <- Week6to8.sol[nrow(Week6to8.sol),34]
  D4_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),35]
  D4_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),36]
  WF2Hog_Count4_int <- Week6to8.sol[nrow(Week6to8.sol),37]
  
  HS_int <-  Week6to8.sol[nrow(Week6to8.sol),38]
  HE_int <-  Week6to8.sol[nrow(Week6to8.sol),39]
  HI_int <-  Week6to8.sol[nrow(Week6to8.sol),40]
  HR_int <-  Week6to8.sol[nrow(Week6to8.sol),41]
  H_DN_int <-  Week6to8.sol[nrow(Week6to8.sol),42]
  H_DI_int <-  Week6to8.sol[nrow(Week6to8.sol),43]
  Hog2WF_Count_int <-  Week6to8.sol[nrow(Week6to8.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 8 to 23
  #====================================================================================================
  Week8to23.t <- seq(56.1,Week8to23,0.1)          
  Week8to23.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                      S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                      S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                      S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                      HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week8to23.sol <- lsoda(Week8to23.init, Week8to23.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 23 TO 25 (Room 2 Exit) ####
  Week23to25 <- 175
  #====================================================================================================
  # Model Initial Values  Weeks 23 to 25
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week8to23.sol[nrow(Week8to23.sol),6]
  Total_Rec1_int <- Week8to23.sol[nrow(Week8to23.sol),7]
  D1_Nat_int <-  Week8to23.sol[nrow(Week8to23.sol),8]
  D1_Inf_int <-  Week8to23.sol[nrow(Week8to23.sol),9]
  WF2Hog_Count1_int <-  Week8to23.sol[nrow(Week8to23.sol),10]
  
  
  S2_int <- Week8to23.sol[nrow(Week8to23.sol),11]
  E2_int <- Week8to23.sol[nrow(Week8to23.sol),12]
  I2_int <- Week8to23.sol[nrow(Week8to23.sol),13]
  R2_int <- Week8to23.sol[nrow(Week8to23.sol),14]
  Total_inf2_int <- Week8to23.sol[nrow(Week8to23.sol),15]
  Total_Rec2_int <- Week8to23.sol[nrow(Week8to23.sol),16]
  D2_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),17]
  D2_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),18]
  WF2Hog_Count2_int <- Week8to23.sol[nrow(Week8to23.sol),19]
  
  S3_int <- Week8to23.sol[nrow(Week8to23.sol),20]
  E3_int <- Week8to23.sol[nrow(Week8to23.sol),21]
  I3_int <- Week8to23.sol[nrow(Week8to23.sol),22]
  R3_int <-  Week8to23.sol[nrow(Week8to23.sol),23]
  Total_inf3_int <- Week8to23.sol[nrow(Week8to23.sol),24]
  Total_Rec3_int <- Week8to23.sol[nrow(Week8to23.sol),25]
  D3_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),26]
  D3_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),27]
  WF2Hog_Count3_int <- Week8to23.sol[nrow(Week8to23.sol),28]
  
  S4_int <- Week8to23.sol[nrow(Week8to23.sol),29]
  E4_int <- Week8to23.sol[nrow(Week8to23.sol),30]
  I4_int <- Week8to23.sol[nrow(Week8to23.sol),31]
  R4_int <- Week8to23.sol[nrow(Week8to23.sol),32]
  Total_inf4_int <- Week8to23.sol[nrow(Week8to23.sol),33]
  Total_Rec4_int <- Week8to23.sol[nrow(Week8to23.sol),34]
  D4_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),35]
  D4_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),36]
  WF2Hog_Count4_int <- Week8to23.sol[nrow(Week8to23.sol),37]
  
  HS_int <-  Week8to23.sol[nrow(Week8to23.sol),38]
  HE_int <-  Week8to23.sol[nrow(Week8to23.sol),39]
  HI_int <-  Week8to23.sol[nrow(Week8to23.sol),40]
  HR_int <-  Week8to23.sol[nrow(Week8to23.sol),41]
  H_DN_int <-  Week8to23.sol[nrow(Week8to23.sol),42]
  H_DI_int <-  Week8to23.sol[nrow(Week8to23.sol),43]
  Hog2WF_Count_int <-  Week8to23.sol[nrow(Week8to23.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 23 to 25
  #====================================================================================================
  Week23to25.t <- seq(161.1,Week23to25,0.1)         
  Week23to25.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week23to25.sol <- lsoda(Week23to25.init, Week23to25.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 25 TO 27 (Room 3 Exit) ####
  Week25to27 <- 189
  #====================================================================================================
  # Model Initial Values  Weeks 25 to 27
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week23to25.sol[nrow(Week23to25.sol),6]
  Total_Rec1_int <- Week23to25.sol[nrow(Week23to25.sol),7]
  D1_Nat_int <-  Week23to25.sol[nrow(Week23to25.sol),8]
  D1_Inf_int <-  Week23to25.sol[nrow(Week23to25.sol),9]
  WF2Hog_Count1_int <-  Week23to25.sol[nrow(Week23to25.sol),10]
  
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- Week23to25.sol[nrow(Week23to25.sol),15]
  Total_Rec2_int <- Week23to25.sol[nrow(Week23to25.sol),16]
  D2_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),17]
  D2_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),18]
  WF2Hog_Count2_int <- Week23to25.sol[nrow(Week23to25.sol),19]
  
  S3_int <- Week23to25.sol[nrow(Week23to25.sol),20]
  E3_int <- Week23to25.sol[nrow(Week23to25.sol),21]
  I3_int <- Week23to25.sol[nrow(Week23to25.sol),22]
  R3_int <-  Week23to25.sol[nrow(Week23to25.sol),23]
  Total_inf3_int <- Week23to25.sol[nrow(Week23to25.sol),24]
  Total_Rec3_int <- Week23to25.sol[nrow(Week23to25.sol),25]
  D3_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),26]
  D3_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),27]
  WF2Hog_Count3_int <- Week23to25.sol[nrow(Week23to25.sol),28]
  
  S4_int <- Week23to25.sol[nrow(Week23to25.sol),29]
  E4_int <- Week23to25.sol[nrow(Week23to25.sol),30]
  I4_int <- Week23to25.sol[nrow(Week23to25.sol),31]
  R4_int <- Week23to25.sol[nrow(Week23to25.sol),32]
  Total_inf4_int <- Week23to25.sol[nrow(Week23to25.sol),33]
  Total_Rec4_int <- Week23to25.sol[nrow(Week23to25.sol),34]
  D4_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),35]
  D4_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),36]
  WF2Hog_Count4_int <- Week23to25.sol[nrow(Week23to25.sol),37]
  
  HS_int <-  Week23to25.sol[nrow(Week23to25.sol),38]
  HE_int <-  Week23to25.sol[nrow(Week23to25.sol),39]
  HI_int <-  Week23to25.sol[nrow(Week23to25.sol),40]
  HR_int <-  Week23to25.sol[nrow(Week23to25.sol),41]
  H_DN_int <-  Week23to25.sol[nrow(Week23to25.sol),42]
  H_DI_int <-  Week23to25.sol[nrow(Week23to25.sol),43]
  Hog2WF_Count_int <-  Week23to25.sol[nrow(Week23to25.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 25 to 27
  #====================================================================================================
  Week25to27.t <- seq(175.1,Week25to27,0.1)          
  Week25to27.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week25to27.sol <- lsoda(Week25to27.init, Week25to27.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 27 TO 29 (Room 4 Exit) ####
  Week27to29 <- 203
  #====================================================================================================
  # Model Initial Values  Weeks 27to 29
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week25to27.sol[nrow(Week25to27.sol),6]
  Total_Rec1_int <- Week25to27.sol[nrow(Week25to27.sol),7]
  D1_Nat_int <-  Week25to27.sol[nrow(Week25to27.sol),8]
  D1_Inf_int <-  Week25to27.sol[nrow(Week25to27.sol),9]
  WF2Hog_Count1_int <-  Week25to27.sol[nrow(Week25to27.sol),10]
  
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- Week25to27.sol[nrow(Week25to27.sol),15]
  Total_Rec2_int <- Week25to27.sol[nrow(Week25to27.sol),16]
  D2_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),17]
  D2_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),18]
  WF2Hog_Count2_int <- Week25to27.sol[nrow(Week25to27.sol),19]
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- Week25to27.sol[nrow(Week25to27.sol),24]
  Total_Rec3_int <- Week25to27.sol[nrow(Week25to27.sol),25]
  D3_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),26]
  D3_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),27]
  WF2Hog_Count3_int <- Week25to27.sol[nrow(Week25to27.sol),28]
  
  S4_int <- Week25to27.sol[nrow(Week25to27.sol),29]
  E4_int <- Week25to27.sol[nrow(Week25to27.sol),30]
  I4_int <- Week25to27.sol[nrow(Week25to27.sol),31]
  R4_int <- Week25to27.sol[nrow(Week25to27.sol),32]
  Total_inf4_int <- Week25to27.sol[nrow(Week25to27.sol),33]
  Total_Rec4_int <- Week25to27.sol[nrow(Week25to27.sol),34]
  D4_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),35]
  D4_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),36]
  WF2Hog_Count4_int <- Week25to27.sol[nrow(Week25to27.sol),37]
  
  HS_int <-  Week25to27.sol[nrow(Week25to27.sol),38]
  HE_int <-  Week25to27.sol[nrow(Week25to27.sol),39]
  HI_int <-  Week25to27.sol[nrow(Week25to27.sol),40]
  HR_int <-  Week25to27.sol[nrow(Week25to27.sol),41]
  H_DN_int <-  Week25to27.sol[nrow(Week25to27.sol),42]
  H_DI_int <-  Week25to27.sol[nrow(Week25to27.sol),43]
  Hog2WF_Count_int <-  Week25to27.sol[nrow(Week25to27.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 27 to 29
  #====================================================================================================
  Week27to29.t <- seq(189.1,Week27to29,0.1)            
  Week27to29.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week27to29.sol <- lsoda(Week27to29.init, Week27to29.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### COMBINING THE DATASETS ####
  Week0to2.sol.df <- as.data.frame(Week0to2.sol)
  colnames(Week0to2.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week2to4.sol.df <- as.data.frame(Week2to4.sol)
  colnames(Week2to4.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week4to6.sol.df <- as.data.frame(Week4to6.sol)
  colnames(Week4to6.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week6to8.sol.df <- as.data.frame(Week6to8.sol)
  colnames(Week6to8.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week8to23.sol.df <- as.data.frame(Week8to23.sol)
  colnames(Week8to23.sol.df) <- c("time",
                                  "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                  "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                  "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                  "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                  "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week23to25.sol.df <- as.data.frame(Week23to25.sol)
  colnames(Week23to25.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week25to27.sol.df <- as.data.frame(Week25to27.sol)
  colnames(Week25to27.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week27to29.sol.df <- as.data.frame(Week27to29.sol)
  colnames(Week27to29.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  
  fulldata.df <- rbind(Week0to2.sol.df, Week2to4.sol.df, Week4to6.sol.df, Week6to8.sol.df, 
                       Week8to23.sol.df, Week23to25.sol.df, Week25to27.sol.df, Week27to29.sol.df)
  
  ## Is there a WF infection
  HI_inf <- fulldata.df %>%
    summarise(total_I = sum(HI)) %>%
    mutate(HI_Inf = ifelse(total_I>1, 1, 0)) %>%
    select(HI_Inf) %>% as.numeric()
  #Time to 1st WF infection
  Time2WFinf0 <- fulldata.df %>%
    arrange(time) %>%
    filter(HI > 0) %>%
    filter(row_number() == 1) %>%
    mutate(TimeSineInfIntro = time - 14) %>%
    select(TimeSineInfIntro) %>% as.numeric()
  
  Time2WFinf.5 <- fulldata.df %>%
    arrange(time) %>%
    filter(HI > 0.5) %>%
    filter(row_number() == 1) %>%
    mutate(TimeSineInfIntro = time - 14) %>%
    select(TimeSineInfIntro) %>% as.numeric()
  ## Time to Second Room Infection
  Timeto2room <- fulldata.df %>% 
    summarise(t_cross = time[which(E1 > 1 | 
                                     E3 > 1 |
                                     E4 > 1)]) %>% 
    summarise(min_t_cross = min(t_cross)) %>% as.numeric()
  ## Total Hogs infected
  TotHogsInf <- fulldata.df %>% 
    #select(t, iter, Total_inf1, Total_inf2, Total_inf3, Total_inf4) %>% 
    mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
           Tot_rec = Total_Rec1 + Total_Rec2 + Total_Rec3 + Total_Rec4) %>% 
    summarise(Tot_inf2 = max(Tot_inf),
              Tot_rec2 = max(Tot_rec))
  ## Minimum Time to Peak Infection
  Time2PeakInf <- fulldata.df %>% 
    mutate(S_tot = S1+S2+S3+S4,                             
           E_tot = E1+E2+E3+E4,
           I_tot = I1+I2+I3+I4,
           R_tot = R1+R2+R3+R4) %>%                                     
    summarise(max_I = max(I_tot),
              t_peak = time[which(I_tot == max(I_tot))]) %>% 
    summarise(min_t_peak = min(t_peak)) %>% as.numeric()
  ##Ro
  Ro <- fulldata.df %>% 
    arrange(time) %>% 
    mutate(Tot_rec = Total_Rec1 + Total_Rec2 + Total_Rec3 + Total_Rec4) %>%  
    filter(row_number() == n()) %>% 
    select(Tot_rec) %>% 
    mutate(Ro = ((log((4000 - Tot_rec)/4000) - log((4000 - (Hog_Vacc_Eff*4000))/4000)) /
                   ((((4000 - Tot_rec) / 4000) - (4000 - (Hog_Vacc_Eff*4000)) / 4000)))) %>% 
    select(Ro) %>% as.numeric()
  
  #Create Summary data frame
  fulldata.df.sum <- data.frame(iter = as.factor(i),
                                HI_inf = HI_inf,
                                Time2WFinf0 = Time2WFinf0,
                                Time2WFinf.5 = Time2WFinf.5,
                                Timeto2room = Timeto2room,
                                TotHogsInf = TotHogsInf$Tot_inf2,
                                TotHogsRec = TotHogsInf$Tot_rec2,
                                Time2PeakInf = Time2PeakInf,
                                Ro = Ro,
                                beta_Hog = beta_Hog,
                                beta_Hog_ind = beta_Hog_ind,
                                beta_S2WF = beta_S2WF,
                                beta_WF2S = beta_WF2S,
                                sigma_Hog = sigma_Hog,
                                delta_Hog = delta_Hog,
                                beta_WF2WF = beta_WF2WF,
                                sigma_WF = sigma_WF,
                                delta_WF = delta_WF)
  
  results <- rbind(results, fulldata.df.sum)
  
  rm(list=ls()[! ls() %in% c("results", "n.itr")])  #Deletes everything from environment expect growing results dataframe
}

saveRDS(results, file = "Current code/Control Meausres/Sensativity Analysis/WFFlowRoom2_Sum_SensAnlaysis.rds")

#### WF Flow Only (Room 1 Intro) ####
results <- NULL
n.itr <- 5000
# seed.num <- 616
seed.nums <- c(100:(100+n.itr))
for (i in c(1:n.itr)) {
  print(i)
  
  #### MODEL CONSTANTS ####
  #Hog constants
  beta_Hog <- rtri(n = 1, min = 0.001, max = 0.1, mode = (10/(1000*5)))
  beta_Hog_ind <- beta_Hog / 178 #(following Etbaigha et al 2018 paper)
  # beta_Hog_ind <- beta_Hog / 500
  Hog_Vacc_Eff <- 0
  u_Hog <- 0.00028     #natural death rate (following Etbaigha et al 2018 paper)
  u_inf_Hog <- 0 #0.1     #infected death rate (0 for now to problem solve)
  w_Hog <- 0 #1/180         #(following Etbaigha et al 2018 paper)  White et al shows a range from 56 - 112 days
  sigma_Hog <- abs(1 / rnorm(n = 1, mean = 2, sd = 1))   #Latent period h1n1 and h3n2
  delta_Hog <- abs(1 / rnorm(n = 1, mean = 5, sd = 1))   #Infectious period h1n1 and h3n2
  
  #Interspecies transmission
  ##applied the beta = Ro / N*D from the modeling book
  # Ro from the following website for swine to human infection 
  # https://bmcmedicine.biomedcentral.com/articles/10.1186/1741-7015-7-30#Sec6   (from super-strain figure) possible the high end of transmisilibty for this parameter
  # beta_S2WF <- (2.3 / (4002*5))  #  #Hogs + #Workforce * Duration of Hogs
  #Alt beta_S2WF from fair outbreak
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3557885/
  #assuming one hog got the first three confirmed cases sick and infectious period of 5 days Ro = 3/5 = 0.6
  beta_S2WF <- runif(1, min = 0.000001, max = 0.0001)
  #when Ro WF2S = Ro S2WF
  # beta_WF2S <- (2.3 / (4002*3))  #  #Hogs + #Workforce * Duration of Hogs
  # WF2S calculated from Swine Outbreak of Pandemic Influenza A Virus on a Canadian Research Farm Supports Human-to-Swine Transmission paper
  beta_WF2S <- runif(1, min = 0.000001, max = 0.0001)    #Hogs + #Workforce * Duration of WF
  #Human constants
  #applied the beta = Ro / N*D from the modeling book
  beta_WF2WF <-   rtri(n = 1, min = 0.0000001, max = 0.64, mode = 0.32) # http://m-hikari.com/ams/ams-2013/ams-41-44-2013/joseAMS41-44-2013.pdf 
  Hum_Vacc_Eff <- 0.0 
  u_WF <- 0     #Natural death rate of human (0 for now to problem solve)
  u_inf_WF <- 0   #Infected death rate of human (0 for now to problem solve)
  w_WF <- 0      #Recovery rate for human (0 for now to problem solve)
  sigma_WF <- abs(1 / rnorm(n = 1, mean = 2, sd = 1))   #Latent period for humans
  delta_WF <- abs(1 / rnorm(n = 1, mean = 3, sd = 1))   #Infectious Period for humans
  
  #### MODEL SETUP ####
  # Model Parameters
  TransmissionModel.par <- c(beta_Hog = beta_Hog, beta_Hog_ind = beta_Hog_ind, beta_S2WF = beta_S2WF, beta_WF2S = beta_WF2S,
                             w_Hog = w_Hog, sigma_Hog = sigma_Hog, delta_Hog = delta_Hog, u_Hog = u_Hog, u_inf_Hog = u_inf_Hog,
                             beta_WF2WF = beta_WF2WF,
                             w_WF = w_WF, sigma_WF = sigma_WF, delta_WF = delta_WF, u_WF = u_WF, u_inf_WF = u_inf_WF)
  #=====================================================================================================
  #Model Equations 
  #=====================================================================================================
  TransmissionModel.dyn <- function(t,var,par) {
    # Rename the parameters 
    beta_Hog <- par[1]
    beta_Hog_ind <- par[2]
    beta_S2WF <- par[3]
    beta_WF2S <- par[4]
    w_Hog <- par[5]
    sigma_Hog <- par[6]
    delta_Hog <- par[7]
    u_Hog <- par[8]
    u_inf_Hog<- par[9]
    beta_WF2WF <- par[10]
    w_WF <- par[11]
    sigma_WF <- par[12]
    delta_WF <- par[13]
    u_WF <- par[14]
    u_inf_WF <- par[15]
    
    S1 <- var[1]
    E1 <- var[2]
    I1 <- var[3]
    R1 <- var[4]
    Total_inf1 <- var[5] 
    True_Rec1 <- var[6]
    D1_Nat <- var[7]
    D1_Inf <- var[8]
    WF2Hog_Count1 <- var[9] 
    
    S2 <- var[10]
    E2 <- var[11]
    I2 <- var[12]
    R2 <- var[13]
    Total_inf2 <- var[14]
    True_Rec2 <- var[15]
    D2_Nat <- var[16]
    D2_Inf <- var[17]
    WF2Hog_Count2 <- var[18]
    
    S3 <- var[19]
    E3 <- var[20]
    I3 <- var[21]
    R3 <- var[22]
    Total_inf3 <- var[23]
    True_Rec3 <- var[24]
    D3_Nat <- var[25]
    D3_Inf <- var[26]
    WF2Hog_Count3 <- var[27]   
    
    S4 <- var[28]
    E4 <- var[29]
    I4 <- var[30]
    R4 <- var[31]
    Total_inf4 <- var[32]
    True_Rec4 <- var[33]
    D4_Nat <- var[34]
    D4_Inf <- var[35]
    WF2Hog_Count4 <- var[36]
    
    HS <- var[37]
    HE <- var[38]
    HI <- var[39] 
    HR <- var[40]
    H_DN <- var[41]
    H_DI <- var[42]
    Hog2WF_Count <- var[43] 
    
    # Calculate the derivatives
    #Room 1
    dS1 <- (w_Hog*R1) - (beta_Hog*I1*S1) - (beta_Hog_ind*(I2+I3+I4)*S1) - (beta_WF2S*HI*S1) - (u_Hog*S1)
    dE1 <-(beta_Hog*I1*S1) + (beta_Hog_ind*(I2+I3+I4)*S1) + (beta_WF2S*HI*S1) - sigma_Hog*E1 - (u_Hog*E1)
    dI1 <- (sigma_Hog*E1) - (delta_Hog*I1) - (u_inf_Hog*I1)
    dR1 <- (delta_Hog*I1) - (w_Hog*R1) - (u_Hog*R1)
    dTotal_inf1 <- (sigma_Hog*E1)
    dTotal_Rec1 <- (delta_Hog*I1)
    dD1_Nat <- (u_Hog*S1) + (u_Hog*E1) + (u_Hog*R1)
    dD1_Inf <- (u_inf_Hog*I1)
    dWF2Hog_Count1 <- (beta_WF2S*HI*S1)
    
    #Room 2
    dS2 <- (w_Hog*R2) - (beta_Hog*I2*S2) - (beta_Hog_ind*(I3+I4)*S2) - (beta_WF2S*HI*S2) - (u_Hog*S2)
    dE2 <-(beta_Hog*I2*S2) + (beta_Hog_ind*(I3+I4)*S2) + (beta_WF2S*HI*S2) - (sigma_Hog*E2) - (u_Hog*E2) 
    dI2 <- (sigma_Hog*E2) - (delta_Hog*I2) - (u_inf_Hog*I2)
    dR2 <- (delta_Hog*I2) - (w_Hog*R2) - (u_Hog*R2)
    dTotal_inf2 <- (sigma_Hog*E2)
    dTotal_Rec2 <- (delta_Hog*I2)
    dD2_Nat <- (u_Hog*S2) + (u_Hog*E2) + (u_Hog*R2)
    dD2_Inf <- (u_inf_Hog*I2)
    dWF2Hog_Count2 <- (beta_WF2S*HI*S2)
    
    #Room 3
    dS3 <- (w_Hog*R3) - (beta_Hog*I3*S3) - (beta_Hog_ind*(I4)*S3) - (beta_WF2S*HI*S3) - (u_Hog*S3)
    dE3 <-(beta_Hog*I3*S3) + (beta_Hog_ind*(I4)*S3) + (beta_WF2S*HI*S3) - (sigma_Hog*E3) - (u_Hog*E3)
    dI3 <- (sigma_Hog*E3) - (delta_Hog*I3) - (u_inf_Hog*I3)
    dR3 <- (delta_Hog*I3) - (w_Hog*R3) - (u_Hog*R3)
    dTotal_inf3 <- (sigma_Hog*E3)
    dTotal_Rec3 <- (delta_Hog*I3)
    dD3_Nat <- (u_Hog*S3) + (u_Hog*E3) + (u_Hog*R3)
    dD3_Inf <- (u_inf_Hog*I3)
    dWF2Hog_Count3 <- (beta_WF2S*HI*S3)
    
    #Room 4
    dS4 <- w_Hog*R4 - (beta_Hog*I4*S4) - (beta_Hog_ind*(0)*S4) - (beta_WF2S*HI*S4) - (u_Hog*S4)
    dE4 <-(beta_Hog*I4*S4) + (beta_Hog_ind*(0)*S4) + (beta_WF2S*HI*S4) - (sigma_Hog*E4) - (u_Hog*E4) 
    dI4 <- (sigma_Hog*E4) - (delta_Hog*I4) - (u_inf_Hog*I4)
    dR4 <- (delta_Hog*I4) - (w_Hog*R4) - (u_Hog*R4)
    dTotal_inf4 <- (sigma_Hog*E4)
    dTotal_Rec4 <- (delta_Hog*I4)
    dD4_Nat <- (u_Hog*S4) + (u_Hog*E4) + (u_Hog*R4)
    dD4_Inf <- (u_inf_Hog*I4)
    dWF2Hog_Count4 <- (beta_WF2S*HI*S4)
    
    #Humans
    dHS <- (w_WF*HR)- (beta_WF2WF*HI*HS) - (beta_S2WF*(I1+I2+I3+I4)*HS) - (u_WF*HS)
    dHE <- (beta_WF2WF*HI*HS) + (beta_S2WF*(I1+I2+I3+I4)*HS) - (sigma_WF*HE) - (u_WF*HE)
    dHI <- (sigma_WF*HE) - (delta_WF*HI) - (u_inf_WF*HI)
    dHR <- (delta_WF*HI) - (w_WF*HR) - u_WF*HR
    dH_DN <- (u_WF*HS) + (u_WF*HE) + (u_WF*HR)
    dH_DI <- (u_inf_WF*HI)
    dHog2WF_Count <- (beta_S2WF*(I1+I2+I3+I4)*HS)
    
    
    # Last instruction: return a list 
    
    return(list(c(dS1, dE1, dI1, dR1, dTotal_inf1, dTotal_Rec1, dD1_Nat, dD1_Inf, dWF2Hog_Count1,
                  dS2, dE2, dI2, dR2, dTotal_inf2, dTotal_Rec2, dD2_Nat, dD2_Inf, dWF2Hog_Count2,
                  dS3, dE3, dI3, dR3, dTotal_inf3, dTotal_Rec3, dD3_Nat, dD3_Inf, dWF2Hog_Count3,
                  dS4, dE4, dI4, dR4, dTotal_inf4, dTotal_Rec4, dD4_Nat, dD4_Inf, dWF2Hog_Count4,
                  dHS, dHE, dHI, dHR, dH_DN, dH_DI, dHog2WF_Count)))
    
  }
  
  #### MODEL ITERATIONS ####
  #### WEEK 0 TO 2 (Room 1 fill) ####
  Week0to2 <- 14
  #====================================================================================================
  # Model Initial Values  Weeks 0 to 2
  #===========================================================================
  S1_int <- (1-Hog_Vacc_Eff)*999
  E1_int <- 0
  I1_int <- 1
  R1_int <- Hog_Vacc_Eff*999
  Total_inf1_int <- 0
  Total_Rec1_int <- 0
  D1_Nat_int <- 0
  D1_Inf_int <- 0
  WF2Hog_Count1_int <- 0
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- 0
  Total_Rec2_int <- 0
  D2_Nat_int <- 0
  D2_Inf_int <- 0
  WF2Hog_Count2_int <- 0
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <- 2
  HE_int <- 0
  HI_int <- 0
  HR_int <- 0
  H_DN_int <- 0
  H_DI_int <- 0
  Hog2WF_Count_int <- 0
  
  #====================================================================================================
  # Model Run  Weeks 0 to 2
  #====================================================================================================
  Week0to2.t <- seq(0,Week0to2,0.1)   
  Week0to2.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week0to2.sol <- lsoda(Week0to2.init, Week0to2.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 2 TO 4 (Room 2 Fill)####
  Week2to4 <- 28
  #====================================================================================================
  # Model Initial Values  Weeks 2 to 4
  #===========================================================================
  S1_int <- Week0to2.sol[nrow(Week0to2.sol),2]
  E1_int <- Week0to2.sol[nrow(Week0to2.sol),3]
  I1_int <- Week0to2.sol[nrow(Week0to2.sol),4]
  R1_int <- Week0to2.sol[nrow(Week0to2.sol),5]
  Total_inf1_int <-  Week0to2.sol[nrow(Week0to2.sol),6]
  Total_Rec1_int <-  Week0to2.sol[nrow(Week0to2.sol),7]
  D1_Nat_int <-  Week0to2.sol[nrow(Week0to2.sol),8]
  D1_Inf_int <-  Week0to2.sol[nrow(Week0to2.sol),9]
  WF2Hog_Count1_int <-  Week0to2.sol[nrow(Week0to2.sol),10]
  
  S2_int <- (1-Hog_Vacc_Eff)*1000
  E2_int <- 0
  I2_int <- 0
  R2_int <- Hog_Vacc_Eff*1000
  Total_inf2_int <- 0
  Total_Rec2_int <- 0
  D2_Nat_int <- 0
  D2_Inf_int <- 0
  WF2Hog_Count2_int <- 0
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week0to2.sol[nrow(Week0to2.sol),38]
  HE_int <-  Week0to2.sol[nrow(Week0to2.sol),39]
  HI_int <-  Week0to2.sol[nrow(Week0to2.sol),40]
  HR_int <-  Week0to2.sol[nrow(Week0to2.sol),41]
  H_DN_int <-  Week0to2.sol[nrow(Week0to2.sol),42]
  H_DI_int <-  Week0to2.sol[nrow(Week0to2.sol),43]
  Hog2WF_Count_int <-  Week0to2.sol[nrow(Week0to2.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 2 to 4
  #====================================================================================================
  Week2to4.t <- seq(14.1,Week2to4,0.1)          
  Week2to4.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week2to4.sol <- lsoda(Week2to4.init, Week2to4.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 4 TO 6 (Room 3 Fill)####
  Week4to6 <- 42
  #====================================================================================================
  # Model Initial Values  Weeks 4 to 6
  #===========================================================================
  S1_int <- Week2to4.sol[nrow(Week2to4.sol),2]
  E1_int <- Week2to4.sol[nrow(Week2to4.sol),3]
  I1_int <- Week2to4.sol[nrow(Week2to4.sol),4]
  R1_int <- Week2to4.sol[nrow(Week2to4.sol),5]
  Total_inf1_int <- Week2to4.sol[nrow(Week2to4.sol),6]
  Total_Rec1_int <- Week2to4.sol[nrow(Week2to4.sol),7]
  D1_Nat_int <-  Week2to4.sol[nrow(Week2to4.sol),8]
  D1_Inf_int <-  Week2to4.sol[nrow(Week2to4.sol),9]
  WF2Hog_Count1_int <-  Week2to4.sol[nrow(Week2to4.sol),10]
  
  S2_int <- Week2to4.sol[nrow(Week2to4.sol),11]
  E2_int <- Week2to4.sol[nrow(Week2to4.sol),12]
  I2_int <- Week2to4.sol[nrow(Week2to4.sol),13]
  R2_int <- Week2to4.sol[nrow(Week2to4.sol),14]
  Total_inf2_int <- Week2to4.sol[nrow(Week2to4.sol),15]
  Total_Rec2_int <- Week2to4.sol[nrow(Week2to4.sol),16]
  D2_Nat_int <- Week2to4.sol[nrow(Week2to4.sol),17]
  D2_Inf_int <- Week2to4.sol[nrow(Week2to4.sol),18]
  WF2Hog_Count2_int <- Week2to4.sol[nrow(Week2to4.sol),19]
  
  S3_int <- (1-Hog_Vacc_Eff)*1000
  E3_int <- 0
  I3_int <- 0
  R3_int <-  Hog_Vacc_Eff*1000
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week2to4.sol[nrow(Week2to4.sol),38]
  HE_int <-  Week2to4.sol[nrow(Week2to4.sol),39]
  HI_int <-  Week2to4.sol[nrow(Week2to4.sol),40]
  HR_int <-  Week2to4.sol[nrow(Week2to4.sol),41]
  H_DN_int <-  Week2to4.sol[nrow(Week2to4.sol),42]
  H_DI_int <-  Week2to4.sol[nrow(Week2to4.sol),43]
  Hog2WF_Count_int <-  Week2to4.sol[nrow(Week2to4.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 4 to 6
  #====================================================================================================
  Week4to6.t <- seq(28.1,Week4to6,0.1)         
  Week4to6.init <-c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                    S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                    S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                    S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                    HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week4to6.sol <- lsoda(Week4to6.init, Week4to6.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 6 TO 8 (Room 4 Fill)####
  Week6to8 <- 56
  #====================================================================================================
  # Model Initial Values  Weeks 6 to 8
  #===========================================================================
  S1_int <- Week4to6.sol[nrow(Week4to6.sol),2]
  E1_int <- Week4to6.sol[nrow(Week4to6.sol),3]
  I1_int <- Week4to6.sol[nrow(Week4to6.sol),4]
  R1_int <- Week4to6.sol[nrow(Week4to6.sol),5]
  Total_inf1_int <- Week4to6.sol[nrow(Week4to6.sol),6]
  Total_Rec1_int <- Week4to6.sol[nrow(Week4to6.sol),7]
  D1_Nat_int <-  Week4to6.sol[nrow(Week4to6.sol),8]
  D1_Inf_int <-  Week4to6.sol[nrow(Week4to6.sol),9]
  WF2Hog_Count1_int <-  Week4to6.sol[nrow(Week4to6.sol),10]
  
  S2_int <- Week4to6.sol[nrow(Week4to6.sol),11]
  E2_int <- Week4to6.sol[nrow(Week4to6.sol),12]
  I2_int <- Week4to6.sol[nrow(Week4to6.sol),13]
  R2_int <- Week4to6.sol[nrow(Week4to6.sol),14]
  Total_inf2_int <- Week4to6.sol[nrow(Week4to6.sol),15]
  Total_Rec2_int <- Week4to6.sol[nrow(Week4to6.sol),16]
  D2_Nat_int <- Week4to6.sol[nrow(Week4to6.sol),17]
  D2_Inf_int <- Week4to6.sol[nrow(Week4to6.sol),18]
  WF2Hog_Count2_int <- Week4to6.sol[nrow(Week4to6.sol),19]
  
  S3_int <- Week4to6.sol[nrow(Week4to6.sol),20]
  E3_int <- Week4to6.sol[nrow(Week4to6.sol),21]
  I3_int <- Week4to6.sol[nrow(Week4to6.sol),22]
  R3_int <-  Week4to6.sol[nrow(Week4to6.sol),23]
  Total_inf3_int <- Week4to6.sol[nrow(Week4to6.sol),24]
  Total_Rec3_int <- Week4to6.sol[nrow(Week4to6.sol),25]
  D3_Nat_int <- Week4to6.sol[nrow(Week4to6.sol),26]
  D3_Inf_int <- Week4to6.sol[nrow(Week4to6.sol),27]
  WF2Hog_Count3_int <- Week4to6.sol[nrow(Week4to6.sol),28]
  
  S4_int <- (1-Hog_Vacc_Eff)*1000
  E4_int <- 0
  I4_int <- 0
  R4_int <- Hog_Vacc_Eff*1000
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week4to6.sol[nrow(Week4to6.sol),38]
  HE_int <-  Week4to6.sol[nrow(Week4to6.sol),39]
  HI_int <-  Week4to6.sol[nrow(Week4to6.sol),40]
  HR_int <-  Week4to6.sol[nrow(Week4to6.sol),41]
  H_DN_int <-  Week4to6.sol[nrow(Week4to6.sol),42]
  H_DI_int <-  Week4to6.sol[nrow(Week4to6.sol),43]
  Hog2WF_Count_int <-  Week4to6.sol[nrow(Week4to6.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 6 to 8
  #====================================================================================================
  Week6to8.t <- seq(42.1,Week6to8,0.1)           
  Week6to8.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week6to8.sol <- lsoda(Week6to8.init, Week6to8.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 8 TO 23 (Room 1 Exit) ####
  Week8to23 <- 161
  #====================================================================================================
  # Model Initial Values  Weeks 8 to 23
  #===========================================================================
  S1_int <- Week6to8.sol[nrow(Week6to8.sol),2]
  E1_int <- Week6to8.sol[nrow(Week6to8.sol),3]
  I1_int <- Week6to8.sol[nrow(Week6to8.sol),4]
  R1_int <- Week6to8.sol[nrow(Week6to8.sol),5]
  Total_inf1_int <- Week6to8.sol[nrow(Week6to8.sol),6]
  Total_Rec1_int <- Week6to8.sol[nrow(Week6to8.sol),7]
  D1_Nat_int <-  Week6to8.sol[nrow(Week6to8.sol),8]
  D1_Inf_int <-  Week6to8.sol[nrow(Week6to8.sol),9]
  WF2Hog_Count1_int <-  Week6to8.sol[nrow(Week6to8.sol),10]
  
  
  S2_int <- Week6to8.sol[nrow(Week6to8.sol),11]
  E2_int <- Week6to8.sol[nrow(Week6to8.sol),12]
  I2_int <- Week6to8.sol[nrow(Week6to8.sol),13]
  R2_int <- Week6to8.sol[nrow(Week6to8.sol),14]
  Total_inf2_int <- Week6to8.sol[nrow(Week6to8.sol),15]
  Total_Rec2_int <- Week6to8.sol[nrow(Week6to8.sol),16]
  D2_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),17]
  D2_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),18]
  WF2Hog_Count2_int <- Week6to8.sol[nrow(Week6to8.sol),19]
  
  S3_int <- Week6to8.sol[nrow(Week6to8.sol),20]
  E3_int <- Week6to8.sol[nrow(Week6to8.sol),21]
  I3_int <- Week6to8.sol[nrow(Week6to8.sol),22]
  R3_int <-  Week6to8.sol[nrow(Week6to8.sol),23]
  Total_inf3_int <- Week6to8.sol[nrow(Week6to8.sol),24]
  Total_Rec3_int <- Week6to8.sol[nrow(Week6to8.sol),25]
  D3_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),26]
  D3_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),27]
  WF2Hog_Count3_int <- Week6to8.sol[nrow(Week6to8.sol),28]
  
  S4_int <- Week6to8.sol[nrow(Week6to8.sol),29]
  E4_int <- Week6to8.sol[nrow(Week6to8.sol),30]
  I4_int <- Week6to8.sol[nrow(Week6to8.sol),31]
  R4_int <- Week6to8.sol[nrow(Week6to8.sol),32]
  Total_inf4_int <- Week6to8.sol[nrow(Week6to8.sol),33]
  Total_Rec4_int <- Week6to8.sol[nrow(Week6to8.sol),34]
  D4_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),35]
  D4_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),36]
  WF2Hog_Count4_int <- Week6to8.sol[nrow(Week6to8.sol),37]
  
  HS_int <-  Week6to8.sol[nrow(Week6to8.sol),38]
  HE_int <-  Week6to8.sol[nrow(Week6to8.sol),39]
  HI_int <-  Week6to8.sol[nrow(Week6to8.sol),40]
  HR_int <-  Week6to8.sol[nrow(Week6to8.sol),41]
  H_DN_int <-  Week6to8.sol[nrow(Week6to8.sol),42]
  H_DI_int <-  Week6to8.sol[nrow(Week6to8.sol),43]
  Hog2WF_Count_int <-  Week6to8.sol[nrow(Week6to8.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 8 to 23
  #====================================================================================================
  Week8to23.t <- seq(56.1,Week8to23,0.1)          
  Week8to23.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                      S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                      S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                      S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                      HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week8to23.sol <- lsoda(Week8to23.init, Week8to23.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 23 TO 25 (Room 2 Exit) ####
  Week23to25 <- 175
  #====================================================================================================
  # Model Initial Values  Weeks 23 to 25
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week8to23.sol[nrow(Week8to23.sol),6]
  Total_Rec1_int <- Week8to23.sol[nrow(Week8to23.sol),7]
  D1_Nat_int <-  Week8to23.sol[nrow(Week8to23.sol),8]
  D1_Inf_int <-  Week8to23.sol[nrow(Week8to23.sol),9]
  WF2Hog_Count1_int <-  Week8to23.sol[nrow(Week8to23.sol),10]
  
  
  S2_int <- Week8to23.sol[nrow(Week8to23.sol),11]
  E2_int <- Week8to23.sol[nrow(Week8to23.sol),12]
  I2_int <- Week8to23.sol[nrow(Week8to23.sol),13]
  R2_int <- Week8to23.sol[nrow(Week8to23.sol),14]
  Total_inf2_int <- Week8to23.sol[nrow(Week8to23.sol),15]
  Total_Rec2_int <- Week8to23.sol[nrow(Week8to23.sol),16]
  D2_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),17]
  D2_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),18]
  WF2Hog_Count2_int <- Week8to23.sol[nrow(Week8to23.sol),19]
  
  S3_int <- Week8to23.sol[nrow(Week8to23.sol),20]
  E3_int <- Week8to23.sol[nrow(Week8to23.sol),21]
  I3_int <- Week8to23.sol[nrow(Week8to23.sol),22]
  R3_int <-  Week8to23.sol[nrow(Week8to23.sol),23]
  Total_inf3_int <- Week8to23.sol[nrow(Week8to23.sol),24]
  Total_Rec3_int <- Week8to23.sol[nrow(Week8to23.sol),25]
  D3_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),26]
  D3_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),27]
  WF2Hog_Count3_int <- Week8to23.sol[nrow(Week8to23.sol),28]
  
  S4_int <- Week8to23.sol[nrow(Week8to23.sol),29]
  E4_int <- Week8to23.sol[nrow(Week8to23.sol),30]
  I4_int <- Week8to23.sol[nrow(Week8to23.sol),31]
  R4_int <- Week8to23.sol[nrow(Week8to23.sol),32]
  Total_inf4_int <- Week8to23.sol[nrow(Week8to23.sol),33]
  Total_Rec4_int <- Week8to23.sol[nrow(Week8to23.sol),34]
  D4_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),35]
  D4_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),36]
  WF2Hog_Count4_int <- Week8to23.sol[nrow(Week8to23.sol),37]
  
  HS_int <-  Week8to23.sol[nrow(Week8to23.sol),38]
  HE_int <-  Week8to23.sol[nrow(Week8to23.sol),39]
  HI_int <-  Week8to23.sol[nrow(Week8to23.sol),40]
  HR_int <-  Week8to23.sol[nrow(Week8to23.sol),41]
  H_DN_int <-  Week8to23.sol[nrow(Week8to23.sol),42]
  H_DI_int <-  Week8to23.sol[nrow(Week8to23.sol),43]
  Hog2WF_Count_int <-  Week8to23.sol[nrow(Week8to23.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 23 to 25
  #====================================================================================================
  Week23to25.t <- seq(161.1,Week23to25,0.1)         
  Week23to25.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week23to25.sol <- lsoda(Week23to25.init, Week23to25.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 25 TO 27 (Room 3 Exit) ####
  Week25to27 <- 189
  #====================================================================================================
  # Model Initial Values  Weeks 25 to 27
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week23to25.sol[nrow(Week23to25.sol),6]
  Total_Rec1_int <- Week23to25.sol[nrow(Week23to25.sol),7]
  D1_Nat_int <-  Week23to25.sol[nrow(Week23to25.sol),8]
  D1_Inf_int <-  Week23to25.sol[nrow(Week23to25.sol),9]
  WF2Hog_Count1_int <-  Week23to25.sol[nrow(Week23to25.sol),10]
  
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- Week23to25.sol[nrow(Week23to25.sol),15]
  Total_Rec2_int <- Week23to25.sol[nrow(Week23to25.sol),16]
  D2_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),17]
  D2_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),18]
  WF2Hog_Count2_int <- Week23to25.sol[nrow(Week23to25.sol),19]
  
  S3_int <- Week23to25.sol[nrow(Week23to25.sol),20]
  E3_int <- Week23to25.sol[nrow(Week23to25.sol),21]
  I3_int <- Week23to25.sol[nrow(Week23to25.sol),22]
  R3_int <-  Week23to25.sol[nrow(Week23to25.sol),23]
  Total_inf3_int <- Week23to25.sol[nrow(Week23to25.sol),24]
  Total_Rec3_int <- Week23to25.sol[nrow(Week23to25.sol),25]
  D3_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),26]
  D3_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),27]
  WF2Hog_Count3_int <- Week23to25.sol[nrow(Week23to25.sol),28]
  
  S4_int <- Week23to25.sol[nrow(Week23to25.sol),29]
  E4_int <- Week23to25.sol[nrow(Week23to25.sol),30]
  I4_int <- Week23to25.sol[nrow(Week23to25.sol),31]
  R4_int <- Week23to25.sol[nrow(Week23to25.sol),32]
  Total_inf4_int <- Week23to25.sol[nrow(Week23to25.sol),33]
  Total_Rec4_int <- Week23to25.sol[nrow(Week23to25.sol),34]
  D4_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),35]
  D4_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),36]
  WF2Hog_Count4_int <- Week23to25.sol[nrow(Week23to25.sol),37]
  
  HS_int <-  Week23to25.sol[nrow(Week23to25.sol),38]
  HE_int <-  Week23to25.sol[nrow(Week23to25.sol),39]
  HI_int <-  Week23to25.sol[nrow(Week23to25.sol),40]
  HR_int <-  Week23to25.sol[nrow(Week23to25.sol),41]
  H_DN_int <-  Week23to25.sol[nrow(Week23to25.sol),42]
  H_DI_int <-  Week23to25.sol[nrow(Week23to25.sol),43]
  Hog2WF_Count_int <-  Week23to25.sol[nrow(Week23to25.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 25 to 27
  #====================================================================================================
  Week25to27.t <- seq(175.1,Week25to27,0.1)          
  Week25to27.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week25to27.sol <- lsoda(Week25to27.init, Week25to27.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 27 TO 29 (Room 4 Exit) ####
  Week27to29 <- 203
  #====================================================================================================
  # Model Initial Values  Weeks 27to 29
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week25to27.sol[nrow(Week25to27.sol),6]
  Total_Rec1_int <- Week25to27.sol[nrow(Week25to27.sol),7]
  D1_Nat_int <-  Week25to27.sol[nrow(Week25to27.sol),8]
  D1_Inf_int <-  Week25to27.sol[nrow(Week25to27.sol),9]
  WF2Hog_Count1_int <-  Week25to27.sol[nrow(Week25to27.sol),10]
  
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- Week25to27.sol[nrow(Week25to27.sol),15]
  Total_Rec2_int <- Week25to27.sol[nrow(Week25to27.sol),16]
  D2_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),17]
  D2_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),18]
  WF2Hog_Count2_int <- Week25to27.sol[nrow(Week25to27.sol),19]
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- Week25to27.sol[nrow(Week25to27.sol),24]
  Total_Rec3_int <- Week25to27.sol[nrow(Week25to27.sol),25]
  D3_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),26]
  D3_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),27]
  WF2Hog_Count3_int <- Week25to27.sol[nrow(Week25to27.sol),28]
  
  S4_int <- Week25to27.sol[nrow(Week25to27.sol),29]
  E4_int <- Week25to27.sol[nrow(Week25to27.sol),30]
  I4_int <- Week25to27.sol[nrow(Week25to27.sol),31]
  R4_int <- Week25to27.sol[nrow(Week25to27.sol),32]
  Total_inf4_int <- Week25to27.sol[nrow(Week25to27.sol),33]
  Total_Rec4_int <- Week25to27.sol[nrow(Week25to27.sol),34]
  D4_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),35]
  D4_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),36]
  WF2Hog_Count4_int <- Week25to27.sol[nrow(Week25to27.sol),37]
  
  HS_int <-  Week25to27.sol[nrow(Week25to27.sol),38]
  HE_int <-  Week25to27.sol[nrow(Week25to27.sol),39]
  HI_int <-  Week25to27.sol[nrow(Week25to27.sol),40]
  HR_int <-  Week25to27.sol[nrow(Week25to27.sol),41]
  H_DN_int <-  Week25to27.sol[nrow(Week25to27.sol),42]
  H_DI_int <-  Week25to27.sol[nrow(Week25to27.sol),43]
  Hog2WF_Count_int <-  Week25to27.sol[nrow(Week25to27.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 27 to 29
  #====================================================================================================
  Week27to29.t <- seq(189.1,Week27to29,0.1)            
  Week27to29.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week27to29.sol <- lsoda(Week27to29.init, Week27to29.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### COMBINING THE DATASETS ####
  Week0to2.sol.df <- as.data.frame(Week0to2.sol)
  colnames(Week0to2.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week2to4.sol.df <- as.data.frame(Week2to4.sol)
  colnames(Week2to4.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week4to6.sol.df <- as.data.frame(Week4to6.sol)
  colnames(Week4to6.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week6to8.sol.df <- as.data.frame(Week6to8.sol)
  colnames(Week6to8.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week8to23.sol.df <- as.data.frame(Week8to23.sol)
  colnames(Week8to23.sol.df) <- c("time",
                                  "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                  "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                  "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                  "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                  "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week23to25.sol.df <- as.data.frame(Week23to25.sol)
  colnames(Week23to25.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week25to27.sol.df <- as.data.frame(Week25to27.sol)
  colnames(Week25to27.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week27to29.sol.df <- as.data.frame(Week27to29.sol)
  colnames(Week27to29.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  
  fulldata.df <- rbind(Week0to2.sol.df, Week2to4.sol.df, Week4to6.sol.df, Week6to8.sol.df, 
                       Week8to23.sol.df, Week23to25.sol.df, Week25to27.sol.df, Week27to29.sol.df)
  
  ## Is there a WF infection
  HI_inf <- fulldata.df %>%
    summarise(total_I = sum(HI)) %>%
    mutate(HI_Inf = ifelse(total_I>1, 1, 0)) %>%
    select(HI_Inf) %>% as.numeric()
  #Time to 1st WF infection
  Time2WFinf0 <- fulldata.df %>%
    arrange(time) %>%
    filter(HI > 0) %>%
    filter(row_number() == 1) %>%
    mutate(TimeSineInfIntro = time - 0) %>%
    select(TimeSineInfIntro) %>% as.numeric()
  
  Time2WFinf.5 <- fulldata.df %>%
    arrange(time) %>%
    filter(HI > 0.5) %>%
    filter(row_number() == 1) %>%
    mutate(TimeSineInfIntro = time - 0) %>%
    select(TimeSineInfIntro) %>% as.numeric()
  ## Time to Second Room Infection
  Timeto2room <- fulldata.df %>% 
    summarise(t_cross = time[which(E4 > 1 | 
                                     E2 > 1 |
                                     E3 > 1)]) %>% 
    summarise(min_t_cross = min(t_cross)) %>% as.numeric()
  ## Total Hogs infected
  TotHogsInf <- fulldata.df %>% 
    #select(t, iter, Total_inf1, Total_inf2, Total_inf3, Total_inf4) %>% 
    mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
           Tot_rec = Total_Rec1 + Total_Rec2 + Total_Rec3 + Total_Rec4) %>% 
    summarise(Tot_inf2 = max(Tot_inf),
              Tot_rec2 = max(Tot_rec))
  ## Minimum Time to Peak Infection
  Time2PeakInf <- fulldata.df %>% 
    mutate(S_tot = S1+S2+S3+S4,                             
           E_tot = E1+E2+E3+E4,
           I_tot = I1+I2+I3+I4,
           R_tot = R1+R2+R3+R4) %>%                                     
    summarise(max_I = max(I_tot),
              t_peak = time[which(I_tot == max(I_tot))]) %>% 
    summarise(min_t_peak = min(t_peak)) %>% as.numeric()
  ##Ro
  Ro <- fulldata.df %>% 
    arrange(time) %>% 
    mutate(Tot_rec = Total_Rec1 + Total_Rec2 + Total_Rec3 + Total_Rec4) %>%  
    filter(row_number() == n()) %>% 
    select(Tot_rec) %>% 
    mutate(Ro = ((log((4000 - Tot_rec)/4000) - log((4000 - (Hog_Vacc_Eff*4000))/4000)) /
                   ((((4000 - Tot_rec) / 4000) - (4000 - (Hog_Vacc_Eff*4000)) / 4000)))) %>% 
    select(Ro) %>% as.numeric()
  
  #Create Summary data frame
  fulldata.df.sum <- data.frame(iter = as.factor(i),
                                HI_inf = HI_inf,
                                Time2WFinf0 = Time2WFinf0,
                                Time2WFinf.5 = Time2WFinf.5,
                                Timeto2room = Timeto2room,
                                TotHogsInf = TotHogsInf$Tot_inf2,
                                TotHogsRec = TotHogsInf$Tot_rec2,
                                Time2PeakInf = Time2PeakInf,
                                Ro = Ro,
                                beta_Hog = beta_Hog,
                                beta_Hog_ind = beta_Hog_ind,
                                beta_S2WF = beta_S2WF,
                                beta_WF2S = beta_WF2S,
                                sigma_Hog = sigma_Hog,
                                delta_Hog = delta_Hog,
                                beta_WF2WF = beta_WF2WF,
                                sigma_WF = sigma_WF,
                                delta_WF = delta_WF)
  
  results <- rbind(results, fulldata.df.sum)
  
  rm(list=ls()[! ls() %in% c("results", "n.itr")])  #Deletes everything from environment expect growing results dataframe
}

saveRDS(results, file = "Current code/Control Meausres/Sensativity Analysis/WFFlowRoom1_Sum_SensAnlaysis.rds")

#### WF intro virus ####
results <- NULL
n.itr <- 5000
# seed.num <- 616
seed.nums <- c(100:(100+n.itr))
for (i in c(1:n.itr)) {
  print(i)
  
  #### MODEL CONSTANTS ####
  #Hog constants
  beta_Hog <- rtri(n = 1, min = 0.001, max = 0.1, mode = (10/(1000*5)))
  beta_Hog_ind <- beta_Hog / 178 #(following Etbaigha et al 2018 paper)
  # beta_Hog_ind <- beta_Hog / 500
  Hog_Vacc_Eff <- 0
  u_Hog <- 0.00028     #natural death rate (following Etbaigha et al 2018 paper)
  u_inf_Hog <- 0 #0.1     #infected death rate (0 for now to problem solve)
  w_Hog <- 0 #1/180         #(following Etbaigha et al 2018 paper)  White et al shows a range from 56 - 112 days
  sigma_Hog <- abs(1 / rnorm(n = 1, mean = 2, sd = 1))   #Latent period h1n1 and h3n2
  delta_Hog <- abs(1 / rnorm(n = 1, mean = 5, sd = 1))   #Infectious period h1n1 and h3n2
  
  #Interspecies transmission
  ##applied the beta = Ro / N*D from the modeling book
  # Ro from the following website for swine to human infection 
  # https://bmcmedicine.biomedcentral.com/articles/10.1186/1741-7015-7-30#Sec6   (from super-strain figure) possible the high end of transmisilibty for this parameter
  # beta_S2WF <- (2.3 / (4002*5))  #  #Hogs + #Workforce * Duration of Hogs
  #Alt beta_S2WF from fair outbreak
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3557885/
  #assuming one hog got the first three confirmed cases sick and infectious period of 5 days Ro = 3/5 = 0.6
  beta_S2WF <- runif(1, min = 0.000001, max = 0.0001)
  #when Ro WF2S = Ro S2WF
  # beta_WF2S <- (2.3 / (4002*3))  #  #Hogs + #Workforce * Duration of Hogs
  # WF2S calculated from Swine Outbreak of Pandemic Influenza A Virus on a Canadian Research Farm Supports Human-to-Swine Transmission paper
  beta_WF2S <- runif(1, min = 0.000001, max = 0.0001)    #Hogs + #Workforce * Duration of WF
  #Human constants
  #applied the beta = Ro / N*D from the modeling book
  beta_WF2WF <-   rtri(n = 1, min = 0.0000001, max = 0.64, mode = 0.32) # http://m-hikari.com/ams/ams-2013/ams-41-44-2013/joseAMS41-44-2013.pdf 
  Hum_Vacc_Eff <- 0.0 
  u_WF <- 0     #Natural death rate of human (0 for now to problem solve)
  u_inf_WF <- 0   #Infected death rate of human (0 for now to problem solve)
  w_WF <- 0      #Recovery rate for human (0 for now to problem solve)
  sigma_WF <- abs(1 / rnorm(n = 1, mean = 2, sd = 1))   #Latent period for humans
  delta_WF <- abs(1 / rnorm(n = 1, mean = 3, sd = 1))   #Infectious Period for humans
  
  #### MODEL SETUP ####
  # Model Parameters
  TransmissionModel.par <- c(beta_Hog = beta_Hog, beta_Hog_ind = beta_Hog_ind, beta_S2WF = beta_S2WF, beta_WF2S = beta_WF2S,
                             w_Hog = w_Hog, sigma_Hog = sigma_Hog, delta_Hog = delta_Hog, u_Hog = u_Hog, u_inf_Hog = u_inf_Hog,
                             beta_WF2WF = beta_WF2WF,
                             w_WF = w_WF, sigma_WF = sigma_WF, delta_WF = delta_WF, u_WF = u_WF, u_inf_WF = u_inf_WF)
  #=====================================================================================================
  #Model Equations 
  #=====================================================================================================
  TransmissionModel.dyn <- function(t,var,par) {
    # Rename the parameters 
    beta_Hog <- par[1]
    beta_Hog_ind <- par[2]
    beta_S2WF <- par[3]
    beta_WF2S <- par[4]
    w_Hog <- par[5]
    sigma_Hog <- par[6]
    delta_Hog <- par[7]
    u_Hog <- par[8]
    u_inf_Hog<- par[9]
    beta_WF2WF <- par[10]
    w_WF <- par[11]
    sigma_WF <- par[12]
    delta_WF <- par[13]
    u_WF <- par[14]
    u_inf_WF <- par[15]
    
    S1 <- var[1]
    E1 <- var[2]
    I1 <- var[3]
    R1 <- var[4]
    Total_inf1 <- var[5] 
    True_Rec1 <- var[6]
    D1_Nat <- var[7]
    D1_Inf <- var[8]
    WF2Hog_Count1 <- var[9] 
    
    S2 <- var[10]
    E2 <- var[11]
    I2 <- var[12]
    R2 <- var[13]
    Total_inf2 <- var[14]
    True_Rec2 <- var[15]
    D2_Nat <- var[16]
    D2_Inf <- var[17]
    WF2Hog_Count2 <- var[18]
    
    S3 <- var[19]
    E3 <- var[20]
    I3 <- var[21]
    R3 <- var[22]
    Total_inf3 <- var[23]
    True_Rec3 <- var[24]
    D3_Nat <- var[25]
    D3_Inf <- var[26]
    WF2Hog_Count3 <- var[27]   
    
    S4 <- var[28]
    E4 <- var[29]
    I4 <- var[30]
    R4 <- var[31]
    Total_inf4 <- var[32]
    True_Rec4 <- var[33]
    D4_Nat <- var[34]
    D4_Inf <- var[35]
    WF2Hog_Count4 <- var[36]
    
    HS <- var[37]
    HE <- var[38]
    HI <- var[39] 
    HR <- var[40]
    H_DN <- var[41]
    H_DI <- var[42]
    Hog2WF_Count <- var[43] 
    
    # Calculate the derivatives
    #Room 1
    dS1 <- (w_Hog*R1) - (beta_Hog*I1*S1) - (beta_Hog_ind*(I2+I3+I4)*S1) - (beta_WF2S*HI*S1) - (u_Hog*S1)
    dE1 <-(beta_Hog*I1*S1) + (beta_Hog_ind*(I2+I3+I4)*S1) + (beta_WF2S*HI*S1) - sigma_Hog*E1 - (u_Hog*E1)
    dI1 <- (sigma_Hog*E1) - (delta_Hog*I1) - (u_inf_Hog*I1)
    dR1 <- (delta_Hog*I1) - (w_Hog*R1) - (u_Hog*R1)
    dTotal_inf1 <- (sigma_Hog*E1)
    dTotal_Rec1 <- (delta_Hog*I1)
    dD1_Nat <- (u_Hog*S1) + (u_Hog*E1) + (u_Hog*R1)
    dD1_Inf <- (u_inf_Hog*I1)
    dWF2Hog_Count1 <- (beta_WF2S*HI*S1)
    
    #Room 2
    dS2 <- (w_Hog*R2) - (beta_Hog*I2*S2) - (beta_Hog_ind*(I1+I3+I4)*S2) - (beta_WF2S*HI*S2) - (u_Hog*S2)
    dE2 <-(beta_Hog*I2*S2) + (beta_Hog_ind*(I1+I3+I4)*S2) + (beta_WF2S*HI*S2) - (sigma_Hog*E2) - (u_Hog*E2) 
    dI2 <- (sigma_Hog*E2) - (delta_Hog*I2) - (u_inf_Hog*I2)
    dR2 <- (delta_Hog*I2) - (w_Hog*R2) - (u_Hog*R2)
    dTotal_inf2 <- (sigma_Hog*E2)
    dTotal_Rec2 <- (delta_Hog*I2)
    dD2_Nat <- (u_Hog*S2) + (u_Hog*E2) + (u_Hog*R2)
    dD2_Inf <- (u_inf_Hog*I2)
    dWF2Hog_Count2 <- (beta_WF2S*HI*S2)
    
    #Room 3
    dS3 <- (w_Hog*R3) - (beta_Hog*I3*S3) - (beta_Hog_ind*(I1+I2+I4)*S3) - (beta_WF2S*HI*S3) - (u_Hog*S3)
    dE3 <-(beta_Hog*I3*S3) + (beta_Hog_ind*(I1+I2+I4)*S3) + (beta_WF2S*HI*S3) - (sigma_Hog*E3) - (u_Hog*E3)
    dI3 <- (sigma_Hog*E3) - (delta_Hog*I3) - (u_inf_Hog*I3)
    dR3 <- (delta_Hog*I3) - (w_Hog*R3) - (u_Hog*R3)
    dTotal_inf3 <- (sigma_Hog*E3)
    dTotal_Rec3 <- (delta_Hog*I3)
    dD3_Nat <- (u_Hog*S3) + (u_Hog*E3) + (u_Hog*R3)
    dD3_Inf <- (u_inf_Hog*I3)
    dWF2Hog_Count3 <- (beta_WF2S*HI*S3)
    
    #Room 4
    dS4 <- w_Hog*R4 - (beta_Hog*I4*S4) - (beta_Hog_ind*(I1+I2+I3)*S4) - (beta_WF2S*HI*S4) - (u_Hog*S4)
    dE4 <-(beta_Hog*I4*S4) + (beta_Hog_ind*(I1+I2+I3)*S4) + (beta_WF2S*HI*S4) - (sigma_Hog*E4) - (u_Hog*E4) 
    dI4 <- (sigma_Hog*E4) - (delta_Hog*I4) - (u_inf_Hog*I4)
    dR4 <- (delta_Hog*I4) - (w_Hog*R4) - (u_Hog*R4)
    dTotal_inf4 <- (sigma_Hog*E4)
    dTotal_Rec4 <- (delta_Hog*I4)
    dD4_Nat <- (u_Hog*S4) + (u_Hog*E4) + (u_Hog*R4)
    dD4_Inf <- (u_inf_Hog*I4)
    dWF2Hog_Count4 <- (beta_WF2S*HI*S4)
    
    #Humans
    dHS <- (w_WF*HR)- (beta_WF2WF*HI*HS) - (beta_S2WF*(I1+I2+I3+I4)*HS) - (u_WF*HS)
    dHE <- (beta_WF2WF*HI*HS) + (beta_S2WF*(I1+I2+I3+I4)*HS) - (sigma_WF*HE) - (u_WF*HE)
    dHI <- (sigma_WF*HE) - (delta_WF*HI) - (u_inf_WF*HI)
    dHR <- (delta_WF*HI) - (w_WF*HR) - u_WF*HR
    dH_DN <- (u_WF*HS) + (u_WF*HE) + (u_WF*HR)
    dH_DI <- (u_inf_WF*HI)
    dHog2WF_Count <- (beta_S2WF*(I1+I2+I3+I4)*HS)
    
    
    # Last instruction: return a list 
    
    return(list(c(dS1, dE1, dI1, dR1, dTotal_inf1, dTotal_Rec1, dD1_Nat, dD1_Inf, dWF2Hog_Count1,
                  dS2, dE2, dI2, dR2, dTotal_inf2, dTotal_Rec2, dD2_Nat, dD2_Inf, dWF2Hog_Count2,
                  dS3, dE3, dI3, dR3, dTotal_inf3, dTotal_Rec3, dD3_Nat, dD3_Inf, dWF2Hog_Count3,
                  dS4, dE4, dI4, dR4, dTotal_inf4, dTotal_Rec4, dD4_Nat, dD4_Inf, dWF2Hog_Count4,
                  dHS, dHE, dHI, dHR, dH_DN, dH_DI, dHog2WF_Count)))
    
  }
  
  #### MODEL ITERATIONS ####
  #### WEEK 0 TO 2 (Room 1 fill) ####
  Week0to2 <- 14
  #====================================================================================================
  # Model Initial Values  Weeks 0 to 2
  #===========================================================================
  S1_int <- (1-Hog_Vacc_Eff)*1000
  E1_int <- 0
  I1_int <- 0 #1
  R1_int <- Hog_Vacc_Eff*1000
  Total_inf1_int <- 0
  Total_Rec1_int <- 0
  D1_Nat_int <- 0
  D1_Inf_int <- 0
  WF2Hog_Count1_int <- 0
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- 0
  Total_Rec2_int <- 0
  D2_Nat_int <- 0
  D2_Inf_int <- 0
  WF2Hog_Count2_int <- 0
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <- 1
  HE_int <- 0
  HI_int <- 1
  HR_int <- 0
  H_DN_int <- 0
  H_DI_int <- 0
  Hog2WF_Count_int <- 0
  
  #====================================================================================================
  # Model Run  Weeks 0 to 2
  #====================================================================================================
  Week0to2.t <- seq(0,Week0to2,0.1)   
  Week0to2.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week0to2.sol <- lsoda(Week0to2.init, Week0to2.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 2 TO 4 (Room 2 Fill)####
  Week2to4 <- 28
  #====================================================================================================
  # Model Initial Values  Weeks 2 to 4
  #===========================================================================
  S1_int <- Week0to2.sol[nrow(Week0to2.sol),2]
  E1_int <- Week0to2.sol[nrow(Week0to2.sol),3]
  I1_int <- Week0to2.sol[nrow(Week0to2.sol),4]
  R1_int <- Week0to2.sol[nrow(Week0to2.sol),5]
  Total_inf1_int <-  Week0to2.sol[nrow(Week0to2.sol),6]
  Total_Rec1_int <-  Week0to2.sol[nrow(Week0to2.sol),7]
  D1_Nat_int <-  Week0to2.sol[nrow(Week0to2.sol),8]
  D1_Inf_int <-  Week0to2.sol[nrow(Week0to2.sol),9]
  WF2Hog_Count1_int <-  Week0to2.sol[nrow(Week0to2.sol),10]
  
  S2_int <- (1-Hog_Vacc_Eff)*1000
  E2_int <- 0
  I2_int <- 0 #1
  R2_int <- Hog_Vacc_Eff*1000
  Total_inf2_int <- 0
  Total_Rec2_int <- 0
  D2_Nat_int <- 0
  D2_Inf_int <- 0
  WF2Hog_Count2_int <- 0
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week0to2.sol[nrow(Week0to2.sol),38]
  HE_int <-  Week0to2.sol[nrow(Week0to2.sol),39]
  HI_int <-  Week0to2.sol[nrow(Week0to2.sol),40]
  HR_int <-  Week0to2.sol[nrow(Week0to2.sol),41]
  H_DN_int <-  Week0to2.sol[nrow(Week0to2.sol),42]
  H_DI_int <-  Week0to2.sol[nrow(Week0to2.sol),43]
  Hog2WF_Count_int <-  Week0to2.sol[nrow(Week0to2.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 2 to 4
  #====================================================================================================
  Week2to4.t <- seq(14.1,Week2to4,0.1)          
  Week2to4.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week2to4.sol <- lsoda(Week2to4.init, Week2to4.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 4 TO 6 (Room 3 Fill)####
  Week4to6 <- 42
  #====================================================================================================
  # Model Initial Values  Weeks 4 to 6
  #===========================================================================
  S1_int <- Week2to4.sol[nrow(Week2to4.sol),2]
  E1_int <- Week2to4.sol[nrow(Week2to4.sol),3]
  I1_int <- Week2to4.sol[nrow(Week2to4.sol),4]
  R1_int <- Week2to4.sol[nrow(Week2to4.sol),5]
  Total_inf1_int <- Week2to4.sol[nrow(Week2to4.sol),6]
  Total_Rec1_int <- Week2to4.sol[nrow(Week2to4.sol),7]
  D1_Nat_int <-  Week2to4.sol[nrow(Week2to4.sol),8]
  D1_Inf_int <-  Week2to4.sol[nrow(Week2to4.sol),9]
  WF2Hog_Count1_int <-  Week2to4.sol[nrow(Week2to4.sol),10]
  
  S2_int <- Week2to4.sol[nrow(Week2to4.sol),11]
  E2_int <- Week2to4.sol[nrow(Week2to4.sol),12]
  I2_int <- Week2to4.sol[nrow(Week2to4.sol),13]
  R2_int <- Week2to4.sol[nrow(Week2to4.sol),14]
  Total_inf2_int <- Week2to4.sol[nrow(Week2to4.sol),15]
  Total_Rec2_int <- Week2to4.sol[nrow(Week2to4.sol),16]
  D2_Nat_int <- Week2to4.sol[nrow(Week2to4.sol),17]
  D2_Inf_int <- Week2to4.sol[nrow(Week2to4.sol),18]
  WF2Hog_Count2_int <- Week2to4.sol[nrow(Week2to4.sol),19]
  
  S3_int <- (1-Hog_Vacc_Eff)*1000
  E3_int <- 0
  I3_int <- 0
  R3_int <-  Hog_Vacc_Eff*1000
  Total_inf3_int <- 0
  Total_Rec3_int <- 0
  D3_Nat_int <- 0
  D3_Inf_int <- 0
  WF2Hog_Count3_int <- 0
  
  S4_int <- 0
  E4_int <- 0
  I4_int <- 0
  R4_int <- 0
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week2to4.sol[nrow(Week2to4.sol),38]
  HE_int <-  Week2to4.sol[nrow(Week2to4.sol),39]
  HI_int <-  Week2to4.sol[nrow(Week2to4.sol),40]
  HR_int <-  Week2to4.sol[nrow(Week2to4.sol),41]
  H_DN_int <-  Week2to4.sol[nrow(Week2to4.sol),42]
  H_DI_int <-  Week2to4.sol[nrow(Week2to4.sol),43]
  Hog2WF_Count_int <-  Week2to4.sol[nrow(Week2to4.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 4 to 6
  #====================================================================================================
  Week4to6.t <- seq(28.1,Week4to6,0.1)         
  Week4to6.init <-c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                    S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                    S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                    S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                    HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week4to6.sol <- lsoda(Week4to6.init, Week4to6.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 6 TO 8 (Room 4 Fill)####
  Week6to8 <- 56
  #====================================================================================================
  # Model Initial Values  Weeks 6 to 8
  #===========================================================================
  S1_int <- Week4to6.sol[nrow(Week4to6.sol),2]
  E1_int <- Week4to6.sol[nrow(Week4to6.sol),3]
  I1_int <- Week4to6.sol[nrow(Week4to6.sol),4]
  R1_int <- Week4to6.sol[nrow(Week4to6.sol),5]
  Total_inf1_int <- Week4to6.sol[nrow(Week4to6.sol),6]
  Total_Rec1_int <- Week4to6.sol[nrow(Week4to6.sol),7]
  D1_Nat_int <-  Week4to6.sol[nrow(Week4to6.sol),8]
  D1_Inf_int <-  Week4to6.sol[nrow(Week4to6.sol),9]
  WF2Hog_Count1_int <-  Week4to6.sol[nrow(Week4to6.sol),10]
  
  S2_int <- Week4to6.sol[nrow(Week4to6.sol),11]
  E2_int <- Week4to6.sol[nrow(Week4to6.sol),12]
  I2_int <- Week4to6.sol[nrow(Week4to6.sol),13]
  R2_int <- Week4to6.sol[nrow(Week4to6.sol),14]
  Total_inf2_int <- Week4to6.sol[nrow(Week4to6.sol),15]
  Total_Rec2_int <- Week4to6.sol[nrow(Week4to6.sol),16]
  D2_Nat_int <- Week4to6.sol[nrow(Week4to6.sol),17]
  D2_Inf_int <- Week4to6.sol[nrow(Week4to6.sol),18]
  WF2Hog_Count2_int <- Week4to6.sol[nrow(Week4to6.sol),19]
  
  S3_int <- Week4to6.sol[nrow(Week4to6.sol),20]
  E3_int <- Week4to6.sol[nrow(Week4to6.sol),21]
  I3_int <- Week4to6.sol[nrow(Week4to6.sol),22]
  R3_int <-  Week4to6.sol[nrow(Week4to6.sol),23]
  Total_inf3_int <- Week4to6.sol[nrow(Week4to6.sol),24]
  Total_Rec3_int <- Week4to6.sol[nrow(Week4to6.sol),25]
  D3_Nat_int <- Week4to6.sol[nrow(Week4to6.sol),26]
  D3_Inf_int <- Week4to6.sol[nrow(Week4to6.sol),27]
  WF2Hog_Count3_int <- Week4to6.sol[nrow(Week4to6.sol),28]
  
  S4_int <- (1-Hog_Vacc_Eff)*1000
  E4_int <- 0
  I4_int <- 0
  R4_int <- Hog_Vacc_Eff*1000
  Total_inf4_int <- 0
  Total_Rec4_int <- 0
  D4_Nat_int <- 0
  D4_Inf_int <- 0
  WF2Hog_Count4_int <- 0
  
  HS_int <-  Week4to6.sol[nrow(Week4to6.sol),38]
  HE_int <-  Week4to6.sol[nrow(Week4to6.sol),39]
  HI_int <-  Week4to6.sol[nrow(Week4to6.sol),40]
  HR_int <-  Week4to6.sol[nrow(Week4to6.sol),41]
  H_DN_int <-  Week4to6.sol[nrow(Week4to6.sol),42]
  H_DI_int <-  Week4to6.sol[nrow(Week4to6.sol),43]
  Hog2WF_Count_int <-  Week4to6.sol[nrow(Week4to6.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 6 to 8
  #====================================================================================================
  Week6to8.t <- seq(42.1,Week6to8,0.1)           
  Week6to8.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                     S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                     S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                     S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                     HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week6to8.sol <- lsoda(Week6to8.init, Week6to8.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 8 TO 23 (Room 1 Exit) ####
  Week8to23 <- 161
  #====================================================================================================
  # Model Initial Values  Weeks 8 to 23
  #===========================================================================
  S1_int <- Week6to8.sol[nrow(Week6to8.sol),2]
  E1_int <- Week6to8.sol[nrow(Week6to8.sol),3]
  I1_int <- Week6to8.sol[nrow(Week6to8.sol),4]
  R1_int <- Week6to8.sol[nrow(Week6to8.sol),5]
  Total_inf1_int <- Week6to8.sol[nrow(Week6to8.sol),6]
  Total_Rec1_int <- Week6to8.sol[nrow(Week6to8.sol),7]
  D1_Nat_int <-  Week6to8.sol[nrow(Week6to8.sol),8]
  D1_Inf_int <-  Week6to8.sol[nrow(Week6to8.sol),9]
  WF2Hog_Count1_int <-  Week6to8.sol[nrow(Week6to8.sol),10]
  
  
  S2_int <- Week6to8.sol[nrow(Week6to8.sol),11]
  E2_int <- Week6to8.sol[nrow(Week6to8.sol),12]
  I2_int <- Week6to8.sol[nrow(Week6to8.sol),13]
  R2_int <- Week6to8.sol[nrow(Week6to8.sol),14]
  Total_inf2_int <- Week6to8.sol[nrow(Week6to8.sol),15]
  Total_Rec2_int <- Week6to8.sol[nrow(Week6to8.sol),16]
  D2_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),17]
  D2_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),18]
  WF2Hog_Count2_int <- Week6to8.sol[nrow(Week6to8.sol),19]
  
  S3_int <- Week6to8.sol[nrow(Week6to8.sol),20]
  E3_int <- Week6to8.sol[nrow(Week6to8.sol),21]
  I3_int <- Week6to8.sol[nrow(Week6to8.sol),22]
  R3_int <-  Week6to8.sol[nrow(Week6to8.sol),23]
  Total_inf3_int <- Week6to8.sol[nrow(Week6to8.sol),24]
  Total_Rec3_int <- Week6to8.sol[nrow(Week6to8.sol),25]
  D3_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),26]
  D3_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),27]
  WF2Hog_Count3_int <- Week6to8.sol[nrow(Week6to8.sol),28]
  
  S4_int <- Week6to8.sol[nrow(Week6to8.sol),29]
  E4_int <- Week6to8.sol[nrow(Week6to8.sol),30]
  I4_int <- Week6to8.sol[nrow(Week6to8.sol),31]
  R4_int <- Week6to8.sol[nrow(Week6to8.sol),32]
  Total_inf4_int <- Week6to8.sol[nrow(Week6to8.sol),33]
  Total_Rec4_int <- Week6to8.sol[nrow(Week6to8.sol),34]
  D4_Nat_int <- Week6to8.sol[nrow(Week6to8.sol),35]
  D4_Inf_int <- Week6to8.sol[nrow(Week6to8.sol),36]
  WF2Hog_Count4_int <- Week6to8.sol[nrow(Week6to8.sol),37]
  
  HS_int <-  Week6to8.sol[nrow(Week6to8.sol),38]
  HE_int <-  Week6to8.sol[nrow(Week6to8.sol),39]
  HI_int <-  Week6to8.sol[nrow(Week6to8.sol),40]
  HR_int <-  Week6to8.sol[nrow(Week6to8.sol),41]
  H_DN_int <-  Week6to8.sol[nrow(Week6to8.sol),42]
  H_DI_int <-  Week6to8.sol[nrow(Week6to8.sol),43]
  Hog2WF_Count_int <-  Week6to8.sol[nrow(Week6to8.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 8 to 23
  #====================================================================================================
  Week8to23.t <- seq(56.1,Week8to23,0.1)          
  Week8to23.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                      S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                      S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                      S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                      HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week8to23.sol <- lsoda(Week8to23.init, Week8to23.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 23 TO 25 (Room 2 Exit) ####
  Week23to25 <- 175
  #====================================================================================================
  # Model Initial Values  Weeks 23 to 25
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week8to23.sol[nrow(Week8to23.sol),6]
  Total_Rec1_int <- Week8to23.sol[nrow(Week8to23.sol),7]
  D1_Nat_int <-  Week8to23.sol[nrow(Week8to23.sol),8]
  D1_Inf_int <-  Week8to23.sol[nrow(Week8to23.sol),9]
  WF2Hog_Count1_int <-  Week8to23.sol[nrow(Week8to23.sol),10]
  
  
  S2_int <- Week8to23.sol[nrow(Week8to23.sol),11]
  E2_int <- Week8to23.sol[nrow(Week8to23.sol),12]
  I2_int <- Week8to23.sol[nrow(Week8to23.sol),13]
  R2_int <- Week8to23.sol[nrow(Week8to23.sol),14]
  Total_inf2_int <- Week8to23.sol[nrow(Week8to23.sol),15]
  Total_Rec2_int <- Week8to23.sol[nrow(Week8to23.sol),16]
  D2_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),17]
  D2_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),18]
  WF2Hog_Count2_int <- Week8to23.sol[nrow(Week8to23.sol),19]
  
  S3_int <- Week8to23.sol[nrow(Week8to23.sol),20]
  E3_int <- Week8to23.sol[nrow(Week8to23.sol),21]
  I3_int <- Week8to23.sol[nrow(Week8to23.sol),22]
  R3_int <-  Week8to23.sol[nrow(Week8to23.sol),23]
  Total_inf3_int <- Week8to23.sol[nrow(Week8to23.sol),24]
  Total_Rec3_int <- Week8to23.sol[nrow(Week8to23.sol),25]
  D3_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),26]
  D3_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),27]
  WF2Hog_Count3_int <- Week8to23.sol[nrow(Week8to23.sol),28]
  
  S4_int <- Week8to23.sol[nrow(Week8to23.sol),29]
  E4_int <- Week8to23.sol[nrow(Week8to23.sol),30]
  I4_int <- Week8to23.sol[nrow(Week8to23.sol),31]
  R4_int <- Week8to23.sol[nrow(Week8to23.sol),32]
  Total_inf4_int <- Week8to23.sol[nrow(Week8to23.sol),33]
  Total_Rec4_int <- Week8to23.sol[nrow(Week8to23.sol),34]
  D4_Nat_int <- Week8to23.sol[nrow(Week8to23.sol),35]
  D4_Inf_int <- Week8to23.sol[nrow(Week8to23.sol),36]
  WF2Hog_Count4_int <- Week8to23.sol[nrow(Week8to23.sol),37]
  
  HS_int <-  Week8to23.sol[nrow(Week8to23.sol),38]
  HE_int <-  Week8to23.sol[nrow(Week8to23.sol),39]
  HI_int <-  Week8to23.sol[nrow(Week8to23.sol),40]
  HR_int <-  Week8to23.sol[nrow(Week8to23.sol),41]
  H_DN_int <-  Week8to23.sol[nrow(Week8to23.sol),42]
  H_DI_int <-  Week8to23.sol[nrow(Week8to23.sol),43]
  Hog2WF_Count_int <-  Week8to23.sol[nrow(Week8to23.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 23 to 25
  #====================================================================================================
  Week23to25.t <- seq(161.1,Week23to25,0.1)         
  Week23to25.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week23to25.sol <- lsoda(Week23to25.init, Week23to25.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 25 TO 27 (Room 3 Exit) ####
  Week25to27 <- 189
  #====================================================================================================
  # Model Initial Values  Weeks 25 to 27
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week23to25.sol[nrow(Week23to25.sol),6]
  Total_Rec1_int <- Week23to25.sol[nrow(Week23to25.sol),7]
  D1_Nat_int <-  Week23to25.sol[nrow(Week23to25.sol),8]
  D1_Inf_int <-  Week23to25.sol[nrow(Week23to25.sol),9]
  WF2Hog_Count1_int <-  Week23to25.sol[nrow(Week23to25.sol),10]
  
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- Week23to25.sol[nrow(Week23to25.sol),15]
  Total_Rec2_int <- Week23to25.sol[nrow(Week23to25.sol),16]
  D2_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),17]
  D2_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),18]
  WF2Hog_Count2_int <- Week23to25.sol[nrow(Week23to25.sol),19]
  
  S3_int <- Week23to25.sol[nrow(Week23to25.sol),20]
  E3_int <- Week23to25.sol[nrow(Week23to25.sol),21]
  I3_int <- Week23to25.sol[nrow(Week23to25.sol),22]
  R3_int <-  Week23to25.sol[nrow(Week23to25.sol),23]
  Total_inf3_int <- Week23to25.sol[nrow(Week23to25.sol),24]
  Total_Rec3_int <- Week23to25.sol[nrow(Week23to25.sol),25]
  D3_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),26]
  D3_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),27]
  WF2Hog_Count3_int <- Week23to25.sol[nrow(Week23to25.sol),28]
  
  S4_int <- Week23to25.sol[nrow(Week23to25.sol),29]
  E4_int <- Week23to25.sol[nrow(Week23to25.sol),30]
  I4_int <- Week23to25.sol[nrow(Week23to25.sol),31]
  R4_int <- Week23to25.sol[nrow(Week23to25.sol),32]
  Total_inf4_int <- Week23to25.sol[nrow(Week23to25.sol),33]
  Total_Rec4_int <- Week23to25.sol[nrow(Week23to25.sol),34]
  D4_Nat_int <- Week23to25.sol[nrow(Week23to25.sol),35]
  D4_Inf_int <- Week23to25.sol[nrow(Week23to25.sol),36]
  WF2Hog_Count4_int <- Week23to25.sol[nrow(Week23to25.sol),37]
  
  HS_int <-  Week23to25.sol[nrow(Week23to25.sol),38]
  HE_int <-  Week23to25.sol[nrow(Week23to25.sol),39]
  HI_int <-  Week23to25.sol[nrow(Week23to25.sol),40]
  HR_int <-  Week23to25.sol[nrow(Week23to25.sol),41]
  H_DN_int <-  Week23to25.sol[nrow(Week23to25.sol),42]
  H_DI_int <-  Week23to25.sol[nrow(Week23to25.sol),43]
  Hog2WF_Count_int <-  Week23to25.sol[nrow(Week23to25.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 25 to 27
  #====================================================================================================
  Week25to27.t <- seq(175.1,Week25to27,0.1)          
  Week25to27.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week25to27.sol <- lsoda(Week25to27.init, Week25to27.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### WEEK 27 TO 29 (Room 4 Exit) ####
  Week27to29 <- 203
  #====================================================================================================
  # Model Initial Values  Weeks 27to 29
  #===========================================================================
  S1_int <- 0
  E1_int <- 0
  I1_int <- 0
  R1_int <- 0
  Total_inf1_int <- Week25to27.sol[nrow(Week25to27.sol),6]
  Total_Rec1_int <- Week25to27.sol[nrow(Week25to27.sol),7]
  D1_Nat_int <-  Week25to27.sol[nrow(Week25to27.sol),8]
  D1_Inf_int <-  Week25to27.sol[nrow(Week25to27.sol),9]
  WF2Hog_Count1_int <-  Week25to27.sol[nrow(Week25to27.sol),10]
  
  
  S2_int <- 0
  E2_int <- 0
  I2_int <- 0
  R2_int <- 0
  Total_inf2_int <- Week25to27.sol[nrow(Week25to27.sol),15]
  Total_Rec2_int <- Week25to27.sol[nrow(Week25to27.sol),16]
  D2_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),17]
  D2_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),18]
  WF2Hog_Count2_int <- Week25to27.sol[nrow(Week25to27.sol),19]
  
  S3_int <- 0
  E3_int <- 0
  I3_int <- 0
  R3_int <- 0
  Total_inf3_int <- Week25to27.sol[nrow(Week25to27.sol),24]
  Total_Rec3_int <- Week25to27.sol[nrow(Week25to27.sol),25]
  D3_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),26]
  D3_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),27]
  WF2Hog_Count3_int <- Week25to27.sol[nrow(Week25to27.sol),28]
  
  S4_int <- Week25to27.sol[nrow(Week25to27.sol),29]
  E4_int <- Week25to27.sol[nrow(Week25to27.sol),30]
  I4_int <- Week25to27.sol[nrow(Week25to27.sol),31]
  R4_int <- Week25to27.sol[nrow(Week25to27.sol),32]
  Total_inf4_int <- Week25to27.sol[nrow(Week25to27.sol),33]
  Total_Rec4_int <- Week25to27.sol[nrow(Week25to27.sol),34]
  D4_Nat_int <- Week25to27.sol[nrow(Week25to27.sol),35]
  D4_Inf_int <- Week25to27.sol[nrow(Week25to27.sol),36]
  WF2Hog_Count4_int <- Week25to27.sol[nrow(Week25to27.sol),37]
  
  HS_int <-  Week25to27.sol[nrow(Week25to27.sol),38]
  HE_int <-  Week25to27.sol[nrow(Week25to27.sol),39]
  HI_int <-  Week25to27.sol[nrow(Week25to27.sol),40]
  HR_int <-  Week25to27.sol[nrow(Week25to27.sol),41]
  H_DN_int <-  Week25to27.sol[nrow(Week25to27.sol),42]
  H_DI_int <-  Week25to27.sol[nrow(Week25to27.sol),43]
  Hog2WF_Count_int <-  Week25to27.sol[nrow(Week25to27.sol),44]
  
  #====================================================================================================
  # Model Run  Weeks 27 to 29
  #====================================================================================================
  Week27to29.t <- seq(189.1,Week27to29,0.1)            
  Week27to29.init <- c(S1_int, E1_int, I1_int, R1_int, Total_inf1_int, Total_Rec1_int, D1_Nat_int, D1_Inf_int, WF2Hog_Count1_int,
                       S2_int, E2_int, I2_int, R2_int, Total_inf2_int, Total_Rec2_int, D2_Nat_int, D2_Inf_int, WF2Hog_Count2_int,
                       S3_int, E3_int, I3_int, R3_int, Total_inf3_int, Total_Rec3_int, D3_Nat_int, D3_Inf_int, WF2Hog_Count3_int,
                       S4_int, E4_int, I4_int, R4_int, Total_inf4_int, Total_Rec4_int, D4_Nat_int, D4_Inf_int, WF2Hog_Count4_int,
                       HS_int, HE_int, HI_int, HR_int, H_DN_int, H_DI_int, Hog2WF_Count_int)
  Week27to29.sol <- lsoda(Week27to29.init, Week27to29.t, TransmissionModel.dyn, TransmissionModel.par)
  
  #### COMBINING THE DATASETS ####
  Week0to2.sol.df <- as.data.frame(Week0to2.sol)
  colnames(Week0to2.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week2to4.sol.df <- as.data.frame(Week2to4.sol)
  colnames(Week2to4.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week4to6.sol.df <- as.data.frame(Week4to6.sol)
  colnames(Week4to6.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week6to8.sol.df <- as.data.frame(Week6to8.sol)
  colnames(Week6to8.sol.df) <- c("time",
                                 "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                 "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                 "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                 "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                 "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week8to23.sol.df <- as.data.frame(Week8to23.sol)
  colnames(Week8to23.sol.df) <- c("time",
                                  "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                  "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                  "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                  "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                  "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week23to25.sol.df <- as.data.frame(Week23to25.sol)
  colnames(Week23to25.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week25to27.sol.df <- as.data.frame(Week25to27.sol)
  colnames(Week25to27.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  Week27to29.sol.df <- as.data.frame(Week27to29.sol)
  colnames(Week27to29.sol.df) <- c("time",
                                   "S1", "E1", "I1", "R1", "Total_inf1", "Total_Rec1", "D1_Nat", "D1_Inf", "WF2Hog_Count1",
                                   "S2", "E2", "I2", "R2", "Total_inf2", "Total_Rec2", "D2_Nat", "D2_Inf", "WF2Hog_Count2",
                                   "S3", "E3", "I3", "R3", "Total_inf3", "Total_Rec3", "D3_Nat", "D3_Inf", "WF2Hog_Count3",
                                   "S4", "E4", "I4", "R4", "Total_inf4", "Total_Rec4", "D4_Nat", "D4_Inf", "WF2Hog_Count4",
                                   "HS", "HE", "HI", "HR", "H_DN", "H_DIf", "Hog2WF_Count")
  
  fulldata.df <- rbind(Week0to2.sol.df, Week2to4.sol.df, Week4to6.sol.df, Week6to8.sol.df, 
                       Week8to23.sol.df, Week23to25.sol.df, Week25to27.sol.df, Week27to29.sol.df)
  
  ## Is there a WF infection
  HI_inf <- fulldata.df %>%
    summarise(total_I = sum(HI)) %>%
    mutate(HI_Inf = ifelse(total_I>0, 1, 0)) %>%
    select(HI_Inf) %>% as.numeric()
  #Time to 1st WF infection
  Time2WFinf0 <- fulldata.df %>%
    arrange(time) %>%
    filter(HI > 0) %>%
    filter(row_number() == 1) %>%
    mutate(TimeSineInfIntro = time) %>%
    select(TimeSineInfIntro) %>% as.numeric()
  
  Time2WFinf.5 <- fulldata.df %>%
    arrange(time) %>%
    filter(HI > 0.5) %>%
    filter(row_number() == 1) %>%
    mutate(TimeSineInfIntro = time) %>%
    select(TimeSineInfIntro) %>% as.numeric()
  ## Time to Second Room Infection
  Timeto2room <- fulldata.df %>% 
    summarise(t_cross = time[which(E1 >= 1 | 
                                     E2 >= 1 |
                                     E3 >= 1)]) %>% 
    summarise(min_t_cross = min(t_cross)) %>% as.numeric()
  ## Total Hogs infected
  TotHogsInf <- fulldata.df %>% 
    #select(t, iter, Total_inf1, Total_inf2, Total_inf3, Total_inf4) %>% 
    mutate(Tot_inf = Total_inf1 + Total_inf2 + Total_inf3 + Total_inf4,
           Tot_rec = Total_Rec1 + Total_Rec2 + Total_Rec3 + Total_Rec4) %>% 
    summarise(Tot_inf2 = max(Tot_inf),
              Tot_rec2 = max(Tot_rec))
  ## Minimum Time to Peak Infection
  Time2PeakInf <- fulldata.df %>% 
    mutate(S_tot = S1+S2+S3+S4,                             
           E_tot = E1+E2+E3+E4,
           I_tot = I1+I2+I3+I4,
           R_tot = R1+R2+R3+R4) %>%                                     
    summarise(max_I = max(I_tot),
              t_peak = time[which(I_tot == max(I_tot))]) %>% 
    summarise(min_t_peak = min(t_peak)) %>% as.numeric()
  ##Ro
  Ro <- fulldata.df %>% 
    arrange(time) %>% 
    mutate(Tot_rec = Total_Rec1 + Total_Rec2 + Total_Rec3 + Total_Rec4) %>%  
    filter(row_number() == n()) %>% 
    select(Tot_rec) %>% 
    mutate(Ro = ((log((4000 - Tot_rec)/4000) - log((4000 - (Hog_Vacc_Eff*4000))/4000)) /
                   ((((4000 - Tot_rec) / 4000) - (4000 - (Hog_Vacc_Eff*4000)) / 4000)))) %>% 
    select(Ro) %>% as.numeric()
  
  #Create Summary data frame
  fulldata.df.sum <- data.frame(iter = as.factor(i),
                                HI_inf = HI_inf,
                                Time2WFinf0 = Time2WFinf0,
                                Time2WFinf.5 = Time2WFinf.5,
                                Timeto2room = Timeto2room,
                                TotHogsInf = TotHogsInf$Tot_inf2,
                                TotHogsRec = TotHogsInf$Tot_rec2,
                                Time2PeakInf = Time2PeakInf,
                                Ro = Ro,
                                beta_Hog = beta_Hog,
                                beta_Hog_ind = beta_Hog_ind,
                                beta_S2WF = beta_S2WF,
                                beta_WF2S = beta_WF2S,
                                sigma_Hog = sigma_Hog,
                                delta_Hog = delta_Hog,
                                beta_WF2WF = beta_WF2WF,
                                sigma_WF = sigma_WF,
                                delta_WF = delta_WF)
  
  results <- rbind(results, fulldata.df.sum)
  
  rm(list=ls()[! ls() %in% c("results", "n.itr")])  #Deletes everything from environment expect growing results dataframe
}

saveRDS(results, file = "Current code/Control Meausres/Sensativity Analysis/WFInfIntro_Sum_SensAnlaysis.rds")

#### SUMMARIZE ####
results <- readRDS(file = "Current code/Control Meausres/Sensativity Analysis/WFFlowRoom1_Sum_SensAnlaysis.rds")

##Iterations with interspecies transmssion
results %>% 
  filter(HI_inf > 0) %>% 
  nrow() / 5000 * 100

##time to first intersepcies infection
results %>% 
  summarize(Med_Time2WFInf0 = round(median(Time2WFinf0, na.rm = T), digits = 2),
            LB_Time2WFInf0 = round(quantile(Time2WFinf0, probs = 0.05, na.rm = T), digits = 2),
            UB_Time2WFInf0 = round(quantile(Time2WFinf0, probs = 0.95, na.rm = T), digits = 2),
            Med_Time2WFinf.5 = round(median(Time2WFinf.5, na.rm = T), digits = 2),
            LB_Time2WFinf.5 = round(quantile(Time2WFinf.5, probs = 0.05, na.rm = T), digits = 2),
            UB_Time2WFinf.5 = round(quantile(Time2WFinf.5, probs = 0.95, na.rm = T), digits = 2))

## total hogs infected
results %>% 
  summarize(Med_TotInf = round(median(TotHogsInf, na.rm = T), digits = 2),
            LB_TotInf = round(quantile(TotHogsInf, probs = 0.05, na.rm = T), digits = 2),
            UB_TotInf = round(quantile(TotHogsInf, probs = 0.95, na.rm = T), digits = 2),
            Med_TotRec = round(median(TotHogsRec, na.rm = T), digits = 2),
            LB_TotRec = round(quantile(TotHogsRec, probs = 0.05, na.rm = T), digits = 2),
            UB_TotRec = round(quantile(TotHogsRec, probs = 0.95, na.rm = T), digits = 2))

## Minimum time to peak infection
results %>% 
  summarise(Med_Time2PeakInf = round(median(Time2PeakInf, na.rm = T), digits = 2),
            LB_Time2PeakInf = round(quantile(Time2PeakInf, probs = 0.05, na.rm = T), digits = 2),
            UB_Time2PeakInf = round(quantile(Time2PeakInf, probs = 0.95, na.rm = T), digits = 2))

## R0
results %>% 
  summarize(Med_Ro = round(median(Ro, na.rm = T), digits = 2),
            LB_Ro = round(quantile(Ro, probs = 0.05, na.rm = T), digits = 2),
            UB_Ro = round(quantile(Ro, probs = 0.95, na.rm = T), digits = 2))
## Beta value
results %>% 
  summarize(Med_BetaHog = median(beta_Hog, na.rm = T),
            LB_BetaHog = quantile(beta_Hog, probs = 0.05, na.rm = T),
            UB_BetaHog = quantile(beta_Hog, probs = 0.95, na.rm = T),
            Med_BetaS2WF = median(beta_S2WF, na.rm = T),
            LB_BetaS2WF = quantile(beta_S2WF, probs = 0.05, na.rm = T),
            UB_BetaS2WF = quantile(beta_S2WF, probs = 0.95, na.rm = T),
            Med_BetaWF2S = median(beta_WF2S, na.rm = T),
            LB_BetaWF2S = quantile(beta_WF2S, probs = 0.05, na.rm = T),
            UB_BetaWF2S = quantile(beta_WF2S, probs = 0.95, na.rm = T))

results %>% 
  summarize(Avg_BetaHog = mean(beta_Hog, na.rm = T),
            LB_BetaHog = quantile(beta_Hog, probs = 0.05, na.rm = T),
            UB_BetaHog = quantile(beta_Hog, probs = 0.95, na.rm = T),
            Avg_BetaS2WF = mean(beta_S2WF, na.rm = T),
            LB_BetaS2WF = quantile(beta_S2WF, probs = 0.05, na.rm = T),
            UB_BetaS2WF = quantile(beta_S2WF, probs = 0.95, na.rm = T),
            Avg_BetaWF2S = mean(beta_WF2S, na.rm = T),
            LB_BetaWF2S = quantile(beta_WF2S, probs = 0.05, na.rm = T),
            UB_BetaWF2S = quantile(beta_WF2S, probs = 0.95, na.rm = T))




