#-------------------------------------------------------#
#                                                       #
# This program estimates first- and second-stage        #
# regressions and simulates the outcomes for the        #
# 2018-drought, for the model with heterogeneity.       #
#                                                       #
# Note: The farm-level data are confidential and cannot # 
# be loaded (line 18). The results at the sample mean   #
# are saved in lines 558, 722, and 805, and the figures #
# in Section 5 can be created.                          #
#                                                       #
#-------------------------------------------------------#

library(dplyr)
library(systemfit)
library(ggplot2)
library(ggpubr)
load("rOutput/farm_ready.Rda")

#-----------------------------------------------------------#
#### 1) Prepare the weather variables for the simulation ####
#-----------------------------------------------------------#

# Load simulated weather

  load("rOutput/sim2018.Rda")

# Load LTA weather and weather in 2018

  load("rOutput/df_weather_lta.Rda")
  load("rOutput/df_weather_2018.Rda")

# Add observed weather to simulated weather

sim2018$gdd_obs <- c(df_weather_lta$gdd_1030, df_weather_2018$gdd_1030, rep(df_weather_lta$gdd_1030,10))
sim2018$prec_obs <- c(df_weather_lta$prec, df_weather_2018$prec, rep(df_weather_lta$prec,10))
sim2018$gddHigh_obs <- c(df_weather_lta$gdd_30, df_weather_2018$gdd_30, rep(df_weather_lta$gdd_30,10))
sim2018$dd_obs <- c(df_weather_lta$dd, df_weather_2018$dd, rep(df_weather_lta$dd,10))

# Add interactions between observed and experienced weather

sim2018$gdd_obs_1to3     <- sim2018$gdd_obs * sim2018$gdd_1to3
sim2018$prec_obs_1to3    <- sim2018$prec_obs * sim2018$prec_1to3
sim2018$gddHigh_obs_1to3 <- sim2018$gddHigh_obs * sim2018$gddHigh_1to3
sim2018$dd_obs_1to3          <- sim2018$dd_obs * sim2018$dd_1to3

sim2018$gdd_obs_4to10     <- sim2018$gdd_obs * sim2018$gdd_4to10
sim2018$prec_obs_4to10    <- sim2018$prec_obs * sim2018$prec_4to10
sim2018$gddHigh_obs_4to10 <- sim2018$gddHigh_obs * sim2018$gddHigh_4to10
sim2018$dd_obs_4to10          <- sim2018$dd_obs * sim2018$dd_4to10

# Add time indicator for simulation

sim2018$t <- c(-1:10)

#----------------------------------#
#### 2) Run probit regressions #####
#----------------------------------#

#Define probit regressions
cereal.eq <- iCereals ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
  k_land + k_labor + k_capital + trend + trend2 + 
  lsh_cereals + lsh_protein + lsh_oilseed + lsh_roots + 
  gdd_1to3 + prec_1to3 + gddHigh_1to3 + dd_1to3 +
  gdd_4to10 + prec_4to10 + gddHigh_4to10 + dd_4to10 +
  gdd_1to3_land + prec_1to3_land + gddHigh_1to3_land + dd_1to3_land + 
  gdd_4to10_land + prec_4to10_land + gddHigh_4to10_land + dd_4to10_land + 
  np_cereals_fm + np_protein_fm + np_oilseed_fm + np_roots_fm + np_corn_fm + nw_fert_fm + 
  k_land_fm + k_labor_fm + k_capital_fm + trend_fm + trend2_fm + 
  lsh_cereals_fm + lsh_protein_fm + lsh_oilseed_fm + lsh_roots_fm + 
  gdd_1to3_fm + prec_1to3_fm + gddHigh_1to3_fm + dd_1to3_fm +
  gdd_4to10_fm + prec_4to10_fm + gddHigh_4to10_fm + dd_4to10_fm  +
  gdd_1to3_land_fm + prec_1to3_land_fm + gddHigh_1to3_land_fm + dd_1to3_land_fm + 
  gdd_4to10_land_fm + prec_4to10_land_fm + gddHigh_4to10_land_fm + dd_4to10_land_fm 

protein.eq <- iProtein ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
  k_land + k_labor + k_capital + trend + trend2 + 
  lsh_cereals + lsh_protein + lsh_oilseed + lsh_roots + 
  gdd_1to3 + prec_1to3 + gddHigh_1to3 + dd_1to3 +
  gdd_4to10 + prec_4to10 + gddHigh_4to10 + dd_4to10 +
  gdd_1to3_land + prec_1to3_land + gddHigh_1to3_land + dd_1to3_land +
  gdd_4to10_land + prec_4to10_land + gddHigh_4to10_land + dd_4to10_land +
  np_cereals_fm + np_protein_fm + np_oilseed_fm + np_roots_fm + np_corn_fm + nw_fert_fm + 
  k_land_fm + k_labor_fm + k_capital_fm + trend_fm + trend2_fm + 
  lsh_cereals_fm + lsh_protein_fm + lsh_oilseed_fm + lsh_roots_fm + 
  gdd_1to3_fm + prec_1to3_fm + gddHigh_1to3_fm + dd_1to3_fm +
  gdd_4to10_fm + prec_4to10_fm + gddHigh_4to10_fm + dd_4to10_fm  +
  gdd_1to3_land_fm + prec_1to3_land_fm + gddHigh_1to3_land_fm + dd_1to3_land_fm +
  gdd_4to10_land_fm + prec_4to10_land_fm + gddHigh_4to10_land_fm + dd_4to10_land_fm

oilseed.eq <- iOilseed ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
  k_land + k_labor + k_capital + trend + trend2 + 
  lsh_cereals + lsh_protein + lsh_oilseed + lsh_roots + 
  gdd_1to3 + prec_1to3 + gddHigh_1to3 + dd_1to3 +
  gdd_4to10 + prec_4to10 + gddHigh_4to10 + dd_4to10 +
  gdd_1to3_land + prec_1to3_land + gddHigh_1to3_land + dd_1to3_land +
  gdd_4to10_land + prec_4to10_land + gddHigh_4to10_land + dd_4to10_land +
  np_cereals_fm + np_protein_fm + np_oilseed_fm + np_roots_fm + np_corn_fm + nw_fert_fm + 
  k_land_fm + k_labor_fm + k_capital_fm + trend_fm + trend2_fm + 
  lsh_cereals_fm + lsh_protein_fm + lsh_oilseed_fm + lsh_roots_fm + 
  gdd_1to3_fm + prec_1to3_fm + gddHigh_1to3_fm + dd_1to3_fm +
  gdd_4to10_fm + prec_4to10_fm + gddHigh_4to10_fm + dd_4to10_fm  +
  gdd_1to3_land_fm + prec_1to3_land_fm + gddHigh_1to3_land_fm + dd_1to3_land_fm +
  gdd_4to10_land_fm + prec_4to10_land_fm + gddHigh_4to10_land_fm + dd_4to10_land_fm

roots.eq <- iRoots ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
  k_land + k_labor + k_capital + trend + trend2 + 
  lsh_cereals + lsh_protein + lsh_oilseed + lsh_roots + 
  gdd_1to3 + prec_1to3 + gddHigh_1to3 + dd_1to3 +
  gdd_4to10 + prec_4to10 + gddHigh_4to10 + dd_4to10 +
  gdd_1to3_land + prec_1to3_land + gddHigh_1to3_land + dd_1to3_land +
  gdd_4to10_land + prec_4to10_land + gddHigh_4to10_land + dd_4to10_land +
  np_cereals_fm + np_protein_fm + np_oilseed_fm + np_roots_fm + np_corn_fm + nw_fert_fm + 
  k_land_fm + k_labor_fm + k_capital_fm + trend_fm + trend2_fm + 
  lsh_cereals_fm + lsh_protein_fm + lsh_oilseed_fm + lsh_roots_fm + 
  gdd_1to3_fm + prec_1to3_fm + gddHigh_1to3_fm + dd_1to3_fm +
  gdd_4to10_fm + prec_4to10_fm + gddHigh_4to10_fm + dd_4to10_fm  +
  gdd_1to3_land_fm + prec_1to3_land_fm + gddHigh_1to3_land_fm + dd_1to3_land_fm +
  gdd_4to10_land_fm + prec_4to10_land_fm + gddHigh_4to10_land_fm + dd_4to10_land_fm

corn.eq <- iCorn ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
  k_land + k_labor + k_capital + trend + trend2 + 
  lsh_cereals + lsh_protein + lsh_oilseed + lsh_roots + 
  gdd_1to3 + prec_1to3 + gddHigh_1to3 + dd_1to3 +
  gdd_4to10 + prec_4to10 + gddHigh_4to10 + dd_4to10 +
  gdd_1to3_land + prec_1to3_land + gddHigh_1to3_land + dd_1to3_land +
  gdd_4to10_land + prec_4to10_land + gddHigh_4to10_land + dd_4to10_land +
  np_cereals_fm + np_protein_fm + np_oilseed_fm + np_roots_fm + np_corn_fm + nw_fert_fm + 
  k_land_fm + k_labor_fm + k_capital_fm + trend_fm + trend2_fm + 
  lsh_cereals_fm + lsh_protein_fm + lsh_oilseed_fm + lsh_roots_fm + 
  gdd_1to3_fm + prec_1to3_fm + gddHigh_1to3_fm + dd_1to3_fm +
  gdd_4to10_fm + prec_4to10_fm + gddHigh_4to10_fm + dd_4to10_fm +
  gdd_1to3_land_fm + prec_1to3_land_fm + gddHigh_1to3_land_fm + dd_1to3_land_fm +
  gdd_4to10_land_fm + prec_4to10_land_fm + gddHigh_4to10_land_fm + dd_4to10_land_fm


#set-up the loop over crops
list_crops <- data.frame(crops=c("Cereals","Protein","Oilseed","Roots","Corn"),
                         probit=c("cereal.eq","protein.eq","oilseed.eq","roots.eq","corn.eq"),
                         abbrev=c("cer","prot","oil","roots","corn"))
n_crops <- list_crops %>% count() %>% as.numeric() # number of crops

#set-up the lists to store the results

  # List for actual results
  list_probit_act <- list(Cereals=list(PHI="", phi=""),
                          Protein=list(PHI="", phi=""),
                          Oilseed=list(PHI="", phi=""),
                          Roots=  list(PHI="", phi=""),
                          Corn=list(PHI="", phi=""))
  
  # List for simulated results (Note: I need t=1 to t=12 because t-1 does not work in lists)
  list_probit_sim <- list(Cereals=list(PHI=list(t1="", t2="", t3="", t4="", t5="",t6="", t7="", t8="", t9="", t10="", t11="", t12=""), 
                                       phi=list(t1="", t2="", t3="", t4="", t5="",t6="", t7="", t8="", t9="", t10="", t11="", t12="")),
                          Protein=list(PHI=list(t1="", t2="", t3="", t4="", t5="",t6="", t7="", t8="", t9="", t10="", t11="", t12=""), 
                                       phi=list(t1="", t2="", t3="", t4="", t5="",t6="", t7="", t8="", t9="", t10="", t11="", t12="")),
                          Oilseed=list(PHI=list(t1="", t2="", t3="", t4="", t5="",t6="", t7="", t8="", t9="", t10="", t11="", t12=""), 
                                       phi=list(t1="", t2="", t3="", t4="", t5="",t6="", t7="", t8="", t9="", t10="", t11="", t12="")),
                          Roots=  list(PHI=list(t1="", t2="", t3="", t4="", t5="",t6="", t7="", t8="", t9="", t10="", t11="", t12=""), 
                                       phi=list(t1="", t2="", t3="", t4="", t5="",t6="", t7="", t8="", t9="", t10="", t11="", t12="")),
                          Corn=list(PHI=list(t1="", t2="", t3="", t4="", t5="",t6="", t7="", t8="", t9="", t10="", t11="", t12=""), 
                                       phi=list(t1="", t2="", t3="", t4="", t5="",t6="", t7="", t8="", t9="", t10="", t11="", t12="")))

# Loop over all crops and estimate probits 
for (i in 1:n_crops) { 
  
  # Probit estimation
  crop.eq <- get(list_crops[i,"probit"])
  iCrop <- glm(crop.eq, family = binomial(link = "probit"), 
               data = df_farm)
  
  # Actual prediction
  
      # At the farm-level
      iCrop_linpred <- predict(iCrop) 
      
      # PHI
      list_probit_act[[i]][[1]] <- pnorm(iCrop_linpred)
      
      # phi
      list_probit_act[[i]][[2]] <- dnorm(iCrop_linpred)
  
  # Simulation
  
      # Data used in probit model
      X <- as_tibble(model.matrix(iCrop)) #note: with as_tibble, I retain the column names
      
      for (t in -1:10)
      {
        
        # Simulate t = -1
        X$gdd_1to3 <- sim2018$gdd_1to3[sim2018$t==t]
        X$prec_1to3 <- sim2018$prec_1to3[sim2018$t==t]
        X$gddHigh_1to3 <- sim2018$gddHigh_1to3[sim2018$t==t]
        X$dd_1to3 <- sim2018$dd_1to3[sim2018$t==t]
        X$gdd_4to10 <- sim2018$gdd_4to10[sim2018$t==t]
        X$prec_4to10 <- sim2018$prec_4to10[sim2018$t==t]
        X$gddHigh_4to10 <- sim2018$gddHigh_4to10[sim2018$t==t]
        X$dd_4to10  <- sim2018$dd_4to10[sim2018$t==t]
        
        X$gdd_1to3_land <- sim2018$gdd_1to3[sim2018$t==t] * X$k_land
        X$prec_1to3_land <- sim2018$prec_1to3[sim2018$t==t] * X$k_land
        X$gddHigh_1to3_land <- sim2018$gddHigh_1to3[sim2018$t==t] * X$k_land
        X$dd_1to3_land <- sim2018$dd_1to3[sim2018$t==t] * X$k_land
        X$gdd_4to10_land <- sim2018$gdd_4to10[sim2018$t==t] * X$k_land
        X$prec_4to10_land <- sim2018$prec_4to10[sim2018$t==t] * X$k_land
        X$gddHigh_4to10_land <- sim2018$gddHigh_4to10[sim2018$t==t] * X$k_land
        X$dd_4to10_land <- sim2018$dd_4to10[sim2018$t==t] * X$k_land
        
        # Store PHI's and phi's
        list_probit_sim[[i]][[1]][[t+2]] <- as.vector ( pnorm( t( as.matrix(as.numeric(iCrop$coefficients) %*%  t(as.matrix(X))) ) ) )
        list_probit_sim[[i]][[2]][[t+2]] <- as.vector ( dnorm( t( as.matrix(as.numeric(iCrop$coefficients) %*%  t(as.matrix(X))) ) ) )
    
  }
}

# Append actual farm-level PHIs and phis to data set
df_farm$PHI_cereals <- as.numeric(do.call(c, list_probit_act[[1]][1]))
df_farm$PHI_protein <- as.numeric(do.call(c, list_probit_act[[2]][1]))
df_farm$PHI_oilseed <- as.numeric(do.call(c, list_probit_act[[3]][1]))
df_farm$PHI_roots <-   as.numeric(do.call(c, list_probit_act[[4]][1]))
df_farm$PHI_corn <- as.numeric(do.call(c, list_probit_act[[5]][1]))

df_farm$phi_cereals <- as.numeric(do.call(c, list_probit_act[[1]][2]))
df_farm$phi_protein <- as.numeric(do.call(c, list_probit_act[[2]][2]))
df_farm$phi_oilseed <- as.numeric(do.call(c, list_probit_act[[3]][2]))
df_farm$phi_roots <-   as.numeric(do.call(c, list_probit_act[[4]][2]))
df_farm$phi_corn <- as.numeric(do.call(c, list_probit_act[[5]][2]))

#--------------------------------------#
#### 3) Run structural regressions #####
#--------------------------------------#

#Define structural equations
eqQQcereals <- qq_cereals ~ PHI_cereals + I(PHI_cereals*np_cereals) + I(PHI_cereals*np_protein) + I(PHI_cereals*np_oilseed) + I(PHI_cereals*np_roots) + I(PHI_cereals*np_corn) + I(PHI_cereals*nw_fert) + 
  I(PHI_cereals*k_land) + I(PHI_cereals*k_labor) + I(PHI_cereals*k_capital) + I(PHI_cereals*trend) + I(PHI_cereals*trend2) +
  I(PHI_cereals*gdd_obs) + I(PHI_cereals*prec_obs) + I(PHI_cereals*gddHigh_obs) + I(PHI_cereals*dd_obs) + 
  I(PHI_cereals*gdd_1to3) + I(PHI_cereals*prec_1to3) + I(PHI_cereals*gddHigh_1to3) + I(PHI_cereals*dd_1to3) +
  I(PHI_cereals*gdd_4to10) + I(PHI_cereals*prec_4to10) + I(PHI_cereals*gddHigh_4to10) + I(PHI_cereals*dd_4to10) +
  I(PHI_cereals*gdd_obs_land) + I(PHI_cereals*prec_obs_land) + I(PHI_cereals*gddHigh_obs_land) + I(PHI_cereals*dd_obs_land) + # Compared to the baseline model, weather vars are interacted with land
  I(PHI_cereals*gdd_1to3_land) + I(PHI_cereals*prec_1to3_land) + I(PHI_cereals*gddHigh_1to3_land) + I(PHI_cereals*dd_1to3_land) + # Compared to the baseline model, weather vars are interacted with land
  I(PHI_cereals*gdd_4to10_land) + I(PHI_cereals*prec_4to10_land) + I(PHI_cereals*gddHigh_4to10_land) + I(PHI_cereals*dd_4to10_land) + # Compared to the baseline model, weather vars are interacted with land
  I(PHI_cereals*gdd_obs_1to3) + I(PHI_cereals*prec_obs_1to3) + I(PHI_cereals*gddHigh_obs_1to3) + I(PHI_cereals*dd_obs_1to3) + # Weather-climate Interactions as in the baseline model
  I(PHI_cereals*gdd_obs_4to10) + I(PHI_cereals*prec_obs_4to10) + I(PHI_cereals*gddHigh_obs_4to10) + I(PHI_cereals*dd_obs_4to10) + # Weather-climate Interactions as in the baseline model
  I(PHI_cereals*np_cereals_fm) + I(PHI_cereals*np_protein_fm) + I(PHI_cereals*np_oilseed_fm) + I(PHI_cereals*np_roots_fm) + I(PHI_cereals*np_corn_fm) + I(PHI_cereals*nw_fert_fm) +
  I(PHI_cereals*k_land_fm) + I(PHI_cereals*k_labor_fm) + I(PHI_cereals*k_capital_fm) + I(PHI_cereals*trend_fm) + I(PHI_cereals*trend2_fm) + 
  I(PHI_cereals*gdd_obs_fm) + I(PHI_cereals*prec_obs_fm) + I(PHI_cereals*gddHigh_obs_fm) + I(PHI_cereals*dd_obs_fm) + 
  I(PHI_cereals*gdd_1to3_fm) + I(PHI_cereals*prec_1to3_fm) + I(PHI_cereals*gddHigh_1to3_fm) + I(PHI_cereals*dd_1to3_fm) +
  I(PHI_cereals*gdd_4to10_fm) + I(PHI_cereals*prec_4to10_fm) + I(PHI_cereals*gddHigh_4to10_fm) + I(PHI_cereals*dd_4to10_fm) + 
  I(PHI_cereals*gdd_obs_land_fm) + I(PHI_cereals*prec_obs_land_fm) + I(PHI_cereals*gddHigh_obs_land_fm) + I(PHI_cereals*dd_obs_land_fm) + # Of course, the interactions are also here
  I(PHI_cereals*gdd_1to3_land_fm) + I(PHI_cereals*prec_1to3_land_fm) + I(PHI_cereals*gddHigh_1to3_land_fm) + I(PHI_cereals*dd_1to3_land_fm) + # Of course, the interactions are also here
  I(PHI_cereals*gdd_4to10_land_fm) + I(PHI_cereals*prec_4to10_land_fm) + I(PHI_cereals*gddHigh_4to10_land_fm) + I(PHI_cereals*dd_4to10_land_fm) + # Of course, the interactions are also here
  I(PHI_cereals*gdd_obs_1to3_fm) + I(PHI_cereals*prec_obs_1to3_fm) + I(PHI_cereals*gddHigh_obs_1to3_fm) + I(PHI_cereals*dd_obs_1to3_fm) + # Weather-climate Interactions as in the baseline model
  I(PHI_cereals*gdd_obs_4to10_fm) + I(PHI_cereals*prec_obs_4to10_fm) + I(PHI_cereals*gddHigh_obs_4to10_fm) + I(PHI_cereals*dd_obs_4to10_fm) + # Weather-climate Interactions as in the baseline model
  phi_cereals -1
eqQQprotein <- qq_protein ~ PHI_protein + I(PHI_protein*np_cereals) + I(PHI_protein*np_protein) + I(PHI_protein*np_oilseed) + I(PHI_protein*np_roots) + I(PHI_protein*np_corn) + I(PHI_protein*nw_fert) + 
  I(PHI_protein*k_land) + I(PHI_protein*k_labor) + I(PHI_protein*k_capital) + I(PHI_protein*trend) + I(PHI_protein*trend2) +
  I(PHI_protein*gdd_obs) + I(PHI_protein*prec_obs) + I(PHI_protein*gddHigh_obs) + I(PHI_protein*dd_obs) + 
  I(PHI_protein*gdd_1to3) + I(PHI_protein*prec_1to3) + I(PHI_protein*gddHigh_1to3) + I(PHI_protein*dd_1to3) +
  I(PHI_protein*gdd_4to10) + I(PHI_protein*prec_4to10) + I(PHI_protein*gddHigh_4to10) + I(PHI_protein*dd_4to10) +
  I(PHI_protein*gdd_obs_land) + I(PHI_protein*prec_obs_land) + I(PHI_protein*gddHigh_obs_land) + I(PHI_protein*dd_obs_land) + 
  I(PHI_protein*gdd_1to3_land) + I(PHI_protein*prec_1to3_land) + I(PHI_protein*gddHigh_1to3_land) + I(PHI_protein*dd_1to3_land) +
  I(PHI_protein*gdd_4to10_land) + I(PHI_protein*prec_4to10_land) + I(PHI_protein*gddHigh_4to10_land) + I(PHI_protein*dd_4to10_land) +
  I(PHI_protein*gdd_obs_1to3) + I(PHI_protein*prec_obs_1to3) + I(PHI_protein*gddHigh_obs_1to3) + I(PHI_protein*dd_obs_1to3) + 
  I(PHI_protein*gdd_obs_4to10) + I(PHI_protein*prec_obs_4to10) + I(PHI_protein*gddHigh_obs_4to10) + I(PHI_protein*dd_obs_4to10) + 
  I(PHI_protein*np_cereals_fm) + I(PHI_protein*np_protein_fm) + I(PHI_protein*np_oilseed_fm) + I(PHI_protein*np_roots_fm) + I(PHI_protein*np_corn_fm) + I(PHI_protein*nw_fert_fm) +
  I(PHI_protein*k_land_fm) + I(PHI_protein*k_labor_fm) + I(PHI_protein*k_capital_fm) + I(PHI_protein*trend_fm) + I(PHI_protein*trend2_fm) + 
  I(PHI_protein*gdd_obs_fm) + I(PHI_protein*prec_obs_fm) + I(PHI_protein*gddHigh_obs_fm) + I(PHI_protein*dd_obs_fm) + 
  I(PHI_protein*gdd_1to3_fm) + I(PHI_protein*prec_1to3_fm) + I(PHI_protein*gddHigh_1to3_fm) + I(PHI_protein*dd_1to3_fm) +
  I(PHI_protein*gdd_4to10_fm) + I(PHI_protein*prec_4to10_fm) + I(PHI_protein*gddHigh_4to10_fm) + I(PHI_protein*dd_4to10_fm) + 
  I(PHI_protein*gdd_obs_land_fm) + I(PHI_protein*prec_obs_land_fm) + I(PHI_protein*gddHigh_obs_land_fm) + I(PHI_protein*dd_obs_land_fm) + 
  I(PHI_protein*gdd_1to3_land_fm) + I(PHI_protein*prec_1to3_land_fm) + I(PHI_protein*gddHigh_1to3_land_fm) + I(PHI_protein*dd_1to3_land_fm) +
  I(PHI_protein*gdd_4to10_land_fm) + I(PHI_protein*prec_4to10_land_fm) + I(PHI_protein*gddHigh_4to10_land_fm) + I(PHI_protein*dd_4to10_land_fm) +
  I(PHI_protein*gdd_obs_1to3_fm) + I(PHI_protein*prec_obs_1to3_fm) + I(PHI_protein*gddHigh_obs_1to3_fm) + I(PHI_protein*dd_obs_1to3_fm) + 
  I(PHI_protein*gdd_obs_4to10_fm) + I(PHI_protein*prec_obs_4to10_fm) + I(PHI_protein*gddHigh_obs_4to10_fm) + I(PHI_protein*dd_obs_4to10_fm) + 
  phi_protein -1
eqQQoilseed <- qq_oilseed ~ PHI_oilseed + I(PHI_oilseed*np_cereals) + I(PHI_oilseed*np_protein) + I(PHI_oilseed*np_oilseed) + I(PHI_oilseed*np_roots) + I(PHI_oilseed*np_corn) + I(PHI_oilseed*nw_fert) + 
  I(PHI_oilseed*k_land) + I(PHI_oilseed*k_labor) + I(PHI_oilseed*k_capital) + I(PHI_oilseed*trend) + I(PHI_oilseed*trend2) +
  I(PHI_oilseed*gdd_obs) + I(PHI_oilseed*prec_obs) + I(PHI_oilseed*gddHigh_obs) + I(PHI_oilseed*dd_obs) + 
  I(PHI_oilseed*gdd_1to3) + I(PHI_oilseed*prec_1to3) + I(PHI_oilseed*gddHigh_1to3) + I(PHI_oilseed*dd_1to3) +
  I(PHI_oilseed*gdd_4to10) + I(PHI_oilseed*prec_4to10) + I(PHI_oilseed*gddHigh_4to10) + I(PHI_oilseed*dd_4to10) +
  I(PHI_oilseed*gdd_obs_land) + I(PHI_oilseed*prec_obs_land) + I(PHI_oilseed*gddHigh_obs_land) + I(PHI_oilseed*dd_obs_land) + 
  I(PHI_oilseed*gdd_1to3_land) + I(PHI_oilseed*prec_1to3_land) + I(PHI_oilseed*gddHigh_1to3_land) + I(PHI_oilseed*dd_1to3_land) +
  I(PHI_oilseed*gdd_4to10_land) + I(PHI_oilseed*prec_4to10_land) + I(PHI_oilseed*gddHigh_4to10_land) + I(PHI_oilseed*dd_4to10_land) +
  I(PHI_oilseed*gdd_obs_1to3) + I(PHI_oilseed*prec_obs_1to3) + I(PHI_oilseed*gddHigh_obs_1to3) + I(PHI_oilseed*dd_obs_1to3) + 
  I(PHI_oilseed*gdd_obs_4to10) + I(PHI_oilseed*prec_obs_4to10) + I(PHI_oilseed*gddHigh_obs_4to10) + I(PHI_oilseed*dd_obs_4to10) + 
  I(PHI_oilseed*np_cereals_fm) + I(PHI_oilseed*np_protein_fm) + I(PHI_oilseed*np_oilseed_fm) + I(PHI_oilseed*np_roots_fm) + I(PHI_oilseed*np_corn_fm) + I(PHI_oilseed*nw_fert_fm) +
  I(PHI_oilseed*k_land_fm) + I(PHI_oilseed*k_labor_fm) + I(PHI_oilseed*k_capital_fm) + I(PHI_oilseed*trend_fm) + I(PHI_oilseed*trend2_fm) + 
  I(PHI_oilseed*gdd_obs_fm) + I(PHI_oilseed*prec_obs_fm) + I(PHI_oilseed*gddHigh_obs_fm) + I(PHI_oilseed*dd_obs_fm) + 
  I(PHI_oilseed*gdd_1to3_fm) + I(PHI_oilseed*prec_1to3_fm) + I(PHI_oilseed*gddHigh_1to3_fm) + I(PHI_oilseed*dd_1to3_fm) +
  I(PHI_oilseed*gdd_4to10_fm) + I(PHI_oilseed*prec_4to10_fm) + I(PHI_oilseed*gddHigh_4to10_fm) + I(PHI_oilseed*dd_4to10_fm) + 
  I(PHI_oilseed*gdd_obs_land_fm) + I(PHI_oilseed*prec_obs_land_fm) + I(PHI_oilseed*gddHigh_obs_land_fm) + I(PHI_oilseed*dd_obs_land_fm) + 
  I(PHI_oilseed*gdd_1to3_land_fm) + I(PHI_oilseed*prec_1to3_land_fm) + I(PHI_oilseed*gddHigh_1to3_land_fm) + I(PHI_oilseed*dd_1to3_land_fm) +
  I(PHI_oilseed*gdd_4to10_land_fm) + I(PHI_oilseed*prec_4to10_land_fm) + I(PHI_oilseed*gddHigh_4to10_land_fm) + I(PHI_oilseed*dd_4to10_land_fm) +
  I(PHI_oilseed*gdd_obs_1to3_fm) + I(PHI_oilseed*prec_obs_1to3_fm) + I(PHI_oilseed*gddHigh_obs_1to3_fm) + I(PHI_oilseed*dd_obs_1to3_fm) + 
  I(PHI_oilseed*gdd_obs_4to10_fm) + I(PHI_oilseed*prec_obs_4to10_fm) + I(PHI_oilseed*gddHigh_obs_4to10_fm) + I(PHI_oilseed*dd_obs_4to10_fm) + 
  phi_oilseed -1
eqQQroots <- qq_roots ~ PHI_roots + I(PHI_roots*np_cereals) + I(PHI_roots*np_protein) + I(PHI_roots*np_oilseed) + I(PHI_roots*np_roots) + I(PHI_roots*np_corn) + I(PHI_roots*nw_fert) + 
  I(PHI_roots*k_land) + I(PHI_roots*k_labor) + I(PHI_roots*k_capital) + I(PHI_roots*trend) + I(PHI_roots*trend2) +
  I(PHI_roots*gdd_obs) + I(PHI_roots*prec_obs) + I(PHI_roots*gddHigh_obs) + I(PHI_roots*dd_obs) + 
  I(PHI_roots*gdd_1to3) + I(PHI_roots*prec_1to3) + I(PHI_roots*gddHigh_1to3) + I(PHI_roots*dd_1to3) +
  I(PHI_roots*gdd_4to10) + I(PHI_roots*prec_4to10) + I(PHI_roots*gddHigh_4to10) + I(PHI_roots*dd_4to10) +
  I(PHI_roots*gdd_obs_land) + I(PHI_roots*prec_obs_land) + I(PHI_roots*gddHigh_obs_land) + I(PHI_roots*dd_obs_land) + 
  I(PHI_roots*gdd_1to3_land) + I(PHI_roots*prec_1to3_land) + I(PHI_roots*gddHigh_1to3_land) + I(PHI_roots*dd_1to3_land) +
  I(PHI_roots*gdd_4to10_land) + I(PHI_roots*prec_4to10_land) + I(PHI_roots*gddHigh_4to10_land) + I(PHI_roots*dd_4to10_land) +
  I(PHI_roots*gdd_obs_1to3) + I(PHI_roots*prec_obs_1to3) + I(PHI_roots*gddHigh_obs_1to3) + I(PHI_roots*dd_obs_1to3) + 
  I(PHI_roots*gdd_obs_4to10) + I(PHI_roots*prec_obs_4to10) + I(PHI_roots*gddHigh_obs_4to10) + I(PHI_roots*dd_obs_4to10) + 
  I(PHI_roots*np_cereals_fm) + I(PHI_roots*np_protein_fm) + I(PHI_roots*np_oilseed_fm) + I(PHI_roots*np_roots_fm) + I(PHI_roots*np_corn_fm) + I(PHI_roots*nw_fert_fm) +
  I(PHI_roots*k_land_fm) + I(PHI_roots*k_labor_fm) + I(PHI_roots*k_capital_fm) + I(PHI_roots*trend_fm) + I(PHI_roots*trend2_fm) + 
  I(PHI_roots*gdd_obs_fm) + I(PHI_roots*prec_obs_fm) + I(PHI_roots*gddHigh_obs_fm) + I(PHI_roots*dd_obs_fm) + 
  I(PHI_roots*gdd_1to3_fm) + I(PHI_roots*prec_1to3_fm) + I(PHI_roots*gddHigh_1to3_fm) + I(PHI_roots*dd_1to3_fm) +
  I(PHI_roots*gdd_4to10_fm) + I(PHI_roots*prec_4to10_fm) + I(PHI_roots*gddHigh_4to10_fm) + I(PHI_roots*dd_4to10_fm) +
  I(PHI_roots*gdd_obs_land_fm) + I(PHI_roots*prec_obs_land_fm) + I(PHI_roots*gddHigh_obs_land_fm) + I(PHI_roots*dd_obs_land_fm) + 
  I(PHI_roots*gdd_1to3_land_fm) + I(PHI_roots*prec_1to3_land_fm) + I(PHI_roots*gddHigh_1to3_land_fm) + I(PHI_roots*dd_1to3_land_fm) +
  I(PHI_roots*gdd_4to10_land_fm) + I(PHI_roots*prec_4to10_land_fm) + I(PHI_roots*gddHigh_4to10_land_fm) + I(PHI_roots*dd_4to10_land_fm) +
  I(PHI_roots*gdd_obs_1to3_fm) + I(PHI_roots*prec_obs_1to3_fm) + I(PHI_roots*gddHigh_obs_1to3_fm) + I(PHI_roots*dd_obs_1to3_fm) + 
  I(PHI_roots*gdd_obs_4to10_fm) + I(PHI_roots*prec_obs_4to10_fm) + I(PHI_roots*gddHigh_obs_4to10_fm) + I(PHI_roots*dd_obs_4to10_fm) + 
  phi_roots -1
eqQQcorn <- qq_corn ~ PHI_corn + I(PHI_corn*np_cereals) + I(PHI_corn*np_protein) + I(PHI_corn*np_oilseed) + I(PHI_corn*np_roots) + I(PHI_corn*np_corn) + I(PHI_corn*nw_fert) + 
  I(PHI_corn*k_land) + I(PHI_corn*k_labor) + I(PHI_corn*k_capital) + I(PHI_corn*trend) + I(PHI_corn*trend2) +
  I(PHI_corn*gdd_obs) + I(PHI_corn*prec_obs) + I(PHI_corn*gddHigh_obs) + I(PHI_corn*dd_obs) + 
  I(PHI_corn*gdd_1to3) + I(PHI_corn*prec_1to3) + I(PHI_corn*gddHigh_1to3) + I(PHI_corn*dd_1to3) +
  I(PHI_corn*gdd_4to10) + I(PHI_corn*prec_4to10) + I(PHI_corn*gddHigh_4to10) + I(PHI_corn*dd_4to10) +
  I(PHI_corn*gdd_obs_land) + I(PHI_corn*prec_obs_land) + I(PHI_corn*gddHigh_obs_land) + I(PHI_corn*dd_obs_land) + 
  I(PHI_corn*gdd_1to3_land) + I(PHI_corn*prec_1to3_land) + I(PHI_corn*gddHigh_1to3_land) + I(PHI_corn*dd_1to3_land) +
  I(PHI_corn*gdd_4to10_land) + I(PHI_corn*prec_4to10_land) + I(PHI_corn*gddHigh_4to10_land) + I(PHI_corn*dd_4to10_land) +
  I(PHI_corn*gdd_obs_1to3) + I(PHI_corn*prec_obs_1to3) + I(PHI_corn*gddHigh_obs_1to3) + I(PHI_corn*dd_obs_1to3) + 
  I(PHI_corn*gdd_obs_4to10) + I(PHI_corn*prec_obs_4to10) + I(PHI_corn*gddHigh_obs_4to10) + I(PHI_corn*dd_obs_4to10) + 
  I(PHI_corn*np_cereals_fm) + I(PHI_corn*np_protein_fm) + I(PHI_corn*np_oilseed_fm) + I(PHI_corn*np_roots_fm) + I(PHI_corn*np_corn_fm) + I(PHI_corn*nw_fert_fm) +
  I(PHI_corn*k_land_fm) + I(PHI_corn*k_labor_fm) + I(PHI_corn*k_capital_fm) + I(PHI_corn*trend_fm) + I(PHI_corn*trend2_fm) + 
  I(PHI_corn*gdd_obs_fm) + I(PHI_corn*prec_obs_fm) + I(PHI_corn*gddHigh_obs_fm) + I(PHI_corn*dd_obs_fm) + 
  I(PHI_corn*gdd_1to3_fm) + I(PHI_corn*prec_1to3_fm) + I(PHI_corn*gddHigh_1to3_fm) + I(PHI_corn*dd_1to3_fm) +
  I(PHI_corn*gdd_4to10_fm) + I(PHI_corn*prec_4to10_fm) + I(PHI_corn*gddHigh_4to10_fm) + I(PHI_corn*dd_4to10_fm) +
  I(PHI_corn*gdd_obs_land_fm) + I(PHI_corn*prec_obs_land_fm) + I(PHI_corn*gddHigh_obs_land_fm) + I(PHI_corn*dd_obs_land_fm) + 
  I(PHI_corn*gdd_1to3_land_fm) + I(PHI_corn*prec_1to3_land_fm) + I(PHI_corn*gddHigh_1to3_land_fm) + I(PHI_corn*dd_1to3_land_fm) +
  I(PHI_corn*gdd_4to10_land_fm) + I(PHI_corn*prec_4to10_land_fm) + I(PHI_corn*gddHigh_4to10_land_fm) + I(PHI_corn*dd_4to10_land_fm) +
  I(PHI_corn*gdd_obs_1to3_fm) + I(PHI_corn*prec_obs_1to3_fm) + I(PHI_corn*gddHigh_obs_1to3_fm) + I(PHI_corn*dd_obs_1to3_fm) + 
  I(PHI_corn*gdd_obs_4to10_fm) + I(PHI_corn*prec_obs_4to10_fm) + I(PHI_corn*gddHigh_obs_4to10_fm) + I(PHI_corn*dd_obs_4to10_fm) + 
  phi_corn -1
eqNXfert <- nx_fert ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
  k_land + k_labor + k_capital + trend + trend2 +
  gdd_obs + prec_obs + gddHigh_obs + dd_obs + 
  gdd_1to3 + prec_1to3 + gddHigh_1to3 + dd_1to3 +
  gdd_4to10 + prec_4to10 + gddHigh_4to10 + dd_4to10 +
  gdd_obs_land + prec_obs_land + gddHigh_obs_land + dd_obs_land + 
  gdd_1to3_land + prec_1to3_land + gddHigh_1to3_land + dd_1to3_land +
  gdd_4to10_land + prec_4to10_land + gddHigh_4to10_land + dd_4to10_land +
  gdd_obs_1to3 + prec_obs_1to3 + gddHigh_obs_1to3 + dd_obs_1to3 + 
  gdd_obs_4to10 + prec_obs_4to10 + gddHigh_obs_4to10 + dd_obs_4to10 + 
  np_cereals_fm + np_protein_fm + np_oilseed_fm + np_roots_fm + np_corn_fm + nw_fert_fm +
  k_land_fm + k_labor_fm + k_capital_fm + trend_fm + trend2_fm + 
  gdd_obs_fm + prec_obs_fm + gddHigh_obs_fm + dd_obs_fm + 
  gdd_1to3_fm + prec_1to3_fm + gddHigh_1to3_fm + dd_1to3_fm +
  gdd_4to10_fm + prec_4to10_fm + gddHigh_4to10_fm + dd_4to10_fm +
  gdd_obs_land_fm + prec_obs_land_fm + gddHigh_obs_land_fm + dd_obs_land_fm + 
  gdd_1to3_land_fm + prec_1to3_land_fm + gddHigh_1to3_land_fm + dd_1to3_land_fm +
  gdd_4to10_land_fm + prec_4to10_land_fm + gddHigh_4to10_land_fm + dd_4to10_land_fm +
  gdd_obs_1to3_fm + prec_obs_1to3_fm + gddHigh_obs_1to3_fm + dd_obs_1to3_fm +
  gdd_obs_4to10_fm + prec_obs_4to10_fm + gddHigh_obs_4to10_fm + dd_obs_4to10_fm

system <- list( QQcereals = eqQQcereals, 
                QQprotein = eqQQprotein,
                QQoilseed = eqQQoilseed,
                QQroots = eqQQroots,
                QQcorn = eqQQcorn,
                NXfert = eqNXfert)

## restrictions
restrict <- c( "QQcereals_I(PHI_cereals * np_protein) - QQprotein_I(PHI_protein * np_cereals) = 0",
               "QQcereals_I(PHI_cereals * np_oilseed) - QQoilseed_I(PHI_oilseed * np_cereals) = 0",
               "QQcereals_I(PHI_cereals * np_roots)   - QQroots_I(PHI_roots * np_cereals) = 0",
               "QQcereals_I(PHI_cereals * np_corn)  - QQcorn_I(PHI_corn * np_cereals) = 0",
               "QQcereals_I(PHI_cereals * nw_fert)    - NXfert_np_cereals = 0",
               "QQprotein_I(PHI_protein * np_oilseed) - QQoilseed_I(PHI_oilseed * np_protein) = 0",
               "QQprotein_I(PHI_protein * np_roots)   - QQroots_I(PHI_roots * np_protein) = 0",
               "QQprotein_I(PHI_protein * np_corn)  - QQcorn_I(PHI_corn * np_protein) = 0",
               "QQprotein_I(PHI_protein * nw_fert)    - NXfert_np_protein = 0",
               "QQoilseed_I(PHI_oilseed * np_roots)   - QQroots_I(PHI_roots * np_oilseed) = 0",
               "QQoilseed_I(PHI_oilseed * np_corn)  - QQcorn_I(PHI_corn * np_oilseed) = 0",
               "QQoilseed_I(PHI_oilseed * nw_fert)    - NXfert_np_oilseed = 0",
               "QQroots_I(PHI_roots * np_corn)      - QQcorn_I(PHI_corn * np_roots) = 0",
               "QQroots_I(PHI_roots * nw_fert)        - NXfert_np_roots = 0",
               "QQcorn_I(PHI_corn * nw_fert)      - NXfert_np_corn = 0")

## Regression with iterated SUR estimation
model_linear <- systemfit( formula = system, method = "SUR", 
                           data = df_farm, restrict.matrix = restrict, 
                           maxit = 100 )

# Extract data used in this estimation for simulation later #

  dat_str_cereals <- as_tibble(model.matrix(model_linear)[grep('^QQcereals', rownames(model.matrix(model_linear))),]) %>% 
    mutate_at(vars(c(1:dim(.)[2])), funs(./df_farm$PHI_cereals)) %>% 
    select(where(~ any(. != 0)))
  
  dat_str_protein <- as_tibble(model.matrix(model_linear)[grep('^QQprotein', rownames(model.matrix(model_linear))),]) %>% 
    mutate_at(vars(c(1:dim(.)[2])), funs(./df_farm$PHI_protein)) %>% 
    select(where(~ any(. != 0)))
  
  dat_str_oilseed <- as_tibble(model.matrix(model_linear)[grep('^QQoilseed', rownames(model.matrix(model_linear))),]) %>% 
    mutate_at(vars(c(1:dim(.)[2])), funs(./df_farm$PHI_oilseed)) %>% 
    select(where(~ any(. != 0)))
  
  dat_str_roots <- as_tibble(model.matrix(model_linear)[grep('^QQroots', rownames(model.matrix(model_linear))),]) %>% 
    mutate_at(vars(c(1:dim(.)[2])), funs(./df_farm$PHI_roots)) %>% 
    select(where(~ any(. != 0)))
  
  dat_str_corn <- as_tibble(model.matrix(model_linear)[grep('^QQcorn', rownames(model.matrix(model_linear))),]) %>% 
    mutate_at(vars(c(1:dim(.)[2])), funs(./df_farm$PHI_corn)) %>% 
    select(where(~ any(. != 0)))
  
  dat_str_fert <- as_tibble(model.matrix(model_linear)[grep('^NXfert', rownames(model.matrix(model_linear))),]) %>% 
    select(where(~ any(. != 0)))

# ------------------------------------ #
#### 4) Simulation of final outcome ####
# ------------------------------------ #

#set-up the loop over crops
list_outcome <- data.frame(outcome=c("cereals","protein","oilseed","roots","corn", "fert"),
                           sur=c("QQcereals","QQprotein","QQoilseed","QQroots","QQcorn", "NXfert"))
n_outcomes <- list_outcome %>% count() %>% as.numeric() # number of crops

#set-up the lists to store the results
list_simresults <- list(Cereals=list(t1="", t2="", t3="", t4="", t5="",t6="", t7="", t8="", t9="", t10="", t11="", t12=""), 
                        Protein=list(t1="", t2="", t3="", t4="", t5="",t6="", t7="", t8="", t9="", t10="", t11="", t12=""), 
                        Oilseed=list(t1="", t2="", t3="", t4="", t5="",t6="", t7="", t8="", t9="", t10="", t11="", t12=""), 
                        Roots=list(t1="", t2="", t3="", t4="", t5="",t6="", t7="", t8="", t9="", t10="", t11="", t12=""), 
                        Corn=list(t1="", t2="", t3="", t4="", t5="",t6="", t7="", t8="", t9="", t10="", t11="", t12=""), 
                        Fertilizer=list(t1="", t2="", t3="", t4="", t5="",t6="", t7="", t8="", t9="", t10="", t11="", t12="")) 

# Loop over all outcomes and estimate probits 
for (i in 1:n_outcomes) { 
  
  # Get data
  X <- as_tibble(get(paste0("dat_str_", list_outcome[i,"outcome"])))
  
  # "Clean" variable names
  names(X) <- c("const","np_cereals","np_protein","np_oilseed","np_roots","np_corn","nw_fert",
                "k_land","k_labor","k_capital","trend","trend2",
                "gdd_obs","prec_obs","gddHigh_obs","dd_obs",
                "gdd_1to3","prec_1to3","gddHigh_1to3","dd_1to3",
                "gdd_4to10","prec_4to10","gddHigh_4to10","dd_4to10",
                "gdd_obs_land", "prec_obs_land", "gddHigh_obs_land", "dd_obs_land", 
                "gdd_1to3_land", "prec_1to3_land", "gddHigh_1to3_land", "dd_1to3_land",
                "gdd_4to10_land", "prec_4to10_land", "gddHigh_4to10_land", "dd_4to10_land",
                "gdd_obs_1to3","prec_obs_1to3","gddHigh_obs_1to3","dd_obs_1to3",
                "gdd_obs_4to10","prec_obs_4to10","gddHigh_obs_4to10","dd_obs_4to10",
                "np_cereals_fm","np_protein_fm","np_oilseed_fm","np_roots_fm","np_corn_fm","nw_fert_fm",
                "k_land_fm","k_labor_fm","k_capital_fm","trend_fm","trend2_fm",
                "gdd_obs_fm","prec_obs_fm","gddHigh_obs_fm","dd_obs_fm",
                "gdd_1to3_fm","prec_1to3_fm","gddHigh_1to3_fm","dd_1to3_fm",
                "gdd_4to10_fm","prec_4to10_fm","gddHigh_4to10_fm","dd_4to10_fm",
                "gdd_obs_land_fm", "prec_obs_land_fm", "gddHigh_obs_land_fm", "dd_obs_land_fm", 
                "gdd_1to3_land_fm", "prec_1to3_land_fm", "gddHigh_1to3_land_fm", "dd_1to3_land_fm", 
                "gdd_4to10_land_fm", "prec_4to10_land_fm", "gddHigh_4to10_land_fm", "dd_4to10_land_fm",
                "gdd_obs_1to3_fm","prec_obs_1to3_fm","gddHigh_obs_1to3_fm","dd_obs_1to3_fm",
                "gdd_obs_4to10_fm","prec_obs_4to10_fm","gddHigh_obs_4to10_fm","dd_obs_4to10_fm",
                "phi")
  
  # Simulate the outcome
  for (t in -1:10) {
    
    X_sim <- X
    
    # Replace the weather variables for the simulation
    X_sim$gdd_obs <- sim2018$gdd_obs[sim2018$t==t]
    X_sim$prec_obs <- sim2018$prec_obs[sim2018$t==t]
    X_sim$gddHigh_obs <- sim2018$gddHigh_obs[sim2018$t==t]
    X_sim$dd_obs <- sim2018$dd_obs[sim2018$t==t]
    
    X_sim$gdd_1to3 <- sim2018$gdd_1to3[sim2018$t==t]
    X_sim$prec_1to3 <- sim2018$prec_1to3[sim2018$t==t]
    X_sim$gddHigh_1to3 <- sim2018$gddHigh_1to3[sim2018$t==t]
    X_sim$dd_1to3 <- sim2018$dd_1to3[sim2018$t==t]
    
    X_sim$gdd_4to10 <- sim2018$gdd_4to10[sim2018$t==t]
    X_sim$prec_4to10 <- sim2018$prec_4to10[sim2018$t==t]
    X_sim$gddHigh_4to10 <- sim2018$gddHigh_4to10[sim2018$t==t]
    X_sim$dd_4to10  <- sim2018$dd_4to10[sim2018$t==t]
    
    X_sim$gdd_obs_land <- sim2018$gdd_obs[sim2018$t==t] * X$k_land
    X_sim$prec_obs_land <- sim2018$prec_obs[sim2018$t==t] * X$k_land
    X_sim$gddHigh_obs_land <- sim2018$gddHigh_obs[sim2018$t==t] * X$k_land
    X_sim$dd_obs_land <- sim2018$dd_obs[sim2018$t==t] * X$k_land
    X_sim$gdd_1to3_land <- sim2018$gdd_1to3[sim2018$t==t] * X$k_land
    X_sim$prec_1to3_land <- sim2018$prec_1to3[sim2018$t==t] * X$k_land
    X_sim$gddHigh_1to3_land <- sim2018$gddHigh_1to3[sim2018$t==t] * X$k_land
    X_sim$dd_1to3_land <- sim2018$dd_1to3[sim2018$t==t] * X$k_land
    X_sim$gdd_4to10_land <- sim2018$gdd_4to10[sim2018$t==t] * X$k_land
    X_sim$prec_4to10_land <- sim2018$prec_4to10[sim2018$t==t] * X$k_land
    X_sim$gddHigh_4to10_land <- sim2018$gddHigh_4to10[sim2018$t==t] * X$k_land
    X_sim$dd_4to10_land <- sim2018$dd_4to10[sim2018$t==t] * X$k_land
    
    X_sim$gdd_obs_1to3 <- sim2018$gdd_obs_1to3[sim2018$t==t]
    X_sim$prec_obs_1to3 <- sim2018$prec_obs_1to3[sim2018$t==t]
    X_sim$gddHigh_obs_1to3 <- sim2018$gddHigh_obs_1to3[sim2018$t==t]
    X_sim$dd_obs_1to3 <- sim2018$dd_obs_1to3[sim2018$t==t]
    X_sim$gdd_obs_4to10 <- sim2018$gdd_obs_4to10[sim2018$t==t]
    X_sim$prec_obs_4to10 <- sim2018$prec_obs_4to10[sim2018$t==t]
    X_sim$gddHigh_obs_4to10 <- sim2018$gddHigh_obs_4to10[sim2018$t==t]
    X_sim$dd_obs_4to10 <- sim2018$dd_obs_4to10[sim2018$t==t]
    
    # Replace phi and PHI for all crops, but not for fert
    if (list_outcome[i,"outcome"] != "fert") {
      X_sim$phi <- do.call(c,list_probit_sim[[i]][[2]][t+2])
      PHI <- do.call(c,list_probit_sim[[i]][[1]][t+2])
    } else {
      X_sim$phi <- NULL
      PHI <- 1
    } 
    
    # Predict output 
    X_sim <- as.matrix(X_sim)
    C <- as.data.frame(t(coef_linear) ) %>% select(starts_with(list_outcome[i,"sur"]))
    C <- t(C) 
    
    list_simresults[[i]][[t+2]] <- PHI * (X_sim %*% C) 
    
  }
}

# Summarize simulation results at sample mean, for large farms, and for small farms

  #--------------------#
  # At the sample mean #
  #--------------------#

  simresults <- data.frame(
    t = c(-1:10),
    Cereals = c(mean(list_simresults$Cereals$t1), mean(list_simresults$Cereals$t2), mean(list_simresults$Cereals$t3),
                mean(list_simresults$Cereals$t4), mean(list_simresults$Cereals$t5), mean(list_simresults$Cereals$t6),
                mean(list_simresults$Cereals$t7), mean(list_simresults$Cereals$t8), mean(list_simresults$Cereals$t9),
                mean(list_simresults$Cereals$t10), mean(list_simresults$Cereals$t11), mean(list_simresults$Cereals$t12)),
    Protein = c(mean(list_simresults$Protein$t1), mean(list_simresults$Protein$t2), mean(list_simresults$Protein$t3),
                mean(list_simresults$Protein$t4), mean(list_simresults$Protein$t5), mean(list_simresults$Protein$t6),
                mean(list_simresults$Protein$t7), mean(list_simresults$Protein$t8), mean(list_simresults$Protein$t9),
                mean(list_simresults$Protein$t10), mean(list_simresults$Protein$t11), mean(list_simresults$Protein$t12)),
    Oilseed = c(mean(list_simresults$Oilseed$t1), mean(list_simresults$Oilseed$t2), mean(list_simresults$Oilseed$t3),
                mean(list_simresults$Oilseed$t4), mean(list_simresults$Oilseed$t5), mean(list_simresults$Oilseed$t6),
                mean(list_simresults$Oilseed$t7), mean(list_simresults$Oilseed$t8), mean(list_simresults$Oilseed$t9),
                mean(list_simresults$Oilseed$t10), mean(list_simresults$Oilseed$t11), mean(list_simresults$Oilseed$t12)),
    Roots = c(mean(list_simresults$Roots$t1), mean(list_simresults$Roots$t2), mean(list_simresults$Roots$t3),
              mean(list_simresults$Roots$t4), mean(list_simresults$Roots$t5), mean(list_simresults$Roots$t6),
              mean(list_simresults$Roots$t7), mean(list_simresults$Roots$t8), mean(list_simresults$Roots$t9),
              mean(list_simresults$Roots$t10), mean(list_simresults$Roots$t11), mean(list_simresults$Roots$t12)),
    Corn = c(mean(list_simresults$Corn$t1), mean(list_simresults$Corn$t2), mean(list_simresults$Corn$t3),
                mean(list_simresults$Corn$t4), mean(list_simresults$Corn$t5), mean(list_simresults$Corn$t6),
                mean(list_simresults$Corn$t7), mean(list_simresults$Corn$t8), mean(list_simresults$Corn$t9),
                mean(list_simresults$Corn$t10), mean(list_simresults$Corn$t11), mean(list_simresults$Corn$t12)),
    Fertilizer = -c(mean(list_simresults$Fertilizer$t1), mean(list_simresults$Fertilizer$t2), mean(list_simresults$Fertilizer$t3),
                    mean(list_simresults$Fertilizer$t4), mean(list_simresults$Fertilizer$t5), mean(list_simresults$Fertilizer$t6),
                    mean(list_simresults$Fertilizer$t7), mean(list_simresults$Fertilizer$t8), mean(list_simresults$Fertilizer$t9),
                    mean(list_simresults$Fertilizer$t10), mean(list_simresults$Fertilizer$t11), mean(list_simresults$Fertilizer$t12))
  )

  # Save simulation result
  save(simresults, file="rOutput/simresults_samplemean.Rda")

  #-----------------#
  # For large farms #
  #-----------------#
  
  # Add results from list to data frame so results can be summarized by size group
  df_farm$cereals_tm1 <- as.numeric(do.call(c, list_simresults[[1]][1]))
  df_farm$cereals_t0 <-  as.numeric(do.call(c, list_simresults[[1]][2]))
  df_farm$cereals_t1 <-  as.numeric(do.call(c, list_simresults[[1]][3]))
  df_farm$cereals_t2 <-  as.numeric(do.call(c, list_simresults[[1]][4]))
  df_farm$cereals_t3 <-  as.numeric(do.call(c, list_simresults[[1]][5]))
  df_farm$cereals_t4 <-  as.numeric(do.call(c, list_simresults[[1]][6]))
  df_farm$cereals_t5 <-  as.numeric(do.call(c, list_simresults[[1]][7]))
  df_farm$cereals_t6 <-  as.numeric(do.call(c, list_simresults[[1]][8]))
  df_farm$cereals_t7 <-  as.numeric(do.call(c, list_simresults[[1]][9]))
  df_farm$cereals_t8 <-  as.numeric(do.call(c, list_simresults[[1]][10]))
  df_farm$cereals_t9 <-  as.numeric(do.call(c, list_simresults[[1]][11]))
  df_farm$cereals_t10 <- as.numeric(do.call(c, list_simresults[[1]][12]))
  
  df_farm$protein_tm1 <- as.numeric(do.call(c, list_simresults[[2]][1]))
  df_farm$protein_t0 <-  as.numeric(do.call(c, list_simresults[[2]][2]))
  df_farm$protein_t1 <-  as.numeric(do.call(c, list_simresults[[2]][3]))
  df_farm$protein_t2 <-  as.numeric(do.call(c, list_simresults[[2]][4]))
  df_farm$protein_t3 <-  as.numeric(do.call(c, list_simresults[[2]][5]))
  df_farm$protein_t4 <-  as.numeric(do.call(c, list_simresults[[2]][6]))
  df_farm$protein_t5 <-  as.numeric(do.call(c, list_simresults[[2]][7]))
  df_farm$protein_t6 <-  as.numeric(do.call(c, list_simresults[[2]][8]))
  df_farm$protein_t7 <-  as.numeric(do.call(c, list_simresults[[2]][9]))
  df_farm$protein_t8 <-  as.numeric(do.call(c, list_simresults[[2]][10]))
  df_farm$protein_t9 <-  as.numeric(do.call(c, list_simresults[[2]][11]))
  df_farm$protein_t10 <- as.numeric(do.call(c, list_simresults[[2]][12]))
  
  df_farm$oilseed_tm1 <- as.numeric(do.call(c, list_simresults[[3]][1]))
  df_farm$oilseed_t0 <-  as.numeric(do.call(c, list_simresults[[3]][2]))
  df_farm$oilseed_t1 <-  as.numeric(do.call(c, list_simresults[[3]][3]))
  df_farm$oilseed_t2 <-  as.numeric(do.call(c, list_simresults[[3]][4]))
  df_farm$oilseed_t3 <-  as.numeric(do.call(c, list_simresults[[3]][5]))
  df_farm$oilseed_t4 <-  as.numeric(do.call(c, list_simresults[[3]][6]))
  df_farm$oilseed_t5 <-  as.numeric(do.call(c, list_simresults[[3]][7]))
  df_farm$oilseed_t6 <-  as.numeric(do.call(c, list_simresults[[3]][8]))
  df_farm$oilseed_t7 <-  as.numeric(do.call(c, list_simresults[[3]][9]))
  df_farm$oilseed_t8 <-  as.numeric(do.call(c, list_simresults[[3]][10]))
  df_farm$oilseed_t9 <-  as.numeric(do.call(c, list_simresults[[3]][11]))
  df_farm$oilseed_t10 <- as.numeric(do.call(c, list_simresults[[3]][12]))
  
  df_farm$roots_tm1 <- as.numeric(do.call(c, list_simresults[[4]][1]))
  df_farm$roots_t0 <-  as.numeric(do.call(c, list_simresults[[4]][2]))
  df_farm$roots_t1 <-  as.numeric(do.call(c, list_simresults[[4]][3]))
  df_farm$roots_t2 <-  as.numeric(do.call(c, list_simresults[[4]][4]))
  df_farm$roots_t3 <-  as.numeric(do.call(c, list_simresults[[4]][5]))
  df_farm$roots_t4 <-  as.numeric(do.call(c, list_simresults[[4]][6]))
  df_farm$roots_t5 <-  as.numeric(do.call(c, list_simresults[[4]][7]))
  df_farm$roots_t6 <-  as.numeric(do.call(c, list_simresults[[4]][8]))
  df_farm$roots_t7 <-  as.numeric(do.call(c, list_simresults[[4]][9]))
  df_farm$roots_t8 <-  as.numeric(do.call(c, list_simresults[[4]][10]))
  df_farm$roots_t9 <-  as.numeric(do.call(c, list_simresults[[4]][11]))
  df_farm$roots_t10 <- as.numeric(do.call(c, list_simresults[[4]][12]))
  
  df_farm$corn_tm1 <- as.numeric(do.call(c, list_simresults[[5]][1]))
  df_farm$corn_t0 <-  as.numeric(do.call(c, list_simresults[[5]][2]))
  df_farm$corn_t1 <-  as.numeric(do.call(c, list_simresults[[5]][3]))
  df_farm$corn_t2 <-  as.numeric(do.call(c, list_simresults[[5]][4]))
  df_farm$corn_t3 <-  as.numeric(do.call(c, list_simresults[[5]][5]))
  df_farm$corn_t4 <-  as.numeric(do.call(c, list_simresults[[5]][6]))
  df_farm$corn_t5 <-  as.numeric(do.call(c, list_simresults[[5]][7]))
  df_farm$corn_t6 <-  as.numeric(do.call(c, list_simresults[[5]][8]))
  df_farm$corn_t7 <-  as.numeric(do.call(c, list_simresults[[5]][9]))
  df_farm$corn_t8 <-  as.numeric(do.call(c, list_simresults[[5]][10]))
  df_farm$corn_t9 <-  as.numeric(do.call(c, list_simresults[[5]][11]))
  df_farm$corn_t10 <- as.numeric(do.call(c, list_simresults[[5]][12]))
  
  df_farm$fert_tm1 <- -as.numeric(do.call(c, list_simresults[[6]][1]))
  df_farm$fert_t0 <-  -as.numeric(do.call(c, list_simresults[[6]][2]))
  df_farm$fert_t1 <-  -as.numeric(do.call(c, list_simresults[[6]][3]))
  df_farm$fert_t2 <-  -as.numeric(do.call(c, list_simresults[[6]][4]))
  df_farm$fert_t3 <-  -as.numeric(do.call(c, list_simresults[[6]][5]))
  df_farm$fert_t4 <-  -as.numeric(do.call(c, list_simresults[[6]][6]))
  df_farm$fert_t5 <-  -as.numeric(do.call(c, list_simresults[[6]][7]))
  df_farm$fert_t6 <-  -as.numeric(do.call(c, list_simresults[[6]][8]))
  df_farm$fert_t7 <-  -as.numeric(do.call(c, list_simresults[[6]][9]))
  df_farm$fert_t8 <-  -as.numeric(do.call(c, list_simresults[[6]][10]))
  df_farm$fert_t9 <-  -as.numeric(do.call(c, list_simresults[[6]][11]))
  df_farm$fert_t10 <- -as.numeric(do.call(c, list_simresults[[6]][12]))
  
  
  simresults <- data.frame(
    t = c(-1:10),
    Cereals = c(mean(df_farm$cereals_tm1[df_farm$k_land > median(df_farm$k_land)]), 
                mean(df_farm$cereals_t0[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$cereals_t1[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$cereals_t2[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$cereals_t3[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$cereals_t4[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$cereals_t5[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$cereals_t6[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$cereals_t7[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$cereals_t8[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$cereals_t9[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$cereals_t10[df_farm$k_land > median(df_farm$k_land)])),
    Protein = c(mean(df_farm$protein_tm1[df_farm$k_land > median(df_farm$k_land)]), 
                mean(df_farm$protein_t0[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$protein_t1[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$protein_t2[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$protein_t3[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$protein_t4[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$protein_t5[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$protein_t6[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$protein_t7[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$protein_t8[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$protein_t9[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$protein_t10[df_farm$k_land > median(df_farm$k_land)])),
    Oilseed = c(mean(df_farm$oilseed_tm1[df_farm$k_land > median(df_farm$k_land)]), 
                mean(df_farm$oilseed_t0[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$oilseed_t1[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$oilseed_t2[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$oilseed_t3[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$oilseed_t4[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$oilseed_t5[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$oilseed_t6[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$oilseed_t7[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$oilseed_t8[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$oilseed_t9[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$oilseed_t10[df_farm$k_land > median(df_farm$k_land)])),
    Roots = c(mean(df_farm$roots_tm1[df_farm$k_land > median(df_farm$k_land)]), 
              mean(df_farm$roots_t0[df_farm$k_land > median(df_farm$k_land)]),
              mean(df_farm$roots_t1[df_farm$k_land > median(df_farm$k_land)]),
              mean(df_farm$roots_t2[df_farm$k_land > median(df_farm$k_land)]),
              mean(df_farm$roots_t3[df_farm$k_land > median(df_farm$k_land)]),
              mean(df_farm$roots_t4[df_farm$k_land > median(df_farm$k_land)]),
              mean(df_farm$roots_t5[df_farm$k_land > median(df_farm$k_land)]),
              mean(df_farm$roots_t6[df_farm$k_land > median(df_farm$k_land)]),
              mean(df_farm$roots_t7[df_farm$k_land > median(df_farm$k_land)]),
              mean(df_farm$roots_t8[df_farm$k_land > median(df_farm$k_land)]),
              mean(df_farm$roots_t9[df_farm$k_land > median(df_farm$k_land)]),
              mean(df_farm$roots_t10[df_farm$k_land > median(df_farm$k_land)])),
    Corn = c(mean(df_farm$corn_tm1[df_farm$k_land > median(df_farm$k_land)]), 
                mean(df_farm$corn_t0[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$corn_t1[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$corn_t2[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$corn_t3[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$corn_t4[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$corn_t5[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$corn_t6[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$corn_t7[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$corn_t8[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$corn_t9[df_farm$k_land > median(df_farm$k_land)]),
                mean(df_farm$corn_t10[df_farm$k_land > median(df_farm$k_land)])),
    Fertilizer = c(mean(df_farm$fert_tm1[df_farm$k_land > median(df_farm$k_land)]), 
                   mean(df_farm$fert_t0[df_farm$k_land > median(df_farm$k_land)]),
                   mean(df_farm$fert_t1[df_farm$k_land > median(df_farm$k_land)]),
                   mean(df_farm$fert_t2[df_farm$k_land > median(df_farm$k_land)]),
                   mean(df_farm$fert_t3[df_farm$k_land > median(df_farm$k_land)]),
                   mean(df_farm$fert_t4[df_farm$k_land > median(df_farm$k_land)]),
                   mean(df_farm$fert_t5[df_farm$k_land > median(df_farm$k_land)]),
                   mean(df_farm$fert_t6[df_farm$k_land > median(df_farm$k_land)]),
                   mean(df_farm$fert_t7[df_farm$k_land > median(df_farm$k_land)]),
                   mean(df_farm$fert_t8[df_farm$k_land > median(df_farm$k_land)]),
                   mean(df_farm$fert_t9[df_farm$k_land > median(df_farm$k_land)]),
                   mean(df_farm$fert_t10[df_farm$k_land > median(df_farm$k_land)]))
  )
  
  
  # Save simulation result
  save(simresults, file="rOutput/simresults_large.Rda")
  
  #------------------#
  # For small farms  #
  #------------------#
  
  simresults <- data.frame(
    t = c(-1:10),
    Cereals = c(mean(df_farm$cereals_tm1[df_farm$k_land < median(df_farm$k_land)]), 
                mean(df_farm$cereals_t0[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$cereals_t1[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$cereals_t2[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$cereals_t3[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$cereals_t4[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$cereals_t5[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$cereals_t6[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$cereals_t7[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$cereals_t8[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$cereals_t9[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$cereals_t10[df_farm$k_land < median(df_farm$k_land)])),
    Protein = c(mean(df_farm$protein_tm1[df_farm$k_land < median(df_farm$k_land)]), 
                mean(df_farm$protein_t0[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$protein_t1[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$protein_t2[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$protein_t3[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$protein_t4[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$protein_t5[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$protein_t6[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$protein_t7[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$protein_t8[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$protein_t9[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$protein_t10[df_farm$k_land < median(df_farm$k_land)])),
    Oilseed = c(mean(df_farm$oilseed_tm1[df_farm$k_land < median(df_farm$k_land)]), 
                mean(df_farm$oilseed_t0[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$oilseed_t1[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$oilseed_t2[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$oilseed_t3[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$oilseed_t4[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$oilseed_t5[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$oilseed_t6[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$oilseed_t7[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$oilseed_t8[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$oilseed_t9[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$oilseed_t10[df_farm$k_land < median(df_farm$k_land)])),
    Roots = c(mean(df_farm$roots_tm1[df_farm$k_land < median(df_farm$k_land)]), 
              mean(df_farm$roots_t0[df_farm$k_land < median(df_farm$k_land)]),
              mean(df_farm$roots_t1[df_farm$k_land < median(df_farm$k_land)]),
              mean(df_farm$roots_t2[df_farm$k_land < median(df_farm$k_land)]),
              mean(df_farm$roots_t3[df_farm$k_land < median(df_farm$k_land)]),
              mean(df_farm$roots_t4[df_farm$k_land < median(df_farm$k_land)]),
              mean(df_farm$roots_t5[df_farm$k_land < median(df_farm$k_land)]),
              mean(df_farm$roots_t6[df_farm$k_land < median(df_farm$k_land)]),
              mean(df_farm$roots_t7[df_farm$k_land < median(df_farm$k_land)]),
              mean(df_farm$roots_t8[df_farm$k_land < median(df_farm$k_land)]),
              mean(df_farm$roots_t9[df_farm$k_land < median(df_farm$k_land)]),
              mean(df_farm$roots_t10[df_farm$k_land < median(df_farm$k_land)])),
    Corn = c(mean(df_farm$corn_tm1[df_farm$k_land < median(df_farm$k_land)]), 
                mean(df_farm$corn_t0[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$corn_t1[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$corn_t2[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$corn_t3[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$corn_t4[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$corn_t5[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$corn_t6[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$corn_t7[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$corn_t8[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$corn_t9[df_farm$k_land < median(df_farm$k_land)]),
                mean(df_farm$corn_t10[df_farm$k_land < median(df_farm$k_land)])),
    Fertilizer = c(mean(df_farm$fert_tm1[df_farm$k_land < median(df_farm$k_land)]), 
                   mean(df_farm$fert_t0[df_farm$k_land < median(df_farm$k_land)]),
                   mean(df_farm$fert_t1[df_farm$k_land < median(df_farm$k_land)]),
                   mean(df_farm$fert_t2[df_farm$k_land < median(df_farm$k_land)]),
                   mean(df_farm$fert_t3[df_farm$k_land < median(df_farm$k_land)]),
                   mean(df_farm$fert_t4[df_farm$k_land < median(df_farm$k_land)]),
                   mean(df_farm$fert_t5[df_farm$k_land < median(df_farm$k_land)]),
                   mean(df_farm$fert_t6[df_farm$k_land < median(df_farm$k_land)]),
                   mean(df_farm$fert_t7[df_farm$k_land < median(df_farm$k_land)]),
                   mean(df_farm$fert_t8[df_farm$k_land < median(df_farm$k_land)]),
                   mean(df_farm$fert_t9[df_farm$k_land < median(df_farm$k_land)]),
                   mean(df_farm$fert_t10[df_farm$k_land < median(df_farm$k_land)]))
  )
  
  # Save simulation result
  save(simresults, file="rOutput/simresults_small.Rda")
  
  
#---------------------------#
#### 5) Plot the results ####
#---------------------------#

#------------------------------------#
##### 5a) ... at the sample mean #####
#------------------------------------#

# Load simulation results
load("rOutput/simresults_samplemean.Rda")

# Normalize the data for plotting
simresults <- simresults %>% 
  mutate_each(funs(./.[1]), setdiff(names(.), c("t")))
  
# Transform data to long data
simresults_long <- reshape2::melt(simresults, id.vars="t")
  
# Plot the results
Fig_het_samplemean <- ggplot(simresults_long, aes(x=t, y=value, group=variable, shape=variable, color=variable)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_manual(values=c("#8c510a", "#bf812d", "#dfc27d", "#80cdc1", "#35978f", "#003c30"),
                     name = "",
                     labels = c("Cereals",
                                "Protein crops",
                                "Oilseeds",
                                "Root crops",
                                "Corn",
                                "Fertilizer")) +
  scale_shape_manual(values=c(10,16,17,18,15,13),
                     name = "",
                     labels = c("Cereals",
                                "Protein crops",
                                "Oilseeds",
                                "Root crops",
                                "Corn",
                                "Fertilizer")) +
  ggtitle("") +
  scale_x_continuous(breaks = seq(-1, 10, by = 1),
                     labels = c("-1","0","1","2","3","4",
                                "5","6","7","8","9","10"),
                     name = "Periods (t) since drought shock") +
  scale_y_continuous(breaks = seq(0.4, 1.4, by = 0.2), limits=c(0.4,1.52)) + 
  theme(text = element_text(size = 11),
        legend.text=element_text(size=9),
        axis.title = element_text(size = 10),
        legend.margin=margin(t = -0.7, unit='cm'),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank()) +
  labs(
    x = "",
    y = "Index",
    linetype = ""
  ) +
  geom_segment(aes(x = 0, y = 0.4, xend = 0, yend = 1.45), linetype="dashed", color="grey") +
  annotate(geom="text", x=0, y=1.52, label="Drought shock",
           color="darkgrey", size=3.5)

Fig_het_samplemean


#---------------------------------#
##### 5b) ... for large farms #####
#---------------------------------#

# Load simulation results
load("rOutput/simresults_large.Rda")

# Normalize the data for plotting
simresults <- simresults %>% 
  mutate_each(funs(./.[1]), setdiff(names(.), c("t")))

# Transform data to long data
simresults_long <- reshape2::melt(simresults, id.vars="t")

# Plot the results
Fig_het_large <- ggplot(simresults_long, aes(x=t, y=value, group=variable, shape=variable, color=variable)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_manual(values=c("#8c510a", "#bf812d", "#dfc27d", "#80cdc1", "#35978f", "#003c30"),
                     name = "",
                     labels = c("Cereals",
                                "Protein crops",
                                "Oilseeds",
                                "Root crops",
                                "Corn",
                                "Fertilizer")) +
  scale_shape_manual(values=c(10,16,17,18,15,13),
                     name = "",
                     labels = c("Cereals",
                                "Protein crops",
                                "Oilseeds",
                                "Root crops",
                                "Corn",
                                "Fertilizer")) +
  ggtitle("") +
  scale_x_continuous(breaks = seq(-1, 10, by = 1),
                     labels = c("-1","0","1","2","3","4",
                                "5","6","7","8","9","10"),
                     name = "Periods (t) since drought shock") +
  scale_y_continuous(breaks = seq(0.4, 1.4, by = 0.2), limits=c(0.4,1.52)) + 
  theme(text = element_text(size = 11),
        legend.text=element_text(size=9),
        axis.title = element_text(size = 10),
        legend.margin=margin(t = -0.7, unit='cm'),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank()) +
  labs(
    x = "",
    y = "Index",
    linetype = ""
  ) +
  geom_segment(aes(x = 0, y = 0.4, xend = 0, yend = 1.45), linetype="dashed", color="grey") +
  annotate(geom="text", x=0, y=1.52, label="Drought shock",
           color="darkgrey", size=3.5)

Fig_het_large

#---------------------------------#
##### 5c) ... for small farms #####
#---------------------------------#

# Load simulation results
load("rOutput/simresults_small.Rda")

# Normalize the data for plotting
simresults <- simresults %>% 
  mutate_each(funs(./.[1]), setdiff(names(.), c("t")))

# Transform data to long data
simresults_long <- reshape2::melt(simresults, id.vars="t")


# Plot the results
Fig_het_small <- ggplot(simresults_long, aes(x=t, y=value, group=variable, shape=variable, color=variable)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_manual(values=c("#8c510a", "#bf812d", "#dfc27d", "#80cdc1", "#35978f", "#003c30"),
                     name = "",
                     labels = c("Cereals",
                                "Protein crops",
                                "Oilseeds",
                                "Root crops",
                                "Corn",
                                "Fertilizer")) +
  scale_shape_manual(values=c(10,16,17,18,15,13),
                     name = "",
                     labels = c("Cereals",
                                "Protein crops",
                                "Oilseeds",
                                "Root crops",
                                "Corn",
                                "Fertilizer")) +
  ggtitle("") +
  scale_x_continuous(breaks = seq(-1, 10, by = 1),
                     labels = c("-1","0","1","2","3","4",
                                "5","6","7","8","9","10"),
                     name = "Periods (t) since drought shock") +
  scale_y_continuous(breaks = seq(0.4, 1.4, by = 0.2), limits=c(0.4,1.52)) + 
  theme(text = element_text(size = 11),
        legend.text=element_text(size=9),
        axis.title = element_text(size = 10),
        legend.margin=margin(t = -0.7, unit='cm'),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank()) +
  labs(
    x = "",
    y = "Index",
    linetype = ""
  ) +
  geom_segment(aes(x = 0, y = 0.4, xend = 0, yend = 1.45), linetype="dashed", color="grey") +
  annotate(geom="text", x=0, y=1.52, label="Drought shock",
           color="darkgrey", size=3.5)

Fig_het_small


# ------------------------ #
# Make three plots in one  #
# -------------------------#

Fig_het_all <- ggarrange(Fig_het_small,
                         Fig_het_samplemean,
                         Fig_het_large,
                         ncol = 1, nrow=3,
                         labels=c("Small farms", "Sample mean", "Large farms"),
                         font.label = list(size = 11, color = "black", face = "italic"),
                         legend = "right",
                         common.legend = TRUE)
Fig_het_all

# Save the plot as .png and .eps files
ggsave("Figures/Figure_4.png", Fig_het_all, device="png", width = 6, height = 9, units = "in", dpi=1200)
ggsave("Figures/Figure_4.eps", Fig_het_all, device="eps", width = 6, height = 9, units = "in", dpi=1200)


