#-----------------------------------------------#
#                                               #
# This program estimates the main model without #
# the year 2019 and provides an out-of-sample   # 
# prediction for the year 2019.                 #
#                                               #
# 1) Some preliminaries                         #
# 2) Selection equations                        #
# 3) Structural equations                       #
# 4) Out-of-sample prediction                   #
#                                               #
#-----------------------------------------------#

library(dplyr)
library(GJRM)
library(systemfit)

#-----------------------------#
#### 0) Some preliminaries ####
#-----------------------------#

load("rOutput/farm_ready.Rda")

# Generate subsample without 2019
df_farm_no2019 <- df_farm %>%
  filter(year < 2019)

# Generate subsample with 2019 only (in this case they already have the "right" Mundlaks)
df_farm_2019 <- df_farm %>%
  filter(year == 2019)

# Generate Mundlak devices again because the sample is different now

allindepvars <- c("np_cereals", "np_oilseed", "np_roots", 
                  "np_protein", "np_corn", "nw_fert", "k_land", "k_labor",
                  "k_capital", "trend", "trend2", 
                  "aa_crops", "p_cereals", "p_oilseed", "p_roots", "p_protein", "p_corn",
                  "sh_cereals", "sh_protein", "sh_oilseed", "sh_roots", "sh_corn",
                  "lsh_cereals", "lsh_protein", "lsh_oilseed", "lsh_roots", "lsh_corn",
                  "gdd_obs", "gdd2_obs", "prec_obs", "prec2_obs", "gddHigh_obs", "dd_obs",
                  "gdd_1to3", "gdd2_1to3", "prec_1to3", "prec2_1to3", "gddHigh_1to3", "dd_1to3",
                  "gdd_4to10","gdd2_4to10", "prec_4to10", "prec2_4to10", "gddHigh_4to10", "dd_4to10",
                  "gdd_1to5", "gdd2_1to5", "prec_1to5", "prec2_1to5", "gddHigh_1to5", "dd_1to5",
                  "gdd_6to10", "gdd2_6to10",  "prec_6to10", "prec2_6to10", "gddHigh_6to10", "dd_6to10",
                  "gdd_6to20", "gdd2_6to20",  "prec_6to20", "prec2_6to20", "gddHigh_6to20", "dd_6to20",
                  "gdd_1", "gdd2_1", "prec_1", "prec2_1", "gddHigh_1", "dd_1",
                  "gdd_2to5", "gdd2_2to5",  "prec_2to5", "prec2_2to5", "gddHigh_2to5", "dd_2to5",
                  "gdd_2to10", "gdd2_2to10", "prec_2to10", "prec2_2to10", "gddHigh_2to10", "dd_2to10",
                  "gdd_1to10", "gdd2_1to10", "prec_1to10", "prec2_1to10", "gddHigh_1to10", "dd_1to10",
                  "gdd_1to30", "gdd2_1to30", "prec_1to30", "prec2_1to30", "gddHigh_1to30", "dd_1to30",
                  "gdd_obs_land" ,"prec_obs_land", "gddHigh_obs_land", "dd_obs_land",
                  "gdd_1to3_land", "prec_1to3_land", "gddHigh_1to3_land", "dd_1to3_land",
                  "gdd_4to10_land", "prec_4to10_land", "gddHigh_4to10_land", "dd_4to10_land",
                  "gdd_1to5_land", "prec_1to5_land", "gddHigh_1to5_land", "dd_1to5_land", 
                  "gdd_6to10_land", "prec_6to10_land", "gddHigh_6to10_land", "dd_6to10_land",
                  "gdd_6to20_land", "prec_6to20_land", "gddHigh_6to20_land", "dd_6to20_land",
                  "gdd_1_land", "prec_1_land", "gddHigh_1_land", "dd_1_land",
                  "gdd_2to5_land", "prec_2to5_land", "gddHigh_2to5_land", "dd_2to5_land",
                  "gdd_2to10_land", "prec_2to10_land", "gddHigh_2to10_land", "dd_2to10_land",
                  "gdd_1to10_land", "prec_1to10_land", "gddHigh_1to10_land", "dd_1to10_land",
                  "gdd_1to30_land", "prec_1to30_land", "gddHigh_1to30_land", "dd_1to30_land",
                  "gdd_obs_1to3", "gdd2_obs_1to3", "prec_obs_1to3", "prec2_obs_1to3", "gddHigh_obs_1to3", "dd_obs_1to3",
                  "gdd_obs_4to10", "gdd2_obs_4to10", "prec_obs_4to10", "prec2_obs_4to10", "gddHigh_obs_4to10", "dd_obs_4to10",
                  "gdd_obs_1to5", "gdd2_obs_1to5", "prec_obs_1to5", "prec2_obs_1to5", "gddHigh_obs_1to5", "dd_obs_1to5", 
                  "gdd_obs_6to10", "gdd2_obs_6to10", "prec_obs_6to10", "prec2_obs_6to10", "gddHigh_obs_6to10", "dd_obs_6to10",
                  "gdd_obs_6to20", "gdd2_obs_6to20", "prec_obs_6to20", "prec2_obs_6to20", "gddHigh_obs_6to20", "dd_obs_6to20",
                  "gdd_obs_1", "gdd2_obs_1", "prec_obs_1", "prec2_obs_1", "gddHigh_obs_1", "dd_obs_1", 
                  "gdd_obs_2to10", "gdd2_obs_2to10", "prec_obs_2to10", "prec2_obs_2to10", "gddHigh_obs_2to10", "dd_obs_2to10")

df_farm_no2019_mundl <- df_farm_no2019 %>%
  group_by(key) %>%
  summarise_at(c(all_of(allindepvars)), list(fm = mean)) 
df_farm_no2019 <- left_join(df_farm_no2019, df_farm_no2019_mundl)

#-------------------------------#
#### 2) Selection equations #####
#-------------------------------#

#Define probit regressions
cereal.eq <- iCereals ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
  k_land + k_labor + k_capital + trend + trend2 + 
  lsh_cereals + lsh_protein + lsh_oilseed + lsh_roots + 
  gdd_1to3 + prec_1to3 + gddHigh_1to3 + dd_1to3 +
  gdd_4to10 + prec_4to10 + gddHigh_4to10 + dd_4to10 + 
  np_cereals_fm + np_protein_fm + np_oilseed_fm + np_roots_fm + np_corn_fm + nw_fert_fm + 
  k_land_fm + k_labor_fm + k_capital_fm + trend_fm + trend2_fm + 
  lsh_cereals_fm + lsh_protein_fm + lsh_oilseed_fm + lsh_roots_fm + 
  gdd_1to3_fm + prec_1to3_fm + gddHigh_1to3_fm + dd_1to3_fm +
  gdd_4to10_fm + prec_4to10_fm + gddHigh_4to10_fm + dd_4to10_fm

protein.eq <- iProtein ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
  k_land + k_labor + k_capital + trend + trend2 + 
  lsh_cereals + lsh_protein + lsh_oilseed + lsh_roots + 
  gdd_1to3 + prec_1to3 + gddHigh_1to3 + dd_1to3 +
  gdd_4to10 + prec_4to10 + gddHigh_4to10 + dd_4to10 + 
  np_cereals_fm + np_protein_fm + np_oilseed_fm + np_roots_fm + np_corn_fm + nw_fert_fm + 
  k_land_fm + k_labor_fm + k_capital_fm + trend_fm + trend2_fm + 
  lsh_cereals_fm + lsh_protein_fm + lsh_oilseed_fm + lsh_roots_fm + 
  gdd_1to3_fm + prec_1to3_fm + gddHigh_1to3_fm + dd_1to3_fm +
  gdd_4to10_fm + prec_4to10_fm + gddHigh_4to10_fm + dd_4to10_fm

oilseed.eq <- iOilseed ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
  k_land + k_labor + k_capital + trend + trend2 + 
  lsh_cereals + lsh_protein + lsh_oilseed + lsh_roots + 
  gdd_1to3 + prec_1to3 + gddHigh_1to3 + dd_1to3 +
  gdd_4to10 + prec_4to10 + gddHigh_4to10 + dd_4to10 + 
  np_cereals_fm + np_protein_fm + np_oilseed_fm + np_roots_fm + np_corn_fm + nw_fert_fm + 
  k_land_fm + k_labor_fm + k_capital_fm + trend_fm + trend2_fm + 
  lsh_cereals_fm + lsh_protein_fm + lsh_oilseed_fm + lsh_roots_fm + 
  gdd_1to3_fm + prec_1to3_fm + gddHigh_1to3_fm + dd_1to3_fm +
  gdd_4to10_fm + prec_4to10_fm + gddHigh_4to10_fm + dd_4to10_fm

roots.eq <- iRoots ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
  k_land + k_labor + k_capital + trend + trend2 + 
  lsh_cereals + lsh_protein + lsh_oilseed + lsh_roots + 
  gdd_1to3 + prec_1to3 + gddHigh_1to3 + dd_1to3 +
  gdd_4to10 + prec_4to10 + gddHigh_4to10 + dd_4to10 + 
  np_cereals_fm + np_protein_fm + np_oilseed_fm + np_roots_fm + np_corn_fm + nw_fert_fm + 
  k_land_fm + k_labor_fm + k_capital_fm + trend_fm + trend2_fm + 
  lsh_cereals_fm + lsh_protein_fm + lsh_oilseed_fm + lsh_roots_fm + 
  gdd_1to3_fm + prec_1to3_fm + gddHigh_1to3_fm + dd_1to3_fm +
  gdd_4to10_fm + prec_4to10_fm + gddHigh_4to10_fm + dd_4to10_fm

corn.eq <- iCorn ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
  k_land + k_labor + k_capital + trend + trend2 + 
  lsh_cereals + lsh_protein + lsh_oilseed + lsh_roots + 
  gdd_1to3 + prec_1to3 + gddHigh_1to3 + dd_1to3 +
  gdd_4to10 + prec_4to10 + gddHigh_4to10 + dd_4to10 + 
  np_cereals_fm + np_protein_fm + np_oilseed_fm + np_roots_fm + np_corn_fm + nw_fert_fm + 
  k_land_fm + k_labor_fm + k_capital_fm + trend_fm + trend2_fm + 
  lsh_cereals_fm + lsh_protein_fm + lsh_oilseed_fm + lsh_roots_fm + 
  gdd_1to3_fm + prec_1to3_fm + gddHigh_1to3_fm + dd_1to3_fm +
  gdd_4to10_fm + prec_4to10_fm + gddHigh_4to10_fm + dd_4to10_fm

# Before I run the probit regressions, I create a data set containing the RHS-variables 
# of the excluded subsample (2019) so I can easily make the predictions.
iCereals <- glm(cereal.eq, family = binomial(link = "probit"), 
                data = df_farm)
vars <- colnames(model.matrix(iCereals))[-1]
model.matrix.2019 <- cbind(1,df_farm_2019[vars])

#-------------------------#
# Selection for q_cereals #
#-------------------------#

# Probit estimation
iCereals <- glm(cereal.eq, family = binomial(link = "probit"), 
                data = df_farm_no2019)

  # Predictions for in-sample (1996-2018)

    # Linear prediction
    df_farm_no2019$iCereals_linpred <- as.numeric(crossprod(matrix(iCereals$coefficients),
                                     t(model.matrix(iCereals))))
    # PHI
    df_farm_no2019$PHI_cereals <- pnorm(df_farm_no2019$iCereals_linpred)
    
    # phi
    df_farm_no2019$phi_cereals <- dnorm(df_farm_no2019$iCereals_linpred)
    
  # Predictions for out-of-sample (2019)
      
    # Linear prediction
    df_farm_2019$iCereals_linpred <- as.numeric(crossprod(matrix(iCereals$coefficients),
                                                    t(model.matrix.2019)))
    # PHI
    df_farm_2019$PHI_cereals <- pnorm(df_farm_2019$iCereals_linpred)
    
    # phi
    df_farm_2019$phi_cereals <- dnorm(df_farm_2019$iCereals_linpred) 
    
#-------------------------#
# Selection for q_protein #
#-------------------------#

# Probit estimation
iProtein <- glm(protein.eq, family = binomial(link = "probit"), 
                data = df_farm_no2019)

# Predictions for in-sample (1996-2018)

    # Linear prediction
    df_farm_no2019$iProtein_linpred <- as.numeric(crossprod(matrix(iProtein$coefficients),
                                                 t(model.matrix(iProtein))))
    # PHI
    df_farm_no2019$PHI_protein <- pnorm(df_farm_no2019$iProtein_linpred)
    
    # phi
    df_farm_no2019$phi_protein <- dnorm(df_farm_no2019$iProtein_linpred)

# Predictions for out-of-sample (2019)

    # Linear prediction
    df_farm_2019$iProtein_linpred <- as.numeric(crossprod(matrix(iProtein$coefficients),
                                                          t(model.matrix.2019)))
    # PHI
    df_farm_2019$PHI_protein <- pnorm(df_farm_2019$iProtein_linpred)
    
    # phi
    df_farm_2019$phi_protein <- dnorm(df_farm_2019$iProtein_linpred) 

#-------------------------#
# Selection for q_oilseed #
#-------------------------#

# Probit estimation
iOilseed <- glm(oilseed.eq, family = binomial(link = "probit"), 
                data = df_farm_no2019)

# Predictions for in-sample (1996-2018)

    # Linear prediction
    df_farm_no2019$iOilseed_linpred <- as.numeric(crossprod(matrix(iOilseed$coefficients),
                                                 t(model.matrix(iOilseed))))
    # PHI
    df_farm_no2019$PHI_oilseed <- pnorm(df_farm_no2019$iOilseed_linpred)
    
    # phi
    df_farm_no2019$phi_oilseed <- dnorm(df_farm_no2019$iOilseed_linpred)

# Predictions for out-of-sample (2019)

    # Linear prediction
    df_farm_2019$iOilseed_linpred <- as.numeric(crossprod(matrix(iOilseed$coefficients),
                                                          t(model.matrix.2019)))
    # PHI
    df_farm_2019$PHI_oilseed <- pnorm(df_farm_2019$iOilseed_linpred)
    
    # phi
    df_farm_2019$phi_oilseed <- dnorm(df_farm_2019$iOilseed_linpred) 

#-----------------------#
# Selection for q_roots #
#-----------------------#

# Probit estimation
iRoots <- glm(roots.eq, family = binomial(link = "probit"), 
              data = df_farm_no2019)

# Predictions for in-sample (1996-2018)

    # Linear prediction
    df_farm_no2019$iRoots_linpred <- as.numeric(crossprod(matrix(iRoots$coefficients),
                                               t(model.matrix(iRoots))))
    # PHI
    df_farm_no2019$PHI_roots <- pnorm(df_farm_no2019$iRoots_linpred)
    
    # phi
    df_farm_no2019$phi_roots <- dnorm(df_farm_no2019$iRoots_linpred)

# Predictions for out-of-sample (2019)

    # Linear prediction
    df_farm_2019$iRoots_linpred <- as.numeric(crossprod(matrix(iRoots$coefficients),
                                                        t(model.matrix.2019)))
    # PHI
    df_farm_2019$PHI_roots <- pnorm(df_farm_2019$iRoots_linpred)
    
    # phi
    df_farm_2019$phi_roots <- dnorm(df_farm_2019$iRoots_linpred) 

#----------------------#
# Selection for q_corn #
#----------------------#

# Probit estimation
iCorn <- glm(corn.eq, family = binomial(link = "probit"), 
                data = df_farm_no2019)

# Predictions for in-sample (1996-2018)

    # Linear prediction
    df_farm_no2019$iCorn_linpred <- as.numeric(crossprod(matrix(iCorn$coefficients),
                                                 t(model.matrix(iCorn))))
    # PHI
    df_farm_no2019$PHI_corn <- pnorm(df_farm_no2019$iCorn_linpred)
    
    # phi
    df_farm_no2019$phi_corn <- dnorm(df_farm_no2019$iCorn_linpred)

# Predictions for out-of-sample (2019)

    # Linear prediction
    df_farm_2019$iCorn_linpred <- as.numeric(crossprod(matrix(iCorn$coefficients),
                                                          t(model.matrix.2019)))
    # PHI
    df_farm_2019$PHI_corn <- pnorm(df_farm_2019$iCorn_linpred)
    
    # phi
    df_farm_2019$phi_corn <- dnorm(df_farm_2019$iCorn_linpred) 

#-------------------------------#
#### 3) Structural equations ####
#-------------------------------#
    
    #Define structural equations
    eqQQcereals <- qq_cereals ~ PHI_cereals + I(PHI_cereals*np_cereals) + I(PHI_cereals*np_protein) + I(PHI_cereals*np_oilseed) + I(PHI_cereals*np_roots) + I(PHI_cereals*np_corn) + I(PHI_cereals*nw_fert) + 
      I(PHI_cereals*k_land) + I(PHI_cereals*k_labor) + I(PHI_cereals*k_capital) + I(PHI_cereals*trend) + I(PHI_cereals*trend2) +
      I(PHI_cereals*gdd_obs) + I(PHI_cereals*prec_obs) + I(PHI_cereals*gddHigh_obs) + I(PHI_cereals*dd_obs) + 
      I(PHI_cereals*gdd_1to3) + I(PHI_cereals*prec_1to3) + I(PHI_cereals*gddHigh_1to3) + I(PHI_cereals*dd_1to3) +
      I(PHI_cereals*gdd_4to10) + I(PHI_cereals*prec_4to10) + I(PHI_cereals*gddHigh_4to10) + I(PHI_cereals*dd_4to10) +
      I(PHI_cereals*gdd_obs_1to3) + I(PHI_cereals*prec_obs_1to3) + I(PHI_cereals*gddHigh_obs_1to3) + I(PHI_cereals*dd_obs_1to3) + 
      I(PHI_cereals*gdd_obs_4to10) + I(PHI_cereals*prec_obs_4to10) + I(PHI_cereals*gddHigh_obs_4to10) + I(PHI_cereals*dd_obs_4to10) + 
      I(PHI_cereals*np_cereals_fm) + I(PHI_cereals*np_protein_fm) + I(PHI_cereals*np_oilseed_fm) + I(PHI_cereals*np_roots_fm) + I(PHI_cereals*np_corn_fm) + I(PHI_cereals*nw_fert_fm) +
      I(PHI_cereals*k_land_fm) + I(PHI_cereals*k_labor_fm) + I(PHI_cereals*k_capital_fm) + I(PHI_cereals*trend_fm) + I(PHI_cereals*trend2_fm) + 
      I(PHI_cereals*gdd_obs_fm) + I(PHI_cereals*prec_obs_fm) + I(PHI_cereals*gddHigh_obs_fm) + I(PHI_cereals*dd_obs_fm) + 
      I(PHI_cereals*gdd_1to3_fm) + I(PHI_cereals*prec_1to3_fm) + I(PHI_cereals*gddHigh_1to3_fm) + I(PHI_cereals*dd_1to3_fm) +
      I(PHI_cereals*gdd_4to10_fm) + I(PHI_cereals*prec_4to10_fm) + I(PHI_cereals*gddHigh_4to10_fm) + I(PHI_cereals*dd_4to10_fm) + 
      I(PHI_cereals*gdd_obs_1to3_fm) + I(PHI_cereals*prec_obs_1to3_fm) + I(PHI_cereals*gddHigh_obs_1to3_fm) + I(PHI_cereals*dd_obs_1to3_fm) + 
      I(PHI_cereals*gdd_obs_4to10_fm) + I(PHI_cereals*prec_obs_4to10_fm) + I(PHI_cereals*gddHigh_obs_4to10_fm) + I(PHI_cereals*dd_obs_4to10_fm) + 
      phi_cereals -1
    eqQQprotein <- qq_protein ~ PHI_protein + I(PHI_protein*np_cereals) + I(PHI_protein*np_protein) + I(PHI_protein*np_oilseed) + I(PHI_protein*np_roots) + I(PHI_protein*np_corn) + I(PHI_protein*nw_fert) + 
      I(PHI_protein*k_land) + I(PHI_protein*k_labor) + I(PHI_protein*k_capital) + I(PHI_protein*trend) + I(PHI_protein*trend2) +
      I(PHI_protein*gdd_obs) + I(PHI_protein*prec_obs) + I(PHI_protein*gddHigh_obs) + I(PHI_protein*dd_obs) + 
      I(PHI_protein*gdd_1to3) + I(PHI_protein*prec_1to3) + I(PHI_protein*gddHigh_1to3) + I(PHI_protein*dd_1to3) +
      I(PHI_protein*gdd_4to10) + I(PHI_protein*prec_4to10) + I(PHI_protein*gddHigh_4to10) + I(PHI_protein*dd_4to10) +
      I(PHI_protein*gdd_obs_1to3) + I(PHI_protein*prec_obs_1to3) + I(PHI_protein*gddHigh_obs_1to3) + I(PHI_protein*dd_obs_1to3) + 
      I(PHI_protein*gdd_obs_4to10) + I(PHI_protein*prec_obs_4to10) + I(PHI_protein*gddHigh_obs_4to10) + I(PHI_protein*dd_obs_4to10) + 
      I(PHI_protein*np_cereals_fm) + I(PHI_protein*np_protein_fm) + I(PHI_protein*np_oilseed_fm) + I(PHI_protein*np_roots_fm) + I(PHI_protein*np_corn_fm) + I(PHI_protein*nw_fert_fm) +
      I(PHI_protein*k_land_fm) + I(PHI_protein*k_labor_fm) + I(PHI_protein*k_capital_fm) + I(PHI_protein*trend_fm) + I(PHI_protein*trend2_fm) + 
      I(PHI_protein*gdd_obs_fm) + I(PHI_protein*prec_obs_fm) + I(PHI_protein*gddHigh_obs_fm) + I(PHI_protein*dd_obs_fm) + 
      I(PHI_protein*gdd_1to3_fm) + I(PHI_protein*prec_1to3_fm) + I(PHI_protein*gddHigh_1to3_fm) + I(PHI_protein*dd_1to3_fm) +
      I(PHI_protein*gdd_4to10_fm) + I(PHI_protein*prec_4to10_fm) + I(PHI_protein*gddHigh_4to10_fm) + I(PHI_protein*dd_4to10_fm) + 
      I(PHI_protein*gdd_obs_1to3_fm) + I(PHI_protein*prec_obs_1to3_fm) + I(PHI_protein*gddHigh_obs_1to3_fm) + I(PHI_protein*dd_obs_1to3_fm) + 
      I(PHI_protein*gdd_obs_4to10_fm) + I(PHI_protein*prec_obs_4to10_fm) + I(PHI_protein*gddHigh_obs_4to10_fm) + I(PHI_protein*dd_obs_4to10_fm) + 
      phi_protein -1
    eqQQoilseed <- qq_oilseed ~ PHI_oilseed + I(PHI_oilseed*np_cereals) + I(PHI_oilseed*np_protein) + I(PHI_oilseed*np_oilseed) + I(PHI_oilseed*np_roots) + I(PHI_oilseed*np_corn) + I(PHI_oilseed*nw_fert) + 
      I(PHI_oilseed*k_land) + I(PHI_oilseed*k_labor) + I(PHI_oilseed*k_capital) + I(PHI_oilseed*trend) + I(PHI_oilseed*trend2) +
      I(PHI_oilseed*gdd_obs) + I(PHI_oilseed*prec_obs) + I(PHI_oilseed*gddHigh_obs) + I(PHI_oilseed*dd_obs) + 
      I(PHI_oilseed*gdd_1to3) + I(PHI_oilseed*prec_1to3) + I(PHI_oilseed*gddHigh_1to3) + I(PHI_oilseed*dd_1to3) +
      I(PHI_oilseed*gdd_4to10) + I(PHI_oilseed*prec_4to10) + I(PHI_oilseed*gddHigh_4to10) + I(PHI_oilseed*dd_4to10) +
      I(PHI_oilseed*gdd_obs_1to3) + I(PHI_oilseed*prec_obs_1to3) + I(PHI_oilseed*gddHigh_obs_1to3) + I(PHI_oilseed*dd_obs_1to3) + 
      I(PHI_oilseed*gdd_obs_4to10) + I(PHI_oilseed*prec_obs_4to10) + I(PHI_oilseed*gddHigh_obs_4to10) + I(PHI_oilseed*dd_obs_4to10) + 
      I(PHI_oilseed*np_cereals_fm) + I(PHI_oilseed*np_protein_fm) + I(PHI_oilseed*np_oilseed_fm) + I(PHI_oilseed*np_roots_fm) + I(PHI_oilseed*np_corn_fm) + I(PHI_oilseed*nw_fert_fm) +
      I(PHI_oilseed*k_land_fm) + I(PHI_oilseed*k_labor_fm) + I(PHI_oilseed*k_capital_fm) + I(PHI_oilseed*trend_fm) + I(PHI_oilseed*trend2_fm) + 
      I(PHI_oilseed*gdd_obs_fm) + I(PHI_oilseed*prec_obs_fm) + I(PHI_oilseed*gddHigh_obs_fm) + I(PHI_oilseed*dd_obs_fm) + 
      I(PHI_oilseed*gdd_1to3_fm) + I(PHI_oilseed*prec_1to3_fm) + I(PHI_oilseed*gddHigh_1to3_fm) + I(PHI_oilseed*dd_1to3_fm) +
      I(PHI_oilseed*gdd_4to10_fm) + I(PHI_oilseed*prec_4to10_fm) + I(PHI_oilseed*gddHigh_4to10_fm) + I(PHI_oilseed*dd_4to10_fm) + 
      I(PHI_oilseed*gdd_obs_1to3_fm) + I(PHI_oilseed*prec_obs_1to3_fm) + I(PHI_oilseed*gddHigh_obs_1to3_fm) + I(PHI_oilseed*dd_obs_1to3_fm) + 
      I(PHI_oilseed*gdd_obs_4to10_fm) + I(PHI_oilseed*prec_obs_4to10_fm) + I(PHI_oilseed*gddHigh_obs_4to10_fm) + I(PHI_oilseed*dd_obs_4to10_fm) + 
      phi_oilseed -1
    eqQQroots <- qq_roots ~ PHI_roots + I(PHI_roots*np_cereals) + I(PHI_roots*np_protein) + I(PHI_roots*np_oilseed) + I(PHI_roots*np_roots) + I(PHI_roots*np_corn) + I(PHI_roots*nw_fert) + 
      I(PHI_roots*k_land) + I(PHI_roots*k_labor) + I(PHI_roots*k_capital) + I(PHI_roots*trend) + I(PHI_roots*trend2) +
      I(PHI_roots*gdd_obs) + I(PHI_roots*prec_obs) + I(PHI_roots*gddHigh_obs) + I(PHI_roots*dd_obs) + 
      I(PHI_roots*gdd_1to3) + I(PHI_roots*prec_1to3) + I(PHI_roots*gddHigh_1to3) + I(PHI_roots*dd_1to3) +
      I(PHI_roots*gdd_4to10) + I(PHI_roots*prec_4to10) + I(PHI_roots*gddHigh_4to10) + I(PHI_roots*dd_4to10) +
      I(PHI_roots*gdd_obs_1to3) + I(PHI_roots*prec_obs_1to3) + I(PHI_roots*gddHigh_obs_1to3) + I(PHI_roots*dd_obs_1to3) + 
      I(PHI_roots*gdd_obs_4to10) + I(PHI_roots*prec_obs_4to10) + I(PHI_roots*gddHigh_obs_4to10) + I(PHI_roots*dd_obs_4to10) + 
      I(PHI_roots*np_cereals_fm) + I(PHI_roots*np_protein_fm) + I(PHI_roots*np_oilseed_fm) + I(PHI_roots*np_roots_fm) + I(PHI_roots*np_corn_fm) + I(PHI_roots*nw_fert_fm) +
      I(PHI_roots*k_land_fm) + I(PHI_roots*k_labor_fm) + I(PHI_roots*k_capital_fm) + I(PHI_roots*trend_fm) + I(PHI_roots*trend2_fm) + 
      I(PHI_roots*gdd_obs_fm) + I(PHI_roots*prec_obs_fm) + I(PHI_roots*gddHigh_obs_fm) + I(PHI_roots*dd_obs_fm) + 
      I(PHI_roots*gdd_1to3_fm) + I(PHI_roots*prec_1to3_fm) + I(PHI_roots*gddHigh_1to3_fm) + I(PHI_roots*dd_1to3_fm) +
      I(PHI_roots*gdd_4to10_fm) + I(PHI_roots*prec_4to10_fm) + I(PHI_roots*gddHigh_4to10_fm) + I(PHI_roots*dd_4to10_fm) +
      I(PHI_roots*gdd_obs_1to3_fm) + I(PHI_roots*prec_obs_1to3_fm) + I(PHI_roots*gddHigh_obs_1to3_fm) + I(PHI_roots*dd_obs_1to3_fm) + 
      I(PHI_roots*gdd_obs_4to10_fm) + I(PHI_roots*prec_obs_4to10_fm) + I(PHI_roots*gddHigh_obs_4to10_fm) + I(PHI_roots*dd_obs_4to10_fm) + 
      phi_roots -1
    eqQQcorn <- qq_corn ~ PHI_corn + I(PHI_corn*np_cereals) + I(PHI_corn*np_protein) + I(PHI_corn*np_oilseed) + I(PHI_corn*np_roots) + I(PHI_corn*np_corn) + I(PHI_corn*nw_fert) + 
      I(PHI_corn*k_land) + I(PHI_corn*k_labor) + I(PHI_corn*k_capital) + I(PHI_corn*trend) + I(PHI_corn*trend2) +
      I(PHI_corn*gdd_obs) + I(PHI_corn*prec_obs) + I(PHI_corn*gddHigh_obs) + I(PHI_corn*dd_obs) + 
      I(PHI_corn*gdd_1to3) + I(PHI_corn*prec_1to3) + I(PHI_corn*gddHigh_1to3) + I(PHI_corn*dd_1to3) +
      I(PHI_corn*gdd_4to10) + I(PHI_corn*prec_4to10) + I(PHI_corn*gddHigh_4to10) + I(PHI_corn*dd_4to10) +
      I(PHI_corn*gdd_obs_1to3) + I(PHI_corn*prec_obs_1to3) + I(PHI_corn*gddHigh_obs_1to3) + I(PHI_corn*dd_obs_1to3) + 
      I(PHI_corn*gdd_obs_4to10) + I(PHI_corn*prec_obs_4to10) + I(PHI_corn*gddHigh_obs_4to10) + I(PHI_corn*dd_obs_4to10) + 
      I(PHI_corn*np_cereals_fm) + I(PHI_corn*np_protein_fm) + I(PHI_corn*np_oilseed_fm) + I(PHI_corn*np_roots_fm) + I(PHI_corn*np_corn_fm) + I(PHI_corn*nw_fert_fm) +
      I(PHI_corn*k_land_fm) + I(PHI_corn*k_labor_fm) + I(PHI_corn*k_capital_fm) + I(PHI_corn*trend_fm) + I(PHI_corn*trend2_fm) + 
      I(PHI_corn*gdd_obs_fm) + I(PHI_corn*prec_obs_fm) + I(PHI_corn*gddHigh_obs_fm) + I(PHI_corn*dd_obs_fm) + 
      I(PHI_corn*gdd_1to3_fm) + I(PHI_corn*prec_1to3_fm) + I(PHI_corn*gddHigh_1to3_fm) + I(PHI_corn*dd_1to3_fm) +
      I(PHI_corn*gdd_4to10_fm) + I(PHI_corn*prec_4to10_fm) + I(PHI_corn*gddHigh_4to10_fm) + I(PHI_corn*dd_4to10_fm) +
      I(PHI_corn*gdd_obs_1to3_fm) + I(PHI_corn*prec_obs_1to3_fm) + I(PHI_corn*gddHigh_obs_1to3_fm) + I(PHI_corn*dd_obs_1to3_fm) + 
      I(PHI_corn*gdd_obs_4to10_fm) + I(PHI_corn*prec_obs_4to10_fm) + I(PHI_corn*gddHigh_obs_4to10_fm) + I(PHI_corn*dd_obs_4to10_fm) + 
      phi_corn -1
    eqNXfert <- nx_fert ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
      k_land + k_labor + k_capital + trend + trend2 +
      gdd_obs + prec_obs + gddHigh_obs + dd_obs + 
      gdd_1to3 + prec_1to3 + gddHigh_1to3 + dd_1to3 +
      gdd_4to10 + prec_4to10 + gddHigh_4to10 + dd_4to10 +
      gdd_obs_1to3 + prec_obs_1to3 + gddHigh_obs_1to3 + dd_obs_1to3 + 
      gdd_obs_4to10 + prec_obs_4to10 + gddHigh_obs_4to10 + dd_obs_4to10 + 
      np_cereals_fm + np_protein_fm + np_oilseed_fm + np_roots_fm + np_corn_fm + nw_fert_fm +
      k_land_fm + k_labor_fm + k_capital_fm + trend_fm + trend2_fm + 
      gdd_obs_fm + prec_obs_fm + gddHigh_obs_fm + dd_obs_fm + 
      gdd_1to3_fm + prec_1to3_fm + gddHigh_1to3_fm + dd_1to3_fm +
      gdd_4to10_fm + prec_4to10_fm + gddHigh_4to10_fm + dd_4to10_fm +
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
    model_nonlinear <- systemfit( formula = system, method = "SUR", 
                                  data = df_farm_no2019, restrict.matrix = restrict, 
                                  maxit = 100 )
    
    # Calculate MAE (to compare with out-of-sample prediction)
    predicted_outcomes <- predict(model_nonlinear)
    caret::postResample(df_farm_no2019$qq_cereals,predicted_outcomes$QQcereals.pred)
    caret::postResample(df_farm_no2019$qq_protein,predicted_outcomes$QQprotein.pred)
    caret::postResample(df_farm_no2019$qq_oilseed,predicted_outcomes$QQoilseed.pred)
    caret::postResample(df_farm_no2019$qq_roots,predicted_outcomes$QQroots.pred)
    caret::postResample(df_farm_no2019$qq_corn,predicted_outcomes$QQcorn.pred)
    caret::postResample(df_farm_no2019$x_fert,-predicted_outcomes$NXfert.pred)
    
    
#-----------------------------------#
#### 4) Out-of-sample prediction ####
#-----------------------------------#    
    
    # ---------------------------------------------------- #
    # Predicted values for farms in the sample (1996-2018) #
    # -----------------------------------------------------#
    
    # coefficients:
    matrix(model_nonlinear$coefficients)
    matrix(model_nonlinear$eq[[1]]$coefficients)
    
    # data
    df_in_sample <- data.frame(x1=rep(1,nrow(df_farm_no2019)))
    df_in_sample$x2 <- df_farm_no2019$np_cereals
    df_in_sample$x3 <- df_farm_no2019$np_protein
    df_in_sample$x4 <- df_farm_no2019$np_oilseed
    df_in_sample$x5 <- df_farm_no2019$np_roots
    df_in_sample$x6 <- df_farm_no2019$np_corn
    df_in_sample$x7 <- df_farm_no2019$nw_fert
    df_in_sample$x8 <- df_farm_no2019$k_land
    df_in_sample$x9 <- df_farm_no2019$k_labor
    df_in_sample$x10 <- df_farm_no2019$k_capital
    df_in_sample$x11 <- df_farm_no2019$trend
    df_in_sample$x12 <- df_farm_no2019$trend2
    df_in_sample$x13 <- df_farm_no2019$gdd_obs
    df_in_sample$x14 <- df_farm_no2019$prec_obs
    df_in_sample$x15 <- df_farm_no2019$gddHigh_obs
    df_in_sample$x16 <- df_farm_no2019$dd_obs
    df_in_sample$x17 <- df_farm_no2019$gdd_1to3
    df_in_sample$x18 <- df_farm_no2019$prec_1to3
    df_in_sample$x19 <- df_farm_no2019$gddHigh_1to3
    df_in_sample$x20 <- df_farm_no2019$dd_1to3
    df_in_sample$x21 <- df_farm_no2019$gdd_4to10
    df_in_sample$x22 <- df_farm_no2019$prec_4to10
    df_in_sample$x23 <- df_farm_no2019$gddHigh_4to10
    df_in_sample$x24 <- df_farm_no2019$dd_4to10
    df_in_sample$x25 <- df_farm_no2019$gdd_obs_1to3
    df_in_sample$x26 <- df_farm_no2019$prec_obs_1to3
    df_in_sample$x27 <- df_farm_no2019$gddHigh_obs_1to3
    df_in_sample$x28 <- df_farm_no2019$dd_obs_1to3
    df_in_sample$x29 <- df_farm_no2019$gdd_obs_4to10
    df_in_sample$x30 <- df_farm_no2019$prec_obs_4to10
    df_in_sample$x31 <- df_farm_no2019$gddHigh_obs_4to10
    df_in_sample$x32 <- df_farm_no2019$dd_obs_4to10
    df_in_sample$x33 <- df_farm_no2019$np_cereals_fm
    df_in_sample$x34 <- df_farm_no2019$np_protein_fm
    df_in_sample$x35 <- df_farm_no2019$np_oilseed_fm
    df_in_sample$x36 <- df_farm_no2019$np_roots_fm
    df_in_sample$x37 <- df_farm_no2019$np_corn_fm
    df_in_sample$x38 <- df_farm_no2019$nw_fert_fm
    df_in_sample$x39 <- df_farm_no2019$k_land_fm
    df_in_sample$x40 <- df_farm_no2019$k_labor_fm
    df_in_sample$x41 <- df_farm_no2019$k_capital_fm
    df_in_sample$x42 <- df_farm_no2019$trend_fm
    df_in_sample$x43 <- df_farm_no2019$trend2_fm
    df_in_sample$x44 <- df_farm_no2019$gdd_obs_fm
    df_in_sample$x45 <- df_farm_no2019$prec_obs_fm
    df_in_sample$x46 <- df_farm_no2019$gddHigh_obs_fm
    df_in_sample$x47 <- df_farm_no2019$dd_obs_fm
    df_in_sample$x48 <- df_farm_no2019$gdd_1to3_fm
    df_in_sample$x49 <- df_farm_no2019$prec_1to3_fm
    df_in_sample$x50 <- df_farm_no2019$gddHigh_1to3_fm
    df_in_sample$x51 <- df_farm_no2019$dd_1to3_fm
    df_in_sample$x52 <- df_farm_no2019$gdd_4to10_fm
    df_in_sample$x53 <- df_farm_no2019$prec_4to10_fm
    df_in_sample$x54 <- df_farm_no2019$gddHigh_4to10_fm
    df_in_sample$x55 <- df_farm_no2019$dd_4to10_fm
    df_in_sample$x56 <- df_farm_no2019$gdd_obs_1to3_fm
    df_in_sample$x57 <- df_farm_no2019$prec_obs_1to3_fm
    df_in_sample$x58 <- df_farm_no2019$gddHigh_obs_1to3_fm
    df_in_sample$x59 <- df_farm_no2019$dd_obs_1to3_fm
    df_in_sample$x60 <- df_farm_no2019$gdd_obs_4to10_fm
    df_in_sample$x61 <- df_farm_no2019$prec_obs_4to10_fm
    df_in_sample$x62 <- df_farm_no2019$gddHigh_obs_4to10_fm
    df_in_sample$x63 <- df_farm_no2019$dd_obs_4to10_fm
    
    # For cereals
    df_in_sample_cereals <- df_in_sample * df_farm_no2019$PHI_cereals
    df_in_sample_cereals$x64 <- df_farm_no2019$phi_cereals
    df_in_sample_cereals$qq_cereals_pred <- t(crossprod(matrix(model_nonlinear$eq[[1]]$coefficients),
                                                        (t(df_in_sample_cereals[1:64]))))
    
    # For protein
    df_in_sample_protein <- df_in_sample * df_farm_no2019$PHI_protein
    df_in_sample_protein$x64 <- df_farm_no2019$phi_protein
    df_in_sample_protein$qq_protein_pred <- t(crossprod(matrix(model_nonlinear$eq[[2]]$coefficients),
                                                        (t(df_in_sample_protein[1:64]))))
    
    # For oilseed
    df_in_sample_oilseed <- df_in_sample * df_farm_no2019$PHI_oilseed
    df_in_sample_oilseed$x64 <- df_farm_no2019$phi_oilseed
    df_in_sample_oilseed$qq_oilseed_pred <- t(crossprod(matrix(model_nonlinear$eq[[3]]$coefficients),
                                                        (t(df_in_sample_oilseed[1:64]))))
    
    # For roots
    df_in_sample_roots <- df_in_sample * df_farm_no2019$PHI_roots
    df_in_sample_roots$x64 <- df_farm_no2019$phi_roots
    df_in_sample_roots$qq_roots_pred <- t(crossprod(matrix(model_nonlinear$eq[[4]]$coefficients),
                                                    (t(df_in_sample_roots[1:64]))))
    
    # For corn
    df_in_sample_corn <- df_in_sample * df_farm_no2019$PHI_corn
    df_in_sample_corn$x64 <- df_farm_no2019$phi_corn
    df_in_sample_corn$qq_corn_pred <- t(crossprod(matrix(model_nonlinear$eq[[5]]$coefficients),
                                                  (t(df_in_sample_corn[1:64]))))
    
    # For fert
    df_in_sample_fert <- df_in_sample 
    df_in_sample_fert$nx_fert_pred <- t(crossprod(matrix(model_nonlinear$eq[[6]]$coefficients),
                                                  (t(df_in_sample_fert[1:63]))))
    
    # ------------------------------------------------#
    # Predicted values for out-of-sample farms (2019) #
    # ------------------------------------------------#
    
    # Data
    df_out_of_sample <- data.frame(x1=rep(1,nrow(df_farm_2019)))
    df_out_of_sample$x2 <- df_farm_2019$np_cereals
    df_out_of_sample$x3 <- df_farm_2019$np_protein
    df_out_of_sample$x4 <- df_farm_2019$np_oilseed
    df_out_of_sample$x5 <- df_farm_2019$np_roots
    df_out_of_sample$x6 <- df_farm_2019$np_corn
    df_out_of_sample$x7 <- df_farm_2019$nw_fert
    df_out_of_sample$x8 <- df_farm_2019$k_land
    df_out_of_sample$x9 <- df_farm_2019$k_labor
    df_out_of_sample$x10 <- df_farm_2019$k_capital
    df_out_of_sample$x11 <- df_farm_2019$trend
    df_out_of_sample$x12 <- df_farm_2019$trend2
    df_out_of_sample$x13 <- df_farm_2019$gdd_obs
    df_out_of_sample$x14 <- df_farm_2019$prec_obs
    df_out_of_sample$x15 <- df_farm_2019$gddHigh_obs
    df_out_of_sample$x16 <- df_farm_2019$dd_obs
    df_out_of_sample$x17 <- df_farm_2019$gdd_1to3
    df_out_of_sample$x18 <- df_farm_2019$prec_1to3
    df_out_of_sample$x19 <- df_farm_2019$gddHigh_1to3
    df_out_of_sample$x20 <- df_farm_2019$dd_1to3
    df_out_of_sample$x21 <- df_farm_2019$gdd_4to10
    df_out_of_sample$x22 <- df_farm_2019$prec_4to10
    df_out_of_sample$x23 <- df_farm_2019$gddHigh_4to10
    df_out_of_sample$x24 <- df_farm_2019$dd_4to10
    df_out_of_sample$x25 <- df_farm_2019$gdd_obs_1to3
    df_out_of_sample$x26 <- df_farm_2019$prec_obs_1to3
    df_out_of_sample$x27 <- df_farm_2019$gddHigh_obs_1to3
    df_out_of_sample$x28 <- df_farm_2019$dd_obs_1to3
    df_out_of_sample$x29 <- df_farm_2019$gdd_obs_4to10
    df_out_of_sample$x30 <- df_farm_2019$prec_obs_4to10
    df_out_of_sample$x31 <- df_farm_2019$gddHigh_obs_4to10
    df_out_of_sample$x32 <- df_farm_2019$dd_obs_4to10
    df_out_of_sample$x33 <- df_farm_2019$np_cereals_fm
    df_out_of_sample$x34 <- df_farm_2019$np_protein_fm
    df_out_of_sample$x35 <- df_farm_2019$np_oilseed_fm
    df_out_of_sample$x36 <- df_farm_2019$np_roots_fm
    df_out_of_sample$x37 <- df_farm_2019$np_corn_fm
    df_out_of_sample$x38 <- df_farm_2019$nw_fert_fm
    df_out_of_sample$x39 <- df_farm_2019$k_land_fm
    df_out_of_sample$x40 <- df_farm_2019$k_labor_fm
    df_out_of_sample$x41 <- df_farm_2019$k_capital_fm
    df_out_of_sample$x42 <- df_farm_2019$trend_fm
    df_out_of_sample$x43 <- df_farm_2019$trend2_fm
    df_out_of_sample$x44 <- df_farm_2019$gdd_obs_fm
    df_out_of_sample$x45 <- df_farm_2019$prec_obs_fm
    df_out_of_sample$x46 <- df_farm_2019$gddHigh_obs_fm
    df_out_of_sample$x47 <- df_farm_2019$dd_obs_fm
    df_out_of_sample$x48 <- df_farm_2019$gdd_1to3_fm
    df_out_of_sample$x49 <- df_farm_2019$prec_1to3_fm
    df_out_of_sample$x50 <- df_farm_2019$gddHigh_1to3_fm
    df_out_of_sample$x51 <- df_farm_2019$dd_1to3_fm
    df_out_of_sample$x52 <- df_farm_2019$gdd_4to10_fm
    df_out_of_sample$x53 <- df_farm_2019$prec_4to10_fm
    df_out_of_sample$x54 <- df_farm_2019$gddHigh_4to10_fm
    df_out_of_sample$x55 <- df_farm_2019$dd_4to10_fm
    df_out_of_sample$x56 <- df_farm_2019$gdd_obs_1to3_fm
    df_out_of_sample$x57 <- df_farm_2019$prec_obs_1to3_fm
    df_out_of_sample$x58 <- df_farm_2019$gddHigh_obs_1to3_fm
    df_out_of_sample$x59 <- df_farm_2019$dd_obs_1to3_fm
    df_out_of_sample$x60 <- df_farm_2019$gdd_obs_4to10_fm
    df_out_of_sample$x61 <- df_farm_2019$prec_obs_4to10_fm
    df_out_of_sample$x62 <- df_farm_2019$gddHigh_obs_4to10_fm
    df_out_of_sample$x63 <- df_farm_2019$dd_obs_4to10_fm
    
    # For cereals
    df_out_of_sample_cereals <- df_out_of_sample * df_farm_2019$PHI_cereals
    df_out_of_sample_cereals$x64 <- df_farm_2019$phi_cereals
    df_out_of_sample_cereals$qq_cereals_pred <- t(crossprod(matrix(model_nonlinear$eq[[1]]$coefficients),
                                                            (t(df_out_of_sample_cereals[1:64]))))
    df_out_of_sample_cereals$key <- df_farm_2019$key
    
    # For protein
    df_out_of_sample_protein <- df_out_of_sample * df_farm_2019$PHI_protein
    df_out_of_sample_protein$x64 <- df_farm_2019$phi_protein
    df_out_of_sample_protein$qq_protein_pred <- t(crossprod(matrix(model_nonlinear$eq[[2]]$coefficients),
                                                            (t(df_out_of_sample_protein[1:64]))))
    df_out_of_sample_protein$key <- df_farm_2019$key
    
    # For oilseed
    df_out_of_sample_oilseed <- df_out_of_sample * df_farm_2019$PHI_oilseed
    df_out_of_sample_oilseed$x64 <- df_farm_2019$phi_oilseed
    df_out_of_sample_oilseed$qq_oilseed_pred <- t(crossprod(matrix(model_nonlinear$eq[[3]]$coefficients),
                                                            (t(df_out_of_sample_oilseed[1:64]))))
    df_out_of_sample_oilseed$key <- df_farm_2019$key
    
    # For roots
    df_out_of_sample_roots <- df_out_of_sample * df_farm_2019$PHI_roots
    df_out_of_sample_roots$x64 <- df_farm_2019$phi_roots
    df_out_of_sample_roots$qq_roots_pred <- t(crossprod(matrix(model_nonlinear$eq[[4]]$coefficients),
                                                        (t(df_out_of_sample_roots[1:64]))))
    df_out_of_sample_roots$key <- df_farm_2019$key
    
    # For corn
    df_out_of_sample_corn <- df_out_of_sample * df_farm_2019$PHI_corn
    df_out_of_sample_corn$x64 <- df_farm_2019$phi_corn
    df_out_of_sample_corn$qq_corn_pred <- t(crossprod(matrix(model_nonlinear$eq[[5]]$coefficients),
                                                      (t(df_out_of_sample_corn[1:64]))))
    df_out_of_sample_corn$key <- df_farm_2019$key
    
    # For fert
    df_out_of_sample_fert <- df_out_of_sample 
    df_out_of_sample_fert$nx_fert_pred <- t(crossprod(matrix(model_nonlinear$eq[[6]]$coefficients),
                                                      (t(df_out_of_sample_fert[1:63]))))
    df_out_of_sample_fert$key <- df_farm_2019$key
    
    
    # ----------------------------------------------------#
    # Compare out-of-sample predicted vs. observed values #
    # ----------------------------------------------------#  
    
    comp <- left_join(df_farm_2019,df_out_of_sample_cereals)
    
    observed <- df_farm_2019[c("qq_cereals", "qq_protein", "qq_oilseed", "qq_roots", "qq_corn", "nx_fert")]
    predicted <- df_out_of_sample_cereals[c("qq_cereals_pred", "key")] %>% 
      full_join(df_out_of_sample_protein[c("qq_protein_pred", "key")]) %>% 
      full_join(df_out_of_sample_oilseed[c("qq_oilseed_pred", "key")]) %>% 
      full_join(df_out_of_sample_roots[c("qq_roots_pred", "key")]) %>% 
      full_join(df_out_of_sample_corn[c("qq_corn_pred", "key")]) %>% 
      full_join(df_out_of_sample_fert[c("nx_fert_pred", "key")])
    
    # correlations                               
    cor_cer <- cor(observed$qq_cereals,predicted$qq_cereals_pred)   
    cor_prot <- cor(observed$qq_protein,predicted$qq_protein_pred)   
    cor_oil <- cor(observed$qq_oilseed,predicted$qq_oilseed_pred)   
    cor_roots <- cor(observed$qq_roots,predicted$qq_roots_pred)   
    cor_corn <- cor(observed$qq_corn,predicted$qq_corn_pred)   
    cor_fert <- cor(observed$nx_fert,predicted$nx_fert_pred)   
    
    # RMSE and MAE
    summary(df_farm_2019[c("qq_cereals","qq_protein","qq_oilseed", "qq_roots", "qq_corn", "nx_fert")])
    caret_cer <- caret::postResample(observed$qq_cereals,predicted$qq_cereals_pred)   
    caret_prot <- caret::postResample(observed$qq_protein,predicted$qq_protein_pred)   
    caret_oil <- caret::postResample(observed$qq_oilseed,predicted$qq_oilseed_pred)   
    caret_roots <- caret::postResample(observed$qq_roots,predicted$qq_roots_pred)   
    caret_corn <- caret::postResample(observed$qq_corn,predicted$qq_corn_pred)   
    caret_fert <- caret::postResample(observed$nx_fert,predicted$nx_fert_pred) 
    
