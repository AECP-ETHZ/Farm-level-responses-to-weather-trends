#-------------------------------------------------------#
#                                                       #
# This program estimates first- and second-stage        #
# regressions and simulates the outcomes for the        #
# 2018-drought, for the model with curvature imposed.   #
#                                                       #
# Note: The farm-level data are confidential and cannot # 
# be loaded (line 17). The results at the sample mean   #
# are saved in line 1234 and the figure in Section 5    # 
# can be created.                                       #
#                                                       #
#-------------------------------------------------------#

library(dplyr)
library(systemfit)
library(ggplot2)
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


#-------------------------------------------------#
# Estimate and store probit results for all crops #
#-------------------------------------------------#

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

#De-mean the data so the regression can converge
vars <- c("qq_cereals", "qq_protein", "qq_oilseed", "qq_roots", "qq_corn", "x_fert", "nx_fert",
          "np_cereals", "np_protein", "np_oilseed", "np_roots", "np_corn", "nw_fert", 
          "k_land", "k_labor", "k_capital", "trend", "trend2",
          "gdd_obs", "prec_obs", "gddHigh_obs", "dd_obs", 
          "gdd_1to3", "prec_1to3", "gddHigh_1to3", "dd_1to3", 
          "gdd_4to10", "prec_4to10", "gddHigh_4to10", "dd_4to10", 
          "gdd_obs_1to3", "prec_obs_1to3", "gddHigh_obs_1to3", "dd_obs_1to3", 
          "gdd_obs_4to10", "prec_obs_4to10", "gddHigh_obs_4to10", "dd_obs_4to10", 
          "np_cereals_fm", "np_protein_fm", "np_oilseed_fm", "np_roots_fm", "np_corn_fm", "nw_fert_fm", 
          "k_land_fm", "k_labor_fm", "k_capital_fm", "trend_fm", "trend2_fm", 
          "gdd_obs_fm", "prec_obs_fm", "gddHigh_obs_fm", "dd_obs_fm", 
          "gdd_1to3_fm", "prec_1to3_fm", "gddHigh_1to3_fm", "dd_1to3_fm", 
          "gdd_4to10_fm", "prec_4to10_fm", "gddHigh_4to10_fm", "dd_4to10_fm", 
          "gdd_obs_1to3_fm", "prec_obs_1to3_fm", "gddHigh_obs_1to3_fm", "dd_obs_1to3_fm", 
          "gdd_obs_4to10_fm", "prec_obs_4to10_fm", "gddHigh_obs_4to10_fm", "dd_obs_4to10_fm", 
          "PHI_cereals", "phi_cereals", "PHI_protein", "phi_protein", "PHI_oilseed", "phi_oilseed", "PHI_roots", "phi_roots", "PHI_corn", "phi_corn")

temp_fn <- function(x) x/mean(x) 
df_farm_dem <- df_farm %>% mutate(across(vars, temp_fn))


#Define structural equations
eqns <-
  list(qq_cereals ~ (b1*PHI_cereals + 
                       (c11*c11)*np_cereals*PHI_cereals +
                       (c11*c12)*np_protein*PHI_cereals +
                       (c11*c13)*np_oilseed*PHI_cereals +
                       (c11*c14)*np_roots*PHI_cereals +
                       (c11*c15)*np_corn*PHI_cereals +
                       (c11*c16)*nw_fert*PHI_cereals +
                       k11*k_land*PHI_cereals +
                       k12*k_labor*PHI_cereals + 
                       k13*k_capital*PHI_cereals +
                       t1*trend*PHI_cereals +
                       tsq1*trend2*PHI_cereals +
                       w11*gdd_obs*PHI_cereals + 
                       w12*prec_obs*PHI_cereals +
                       w13*gddHigh_obs*PHI_cereals +
                       w14*dd_obs*PHI_cereals +
                       w15*gdd_1to3*PHI_cereals + 
                       w16*prec_1to3*PHI_cereals +
                       w17*gddHigh_1to3*PHI_cereals +
                       w18*dd_1to3*PHI_cereals +
                       w19*gdd_4to10*PHI_cereals + 
                       w110*prec_4to10*PHI_cereals +
                       w111*gddHigh_4to10*PHI_cereals +
                       w112*dd_4to10*PHI_cereals +
                       m11*np_cereals_fm*PHI_cereals +
                       m12*np_protein_fm*PHI_cereals +
                       m13*np_oilseed_fm*PHI_cereals +
                       m14*np_roots_fm*PHI_cereals +
                       m15*np_corn_fm*PHI_cereals +
                       m16*nw_fert_fm*PHI_cereals +
                       m17*k_land_fm*PHI_cereals +
                       m18*k_labor_fm*PHI_cereals +
                       m19*k_capital_fm*PHI_cereals + 
                       t1_fm*trend_fm*PHI_cereals +
                       tsq1_fm*trend2_fm*PHI_cereals +
                       w11_fm*gdd_obs_fm*PHI_cereals + 
                       w12_fm*prec_obs_fm*PHI_cereals +
                       w13_fm*gddHigh_obs_fm*PHI_cereals +
                       w14_fm*dd_obs_fm*PHI_cereals +
                       w15_fm*gdd_1to3_fm*PHI_cereals + 
                       w16_fm*prec_1to3_fm*PHI_cereals +
                       w17_fm*gddHigh_1to3_fm*PHI_cereals +
                       w18_fm*dd_1to3_fm*PHI_cereals +
                       w19_fm*gdd_4to10_fm*PHI_cereals + 
                       w110_fm*prec_4to10_fm*PHI_cereals +
                       w111_fm*gddHigh_4to10_fm*PHI_cereals +
                       w112_fm*dd_4to10_fm*PHI_cereals) +
         phi11*phi_cereals,
       qq_protein ~  (b2*PHI_protein +  
                        (c12*c11)*np_cereals*PHI_protein +
                        (c12*c12+c22*c22)*np_protein*PHI_protein +
                        (c12*c13+c22*c23)*np_oilseed*PHI_protein +
                        (c12*c14+c22*c24)*np_roots*PHI_protein +
                        (c12*c15+c22*c25)*np_corn*PHI_protein + 
                        (c12*c16+c22*c26)*nw_fert*PHI_protein +
                        k21*k_land*PHI_protein +
                        k22*k_labor*PHI_protein + 
                        k23*k_capital*PHI_protein +
                        t2*trend*PHI_protein +
                        tsq2*trend2*PHI_protein +
                        w21*gdd_obs*PHI_protein + 
                        w22*prec_obs*PHI_protein +
                        w23*gddHigh_obs*PHI_protein +
                        w24*dd_obs*PHI_protein +
                        w25*gdd_1to3*PHI_protein + 
                        w26*prec_1to3*PHI_protein +
                        w27*gddHigh_1to3*PHI_protein +
                        w28*dd_1to3*PHI_protein +
                        w29*gdd_4to10*PHI_protein + 
                        w210*prec_4to10*PHI_protein +
                        w211*gddHigh_4to10*PHI_protein +
                        w212*dd_4to10*PHI_protein +
                        m21*np_cereals_fm*PHI_protein +
                        m22*np_protein_fm*PHI_protein +
                        m23*np_oilseed_fm*PHI_protein +
                        m24*np_roots_fm*PHI_protein +
                        m25*np_corn_fm*PHI_protein +
                        m26*nw_fert_fm*PHI_protein +
                        m27*k_land_fm*PHI_protein +
                        m28*k_labor_fm*PHI_protein +
                        m29*k_capital_fm*PHI_protein + 
                        t2_fm*trend_fm*PHI_protein +
                        tsq2_fm*trend2_fm*PHI_protein +
                        w21_fm*gdd_obs_fm*PHI_protein + 
                        w22_fm*prec_obs_fm*PHI_protein +
                        w23_fm*gddHigh_obs_fm*PHI_protein +
                        w24_fm*dd_obs_fm*PHI_protein +
                        w25_fm*gdd_1to3_fm*PHI_protein + 
                        w26_fm*prec_1to3_fm*PHI_protein +
                        w27_fm*gddHigh_1to3_fm*PHI_protein +
                        w28_fm*dd_1to3_fm*PHI_protein +
                        w29_fm*gdd_4to10_fm*PHI_protein + 
                        w210_fm*prec_4to10_fm*PHI_protein +
                        w211_fm*gddHigh_4to10_fm*PHI_protein +
                        w212_fm*dd_4to10_fm*PHI_protein) +
         phi22*phi_protein,
       qq_oilseed   ~ (b3*PHI_oilseed + 
                         (c13*c11)*np_cereals*PHI_oilseed +
                         (c13*c12+c23*c22)*np_protein*PHI_oilseed +
                         (c13*c13+c23*c23+c33*c33)*np_oilseed*PHI_oilseed +
                         (c13*c14+c23*c24+c33*c34)*np_roots*PHI_oilseed +
                         (c13*c15+c23*c25+c33*c35)*np_corn*PHI_oilseed + 
                         (c13*c16+c23*c26+c33*c36)*nw_fert*PHI_oilseed +
                         k31*k_land*PHI_oilseed +
                         k32*k_labor*PHI_oilseed + 
                         k33*k_capital*PHI_oilseed +
                         t3*trend*PHI_oilseed +
                         tsq3*trend2*PHI_oilseed +
                         w31*gdd_obs*PHI_oilseed + 
                         w32*prec_obs*PHI_oilseed +
                         w33*gddHigh_obs*PHI_oilseed +
                         w34*dd_obs*PHI_oilseed +
                         w35*gdd_1to3*PHI_oilseed + 
                         w36*prec_1to3*PHI_oilseed +
                         w37*gddHigh_1to3*PHI_oilseed +
                         w38*dd_1to3*PHI_oilseed +
                         w39*gdd_4to10*PHI_oilseed + 
                         w310*prec_4to10*PHI_oilseed +
                         w311*gddHigh_4to10*PHI_oilseed +
                         w312*dd_4to10*PHI_oilseed +
                         m31*np_cereals_fm*PHI_oilseed +
                         m32*np_protein_fm*PHI_oilseed +
                         m33*np_oilseed_fm*PHI_oilseed +
                         m34*np_roots_fm*PHI_oilseed +
                         m35*np_corn_fm*PHI_oilseed +
                         m36*nw_fert_fm*PHI_oilseed +
                         m37*k_land_fm*PHI_oilseed +
                         m38*k_labor_fm*PHI_oilseed +
                         m39*k_capital_fm*PHI_oilseed + 
                         t3_fm*trend_fm*PHI_oilseed +
                         tsq3_fm*trend2_fm*PHI_oilseed +
                         w31_fm*gdd_obs_fm*PHI_oilseed + 
                         w32_fm*prec_obs_fm*PHI_oilseed +
                         w33_fm*gddHigh_obs_fm*PHI_oilseed +
                         w34_fm*dd_obs_fm*PHI_oilseed +
                         w35_fm*gdd_1to3_fm*PHI_oilseed + 
                         w36_fm*prec_1to3_fm*PHI_oilseed +
                         w37_fm*gddHigh_1to3_fm*PHI_oilseed +
                         w38_fm*dd_1to3_fm*PHI_oilseed +
                         w39_fm*gdd_4to10_fm*PHI_oilseed + 
                         w310_fm*prec_4to10_fm*PHI_oilseed +
                         w311_fm*gddHigh_4to10_fm*PHI_oilseed +
                         w312_fm*dd_4to10_fm*PHI_oilseed) +
         phi33*phi_oilseed,
       qq_roots ~   (b4*PHI_roots +
                       (c14*c11)*np_cereals*PHI_roots +
                       (c14*c12+c24*c22)*np_protein*PHI_roots +
                       (c14*c13+c24*c23+c34*c33)*np_oilseed*PHI_roots +
                       (c14*c14+c24*c24+c34*c34+c44*c44)*np_roots*PHI_roots +
                       (c14*c15+c24*c25+c34*c35+c44*c45)*np_corn*PHI_roots + 
                       (c14*c16+c24*c26+c34*c36+c44*c46)*nw_fert*PHI_roots +
                       k41*k_land*PHI_roots +
                       k42*k_labor*PHI_roots + 
                       k43*k_capital*PHI_roots +
                       t4*trend*PHI_roots +
                       tsq4*trend2*PHI_roots +
                       w41*gdd_obs*PHI_roots + 
                       w42*prec_obs*PHI_roots +
                       w43*gddHigh_obs*PHI_roots +
                       w44*dd_obs*PHI_roots +
                       w45*gdd_1to3*PHI_roots + 
                       w46*prec_1to3*PHI_roots +
                       w47*gddHigh_1to3*PHI_roots +
                       w48*dd_1to3*PHI_roots +
                       w49*gdd_4to10*PHI_roots + 
                       w410*prec_4to10*PHI_roots +
                       w411*gddHigh_4to10*PHI_roots +
                       w412*dd_4to10*PHI_roots +
                       m41*np_cereals_fm*PHI_roots +
                       m42*np_protein_fm*PHI_roots +
                       m43*np_oilseed_fm*PHI_roots +
                       m44*np_roots_fm*PHI_roots +
                       m45*np_corn_fm*PHI_roots +
                       m46*nw_fert_fm*PHI_roots +
                       m47*k_land_fm*PHI_roots +
                       m48*k_labor_fm*PHI_roots +
                       m49*k_capital_fm*PHI_roots + 
                       t4_fm*trend_fm*PHI_roots +
                       tsq4_fm*trend2_fm*PHI_roots +
                       w41_fm*gdd_obs_fm*PHI_roots + 
                       w42_fm*prec_obs_fm*PHI_roots +
                       w43_fm*gddHigh_obs_fm*PHI_roots +
                       w44_fm*dd_obs_fm*PHI_roots +
                       w45_fm*gdd_1to3_fm*PHI_roots + 
                       w46_fm*prec_1to3_fm*PHI_roots +
                       w47_fm*gddHigh_1to3_fm*PHI_roots +
                       w48_fm*dd_1to3_fm*PHI_roots +
                       w49_fm*gdd_4to10_fm*PHI_roots + 
                       w410_fm*prec_4to10_fm*PHI_roots +
                       w411_fm*gddHigh_4to10_fm*PHI_roots +
                       w412_fm*dd_4to10_fm*PHI_roots) +
         phi44*phi_roots,
       qq_corn ~  (b5*PHI_corn +
                        (c15*c11)*np_cereals*PHI_corn +
                        (c15*c12+c25*c22)*np_protein*PHI_corn +
                        (c15*c13+c25*c23+c35*c33)*np_oilseed*PHI_corn +
                        (c15*c14+c25*c24+c35*c34+c45*c44)*np_roots*PHI_corn +
                        (c15*c15+c25*c25+c35*c35+c45*c45)*np_corn*PHI_corn + 
                        (c15*c16+c25*c26+c35*c36+c45*c46)*nw_fert*PHI_corn +		
                        k51*k_land*PHI_corn +
                        k52*k_labor*PHI_corn + 
                        k53*k_capital*PHI_corn +
                        t5*trend*PHI_corn +
                        tsq5*trend2*PHI_corn +
                        w51*gdd_obs*PHI_corn + 
                        w52*prec_obs*PHI_corn +
                        w53*gddHigh_obs*PHI_corn +
                        w54*dd_obs*PHI_corn +
                        w55*gdd_1to3*PHI_corn + 
                        w56*prec_1to3*PHI_corn +
                        w57*gddHigh_1to3*PHI_corn +
                        w58*dd_1to3*PHI_corn +
                        w59*gdd_4to10*PHI_corn + 
                        w510*prec_4to10*PHI_corn +
                        w511*gddHigh_4to10*PHI_corn +
                        w512*dd_4to10*PHI_corn +
                        m51*np_cereals_fm*PHI_corn +
                        m52*np_protein_fm*PHI_corn +
                        m53*np_oilseed_fm*PHI_corn +
                        m54*np_roots_fm*PHI_corn +
                        m55*np_corn_fm*PHI_corn +
                        m56*nw_fert_fm*PHI_corn +
                        m57*k_land_fm*PHI_corn +
                        m58*k_labor_fm*PHI_corn +
                        m59*k_capital_fm*PHI_corn + 
                        t5_fm*trend_fm*PHI_corn +
                        tsq5_fm*trend2_fm*PHI_corn +
                        w51_fm*gdd_obs_fm*PHI_corn + 
                        w52_fm*prec_obs_fm*PHI_corn +
                        w53_fm*gddHigh_obs_fm*PHI_corn +
                        w54_fm*dd_obs_fm*PHI_corn +
                        w55_fm*gdd_1to3_fm*PHI_corn + 
                        w56_fm*prec_1to3_fm*PHI_corn +
                        w57_fm*gddHigh_1to3_fm*PHI_corn +
                        w58_fm*dd_1to3_fm*PHI_corn +
                        w59_fm*gdd_4to10_fm*PHI_corn + 
                        w510_fm*prec_4to10_fm*PHI_corn +
                        w511_fm*gddHigh_4to10_fm*PHI_corn +
                        w512_fm*dd_4to10_fm*PHI_corn) +
         phi55*phi_corn,
       nx_fert ~     b6 +
         (c16*c11)*np_cereals +
         (c16*c12+c26*c22)*np_protein +
         (c16*c13+c26*c23+c36*c33)*np_oilseed +
         (c16*c14+c26*c24+c36*c34+c46*c44)*np_roots +
         (c16*c15+c26*c25+c36*c35+c46*c45)*np_corn + 
         (c16*c16+c26*c26+c36*c36+c46*c46)*nw_fert +		
         k61*k_land +
         k62*k_labor + 
         k63*k_capital +
         t6*trend +
         tsq6*trend2 +
         w61*gdd_obs + 
         w62*prec_obs +
         w63*gddHigh_obs +
         w64*dd_obs +
         w65*gdd_1to3 + 
         w66*prec_1to3 +
         w67*gddHigh_1to3 +
         w68*dd_1to3 +
         w69*gdd_4to10 + 
         w610*prec_4to10 +
         w611*gddHigh_4to10 +
         w612*dd_4to10 +
         m61*np_cereals_fm +
         m62*np_protein_fm +
         m63*np_oilseed_fm +
         m64*np_roots_fm +
         m65*np_corn_fm +
         m66*nw_fert_fm +
         m67*k_land_fm +
         m68*k_labor_fm +
         m69*k_capital_fm + 
         t6_fm*trend_fm +
         tsq6_fm*trend2_fm +
         w61_fm*gdd_obs_fm + 
         w62_fm*prec_obs_fm +
         w63_fm*gddHigh_obs_fm +
         w64_fm*dd_obs_fm +
         w65_fm*gdd_1to3_fm + 
         w66_fm*prec_1to3_fm +
         w67_fm*gddHigh_1to3_fm +
         w68_fm*dd_1to3_fm +
         w69_fm*gdd_4to10_fm + 
         w610_fm*prec_4to10_fm +
         w611_fm*gddHigh_4to10_fm +
         w612_fm*dd_4to10_fm)


# Starting values
strtvls <- c(
          b1 = 4.744e+00,
          c11 = 6.026e-01,
          c12 = 3.133e-02,
          c13 = 2.506e-01,
          c14 = 4.535e-02,
          c15 = 1.821e-04,
          c16 = 1.714e-01,
          b2 = 2.370e+00,
          c22 = 2.649e-01,
          c23 = 5.500e-01,
          c24 = 2.463e-01,
          c25 = 6.245e-03,
          c26 = -2.012e-01,
          b3 = -7.524e+00,
          c33 = 8.229e-01,
          c34 = -2.346e-01,
          c35 = -4.137e-03,
          c36 = -1.313e-01,
          b4 = -3.511e+00,
          c44 = 1.060e-01,
          c45 = -1.791e-02,
          c46 = -4.182e-02,
          b5 = -2.147e+01,
          b6 = -1.637e+00,
          k11 = 8.506e-01,
          k12 = 1.971e-02,
          k13 = 2.880e-02,
          k21 = 4.423e-01,
          k22 = -1.504e-01,
          k23 = -3.305e-03,
          k31 = 8.172e-01,
          k32 = -1.141e-01,
          k33 = -7.944e-03,
          k41 = 2.216e-01,
          k42 = -2.688e-02,
          k43 = 9.322e-02,
          k51 = -4.386e-01,
          k52 = -6.078e-03,
          k53 = -1.145e-01,
          k61 = 9.240e-01,
          k62 = 5.283e-02,
          k63 = 1.645e-02,
          t1 = 5.591e-01,
          t2 = 7.405e-01,
          t3 = 1.764e+00,
          t4 = -4.822e-01,
          t5 = 1.258e+00,
          t6 = 1.471e+00,
          tsq1 = -2.616e-01,
          tsq2 = -3.847e-01,
          tsq3 = -8.903e-01,
          tsq4 = 4.978e-01,
          tsq5 = 6.851e-03,
          tsq6 = -9.430e-01,
          m11 = -3.285e-01,
          m12 = -4.149e-02,
          m13 = 6.116e-01,
          m14 = -2.290e-01,
          m15 = 7.683e-04,
          m16 = -1.566e+00,
          m17 = 1.349e-01,
          m18 = -9.083e-02,
          m19 = 1.771e-01,
          m21 = -3.537e-01,
          m22 = 1.381e-01,
          m23 = -6.508e-01,
          m24 = -6.999e-01,
          m25 = -8.149e-03,
          m26 = -2.438e+00,
          m27 = -1.030e-01,
          m28 = 1.211e-01,
          m29 = 1.302e-01,
          m31 = 1.282e-01,
          m32 = -4.534e-01,
          m33 = 1.160e+00,
          m34 = -4.210e-02,
          m35 = 5.258e-03,
          m36 = 1.193e+00,
          m37 = -1.690e-01,
          m38 = 9.611e-02,
          m39 = 4.714e-02,
          m41 = -5.264e-01,
          m42 = 2.478e-01,
          m43 = -3.216e-01,
          m44 = 2.055e-01,
          m45 = -4.556e-03,
          m46 = 1.057e+00,
          m47 = 2.963e-02,
          m48 = 8.101e-02,
          m49 = 1.212e-01,
          m51 = 2.782e+00,
          m52 = 2.398e-04,
          m53 = -6.510e-01,
          m54 = 6.485e-01,
          m55 = 1.062e-02,
          m56 = 4.782e+00,
          m57 = 1.025e+00,
          m58 = -1.659e-01,
          m59 = 1.113e-01,
          m61 = 3.381e-01,
          m62 = -7.972e-02,
          m63 = 2.413e-01,
          m64 = -3.678e-03,
          m65 = 1.031e-03,
          m66 = -4.943e-03,
          m67 = 8.042e-02,
          m68 = -4.373e-02,
          m69 = -1.008e-02,
          w11 = -5.737e-01,
          w12 = -1.041e-01,
          w13 = 1.269e-02,
          w14 = -4.374e-01,
          w15 = -2.429e-01,
          w16 = 1.799e-01,
          w17 = -9.369e-03,
          w18 = 3.181e-01,
          w19 = -2.698e+00,
          w110 = -7.678e-02,
          w111 = 5.092e-02,
          w112 = 5.239e-01,
          w21 = -1.146e+00,
          w22 = -1.995e-01,
          w23 = 6.001e-02,
          w24 = -9.992e-01,
          w25 = -2.017e+00,
          w26 = 2.651e-01,
          w27 = 3.794e-02,
          w28 = 4.538e-01,
          w29 = -5.314e+00,
          w210 = 2.795e+00,
          w211 = 2.413e-01,
          w212 = 1.594e+00,
          w31 = -1.071e+00,
          w32 = -3.907e-01,
          w33 = 2.911e-02,
          w34 = -7.920e-01,
          w35 = -1.453e+00,
          w36 = 7.968e-03,
          w37 = 2.429e-02,
          w38 = -5.488e-01,
          w39 = -3.730e+00,
          w310 = -2.264e+00,
          w311 = 2.021e-02,
          w312 = -3.697e+00,
          w41 = 1.993e-01,
          w42 = -3.278e-02,
          w43 = -7.838e-02,
          w44 = 1.265e-01,
          w45 = 1.625e+00,
          w46 = 4.335e-01,
          w47 = -1.508e-01,
          w48 = 1.729e-01,
          w49 = 6.157e+00,
          w410 = 1.079e+00,
          w411 = -4.365e-01,
          w412 = 3.608e-01,
          w51 = -1.009e+00,
          w52 = -4.142e-01,
          w53 = -9.755e-02,
          w54 = -7.066e-01,
          w55 = -5.234e-01,
          w56 = -9.629e-01,
          w57 = -2.581e-01,
          w58 = -2.359e+00,
          w59 = -1.210e+00,
          w510 = -3.034e+00,
          w511 = -4.900e-01,
          w512 = -1.239e+01,
          w61 = -1.469e+00,
          w62 = 5.296e-02,
          w63 = 5.357e-02,
          w64 = 1.269e+00,
          w65 = -2.039e+00,
          w66 = -1.499e-01,
          w67 = -8.245e-02,
          w68 = 4.067e-02,
          w69 = -3.168e-01,
          w610 = -3.940e-01,
          w611 = -6.995e-02,
          w612 = 2.666e-01,
          phi11 = -2.183e-03,
          phi22 = -4.959e-01,
          phi33 = -9.127e-02,
          phi44 = -2.256e-01,
          phi55 = -5.720e-01,
          t1_fm = -1.390e-01,
          t2_fm = 4.151e-01,
          t3_fm = -1.513e+00,
          t4_fm = 3.726e-01,
          t5_fm = 8.685e-01,
          t6_fm = -5.762e-01,
          tsq1_fm = -1.736e-01,
          tsq2_fm = -4.534e-01,
          tsq3_fm = 7.420e-01,
          tsq4_fm = 1.914e-01,
          tsq5_fm = -3.968e-01,
          tsq6_fm = 3.613e-01,
          w11_fm = -4.496e+00,
          w12_fm = 8.577e-02,
          w13_fm = 2.534e-01,
          w14_fm = 2.011e-01,
          w15_fm = 1.905e+00,
          w16_fm = -6.720e-01,
          w17_fm = -1.322e-02,
          w18_fm = -2.990e+00,
          w19_fm = 4.685e+00,
          w110_fm = 3.109e-01,
          w111_fm = -8.961e-02,
          w112_fm = -9.329e-02,
          w21_fm = -6.051e+00,
          w22_fm = -2.119e+00,
          w23_fm = 1.082e-01,
          w24_fm = 4.192e+00,
          w25_fm = 1.762e+00,
          w26_fm = 6.114e-01,
          w27_fm = 5.651e-01,
          w28_fm = -4.839e+00,
          w29_fm = 1.139e+01,
          w210_fm = -7.979e-01,
          w211_fm = -9.319e-01,
          w212_fm = 1.838e+00,
          w31_fm = -9.971e-01,
          w32_fm = 2.136e+00,
          w33_fm = 1.551e-01,
          w34_fm = 2.807e+00,
          w35_fm = 5.321e+00,
          w36_fm = -8.559e-01,
          w37_fm = -2.256e-01,
          w38_fm = -1.902e+00,
          w39_fm = 2.518e+00,
          w310_fm = 2.041e+00,
          w311_fm = -3.737e-01,
          w312_fm = 7.639e+00,
          w41_fm = -2.341e-01,
          w42_fm = 6.146e-01,
          w43_fm = -3.271e-02,
          w44_fm = 2.645e+00,
          w45_fm = -4.284e+00,
          w46_fm = -1.399e+00,
          w47_fm = 8.840e-02,
          w48_fm = -2.255e+00,
          w49_fm = -3.481e-01,
          w410_fm = -1.055e+00,
          w411_fm = -1.978e-01,
          w412_fm = -1.842e-01,
          w51_fm = 1.118e+01,
          w52_fm = -3.170e-01,
          w53_fm = -9.601e-01,
          w54_fm = 4.075e+00,
          w55_fm = -1.602e+01,
          w56_fm = 1.826e+00,
          w57_fm = 1.457e+00,
          w58_fm = -5.088e+00,
          w59_fm = 1.145e+01,
          w510_fm = 4.361e+00,
          w511_fm = -1.073e+00,
          w512_fm = 2.578e+01,
          w61_fm = -2.735e+00,
          w62_fm = 1.756e+00,
          w63_fm = 5.111e-02,
          w64_fm = 2.907e+00,
          w65_fm = 2.011e+00,
          w66_fm = -1.781e+00,
          w67_fm = 9.353e-02,
          w68_fm = -6.691e+00,
          w69_fm = 3.138e+00,
          w610_fm = 1.207e+00,
          w611_fm = -2.754e-02,
          w612_fm = 3.767e+00)

## Regression with nonlinear SUR 
system.time(model_linear <- nlsur(eqns = eqns, data = df_farm_dem, startvalues = strtvls, type = "FGNLS",
                                  trace = TRUE, eps = 1e-03, tol = 1e-07, maxiter=1000))



# Extract data used in this estimation for simulation later

dat_str_all <- df_farm %>% 
  select(c(np_cereals,np_protein,np_oilseed,np_roots,np_corn,nw_fert,
           k_land,k_labor,k_capital,trend,trend2,
           gdd_obs,prec_obs,gddHigh_obs,dd_obs,
           gdd_1to3,prec_1to3,gddHigh_1to3,dd_1to3,
           gdd_4to10,prec_4to10,gddHigh_4to10,dd_4to10,
           np_cereals_fm,np_protein_fm,np_oilseed_fm,np_roots_fm,np_corn_fm,nw_fert_fm,
           k_land_fm,k_labor_fm,k_capital_fm,trend_fm,trend2_fm,
           gdd_obs_fm,prec_obs_fm,gddHigh_obs_fm,dd_obs_fm,
           gdd_1to3_fm,prec_1to3_fm,gddHigh_1to3_fm,dd_1to3_fm,
           gdd_4to10_fm,prec_4to10_fm,gddHigh_4to10_fm,dd_4to10_fm))

dat_str_cereals <- cbind(1,dat_str_all,df_farm$phi_cereals) 
dat_str_protein <- cbind(1,dat_str_all,df_farm$phi_protein)
dat_str_oilseed <- cbind(1,dat_str_all,df_farm$phi_oilseed)
dat_str_roots <- cbind(1,dat_str_all,df_farm$phi_roots)
dat_str_corn <- cbind(1,dat_str_all,df_farm$phi_corn)
dat_str_fert <- cbind(1,dat_str_all)

#------------------------------------------------#
##### 3b) Restore "non-Cholesky" parameters ######
#------------------------------------------------#
  
  coef_cereals <- c(b1=model_linear$coefficients["b1"] * mean(df_farm$qq_cereals) / mean(df_farm$PHI_cereals),
                    b11=model_linear$coefficients["c11"]*model_linear$coefficients["c11"] * mean(df_farm$qq_cereals) / mean(df_farm$np_cereals) / mean(df_farm$PHI_cereals),
                    b12=model_linear$coefficients["c11"]*model_linear$coefficients["c12"] * mean(df_farm$qq_cereals) / mean(df_farm$np_protein) / mean(df_farm$PHI_cereals),
                    b13=model_linear$coefficients["c11"]*model_linear$coefficients["c13"] * mean(df_farm$qq_cereals) / mean(df_farm$np_oilseed) / mean(df_farm$PHI_cereals),
                    b14=model_linear$coefficients["c11"]*model_linear$coefficients["c14"] * mean(df_farm$qq_cereals) / mean(df_farm$np_roots) / mean(df_farm$PHI_cereals),
                    b15=model_linear$coefficients["c11"]*model_linear$coefficients["c15"] * mean(df_farm$qq_cereals) / mean(df_farm$np_corn) / mean(df_farm$PHI_cereals),
                    b16=model_linear$coefficients["c11"]*model_linear$coefficients["c16"] * mean(df_farm$qq_cereals) / mean(df_farm$nw_fert) / mean(df_farm$PHI_cereals),
                    k11 <- model_linear$coefficients["k11"] * mean(df_farm$qq_cereals) / mean(df_farm$k_land) / mean(df_farm$PHI_cereals),
                    k12 <- model_linear$coefficients["k12"] * mean(df_farm$qq_cereals) / mean(df_farm$k_labor) / mean(df_farm$PHI_cereals),
                    k13 <- model_linear$coefficients["k13"] * mean(df_farm$qq_cereals) / mean(df_farm$k_capital) / mean(df_farm$PHI_cereals),
                    t1 <- model_linear$coefficients["t1"] * mean(df_farm$qq_cereals) / mean(df_farm$trend) / mean(df_farm$PHI_cereals),
                    tsq1 <- model_linear$coefficients["tsq1"] * mean(df_farm$qq_cereals) / mean(df_farm$trend2) / mean(df_farm$PHI_cereals),
                    w11 <- model_linear$coefficients["w11"] * mean(df_farm$qq_cereals) / mean(df_farm$gdd_obs) / mean(df_farm$PHI_cereals),
                    w12 <- model_linear$coefficients["w12"] * mean(df_farm$qq_cereals) / mean(df_farm$prec_obs) / mean(df_farm$PHI_cereals),
                    w13 <- model_linear$coefficients["w13"] * mean(df_farm$qq_cereals) / mean(df_farm$gddHigh_obs) / mean(df_farm$PHI_cereals),
                    w14 <- model_linear$coefficients["w14"] * mean(df_farm$qq_cereals) / mean(df_farm$dd_obs) / mean(df_farm$PHI_cereals),
                    w15 <- model_linear$coefficients["w15"] * mean(df_farm$qq_cereals) / mean(df_farm$gdd_1to3) / mean(df_farm$PHI_cereals),
                    w16 <- model_linear$coefficients["w16"] * mean(df_farm$qq_cereals) / mean(df_farm$prec_1to3) / mean(df_farm$PHI_cereals),
                    w17 <- model_linear$coefficients["w17"] * mean(df_farm$qq_cereals) / mean(df_farm$gddHigh_1to3) / mean(df_farm$PHI_cereals),
                    w18 <- model_linear$coefficients["w18"] * mean(df_farm$qq_cereals) / mean(df_farm$dd_1to3) / mean(df_farm$PHI_cereals),
                    w19 <- model_linear$coefficients["w19"] * mean(df_farm$qq_cereals) / mean(df_farm$gdd_4to10) / mean(df_farm$PHI_cereals),
                    w110 <- model_linear$coefficients["w110"] * mean(df_farm$qq_cereals) / mean(df_farm$prec_4to10) / mean(df_farm$PHI_cereals),
                    w111 <- model_linear$coefficients["w111"] * mean(df_farm$qq_cereals) / mean(df_farm$gddHigh_4to10) / mean(df_farm$PHI_cereals),
                    w112 <- model_linear$coefficients["w112"] * mean(df_farm$qq_cereals) / mean(df_farm$dd_4to10) / mean(df_farm$PHI_cereals),
                    m11 <- model_linear$coefficients["m11"] * mean(df_farm$qq_cereals) / mean(df_farm$np_cereals_fm) / mean(df_farm$PHI_cereals),
                    m12 <- model_linear$coefficients["m12"] * mean(df_farm$qq_cereals) / mean(df_farm$np_protein_fm) / mean(df_farm$PHI_cereals),
                    m13 <- model_linear$coefficients["m13"] * mean(df_farm$qq_cereals) / mean(df_farm$np_oilseed_fm) / mean(df_farm$PHI_cereals),
                    m14 <- model_linear$coefficients["m14"] * mean(df_farm$qq_cereals) / mean(df_farm$np_roots_fm) / mean(df_farm$PHI_cereals),
                    m15 <- model_linear$coefficients["m15"] * mean(df_farm$qq_cereals) / mean(df_farm$np_corn_fm) / mean(df_farm$PHI_cereals),
                    m16 <- model_linear$coefficients["m16"] * mean(df_farm$qq_cereals) / mean(df_farm$nw_fert_fm) / mean(df_farm$PHI_cereals),
                    m17 <- model_linear$coefficients["m17"] * mean(df_farm$qq_cereals) / mean(df_farm$k_land_fm) / mean(df_farm$PHI_cereals),
                    m18 <- model_linear$coefficients["m18"] * mean(df_farm$qq_cereals) / mean(df_farm$k_labor_fm) / mean(df_farm$PHI_cereals),
                    m19 <- model_linear$coefficients["m19"] * mean(df_farm$qq_cereals) / mean(df_farm$k_capital_fm) / mean(df_farm$PHI_cereals),
                    t1_fm <- model_linear$coefficients["t1_fm"] * mean(df_farm$qq_cereals) / mean(df_farm$trend_fm) / mean(df_farm$PHI_cereals),
                    tsq1_fm <- model_linear$coefficients["tsq1_fm"] * mean(df_farm$qq_cereals) / mean(df_farm$trend2_fm) / mean(df_farm$PHI_cereals),
                    w11_fm <- model_linear$coefficients["w11_fm"] * mean(df_farm$qq_cereals) / mean(df_farm$gdd_obs_fm) / mean(df_farm$PHI_cereals),
                    w12_fm <- model_linear$coefficients["w12_fm"] * mean(df_farm$qq_cereals) / mean(df_farm$prec_obs_fm) / mean(df_farm$PHI_cereals),
                    w13_fm <- model_linear$coefficients["w13_fm"] * mean(df_farm$qq_cereals) / mean(df_farm$gddHigh_obs_fm) / mean(df_farm$PHI_cereals),
                    w14_fm <- model_linear$coefficients["w14_fm"] * mean(df_farm$qq_cereals) / mean(df_farm$dd_obs_fm) / mean(df_farm$PHI_cereals),
                    w15_fm <- model_linear$coefficients["w15_fm"] * mean(df_farm$qq_cereals) / mean(df_farm$gdd_1to3_fm) / mean(df_farm$PHI_cereals),
                    w16_fm <- model_linear$coefficients["w16_fm"] * mean(df_farm$qq_cereals) / mean(df_farm$prec_1to3_fm) / mean(df_farm$PHI_cereals),
                    w17_fm <- model_linear$coefficients["w17_fm"] * mean(df_farm$qq_cereals) / mean(df_farm$gddHigh_1to3_fm) / mean(df_farm$PHI_cereals),
                    w18_fm <- model_linear$coefficients["w18_fm"] * mean(df_farm$qq_cereals) / mean(df_farm$dd_1to3_fm) / mean(df_farm$PHI_cereals),
                    w19_fm <- model_linear$coefficients["w19_fm"] * mean(df_farm$qq_cereals) / mean(df_farm$gdd_4to10_fm) / mean(df_farm$PHI_cereals),
                    w110_fm <- model_linear$coefficients["w110_fm"] * mean(df_farm$qq_cereals) / mean(df_farm$prec_4to10_fm) / mean(df_farm$PHI_cereals),
                    w111_fm <- model_linear$coefficients["w111_fm"] * mean(df_farm$qq_cereals) / mean(df_farm$gddHigh_4to10_fm) / mean(df_farm$PHI_cereals),
                    w112_fm <- model_linear$coefficients["w112_fm"] * mean(df_farm$qq_cereals) / mean(df_farm$dd_4to10_fm) / mean(df_farm$PHI_cereals),
                    phi1 <- model_linear$coefficients["phi11"] * mean(df_farm$qq_cereals) / mean(df_farm$phi_cereals))
  
  coef_protein <- c(
                    b2 <- model_linear$coefficients["b2"] * mean(df_farm$qq_protein) / mean(df_farm$PHI_protein),
                    b21 <- (model_linear$coefficients["c12"]*model_linear$coefficients["c11"]) * mean(df_farm$qq_protein) / mean(df_farm$np_cereals) / mean(df_farm$PHI_protein),
                    b22 <- (model_linear$coefficients["c12"]*model_linear$coefficients["c12"]+model_linear$coefficients["c22"]*model_linear$coefficients["c22"]) * mean(df_farm$qq_protein) / mean(df_farm$np_protein) / mean(df_farm$PHI_protein),
                    b23 <- (model_linear$coefficients["c12"]*model_linear$coefficients["c13"]+model_linear$coefficients["c22"]*model_linear$coefficients["c23"]) * mean(df_farm$qq_protein) / mean(df_farm$np_oilseed) / mean(df_farm$PHI_protein),
                    b24 <- (model_linear$coefficients["c12"]*model_linear$coefficients["c14"]+model_linear$coefficients["c22"]*model_linear$coefficients["c24"]) * mean(df_farm$qq_protein) / mean(df_farm$np_roots) / mean(df_farm$PHI_protein),
                    b25 <- (model_linear$coefficients["c12"]*model_linear$coefficients["c15"]+model_linear$coefficients["c22"]*model_linear$coefficients["c25"]) * mean(df_farm$qq_protein) / mean(df_farm$np_corn) / mean(df_farm$PHI_protein),
                    b26 <- (model_linear$coefficients["c12"]*model_linear$coefficients["c16"]+model_linear$coefficients["c22"]*model_linear$coefficients["c26"]) * mean(df_farm$qq_protein) / mean(df_farm$nw_fert) / mean(df_farm$PHI_protein),
                    k21 <- model_linear$coefficients["k21"] * mean(df_farm$qq_protein) / mean(df_farm$k_land) / mean(df_farm$PHI_protein),
                    k22 <- model_linear$coefficients["k22"] * mean(df_farm$qq_protein) / mean(df_farm$k_labor) / mean(df_farm$PHI_protein),
                    k23 <- model_linear$coefficients["k23"] * mean(df_farm$qq_protein) / mean(df_farm$k_capital) / mean(df_farm$PHI_protein),
                    t2 <- model_linear$coefficients["t2"] * mean(df_farm$qq_protein) / mean(df_farm$trend) / mean(df_farm$PHI_protein),
                    tsq2 <- model_linear$coefficients["tsq2"] * mean(df_farm$qq_protein) / mean(df_farm$trend2) / mean(df_farm$PHI_protein),
                    w21 <- model_linear$coefficients["w21"] * mean(df_farm$qq_protein) / mean(df_farm$gdd_obs) / mean(df_farm$PHI_protein),
                    w22 <- model_linear$coefficients["w22"] * mean(df_farm$qq_protein) / mean(df_farm$prec_obs) / mean(df_farm$PHI_protein),
                    w23 <- model_linear$coefficients["w23"] * mean(df_farm$qq_protein) / mean(df_farm$gddHigh_obs) / mean(df_farm$PHI_protein),
                    w24 <- model_linear$coefficients["w24"] * mean(df_farm$qq_protein) / mean(df_farm$dd_obs) / mean(df_farm$PHI_protein),
                    w25 <- model_linear$coefficients["w25"] * mean(df_farm$qq_protein) / mean(df_farm$gdd_1to3) / mean(df_farm$PHI_protein),
                    w26 <- model_linear$coefficients["w26"] * mean(df_farm$qq_protein) / mean(df_farm$prec_1to3) / mean(df_farm$PHI_protein),
                    w27 <- model_linear$coefficients["w27"] * mean(df_farm$qq_protein) / mean(df_farm$gddHigh_1to3) / mean(df_farm$PHI_protein),
                    w28 <- model_linear$coefficients["w28"] * mean(df_farm$qq_protein) / mean(df_farm$dd_1to3) / mean(df_farm$PHI_protein),
                    w29 <- model_linear$coefficients["w29"] * mean(df_farm$qq_protein) / mean(df_farm$gdd_4to10) / mean(df_farm$PHI_protein),
                    w210 <- model_linear$coefficients["w210"] * mean(df_farm$qq_protein) / mean(df_farm$prec_4to10) / mean(df_farm$PHI_protein),
                    w211 <- model_linear$coefficients["w211"] * mean(df_farm$qq_protein) / mean(df_farm$gddHigh_4to10) / mean(df_farm$PHI_protein),
                    w212 <- model_linear$coefficients["w212"] * mean(df_farm$qq_protein) / mean(df_farm$dd_4to10) / mean(df_farm$PHI_protein),
                    m21 <- model_linear$coefficients["m21"] * mean(df_farm$qq_protein) / mean(df_farm$np_cereals_fm) / mean(df_farm$PHI_protein),
                    m22 <- model_linear$coefficients["m22"] * mean(df_farm$qq_protein) / mean(df_farm$np_protein_fm) / mean(df_farm$PHI_protein),
                    m23 <- model_linear$coefficients["m23"] * mean(df_farm$qq_protein) / mean(df_farm$np_oilseed_fm) / mean(df_farm$PHI_protein),
                    m24 <- model_linear$coefficients["m24"] * mean(df_farm$qq_protein) / mean(df_farm$np_roots_fm) / mean(df_farm$PHI_protein),
                    m25 <- model_linear$coefficients["m25"] * mean(df_farm$qq_protein) / mean(df_farm$np_corn_fm) / mean(df_farm$PHI_protein),
                    m26 <- model_linear$coefficients["m26"] * mean(df_farm$qq_protein) / mean(df_farm$nw_fert_fm) / mean(df_farm$PHI_protein),
                    m27 <- model_linear$coefficients["m27"] * mean(df_farm$qq_protein) / mean(df_farm$k_land_fm) / mean(df_farm$PHI_protein),
                    m28 <- model_linear$coefficients["m28"] * mean(df_farm$qq_protein) / mean(df_farm$k_labor_fm) / mean(df_farm$PHI_protein),
                    m29 <- model_linear$coefficients["m29"] * mean(df_farm$qq_protein) / mean(df_farm$k_capital_fm) / mean(df_farm$PHI_protein),
                    t2_fm <- model_linear$coefficients["t2_fm"] * mean(df_farm$qq_protein) / mean(df_farm$trend_fm) / mean(df_farm$PHI_protein),
                    tsq2_fm <- model_linear$coefficients["tsq2_fm"] * mean(df_farm$qq_protein) / mean(df_farm$trend2_fm) / mean(df_farm$PHI_protein),
                    w21_fm <- model_linear$coefficients["w21_fm"] * mean(df_farm$qq_protein) / mean(df_farm$gdd_obs_fm) / mean(df_farm$PHI_protein),
                    w22_fm <- model_linear$coefficients["w22_fm"] * mean(df_farm$qq_protein) / mean(df_farm$prec_obs_fm) / mean(df_farm$PHI_protein),
                    w23_fm <- model_linear$coefficients["w23_fm"] * mean(df_farm$qq_protein) / mean(df_farm$gddHigh_obs_fm) / mean(df_farm$PHI_protein),
                    w24_fm <- model_linear$coefficients["w24_fm"] * mean(df_farm$qq_protein) / mean(df_farm$dd_obs_fm) / mean(df_farm$PHI_protein),
                    w25_fm <- model_linear$coefficients["w25_fm"] * mean(df_farm$qq_protein) / mean(df_farm$gdd_1to3_fm) / mean(df_farm$PHI_protein),
                    w26_fm <- model_linear$coefficients["w26_fm"] * mean(df_farm$qq_protein) / mean(df_farm$prec_1to3_fm) / mean(df_farm$PHI_protein),
                    w27_fm <- model_linear$coefficients["w27_fm"] * mean(df_farm$qq_protein) / mean(df_farm$gddHigh_1to3_fm) / mean(df_farm$PHI_protein),
                    w28_fm <- model_linear$coefficients["w28_fm"] * mean(df_farm$qq_protein) / mean(df_farm$dd_1to3_fm) / mean(df_farm$PHI_protein),
                    w29_fm <- model_linear$coefficients["w29_fm"] * mean(df_farm$qq_protein) / mean(df_farm$gdd_4to10_fm) / mean(df_farm$PHI_protein),
                    w210_fm <- model_linear$coefficients["w210_fm"] * mean(df_farm$qq_protein) / mean(df_farm$prec_4to10_fm) / mean(df_farm$PHI_protein),
                    w211_fm <- model_linear$coefficients["w211_fm"] * mean(df_farm$qq_protein) / mean(df_farm$gddHigh_4to10_fm) / mean(df_farm$PHI_protein),
                    w212_fm <- model_linear$coefficients["w212_fm"] * mean(df_farm$qq_protein) / mean(df_farm$dd_4to10_fm) / mean(df_farm$PHI_protein),
                    phi2 <- model_linear$coefficients["phi22"] * mean(df_farm$qq_protein) / mean(df_farm$phi_protein)
  )
  
  coef_oilseed <- c(
                    b3 <- model_linear$coefficients["b3"] * mean(df_farm$qq_oilseed) / mean(df_farm$PHI_oilseed),
                    b31 <- (model_linear$coefficients["c13"]*model_linear$coefficients["c11"]) * mean(df_farm$qq_oilseed) / mean(df_farm$np_cereals) / mean(df_farm$PHI_oilseed),
                    b32 <- (model_linear$coefficients["c13"]*model_linear$coefficients["c12"]+model_linear$coefficients["c23"]*model_linear$coefficients["c22"]) * mean(df_farm$qq_oilseed) / mean(df_farm$np_protein) / mean(df_farm$PHI_oilseed),
                    b33 <- (model_linear$coefficients["c13"]*model_linear$coefficients["c13"]+model_linear$coefficients["c23"]*model_linear$coefficients["c23"]+model_linear$coefficients["c33"]*model_linear$coefficients["c33"]) * mean(df_farm$qq_oilseed) / mean(df_farm$np_oilseed) / mean(df_farm$PHI_oilseed),
                    b34 <- (model_linear$coefficients["c13"]*model_linear$coefficients["c14"]+model_linear$coefficients["c23"]*model_linear$coefficients["c24"]+model_linear$coefficients["c33"]*model_linear$coefficients["c34"]) * mean(df_farm$qq_oilseed) / mean(df_farm$np_roots) / mean(df_farm$PHI_oilseed),
                    b35 <- (model_linear$coefficients["c13"]*model_linear$coefficients["c15"]+model_linear$coefficients["c23"]*model_linear$coefficients["c25"]+model_linear$coefficients["c33"]*model_linear$coefficients["c35"]) * mean(df_farm$qq_oilseed) / mean(df_farm$np_corn) / mean(df_farm$PHI_oilseed),
                    b36 <- (model_linear$coefficients["c13"]*model_linear$coefficients["c16"]+model_linear$coefficients["c23"]*model_linear$coefficients["c26"]+model_linear$coefficients["c33"]*model_linear$coefficients["c36"]) * mean(df_farm$qq_oilseed) / mean(df_farm$nw_fert) / mean(df_farm$PHI_oilseed),
                    k31 <- model_linear$coefficients["k31"] * mean(df_farm$qq_oilseed) / mean(df_farm$k_land) / mean(df_farm$PHI_oilseed),
                    k32 <- model_linear$coefficients["k32"] * mean(df_farm$qq_oilseed) / mean(df_farm$k_labor) / mean(df_farm$PHI_oilseed),
                    k33 <- model_linear$coefficients["k33"] * mean(df_farm$qq_oilseed) / mean(df_farm$k_capital) / mean(df_farm$PHI_oilseed),
                    t3 <- model_linear$coefficients["t3"] * mean(df_farm$qq_oilseed) / mean(df_farm$trend) / mean(df_farm$PHI_oilseed),
                    tsq3 <- model_linear$coefficients["tsq3"] * mean(df_farm$qq_oilseed) / mean(df_farm$trend2) / mean(df_farm$PHI_oilseed),
                    w31 <- model_linear$coefficients["w31"] * mean(df_farm$qq_oilseed) / mean(df_farm$gdd_obs) / mean(df_farm$PHI_oilseed),
                    w32 <- model_linear$coefficients["w32"] * mean(df_farm$qq_oilseed) / mean(df_farm$prec_obs) / mean(df_farm$PHI_oilseed),
                    w33 <- model_linear$coefficients["w33"] * mean(df_farm$qq_oilseed) / mean(df_farm$gddHigh_obs) / mean(df_farm$PHI_oilseed),
                    w34 <- model_linear$coefficients["w34"] * mean(df_farm$qq_oilseed) / mean(df_farm$dd_obs) / mean(df_farm$PHI_oilseed),
                    w35 <- model_linear$coefficients["w35"] * mean(df_farm$qq_oilseed) / mean(df_farm$gdd_1to3) / mean(df_farm$PHI_oilseed),
                    w36 <- model_linear$coefficients["w36"] * mean(df_farm$qq_oilseed) / mean(df_farm$prec_1to3) / mean(df_farm$PHI_oilseed),
                    w37 <- model_linear$coefficients["w37"] * mean(df_farm$qq_oilseed) / mean(df_farm$gddHigh_1to3) / mean(df_farm$PHI_oilseed),
                    w38 <- model_linear$coefficients["w38"] * mean(df_farm$qq_oilseed) / mean(df_farm$dd_1to3) / mean(df_farm$PHI_oilseed),
                    w39 <- model_linear$coefficients["w39"] * mean(df_farm$qq_oilseed) / mean(df_farm$gdd_4to10) / mean(df_farm$PHI_oilseed),
                    w310 <- model_linear$coefficients["w310"] * mean(df_farm$qq_oilseed) / mean(df_farm$prec_4to10) / mean(df_farm$PHI_oilseed),
                    w311 <- model_linear$coefficients["w311"] * mean(df_farm$qq_oilseed) / mean(df_farm$gddHigh_4to10) / mean(df_farm$PHI_oilseed),
                    w312 <- model_linear$coefficients["w312"] * mean(df_farm$qq_oilseed) / mean(df_farm$dd_4to10) / mean(df_farm$PHI_oilseed),
                    m31 <- model_linear$coefficients["m31"] * mean(df_farm$qq_oilseed) / mean(df_farm$np_cereals_fm) / mean(df_farm$PHI_oilseed),
                    m32 <- model_linear$coefficients["m32"] * mean(df_farm$qq_oilseed) / mean(df_farm$np_protein_fm) / mean(df_farm$PHI_oilseed),
                    m33 <- model_linear$coefficients["m33"] * mean(df_farm$qq_oilseed) / mean(df_farm$np_oilseed_fm) / mean(df_farm$PHI_oilseed),
                    m34 <- model_linear$coefficients["m34"] * mean(df_farm$qq_oilseed) / mean(df_farm$np_roots_fm) / mean(df_farm$PHI_oilseed),
                    m35 <- model_linear$coefficients["m35"] * mean(df_farm$qq_oilseed) / mean(df_farm$np_corn_fm) / mean(df_farm$PHI_oilseed),
                    m36 <- model_linear$coefficients["m36"] * mean(df_farm$qq_oilseed) / mean(df_farm$nw_fert_fm) / mean(df_farm$PHI_oilseed),
                    m37 <- model_linear$coefficients["m37"] * mean(df_farm$qq_oilseed) / mean(df_farm$k_land_fm) / mean(df_farm$PHI_oilseed),
                    m38 <- model_linear$coefficients["m38"] * mean(df_farm$qq_oilseed) / mean(df_farm$k_labor_fm) / mean(df_farm$PHI_oilseed),
                    m39 <- model_linear$coefficients["m39"] * mean(df_farm$qq_oilseed) / mean(df_farm$k_capital_fm) / mean(df_farm$PHI_oilseed),
                    t3_fm <- model_linear$coefficients["t3_fm"] * mean(df_farm$qq_oilseed) / mean(df_farm$trend_fm) / mean(df_farm$PHI_oilseed),
                    tsq3_fm <- model_linear$coefficients["tsq3_fm"] * mean(df_farm$qq_oilseed) / mean(df_farm$trend2_fm) / mean(df_farm$PHI_oilseed),
                    w31_fm <- model_linear$coefficients["w31_fm"] * mean(df_farm$qq_oilseed) / mean(df_farm$gdd_obs_fm) / mean(df_farm$PHI_oilseed),
                    w32_fm <- model_linear$coefficients["w32_fm"] * mean(df_farm$qq_oilseed) / mean(df_farm$prec_obs_fm) / mean(df_farm$PHI_oilseed),
                    w33_fm <- model_linear$coefficients["w33_fm"] * mean(df_farm$qq_oilseed) / mean(df_farm$gddHigh_obs_fm) / mean(df_farm$PHI_oilseed),
                    w34_fm <- model_linear$coefficients["w34_fm"] * mean(df_farm$qq_oilseed) / mean(df_farm$dd_obs_fm) / mean(df_farm$PHI_oilseed),
                    w35_fm <- model_linear$coefficients["w35_fm"] * mean(df_farm$qq_oilseed) / mean(df_farm$gdd_1to3_fm) / mean(df_farm$PHI_oilseed),
                    w36_fm <- model_linear$coefficients["w36_fm"] * mean(df_farm$qq_oilseed) / mean(df_farm$prec_1to3_fm) / mean(df_farm$PHI_oilseed),
                    w37_fm <- model_linear$coefficients["w37_fm"] * mean(df_farm$qq_oilseed) / mean(df_farm$gddHigh_1to3_fm) / mean(df_farm$PHI_oilseed),
                    w38_fm <- model_linear$coefficients["w38_fm"] * mean(df_farm$qq_oilseed) / mean(df_farm$dd_1to3_fm) / mean(df_farm$PHI_oilseed),
                    w39_fm <- model_linear$coefficients["w39_fm"] * mean(df_farm$qq_oilseed) / mean(df_farm$gdd_4to10_fm) / mean(df_farm$PHI_oilseed),
                    w310_fm <- model_linear$coefficients["w310_fm"] * mean(df_farm$qq_oilseed) / mean(df_farm$prec_4to10_fm) / mean(df_farm$PHI_oilseed),
                    w311_fm <- model_linear$coefficients["w311_fm"] * mean(df_farm$qq_oilseed) / mean(df_farm$gddHigh_4to10_fm) / mean(df_farm$PHI_oilseed),
                    w312_fm <- model_linear$coefficients["w312_fm"] * mean(df_farm$qq_oilseed) / mean(df_farm$dd_4to10_fm) / mean(df_farm$PHI_oilseed),
                    phi3 <- model_linear$coefficients["phi33"] * mean(df_farm$qq_oilseed) / mean(df_farm$phi_oilseed)
  )
  
  coef_roots <- c(
                    b4 <- model_linear$coefficients["b4"] * mean(df_farm$qq_roots) / mean(df_farm$PHI_roots),
                    b41 <- (model_linear$coefficients["c14"]*model_linear$coefficients["c11"]) * mean(df_farm$qq_roots) / mean(df_farm$np_cereals) / mean(df_farm$PHI_roots),
                    b42 <- (model_linear$coefficients["c14"]*model_linear$coefficients["c12"]+model_linear$coefficients["c24"]*model_linear$coefficients["c22"]) * mean(df_farm$qq_roots) / mean(df_farm$np_protein) / mean(df_farm$PHI_roots),
                    b43 <- (model_linear$coefficients["c14"]*model_linear$coefficients["c13"]+model_linear$coefficients["c24"]*model_linear$coefficients["c23"]+model_linear$coefficients["c34"]*model_linear$coefficients["c33"]) * mean(df_farm$qq_roots) / mean(df_farm$np_oilseed) / mean(df_farm$PHI_roots),
                    b44 <- (model_linear$coefficients["c14"]*model_linear$coefficients["c14"]+model_linear$coefficients["c24"]*model_linear$coefficients["c24"]+model_linear$coefficients["c34"]*model_linear$coefficients["c34"]+model_linear$coefficients["c44"]*model_linear$coefficients["c44"]) * mean(df_farm$qq_roots) / mean(df_farm$np_roots) / mean(df_farm$PHI_roots),
                    b45 <- (model_linear$coefficients["c14"]*model_linear$coefficients["c15"]+model_linear$coefficients["c24"]*model_linear$coefficients["c25"]+model_linear$coefficients["c34"]*model_linear$coefficients["c35"]+model_linear$coefficients["c44"]*model_linear$coefficients["c45"]) * mean(df_farm$qq_roots) / mean(df_farm$np_corn) / mean(df_farm$PHI_roots),
                    b46 <- (model_linear$coefficients["c14"]*model_linear$coefficients["c16"]+model_linear$coefficients["c24"]*model_linear$coefficients["c26"]+model_linear$coefficients["c34"]*model_linear$coefficients["c36"]+model_linear$coefficients["c44"]*model_linear$coefficients["c46"]) * mean(df_farm$qq_roots) / mean(df_farm$nw_fert) / mean(df_farm$PHI_roots),
                    k41 <- model_linear$coefficients["k41"] * mean(df_farm$qq_roots) / mean(df_farm$k_land) / mean(df_farm$PHI_roots),
                    k42 <- model_linear$coefficients["k42"] * mean(df_farm$qq_roots) / mean(df_farm$k_labor) / mean(df_farm$PHI_roots),
                    k43 <- model_linear$coefficients["k43"] * mean(df_farm$qq_roots) / mean(df_farm$k_capital) / mean(df_farm$PHI_roots),
                    t4 <- model_linear$coefficients["t4"] * mean(df_farm$qq_roots) / mean(df_farm$trend) / mean(df_farm$PHI_roots),
                    tsq4 <- model_linear$coefficients["tsq4"] * mean(df_farm$qq_roots) / mean(df_farm$trend2) / mean(df_farm$PHI_roots),
                    w41 <- model_linear$coefficients["w41"] * mean(df_farm$qq_roots) / mean(df_farm$gdd_obs) / mean(df_farm$PHI_roots),
                    w42 <- model_linear$coefficients["w42"] * mean(df_farm$qq_roots) / mean(df_farm$prec_obs) / mean(df_farm$PHI_roots),
                    w43 <- model_linear$coefficients["w43"] * mean(df_farm$qq_roots) / mean(df_farm$gddHigh_obs) / mean(df_farm$PHI_roots),
                    w44 <- model_linear$coefficients["w44"] * mean(df_farm$qq_roots) / mean(df_farm$dd_obs) / mean(df_farm$PHI_roots),
                    w45 <- model_linear$coefficients["w45"] * mean(df_farm$qq_roots) / mean(df_farm$gdd_1to3) / mean(df_farm$PHI_roots),
                    w46 <- model_linear$coefficients["w46"] * mean(df_farm$qq_roots) / mean(df_farm$prec_1to3) / mean(df_farm$PHI_roots),
                    w47 <- model_linear$coefficients["w47"] * mean(df_farm$qq_roots) / mean(df_farm$gddHigh_1to3) / mean(df_farm$PHI_roots),
                    w48 <- model_linear$coefficients["w48"] * mean(df_farm$qq_roots) / mean(df_farm$dd_1to3) / mean(df_farm$PHI_roots),
                    w49 <- model_linear$coefficients["w49"] * mean(df_farm$qq_roots) / mean(df_farm$gdd_4to10) / mean(df_farm$PHI_roots),
                    w410 <- model_linear$coefficients["w410"] * mean(df_farm$qq_roots) / mean(df_farm$prec_4to10) / mean(df_farm$PHI_roots),
                    w411 <- model_linear$coefficients["w411"] * mean(df_farm$qq_roots) / mean(df_farm$gddHigh_4to10) / mean(df_farm$PHI_roots),
                    w412 <- model_linear$coefficients["w412"] * mean(df_farm$qq_roots) / mean(df_farm$dd_4to10) / mean(df_farm$PHI_roots),
                    m41 <- model_linear$coefficients["m41"] * mean(df_farm$qq_roots) / mean(df_farm$np_cereals_fm) / mean(df_farm$PHI_roots),
                    m42 <- model_linear$coefficients["m42"] * mean(df_farm$qq_roots) / mean(df_farm$np_protein_fm) / mean(df_farm$PHI_roots),
                    m43 <- model_linear$coefficients["m43"] * mean(df_farm$qq_roots) / mean(df_farm$np_oilseed_fm) / mean(df_farm$PHI_roots),
                    m44 <- model_linear$coefficients["m44"] * mean(df_farm$qq_roots) / mean(df_farm$np_roots_fm) / mean(df_farm$PHI_roots),
                    m45 <- model_linear$coefficients["m45"] * mean(df_farm$qq_roots) / mean(df_farm$np_corn_fm) / mean(df_farm$PHI_roots),
                    m46 <- model_linear$coefficients["m46"] * mean(df_farm$qq_roots) / mean(df_farm$nw_fert_fm) / mean(df_farm$PHI_roots),
                    m47 <- model_linear$coefficients["m47"] * mean(df_farm$qq_roots) / mean(df_farm$k_land_fm) / mean(df_farm$PHI_roots),
                    m48 <- model_linear$coefficients["m48"] * mean(df_farm$qq_roots) / mean(df_farm$k_labor_fm) / mean(df_farm$PHI_roots),
                    m49 <- model_linear$coefficients["m49"] * mean(df_farm$qq_roots) / mean(df_farm$k_capital_fm) / mean(df_farm$PHI_roots),
                    t4_fm <- model_linear$coefficients["t4_fm"] * mean(df_farm$qq_roots) / mean(df_farm$trend_fm) / mean(df_farm$PHI_roots),
                    tsq4_fm <- model_linear$coefficients["tsq4_fm"] * mean(df_farm$qq_roots) / mean(df_farm$trend2_fm) / mean(df_farm$PHI_roots),
                    w41_fm <- model_linear$coefficients["w41_fm"] * mean(df_farm$qq_roots) / mean(df_farm$gdd_obs_fm) / mean(df_farm$PHI_roots),
                    w42_fm <- model_linear$coefficients["w42_fm"] * mean(df_farm$qq_roots) / mean(df_farm$prec_obs_fm) / mean(df_farm$PHI_roots),
                    w43_fm <- model_linear$coefficients["w43_fm"] * mean(df_farm$qq_roots) / mean(df_farm$gddHigh_obs_fm) / mean(df_farm$PHI_roots),
                    w44_fm <- model_linear$coefficients["w44_fm"] * mean(df_farm$qq_roots) / mean(df_farm$dd_obs_fm) / mean(df_farm$PHI_roots),
                    w45_fm <- model_linear$coefficients["w45_fm"] * mean(df_farm$qq_roots) / mean(df_farm$gdd_1to3_fm) / mean(df_farm$PHI_roots),
                    w46_fm <- model_linear$coefficients["w46_fm"] * mean(df_farm$qq_roots) / mean(df_farm$prec_1to3_fm) / mean(df_farm$PHI_roots),
                    w47_fm <- model_linear$coefficients["w47_fm"] * mean(df_farm$qq_roots) / mean(df_farm$gddHigh_1to3_fm) / mean(df_farm$PHI_roots),
                    w48_fm <- model_linear$coefficients["w48_fm"] * mean(df_farm$qq_roots) / mean(df_farm$dd_1to3_fm) / mean(df_farm$PHI_roots),
                    w49_fm <- model_linear$coefficients["w49_fm"] * mean(df_farm$qq_roots) / mean(df_farm$gdd_4to10_fm) / mean(df_farm$PHI_roots),
                    w410_fm <- model_linear$coefficients["w410_fm"] * mean(df_farm$qq_roots) / mean(df_farm$prec_4to10_fm) / mean(df_farm$PHI_roots),
                    w411_fm <- model_linear$coefficients["w411_fm"] * mean(df_farm$qq_roots) / mean(df_farm$gddHigh_4to10_fm) / mean(df_farm$PHI_roots),
                    w412_fm <- model_linear$coefficients["w412_fm"] * mean(df_farm$qq_roots) / mean(df_farm$dd_4to10_fm) / mean(df_farm$PHI_roots),
                    phi4 <- model_linear$coefficients["phi44"] * mean(df_farm$qq_roots) / mean(df_farm$phi_roots)
  )
  
  coef_corn <- c(
                    b5 <- model_linear$coefficients["b5"] * mean(df_farm$qq_corn) / mean(df_farm$PHI_corn),
                    b51 <- (model_linear$coefficients["c15"]*model_linear$coefficients["c11"]) * mean(df_farm$qq_corn) / mean(df_farm$np_cereals) / mean(df_farm$PHI_corn),
                    b52 <- (model_linear$coefficients["c15"]*model_linear$coefficients["c12"]+model_linear$coefficients["c25"]*model_linear$coefficients["c22"]) * mean(df_farm$qq_corn) / mean(df_farm$np_protein) / mean(df_farm$PHI_corn),
                    b53 <- (model_linear$coefficients["c15"]*model_linear$coefficients["c13"]+model_linear$coefficients["c25"]*model_linear$coefficients["c23"]+model_linear$coefficients["c35"]*model_linear$coefficients["c33"]) * mean(df_farm$qq_corn) / mean(df_farm$np_oilseed) / mean(df_farm$PHI_corn),
                    b54 <- (model_linear$coefficients["c15"]*model_linear$coefficients["c14"]+model_linear$coefficients["c25"]*model_linear$coefficients["c24"]+model_linear$coefficients["c35"]*model_linear$coefficients["c34"]+model_linear$coefficients["c45"]*model_linear$coefficients["c44"]) * mean(df_farm$qq_corn) / mean(df_farm$np_roots) / mean(df_farm$PHI_corn),
                    b55 <- (model_linear$coefficients["c15"]*model_linear$coefficients["c15"]+model_linear$coefficients["c25"]*model_linear$coefficients["c25"]+model_linear$coefficients["c35"]*model_linear$coefficients["c35"]+model_linear$coefficients["c45"]*model_linear$coefficients["c45"]) * mean(df_farm$qq_corn) / mean(df_farm$np_corn) / mean(df_farm$PHI_corn),
                    b56 <- (model_linear$coefficients["c15"]*model_linear$coefficients["c16"]+model_linear$coefficients["c25"]*model_linear$coefficients["c26"]+model_linear$coefficients["c35"]*model_linear$coefficients["c36"]+model_linear$coefficients["c45"]*model_linear$coefficients["c46"]) * mean(df_farm$qq_corn) / mean(df_farm$nw_fert) / mean(df_farm$PHI_corn),
                    k51 <- model_linear$coefficients["k51"] * mean(df_farm$qq_corn) / mean(df_farm$k_land) / mean(df_farm$PHI_corn),
                    k52 <- model_linear$coefficients["k52"] * mean(df_farm$qq_corn) / mean(df_farm$k_labor) / mean(df_farm$PHI_corn),
                    k53 <- model_linear$coefficients["k53"] * mean(df_farm$qq_corn) / mean(df_farm$k_capital) / mean(df_farm$PHI_corn),
                    t5 <- model_linear$coefficients["t5"] * mean(df_farm$qq_corn) / mean(df_farm$trend) / mean(df_farm$PHI_corn),
                    tsq5 <- model_linear$coefficients["tsq5"] * mean(df_farm$qq_corn) / mean(df_farm$trend2) / mean(df_farm$PHI_corn),
                    w51 <- model_linear$coefficients["w51"] * mean(df_farm$qq_corn) / mean(df_farm$gdd_obs) / mean(df_farm$PHI_corn),
                    w52 <- model_linear$coefficients["w52"] * mean(df_farm$qq_corn) / mean(df_farm$prec_obs) / mean(df_farm$PHI_corn),
                    w53 <- model_linear$coefficients["w53"] * mean(df_farm$qq_corn) / mean(df_farm$gddHigh_obs) / mean(df_farm$PHI_corn),
                    w54 <- model_linear$coefficients["w54"] * mean(df_farm$qq_corn) / mean(df_farm$dd_obs) / mean(df_farm$PHI_corn),
                    w55 <- model_linear$coefficients["w55"] * mean(df_farm$qq_corn) / mean(df_farm$gdd_1to3) / mean(df_farm$PHI_corn),
                    w56 <- model_linear$coefficients["w56"] * mean(df_farm$qq_corn) / mean(df_farm$prec_1to3) / mean(df_farm$PHI_corn),
                    w57 <- model_linear$coefficients["w57"] * mean(df_farm$qq_corn) / mean(df_farm$gddHigh_1to3) / mean(df_farm$PHI_corn),
                    w58 <- model_linear$coefficients["w58"] * mean(df_farm$qq_corn) / mean(df_farm$dd_1to3) / mean(df_farm$PHI_corn),
                    w59 <- model_linear$coefficients["w59"] * mean(df_farm$qq_corn) / mean(df_farm$gdd_4to10) / mean(df_farm$PHI_corn),
                    w510 <- model_linear$coefficients["w510"] * mean(df_farm$qq_corn) / mean(df_farm$prec_4to10) / mean(df_farm$PHI_corn),
                    w511 <- model_linear$coefficients["w511"] * mean(df_farm$qq_corn) / mean(df_farm$gddHigh_4to10) / mean(df_farm$PHI_corn),
                    w512 <- model_linear$coefficients["w512"] * mean(df_farm$qq_corn) / mean(df_farm$dd_4to10) / mean(df_farm$PHI_corn),
                    m51 <- model_linear$coefficients["m51"] * mean(df_farm$qq_corn) / mean(df_farm$np_cereals_fm) / mean(df_farm$PHI_corn),
                    m52 <- model_linear$coefficients["m52"] * mean(df_farm$qq_corn) / mean(df_farm$np_protein_fm) / mean(df_farm$PHI_corn),
                    m53 <- model_linear$coefficients["m53"] * mean(df_farm$qq_corn) / mean(df_farm$np_oilseed_fm) / mean(df_farm$PHI_corn),
                    m54 <- model_linear$coefficients["m54"] * mean(df_farm$qq_corn) / mean(df_farm$np_roots_fm) / mean(df_farm$PHI_corn),
                    m55 <- model_linear$coefficients["m55"] * mean(df_farm$qq_corn) / mean(df_farm$np_corn_fm) / mean(df_farm$PHI_corn),
                    m56 <- model_linear$coefficients["m56"] * mean(df_farm$qq_corn) / mean(df_farm$nw_fert_fm) / mean(df_farm$PHI_corn),
                    m57 <- model_linear$coefficients["m57"] * mean(df_farm$qq_corn) / mean(df_farm$k_land_fm) / mean(df_farm$PHI_corn),
                    m58 <- model_linear$coefficients["m58"] * mean(df_farm$qq_corn) / mean(df_farm$k_labor_fm) / mean(df_farm$PHI_corn),
                    m59 <- model_linear$coefficients["m59"] * mean(df_farm$qq_corn) / mean(df_farm$k_capital_fm) / mean(df_farm$PHI_corn),
                    t5_fm <- model_linear$coefficients["t5_fm"] * mean(df_farm$qq_corn) / mean(df_farm$trend_fm) / mean(df_farm$PHI_corn),
                    tsq5_fm <- model_linear$coefficients["tsq5_fm"] * mean(df_farm$qq_corn) / mean(df_farm$trend2_fm) / mean(df_farm$PHI_corn),
                    w51_fm <- model_linear$coefficients["w51_fm"] * mean(df_farm$qq_corn) / mean(df_farm$gdd_obs_fm) / mean(df_farm$PHI_corn),
                    w52_fm <- model_linear$coefficients["w52_fm"] * mean(df_farm$qq_corn) / mean(df_farm$prec_obs_fm) / mean(df_farm$PHI_corn),
                    w53_fm <- model_linear$coefficients["w53_fm"] * mean(df_farm$qq_corn) / mean(df_farm$gddHigh_obs_fm) / mean(df_farm$PHI_corn),
                    w54_fm <- model_linear$coefficients["w54_fm"] * mean(df_farm$qq_corn) / mean(df_farm$dd_obs_fm) / mean(df_farm$PHI_corn),
                    w55_fm <- model_linear$coefficients["w55_fm"] * mean(df_farm$qq_corn) / mean(df_farm$gdd_1to3_fm) / mean(df_farm$PHI_corn),
                    w56_fm <- model_linear$coefficients["w56_fm"] * mean(df_farm$qq_corn) / mean(df_farm$prec_1to3_fm) / mean(df_farm$PHI_corn),
                    w57_fm <- model_linear$coefficients["w57_fm"] * mean(df_farm$qq_corn) / mean(df_farm$gddHigh_1to3_fm) / mean(df_farm$PHI_corn),
                    w58_fm <- model_linear$coefficients["w58_fm"] * mean(df_farm$qq_corn) / mean(df_farm$dd_1to3_fm) / mean(df_farm$PHI_corn),
                    w59_fm <- model_linear$coefficients["w59_fm"] * mean(df_farm$qq_corn) / mean(df_farm$gdd_4to10_fm) / mean(df_farm$PHI_corn),
                    w510_fm <- model_linear$coefficients["w510_fm"] * mean(df_farm$qq_corn) / mean(df_farm$prec_4to10_fm) / mean(df_farm$PHI_corn),
                    w511_fm <- model_linear$coefficients["w511_fm"] * mean(df_farm$qq_corn) / mean(df_farm$gddHigh_4to10_fm) / mean(df_farm$PHI_corn),
                    w512_fm <- model_linear$coefficients["w512_fm"] * mean(df_farm$qq_corn) / mean(df_farm$dd_4to10_fm) / mean(df_farm$PHI_corn),
                    phi5 <- model_linear$coefficients["phi55"] * mean(df_farm$qq_corn) / mean(df_farm$phi_corn)
  )
  
  coef_fert <- c(
                    b6 <- model_linear$coefficients["b6"] * mean(df_farm$x_fert),
                    b61 <- (model_linear$coefficients["c16"]*model_linear$coefficients["c11"]) * mean(df_farm$x_fert) / mean(df_farm$np_cereals),
                    b62 <- (model_linear$coefficients["c16"]*model_linear$coefficients["c12"]+model_linear$coefficients["c26"]*model_linear$coefficients["c22"]) * mean(df_farm$x_fert) / mean(df_farm$np_protein),
                    b63 <- (model_linear$coefficients["c16"]*model_linear$coefficients["c13"]+model_linear$coefficients["c26"]*model_linear$coefficients["c23"]+model_linear$coefficients["c36"]*model_linear$coefficients["c33"]) * mean(df_farm$x_fert) / mean(df_farm$np_oilseed),
                    b64 <- (model_linear$coefficients["c16"]*model_linear$coefficients["c14"]+model_linear$coefficients["c26"]*model_linear$coefficients["c24"]+model_linear$coefficients["c36"]*model_linear$coefficients["c34"]+model_linear$coefficients["c46"]*model_linear$coefficients["c44"]) * mean(df_farm$x_fert) / mean(df_farm$np_roots),
                    b65 <- (model_linear$coefficients["c16"]*model_linear$coefficients["c15"]+model_linear$coefficients["c26"]*model_linear$coefficients["c25"]+model_linear$coefficients["c36"]*model_linear$coefficients["c35"]+model_linear$coefficients["c46"]*model_linear$coefficients["c45"]) * mean(df_farm$x_fert) / mean(df_farm$np_corn),
                    b66 <- (model_linear$coefficients["c16"]*model_linear$coefficients["c16"]+model_linear$coefficients["c26"]*model_linear$coefficients["c26"]+model_linear$coefficients["c36"]*model_linear$coefficients["c36"]+model_linear$coefficients["c46"]*model_linear$coefficients["c46"]) * mean(df_farm$x_fert) / mean(df_farm$nw_fert),
                    k61 <- model_linear$coefficients["k61"] * mean(df_farm$x_fert) / mean(df_farm$k_land),
                    k62 <- model_linear$coefficients["k62"] * mean(df_farm$x_fert) / mean(df_farm$k_labor),
                    k63 <- model_linear$coefficients["k63"] * mean(df_farm$x_fert) / mean(df_farm$k_capital),
                    t6 <- model_linear$coefficients["t6"] * mean(df_farm$x_fert) / mean(df_farm$trend) ,
                    tsq6 <- model_linear$coefficients["tsq6"] * mean(df_farm$x_fert) / mean(df_farm$trend2) ,
                    w61 <- model_linear$coefficients["w61"] * mean(df_farm$x_fert) / mean(df_farm$gdd_obs) ,
                    w62 <- model_linear$coefficients["w62"] * mean(df_farm$x_fert) / mean(df_farm$prec_obs) ,
                    w63 <- model_linear$coefficients["w63"] * mean(df_farm$x_fert) / mean(df_farm$gddHigh_obs) ,
                    w64 <- model_linear$coefficients["w64"] * mean(df_farm$x_fert) / mean(df_farm$dd_obs) ,
                    w65 <- model_linear$coefficients["w65"] * mean(df_farm$x_fert) / mean(df_farm$gdd_1to3) ,
                    w66 <- model_linear$coefficients["w66"] * mean(df_farm$x_fert) / mean(df_farm$prec_1to3) ,
                    w67 <- model_linear$coefficients["w67"] * mean(df_farm$x_fert) / mean(df_farm$gddHigh_1to3) ,
                    w68 <- model_linear$coefficients["w68"] * mean(df_farm$x_fert) / mean(df_farm$dd_1to3) ,
                    w69 <- model_linear$coefficients["w69"] * mean(df_farm$x_fert) / mean(df_farm$gdd_4to10) ,
                    w610 <- model_linear$coefficients["w610"] * mean(df_farm$x_fert) / mean(df_farm$prec_4to10) ,
                    w611 <- model_linear$coefficients["w611"] * mean(df_farm$x_fert) / mean(df_farm$gddHigh_4to10) ,
                    w612 <- model_linear$coefficients["w612"] * mean(df_farm$x_fert) / mean(df_farm$dd_4to10),
                    m61 <- model_linear$coefficients["m61"] * mean(df_farm$x_fert) / mean(df_farm$np_cereals_fm),
                    m62 <- model_linear$coefficients["m62"] * mean(df_farm$x_fert) / mean(df_farm$np_protein_fm),
                    m63 <- model_linear$coefficients["m63"] * mean(df_farm$x_fert) / mean(df_farm$np_oilseed_fm),
                    m64 <- model_linear$coefficients["m64"] * mean(df_farm$x_fert) / mean(df_farm$np_roots_fm),
                    m65 <- model_linear$coefficients["m65"] * mean(df_farm$x_fert) / mean(df_farm$np_corn_fm),
                    m66 <- model_linear$coefficients["m66"] * mean(df_farm$x_fert) / mean(df_farm$nw_fert_fm),
                    m67 <- model_linear$coefficients["m67"] * mean(df_farm$x_fert) / mean(df_farm$k_land_fm),
                    m68 <- model_linear$coefficients["m68"] * mean(df_farm$x_fert) / mean(df_farm$k_labor_fm),
                    m69 <- model_linear$coefficients["m69"] * mean(df_farm$x_fert) / mean(df_farm$k_capital_fm),
                    t6_fm <- model_linear$coefficients["t6_fm"] * mean(df_farm$x_fert) / mean(df_farm$trend_fm) ,
                    tsq6_fm <- model_linear$coefficients["tsq6_fm"] * mean(df_farm$x_fert) / mean(df_farm$trend2_fm), 
                    w61_fm <- model_linear$coefficients["w61_fm"] * mean(df_farm$x_fert) / mean(df_farm$gdd_obs_fm),
                    w62_fm <- model_linear$coefficients["w62_fm"] * mean(df_farm$x_fert) / mean(df_farm$prec_obs_fm),
                    w63_fm <- model_linear$coefficients["w63_fm"] * mean(df_farm$x_fert) / mean(df_farm$gddHigh_obs_fm),
                    w64_fm <- model_linear$coefficients["w64_fm"] * mean(df_farm$x_fert) / mean(df_farm$dd_obs_fm),
                    w65_fm <- model_linear$coefficients["w65_fm"] * mean(df_farm$x_fert) / mean(df_farm$gdd_1to3_fm),
                    w66_fm <- model_linear$coefficients["w66_fm"] * mean(df_farm$x_fert) / mean(df_farm$prec_1to3_fm),
                    w67_fm <- model_linear$coefficients["w67_fm"] * mean(df_farm$x_fert) / mean(df_farm$gddHigh_1to3_fm),
                    w68_fm <- model_linear$coefficients["w68_fm"] * mean(df_farm$x_fert) / mean(df_farm$dd_1to3_fm),
                    w69_fm <- model_linear$coefficients["w69_fm"] * mean(df_farm$x_fert) / mean(df_farm$gdd_4to10_fm),
                    w610_fm <- model_linear$coefficients["w610_fm"] * mean(df_farm$x_fert) / mean(df_farm$prec_4to10_fm),
                    w611_fm <- model_linear$coefficients["w611_fm"] * mean(df_farm$x_fert) / mean(df_farm$gddHigh_4to10_fm),
                    w612_fm <- model_linear$coefficients["w612_fm"] * mean(df_farm$x_fert) / mean(df_farm$dd_4to10_fm))

  
# ------------------------------------ #
#### 4) Simulation of drought shock ####
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
                  "np_cereals_fm","np_protein_fm","np_oilseed_fm","np_roots_fm","np_corn_fm","nw_fert_fm",
                  "k_land_fm","k_labor_fm","k_capital_fm","trend_fm","trend2_fm",
                  "gdd_obs_fm","prec_obs_fm","gddHigh_obs_fm","dd_obs_fm",
                  "gdd_1to3_fm","prec_1to3_fm","gddHigh_1to3_fm","dd_1to3_fm",
                  "gdd_4to10_fm","prec_4to10_fm","gddHigh_4to10_fm","dd_4to10_fm",
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
      C <- as.data.frame(t(get(paste0("coef_",list_outcome[i,"outcome"]))) ) ## note: need to check if this works xxx
      C <- t(C) 
      
      list_simresults[[i]][[t+2]] <- PHI * (X_sim %*% C) 
            
    }
  }
    
  # Calculate sample mean of simulations
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
    Fertilizer = c(mean(list_simresults$Fertilizer$t1), mean(list_simresults$Fertilizer$t2), mean(list_simresults$Fertilizer$t3), 
                   mean(list_simresults$Fertilizer$t4), mean(list_simresults$Fertilizer$t5), mean(list_simresults$Fertilizer$t6),
                   mean(list_simresults$Fertilizer$t7), mean(list_simresults$Fertilizer$t8), mean(list_simresults$Fertilizer$t9),
                   mean(list_simresults$Fertilizer$t10), mean(list_simresults$Fertilizer$t11), mean(list_simresults$Fertilizer$t12))
  )
  
  # Save simulation result
  save(simresults, file="rOutput/simresults_restr.Rda")
    
    
#----------------------------------------------#
#### 5) Plot the results at the sample mean ####
#----------------------------------------------#

  # Load results
  load("rOutput/simresults_restr.Rda")
  
  # Normalize the data for plotting
  simresults <- simresults %>% 
    mutate_each(funs(./.[1]), setdiff(names(.), c("t")))
  
  # Transform data to long data
  simresults_long <- reshape2::melt(simresults, id.vars="t")
    
  # Plot the results
  Fig_restr <- ggplot(simresults_long, aes(x=t, y=value, group=variable, shape=variable, color=variable)) +
    geom_line()+
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
    scale_y_continuous(breaks = seq(0, 1.4, by = 0.2), limits=c(0,1.52)) +  
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
    geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1.45), linetype="dashed", color="grey") +
    annotate(geom="text", x=0, y=1.52, label="Drought shock",
             color="darkgrey", size=3.5)
  
  Fig_restr
  
  # Save the plot as .png and .eps files
  ggsave("Figures/Figure_S6.5.1.png", Fig_restr, device="png", width = 6, height = 3, units = "in", dpi=1200)
  ggsave("Figures/Figure_S6.5.1.eps", Fig_restr, device="eps", width = 6, height = 3, units = "in", dpi=1200)
  
  