#----------------------------------------------------------#
#                                                          #
# This program estimates first- and second-stage           #
# regressions and simulates the outcomes for the 2018-     #
# drought, for the model without non-random crop selection #
#                                                          #
# Note: The farm-level data are confidential and cannot    # 
# be loaded (line 17). The results at the sample mean      #
# are saved in line 323 and the figure in Section 5        # 
# can be created.                                          #
#                                                          #
#----------------------------------------------------------#

library(dplyr)
library(systemfit)
library(ggplot2)
load("rOutput/farm_ready.Rda")

#-----------------------------------------------------------#
#### 1) Prepare the weather variables for the simulation ####
#-----------------------------------------------------------#

# Load simulated weather
  
  load("rOutput/sim2018.Rda")
  sim2018 <- sim2018

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

  # --> not needed here because selection is ignored
  
#--------------------------------------#
#### 3) Run structural regressions #####
#--------------------------------------#

  #Define structural equations
  eqQQcereals <- qq_cereals ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
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
  
  eqQQprotein <- qq_protein ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
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
  
  eqQQoilseed <- qq_oilseed ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
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
  
  eqQQroots <- qq_roots ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
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
  
  eqQQcorn <- qq_corn ~ np_cereals + np_protein + np_oilseed + np_roots + np_corn + nw_fert + 
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
restrict <- c( "QQcereals_np_protein - QQprotein_np_cereals = 0",
               "QQcereals_np_oilseed - QQoilseed_np_cereals = 0",
               "QQcereals_np_roots   - QQroots_np_cereals = 0",
               "QQcereals_np_corn - QQcorn_np_cereals = 0",
               "QQcereals_nw_fert    - NXfert_np_cereals = 0",
               "QQprotein_np_oilseed - QQoilseed_np_protein = 0",
               "QQprotein_np_roots   - QQroots_np_protein = 0",
               "QQprotein_np_corn - QQcorn_np_protein = 0",
               "QQprotein_nw_fert    - NXfert_np_protein = 0",
               "QQoilseed_np_roots   - QQroots_np_oilseed = 0",
               "QQoilseed_np_corn - QQcorn_np_oilseed = 0",
               "QQoilseed_nw_fert    - NXfert_np_oilseed = 0",
               "QQroots_np_corn   - QQcorn_np_roots = 0",
               "QQroots_nw_fert      - NXfert_np_roots = 0",
               "QQcorn_nw_fert    - NXfert_np_corn = 0")

## Regression with iterated SUR estimation
model_linear <- systemfit( formula = system, method = "SUR", 
                           data = df_farm, restrict.matrix = restrict, 
                           maxit = 100 )


# Extract data used in this estimation for simulation later

    dat_str_cereals <- as_tibble(model.matrix(model_linear)[grep('^QQcereals', rownames(model.matrix(model_linear))),]) %>% 
      select(where(~ any(. != 0)))
    
    dat_str_protein <- as_tibble(model.matrix(model_linear)[grep('^QQprotein', rownames(model.matrix(model_linear))),]) %>% 
      select(where(~ any(. != 0)))
    
    dat_str_oilseed <- as_tibble(model.matrix(model_linear)[grep('^QQoilseed', rownames(model.matrix(model_linear))),]) %>% 
      select(where(~ any(. != 0)))
    
    dat_str_roots <- as_tibble(model.matrix(model_linear)[grep('^QQroots', rownames(model.matrix(model_linear))),]) %>% 
      select(where(~ any(. != 0)))
    
    dat_str_corn <- as_tibble(model.matrix(model_linear)[grep('^QQcorn', rownames(model.matrix(model_linear))),]) %>% 
      select(where(~ any(. != 0)))
    
    dat_str_fert <- as_tibble(model.matrix(model_linear)[grep('^NXfert', rownames(model.matrix(model_linear))),]) %>% 
      select(where(~ any(. != 0)))
    

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
                  "gdd_obs_1to3","prec_obs_1to3","gddHigh_obs_1to3","dd_obs_1to3",
                  "gdd_obs_4to10","prec_obs_4to10","gddHigh_obs_4to10","dd_obs_4to10",
                  "np_cereals_fm","np_protein_fm","np_oilseed_fm","np_roots_fm","np_corn_fm","nw_fert_fm",
                  "k_land_fm","k_labor_fm","k_capital_fm","trend_fm","trend2_fm",
                  "gdd_obs_fm","prec_obs_fm","gddHigh_obs_fm","dd_obs_fm",
                  "gdd_1to3_fm","prec_1to3_fm","gddHigh_1to3_fm","dd_1to3_fm",
                  "gdd_4to10_fm","prec_4to10_fm","gddHigh_4to10_fm","dd_4to10_fm",
                  "gdd_obs_1to3_fm","prec_obs_1to3_fm","gddHigh_obs_1to3_fm","dd_obs_1to3_fm",
                  "gdd_obs_4to10_fm","prec_obs_4to10_fm","gddHigh_obs_4to10_fm","dd_obs_4to10_fm")
    ## note: need to check if X includes the constant term xxx
      
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
            
      X_sim$gdd_obs_1to3 <- sim2018$gdd_obs_1to3[sim2018$t==t]
      X_sim$prec_obs_1to3 <- sim2018$prec_obs_1to3[sim2018$t==t]
      X_sim$gddHigh_obs_1to3 <- sim2018$gddHigh_obs_1to3[sim2018$t==t]
      X_sim$dd_obs_1to3 <- sim2018$dd_obs_1to3[sim2018$t==t]
      X_sim$gdd_obs_4to10 <- sim2018$gdd_obs_4to10[sim2018$t==t]
      X_sim$prec_obs_4to10 <- sim2018$prec_obs_4to10[sim2018$t==t]
      X_sim$gddHigh_obs_4to10 <- sim2018$gddHigh_obs_4to10[sim2018$t==t]
      X_sim$dd_obs_4to10 <- sim2018$dd_obs_4to10[sim2018$t==t]
        
      # Replace phi and PHI for all crops, but not for fert
      
          # --> not needed here because we ignore selection
          
      # Predict output 
      X_sim <- as.matrix(X_sim)
      C <- as.data.frame(t(coef_linear) ) %>% select(starts_with(list_outcome[i,"sur"])) ## note: check if this works xxx
      C <- t(C) 
      
      list_simresults[[i]][[t+2]] <- (X_sim %*% C) 
            
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
    Fertilizer = -c(mean(list_simresults$Fertilizer$t1), mean(list_simresults$Fertilizer$t2), mean(list_simresults$Fertilizer$t3),
                    mean(list_simresults$Fertilizer$t4), mean(list_simresults$Fertilizer$t5), mean(list_simresults$Fertilizer$t6),
                    mean(list_simresults$Fertilizer$t7), mean(list_simresults$Fertilizer$t8), mean(list_simresults$Fertilizer$t9),
                    mean(list_simresults$Fertilizer$t10), mean(list_simresults$Fertilizer$t11), mean(list_simresults$Fertilizer$t12))
  )
    
  # Save simulation result
  save(simresults, file="rOutput/simresults_wosel.Rda")
    
    
#----------------------------------------------#
#### 5) Plot the results at the sample mean ####
#----------------------------------------------#

  # Load results
  load("rOutput/simresults_wosel.Rda")
  
  # Normalize the data for plotting
  simresults <- simresults %>% 
    mutate_each(funs(./.[1]), setdiff(names(.), c("t")))
  
  # Transform data to long data
  simresults_long <- reshape2::melt(simresults, id.vars="t")
    
  # Plot the results
  Fig_wosel <- ggplot(simresults_long, aes(x=t, y=value, group=variable, shape=variable, color=variable)) +
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
  
  Fig_wosel
  
  # Save the plot as .png and .eps files
  ggsave("Figures/Figure_S8.1.png", Fig_wosel, device="png", width = 6, height = 3, units = "in", dpi=1200)
  ggsave("Figures/Figure_S8.1.eps", Fig_wosel, device="eps", width = 6, height = 3, units = "in", dpi=1200)
  