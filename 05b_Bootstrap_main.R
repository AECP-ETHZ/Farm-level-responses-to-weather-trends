#-----------------------------------------------------#
#                                                     #
# This program estimates the main model with boot-    #
# strapped standard errors.                           #
#                                                     #
# Note: The bootstrap program (Section 1: Bootstrap   #
# program) cannot be run because farm-level data are  #
# confidential. The results of the bootstrap  are     #
# stored in lines 916, 931, and 946, allowing to      # 
# create the tables in Section 2 (Presentation of     #
# results).                                           # 
#                                                     #
#-----------------------------------------------------#

# load packages
library(dplyr)
library(writexl)
library(ggplot2)
library(ggpattern)
library(systemfit)
library(doFuture) 
library(doRNG) 

# load data
load("rOutput/farm_ready.Rda")

# keep only variables needed
df_farm <- df_farm %>% select(key, year, nuts2, 
                              iCereals, iProtein, iOilseed, iRoots, iCorn, 
                              np_cereals, np_protein, np_oilseed, np_roots, np_corn, nw_fert, 
                              k_land, k_labor, k_capital, trend, trend2, 
                              lsh_cereals, lsh_protein, lsh_oilseed, lsh_roots, 
                              gdd_1to3, prec_1to3, gddHigh_1to3, dd_1to3, 
                              gdd_4to10, prec_4to10, gddHigh_4to10, dd_4to10,
                              np_cereals_fm, np_protein_fm, np_oilseed_fm, np_roots_fm, np_corn_fm, nw_fert_fm, 
                              k_land_fm, k_labor_fm, k_capital_fm, trend_fm, trend2_fm, 
                              lsh_cereals_fm, lsh_protein_fm, lsh_oilseed_fm, lsh_roots_fm, 
                              gdd_1to3_fm, prec_1to3_fm, gddHigh_1to3_fm, dd_1to3_fm, 
                              gdd_4to10_fm, prec_4to10_fm, gddHigh_4to10_fm, dd_4to10_fm,
                              qq_cereals, qq_protein, qq_oilseed, qq_roots, qq_corn, nx_fert, x_fert, x_otherinp,
                              gdd_obs, prec_obs, gddHigh_obs, dd_obs,
                              gdd_obs_1to3, prec_obs_1to3, gddHigh_obs_1to3, dd_obs_1to3,
                              gdd_obs_4to10, prec_obs_4to10, gddHigh_obs_4to10, dd_obs_4to10,
                              gdd_obs_fm, prec_obs_fm, gddHigh_obs_fm, dd_obs_fm,
                              gdd_obs_1to3_fm, prec_obs_1to3_fm, gddHigh_obs_1to3_fm, dd_obs_1to3_fm,
                              gdd_obs_4to10_fm, prec_obs_4to10_fm, gddHigh_obs_4to10_fm, dd_obs_4to10_fm)
                
#----------------------------#
#### 1) Bootstrap program ####
#----------------------------#

my.boot <- function(data, nrep, cluster, nCores, actual) {
# Set up the bootstrap 
B <- nrep #Number of draws
registerDoFuture() # Initialize parallel computing
plan(multisession, workers = nCores) # Define parallel computing plan with number of cores 

res <- foreach(1:B, .combine = rbind) %dorng% {
  
  if (actual==FALSE) {
    
    cluster_id <- unique(cluster)
    sb_ID <- sample(cluster_id, replace = TRUE) # For clustering at farm-level, select "ID_farm". Otherwise, select "ID_nuts2"
    data_boot <- list()
    for (j in 1:length(sb_ID)) {
      data_boot[[j]] <- df_farm[which(cluster == sb_ID[j]), ] 
    }
    data_boot <- do.call(rbind, data_boot)
    
  } else if (actual==TRUE) {
    data_boot <- df_farm
  }
    
    tryCatch( { # in case there is an optimization error in one of the draws
      
# -------------------- #
# START OF ESTIMATIONS # 
# -------------------- #
  
  library(dplyr)

  #calculate sample mean for elasticity evaluation
  dat_sm <- dplyr::summarise_all(data_boot, mean) 
  
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
    
  # Define system of equations
  system <- list( QQcereals = eqQQcereals, 
                  QQprotein = eqQQprotein,
                  QQoilseed = eqQQoilseed,
                  QQroots = eqQQroots,
                  QQcorn = eqQQcorn,
                  NXfert = eqNXfert)
    
  # Define restrictions
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

#----------------------------------#
#### Step 1: Probit regressions ####
#----------------------------------#
  
  # Obtain data mean for probit-variables (same for all probit regressions) 
  
    iCrop <- glm(cereal.eq, family = binomial(link = "probit"), 
                    data = data_boot)
        
    dat_prob <- dplyr::as_tibble(model.matrix(iCrop))
    dat_sm_prob <- dat_prob %>% 
      dplyr::summarise_all(mean)
  
  #set-up the loop over crops
  list_crops <- data.frame(crops=c("Cereals","Protein","Oilseed","Roots","Corn"),
                           probit=c("cereal.eq","protein.eq","oilseed.eq","roots.eq","corn.eq"),
                           abbrev=c("cer","prot","oil","roots","corn"))
  n_crops <- list_crops %>% count() %>% as.numeric() # number of crops
    
  # Create list where the results will be stored
  list_probit_act <- list(Cereals=list(linpred_sm="", PHI="", PHI_sm="", phi="", phi_sm="", probit_me=""),
                          Protein=list(linpred_sm="", PHI="", PHI_sm="", phi="", phi_sm="", probit_me=""),
                          Oilseed=list(linpred_sm="", PHI="", PHI_sm="", phi="", phi_sm="", probit_me=""),
                          Roots=  list(linpred_sm="", PHI="", PHI_sm="", phi="", phi_sm="", probit_me=""),
                          Corn=list(linpred_sm="", PHI="", PHI_sm="", phi="", phi_sm="", probit_me=""))
  
  # Loop over all crops and estimate probits 
    
    for (i in 1:n_crops) { 
      
      # Probit estimation
      crop.eq <- get(list_crops[i,"probit"])
      iCrop <- glm(crop.eq, family = binomial(link = "probit"), 
                   data = data_boot)
  
      # Linear prediction at the sample mean (=linpred_sm)
      list_probit_act[[i]][[1]] <- iCrop_linpred_sm <- crossprod(matrix(iCrop$coefficients),
                                  t(dat_sm_prob))
      
      # PHI at the farm-level
      list_probit_act[[i]][[2]] <- pnorm(predict(iCrop)) #note: predict(iCrop) is iCrop_linpred
      
      # PHI at the sample mean
      list_probit_act[[i]][[3]] <- pnorm(iCrop_linpred_sm)
      
      # phi at the farm-level
      list_probit_act[[i]][[4]] <- dnorm(predict(iCrop)) #note: predict(iCrop) is iCrop_linpred
      
      # phi at the sample mean
      list_probit_act[[i]][[5]] <- dnorm(iCrop_linpred_sm)
      
      # Probit coefficients
      if (i == 1) {
        coef_iCereals <- iCrop$coefficients
      } else if (i == 2) {
        coef_iProtein <- iCrop$coefficients
      } else if (i == 3) {
        coef_iOilseed <- iCrop$coefficients
      } else if (i == 4) {
        coef_iRoots <- iCrop$coefficients
      } else if (i == 5) {
        coef_iCorn <- iCrop$coefficients
      }
      
      # Marginal effects (iCrop_me)
      list_probit_act[[i]][[6]] <- dnorm(iCrop_linpred_sm)*coef(iCrop) #note: dnorm(iCrop_linpred_sm) is phi_crop_sm
    }
  
  # Restore the results estimated in the loop
    
    # Probit linear prediction at the sample mean (for elasticity evaluation)
    iCereals_linpred_sm <- as.numeric(do.call(c, list_probit_act[[1]][1]))
    iProtein_linpred_sm <- as.numeric(do.call(c, list_probit_act[[2]][1]))
    iOilseed_linpred_sm <- as.numeric(do.call(c, list_probit_act[[3]][1]))
    iRoots_linpred_sm   <- as.numeric(do.call(c, list_probit_act[[4]][1]))
    iCorn_linpred_sm    <- as.numeric(do.call(c, list_probit_act[[5]][1]))
    
    # PHIs at the farm-level (for estimation of second stage)
    data_boot$PHI_cereals <- as.numeric(do.call(c, list_probit_act[[1]][2]))
    data_boot$PHI_protein <- as.numeric(do.call(c, list_probit_act[[2]][2]))
    data_boot$PHI_oilseed <- as.numeric(do.call(c, list_probit_act[[3]][2]))
    data_boot$PHI_roots   <- as.numeric(do.call(c, list_probit_act[[4]][2]))
    data_boot$PHI_corn <- as.numeric(do.call(c, list_probit_act[[5]][2]))
    
    # PHIs at the sample mean (for elasticity evaluation)
    PHI_cereals_sm  <- as.numeric(do.call(c, list_probit_act[[1]][3]))
    PHI_protein_sm  <- as.numeric(do.call(c, list_probit_act[[2]][3]))
    PHI_oilseed_sm  <- as.numeric(do.call(c, list_probit_act[[3]][3]))
    PHI_roots_sm    <- as.numeric(do.call(c, list_probit_act[[4]][3]))
    PHI_corn_sm  <- as.numeric(do.call(c, list_probit_act[[5]][3]))
    
    # phis at the farm-level (for estimation of second stage)
    data_boot$phi_cereals <- as.numeric(do.call(c, list_probit_act[[1]][4]))
    data_boot$phi_protein <- as.numeric(do.call(c, list_probit_act[[2]][4]))
    data_boot$phi_oilseed <- as.numeric(do.call(c, list_probit_act[[3]][4]))
    data_boot$phi_roots   <- as.numeric(do.call(c, list_probit_act[[4]][4]))
    data_boot$phi_corn <- as.numeric(do.call(c, list_probit_act[[5]][4]))
    
    # phis at the sample mean (for elasticity evaluation)
    phi_cereals_sm  <- as.numeric(do.call(c, list_probit_act[[1]][5]))
    phi_protein_sm  <- as.numeric(do.call(c, list_probit_act[[2]][5]))
    phi_oilseed_sm  <- as.numeric(do.call(c, list_probit_act[[3]][5]))
    phi_roots_sm    <- as.numeric(do.call(c, list_probit_act[[4]][5]))
    phi_corn_sm  <- as.numeric(do.call(c, list_probit_act[[5]][5]))
    
    # Probit marginal effects (for results presentation)
    iCereals_me <- as.numeric(do.call(c, list_probit_act[[1]][6]))
    iProtein_me <- as.numeric(do.call(c, list_probit_act[[2]][6]))
    iOilseed_me <- as.numeric(do.call(c, list_probit_act[[3]][6]))
    iRoots_me   <- as.numeric(do.call(c, list_probit_act[[4]][6]))
    iCorn_me    <- as.numeric(do.call(c, list_probit_act[[5]][6]))
    
#----------------------------------------#
#### Step 2: Run Structural equations ####
#----------------------------------------#
  
  ## Regression with iterated SUR estimation
  model_linear <- systemfit::systemfit( formula = system, method = "SUR", 
                         data = data_boot, restrict.matrix = restrict, 
                         maxit = 100 )
    
  coef_linear <- coef(model_linear)  
  
  # "x times \beta" (needed for elasticities)
  qx_pred <- predict(model_linear)
  
  qq_cereals_pred_b <- mean(qx_pred$QQcereals.pred) - model_linear$coefficients["QQcereals_phi_cereals"]*mean(data_boot$phi_cereals) 
  qq_cereals_pred <- (qq_cereals_pred_b) / mean(data_boot$PHI_cereals)
  
  qq_protein_pred_b <- mean(qx_pred$QQprotein.pred) - model_linear$coefficients["QQprotein_phi_protein"]*mean(data_boot$phi_protein) 
  qq_protein_pred <- (qq_protein_pred_b) / mean(data_boot$PHI_protein)
  
  qq_oilseed_pred_b <- mean(qx_pred$QQoilseed.pred) - model_linear$coefficients["QQoilseed_phi_oilseed"]*mean(data_boot$phi_oilseed) 
  qq_oilseed_pred <- (qq_oilseed_pred_b) / mean(data_boot$PHI_oilseed)
  
  qq_roots_pred_b <- mean(qx_pred$QQroots.pred) - model_linear$coefficients["QQroots_phi_roots"]*mean(data_boot$phi_roots) 
  qq_roots_pred <- (qq_roots_pred_b) / mean(data_boot$PHI_roots)
  
  qq_corn_pred_b <- mean(qx_pred$QQcorn.pred) - model_linear$coefficients["QQcorn_phi_corn"]*mean(data_boot$phi_corn) 
  qq_corn_pred <- (qq_corn_pred_b) / mean(data_boot$PHI_corn)
  
  nx_fert_pred <- mean(qx_pred$NXfert.pred)
  
  # -------------------------------- #
  #### Compute price elasticities ####
  # ---------------------------------#
  
    #set-up a data frame to store the results
    elast_prices <- data.frame("Variables"=c("P Cereals", "P Protein", "P Oilseed", "P Roots", "P Corn", "W Fertilizer", "W Others"),
                               "Q Cereals"=rep(NA,7),
                               "Q Protein"=rep(NA,7),
                               "Q Oilseed"=rep(NA,7),
                               "Q Roots"=rep(NA,7),
                               "Q Corn"=rep(NA,7),
                               "X Fertilizer"=rep(NA,7),
                               "X Others"=rep(NA,7))  
    
    # ------------------------------ #
    # Price elasticities for outputs #
    # ------------------------------ #
    
    #set-up the loop over crops
    list_qcrops <- data.frame(crops=c("cereals","protein","oilseed","roots","corn"),
                              probit=c("iCereals","iProtein","iOilseed","iRoots","iCorn"),
                              struct=c("QQcereals","QQprotein", "QQoilseed", "QQroots", "QQcorn"))
    n_qcrops <- list_qcrops %>% count() %>% as.numeric() # number of crops
    
    list_prices <- data.frame(prices=c("np_cereals", "np_protein", "np_oilseed", "np_roots", "np_corn", "nw_fert"),
                              endswith=c("np_cereals)", "np_protein)", "np_oilseed)", "np_roots)", "np_corn)", "nw_fert)"))
    n_prices <- list_prices %>% count() %>% as.numeric()
    
    
    # Loop over all crops 
    for (i in 1:n_qcrops) { 
      
      coefs_struct_crop <- dplyr::as_tibble(t(coef_linear)) %>% select(starts_with(list_qcrops[i,"struct"])) 
      coefs_probit_crop <- dplyr::as_tibble(t(get(paste0("coef_",list_qcrops[i,"probit"]))))
      coef_phi_struct <- coefs_struct_crop %>% select(ends_with(paste0("phi_",list_qcrops[i,"crops"]), ignore.case = FALSE)) # phi_cereals etc.
      PHI <- get(paste0("PHI_",list_qcrops[i,"crops"],"_sm"))
      phi <- get(paste0("phi_",list_qcrops[i,"crops"],"_sm"))
      linpred <- get(paste0(list_qcrops[i,"probit"],"_linpred_sm"))
      pred <- get(paste0("qq_",list_qcrops[i,"crops"],"_pred"))
      mean_quantity <- dat_sm %>% select(paste0("qq_",list_qcrops[i,"crops"]))
      
      for (j in 1:n_prices) {
        
        coef_probit_p <- coefs_probit_crop %>% select(list_prices[j,"prices"])
        coef_struct_p <- coefs_struct_crop %>% select(ends_with(list_prices[j,"endswith"]))
        mean_price <- dat_sm %>% select(list_prices[j,"prices"])
        
        elast_prices[j,1+i] <- (PHI * coef_struct_p + phi * pred * coef_probit_p - 
                                  coef_phi_struct * linpred * coef_probit_p * phi ) * (mean_price / mean_quantity )  
      }
      
      # Price elasticity for numeraire adds up to 0
      elast_prices[j+1,1+i] <- 0 - elast_prices[1,1+i] - elast_prices[2,1+i] - elast_prices[3,1+i] - elast_prices[4,1+i] - elast_prices[5,1+i] - elast_prices[6,1+i]
    }
    
    
    # ----------------------------- #
    # Price elasticities for inputs #
    # ----------------------------- #
    
    # Fertilizer demand (note: There is no selection equation for fertilizer --> use regular derivative)
    coef_fert <- dplyr::as_tibble(t(coef_linear)) %>% select(starts_with("NXfert")) 
    
    for (j in 1:n_prices) {
      
      coef_struct_p <- coef_fert %>% select(ends_with(list_prices[j,"prices"]))
      mean_price <- dat_sm %>% select(list_prices[j,"prices"])
      
      elast_prices[j,7] <- coef_struct_p * (mean_price/dat_sm$nx_fert)
    }
    
    elast_prices[j+1,7] <- 0 - elast_prices[1,7] - elast_prices[2,7] - elast_prices[3,7] - elast_prices[4,7] - elast_prices[5,7] - elast_prices[6,7]
    
    
    # Other input demand
    list_prices <- data.frame(prices=c("np_cereals", "np_protein", "np_oilseed", "np_roots", "np_corn", "nw_fert"),
                              endswith=c("np_cereals)", "np_protein)", "np_oilseed)", "np_roots)", "np_corn)", "nw_fert)"),
                              quantities=c("qq_cereals", "qq_protein", "qq_oilseed", "qq_roots", "qq_corn", "x_fert")) # Note: Here I use x_fert instead of nx_fert so that I do not have to add a negative sign to X.Others_W.Fert
    n_prices <- list_prices %>% count() %>% as.numeric()
    
    for (j in 1:n_prices) {
      
      elast_sum <- sum(elast_prices[1:n_prices,j+1]) # sum of all elasticities except for w_others
      mean_quantity <- dat_sm %>% select(list_prices[j,"quantities"]) 
      mean_price <- dat_sm %>% select(list_prices[j,"prices"]) 
      
      elast_prices[j,8] <- (mean_quantity * mean_price) / dat_sm$x_otherinp * 
        (elast_sum)
    }
    
    elast_prices[j+1,8] <- 0 - elast_prices[1,8] - elast_prices[2,8] - elast_prices[3,8] - elast_prices[4,8] - elast_prices[5,8] - elast_prices[6,8]
    
    
    # ------------------------------------- #
    # Store price elasticities in bootstrap #
    # --------------------------------------#
    
    #q_cer
    el_qcer_pcer <- c('el_qcer_pcer'            = as.numeric(elast_prices[1,2]))
    el_qcer_pprot <- c('el_qcer_pprot'          = as.numeric(elast_prices[2,2]))
    el_qcer_poil <- c('el_qcer_poil'            = as.numeric(elast_prices[3,2]))
    el_qcer_proots <- c('el_qcer_proots'        = as.numeric(elast_prices[4,2]))
    el_qcer_pcorn <- c('el_qcer_pcorn'    = as.numeric(elast_prices[5,2]))
    el_qcer_wfert <- c('el_qcer_wfert'          = as.numeric(elast_prices[6,2]))
    el_qcer_wotherinp <- c('el_qcer_wotherinp'  = as.numeric(elast_prices[7,2]))
    
    #q_prot
    el_qprot_pcer <- c('el_qprot_pcer'            = as.numeric(elast_prices[1,3]))
    el_qprot_pprot <- c('el_qprot_pprot'          = as.numeric(elast_prices[2,3]))
    el_qprot_poil <- c('el_qprot_poil'            = as.numeric(elast_prices[3,3]))
    el_qprot_proots <- c('el_qprot_proots'        = as.numeric(elast_prices[4,3]))
    el_qprot_pcorn <- c('el_qprot_pcorn'    = as.numeric(elast_prices[5,3]))
    el_qprot_wfert <- c('el_qprot_wfert'          = as.numeric(elast_prices[6,3]))
    el_qprot_wotherinp <- c('el_qprot_wotherinp'  = as.numeric(elast_prices[7,3]))
    
    #q_oil
    el_qoil_pcer <- c('el_qoil_pcer'            = as.numeric(elast_prices[1,4]))
    el_qoil_pprot <- c('el_qoil_pprot'          = as.numeric(elast_prices[2,4]))
    el_qoil_poil <- c('el_qoil_poil'            = as.numeric(elast_prices[3,4]))
    el_qoil_proots <- c('el_qoil_proots'        = as.numeric(elast_prices[4,4]))
    el_qoil_pcorn <- c('el_qoil_pcorn'    = as.numeric(elast_prices[5,4]))
    el_qoil_wfert <- c('el_qoil_wfert'          = as.numeric(elast_prices[6,4]))
    el_qoil_wotherinp <- c('el_qoil_wotherinp'  = as.numeric(elast_prices[7,4]))
    
    #qroots
    el_qroots_pcer <- c('el_qroots_pcer'            = as.numeric(elast_prices[1,5]))
    el_qroots_pprot <- c('el_qroots_pprot'          = as.numeric(elast_prices[2,5]))
    el_qroots_poil <- c('el_qroots_poil'            = as.numeric(elast_prices[3,5]))
    el_qroots_proots <- c('el_qroots_proots'        = as.numeric(elast_prices[4,5]))
    el_qroots_pcorn <- c('el_qroots_pcorn'    = as.numeric(elast_prices[5,5]))
    el_qroots_wfert <- c('el_qroots_wfert'          = as.numeric(elast_prices[6,5]))
    el_qroots_wotherinp <- c('el_qroots_wotherinp'  = as.numeric(elast_prices[7,5]))
    
    #qcorn
    el_qcorn_pcer <- c('el_qcorn_pcer'            = as.numeric(elast_prices[1,6]))
    el_qcorn_pprot <- c('el_qcorn_pprot'          = as.numeric(elast_prices[2,6]))
    el_qcorn_poil <- c('el_qcorn_poil'            = as.numeric(elast_prices[3,6]))
    el_qcorn_proots <- c('el_qcorn_proots'        = as.numeric(elast_prices[4,6]))
    el_qcorn_pcorn <- c('el_qcorn_pcorn'    = as.numeric(elast_prices[5,6]))
    el_qcorn_wfert <- c('el_qcorn_wfert'          = as.numeric(elast_prices[6,6]))
    el_qcorn_wotherinp <- c('el_qcorn_wotherinp'  = as.numeric(elast_prices[7,6]))
    
    #xfert
    el_xfert_pcer <- c('el_xfert_pcer'            = as.numeric(elast_prices[1,7]))
    el_xfert_pprot <- c('el_xfert_pprot'          = as.numeric(elast_prices[2,7]))
    el_xfert_poil <- c('el_xfert_poils'           = as.numeric(elast_prices[3,7]))
    el_xfert_proots <- c('el_xfert_proots'        = as.numeric(elast_prices[4,7]))
    el_xfert_pcorn <- c('el_xfert_pcorn'    = as.numeric(elast_prices[5,7]))
    el_xfert_wfert <- c('el_xfert_wfert'          = as.numeric(elast_prices[6,7]))
    el_xfert_wotherinp <- c('el_xfert_woth'       = as.numeric(elast_prices[7,7]))
    
    #xother
    el_xotherinp_pcer <- c('el_xotherinp_pcer'            = as.numeric(elast_prices[1,8]))
    el_xotherinp_pprot <- c('el_xotherinp_pprot'          = as.numeric(elast_prices[2,8]))
    el_xotherinp_poil <- c('el_xotherinp_poil'            = as.numeric(elast_prices[3,8]))
    el_xotherinp_proots <- c('el_xotherinp_proots'        = as.numeric(elast_prices[4,8]))
    el_xotherinp_pcorn <- c('el_xotherinp_pcorn'    = as.numeric(elast_prices[5,8]))
    el_xotherinp_wfert <- c('el_xotherinp_wfert'          = as.numeric(elast_prices[6,8]))
    el_xotherinp_wotherinp <- c('el_xotherinp_wotherinp'  = as.numeric(elast_prices[7,8]))
    

  # ---------------------------------- #
  #### Compute weather elasticities ####
  # ---------------------------------- #  
  
    # ------------------------------ #
    # Weather elasticities for crops #
    # ------------------------------ #
    
    #set-up a data frame to store the results
    elast_weather <- data.frame("Variables"=c("GDD_obs", "PREC_obs", "GDDHigh_obs", "DD_obs", "GDD_past", "PREC_past", "GDDHigh_past", "DD_past"),
                                "Q Cereals"=rep(NA,8),
                                "Q Protein"=rep(NA,8),
                                "Q Oilseed"=rep(NA,8),
                                "Q Roots"=rep(NA,8),
                                "Q Corn"=rep(NA,8),
                                "X Fertilizer"=rep(NA,8))  
    
    # Weather vars
    list_weather <- data.frame(weather_obs=  c("gdd_obs","prec_obs","gddHigh_obs","dd_obs"),
                               weather_1to3= c("gdd_1to3","prec_1to3","gddHigh_1to3","dd_1to3"),
                               weather_4to10=c("gdd_4to10","prec_4to10","gddHigh_4to10","dd_4to10"),
                               weather_obs_1to3 = c("gdd_obs_1to3","prec_obs_1to3","gddHigh_obs_1to3","dd_obs_1to3"),
                               weather_obs_4to10 = c("gdd_obs_4to10","prec_obs_4to10","gddHigh_obs_4to10","dd_obs_4to10"),
                               weather_obs_endwith=  c("* gdd_obs)","* prec_obs)","* gddHigh_obs)","* dd_obs)"),
                               weather_1to3_endwith= c("* gdd_1to3)","* prec_1to3)","* gddHigh_1to3)","* dd_1to3)"),
                               weather_4to10_endwith=c("* gdd_4to10)","* prec_4to10)","* gddHigh_4to10)","* dd_4to10)"),
                               weather_obs_1to3_endwith = c("* gdd_obs_1to3)","* prec_obs_1to3)","* gddHigh_obs_1to3)","* dd_obs_1to3)"),
                               weather_obs_4to10_endwith = c("* gdd_obs_4to10)","* prec_obs_4to10)","* gddHigh_obs_4to10)","* dd_obs_4to10)"))
    
    n_weather <- list_weather %>% count() %>% as.numeric() 
    
    
    # Loop over all crops 
    for (i in 1:n_qcrops) { 
      
      coefs_struct_crop <- dplyr::as_tibble(t(coef_linear)) %>% select(starts_with(list_qcrops[i,"struct"])) 
      coefs_probit_crop <- dplyr::as_tibble(t(get(paste0("coef_",list_qcrops[i,"probit"]))))
      coef_phi_struct <- coefs_struct_crop %>% select(ends_with(paste0("phi_",list_qcrops[i,"crops"]), ignore.case = FALSE)) # phi_cereals etc.
      PHI <- get(paste0("PHI_",list_qcrops[i,"crops"],"_sm"))
      phi <- get(paste0("phi_",list_qcrops[i,"crops"],"_sm"))
      linpred <- get(paste0(list_qcrops[i,"probit"],"_linpred_sm"))
      pred <- get(paste0("qq_",list_qcrops[i,"crops"],"_pred"))
      mean_quantity <- dat_sm %>% select(paste0("qq_",list_qcrops[i,"crops"]))
      
      # Loop over all weather variables
      for (j in 1:n_weather) {
        
        # Weather-interaction terms
        coef_struct_weather_obs_1to3 <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_obs_1to3_endwith"]))
        coef_struct_weather_obs_4to10 <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_obs_4to10_endwith"]))
        
        # Mean weather for derivatives 
        mean_weather_obs <- dat_sm %>% select(list_weather[j,"weather_obs"])
        mean_weather_1to3 <- dat_sm %>% select(list_weather[j,"weather_1to3"])
        mean_weather_4to10 <- dat_sm %>% select(list_weather[j,"weather_4to10"])
        
        # Observed weather (note: observed weather only affects structural equations, not crop selection; also note that PHI_* is treated like a constant term )
        coef_struct_weather <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_obs_endwith"]))
        elast_weather[j,i+1] <- PHI * (coef_struct_weather + coef_struct_weather_obs_1to3*mean_weather_1to3 + coef_struct_weather_obs_4to10*mean_weather_4to10)  * (1/mean_quantity)
        
        # Weather_1to3
        coef_probit_weather <- coefs_probit_crop %>% select(list_weather[j,"weather_1to3"])
        coef_struct_weather <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_1to3_endwith"]))
        
        elast_weather_1to3 <- (PHI * (coef_struct_weather + coef_struct_weather_obs_1to3*mean_weather_obs) + phi * pred * coef_probit_weather -
                                 coef_phi_struct * linpred * coef_probit_weather * phi) * (1 / mean_quantity )
        
        # Weather_4to10
        coef_probit_weather <- coefs_probit_crop %>% select(list_weather[j,"weather_4to10"])
        coef_struct_weather <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_4to10_endwith"]))
        
        elast_weather_4to10 <- (PHI * (coef_struct_weather + coef_struct_weather_obs_4to10*mean_weather_obs) + phi * pred * coef_probit_weather -
                                  coef_phi_struct * linpred * coef_probit_weather * phi) * (1 / mean_quantity )
        
        # Weather_past
        elast_weather[j+4,i+1] <- elast_weather_1to3 + elast_weather_4to10
        
      }
      
    }
    
    # ----------------------------------- #
    # Weather elasticities for fertilizer #
    # ----------------------------------- #      
    
    # Weather vars
    list_weather <- data.frame(weather_obs=  c("gdd_obs","prec_obs","gddHigh_obs","dd_obs"),
                               weather_1to3= c("gdd_1to3","prec_1to3","gddHigh_1to3","dd_1to3"),
                               weather_4to10=c("gdd_4to10","prec_4to10","gddHigh_4to10","dd_4to10"),
                               weather_obs_1to3 = c("gdd_obs_1to3","prec_obs_1to3","gddHigh_obs_1to3","dd_obs_1to3"),
                               weather_obs_4to10 = c("gdd_obs_4to10","prec_obs_4to10","gddHigh_obs_4to10","dd_obs_4to10"),
                               weather_obs_endwith=  c("_gdd_obs","_prec_obs","_gddHigh_obs","_dd_obs"),
                               weather_1to3_endwith= c("_gdd_1to3","_prec_1to3","_gddHigh_1to3","_dd_1to3"),
                               weather_4to10_endwith=c("_gdd_4to10","_prec_4to10","_gddHigh_4to10","_dd_4to10"),
                               weather_obs_1to3_endwith = c("_gdd_obs_1to3","_prec_obs_1to3","_gddHigh_obs_1to3","_dd_obs_1to3"),
                               weather_obs_4to10_endwith = c("_gdd_obs_4to10","_prec_obs_4to10","_gddHigh_obs_4to10","_dd_obs_4to10"))
    
    n_weather <- list_weather %>% count() %>% as.numeric() 
    
    coefs_fert <- dplyr::as_tibble(t(coef_linear)) %>% select(starts_with("NXfert")) 
    
    
    # Loop over all weather variables
    for (j in 1:n_weather) {
      
      # Weather-interaction terms
      coef_fert_weather_obs_1to3 <- coefs_fert %>% select(ends_with(list_weather[j,"weather_obs_1to3_endwith"]))
      coef_fert_weather_obs_4to10 <- coefs_fert %>% select(ends_with(list_weather[j,"weather_obs_4to10_endwith"]))
      
      # Mean weather for derivatives 
      mean_weather_obs <- dat_sm %>% select(list_weather[j,"weather_obs"])
      mean_weather_1to3 <- dat_sm %>% select(list_weather[j,"weather_1to3"])
      mean_weather_4to10 <- dat_sm %>% select(list_weather[j,"weather_4to10"])
      
      # Observed weather (note: observed weather only affects structural equations, not crop selection; also note that PHI_* is treated like a constant term )
      coef_fert_weather <- coefs_fert %>% select(ends_with(list_weather[j,"weather_obs_endwith"]))
      elast_weather[j,7] <- (coef_fert_weather + coef_fert_weather_obs_1to3*mean_weather_1to3 + coef_fert_weather_obs_4to10*mean_weather_4to10) * (1/dat_sm$nx_fert)
      
      # Weather_1to3
      coef_fert_weather <- coefs_fert %>% select(ends_with(list_weather[j,"weather_1to3_endwith"]))
      elast_weather_1to3 <- (coef_fert_weather + coef_fert_weather_obs_1to3*mean_weather_obs) * (1/dat_sm$nx_fert)
      
      # Weather_4to10
      coef_fert_weather <- coefs_fert %>% select(ends_with(list_weather[j,"weather_4to10_endwith"]))
      elast_weather_4to10 <- (coef_fert_weather + coef_fert_weather_obs_4to10*mean_weather_obs) * (1/dat_sm$nx_fert)
      
      # Weather_past
      elast_weather[j+4,7] <- elast_weather_1to3 + elast_weather_4to10
      
    } 
    
    # --------------------------------------- #
    # Store weather elasticities in bootstrap #
    # ----------------------------------------#
    
      # Rescale results to obtain the units days, mm, days, days
      elast_weather[,2:7] <- elast_weather[,2:7]*100  
    
      #qcer
      
        # obs
        sel_qcer_gdd_obs     <- c('sel_qcer_gdd_obs'      = as.numeric(elast_weather[1,2]))
        sel_qcer_prec_obs    <- c('sel_qcer_prec_obs'     = as.numeric(elast_weather[2,2]))
        sel_qcer_gddHigh_obs <- c('sel_qcer_gddHigh_obs'  = as.numeric(elast_weather[3,2]))
        sel_qcer_dd_obs      <- c('sel_qcer_dd_obs'       = as.numeric(elast_weather[4,2]))
      
        #past
        sel_qcer_gdd_past <- c('sel_qcer_gdd_past'          = as.numeric(elast_weather[5,2]))
        sel_qcer_prec_past  <- c('sel_qcer_prec_past'       = as.numeric(elast_weather[6,2]))
        sel_qcer_gddHigh_past  <- c('sel_qcer_gddHigh_past' = as.numeric(elast_weather[7,2]))
        sel_qcer_dd_past  <- c('sel_qcer_dd_past'           = as.numeric(elast_weather[8,2]))
        
      #qprot
          
        # obs
        sel_qprot_gdd_obs     <- c('sel_qprot_gdd_obs'      = as.numeric(elast_weather[1,3]))
        sel_qprot_prec_obs    <- c('sel_qprot_prec_obs'     = as.numeric(elast_weather[2,3]))
        sel_qprot_gddHigh_obs <- c('sel_qprot_gddHigh_obs'  = as.numeric(elast_weather[3,3]))
        sel_qprot_dd_obs      <- c('sel_qprot_dd_obs'       = as.numeric(elast_weather[4,3]))
        
        #past
        sel_qprot_gdd_past <- c('sel_qprot_gdd_past'          = as.numeric(elast_weather[5,3]))
        sel_qprot_prec_past  <- c('sel_qprot_prec_past'       = as.numeric(elast_weather[6,3]))
        sel_qprot_gddHigh_past  <- c('sel_qprot_gddHigh_past' = as.numeric(elast_weather[7,3]))
        sel_qprot_dd_past  <- c('sel_qprot_dd_past'           = as.numeric(elast_weather[8,3]))
        
      #qoil
      
        # obs
        sel_qoil_gdd_obs     <- c('sel_qoil_gdd_obs'      = as.numeric(elast_weather[1,4]))
        sel_qoil_prec_obs    <- c('sel_qoil_prec_obs'     = as.numeric(elast_weather[2,4]))
        sel_qoil_gddHigh_obs <- c('sel_qoil_gddHigh_obs'  = as.numeric(elast_weather[3,4]))
        sel_qoil_dd_obs      <- c('sel_qoil_dd_obs'       = as.numeric(elast_weather[4,4]))
        
        #past
        sel_qoil_gdd_past <- c('sel_qoil_gdd_past'          = as.numeric(elast_weather[5,4]))
        sel_qoil_prec_past  <- c('sel_qoil_prec_past'       = as.numeric(elast_weather[6,4]))
        sel_qoil_gddHigh_past  <- c('sel_qoil_gddHigh_past' = as.numeric(elast_weather[7,4]))
        sel_qoil_dd_past  <- c('sel_qoil_dd_past'           = as.numeric(elast_weather[8,4]))
        
    
      #qroots
      
        # obs
        sel_qroots_gdd_obs     <- c('sel_qroots_gdd_obs'      = as.numeric(elast_weather[1,5]))
        sel_qroots_prec_obs    <- c('sel_qroots_prec_obs'     = as.numeric(elast_weather[2,5]))
        sel_qroots_gddHigh_obs <- c('sel_qroots_gddHigh_obs'  = as.numeric(elast_weather[3,5]))
        sel_qroots_dd_obs      <- c('sel_qroots_dd_obs'       = as.numeric(elast_weather[4,5]))
        
        #past
        sel_qroots_gdd_past <- c('sel_qroots_gdd_past'          = as.numeric(elast_weather[5,5]))
        sel_qroots_prec_past  <- c('sel_qroots_prec_past'       = as.numeric(elast_weather[6,5]))
        sel_qroots_gddHigh_past  <- c('sel_qroots_gddHigh_past' = as.numeric(elast_weather[7,5]))
        sel_qroots_dd_past  <- c('sel_qroots_dd_past'           = as.numeric(elast_weather[8,5]))
        
      #corn
          
        # obs
        sel_qcorn_gdd_obs     <- c('sel_qcorn_gdd_obs'      = as.numeric(elast_weather[1,6]))
        sel_qcorn_prec_obs    <- c('sel_qcorn_prec_obs'     = as.numeric(elast_weather[2,6]))
        sel_qcorn_gddHigh_obs <- c('sel_qcorn_gddHigh_obs'  = as.numeric(elast_weather[3,6]))
        sel_qcorn_dd_obs      <- c('sel_qcorn_dd_obs'       = as.numeric(elast_weather[4,6]))
    
        #past
        sel_qcorn_gdd_past <- c('sel_qcorn_gdd_past'          = as.numeric(elast_weather[5,6]))
        sel_qcorn_prec_past  <- c('sel_qcorn_prec_past'       = as.numeric(elast_weather[6,6]))
        sel_qcorn_gddHigh_past  <- c('sel_qcorn_gddHigh_past' = as.numeric(elast_weather[7,6]))
        sel_qcorn_dd_past  <- c('sel_qcorn_dd_past'           = as.numeric(elast_weather[8,6]))
        
      #xfert
          
        # obs
        sel_xfert_gdd_obs     <- c('sel_xfert_gdd_obs'      = as.numeric(elast_weather[1,7]))
        sel_xfert_prec_obs    <- c('sel_xfert_prec_obs'     = as.numeric(elast_weather[2,7]))
        sel_xfert_gddHigh_obs <- c('sel_xfert_gddHigh_obs'  = as.numeric(elast_weather[3,7]))
        sel_xfert_dd_obs      <- c('sel_xfert_dd_obs'       = as.numeric(elast_weather[4,7]))
        
        #past
        sel_xfert_gdd_past     <- c('sel_xfert_gdd_past'     = as.numeric(elast_weather[5,7]))
        sel_xfert_prec_past    <- c('sel_xfert_prec_past'    = as.numeric(elast_weather[6,7]))
        sel_xfert_gddHigh_past <- c('sel_xfert_gddHigh_past' = as.numeric(elast_weather[7,7]))
        sel_xfert_dd_past      <- c('sel_xfert_dd_past'      = as.numeric(elast_weather[8,7]))
        
        
  # ------------------------------- #
  #### Store results of interest ####
  # ------------------------------- #    
      
    result <- c(coef_iCereals, # Probit coefficients
                coef_iProtein,
                coef_iOilseed,
                coef_iRoots,
                coef_iCorn,
                iCereals_me, # Probit marginal effects
                iProtein_me,
                iOilseed_me,
                iRoots_me,
                iCorn_me,
                model_linear$coefficients, #structural coefficients
                el_qcer_pcer, #price elasticities
                el_qcer_pprot,
                el_qcer_poil,
                el_qcer_proots,
                el_qcer_pcorn,
                el_qcer_wfert,
                el_qcer_wotherinp,
                el_qprot_pcer,
                el_qprot_pprot,
                el_qprot_poil,
                el_qprot_proots,
                el_qprot_pcorn,
                el_qprot_wfert,
                el_qprot_wotherinp,
                el_qoil_pcer,
                el_qoil_pprot,
                el_qoil_poil,
                el_qoil_proots,
                el_qoil_pcorn,
                el_qoil_wfert,
                el_qoil_wotherinp,
                el_qroots_pcer,
                el_qroots_pprot,
                el_qroots_poil,
                el_qroots_proots,
                el_qroots_pcorn,
                el_qroots_wfert,
                el_qroots_wotherinp,
                el_qcorn_pcer,
                el_qcorn_pprot,
                el_qcorn_poil,
                el_qcorn_proots,
                el_qcorn_pcorn,
                el_qcorn_wfert,
                el_qcorn_wotherinp,
                el_xfert_pcer,
                el_xfert_pprot,
                el_xfert_poil,
                el_xfert_proots,
                el_xfert_pcorn,
                el_xfert_wfert,
                el_xfert_wotherinp,
                el_xotherinp_pcer,
                el_xotherinp_pprot,
                el_xotherinp_poil,
                el_xotherinp_proots,
                el_xotherinp_pcorn,
                el_xotherinp_wfert,
                el_xotherinp_wotherinp,
                sel_qcer_gdd_obs, #weather elasticities
                sel_qcer_prec_obs,
                sel_qcer_gddHigh_obs,
                sel_qcer_dd_obs,
                sel_qprot_gdd_obs,
                sel_qprot_prec_obs,
                sel_qprot_gddHigh_obs,
                sel_qprot_dd_obs,
                sel_qoil_gdd_obs,
                sel_qoil_prec_obs,
                sel_qoil_gddHigh_obs,
                sel_qoil_dd_obs,
                sel_qroots_gdd_obs,
                sel_qroots_prec_obs,
                sel_qroots_gddHigh_obs,
                sel_qroots_dd_obs,
                sel_qcorn_gdd_obs,
                sel_qcorn_prec_obs,
                sel_qcorn_gddHigh_obs,
                sel_qcorn_dd_obs,
                sel_xfert_gdd_obs,
                sel_xfert_prec_obs,
                sel_xfert_gddHigh_obs,
                sel_xfert_dd_obs,
                sel_qcer_gdd_past,
                sel_qcer_prec_past,
                sel_qcer_gddHigh_past,
                sel_qcer_dd_past,
                sel_qprot_gdd_past,
                sel_qprot_prec_past,
                sel_qprot_gddHigh_past,
                sel_qprot_dd_past,
                sel_qoil_gdd_past,
                sel_qoil_prec_past,
                sel_qoil_gddHigh_past,
                sel_qoil_dd_past,
                sel_qroots_gdd_past,
                sel_qroots_prec_past,
                sel_qroots_gddHigh_past,
                sel_qroots_dd_past,
                sel_qcorn_gdd_past,
                sel_qcorn_prec_past,
                sel_qcorn_gddHigh_past,
                sel_qcorn_dd_past,
                sel_xfert_gdd_past,
                sel_xfert_prec_past,
                sel_xfert_gddHigh_past,
                sel_xfert_dd_past) 
  
  result 

  }, error=function(err) {rep(NA,1185)} ) # insert NA's if there is an optimization error
}

}


#-------------------------------------------#
#### Obtain and store the actual results ####
#-------------------------------------------#

# run one draw with the original data
system.time(boot_actresults <- my.boot(data=df_farm, nrep=1, cluster=df_farm$key, nCores=1, actual=TRUE))

# save actual result
save(boot_actresults, file="rOutput/bootresults_main_actual.Rda")
boot_actresults <- NULL

#------------------------------------------------#
#### Run bootstrap: Cluster at the farm-level ####
#------------------------------------------------#

# find number of cores
nCores <- parallel::detectCores() -1 

# run the bootstrap
set.seed(1234)
system.time(boot_results <- my.boot(data=df_farm, nrep=1000, cluster=df_farm$key, nCores=nCores, actual=FALSE))

# save
save(boot_results, file="rOutput/bootresults_main_clustFarm.Rda")
boot_results <- NULL

#-------------------------------------------------#
#### Run bootstrap: Cluster at the nuts2-level ####
#-------------------------------------------------#

# find number of cores
nCores <- parallel::detectCores() -1 

# run
set.seed(1234)
system.time(boot_results <- my.boot(data=df_farm, nrep=1000, cluster=df_farm$nuts2, nCores=nCores, actual=FALSE))

# save
save(boot_results, file="rOutput/bootresults_main_clustNuts2.Rda")
boot_results <- NULL



#----------------------------------#
#### 2) Presentation of results ####
#----------------------------------#

# Pick cluster at farm or nuts2 level
cluster <- "clustFarm"
#cluster <- "clustNuts2"

# Load bootstrap results
load(paste0("rOutput/bootresults_main_",cluster,".Rda"))
load("rOutput/bootresults_main_actual.Rda")

# Turn matrix into data frame
boot_results <- as.data.frame(boot_results)
summary(boot_results$sel_qcer_gdd_obs)

# --------------------------------- #
# Create list with original results #
# --------------------------------- #

names(boot_results)
boot_summary <- data.frame("Variable"=c(names(boot_results)),
                           "Original"=boot_actresults) 

# ------------------------ #
# Add confidence intervals #
# ------------------------ #

# Confidence intervals
CIs <- t(sapply(boot_results,function(i) quantile(i,c(0.025, 0.975), na.rm=TRUE)))
CIs <- as.data.frame(CIs)
colnames(CIs) <- c("L95","U95")

# Add CIs to boot_summary
boot_summary <- cbind(boot_summary, CIs)

# Organize
coef_iCereals <- boot_summary[c(1:47),]
coef_iProtein <- boot_summary[c(48:94),]
coef_iOilseed <- boot_summary[c(95:141),]
coef_iRoots <- boot_summary[c(142:188),]
coef_iCorn <- boot_summary[c(189:235),]
iCereals_me <- boot_summary[c(236:282),]
iProtein_me <- boot_summary[c(283:329),]
iOilseed_me <- boot_summary[c(330:376),]
iRoots_me <- boot_summary[c(377:423),]
iCorn_me <- boot_summary[c(424:470),]

iCereals_me_extr <- boot_summary[c(471:517),]
iProtein_me_extr <- boot_summary[c(518:564),]
iOilseed_me_extr <- boot_summary[c(565:611),]
iRoots_me_extr <- boot_summary[c(612:658),]
iCorn_me_extr <- boot_summary[c(659:705),]

qqCereals <- boot_summary[c(471:534),]
qqProtein <- boot_summary[c(535:598),]
qqOilseed <- boot_summary[c(599:662),]
qqRoots   <- boot_summary[c(663:726),]
qqCorn    <- boot_summary[c(727:790),]
nxFert    <- boot_summary[c(791:853),]

PriceElast_cereals <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "el_qcer_")),]
PriceElast_protein <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "el_qprot_")),]
PriceElast_oilseed <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "el_qoil_")),]
PriceElast_roots   <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "el_qroots_")),]
PriceElast_corn <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "el_qcorn_")),]
PriceElast_fert    <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "el_xfert_")),]
PriceElast_others  <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "el_xotherinp_")),]

WeatherElast_cereals <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qcer_")),]
WeatherElast_protein <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qprot_")),]
WeatherElast_oilseed <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qoil_")),]
WeatherElast_roots   <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qroots_")),]
WeatherElast_corn <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qcorn_")),]
WeatherElast_fert    <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_xfert_")),]


#-------------------------------------------------------------------#
#### Figure 2: Average partial effects of weather on crop choice ####
#-------------------------------------------------------------------#

df.cropchoice <- data.frame(matrix(ncol = 5, nrow = 20))
colnames(df.cropchoice) <- c("Weather", "Crop", "ME", "L95", "U95")

df.cropchoice$Weather <- rep(c("gdd","prec","gddHigh","dd"),5)


#Cereals
df.cropchoice[1,] = c("gdd","Cereals",10*iCereals_me[17,c("Original","L95", "U95")] + 
                        10*iCereals_me[21,c("Original","L95", "U95")])

df.cropchoice[2,] = c("prec","Cereals",10*iCereals_me[18,c("Original","L95", "U95")] + 
                        10*iCereals_me[22,c("Original","L95", "U95")])

df.cropchoice[3,] = c("gddHigh","Cereals",iCereals_me[19,c("Original","L95", "U95")] + 
                        iCereals_me[23,c("Original","L95", "U95")])

df.cropchoice[4,] = c("dd","Cereals",iCereals_me[20,c("Original","L95", "U95")] + 
                        iCereals_me[24,c("Original","L95", "U95")])

#Protein
df.cropchoice[5,] = c("gdd","Protein crops",10*iProtein_me[17,c("Original","L95", "U95")] + 
                        10*iProtein_me[21,c("Original","L95", "U95")])

df.cropchoice[6,] = c("prec","Protein crops",10*iProtein_me[18,c("Original","L95", "U95")] + 
                        10*iProtein_me[22,c("Original","L95", "U95")])

df.cropchoice[7,] = c("gddHigh","Protein crops",iProtein_me[19,c("Original","L95", "U95")] + 
                        iProtein_me[23,c("Original","L95", "U95")])

df.cropchoice[8,] = c("dd","Protein crops",iProtein_me[20,c("Original","L95", "U95")] + 
                        iProtein_me[24,c("Original","L95", "U95")])


#Oilseed
df.cropchoice[9,] = c("gdd","Oilseeds",10*iOilseed_me[17,c("Original","L95", "U95")] + 
                        10*iOilseed_me[21,c("Original","L95", "U95")])

df.cropchoice[10,] = c("prec","Oilseeds",10*iOilseed_me[18,c("Original","L95", "U95")] + 
                         10*iOilseed_me[22,c("Original","L95", "U95")])

df.cropchoice[11,] = c("gddHigh","Oilseeds",iOilseed_me[19,c("Original","L95", "U95")] + 
                         iOilseed_me[23,c("Original","L95", "U95")])

df.cropchoice[12,] = c("dd","Oilseeds",iOilseed_me[20,c("Original","L95", "U95")] + 
                         iOilseed_me[24,c("Original","L95", "U95")])

#Roots
df.cropchoice[13,] = c("gdd","Root crops",10*iRoots_me[17,c("Original","L95", "U95")] + 
                         10*iRoots_me[21,c("Original","L95", "U95")])

df.cropchoice[14,] = c("prec","Root crops",10*iRoots_me[18,c("Original","L95", "U95")] + 
                         10*iRoots_me[22,c("Original","L95", "U95")])

df.cropchoice[15,] = c("gddHigh","Root crops",iRoots_me[19,c("Original","L95", "U95")] + 
                         iRoots_me[23,c("Original","L95", "U95")])

df.cropchoice[16,] = c("dd","Root crops",iRoots_me[20,c("Original","L95", "U95")] + 
                         iRoots_me[24,c("Original","L95", "U95")])

#Corn
df.cropchoice[17,] = c("gdd","Corn",10*iCorn_me[17,c("Original","L95", "U95")] + 
                         10*iCorn_me[21,c("Original","L95", "U95")])

df.cropchoice[18,] = c("prec","Corn",10*iCorn_me[18,c("Original","L95", "U95")] + 
                         10*iCorn_me[22,c("Original","L95", "U95")])

df.cropchoice[19,] = c("gddHigh","Corn",iCorn_me[19,c("Original","L95", "U95")] + 
                         iCorn_me[23,c("Original","L95", "U95")])

df.cropchoice[20,] = c("dd","Corn",iCorn_me[20,c("Original","L95", "U95")] + 
                         iCorn_me[24,c("Original","L95", "U95")])


# Create Plot
df.cropchoice$Weather <- factor(df.cropchoice$Weather,
                                levels = c("gdd", "prec", "gddHigh", "dd")) # to keep the order
df.cropchoice$Crop <- factor(df.cropchoice$Crop,
                             levels = c("Cereals", "Protein crops", "Oilseeds", "Root crops", "Corn")) # to keep the order

fig_cropchoice <- df.cropchoice %>% 
  ggplot() + 
  geom_col_pattern(
    aes(x=Weather, y=ME, pattern_angle=Weather, fill=Weather),
    pattern="stripe",
    pattern_color = "white",
    color = "black", # Color of the borders
    pattern_key_scale_factor = 0.2,
    pattern_density = 0.005, #thickness of lines (default 0.2)
    pattern_spacing = 0.05 #distance between lines (default 0.05)
  ) +
  geom_errorbar(aes(x=Weather, y=ME, ymin=L95, ymax=U95),
                size=.3, # Thinner lines
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  scale_fill_manual(name = "",
                    breaks = c("gdd","prec","gddHigh","dd"),
                    values = c("#8c510a", "#d8b365", "#5ab4ac", "#01665e"),
                    labels = c("Growing degree days (10 days)", 
                               "Precipitation (cm)",
                               "High degree days",
                               "Dry days")) +
  scale_pattern_angle_manual(name = "",
                             breaks = c("gdd","prec","gddHigh","dd"),
                             values=c(45,0,135,90),
                             labels = c("Growing degree days (10 days)", 
                                        "Precipitation (cm)",
                                        "High degree days",
                                        "Dry days")) + 
  scale_y_continuous(breaks = seq(-0.10, 0.10, by = 0.025), limits=c(-0.10,0.10)) +  
  labs(x="", y="Average partial effect") +
  ggtitle("") +
  theme_bw() + 
  theme(text = element_text(size = 12),
        axis.title.y = element_text(size = 10),
        legend.position="bottom",
        legend.key.size = unit(.25, 'cm'),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.box.spacing = unit(0.0, "cm"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "cm"),
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  facet_grid(~Crop, switch = "x") 
fig_cropchoice

#save plot
if (cluster=="clustFarm") {
  
  ggsave("Figures/Figure_2.png", fig_cropchoice, width = 6.5, height = 3, units = "in", dpi=1200)
  ggsave("Figures/Figure_2.eps", fig_cropchoice, width = 6.5, height = 3, units = "in", dpi=1200)
  
} 

#--------------------------------------------------------#
#### Table A.1: Parameter estimates for probit models ####
#--------------------------------------------------------#

#combine probit results
tab_parprobit <- cbind(coef_iCereals[,c("Variable","Original","L95","U95")],
                       coef_iProtein[,c("Original","L95","U95")],
                       coef_iOilseed[,c("Original","L95","U95")],
                       coef_iRoots[,c("Original","L95","U95")],
                       coef_iCorn[,c("Original","L95","U95")])

tab_parprobit <- tab_parprobit[1:24,]

#name columns
names(tab_parprobit) <- c("Variable", 
                          "Cereals", "L95_cereals", "U95_cereals", 
                          "Protein", "L95_protein", "U95_protein", 
                          "Oilseed", "L95_oilseed", "U95_oilseed", 
                          "Roots", "L95_roots", "U95_roots", 
                          "Corn", "L95_corn", "U95_corn")

#name rows
tab_parprobit$Variable <- c("Intercept", "P Cereals", "P Protein", "P Oilseed crops", "P Root crops", "P Corn", "R Fert", 
                            "K Land", "K Labor", "K Capital", "Trend", "Trend^2", 
                            "L.Share Cereals", "L.Share Protein cr.", "L.Share Oilseeds", "L.Share Root cr.", 
                            "Gdd_1to3", "Prec_1to3", "GddHigh_1to3", "DD_1to3",
                            "Gdd_4to10", "Prec_4to10", "GddHigh_4to10", "DD_4to10")

#Numbers smaller than 0.0005 should be in scientific format; others in 0.001 etc.
tab_parprobit <- as.data.frame(cbind(tab_parprobit[,1], sapply(tab_parprobit[,2:ncol(tab_parprobit)], function(x) ifelse(abs(x)<0.0005&is.numeric(x),format(x,scientific=TRUE,digits=1),format(round(x, 3), nsmall = 3, scientific=FALSE)) )) )

# remove white space
tab_parprobit <- as.data.frame(sapply(tab_parprobit, function(x) (gsub(" ", "", x))))

# Merge CIs
tab_parprobit$L95_cereals <- paste0(paste0("(",tab_parprobit$L95_cereals,";"),
                                    paste0(tab_parprobit$U95_cereals,")"))
tab_parprobit$U95_cereals <- NULL
tab_parprobit <- dplyr::rename(tab_parprobit,"CI_cereals" = "L95_cereals")

tab_parprobit$L95_protein <- paste0(paste0("(",tab_parprobit$L95_protein,";"),
                                    paste0(tab_parprobit$U95_protein,")"))
tab_parprobit$U95_protein <- NULL
tab_parprobit <- dplyr::rename(tab_parprobit,"CI_protein" = "L95_protein")

tab_parprobit$L95_oilseed <- paste0(paste0("(",tab_parprobit$L95_oilseed,";"),
                                    paste0(tab_parprobit$U95_oilseed,")"))
tab_parprobit$U95_oilseed <- NULL
tab_parprobit <- dplyr::rename(tab_parprobit,"CI_oilseed" = "L95_oilseed")

tab_parprobit$L95_roots <- paste0(paste0("(",tab_parprobit$L95_roots,";"),
                                  paste0(tab_parprobit$U95_roots,")"))
tab_parprobit$U95_roots <- NULL
tab_parprobit <- dplyr::rename(tab_parprobit,"CI_roots" = "L95_roots")

tab_parprobit$L95_corn <- paste0(paste0("(",tab_parprobit$L95_corn,";"),
                                    paste0(tab_parprobit$U95_corn,")"))
tab_parprobit$U95_corn <- NULL
tab_parprobit <- dplyr::rename(tab_parprobit,"CI_corn" = "L95_corn")

#Insert a blank row to every other row
na.df <- data.frame("na" = NA, "Cereals" = NA, "CI_cereals" = NA, "Protein" = NA,
                    "CI_protein" = NA, "Oilseed" = NA, "CI_oilseed" = NA, "Roots" = NA,
                    "CI_roots" = NA, "Corn" = NA, "CI_corn" = NA)
tab_parprobit <- do.call(rbind, apply(tab_parprobit, 1, function(x) {rbind(x, na.df)}))

## Shift CIs down one line
tab_parprobit$CI_cereals <- c(NA,tab_parprobit$CI_cereals[1:nrow(tab_parprobit)-1])
tab_parprobit$CI_protein <- c(NA,tab_parprobit$CI_protein[1:nrow(tab_parprobit)-1])
tab_parprobit$CI_oilseed <- c(NA,tab_parprobit$CI_oilseed[1:nrow(tab_parprobit)-1])
tab_parprobit$CI_roots <- c(NA,tab_parprobit$CI_roots[1:nrow(tab_parprobit)-1])
tab_parprobit$CI_corn <- c(NA,tab_parprobit$CI_corn[1:nrow(tab_parprobit)-1])

# Shifts CIs to the left
tab_parprobit$Cereals[is.na(tab_parprobit$Cereals)] <- tab_parprobit$CI_cereals[is.na(tab_parprobit$Cereals)]
tab_parprobit$CI_cereals <- NULL

tab_parprobit$Protein[is.na(tab_parprobit$Protein)] <- tab_parprobit$CI_protein[is.na(tab_parprobit$Protein)]
tab_parprobit$CI_protein <- NULL

tab_parprobit$Oilseed[is.na(tab_parprobit$Oilseed)] <- tab_parprobit$CI_oilseed[is.na(tab_parprobit$Oilseed)]
tab_parprobit$CI_oilseed <- NULL

tab_parprobit$Roots[is.na(tab_parprobit$Roots)] <- tab_parprobit$CI_roots[is.na(tab_parprobit$Roots)]
tab_parprobit$CI_roots <- NULL

tab_parprobit$Corn[is.na(tab_parprobit$Corn)] <- tab_parprobit$CI_corn[is.na(tab_parprobit$Corn)]
tab_parprobit$CI_corn <- NULL

# Print
if (cluster=="clustFarm") {
  write_xlsx(tab_parprobit,path="Tables/Table_A1.xlsx")
} else if (cluster=="clustNuts2") {
  write_xlsx(tab_parprobit,path="Tables/Table_S4.1.xlsx")
}

#----------------------------------------------------------#
#### Table A.3: Regression results structural equations ####
#----------------------------------------------------------#

#Take original tables without Mundlak terms
tab_structural <- rbind(qqCereals,qqProtein,qqOilseed,qqRoots,qqCorn,nxFert) %>% 
  filter(!grepl("_fm",Variable)) %>% 
  select("Variable","Coeff"="Original", "L95", "U95")

#Numbers smaller than 0.005 should be in scientific format; others in 0.01 etc.
tab_structural <- as.data.frame(cbind(tab_structural[,1], sapply(tab_structural[,2:ncol(tab_structural)], function(x) ifelse(abs(x)<0.005&is.numeric(x),format(x,scientific=TRUE,digits=1),format(round(x, 2), nsmall = 2, scientific=FALSE)) )) )
tab_structural <- dplyr::rename(tab_structural, "Variable" = "V1")

# remove white space
tab_structural <- as.data.frame(sapply(tab_structural, function(x) (gsub(" ", "", x))))

#Put CIs together
tab_structural$L95 <- paste0(paste0("(",tab_structural$L95,";"),
                             paste0(tab_structural$U95,")"))
tab_structural$U95 <- NULL
tab_structural <- dplyr::rename(tab_structural,"CI" = "L95")

#Use nice variable names
tab_structural$Variable <- c("Intercept", "P Cereals", "P Protein", "P Oilseed crops", "P Root crops", "P Corn", "R Fert", 
                             "K Land", "K Labor", "K Capital", "Trend", "Trend^2", 
                             "Gdd", "Prec", "GddHigh", "DD",
                             "Gdd_1to3", "Prec_1to3", "GddHigh_1to3", "DD_1to3",
                             "Gdd_4to10", "Prec_4to10", "GddHigh_4to10", "DD_4to10",
                             "Gdd x Gdd_1to3", "Prec x Prec_1to3", "GddHigh x GddHigh_1to3", "DD x DD_1to3",
                             "Gdd x Gdd_4to10", "Prec x Prec_4to10", "GddHigh x GddHigh_4to10", "DD x DD_4to10",
                             "phi Cereals",
                             "Intercept", "P Cereals", "P Protein", "P Oilseed crops", "P Root crops", "P Corn", "R Fert", 
                             "K Land", "K Labor", "K Capital", "Trend", "Trend^2", 
                             "Gdd", "Prec", "GddHigh", "DD",
                             "Gdd_1to3", "Prec_1to3", "GddHigh_1to3", "DD_1to3",
                             "Gdd_4to10", "Prec_4to10", "GddHigh_4to10", "DD_4to10",
                             "Gdd x Gdd_1to3", "Prec x Prec_1to3", "GddHigh x GddHigh_1to3", "DD x DD_1to3",
                             "Gdd x Gdd_4to10", "Prec x Prec_4to10", "GddHigh x GddHigh_4to10", "DD x DD_4to10",
                             "phi Protein",
                             "Intercept", "P Cereals", "P Protein", "P Oilseed crops", "P Root crops", "P Corn", "R Fert", 
                             "K Land", "K Labor", "K Capital", "Trend", "Trend^2", 
                             "Gdd", "Prec", "GddHigh", "DD",
                             "Gdd_1to3", "Prec_1to3", "GddHigh_1to3", "DD_1to3",
                             "Gdd_4to10", "Prec_4to10", "GddHigh_4to10", "DD_4to10",
                             "Gdd x Gdd_1to3", "Prec x Prec_1to3", "GddHigh x GddHigh_1to3", "DD x DD_1to3",
                             "Gdd x Gdd_4to10", "Prec x Prec_4to10", "GddHigh x GddHigh_4to10", "DD x DD_4to10",
                             "phi Oilseed",
                             "Intercept", "P Cereals", "P Protein", "P Oilseed crops", "P Root crops", "P Corn", "R Fert", 
                             "K Land", "K Labor", "K Capital", "Trend", "Trend^2", 
                             "Gdd", "Prec", "GddHigh", "DD",
                             "Gdd_1to3", "Prec_1to3", "GddHigh_1to3", "DD_1to3",
                             "Gdd_4to10", "Prec_4to10", "GddHigh_4to10", "DD_4to10",
                             "Gdd x Gdd_1to3", "Prec x Prec_1to3", "GddHigh x GddHigh_1to3", "DD x DD_1to3",
                             "Gdd x Gdd_4to10", "Prec x Prec_4to10", "GddHigh x GddHigh_4to10", "DD x DD_4to10",
                             "phi Roots",
                             "Intercept", "P Cereals", "P Protein", "P Oilseed crops", "P Root crops", "P Corn", "R Fert", 
                             "K Land", "K Labor", "K Capital", "Trend", "Trend^2", 
                             "Gdd", "Prec", "GddHigh", "DD",
                             "Gdd_1to3", "Prec_1to3", "GddHigh_1to3", "DD_1to3",
                             "Gdd_4to10", "Prec_4to10", "GddHigh_4to10", "DD_4to10",
                             "Gdd x Gdd_1to3", "Prec x Prec_1to3", "GddHigh x GddHigh_1to3", "DD x DD_1to3",
                             "Gdd x Gdd_4to10", "Prec x Prec_4to10", "GddHigh x GddHigh_4to10", "DD x DD_4to10",
                             "phi Corn",
                             "Intercept", "P Cereals", "P Protein", "P Oilseed crops", "P Root crops", "P Corn", "R Fert", 
                             "K Land", "K Labor", "K Capital", "Trend", "Trend^2", 
                             "Gdd", "Prec", "GddHigh", "DD",
                             "Gdd_1to3", "Prec_1to3", "GddHigh_1to3", "DD_1to3",
                             "Gdd_4to10", "Prec_4to10", "GddHigh_4to10", "DD_4to10",
                             "Gdd x Gdd_1to3", "Prec x Prec_1to3", "GddHigh x GddHigh_1to3", "DD x DD_1to3",
                             "Gdd x Gdd_4to10", "Prec x Prec_4to10", "GddHigh x GddHigh_4to10", "DD x DD_4to10")


# Insert LHS indicators
tab_structural <- berryFunctions::insertRows(tab_structural, 1 , new = c("Cereal supply","",""))
tab_structural <- berryFunctions::insertRows(tab_structural, 2+nrow(qqCereals %>% filter(!grepl("_fm",Variable))), new = c("Protein supply","",""))
tab_structural <- berryFunctions::insertRows(tab_structural, 3+2*nrow(qqCereals %>% filter(!grepl("_fm",Variable))) , new = c("Oilseed supply","",""))
tab_structural <- berryFunctions::insertRows(tab_structural, 4+3*+nrow(qqCereals %>% filter(!grepl("_fm",Variable))) , new = c("Roots supply","",""))
tab_structural <- berryFunctions::insertRows(tab_structural, 5+4*nrow(qqCereals %>% filter(!grepl("_fm",Variable))) , new = c("Corn supply","",""))
tab_structural <- berryFunctions::insertRows(tab_structural, 6+5*nrow(qqCereals %>% filter(!grepl("_fm",Variable))) , new = c("Fertilizer demand","",""))

# Print
if (cluster=="clustFarm") {
  write_xlsx(tab_structural,path="Tables/Table_A3.xlsx")
} else if (cluster=="clustNuts2") {
  write_xlsx(tab_structural,path="Tables/Table_S4.2.xlsx")
}

#--------------------------------------------------#
#### Table 2: Own- and cross-price elasticities ####
#--------------------------------------------------#

#combine qs and xs
tab_PriceElast <- cbind(PriceElast_cereals[,c("Variable","Original","L95","U95")],
                        PriceElast_protein[,c("Original","L95","U95")],
                        PriceElast_oilseed[,c("Original","L95","U95")],
                        PriceElast_roots[,c("Original","L95","U95")],
                        PriceElast_corn[,c("Original","L95","U95")],
                        PriceElast_fert[,c("Original","L95","U95")],
                        PriceElast_others[,c("Original","L95","U95")])

#name columns
names(tab_PriceElast) <- c("na", 
                           "Q_Cereals", "L95_cereals", "U95_cereals", 
                           "Q_Protein", "L95_protein", "U95_protein", 
                           "Q_Oilseed", "L95_oilseed", "U95_oilseed", 
                           "Q_Roots", "L95_roots", "U95_roots", 
                           "Q_Corn", "L95_corn", "U95_corn", 
                           "X_Fertilizer", "L95_fert", "U95_fert", 
                           "X_Others", "L95_others", "U95_others")

#name rows
tab_PriceElast$na <- c("P_Cereals","P_Protein","P_Oilseed","P_Roots","P_Corn",
                       "X_Fertilizer", "X_Others")

#Numbers smaller than 0.0005 should be in scientific format; others in 0.001 etc.
tab_PriceElast <- as.data.frame(cbind(tab_PriceElast[,1], sapply(tab_PriceElast[,-1], function(x) ifelse(abs(x)<0.0005&is.numeric(x),format(x,scientific=TRUE,digits=1),format(round(x, 3), nsmall = 3, scientific=FALSE)) )) )

# remove white space
tab_PriceElast <- as.data.frame(sapply(tab_PriceElast, function(x) (gsub(" ", "", x))))

# Merge CIs
tab_PriceElast$L95_cereals <- paste0(paste0("(",tab_PriceElast$L95_cereals,";"),
                                     paste0(tab_PriceElast$U95_cereals,")"))
tab_PriceElast$U95_cereals <- NULL
tab_PriceElast <- dplyr::rename(tab_PriceElast,"CI_cereals" = "L95_cereals")

tab_PriceElast$L95_protein <- paste0(paste0("(",tab_PriceElast$L95_protein,";"),
                                     paste0(tab_PriceElast$U95_protein,")"))
tab_PriceElast$U95_protein <- NULL
tab_PriceElast <- dplyr::rename(tab_PriceElast,"CI_protein" = "L95_protein")

tab_PriceElast$L95_oilseed <- paste0(paste0("(",tab_PriceElast$L95_oilseed,";"),
                                     paste0(tab_PriceElast$U95_oilseed,")"))
tab_PriceElast$U95_oilseed <- NULL
tab_PriceElast <- dplyr::rename(tab_PriceElast,"CI_oilseed" = "L95_oilseed")

tab_PriceElast$L95_roots <- paste0(paste0("(",tab_PriceElast$L95_roots,";"),
                                   paste0(tab_PriceElast$U95_roots,")"))
tab_PriceElast$U95_roots <- NULL
tab_PriceElast <- dplyr::rename(tab_PriceElast,"CI_roots" = "L95_roots")

tab_PriceElast$L95_corn <- paste0(paste0("(",tab_PriceElast$L95_corn,";"),
                                     paste0(tab_PriceElast$U95_corn,")"))
tab_PriceElast$U95_corn <- NULL
tab_PriceElast <- dplyr::rename(tab_PriceElast,"CI_corn" = "L95_corn")

tab_PriceElast$L95_fert <- paste0(paste0("(",tab_PriceElast$L95_fert,";"),
                                  paste0(tab_PriceElast$U95_fert,")"))
tab_PriceElast$U95_fert <- NULL
tab_PriceElast <- dplyr::rename(tab_PriceElast,"CI_fert" = "L95_fert")

tab_PriceElast$L95_others <- paste0(paste0("(",tab_PriceElast$L95_others,";"),
                                    paste0(tab_PriceElast$U95_others,")"))
tab_PriceElast$U95_others <- NULL
tab_PriceElast <- dplyr::rename(tab_PriceElast,"CI_others" = "L95_others")

#Insert a blank row to every other row
na.df <- data.frame("na" = NA, "Q_Cereals" = NA, "CI_cereals" = NA, "Q_Protein" = NA,
                    "CI_protein" = NA, "Q_Oilseed" = NA, "CI_oilseed" = NA, "Q_Roots" = NA,
                    "CI_roots" = NA, "Q_Corn" = NA, "CI_corn" = NA, "X_Fertilizer" = NA,
                    "CI_fert" = NA, "X_Others" = NA, "CI_others" = NA)
tab_PriceElast <- do.call(rbind, apply(tab_PriceElast, 1, function(x) {rbind(x, na.df)}))

## Shift CIs down one line
tab_PriceElast$CI_cereals <- c(NA,tab_PriceElast$CI_cereals[1:nrow(tab_PriceElast)-1])
tab_PriceElast$CI_protein <- c(NA,tab_PriceElast$CI_protein[1:nrow(tab_PriceElast)-1])
tab_PriceElast$CI_oilseed <- c(NA,tab_PriceElast$CI_oilseed[1:nrow(tab_PriceElast)-1])
tab_PriceElast$CI_roots <- c(NA,tab_PriceElast$CI_roots[1:nrow(tab_PriceElast)-1])
tab_PriceElast$CI_corn <- c(NA,tab_PriceElast$CI_corn[1:nrow(tab_PriceElast)-1])
tab_PriceElast$CI_fert <- c(NA,tab_PriceElast$CI_fert[1:nrow(tab_PriceElast)-1])
tab_PriceElast$CI_others <- c(NA,tab_PriceElast$CI_others[1:nrow(tab_PriceElast)-1])

# Shifts CIs to the left
tab_PriceElast$Q_Cereals[is.na(tab_PriceElast$Q_Cereals)] <- tab_PriceElast$CI_cereals[is.na(tab_PriceElast$Q_Cereals)]
tab_PriceElast$CI_cereals <- NULL

tab_PriceElast$Q_Protein[is.na(tab_PriceElast$Q_Protein)] <- tab_PriceElast$CI_protein[is.na(tab_PriceElast$Q_Protein)]
tab_PriceElast$CI_protein <- NULL

tab_PriceElast$Q_Oilseed[is.na(tab_PriceElast$Q_Oilseed)] <- tab_PriceElast$CI_oilseed[is.na(tab_PriceElast$Q_Oilseed)]
tab_PriceElast$CI_oilseed <- NULL

tab_PriceElast$Q_Roots[is.na(tab_PriceElast$Q_Roots)] <- tab_PriceElast$CI_roots[is.na(tab_PriceElast$Q_Roots)]
tab_PriceElast$CI_roots <- NULL

tab_PriceElast$Q_Corn[is.na(tab_PriceElast$Q_Corn)] <- tab_PriceElast$CI_corn[is.na(tab_PriceElast$Q_Corn)]
tab_PriceElast$CI_corn <- NULL

tab_PriceElast$X_Fertilizer[is.na(tab_PriceElast$X_Fertilizer)] <- tab_PriceElast$CI_fert[is.na(tab_PriceElast$X_Fertilizer)]
tab_PriceElast$CI_fert <- NULL

tab_PriceElast$X_Others[is.na(tab_PriceElast$X_Others)] <- tab_PriceElast$CI_others[is.na(tab_PriceElast$X_Others)]
tab_PriceElast$CI_others <- NULL

# Print
if (cluster=="clustFarm") {
  write_xlsx(tab_PriceElast,path="Tables/Table_2.xlsx")
} else if (cluster=="clustNuts2") {
  write_xlsx(tab_PriceElast,path="Tables/Table_S4.3.xlsx")
}

#-----------------------------------------------------#
#### Table 3: Weather effects on output and inputs ####
#-----------------------------------------------------#
#(compare obs vs. past)

#combine qs and xs
tab_WeatherElast <- cbind(WeatherElast_cereals[,c("Variable","Original","L95","U95")],
                          WeatherElast_protein[,c("Original","L95","U95")],
                          WeatherElast_oilseed[,c("Original","L95","U95")],
                          WeatherElast_roots[,c("Original","L95","U95")],
                          WeatherElast_corn[,c("Original","L95","U95")],
                          WeatherElast_fert[,c("Original","L95","U95")])

#name columns
names(tab_WeatherElast) <- c("na", 
                             "Q_Cereals", "L95_cereals", "U95_cereals", 
                             "Q_Protein", "L95_protein", "U95_protein", 
                             "Q_Oilseed", "L95_oilseed", "U95_oilseed", 
                             "Q_Roots", "L95_roots", "U95_roots", 
                             "Q_Corn", "L95_corn", "U95_corn", 
                             "X_Fertilizer", "L95_fert", "U95_fert")

#name rows
tab_WeatherElast$na <- c("gdd_obs","prec_obs","gddHigh_obs","dd_obs",
                         "gdd_past","prec_past","gddHigh_past","dd_past")


#Numbers smaller than 0.0005 should be in scientific format; others in 0.001 etc.
tab_WeatherElast <- as.data.frame(cbind(tab_WeatherElast[,1], sapply(tab_WeatherElast[,-1], function(x) ifelse(abs(x)<0.0005&is.numeric(x),format(x,scientific=TRUE,digits=1),format(round(x, 3), nsmall = 3, scientific=FALSE)) )) )

# remove white space
tab_WeatherElast <- as.data.frame(sapply(tab_WeatherElast, function(x) (gsub(" ", "", x))))

# Merge CIs
tab_WeatherElast$L95_cereals <- paste0(paste0("(",tab_WeatherElast$L95_cereals,";"),
                                       paste0(tab_WeatherElast$U95_cereals,")"))
tab_WeatherElast$U95_cereals <- NULL
tab_WeatherElast <- dplyr::rename(tab_WeatherElast,"CI_cereals" = "L95_cereals")

tab_WeatherElast$L95_protein <- paste0(paste0("(",tab_WeatherElast$L95_protein,";"),
                                       paste0(tab_WeatherElast$U95_protein,")"))
tab_WeatherElast$U95_protein <- NULL
tab_WeatherElast <- dplyr::rename(tab_WeatherElast,"CI_protein" = "L95_protein")

tab_WeatherElast$L95_oilseed <- paste0(paste0("(",tab_WeatherElast$L95_oilseed,";"),
                                       paste0(tab_WeatherElast$U95_oilseed,")"))
tab_WeatherElast$U95_oilseed <- NULL
tab_WeatherElast <- dplyr::rename(tab_WeatherElast,"CI_oilseed" = "L95_oilseed")

tab_WeatherElast$L95_roots <- paste0(paste0("(",tab_WeatherElast$L95_roots,";"),
                                     paste0(tab_WeatherElast$U95_roots,")"))
tab_WeatherElast$U95_roots <- NULL
tab_WeatherElast <- dplyr::rename(tab_WeatherElast,"CI_roots" = "L95_roots")

tab_WeatherElast$L95_corn <- paste0(paste0("(",tab_WeatherElast$L95_corn,";"),
                                       paste0(tab_WeatherElast$U95_corn,")"))
tab_WeatherElast$U95_corn <- NULL
tab_WeatherElast <- dplyr::rename(tab_WeatherElast,"CI_corn" = "L95_corn")

tab_WeatherElast$L95_fert <- paste0(paste0("(",tab_WeatherElast$L95_fert,";"),
                                    paste0(tab_WeatherElast$U95_fert,")"))
tab_WeatherElast$U95_fert <- NULL
tab_WeatherElast <- dplyr::rename(tab_WeatherElast,"CI_fert" = "L95_fert")


#Insert a blank row to every other row
na.df <- data.frame("na" = NA, "Q_Cereals" = NA, "CI_cereals" = NA, "Q_Protein" = NA,
                    "CI_protein" = NA, "Q_Oilseed" = NA, "CI_oilseed" = NA, "Q_Roots" = NA,
                    "CI_roots" = NA, "Q_Corn" = NA, "CI_corn" = NA, "X_Fertilizer" = NA,
                    "CI_fert" = NA)
tab_WeatherElast <- do.call(rbind, apply(tab_WeatherElast, 1, function(x) {rbind(x, na.df)}))

## Shift CIs down one line
tab_WeatherElast$CI_cereals <- c(NA,tab_WeatherElast$CI_cereals[1:nrow(tab_WeatherElast)-1])
tab_WeatherElast$CI_protein <- c(NA,tab_WeatherElast$CI_protein[1:nrow(tab_WeatherElast)-1])
tab_WeatherElast$CI_oilseed <- c(NA,tab_WeatherElast$CI_oilseed[1:nrow(tab_WeatherElast)-1])
tab_WeatherElast$CI_roots <- c(NA,tab_WeatherElast$CI_roots[1:nrow(tab_WeatherElast)-1])
tab_WeatherElast$CI_corn <- c(NA,tab_WeatherElast$CI_corn[1:nrow(tab_WeatherElast)-1])
tab_WeatherElast$CI_fert <- c(NA,tab_WeatherElast$CI_fert[1:nrow(tab_WeatherElast)-1])

# Shifts CIs to the left
tab_WeatherElast$Q_Cereals[is.na(tab_WeatherElast$Q_Cereals)] <- tab_WeatherElast$CI_cereals[is.na(tab_WeatherElast$Q_Cereals)]
tab_WeatherElast$CI_cereals <- NULL

tab_WeatherElast$Q_Protein[is.na(tab_WeatherElast$Q_Protein)] <- tab_WeatherElast$CI_protein[is.na(tab_WeatherElast$Q_Protein)]
tab_WeatherElast$CI_protein <- NULL

tab_WeatherElast$Q_Oilseed[is.na(tab_WeatherElast$Q_Oilseed)] <- tab_WeatherElast$CI_oilseed[is.na(tab_WeatherElast$Q_Oilseed)]
tab_WeatherElast$CI_oilseed <- NULL

tab_WeatherElast$Q_Roots[is.na(tab_WeatherElast$Q_Roots)] <- tab_WeatherElast$CI_roots[is.na(tab_WeatherElast$Q_Roots)]
tab_WeatherElast$CI_roots <- NULL

tab_WeatherElast$Q_Corn[is.na(tab_WeatherElast$Q_Corn)] <- tab_WeatherElast$CI_corn[is.na(tab_WeatherElast$Q_Corn)]
tab_WeatherElast$CI_corn <- NULL

tab_WeatherElast$X_Fertilizer[is.na(tab_WeatherElast$X_Fertilizer)] <- tab_WeatherElast$CI_fert[is.na(tab_WeatherElast$X_Fertilizer)]
tab_WeatherElast$CI_fert <- NULL

# Insert observed and past weather indicators
tab_WeatherElast <- berryFunctions::insertRows(tab_WeatherElast, 1 , new = c("Observed weather","","","","","",""))
tab_WeatherElast <- berryFunctions::insertRows(tab_WeatherElast, nrow(WeatherElast_cereals)+2 , new = c("Past weather","","","","","",""))

# Print
if (cluster=="clustFarm") {
  write_xlsx(tab_WeatherElast,path="Tables/Table_3.xlsx")
} else if (cluster=="clustNuts2") {
  write_xlsx(tab_WeatherElast,path="Tables/Table_S4.4.xlsx")
}


