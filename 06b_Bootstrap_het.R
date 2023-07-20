#-----------------------------------------------------#
#                                                     #
# This program estimates the model with heterogeneity # 
# with bootstrapped standard errors.                  #
#                                                     #
# Note: The bootstrap program (Section 1: Bootstrap   #
# program) cannot be run because farm-level data are  #
# confidential. The results of the bootstrap  are     #
# stored in lines 1570 and 1585, allowing to          # 
# create the tables in Section 2 (Presentation of     #
# results).                                           # 
#                                                     #
#-----------------------------------------------------#

# load packages
library(dplyr)
library(writexl)
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
                              gdd_1to3_land, prec_1to3_land, gddHigh_1to3_land, dd_1to3_land,
                              gdd_4to10_land, prec_4to10_land, gddHigh_4to10_land, dd_4to10_land,
                              np_cereals_fm, np_protein_fm, np_oilseed_fm, np_roots_fm, np_corn_fm, nw_fert_fm, 
                              k_land_fm, k_labor_fm, k_capital_fm, trend_fm, trend2_fm, 
                              lsh_cereals_fm, lsh_protein_fm, lsh_oilseed_fm, lsh_roots_fm, 
                              gdd_1to3_fm, prec_1to3_fm, gddHigh_1to3_fm, dd_1to3_fm, 
                              gdd_4to10_fm, prec_4to10_fm, gddHigh_4to10_fm, dd_4to10_fm,
                              gdd_1to3_land_fm, prec_1to3_land_fm, gddHigh_1to3_land_fm, dd_1to3_land_fm,
                              gdd_4to10_land_fm, prec_4to10_land_fm, gddHigh_4to10_land_fm, dd_4to10_land_fm,
                              qq_cereals, qq_protein, qq_oilseed, qq_roots, qq_corn, nx_fert, x_fert, x_otherinp,
                              gdd_obs, prec_obs, gddHigh_obs, dd_obs,
                              gdd_obs_1to3, prec_obs_1to3, gddHigh_obs_1to3, dd_obs_1to3,
                              gdd_obs_4to10, prec_obs_4to10, gddHigh_obs_4to10, dd_obs_4to10,
                              gdd_obs_land, prec_obs_land, gddHigh_obs_land, dd_obs_land, 
                              gdd_obs_fm, prec_obs_fm, gddHigh_obs_fm, dd_obs_fm,
                              gdd_obs_1to3_fm, prec_obs_1to3_fm, gddHigh_obs_1to3_fm, dd_obs_1to3_fm,
                              gdd_obs_4to10_fm, prec_obs_4to10_fm, gddHigh_obs_4to10_fm, dd_obs_4to10_fm,
                              gdd_obs_land_fm, prec_obs_land_fm, gddHigh_obs_land_fm, dd_obs_land_fm,
                              k_land)


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
    
    tryCatch( { 
      
  # -------------------- #
  # START OF ESTIMATIONS # 
  # -------------------- #
      
      library(dplyr)
      
  #calculate sample mean for elasticity evaluation
  dat_sm <- dplyr::summarise_all(data_boot, mean) 
 
  dat_smallsm <- data_boot %>% 
    dplyr::filter(k_land<median(k_land)) %>% 
    dplyr::summarise_all(mean)
      
  dat_largesm <- data_boot %>% 
    dplyr::filter(k_land>median(k_land)) %>% 
    dplyr::summarise_all(mean)
  
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
    
  #Define structural equations
    
    eqQQcereals <- qq_cereals ~ PHI_cereals + I(PHI_cereals*np_cereals) + I(PHI_cereals*np_protein) + I(PHI_cereals*np_oilseed) + I(PHI_cereals*np_roots) + I(PHI_cereals*np_corn) + I(PHI_cereals*nw_fert) + 
      I(PHI_cereals*k_land) + I(PHI_cereals*k_labor) + I(PHI_cereals*k_capital) + I(PHI_cereals*trend) + I(PHI_cereals*trend2) +
      I(PHI_cereals*gdd_obs) + I(PHI_cereals*prec_obs) + I(PHI_cereals*gddHigh_obs) + I(PHI_cereals*dd_obs) + 
      I(PHI_cereals*gdd_1to3) + I(PHI_cereals*prec_1to3) + I(PHI_cereals*gddHigh_1to3) + I(PHI_cereals*dd_1to3) +
      I(PHI_cereals*gdd_4to10) + I(PHI_cereals*prec_4to10) + I(PHI_cereals*gddHigh_4to10) + I(PHI_cereals*dd_4to10) +
      I(PHI_cereals*gdd_obs_land) + I(PHI_cereals*prec_obs_land) + I(PHI_cereals*gddHigh_obs_land) + I(PHI_cereals*dd_obs_land) + 
      I(PHI_cereals*gdd_1to3_land) + I(PHI_cereals*prec_1to3_land) + I(PHI_cereals*gddHigh_1to3_land) + I(PHI_cereals*dd_1to3_land) + 
      I(PHI_cereals*gdd_4to10_land) + I(PHI_cereals*prec_4to10_land) + I(PHI_cereals*gddHigh_4to10_land) + I(PHI_cereals*dd_4to10_land) + 
      I(PHI_cereals*gdd_obs_1to3) + I(PHI_cereals*prec_obs_1to3) + I(PHI_cereals*gddHigh_obs_1to3) + I(PHI_cereals*dd_obs_1to3) + 
      I(PHI_cereals*gdd_obs_4to10) + I(PHI_cereals*prec_obs_4to10) + I(PHI_cereals*gddHigh_obs_4to10) + I(PHI_cereals*dd_obs_4to10) + 
      I(PHI_cereals*np_cereals_fm) + I(PHI_cereals*np_protein_fm) + I(PHI_cereals*np_oilseed_fm) + I(PHI_cereals*np_roots_fm) + I(PHI_cereals*np_corn_fm) + I(PHI_cereals*nw_fert_fm) +
      I(PHI_cereals*k_land_fm) + I(PHI_cereals*k_labor_fm) + I(PHI_cereals*k_capital_fm) + I(PHI_cereals*trend_fm) + I(PHI_cereals*trend2_fm) + 
      I(PHI_cereals*gdd_obs_fm) + I(PHI_cereals*prec_obs_fm) + I(PHI_cereals*gddHigh_obs_fm) + I(PHI_cereals*dd_obs_fm) + 
      I(PHI_cereals*gdd_1to3_fm) + I(PHI_cereals*prec_1to3_fm) + I(PHI_cereals*gddHigh_1to3_fm) + I(PHI_cereals*dd_1to3_fm) +
      I(PHI_cereals*gdd_4to10_fm) + I(PHI_cereals*prec_4to10_fm) + I(PHI_cereals*gddHigh_4to10_fm) + I(PHI_cereals*dd_4to10_fm) + 
      I(PHI_cereals*gdd_obs_land_fm) + I(PHI_cereals*prec_obs_land_fm) + I(PHI_cereals*gddHigh_obs_land_fm) + I(PHI_cereals*dd_obs_land_fm) + 
      I(PHI_cereals*gdd_1to3_land_fm) + I(PHI_cereals*prec_1to3_land_fm) + I(PHI_cereals*gddHigh_1to3_land_fm) + I(PHI_cereals*dd_1to3_land_fm) + 
      I(PHI_cereals*gdd_4to10_land_fm) + I(PHI_cereals*prec_4to10_land_fm) + I(PHI_cereals*gddHigh_4to10_land_fm) + I(PHI_cereals*dd_4to10_land_fm) + 
      I(PHI_cereals*gdd_obs_1to3_fm) + I(PHI_cereals*prec_obs_1to3_fm) + I(PHI_cereals*gddHigh_obs_1to3_fm) + I(PHI_cereals*dd_obs_1to3_fm) + 
      I(PHI_cereals*gdd_obs_4to10_fm) + I(PHI_cereals*prec_obs_4to10_fm) + I(PHI_cereals*gddHigh_obs_4to10_fm) + I(PHI_cereals*dd_obs_4to10_fm) + 
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
    
    # All farms
    dat_sm_prob <- dat_prob %>% 
      summarise_all(mean)
    
    # Small farms
    dat_smallsm_prob <- dat_prob %>% 
      filter(k_land<median(k_land)) %>% 
      summarise_all(mean)
    
    # Large farms
    dat_largesm_prob <- dat_prob %>% 
      filter(k_land>median(k_land)) %>% 
      summarise_all(mean)
  
  #set-up the loop over crops
  list_crops <- data.frame(crops=c("Cereals","Protein","Oilseed","Roots","Corn"),
                           probit=c("cereal.eq","protein.eq","oilseed.eq","roots.eq","corn.eq"),
                           abbrev=c("cer","prot","oil","roots","corn"))
  n_crops <- list_crops %>% count() %>% as.numeric() # number of crops
    
  # Create list where the results will be stored
  list_probit_act <- list(Cereals=list(PHI="", PHI_sm="", phi="", phi_sm="", probit_me="", linpred_sm="",
                                       PHI_sm_small="",phi_sm_small="", linpred_small ="",
                                       PHI_sm_large="",phi_sm_large="", linpred_large =""),
                          Protein=list(PHI="", PHI_sm="", phi="", phi_sm="", probit_me="", linpred_sm="",
                                       PHI_sm_small="",phi_sm_small="", linpred_small ="",
                                       PHI_sm_large="",phi_sm_large="", linpred_large =""),
                          Oilseed=list(PHI="", PHI_sm="", phi="", phi_sm="", probit_me="", linpred_sm="",
                                       PHI_sm_small="",phi_sm_small="", linpred_small ="",
                                       PHI_sm_large="",phi_sm_large="", linpred_large =""),
                          Roots=  list(PHI="", PHI_sm="", phi="", phi_sm="", probit_me="", linpred_sm="",
                                       PHI_sm_small="",phi_sm_small="", linpred_small ="",
                                       PHI_sm_large="",phi_sm_large="", linpred_large =""),
                          Corn=list(PHI="", PHI_sm="", phi="", phi_sm="", probit_me="", linpred_sm="",
                                       PHI_sm_small="",phi_sm_small="", linpred_small ="",
                                       PHI_sm_large="",phi_sm_large="", linpred_large =""))
  
  # Loop over all crops and estimate probits 
    
    for (i in 1:n_crops) { 
      
      # Probit estimation
      crop.eq <- get(list_crops[i,"probit"])
      iCrop <- glm(crop.eq, family = binomial(link = "probit"), 
                   data = data_boot)
  
      # At the farm-level
      iCrop_linpred <- predict(iCrop) 
      
      # At the sample mean (i.e., take crossproduct of coefficients and sample means)
      iCrop_linpred_sm <- crossprod(matrix(iCrop$coefficients),
                                    t(dat_sm_prob))
      list_probit_act[[i]][[6]] <- iCrop_linpred_sm
      
      iCrop_linpred_sm_small <- crossprod(matrix(iCrop$coefficients),
                                          t(dat_smallsm_prob))
      list_probit_act[[i]][[9]] <- iCrop_linpred_sm_small
      
      iCrop_linpred_sm_large <- crossprod(matrix(iCrop$coefficients),
                                          t(dat_largesm_prob))
      list_probit_act[[i]][[12]] <- iCrop_linpred_sm_large
      
      # PHI
      list_probit_act[[i]][[1]] <- pnorm(iCrop_linpred)
      list_probit_act[[i]][[2]] <- pnorm(iCrop_linpred_sm)
      list_probit_act[[i]][[7]] <- pnorm(iCrop_linpred_sm_small)
      list_probit_act[[i]][[10]] <- pnorm(iCrop_linpred_sm_large)
      
      # phi
      list_probit_act[[i]][[3]] <- dnorm(iCrop_linpred)
      list_probit_act[[i]][[4]] <- dnorm(iCrop_linpred_sm)
      
      list_probit_act[[i]][[8]] <- dnorm(iCrop_linpred_sm_small)
      list_probit_act[[i]][[11]] <- dnorm(iCrop_linpred_sm_large)
      
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
      
      # Probit marginal effects
      list_probit_act[[i]][[5]] <- dnorm(iCrop_linpred_sm)*coef(iCrop) #note: dnorm(iCrop_linpred_sm) is phi_crop_sm
  
    }
  
  # Restore the results estimated in the loop
    
  # PHIs at the farm-level (for estimation of second stage)
  data_boot$PHI_cereals <- as.numeric(do.call(c, list_probit_act[[1]][1]))
  data_boot$PHI_protein <- as.numeric(do.call(c, list_probit_act[[2]][1]))
  data_boot$PHI_oilseed <- as.numeric(do.call(c, list_probit_act[[3]][1]))
  data_boot$PHI_roots <-   as.numeric(do.call(c, list_probit_act[[4]][1]))
  data_boot$PHI_corn <- as.numeric(do.call(c, list_probit_act[[5]][1]))
  
  # phis at the farm-level (for estimation of second stage)
  data_boot$phi_cereals <- as.numeric(do.call(c, list_probit_act[[1]][3]))
  data_boot$phi_protein <- as.numeric(do.call(c, list_probit_act[[2]][3]))
  data_boot$phi_oilseed <- as.numeric(do.call(c, list_probit_act[[3]][3]))
  data_boot$phi_roots <-   as.numeric(do.call(c, list_probit_act[[4]][3]))
  data_boot$phi_corn <- as.numeric(do.call(c, list_probit_act[[5]][3]))
  
  # PHIs, phis, and linpreds at sample mean, small farms, and large farms (for elasticity evaluation)
  PHI_cereals_sm <- as.numeric(do.call(c, list_probit_act[[1]][2]))
  phi_cereals_sm <- as.numeric(do.call(c, list_probit_act[[1]][4]))
  iCereals_linpred_sm <- as.numeric(do.call(c, list_probit_act[[1]][6]))
  PHI_cereals_sm_small <- as.numeric(do.call(c, list_probit_act[[1]][7]))
  phi_cereals_sm_small <- as.numeric(do.call(c, list_probit_act[[1]][8]))
  iCereals_linpred_sm_small <- as.numeric(do.call(c, list_probit_act[[1]][9]))
  PHI_cereals_sm_large <- as.numeric(do.call(c, list_probit_act[[1]][10]))
  phi_cereals_sm_large <- as.numeric(do.call(c, list_probit_act[[1]][11]))
  iCereals_linpred_sm_large <- as.numeric(do.call(c, list_probit_act[[1]][12]))
  
  PHI_protein_sm <- as.numeric(do.call(c, list_probit_act[[2]][2]))
  phi_protein_sm <- as.numeric(do.call(c, list_probit_act[[2]][4]))
  iProtein_linpred_sm <- as.numeric(do.call(c, list_probit_act[[2]][6]))
  PHI_protein_sm_small <- as.numeric(do.call(c, list_probit_act[[2]][7]))
  phi_protein_sm_small <- as.numeric(do.call(c, list_probit_act[[2]][8]))
  iProtein_linpred_sm_small <- as.numeric(do.call(c, list_probit_act[[2]][9]))
  PHI_protein_sm_large <- as.numeric(do.call(c, list_probit_act[[2]][10]))
  phi_protein_sm_large <- as.numeric(do.call(c, list_probit_act[[2]][11]))
  iProtein_linpred_sm_large <- as.numeric(do.call(c, list_probit_act[[2]][12]))
  
  PHI_oilseed_sm <- as.numeric(do.call(c, list_probit_act[[3]][2]))
  phi_oilseed_sm <- as.numeric(do.call(c, list_probit_act[[3]][5]))
  iOilseed_linpred_sm <- as.numeric(do.call(c, list_probit_act[[3]][6]))
  PHI_oilseed_sm_small <- as.numeric(do.call(c, list_probit_act[[3]][7]))
  phi_oilseed_sm_small <- as.numeric(do.call(c, list_probit_act[[3]][8]))
  iOilseed_linpred_sm_small <- as.numeric(do.call(c, list_probit_act[[3]][9]))
  PHI_oilseed_sm_large <- as.numeric(do.call(c, list_probit_act[[3]][10]))
  phi_oilseed_sm_large <- as.numeric(do.call(c, list_probit_act[[3]][11]))
  iOilseed_linpred_sm_large <- as.numeric(do.call(c, list_probit_act[[3]][12]))
  
  PHI_roots_sm <- as.numeric(do.call(c, list_probit_act[[4]][2]))
  phi_roots_sm <- as.numeric(do.call(c, list_probit_act[[4]][5]))
  iRoots_linpred_sm <- as.numeric(do.call(c, list_probit_act[[4]][6]))
  PHI_roots_sm_small <- as.numeric(do.call(c, list_probit_act[[4]][7]))
  phi_roots_sm_small <- as.numeric(do.call(c, list_probit_act[[4]][8]))
  iRoots_linpred_sm_small <- as.numeric(do.call(c, list_probit_act[[4]][9]))
  PHI_roots_sm_large <- as.numeric(do.call(c, list_probit_act[[4]][10]))
  phi_roots_sm_large <- as.numeric(do.call(c, list_probit_act[[4]][11]))
  iRoots_linpred_sm_large <- as.numeric(do.call(c, list_probit_act[[4]][12]))
  
  PHI_corn_sm <- as.numeric(do.call(c, list_probit_act[[5]][2]))
  phi_corn_sm <- as.numeric(do.call(c, list_probit_act[[5]][5]))
  iCorn_linpred_sm <- as.numeric(do.call(c, list_probit_act[[5]][6]))
  PHI_corn_sm_small <- as.numeric(do.call(c, list_probit_act[[5]][7]))
  phi_corn_sm_small <- as.numeric(do.call(c, list_probit_act[[5]][8]))
  iCorn_linpred_sm_small <- as.numeric(do.call(c, list_probit_act[[5]][9]))
  PHI_corn_sm_large <- as.numeric(do.call(c, list_probit_act[[5]][10]))
  phi_corn_sm_large <- as.numeric(do.call(c, list_probit_act[[5]][11]))
  iCorn_linpred_sm_large <- as.numeric(do.call(c, list_probit_act[[5]][12]))
  
  # Probit marginal effects (for results presentation)
  iCereals_me <- as.numeric(do.call(c, list_probit_act[[1]][5]))
  iProtein_me <- as.numeric(do.call(c, list_probit_act[[2]][5]))
  iOilseed_me <- as.numeric(do.call(c, list_probit_act[[3]][5]))
  iRoots_me   <- as.numeric(do.call(c, list_probit_act[[4]][5]))
  iCorn_me    <- as.numeric(do.call(c, list_probit_act[[5]][5]))
  
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
  
  # small farms 
  
  qx_pred$k_land <- data_boot$k_land
  
  qq_cereals_pred_b_small <- mean(qx_pred$QQcereals.pred[qx_pred$k_land<median(qx_pred$k_land)]) - 
    model_linear$coefficients["QQcereals_phi_cereals"]*mean(data_boot$phi_cereals[data_boot$k_land<median(data_boot$k_land)]) 
  qq_cereals_pred_small <- (qq_cereals_pred_b_small) / mean(data_boot$PHI_cereals[data_boot$k_land<median(data_boot$k_land)])
  
  qq_protein_pred_b_small <- mean(qx_pred$QQprotein.pred[qx_pred$k_land<median(qx_pred$k_land)]) - 
    model_linear$coefficients["QQprotein_phi_protein"]*mean(data_boot$phi_protein[data_boot$k_land<median(data_boot$k_land)]) 
  qq_protein_pred_small <- (qq_protein_pred_b_small) / mean(data_boot$PHI_protein[data_boot$k_land<median(data_boot$k_land)])
  
  qq_oilseed_pred_b_small <- mean(qx_pred$QQoilseed.pred[qx_pred$k_land<median(qx_pred$k_land)]) - 
    model_linear$coefficients["QQoilseed_phi_oilseed"]*mean(data_boot$phi_oilseed[data_boot$k_land<median(data_boot$k_land)]) 
  qq_oilseed_pred_small <- (qq_oilseed_pred_b_small) / mean(data_boot$PHI_oilseed[data_boot$k_land<median(data_boot$k_land)])
  
  qq_roots_pred_b_small <- mean(qx_pred$QQroots.pred[qx_pred$k_land<median(qx_pred$k_land)]) - 
    model_linear$coefficients["QQroots_phi_roots"]*mean(data_boot$phi_roots[data_boot$k_land<median(data_boot$k_land)]) 
  qq_roots_pred_small <- (qq_roots_pred_b_small) / mean(data_boot$PHI_roots[data_boot$k_land<median(data_boot$k_land)])
  
  qq_corn_pred_b_small <- mean(qx_pred$QQcorn.pred[qx_pred$k_land<median(qx_pred$k_land)]) - 
    model_linear$coefficients["QQcorn_phi_corn"]*mean(data_boot$phi_corn[data_boot$k_land<median(data_boot$k_land)]) 
  qq_corn_pred_small <- (qq_corn_pred_b_small) / mean(data_boot$PHI_corn[data_boot$k_land<median(data_boot$k_land)])
  
  nx_fert_pred_small <- mean(qx_pred$NXfert.pred[qx_pred$k_land<median(qx_pred$k_land)])
  
  # large farms
  
  qq_cereals_pred_b_large <- mean(qx_pred$QQcereals.pred[qx_pred$k_land>median(qx_pred$k_land)]) - 
    model_linear$coefficients["QQcereals_phi_cereals"]*mean(data_boot$phi_cereals[data_boot$k_land>median(data_boot$k_land)]) 
  qq_cereals_pred_large <- (qq_cereals_pred_b_large) / mean(data_boot$PHI_cereals[data_boot$k_land>median(data_boot$k_land)])
  
  qq_protein_pred_b_large <- mean(qx_pred$QQprotein.pred[qx_pred$k_land>median(qx_pred$k_land)]) - 
    model_linear$coefficients["QQprotein_phi_protein"]*mean(data_boot$phi_protein[data_boot$k_land>median(data_boot$k_land)]) 
  qq_protein_pred_large <- (qq_protein_pred_b_large) / mean(data_boot$PHI_protein[data_boot$k_land>median(data_boot$k_land)])
  
  qq_oilseed_pred_b_large <- mean(qx_pred$QQoilseed.pred[qx_pred$k_land>median(qx_pred$k_land)]) - 
    model_linear$coefficients["QQoilseed_phi_oilseed"]*mean(data_boot$phi_oilseed[data_boot$k_land>median(data_boot$k_land)]) 
  qq_oilseed_pred_large <- (qq_oilseed_pred_b_large) / mean(data_boot$PHI_oilseed[data_boot$k_land>median(data_boot$k_land)])
  
  qq_roots_pred_b_large <- mean(qx_pred$QQroots.pred[qx_pred$k_land>median(qx_pred$k_land)]) - 
    model_linear$coefficients["QQroots_phi_roots"]*mean(data_boot$phi_roots[data_boot$k_land>median(data_boot$k_land)]) 
  qq_roots_pred_large <- (qq_roots_pred_b_large) / mean(data_boot$PHI_roots[data_boot$k_land>median(data_boot$k_land)])
  
  qq_corn_pred_b_large <- mean(qx_pred$QQcorn.pred[qx_pred$k_land>median(qx_pred$k_land)]) - 
    model_linear$coefficients["QQcorn_phi_corn"]*mean(data_boot$phi_corn[data_boot$k_land>median(data_boot$k_land)]) 
  qq_corn_pred_large <- (qq_corn_pred_b_large) / mean(data_boot$PHI_corn[data_boot$k_land>median(data_boot$k_land)])
  
  nx_fert_pred_large <- mean(qx_pred$NXfert.pred[qx_pred$k_land>median(qx_pred$k_land)])
  
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
    
    elast_weather_p25 <- data.frame("Variables"=c("GDD_obs", "PREC_obs", "GDDHigh_obs", "DD_obs", "GDD_past", "PREC_past", "GDDHigh_past", "DD_past"),
                                    "Q Cereals"=rep(NA,8),
                                    "Q Protein"=rep(NA,8),
                                    "Q Oilseed"=rep(NA,8),
                                    "Q Roots"=rep(NA,8),
                                    "Q Corn"=rep(NA,8),
                                    "X Fertilizer"=rep(NA,8))  
    
    elast_weather_p75 <- data.frame("Variables"=c("GDD_obs", "PREC_obs", "GDDHigh_obs", "DD_obs", "GDD_past", "PREC_past", "GDDHigh_past", "DD_past"),
                                    "Q Cereals"=rep(NA,8),
                                    "Q Protein"=rep(NA,8),
                                    "Q Oilseed"=rep(NA,8),
                                    "Q Roots"=rep(NA,8),
                                    "Q Corn"=rep(NA,8),
                                    "X Fertilizer"=rep(NA,8)) 
    
    # Weather vars
    list_weather <- data.frame(weather_obs=  c("gdd_obs","prec_obs","gddHigh_obs","dd_obs"),
                               weather_obs_land = c("gdd_obs_land","prec_obs_land","gddHigh_obs_land","dd_obs_land"),
                               weather_1to3= c("gdd_1to3","prec_1to3","gddHigh_1to3","dd_1to3"),
                               weather_1to3_land = c("gdd_1to3_land","prec_1to3_land","gddHigh_1to3_land","dd_1to3_land"),
                               weather_4to10=c("gdd_4to10","prec_4to10","gddHigh_4to10","dd_4to10"),
                               weather_4to10_land = c("gdd_4to10_land","prec_4to10_land","gddHigh_4to10_land","dd_4to10_land"),
                               weather_obs_1to3 = c("gdd_obs_1to3","prec_obs_1to3","gddHigh_obs_1to3","dd_obs_1to3"),
                               weather_obs_4to10 = c("gdd_obs_4to10","prec_obs_4to10","gddHigh_obs_4to10","dd_obs_4to10"),
                               weather_obs_endwith=  c("* gdd_obs)","* prec_obs)","* gddHigh_obs)","* dd_obs)"),
                               weather_obs_land_endwith = c("* gdd_obs_land)","* prec_obs_land)","* gddHigh_obs_land)","* dd_obs_land)"),
                               weather_1to3_endwith= c("* gdd_1to3)","* prec_1to3)","* gddHigh_1to3)","* dd_1to3)"),
                               weather_1to3_land_endwith = c("* gdd_1to3_land)","* prec_1to3_land)","* gddHigh_1to3_land)","* dd_1to3_land)"),
                               weather_4to10_endwith=c("* gdd_4to10)","* prec_4to10)","* gddHigh_4to10)","* dd_4to10)"),
                               weather_4to10_land_endwith = c("* gdd_4to10_land)","* prec_4to10_land)","* gddHigh_4to10_land)","* dd_4to10_land)"),
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
      
      PHI_p25 <- get(paste0("PHI_",list_qcrops[i,"crops"],"_sm_small"))
      phi_p25 <- get(paste0("phi_",list_qcrops[i,"crops"],"_sm_small"))
      linpred_p25 <- get(paste0(list_qcrops[i,"probit"],"_linpred_sm_small"))
      pred_p25 <- get(paste0("qq_",list_qcrops[i,"crops"],"_pred_small"))
      mean_quantity_p25 <- dat_smallsm %>% select(paste0("qq_",list_qcrops[i,"crops"]))
      
      PHI_p75 <- get(paste0("PHI_",list_qcrops[i,"crops"],"_sm_large"))
      phi_p75 <- get(paste0("phi_",list_qcrops[i,"crops"],"_sm_large"))
      linpred_p75 <- get(paste0(list_qcrops[i,"probit"],"_linpred_sm_large"))
      pred_p75 <- get(paste0("qq_",list_qcrops[i,"crops"],"_pred_large")) 
      mean_quantity_p75 <- dat_largesm %>% select(paste0("qq_",list_qcrops[i,"crops"]))
      
      # Loop over all weather variables
      for (j in 1:n_weather) {
        
        # Weather-land terms
        coef_probit_weather_land_1to3 <- coefs_probit_crop %>% select(list_weather[j,"weather_1to3_land"])
        coef_probit_weather_land_4to10 <-  coefs_probit_crop %>% select(list_weather[j,"weather_4to10_land"])
        
        coef_struct_weather_land_obs <-  coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_obs_land_endwith"]))
        coef_struct_weather_land_1to3 <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_1to3_land_endwith"]))
        coef_struct_weather_land_4to10 <-  coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_4to10_land_endwith"]))
        
        # Weather-interaction terms
        coef_struct_weather_obs_1to3 <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_obs_1to3_endwith"]))
        coef_struct_weather_obs_4to10 <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_obs_4to10_endwith"]))
        
        # Mean weather for derivatives 
        mean_weather_obs <- dat_sm %>% select(list_weather[j,"weather_obs"])
        mean_weather_1to3 <- dat_sm %>% select(list_weather[j,"weather_1to3"])
        mean_weather_4to10 <- dat_sm %>% select(list_weather[j,"weather_4to10"])
        
        # Mean land for derivatives (note: can be changed to percentiles to compute heterogeneous elasticities)
        mean_land <- dat_sm %>% select("k_land")
        p25_land <- dat_smallsm %>% select("k_land")
        p75_land <- dat_largesm %>% select("k_land")
        
        # ------------------ #
        # Land = sample mean #
        # ------------------ #
        
        # Observed weather (note: observed weather only affects structural equations, not crop selection; also note that PHI_* is treated like a constant term )
        coef_struct_weather <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_obs_endwith"]))
        elast_weather[j,i+1] <- PHI * (coef_struct_weather + coef_struct_weather_obs_1to3*mean_weather_1to3 + coef_struct_weather_obs_4to10*mean_weather_4to10 + coef_struct_weather_land_obs*mean_land)  * (1/mean_quantity)
        
        # Weather_1to3
        coef_probit_weather <- coefs_probit_crop %>% select(list_weather[j,"weather_1to3"])
        coef_struct_weather <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_1to3_endwith"]))
        
        elast_weather_1to3 <- (PHI * (coef_struct_weather + coef_struct_weather_obs_1to3*mean_weather_obs + coef_struct_weather_land_1to3*mean_land) + phi * pred * 
                                 (coef_probit_weather + coef_probit_weather_land_1to3*mean_land) -
                                 coef_phi_struct * linpred * (coef_probit_weather + coef_probit_weather_land_1to3*mean_land) * phi) * (1 / mean_quantity )
        
        # Weather_4to10
        coef_probit_weather <- coefs_probit_crop %>% select(list_weather[j,"weather_4to10"])
        coef_struct_weather <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_4to10_endwith"]))
        
        elast_weather_4to10 <- (PHI * (coef_struct_weather + coef_struct_weather_obs_4to10*mean_weather_obs + coef_struct_weather_land_4to10*mean_land) + phi * pred * 
                                  (coef_probit_weather + coef_probit_weather_land_4to10*mean_land) -
                                  coef_phi_struct * linpred * (coef_probit_weather + coef_probit_weather_land_4to10*mean_land) * phi) * (1 / mean_quantity )
        
        # Weather_past
        elast_weather[j+4,i+1] <- elast_weather_1to3 + elast_weather_4to10
        
        # ---------- #
        # Land = p25 #
        # ---------- #
        
        # Old codes: PHI is the same as at the sample mean
        
        # Observed weather (note: observed weather only affects structural equations, not crop selection; also note that PHI_* is treated like a constant term )
        coef_struct_weather <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_obs_endwith"]))
        elast_weather_p25[j,i+1] <- PHI * (coef_struct_weather + coef_struct_weather_obs_1to3*mean_weather_1to3 + coef_struct_weather_obs_4to10*mean_weather_4to10 + coef_struct_weather_land_obs*p25_land)  * (1/mean_quantity)
        
        # Weather_1to3
        coef_probit_weather <- coefs_probit_crop %>% select(list_weather[j,"weather_1to3"])
        coef_struct_weather <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_1to3_endwith"]))
        
        elast_weather_1to3 <- (PHI * (coef_struct_weather + coef_struct_weather_obs_1to3*mean_weather_obs + coef_struct_weather_land_1to3*p25_land) + phi * pred * 
                                 (coef_probit_weather + coef_probit_weather_land_1to3*p25_land) -
                                 coef_phi_struct * linpred * (coef_probit_weather + coef_probit_weather_land_1to3*p25_land) * phi) * (1 / mean_quantity )
        
        # Weather_4to10
        coef_probit_weather <- coefs_probit_crop %>% select(list_weather[j,"weather_4to10"])
        coef_struct_weather <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_4to10_endwith"]))
        
        elast_weather_4to10 <- (PHI * (coef_struct_weather + coef_struct_weather_obs_4to10*mean_weather_obs + coef_struct_weather_land_4to10*p25_land) + phi * pred * 
                                  (coef_probit_weather + coef_probit_weather_land_4to10*p25_land) -
                                  coef_phi_struct * linpred * (coef_probit_weather + coef_probit_weather_land_4to10*p25_land) * phi) * (1 / mean_quantity )
        
        # Weather_past
        elast_weather_p25[j+4,i+1] <- elast_weather_1to3 + elast_weather_4to10
        
        # New codes: PHI is specific for small farms
        
        # Observed weather (note: observed weather only affects structural equations, not crop selection; also note that PHI_* is treated like a constant term )
        coef_struct_weather <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_obs_endwith"]))
        elast_weather_p25[j,i+1] <- PHI_p25 * (coef_struct_weather + coef_struct_weather_obs_1to3*mean_weather_1to3 + coef_struct_weather_obs_4to10*mean_weather_4to10 + coef_struct_weather_land_obs*p25_land)  * (1/mean_quantity_p25)
        
        # Weather_1to3
        coef_probit_weather <- coefs_probit_crop %>% select(list_weather[j,"weather_1to3"])
        coef_struct_weather <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_1to3_endwith"]))
        
        elast_weather_1to3 <- (PHI_p25 * (coef_struct_weather + coef_struct_weather_obs_1to3*mean_weather_obs + coef_struct_weather_land_1to3*p25_land) + phi_p25 * pred_p25 * 
                                 (coef_probit_weather + coef_probit_weather_land_1to3*p25_land) -
                                 coef_phi_struct * linpred_p25 * (coef_probit_weather + coef_probit_weather_land_1to3*p25_land) * phi_p25) * (1 / mean_quantity_p25 )
        
        # Weather_4to10
        coef_probit_weather <- coefs_probit_crop %>% select(list_weather[j,"weather_4to10"])
        coef_struct_weather <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_4to10_endwith"]))
        
        elast_weather_4to10 <- (PHI_p25 * (coef_struct_weather + coef_struct_weather_obs_4to10*mean_weather_obs + coef_struct_weather_land_4to10*p25_land) + phi_p25 * pred_p25 * 
                                  (coef_probit_weather + coef_probit_weather_land_4to10*p25_land) -
                                  coef_phi_struct * linpred_p25 * (coef_probit_weather + coef_probit_weather_land_4to10*p25_land) * phi_p25) * (1 / mean_quantity_p25 )
        
        # Weather_past
        elast_weather_p25[j+4,i+1] <- elast_weather_1to3 + elast_weather_4to10
        
        # ---------- #
        # Land = p75 #
        # ---------- #
        
        # Observed weather (note: observed weather only affects structural equations, not crop selection; also note that PHI_* is treated like a constant term )
        coef_struct_weather <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_obs_endwith"]))
        elast_weather_p75[j,i+1] <- PHI_p75 * (coef_struct_weather + coef_struct_weather_obs_1to3*mean_weather_1to3 + coef_struct_weather_obs_4to10*mean_weather_4to10 + coef_struct_weather_land_obs*p75_land)  * (1/mean_quantity)
        
        # Weather_1to3
        coef_probit_weather <- coefs_probit_crop %>% select(list_weather[j,"weather_1to3"])
        coef_struct_weather <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_1to3_endwith"]))
        
        elast_weather_1to3 <- (PHI_p75 * (coef_struct_weather + coef_struct_weather_obs_1to3*mean_weather_obs + coef_struct_weather_land_1to3*p75_land) + phi_p75 * pred_p75 * 
                                 (coef_probit_weather + coef_probit_weather_land_1to3*p75_land) -
                                 coef_phi_struct * linpred_p75 * (coef_probit_weather + coef_probit_weather_land_1to3*p75_land) * phi_p75) * (1 / mean_quantity_p75 )
        
        # Weather_4to10
        coef_probit_weather <- coefs_probit_crop %>% select(list_weather[j,"weather_4to10"])
        coef_struct_weather <- coefs_struct_crop %>% select(ends_with(list_weather[j,"weather_4to10_endwith"]))
        
        elast_weather_4to10 <- (PHI_p75 * (coef_struct_weather + coef_struct_weather_obs_4to10*mean_weather_obs + coef_struct_weather_land_4to10*p75_land) + phi_p75 * pred_p75 * 
                                  (coef_probit_weather + coef_probit_weather_land_4to10*p75_land) -
                                  coef_phi_struct * linpred_p75 * (coef_probit_weather + coef_probit_weather_land_4to10*p75_land) * phi_p75) * (1 / mean_quantity_p75 )
        
        # Weather_past
        elast_weather_p75[j+4,i+1] <- elast_weather_1to3 + elast_weather_4to10
        
      }
      
    }
    
    # ----------------------------------- #
    # Weather elasticities for fertilizer #
    # ----------------------------------- #      
    
    # Weather vars
    list_weather <- data.frame(weather_obs=  c("gdd_obs","prec_obs","gddHigh_obs","dd_obs"),
                               weather_obs_land = c("gdd_obs_land","prec_obs_land","gddHigh_obs_land","dd_obs_land"),
                               weather_1to3= c("gdd_1to3","prec_1to3","gddHigh_1to3","dd_1to3"),
                               weather_1to3_land = c("gdd_1to3_land","prec_1to3_land","gddHigh_1to3_land","dd_1to3_land"),
                               weather_4to10=c("gdd_4to10","prec_4to10","gddHigh_4to10","dd_4to10"),
                               weather_4to10_land = c("gdd_4to10_land","prec_4to10_land","gddHigh_4to10_land","dd_4to10_land"),
                               weather_obs_1to3 = c("gdd_obs_1to3","prec_obs_1to3","gddHigh_obs_1to3","dd_obs_1to3"),
                               weather_obs_4to10 = c("gdd_obs_4to10","prec_obs_4to10","gddHigh_obs_4to10","dd_obs_4to10"),
                               weather_obs_endwith=  c("_gdd_obs","_prec_obs","_gddHigh_obs","_dd_obs"),
                               weather_obs_land_endwith = c("_gdd_obs_land","_prec_obs_land","_gddHigh_obs_land","_dd_obs_land"),
                               weather_1to3_endwith= c("_gdd_1to3","_prec_1to3","_gddHigh_1to3","_dd_1to3"),
                               weather_1to3_land_endwith = c("_gdd_1to3_land","_prec_1to3_land","_gddHigh_1to3_land","_dd_1to3_land"),
                               weather_4to10_endwith=c("_gdd_4to10","_prec_4to10","_gddHigh_4to10","_dd_4to10"),
                               weather_4to10_land_endwith = c("_gdd_4to10_land","_prec_4to10_land","_gddHigh_4to10_land","_dd_4to10_land"),
                               weather_obs_1to3_endwith = c("_gdd_obs_1to3","_prec_obs_1to3","_gddHigh_obs_1to3","_dd_obs_1to3"),
                               weather_obs_4to10_endwith = c("_gdd_obs_4to10","_prec_obs_4to10","_gddHigh_obs_4to10","_dd_obs_4to10"))
    
    
    n_weather <- list_weather %>% count() %>% as.numeric() 
    
    coefs_fert <- dplyr::as_tibble(t(coef_linear)) %>% select(starts_with("NXfert")) 
    
    
    # Loop over all weather variables
    for (j in 1:n_weather) {
      
      # Weather-land terms
      coef_fert_weather_land_obs <-  coefs_fert %>% select(ends_with(list_weather[j,"weather_obs_land_endwith"]))
      coef_fert_weather_land_1to3 <- coefs_fert %>% select(ends_with(list_weather[j,"weather_1to3_land_endwith"]))
      coef_fert_weather_land_4to10 <-  coefs_fert %>% select(ends_with(list_weather[j,"weather_4to10_land_endwith"]))
      
      # Weather-interaction terms
      coef_fert_weather_obs_1to3 <- coefs_fert %>% select(ends_with(list_weather[j,"weather_obs_1to3_endwith"]))
      coef_fert_weather_obs_4to10 <- coefs_fert %>% select(ends_with(list_weather[j,"weather_obs_4to10_endwith"]))
      
      # Mean land for derivatives (note: can be changed to percentiles to compute heterogeneous elasticities)
      mean_land <- dat_sm %>% select("k_land")
      p25_land <- dat_smallsm %>% select("k_land")
      p75_land <- dat_largesm %>% select("k_land")
      
      # Mean weather for derivatives 
      mean_weather_obs <- dat_sm %>% select(list_weather[j,"weather_obs"])
      mean_weather_1to3 <- dat_sm %>% select(list_weather[j,"weather_1to3"])
      mean_weather_4to10 <- dat_sm %>% select(list_weather[j,"weather_4to10"])
      
      # ------------------ #
      # land = sample mean #
      # ------------------ #
      
      # Observed weather (note: observed weather only affects structural equations, not crop selection; also note that PHI_* is treated like a constant term )
      coef_fert_weather <- coefs_fert %>% select(ends_with(list_weather[j,"weather_obs_endwith"]))
      elast_weather[j,7] <- (coef_fert_weather + coef_fert_weather_obs_1to3*mean_weather_1to3 + coef_fert_weather_obs_4to10*mean_weather_4to10 + coef_fert_weather_land_obs*mean_land) * (1/dat_sm$nx_fert)
      
      # Weather_1to3
      coef_fert_weather <- coefs_fert %>% select(ends_with(list_weather[j,"weather_1to3_endwith"]))
      elast_weather_1to3 <- (coef_fert_weather + coef_fert_weather_obs_1to3*mean_weather_obs + coef_fert_weather_land_1to3*mean_land) * (1/dat_sm$nx_fert)
      
      # Weather_4to10
      coef_fert_weather <- coefs_fert %>% select(ends_with(list_weather[j,"weather_4to10_endwith"]))
      elast_weather_4to10 <- (coef_fert_weather + coef_fert_weather_obs_4to10*mean_weather_obs + coef_fert_weather_land_4to10*mean_land) * (1/dat_sm$nx_fert)
      
      # Weather_past
      elast_weather[j+4,7] <- elast_weather_1to3 + elast_weather_4to10
      
      # ---------- #
      # land = p25 #
      # ---------- #
      
      # Observed weather (note: observed weather only affects structural equations, not crop selection; also note that PHI_* is treated like a constant term )
      coef_fert_weather <- coefs_fert %>% select(ends_with(list_weather[j,"weather_obs_endwith"]))
      elast_weather_p25[j,7] <- (coef_fert_weather + coef_fert_weather_obs_1to3*mean_weather_1to3 + coef_fert_weather_obs_4to10*mean_weather_4to10 + coef_fert_weather_land_obs*p25_land) * (1/dat_smallsm$nx_fert)
      
      # Weather_1to3
      coef_fert_weather <- coefs_fert %>% select(ends_with(list_weather[j,"weather_1to3_endwith"]))
      elast_weather_1to3 <- (coef_fert_weather + coef_fert_weather_obs_1to3*mean_weather_obs + coef_fert_weather_land_1to3*p25_land) * (1/dat_smallsm$nx_fert)
      
      # Weather_4to10
      coef_fert_weather <- coefs_fert %>% select(ends_with(list_weather[j,"weather_4to10_endwith"]))
      elast_weather_4to10 <- (coef_fert_weather + coef_fert_weather_obs_4to10*mean_weather_obs + coef_fert_weather_land_4to10*p25_land) * (1/dat_smallsm$nx_fert)
      
      # Weather_past
      elast_weather_p25[j+4,7] <- elast_weather_1to3 + elast_weather_4to10
      
      # ---------- #
      # land = p75 #
      # ---------- #
      
      # Observed weather (note: observed weather only affects structural equations, not crop selection; also note that PHI_* is treated like a constant term )
      coef_fert_weather <- coefs_fert %>% select(ends_with(list_weather[j,"weather_obs_endwith"]))
      elast_weather_p75[j,7] <- (coef_fert_weather + coef_fert_weather_obs_1to3*mean_weather_1to3 + coef_fert_weather_obs_4to10*mean_weather_4to10 + coef_fert_weather_land_obs*p75_land) * (1/dat_largesm$nx_fert)
      
      # Weather_1to3
      coef_fert_weather <- coefs_fert %>% select(ends_with(list_weather[j,"weather_1to3_endwith"]))
      elast_weather_1to3 <- (coef_fert_weather + coef_fert_weather_obs_1to3*mean_weather_obs + coef_fert_weather_land_1to3*p75_land) * (1/dat_largesm$nx_fert)
      
      # Weather_4to10
      coef_fert_weather <- coefs_fert %>% select(ends_with(list_weather[j,"weather_4to10_endwith"]))
      elast_weather_4to10 <- (coef_fert_weather + coef_fert_weather_obs_4to10*mean_weather_obs + coef_fert_weather_land_4to10*p75_land) * (1/dat_largesm$nx_fert)
      
      # Weather_past
      elast_weather_p75[j+4,7] <- elast_weather_1to3 + elast_weather_4to10
      
    } 
    
    # --------------------------------------- #
    # Store weather elasticities in bootstrap #
    # ----------------------------------------#
    
      # Rescale results to obtain the units days, mm, days, days
      elast_weather[,2:7] <- elast_weather[,2:7]*100  
      elast_weather_p25[,2:7] <- elast_weather_p25[,2:7]*100  
      elast_weather_p75[,2:7] <- elast_weather_p75[,2:7]*100  
      
      # WEATHER AT THE SAMPLE MEAN
      
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
              
        # WEATHER AT P25
              
              #qcer
              
              # obs
              sel_qcer_gdd_obs_p25     <- c('sel_qcer_gdd_obs_p25'      = as.numeric(elast_weather_p25[1,2]))
              sel_qcer_prec_obs_p25    <- c('sel_qcer_prec_obs_p25'     = as.numeric(elast_weather_p25[2,2]))
              sel_qcer_gddHigh_obs_p25 <- c('sel_qcer_gddHigh_obs_p25'  = as.numeric(elast_weather_p25[3,2]))
              sel_qcer_dd_obs_p25      <- c('sel_qcer_dd_obs_p25'       = as.numeric(elast_weather_p25[4,2]))
              
              #past
              sel_qcer_gdd_past_p25 <- c('sel_qcer_gdd_past_p25'          = as.numeric(elast_weather_p25[5,2]))
              sel_qcer_prec_past_p25  <- c('sel_qcer_prec_past_p25'       = as.numeric(elast_weather_p25[6,2]))
              sel_qcer_gddHigh_past_p25  <- c('sel_qcer_gddHigh_past_p25' = as.numeric(elast_weather_p25[7,2]))
              sel_qcer_dd_past_p25  <- c('sel_qcer_dd_past_p25'           = as.numeric(elast_weather_p25[8,2]))
              
              #qprot
              
              # obs
              sel_qprot_gdd_obs_p25     <- c('sel_qprot_gdd_obs_p25'      = as.numeric(elast_weather_p25[1,3]))
              sel_qprot_prec_obs_p25    <- c('sel_qprot_prec_obs_p25'     = as.numeric(elast_weather_p25[2,3]))
              sel_qprot_gddHigh_obs_p25 <- c('sel_qprot_gddHigh_obs_p25'  = as.numeric(elast_weather_p25[3,3]))
              sel_qprot_dd_obs_p25      <- c('sel_qprot_dd_obs_p25'       = as.numeric(elast_weather_p25[4,3]))
              
              #past
              sel_qprot_gdd_past_p25 <- c('sel_qprot_gdd_past_p25'          = as.numeric(elast_weather_p25[5,3]))
              sel_qprot_prec_past_p25  <- c('sel_qprot_prec_past_p25'       = as.numeric(elast_weather_p25[6,3]))
              sel_qprot_gddHigh_past_p25  <- c('sel_qprot_gddHigh_past_p25' = as.numeric(elast_weather_p25[7,3]))
              sel_qprot_dd_past_p25  <- c('sel_qprot_dd_past_p25'           = as.numeric(elast_weather_p25[8,3]))
              
              #qoil
              
              # obs
              sel_qoil_gdd_obs_p25     <- c('sel_qoil_gdd_obs_p25'      = as.numeric(elast_weather_p25[1,4]))
              sel_qoil_prec_obs_p25    <- c('sel_qoil_prec_obs_p25'     = as.numeric(elast_weather_p25[2,4]))
              sel_qoil_gddHigh_obs_p25 <- c('sel_qoil_gddHigh_obs_p25'  = as.numeric(elast_weather_p25[3,4]))
              sel_qoil_dd_obs_p25      <- c('sel_qoil_dd_obs_p25'       = as.numeric(elast_weather_p25[4,4]))
              
              #past
              sel_qoil_gdd_past_p25 <- c('sel_qoil_gdd_past_p25'          = as.numeric(elast_weather_p25[5,4]))
              sel_qoil_prec_past_p25  <- c('sel_qoil_prec_past_p25'       = as.numeric(elast_weather_p25[6,4]))
              sel_qoil_gddHigh_past_p25  <- c('sel_qoil_gddHigh_past_p25' = as.numeric(elast_weather_p25[7,4]))
              sel_qoil_dd_past_p25  <- c('sel_qoil_dd_past_p25'           = as.numeric(elast_weather_p25[8,4]))
              
              
              #qroots
              
              # obs
              sel_qroots_gdd_obs_p25     <- c('sel_qroots_gdd_obs_p25'      = as.numeric(elast_weather_p25[1,5]))
              sel_qroots_prec_obs_p25    <- c('sel_qroots_prec_obs_p25'     = as.numeric(elast_weather_p25[2,5]))
              sel_qroots_gddHigh_obs_p25 <- c('sel_qroots_gddHigh_obs_p25'  = as.numeric(elast_weather_p25[3,5]))
              sel_qroots_dd_obs_p25      <- c('sel_qroots_dd_obs_p25'       = as.numeric(elast_weather_p25[4,5]))
              
              #past
              sel_qroots_gdd_past_p25 <- c('sel_qroots_gdd_past_p25'          = as.numeric(elast_weather_p25[5,5]))
              sel_qroots_prec_past_p25  <- c('sel_qroots_prec_past_p25'       = as.numeric(elast_weather_p25[6,5]))
              sel_qroots_gddHigh_past_p25  <- c('sel_qroots_gddHigh_past_p25' = as.numeric(elast_weather_p25[7,5]))
              sel_qroots_dd_past_p25  <- c('sel_qroots_dd_past_p25'           = as.numeric(elast_weather_p25[8,5]))
              
              #corn
              
              # obs
              sel_qcorn_gdd_obs_p25     <- c('sel_qcorn_gdd_obs_p25'      = as.numeric(elast_weather_p25[1,6]))
              sel_qcorn_prec_obs_p25    <- c('sel_qcorn_prec_obs_p25'     = as.numeric(elast_weather_p25[2,6]))
              sel_qcorn_gddHigh_obs_p25 <- c('sel_qcorn_gddHigh_obs_p25'  = as.numeric(elast_weather_p25[3,6]))
              sel_qcorn_dd_obs_p25      <- c('sel_qcorn_dd_obs_p25'       = as.numeric(elast_weather_p25[4,6]))
              
              #past
              sel_qcorn_gdd_past_p25 <- c('sel_qcorn_gdd_past_p25'          = as.numeric(elast_weather_p25[5,6]))
              sel_qcorn_prec_past_p25  <- c('sel_qcorn_prec_past_p25'       = as.numeric(elast_weather_p25[6,6]))
              sel_qcorn_gddHigh_past_p25  <- c('sel_qcorn_gddHigh_past_p25' = as.numeric(elast_weather_p25[7,6]))
              sel_qcorn_dd_past_p25  <- c('sel_qcorn_dd_past_p25'           = as.numeric(elast_weather_p25[8,6]))
              
              #xfert
              
              # obs
              sel_xfert_gdd_obs_p25     <- c('sel_xfert_gdd_obs_p25'      = as.numeric(elast_weather_p25[1,7]))
              sel_xfert_prec_obs_p25    <- c('sel_xfert_prec_obs_p25'     = as.numeric(elast_weather_p25[2,7]))
              sel_xfert_gddHigh_obs_p25 <- c('sel_xfert_gddHigh_obs_p25'  = as.numeric(elast_weather_p25[3,7]))
              sel_xfert_dd_obs_p25      <- c('sel_xfert_dd_obs_p25'       = as.numeric(elast_weather_p25[4,7]))
              
              #past
              sel_xfert_gdd_past_p25     <- c('sel_xfert_gdd_past_p25'     = as.numeric(elast_weather_p25[5,7]))
              sel_xfert_prec_past_p25    <- c('sel_xfert_prec_past_p25'    = as.numeric(elast_weather_p25[6,7]))
              sel_xfert_gddHigh_past_p25 <- c('sel_xfert_gddHigh_past_p25' = as.numeric(elast_weather_p25[7,7]))
              sel_xfert_dd_past_p25      <- c('sel_xfert_dd_past_p25'      = as.numeric(elast_weather_p25[8,7]))
              
          # WEATHER AT P25
              
              #qcer
              
              # obs
              sel_qcer_gdd_obs_p75     <- c('sel_qcer_gdd_obs_p75'      = as.numeric(elast_weather_p75[1,2]))
              sel_qcer_prec_obs_p75    <- c('sel_qcer_prec_obs_p75'     = as.numeric(elast_weather_p75[2,2]))
              sel_qcer_gddHigh_obs_p75 <- c('sel_qcer_gddHigh_obs_p75'  = as.numeric(elast_weather_p75[3,2]))
              sel_qcer_dd_obs_p75      <- c('sel_qcer_dd_obs_p75'       = as.numeric(elast_weather_p75[4,2]))
              
              #past
              sel_qcer_gdd_past_p75 <- c('sel_qcer_gdd_past_p75'          = as.numeric(elast_weather_p75[5,2]))
              sel_qcer_prec_past_p75  <- c('sel_qcer_prec_past_p75'       = as.numeric(elast_weather_p75[6,2]))
              sel_qcer_gddHigh_past_p75  <- c('sel_qcer_gddHigh_past_p75' = as.numeric(elast_weather_p75[7,2]))
              sel_qcer_dd_past_p75  <- c('sel_qcer_dd_past_p75'           = as.numeric(elast_weather_p75[8,2]))
              
              #qprot
              
              # obs
              sel_qprot_gdd_obs_p75     <- c('sel_qprot_gdd_obs_p75'      = as.numeric(elast_weather_p75[1,3]))
              sel_qprot_prec_obs_p75    <- c('sel_qprot_prec_obs_p75'     = as.numeric(elast_weather_p75[2,3]))
              sel_qprot_gddHigh_obs_p75 <- c('sel_qprot_gddHigh_obs_p75'  = as.numeric(elast_weather_p75[3,3]))
              sel_qprot_dd_obs_p75      <- c('sel_qprot_dd_obs_p75'       = as.numeric(elast_weather_p75[4,3]))
              
              #past
              sel_qprot_gdd_past_p75 <- c('sel_qprot_gdd_past_p75'          = as.numeric(elast_weather_p75[5,3]))
              sel_qprot_prec_past_p75  <- c('sel_qprot_prec_past_p75'       = as.numeric(elast_weather_p75[6,3]))
              sel_qprot_gddHigh_past_p75  <- c('sel_qprot_gddHigh_past_p75' = as.numeric(elast_weather_p75[7,3]))
              sel_qprot_dd_past_p75  <- c('sel_qprot_dd_past_p75'           = as.numeric(elast_weather_p75[8,3]))
              
              #qoil
              
              # obs
              sel_qoil_gdd_obs_p75     <- c('sel_qoil_gdd_obs_p75'      = as.numeric(elast_weather_p75[1,4]))
              sel_qoil_prec_obs_p75    <- c('sel_qoil_prec_obs_p75'     = as.numeric(elast_weather_p75[2,4]))
              sel_qoil_gddHigh_obs_p75 <- c('sel_qoil_gddHigh_obs_p75'  = as.numeric(elast_weather_p75[3,4]))
              sel_qoil_dd_obs_p75      <- c('sel_qoil_dd_obs_p75'       = as.numeric(elast_weather_p75[4,4]))
              
              #past
              sel_qoil_gdd_past_p75 <- c('sel_qoil_gdd_past_p75'          = as.numeric(elast_weather_p75[5,4]))
              sel_qoil_prec_past_p75  <- c('sel_qoil_prec_past_p75'       = as.numeric(elast_weather_p75[6,4]))
              sel_qoil_gddHigh_past_p75  <- c('sel_qoil_gddHigh_past_p75' = as.numeric(elast_weather_p75[7,4]))
              sel_qoil_dd_past_p75  <- c('sel_qoil_dd_past_p75'           = as.numeric(elast_weather_p75[8,4]))
              
              
              #qroots
              
              # obs
              sel_qroots_gdd_obs_p75     <- c('sel_qroots_gdd_obs_p75'      = as.numeric(elast_weather_p75[1,5]))
              sel_qroots_prec_obs_p75    <- c('sel_qroots_prec_obs_p75'     = as.numeric(elast_weather_p75[2,5]))
              sel_qroots_gddHigh_obs_p75 <- c('sel_qroots_gddHigh_obs_p75'  = as.numeric(elast_weather_p75[3,5]))
              sel_qroots_dd_obs_p75      <- c('sel_qroots_dd_obs_p75'       = as.numeric(elast_weather_p75[4,5]))
              
              #past
              sel_qroots_gdd_past_p75 <- c('sel_qroots_gdd_past_p75'          = as.numeric(elast_weather_p75[5,5]))
              sel_qroots_prec_past_p75  <- c('sel_qroots_prec_past_p75'       = as.numeric(elast_weather_p75[6,5]))
              sel_qroots_gddHigh_past_p75  <- c('sel_qroots_gddHigh_past_p75' = as.numeric(elast_weather_p75[7,5]))
              sel_qroots_dd_past_p75  <- c('sel_qroots_dd_past_p75'           = as.numeric(elast_weather_p75[8,5]))
              
              #corn
              
              # obs
              sel_qcorn_gdd_obs_p75     <- c('sel_qcorn_gdd_obs_p75'      = as.numeric(elast_weather_p75[1,6]))
              sel_qcorn_prec_obs_p75    <- c('sel_qcorn_prec_obs_p75'     = as.numeric(elast_weather_p75[2,6]))
              sel_qcorn_gddHigh_obs_p75 <- c('sel_qcorn_gddHigh_obs_p75'  = as.numeric(elast_weather_p75[3,6]))
              sel_qcorn_dd_obs_p75      <- c('sel_qcorn_dd_obs_p75'       = as.numeric(elast_weather_p75[4,6]))
              
              #past
              sel_qcorn_gdd_past_p75 <- c('sel_qcorn_gdd_past_p75'          = as.numeric(elast_weather_p75[5,6]))
              sel_qcorn_prec_past_p75  <- c('sel_qcorn_prec_past_p75'       = as.numeric(elast_weather_p75[6,6]))
              sel_qcorn_gddHigh_past_p75  <- c('sel_qcorn_gddHigh_past_p75' = as.numeric(elast_weather_p75[7,6]))
              sel_qcorn_dd_past_p75  <- c('sel_qcorn_dd_past_p75'           = as.numeric(elast_weather_p75[8,6]))
              
              #xfert
              
              # obs
              sel_xfert_gdd_obs_p75     <- c('sel_xfert_gdd_obs_p75'      = as.numeric(elast_weather_p75[1,7]))
              sel_xfert_prec_obs_p75    <- c('sel_xfert_prec_obs_p75'     = as.numeric(elast_weather_p75[2,7]))
              sel_xfert_gddHigh_obs_p75 <- c('sel_xfert_gddHigh_obs_p75'  = as.numeric(elast_weather_p75[3,7]))
              sel_xfert_dd_obs_p75      <- c('sel_xfert_dd_obs_p75'       = as.numeric(elast_weather_p75[4,7]))
              
              #past
              sel_xfert_gdd_past_p75     <- c('sel_xfert_gdd_past_p75'     = as.numeric(elast_weather_p75[5,7]))
              sel_xfert_prec_past_p75    <- c('sel_xfert_prec_past_p75'    = as.numeric(elast_weather_p75[6,7]))
              sel_xfert_gddHigh_past_p75 <- c('sel_xfert_gddHigh_past_p75' = as.numeric(elast_weather_p75[7,7]))
              sel_xfert_dd_past_p75      <- c('sel_xfert_dd_past_p75'      = as.numeric(elast_weather_p75[8,7]))         
      
      
    
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
                sel_xfert_dd_past,
                sel_qcer_gdd_obs_p25, #weather elasticities p25
                sel_qcer_prec_obs_p25,
                sel_qcer_gddHigh_obs_p25,
                sel_qcer_dd_obs_p25,
                sel_qprot_gdd_obs_p25,
                sel_qprot_prec_obs_p25,
                sel_qprot_gddHigh_obs_p25,
                sel_qprot_dd_obs_p25,
                sel_qoil_gdd_obs_p25,
                sel_qoil_prec_obs_p25,
                sel_qoil_gddHigh_obs_p25,
                sel_qoil_dd_obs_p25,
                sel_qroots_gdd_obs_p25,
                sel_qroots_prec_obs_p25,
                sel_qroots_gddHigh_obs_p25,
                sel_qroots_dd_obs_p25,
                sel_qcorn_gdd_obs_p25,
                sel_qcorn_prec_obs_p25,
                sel_qcorn_gddHigh_obs_p25,
                sel_qcorn_dd_obs_p25,
                sel_xfert_gdd_obs_p25,
                sel_xfert_prec_obs_p25,
                sel_xfert_gddHigh_obs_p25,
                sel_xfert_dd_obs_p25,
                sel_qcer_gdd_past_p25,
                sel_qcer_prec_past_p25,
                sel_qcer_gddHigh_past_p25,
                sel_qcer_dd_past_p25,
                sel_qprot_gdd_past_p25,
                sel_qprot_prec_past_p25,
                sel_qprot_gddHigh_past_p25,
                sel_qprot_dd_past_p25,
                sel_qoil_gdd_past_p25,
                sel_qoil_prec_past_p25,
                sel_qoil_gddHigh_past_p25,
                sel_qoil_dd_past_p25,
                sel_qroots_gdd_past_p25,
                sel_qroots_prec_past_p25,
                sel_qroots_gddHigh_past_p25,
                sel_qroots_dd_past_p25,
                sel_qcorn_gdd_past_p25,
                sel_qcorn_prec_past_p25,
                sel_qcorn_gddHigh_past_p25,
                sel_qcorn_dd_past_p25,
                sel_xfert_gdd_past_p25,
                sel_xfert_prec_past_p25,
                sel_xfert_gddHigh_past_p25,
                sel_xfert_dd_past_p25,
                sel_qcer_gdd_obs_p75, #weather elasticities p75
                sel_qcer_prec_obs_p75,
                sel_qcer_gddHigh_obs_p75,
                sel_qcer_dd_obs_p75,
                sel_qprot_gdd_obs_p75,
                sel_qprot_prec_obs_p75,
                sel_qprot_gddHigh_obs_p75,
                sel_qprot_dd_obs_p75,
                sel_qoil_gdd_obs_p75,
                sel_qoil_prec_obs_p75,
                sel_qoil_gddHigh_obs_p75,
                sel_qoil_dd_obs_p75,
                sel_qroots_gdd_obs_p75,
                sel_qroots_prec_obs_p75,
                sel_qroots_gddHigh_obs_p75,
                sel_qroots_dd_obs_p75,
                sel_qcorn_gdd_obs_p75,
                sel_qcorn_prec_obs_p75,
                sel_qcorn_gddHigh_obs_p75,
                sel_qcorn_dd_obs_p75,
                sel_xfert_gdd_obs_p75,
                sel_xfert_prec_obs_p75,
                sel_xfert_gddHigh_obs_p75,
                sel_xfert_dd_obs_p75,
                sel_qcer_gdd_past_p75,
                sel_qcer_prec_past_p75,
                sel_qcer_gddHigh_past_p75,
                sel_qcer_dd_past_p75,
                sel_qprot_gdd_past_p75,
                sel_qprot_prec_past_p75,
                sel_qprot_gddHigh_past_p75,
                sel_qprot_dd_past_p75,
                sel_qoil_gdd_past_p75,
                sel_qoil_prec_past_p75,
                sel_qoil_gddHigh_past_p75,
                sel_qoil_dd_past_p75,
                sel_qroots_gdd_past_p75,
                sel_qroots_prec_past_p75,
                sel_qroots_gddHigh_past_p75,
                sel_qroots_dd_past_p75,
                sel_qcorn_gdd_past_p75,
                sel_qcorn_prec_past_p75,
                sel_qcorn_gddHigh_past_p75,
                sel_qcorn_dd_past_p75,
                sel_xfert_gdd_past_p75,
                sel_xfert_prec_past_p75,
                sel_xfert_gddHigh_past_p75,
                sel_xfert_dd_past_p75) 
  
        result 
        
    }, error=function(err) {rep(NA, 1350)} ) # insert NA's if there is an optimization error
  }
  
}


#-------------------------------------------#
#### Obtain and store the actual results ####
#-------------------------------------------#

# run one draw with the original data
system.time(boot_actresults <- my.boot(data=df_farm, nrep=1, cluster=df_farm$key, nCores=1, actual=TRUE))

# save actual result
save(boot_actresults, file="rOutput/bootresults_het_actual.Rda")
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
save(boot_results, file="rOutput/bootresults_het.Rda")
boot_results <- NULL


#----------------------------------#
#### 2) Presentation of results ####
#----------------------------------#

# Load bootstrap results
load("rOutput/bootresults_het.Rda")
load("rOutput/bootresults_het_actual.Rda")

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
coef_iCereals <- boot_summary[c(1:63),]
coef_iProtein <- boot_summary[c(64:126),]
coef_iOilseed <- boot_summary[c(127:189),]
coef_iRoots <- boot_summary[c(190:252),]
coef_iCorn <- boot_summary[c(253:315),]
iCereals_me <- boot_summary[c(316:378),]
iProtein_me <- boot_summary[c(379:441),]
iOilseed_me <- boot_summary[c(442:504),]
iRoots_me <- boot_summary[c(505:567),]
iCorn_me <- boot_summary[c(568:630),]

qqCereals <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "QQcereals_")),]
qqProtein <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "QQprotein_")),]
qqOilseed <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "QQoilseed_")),]
qqRoots   <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "QQroots_")),]
qqCorn    <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "QQcorn_")),]
nxFert    <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "NXfert_")),]

PriceElast_cereals <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "el_qcer_")),]
PriceElast_protein <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "el_qprot_")),]
PriceElast_oilseed <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "el_qoil_")),]
PriceElast_roots   <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "el_qroots_")),]
PriceElast_corn <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "el_qcorn_")),]
PriceElast_fert    <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "el_xfert_")),]
PriceElast_others  <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "el_xotherinp_")),]

WeatherElast_cereals_sm  <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qcer_") &
                                            (endsWith(x = boot_summary$Variable, suffix = "obs") |
                                               endsWith(x = boot_summary$Variable, suffix = "past"))),]
WeatherElast_cereals_p25 <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qcer_") &
                                            endsWith(x = boot_summary$Variable, suffix = "p25")),]
WeatherElast_cereals_p75 <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qcer_") &
                                            endsWith(x = boot_summary$Variable, suffix = "p75")),]

WeatherElast_protein_sm <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qprot_") &
                                           (endsWith(x = boot_summary$Variable, suffix = "obs") |
                                              endsWith(x = boot_summary$Variable, suffix = "past"))),]
WeatherElast_protein_p25 <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qprot_") &
                                            endsWith(x = boot_summary$Variable, suffix = "p25")),]
WeatherElast_protein_p75 <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qprot_") &
                                            endsWith(x = boot_summary$Variable, suffix = "p75")),]


WeatherElast_oilseed_sm <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qoil_") &
                                           (endsWith(x = boot_summary$Variable, suffix = "obs") |
                                              endsWith(x = boot_summary$Variable, suffix = "past"))),]
WeatherElast_oilseed_p25 <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qoil_") &
                                            endsWith(x = boot_summary$Variable, suffix = "p25")),]
WeatherElast_oilseed_p75 <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qoil_") &
                                            endsWith(x = boot_summary$Variable, suffix = "p75")),]


WeatherElast_roots_sm   <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qroots_") &
                                           (endsWith(x = boot_summary$Variable, suffix = "obs") |
                                              endsWith(x = boot_summary$Variable, suffix = "past"))),]
WeatherElast_roots_p25   <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qroots_") &
                                            endsWith(x = boot_summary$Variable, suffix = "p25")),]
WeatherElast_roots_p75   <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qroots_") &
                                            endsWith(x = boot_summary$Variable, suffix = "p75")),]


WeatherElast_corn_sm <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qcorn_") &
                                           (endsWith(x = boot_summary$Variable, suffix = "obs") |
                                              endsWith(x = boot_summary$Variable, suffix = "past"))),]
WeatherElast_corn_p25 <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qcorn_") &
                                            endsWith(x = boot_summary$Variable, suffix = "p25")),]
WeatherElast_corn_p75 <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_qcorn_") &
                                            endsWith(x = boot_summary$Variable, suffix = "p75")),]


WeatherElast_fert_sm    <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_xfert_") &
                                           (endsWith(x = boot_summary$Variable, suffix = "obs") |
                                              endsWith(x = boot_summary$Variable, suffix = "past"))),]
WeatherElast_fert_p25    <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_xfert_") &
                                            endsWith(x = boot_summary$Variable, suffix = "p25")),]
WeatherElast_fert_p75    <- boot_summary[(startsWith(x = boot_summary$Variable, prefix = "sel_xfert_") &
                                            endsWith(x = boot_summary$Variable, suffix = "p75")),]


#----------------------------------------------------#
#### Table: Parameter estimates for probit models ####
#----------------------------------------------------#

#combine probit results
tab_parprobit <- cbind(coef_iCereals[,c("Variable","Original","L95","U95")],
                       coef_iProtein[,c("Original","L95","U95")],
                       coef_iOilseed[,c("Original","L95","U95")],
                       coef_iRoots[,c("Original","L95","U95")],
                       coef_iCorn[,c("Original","L95","U95")])

tab_parprobit <- tab_parprobit[1:32,]

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
                            "Gdd_4to10", "Prec_4to10", "GddHigh_4to10", "DD_4to10",
                            "Gdd_1to3 x Land", "Prec_1to3 x Land", "GddHigh_1to3 x Land", "DD_1to3 x Land",
                            "Gdd_4to10 x Land", "Prec_4to10 x Land", "GddHigh_4to10 x Land", "DD_4to10 x Land")

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
write_xlsx(tab_parprobit,"Tables/Table_S5.1.xlsx")


#------------------------------------------------------#
#### Table: Regression results structural equations ####
#------------------------------------------------------#

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
                             "Gdd x Land", "Prec x Land", "GddHigh x Land", "DD x Land",
                             "Gdd_1to3 x Land", "Prec_1to3 x Land", "GddHigh_1to3 x Land", "DD_1to3 x Land",
                             "Gdd_4to10 x Land", "Prec_4to10 x Land", "GddHigh_4to10 x Land", "DD_4to10 x Land",
                             "Gdd x Gdd_1to3", "Prec x Prec_1to3", "GddHigh x GddHigh_1to3", "DD x DD_1to3",
                             "Gdd x Gdd_4to10", "Prec x Prec_4to10", "GddHigh x GddHigh_4to10", "DD x DD_4to10",
                             "phi Cereals",
                             "Intercept", "P Cereals", "P Protein", "P Oilseed crops", "P Root crops", "P Corn", "R Fert", 
                             "K Land", "K Labor", "K Capital", "Trend", "Trend^2", 
                             "Gdd", "Prec", "GddHigh", "DD",
                             "Gdd_1to3", "Prec_1to3", "GddHigh_1to3", "DD_1to3",
                             "Gdd_4to10", "Prec_4to10", "GddHigh_4to10", "DD_4to10",
                             "Gdd x Land", "Prec x Land", "GddHigh x Land", "DD x Land",
                             "Gdd_1to3 x Land", "Prec_1to3 x Land", "GddHigh_1to3 x Land", "DD_1to3 x Land",
                             "Gdd_4to10 x Land", "Prec_4to10 x Land", "GddHigh_4to10 x Land", "DD_4to10 x Land",
                             "Gdd x Gdd_1to3", "Prec x Prec_1to3", "GddHigh x GddHigh_1to3", "DD x DD_1to3",
                             "Gdd x Gdd_4to10", "Prec x Prec_4to10", "GddHigh x GddHigh_4to10", "DD x DD_4to10",
                             "phi Protein",
                             "Intercept", "P Cereals", "P Protein", "P Oilseed crops", "P Root crops", "P Corn", "R Fert", 
                             "K Land", "K Labor", "K Capital", "Trend", "Trend^2", 
                             "Gdd", "Prec", "GddHigh", "DD",
                             "Gdd_1to3", "Prec_1to3", "GddHigh_1to3", "DD_1to3",
                             "Gdd_4to10", "Prec_4to10", "GddHigh_4to10", "DD_4to10",
                             "Gdd x Land", "Prec x Land", "GddHigh x Land", "DD x Land",
                             "Gdd_1to3 x Land", "Prec_1to3 x Land", "GddHigh_1to3 x Land", "DD_1to3 x Land",
                             "Gdd_4to10 x Land", "Prec_4to10 x Land", "GddHigh_4to10 x Land", "DD_4to10 x Land",
                             "Gdd x Gdd_1to3", "Prec x Prec_1to3", "GddHigh x GddHigh_1to3", "DD x DD_1to3",
                             "Gdd x Gdd_4to10", "Prec x Prec_4to10", "GddHigh x GddHigh_4to10", "DD x DD_4to10",
                             "phi Oilseed",
                             "Intercept", "P Cereals", "P Protein", "P Oilseed crops", "P Root crops", "P Corn", "R Fert", 
                             "K Land", "K Labor", "K Capital", "Trend", "Trend^2", 
                             "Gdd", "Prec", "GddHigh", "DD",
                             "Gdd_1to3", "Prec_1to3", "GddHigh_1to3", "DD_1to3",
                             "Gdd_4to10", "Prec_4to10", "GddHigh_4to10", "DD_4to10",
                             "Gdd x Land", "Prec x Land", "GddHigh x Land", "DD x Land",
                             "Gdd_1to3 x Land", "Prec_1to3 x Land", "GddHigh_1to3 x Land", "DD_1to3 x Land",
                             "Gdd_4to10 x Land", "Prec_4to10 x Land", "GddHigh_4to10 x Land", "DD_4to10 x Land",
                             "Gdd x Gdd_1to3", "Prec x Prec_1to3", "GddHigh x GddHigh_1to3", "DD x DD_1to3",
                             "Gdd x Gdd_4to10", "Prec x Prec_4to10", "GddHigh x GddHigh_4to10", "DD x DD_4to10",
                             "phi Roots",
                             "Intercept", "P Cereals", "P Protein", "P Oilseed crops", "P Root crops", "P Corn", "R Fert", 
                             "K Land", "K Labor", "K Capital", "Trend", "Trend^2", 
                             "Gdd", "Prec", "GddHigh", "DD",
                             "Gdd_1to3", "Prec_1to3", "GddHigh_1to3", "DD_1to3",
                             "Gdd_4to10", "Prec_4to10", "GddHigh_4to10", "DD_4to10",
                             "Gdd x Land", "Prec x Land", "GddHigh x Land", "DD x Land",
                             "Gdd_1to3 x Land", "Prec_1to3 x Land", "GddHigh_1to3 x Land", "DD_1to3 x Land",
                             "Gdd_4to10 x Land", "Prec_4to10 x Land", "GddHigh_4to10 x Land", "DD_4to10 x Land",
                             "Gdd x Gdd_1to3", "Prec x Prec_1to3", "GddHigh x GddHigh_1to3", "DD x DD_1to3",
                             "Gdd x Gdd_4to10", "Prec x Prec_4to10", "GddHigh x GddHigh_4to10", "DD x DD_4to10",
                             "phi Corn",
                             "Intercept", "P Cereals", "P Protein", "P Oilseed crops", "P Root crops", "P Corn", "R Fert", 
                             "K Land", "K Labor", "K Capital", "Trend", "Trend^2", 
                             "Gdd", "Prec", "GddHigh", "DD",
                             "Gdd_1to3", "Prec_1to3", "GddHigh_1to3", "DD_1to3",
                             "Gdd_4to10", "Prec_4to10", "GddHigh_4to10", "DD_4to10",
                             "Gdd x Land", "Prec x Land", "GddHigh x Land", "DD x Land",
                             "Gdd_1to3 x Land", "Prec_1to3 x Land", "GddHigh_1to3 x Land", "DD_1to3 x Land",
                             "Gdd_4to10 x Land", "Prec_4to10 x Land", "GddHigh_4to10 x Land", "DD_4to10 x Land",
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
write_xlsx(tab_structural,"Tables/Table_S5.2.xlsx")

#------------------------------------------------#
#### Table: Own- and cross-price elasticities ####
#------------------------------------------------#

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
write_xlsx(tab_structural,"Tables/Table_S5.3.xlsx")

#---------------------------------------------------#
#### Table: Weather effects on output and inputs ####
#---------------------------------------------------#
#(compare obs vs. past)

# a) AT THE SAMPLE MEAN

#combine qs and xs
tab_WeatherElast_sm <- cbind(WeatherElast_cereals_sm[,c("Variable","Original","L95","U95")],
                             WeatherElast_protein_sm[,c("Original","L95","U95")],
                             WeatherElast_oilseed_sm[,c("Original","L95","U95")],
                             WeatherElast_roots_sm[,c("Original","L95","U95")],
                             WeatherElast_corn_sm[,c("Original","L95","U95")],
                             WeatherElast_fert_sm[,c("Original","L95","U95")])

#name columns
names(tab_WeatherElast_sm) <- c("na", 
                                "Q_Cereals", "L95_cereals", "U95_cereals", 
                                "Q_Protein", "L95_protein", "U95_protein", 
                                "Q_Oilseed", "L95_oilseed", "U95_oilseed", 
                                "Q_Roots", "L95_roots", "U95_roots", 
                                "Q_Corn", "L95_corn", "U95_corn", 
                                "X_Fertilizer", "L95_fert", "U95_fert")

#name rows
tab_WeatherElast_sm$na <- c("gdd_obs","prec_obs","gddHigh_obs","dd_obs",
                            "gdd_past","prec_past","gddHigh_past","dd_past")


#Numbers smaller than 0.0005 should be in scientific format; others in 0.001 etc.
tab_WeatherElast_sm <- as.data.frame(cbind(tab_WeatherElast_sm[,1], sapply(tab_WeatherElast_sm[,-1], function(x) ifelse(abs(x)<0.0005&is.numeric(x),format(x,scientific=TRUE,digits=1),format(round(x, 3), nsmall = 3, scientific=FALSE)) )) )

# remove white space
tab_WeatherElast_sm <- as.data.frame(sapply(tab_WeatherElast_sm, function(x) (gsub(" ", "", x))))

# Merge CIs
tab_WeatherElast_sm$L95_cereals <- paste0(paste0("(",tab_WeatherElast_sm$L95_cereals,";"),
                                          paste0(tab_WeatherElast_sm$U95_cereals,")"))
tab_WeatherElast_sm$U95_cereals <- NULL
tab_WeatherElast_sm <- dplyr::rename(tab_WeatherElast_sm,"CI_cereals" = "L95_cereals")

tab_WeatherElast_sm$L95_protein <- paste0(paste0("(",tab_WeatherElast_sm$L95_protein,";"),
                                          paste0(tab_WeatherElast_sm$U95_protein,")"))
tab_WeatherElast_sm$U95_protein <- NULL
tab_WeatherElast_sm <- dplyr::rename(tab_WeatherElast_sm,"CI_protein" = "L95_protein")

tab_WeatherElast_sm$L95_oilseed <- paste0(paste0("(",tab_WeatherElast_sm$L95_oilseed,";"),
                                          paste0(tab_WeatherElast_sm$U95_oilseed,")"))
tab_WeatherElast_sm$U95_oilseed <- NULL
tab_WeatherElast_sm <- dplyr::rename(tab_WeatherElast_sm,"CI_oilseed" = "L95_oilseed")

tab_WeatherElast_sm$L95_roots <- paste0(paste0("(",tab_WeatherElast_sm$L95_roots,";"),
                                        paste0(tab_WeatherElast_sm$U95_roots,")"))
tab_WeatherElast_sm$U95_roots <- NULL
tab_WeatherElast_sm <- dplyr::rename(tab_WeatherElast_sm,"CI_roots" = "L95_roots")

tab_WeatherElast_sm$L95_corn <- paste0(paste0("(",tab_WeatherElast_sm$L95_corn,";"),
                                          paste0(tab_WeatherElast_sm$U95_corn,")"))
tab_WeatherElast_sm$U95_corn <- NULL
tab_WeatherElast_sm <- dplyr::rename(tab_WeatherElast_sm,"CI_corn" = "L95_corn")

tab_WeatherElast_sm$L95_fert <- paste0(paste0("(",tab_WeatherElast_sm$L95_fert,";"),
                                       paste0(tab_WeatherElast_sm$U95_fert,")"))
tab_WeatherElast_sm$U95_fert <- NULL
tab_WeatherElast_sm <- dplyr::rename(tab_WeatherElast_sm,"CI_fert" = "L95_fert")


#Insert a blank row to every other row
na.df <- data.frame("na" = NA, "Q_Cereals" = NA, "CI_cereals" = NA, "Q_Protein" = NA,
                    "CI_protein" = NA, "Q_Oilseed" = NA, "CI_oilseed" = NA, "Q_Roots" = NA,
                    "CI_roots" = NA, "Q_Corn" = NA, "CI_corn" = NA, "X_Fertilizer" = NA,
                    "CI_fert" = NA)
tab_WeatherElast_sm <- do.call(rbind, apply(tab_WeatherElast_sm, 1, function(x) {rbind(x, na.df)}))

## Shift CIs down one line
tab_WeatherElast_sm$CI_cereals <- c(NA,tab_WeatherElast_sm$CI_cereals[1:nrow(tab_WeatherElast_sm)-1])
tab_WeatherElast_sm$CI_protein <- c(NA,tab_WeatherElast_sm$CI_protein[1:nrow(tab_WeatherElast_sm)-1])
tab_WeatherElast_sm$CI_oilseed <- c(NA,tab_WeatherElast_sm$CI_oilseed[1:nrow(tab_WeatherElast_sm)-1])
tab_WeatherElast_sm$CI_roots <- c(NA,tab_WeatherElast_sm$CI_roots[1:nrow(tab_WeatherElast_sm)-1])
tab_WeatherElast_sm$CI_corn <- c(NA,tab_WeatherElast_sm$CI_corn[1:nrow(tab_WeatherElast_sm)-1])
tab_WeatherElast_sm$CI_fert <- c(NA,tab_WeatherElast_sm$CI_fert[1:nrow(tab_WeatherElast_sm)-1])

# Shifts CIs to the left
tab_WeatherElast_sm$Q_Cereals[is.na(tab_WeatherElast_sm$Q_Cereals)] <- tab_WeatherElast_sm$CI_cereals[is.na(tab_WeatherElast_sm$Q_Cereals)]
tab_WeatherElast_sm$CI_cereals <- NULL

tab_WeatherElast_sm$Q_Protein[is.na(tab_WeatherElast_sm$Q_Protein)] <- tab_WeatherElast_sm$CI_protein[is.na(tab_WeatherElast_sm$Q_Protein)]
tab_WeatherElast_sm$CI_protein <- NULL

tab_WeatherElast_sm$Q_Oilseed[is.na(tab_WeatherElast_sm$Q_Oilseed)] <- tab_WeatherElast_sm$CI_oilseed[is.na(tab_WeatherElast_sm$Q_Oilseed)]
tab_WeatherElast_sm$CI_oilseed <- NULL

tab_WeatherElast_sm$Q_Roots[is.na(tab_WeatherElast_sm$Q_Roots)] <- tab_WeatherElast_sm$CI_roots[is.na(tab_WeatherElast_sm$Q_Roots)]
tab_WeatherElast_sm$CI_roots <- NULL

tab_WeatherElast_sm$Q_Corn[is.na(tab_WeatherElast_sm$Q_Corn)] <- tab_WeatherElast_sm$CI_corn[is.na(tab_WeatherElast_sm$Q_Corn)]
tab_WeatherElast_sm$CI_corn <- NULL

tab_WeatherElast_sm$X_Fertilizer[is.na(tab_WeatherElast_sm$X_Fertilizer)] <- tab_WeatherElast_sm$CI_fert[is.na(tab_WeatherElast_sm$X_Fertilizer)]
tab_WeatherElast_sm$CI_fert <- NULL

# Insert observed and past weather indicators
tab_WeatherElast_sm <- berryFunctions::insertRows(tab_WeatherElast_sm, 1 , new = c("Observed weather","","","","","",""))
tab_WeatherElast_sm <- berryFunctions::insertRows(tab_WeatherElast_sm, nrow(WeatherElast_cereals_sm)+2 , new = c("Past weather","","","","","",""))

# Print
write_xlsx(tab_structural,"Tables/Table_S5.4.xlsx")

# b) FOR SMALL FARMS

#combine qs and xs
tab_WeatherElast_p25 <- cbind(WeatherElast_cereals_p25[,c("Variable","Original","L95","U95")],
                              WeatherElast_protein_p25[,c("Original","L95","U95")],
                              WeatherElast_oilseed_p25[,c("Original","L95","U95")],
                              WeatherElast_roots_p25[,c("Original","L95","U95")],
                              WeatherElast_corn_p25[,c("Original","L95","U95")],
                              WeatherElast_fert_p25[,c("Original","L95","U95")])

#name columns
names(tab_WeatherElast_p25) <- c("na", 
                                 "Q_Cereals", "L95_cereals", "U95_cereals", 
                                 "Q_Protein", "L95_protein", "U95_protein", 
                                 "Q_Oilseed", "L95_oilseed", "U95_oilseed", 
                                 "Q_Roots", "L95_roots", "U95_roots", 
                                 "Q_Corn", "L95_corn", "U95_corn", 
                                 "X_Fertilizer", "L95_fert", "U95_fert")

#name rows
tab_WeatherElast_p25$na <- c("gdd_obs","prec_obs","gddHigh_obs","dd_obs",
                             "gdd_past","prec_past","gddHigh_past","dd_past")


#Numbers smaller than 0.0005 should be in scientific format; others in 0.001 etc.
tab_WeatherElast_p25 <- as.data.frame(cbind(tab_WeatherElast_p25[,1], sapply(tab_WeatherElast_p25[,-1], function(x) ifelse(abs(x)<0.0005&is.numeric(x),format(x,scientific=TRUE,digits=1),format(round(x, 3), nsmall = 3, scientific=FALSE)) )) )

# remove white space
tab_WeatherElast_p25 <- as.data.frame(sapply(tab_WeatherElast_p25, function(x) (gsub(" ", "", x))))

# Merge CIs
tab_WeatherElast_p25$L95_cereals <- paste0(paste0("(",tab_WeatherElast_p25$L95_cereals,";"),
                                           paste0(tab_WeatherElast_p25$U95_cereals,")"))
tab_WeatherElast_p25$U95_cereals <- NULL
tab_WeatherElast_p25 <- dplyr::rename(tab_WeatherElast_p25,"CI_cereals" = "L95_cereals")

tab_WeatherElast_p25$L95_protein <- paste0(paste0("(",tab_WeatherElast_p25$L95_protein,";"),
                                           paste0(tab_WeatherElast_p25$U95_protein,")"))
tab_WeatherElast_p25$U95_protein <- NULL
tab_WeatherElast_p25 <- dplyr::rename(tab_WeatherElast_p25,"CI_protein" = "L95_protein")

tab_WeatherElast_p25$L95_oilseed <- paste0(paste0("(",tab_WeatherElast_p25$L95_oilseed,";"),
                                           paste0(tab_WeatherElast_p25$U95_oilseed,")"))
tab_WeatherElast_p25$U95_oilseed <- NULL
tab_WeatherElast_p25 <- dplyr::rename(tab_WeatherElast_p25,"CI_oilseed" = "L95_oilseed")

tab_WeatherElast_p25$L95_roots <- paste0(paste0("(",tab_WeatherElast_p25$L95_roots,";"),
                                         paste0(tab_WeatherElast_p25$U95_roots,")"))
tab_WeatherElast_p25$U95_roots <- NULL
tab_WeatherElast_p25 <- dplyr::rename(tab_WeatherElast_p25,"CI_roots" = "L95_roots")

tab_WeatherElast_p25$L95_corn <- paste0(paste0("(",tab_WeatherElast_p25$L95_corn,";"),
                                           paste0(tab_WeatherElast_p25$U95_corn,")"))
tab_WeatherElast_p25$U95_corn <- NULL
tab_WeatherElast_p25 <- dplyr::rename(tab_WeatherElast_p25,"CI_corn" = "L95_corn")

tab_WeatherElast_p25$L95_fert <- paste0(paste0("(",tab_WeatherElast_p25$L95_fert,";"),
                                        paste0(tab_WeatherElast_p25$U95_fert,")"))
tab_WeatherElast_p25$U95_fert <- NULL
tab_WeatherElast_p25 <- dplyr::rename(tab_WeatherElast_p25,"CI_fert" = "L95_fert")


#Insert a blank row to every other row
na.df <- data.frame("na" = NA, "Q_Cereals" = NA, "CI_cereals" = NA, "Q_Protein" = NA,
                    "CI_protein" = NA, "Q_Oilseed" = NA, "CI_oilseed" = NA, "Q_Roots" = NA,
                    "CI_roots" = NA, "Q_Corn" = NA, "CI_corn" = NA, "X_Fertilizer" = NA,
                    "CI_fert" = NA)
tab_WeatherElast_p25 <- do.call(rbind, apply(tab_WeatherElast_p25, 1, function(x) {rbind(x, na.df)}))

## Shift CIs down one line
tab_WeatherElast_p25$CI_cereals <- c(NA,tab_WeatherElast_p25$CI_cereals[1:nrow(tab_WeatherElast_p25)-1])
tab_WeatherElast_p25$CI_protein <- c(NA,tab_WeatherElast_p25$CI_protein[1:nrow(tab_WeatherElast_p25)-1])
tab_WeatherElast_p25$CI_oilseed <- c(NA,tab_WeatherElast_p25$CI_oilseed[1:nrow(tab_WeatherElast_p25)-1])
tab_WeatherElast_p25$CI_roots <- c(NA,tab_WeatherElast_p25$CI_roots[1:nrow(tab_WeatherElast_p25)-1])
tab_WeatherElast_p25$CI_corn <- c(NA,tab_WeatherElast_p25$CI_corn[1:nrow(tab_WeatherElast_p25)-1])
tab_WeatherElast_p25$CI_fert <- c(NA,tab_WeatherElast_p25$CI_fert[1:nrow(tab_WeatherElast_p25)-1])

# Shifts CIs to the left
tab_WeatherElast_p25$Q_Cereals[is.na(tab_WeatherElast_p25$Q_Cereals)] <- tab_WeatherElast_p25$CI_cereals[is.na(tab_WeatherElast_p25$Q_Cereals)]
tab_WeatherElast_p25$CI_cereals <- NULL

tab_WeatherElast_p25$Q_Protein[is.na(tab_WeatherElast_p25$Q_Protein)] <- tab_WeatherElast_p25$CI_protein[is.na(tab_WeatherElast_p25$Q_Protein)]
tab_WeatherElast_p25$CI_protein <- NULL

tab_WeatherElast_p25$Q_Oilseed[is.na(tab_WeatherElast_p25$Q_Oilseed)] <- tab_WeatherElast_p25$CI_oilseed[is.na(tab_WeatherElast_p25$Q_Oilseed)]
tab_WeatherElast_p25$CI_oilseed <- NULL

tab_WeatherElast_p25$Q_Roots[is.na(tab_WeatherElast_p25$Q_Roots)] <- tab_WeatherElast_p25$CI_roots[is.na(tab_WeatherElast_p25$Q_Roots)]
tab_WeatherElast_p25$CI_roots <- NULL

tab_WeatherElast_p25$Q_Corn[is.na(tab_WeatherElast_p25$Q_Corn)] <- tab_WeatherElast_p25$CI_corn[is.na(tab_WeatherElast_p25$Q_Corn)]
tab_WeatherElast_p25$CI_corn <- NULL

tab_WeatherElast_p25$X_Fertilizer[is.na(tab_WeatherElast_p25$X_Fertilizer)] <- tab_WeatherElast_p25$CI_fert[is.na(tab_WeatherElast_p25$X_Fertilizer)]
tab_WeatherElast_p25$CI_fert <- NULL

# Insert observed and past weather indicators
tab_WeatherElast_p25 <- berryFunctions::insertRows(tab_WeatherElast_p25, 1 , new = c("Observed weather","","","","","",""))
tab_WeatherElast_p25 <- berryFunctions::insertRows(tab_WeatherElast_p25, nrow(WeatherElast_cereals_p25)+2 , new = c("Past weather","","","","","",""))

# Print
write_xlsx(tab_structural,"Tables/Table_S5.5.xlsx")


# b) FOR LARGE FARMS

#combine qs and xs
tab_WeatherElast_p75 <- cbind(WeatherElast_cereals_p75[,c("Variable","Original","L95","U95")],
                              WeatherElast_protein_p75[,c("Original","L95","U95")],
                              WeatherElast_oilseed_p75[,c("Original","L95","U95")],
                              WeatherElast_roots_p75[,c("Original","L95","U95")],
                              WeatherElast_corn_p75[,c("Original","L95","U95")],
                              WeatherElast_fert_p75[,c("Original","L95","U95")])

#name columns
names(tab_WeatherElast_p75) <- c("na", 
                                 "Q_Cereals", "L95_cereals", "U95_cereals", 
                                 "Q_Protein", "L95_protein", "U95_protein", 
                                 "Q_Oilseed", "L95_oilseed", "U95_oilseed", 
                                 "Q_Roots", "L95_roots", "U95_roots", 
                                 "Q_Corn", "L95_corn", "U95_corn", 
                                 "X_Fertilizer", "L95_fert", "U95_fert")

#name rows
tab_WeatherElast_p75$na <- c("gdd_obs","prec_obs","gddHigh_obs","dd_obs",
                             "gdd_past","prec_past","gddHigh_past","dd_past")


#Numbers smaller than 0.0005 should be in scientific format; others in 0.001 etc.
tab_WeatherElast_p75 <- as.data.frame(cbind(tab_WeatherElast_p75[,1], sapply(tab_WeatherElast_p75[,-1], function(x) ifelse(abs(x)<0.0005&is.numeric(x),format(x,scientific=TRUE,digits=1),format(round(x, 3), nsmall = 3, scientific=FALSE)) )) )

# remove white space
tab_WeatherElast_p75 <- as.data.frame(sapply(tab_WeatherElast_p75, function(x) (gsub(" ", "", x))))

# Merge CIs
tab_WeatherElast_p75$L95_cereals <- paste0(paste0("(",tab_WeatherElast_p75$L95_cereals,";"),
                                           paste0(tab_WeatherElast_p75$U95_cereals,")"))
tab_WeatherElast_p75$U95_cereals <- NULL
tab_WeatherElast_p75 <- dplyr::rename(tab_WeatherElast_p75,"CI_cereals" = "L95_cereals")

tab_WeatherElast_p75$L95_protein <- paste0(paste0("(",tab_WeatherElast_p75$L95_protein,";"),
                                           paste0(tab_WeatherElast_p75$U95_protein,")"))
tab_WeatherElast_p75$U95_protein <- NULL
tab_WeatherElast_p75 <- dplyr::rename(tab_WeatherElast_p75,"CI_protein" = "L95_protein")

tab_WeatherElast_p75$L95_oilseed <- paste0(paste0("(",tab_WeatherElast_p75$L95_oilseed,";"),
                                           paste0(tab_WeatherElast_p75$U95_oilseed,")"))
tab_WeatherElast_p75$U95_oilseed <- NULL
tab_WeatherElast_p75 <- dplyr::rename(tab_WeatherElast_p75,"CI_oilseed" = "L95_oilseed")

tab_WeatherElast_p75$L95_roots <- paste0(paste0("(",tab_WeatherElast_p75$L95_roots,";"),
                                         paste0(tab_WeatherElast_p75$U95_roots,")"))
tab_WeatherElast_p75$U95_roots <- NULL
tab_WeatherElast_p75 <- dplyr::rename(tab_WeatherElast_p75,"CI_roots" = "L95_roots")

tab_WeatherElast_p75$L95_corn <- paste0(paste0("(",tab_WeatherElast_p75$L95_corn,";"),
                                           paste0(tab_WeatherElast_p75$U95_corn,")"))
tab_WeatherElast_p75$U95_corn <- NULL
tab_WeatherElast_p75 <- dplyr::rename(tab_WeatherElast_p75,"CI_corn" = "L95_corn")

tab_WeatherElast_p75$L95_fert <- paste0(paste0("(",tab_WeatherElast_p75$L95_fert,";"),
                                        paste0(tab_WeatherElast_p75$U95_fert,")"))
tab_WeatherElast_p75$U95_fert <- NULL
tab_WeatherElast_p75 <- dplyr::rename(tab_WeatherElast_p75,"CI_fert" = "L95_fert")


#Insert a blank row to every other row
na.df <- data.frame("na" = NA, "Q_Cereals" = NA, "CI_cereals" = NA, "Q_Protein" = NA,
                    "CI_protein" = NA, "Q_Oilseed" = NA, "CI_oilseed" = NA, "Q_Roots" = NA,
                    "CI_roots" = NA, "Q_Corn" = NA, "CI_corn" = NA, "X_Fertilizer" = NA,
                    "CI_fert" = NA)
tab_WeatherElast_p75 <- do.call(rbind, apply(tab_WeatherElast_p75, 1, function(x) {rbind(x, na.df)}))

## Shift CIs down one line
tab_WeatherElast_p75$CI_cereals <- c(NA,tab_WeatherElast_p75$CI_cereals[1:nrow(tab_WeatherElast_p75)-1])
tab_WeatherElast_p75$CI_protein <- c(NA,tab_WeatherElast_p75$CI_protein[1:nrow(tab_WeatherElast_p75)-1])
tab_WeatherElast_p75$CI_oilseed <- c(NA,tab_WeatherElast_p75$CI_oilseed[1:nrow(tab_WeatherElast_p75)-1])
tab_WeatherElast_p75$CI_roots <- c(NA,tab_WeatherElast_p75$CI_roots[1:nrow(tab_WeatherElast_p75)-1])
tab_WeatherElast_p75$CI_corn <- c(NA,tab_WeatherElast_p75$CI_corn[1:nrow(tab_WeatherElast_p75)-1])
tab_WeatherElast_p75$CI_fert <- c(NA,tab_WeatherElast_p75$CI_fert[1:nrow(tab_WeatherElast_p75)-1])

# Shifts CIs to the left
tab_WeatherElast_p75$Q_Cereals[is.na(tab_WeatherElast_p75$Q_Cereals)] <- tab_WeatherElast_p75$CI_cereals[is.na(tab_WeatherElast_p75$Q_Cereals)]
tab_WeatherElast_p75$CI_cereals <- NULL

tab_WeatherElast_p75$Q_Protein[is.na(tab_WeatherElast_p75$Q_Protein)] <- tab_WeatherElast_p75$CI_protein[is.na(tab_WeatherElast_p75$Q_Protein)]
tab_WeatherElast_p75$CI_protein <- NULL

tab_WeatherElast_p75$Q_Oilseed[is.na(tab_WeatherElast_p75$Q_Oilseed)] <- tab_WeatherElast_p75$CI_oilseed[is.na(tab_WeatherElast_p75$Q_Oilseed)]
tab_WeatherElast_p75$CI_oilseed <- NULL

tab_WeatherElast_p75$Q_Roots[is.na(tab_WeatherElast_p75$Q_Roots)] <- tab_WeatherElast_p75$CI_roots[is.na(tab_WeatherElast_p75$Q_Roots)]
tab_WeatherElast_p75$CI_roots <- NULL

tab_WeatherElast_p75$Q_Corn[is.na(tab_WeatherElast_p75$Q_Corn)] <- tab_WeatherElast_p75$CI_corn[is.na(tab_WeatherElast_p75$Q_Corn)]
tab_WeatherElast_p75$CI_corn <- NULL

tab_WeatherElast_p75$X_Fertilizer[is.na(tab_WeatherElast_p75$X_Fertilizer)] <- tab_WeatherElast_p75$CI_fert[is.na(tab_WeatherElast_p75$X_Fertilizer)]
tab_WeatherElast_p75$CI_fert <- NULL

# Insert observed and past weather indicators
tab_WeatherElast_p75 <- berryFunctions::insertRows(tab_WeatherElast_p75, 1 , new = c("Observed weather","","","","","",""))
tab_WeatherElast_p75 <- berryFunctions::insertRows(tab_WeatherElast_p75, nrow(WeatherElast_cereals_p75)+2 , new = c("Past weather","","","","","",""))

# Print
write_xlsx(tab_structural,"Tables/Table_S5.6.xlsx")

