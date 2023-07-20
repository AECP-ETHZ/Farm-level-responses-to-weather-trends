#-----------------------------------------------------#
#                                                     #
# This program estimates the model with curvature im- #
#                                                     #
# Note: The bootstrap program (Section 1: Bootstrap   #
# program) cannot be run because farm-level data are  #
# confidential. The results of the bootstrap  are     #
# stored in lines 1716 and 1731, allowing to          # 
# create the tables in Section 2 (Presentation of     #
# results).                                           # 
#                                                     #
#-----------------------------------------------------#

# load packages
library(dplyr)
library(writexl)
library(nlsur)
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
                              gdd_obs_fm, prec_obs_fm, gddHigh_obs_fm, dd_obs_fm)


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
      iCorn_me <- as.numeric(do.call(c, list_probit_act[[5]][6]))
      
      #----------------------------------------#
      #### Step 2: Run Structural equations ####
      #----------------------------------------#
      
      #De-mean the data so the regression can converge
      vars <- c("qq_cereals", "qq_protein", "qq_oilseed", "qq_roots", "qq_corn", "x_fert", "nx_fert",
                "np_cereals", "np_protein", "np_oilseed", "np_roots", "np_corn", "nw_fert", 
                "k_land", "k_labor", "k_capital", "trend", "trend2",
                "gdd_obs", "prec_obs", "gddHigh_obs", "dd_obs", 
                "gdd_1to3", "prec_1to3", "gddHigh_1to3", "dd_1to3", 
                "gdd_4to10", "prec_4to10", "gddHigh_4to10", "dd_4to10", 
                "np_cereals_fm", "np_protein_fm", "np_oilseed_fm", "np_roots_fm", "np_corn_fm", "nw_fert_fm", 
                "k_land_fm", "k_labor_fm", "k_capital_fm", "trend_fm", "trend2_fm", 
                "gdd_obs_fm", "prec_obs_fm", "gddHigh_obs_fm", "dd_obs_fm", 
                "gdd_1to3_fm", "prec_1to3_fm", "gddHigh_1to3_fm", "dd_1to3_fm", 
                "gdd_4to10_fm", "prec_4to10_fm", "gddHigh_4to10_fm", "dd_4to10_fm", 
                "PHI_cereals", "phi_cereals", "PHI_protein", "phi_protein", "PHI_oilseed", "phi_oilseed", "PHI_roots", "phi_roots", "PHI_corn", "phi_corn")
      
      temp_fn <- function(x) x/mean(x) 
      data_boot_dem <- data_boot %>% mutate(across(vars, temp_fn))
      
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
      
      #calculate mean of the (original, i.e. non-demeaned) data for elasticities
      dat_sm <- dplyr::summarise_all(data_boot, mean) 
      
      ## Regression with nonlinear SUR 
      system.time(model_linear <- nlsur(eqns = eqns, data = data_boot_dem, startvalues = strtvls, type = "FGNLS",
                                        trace = TRUE, eps = 1e-03, tol = 1e-07, maxiter=50)) #Note: If it has not converged after 50 iters, go to the next bootstrap draw
      
      coef_linear_restr <- coef(model_linear)  
      
      # ---------------------------------------------------------------------------- #
      # Calculate "x times \beta" needed for elasticities that account for selection #
      # -----------------------------------------------------------------------------#
      
      qq_cereals_pred_b <- ( mean(model_linear$fitted$qq_cereals) - coef_linear_restr["phi11"]*mean(data_boot_dem$phi_cereals) )
      qq_cereals_pred <- (qq_cereals_pred_b * mean(data_boot$qq_cereals)) / mean(data_boot$PHI_cereals)
      
      qq_protein_pred_b <- mean(model_linear$fitted$qq_protein) - coef_linear_restr["phi22"]*mean(data_boot_dem$phi_protein)
      qq_protein_pred <- (qq_protein_pred_b * mean(data_boot$qq_protein)) / mean(data_boot$PHI_protein)
      
      qq_oilseed_pred_b <- mean(model_linear$fitted$qq_oilseed) - coef_linear_restr["phi33"]*mean(data_boot_dem$phi_oilseed)
      qq_oilseed_pred <- (qq_oilseed_pred_b * mean(data_boot$qq_oilseed)) / mean(data_boot$PHI_oilseed)
      
      qq_roots_pred_b <- mean(model_linear$fitted$qq_roots) - coef_linear_restr["phi44"]*mean(data_boot_dem$phi_roots)
      qq_roots_pred <- (qq_roots_pred_b * mean(data_boot$qq_roots)) / mean(data_boot$PHI_roots)
      
      qq_corn_pred_b <- mean(model_linear$fitted$qq_corn) - coef_linear_restr["phi55"]*mean(data_boot_dem$phi_corn)
      qq_corn_pred <- (qq_corn_pred_b * mean(data_boot$qq_corn)) / mean(data_boot$PHI_corn)
      
      nx_fert_pred <- mean(model_linear$fitted$nx_fert) * mean(data_boot$x_fert)
      
      
      #----------------------------------------------------------#
      ##### Recover original (i.e., non-demeaned) parameters #####
      #----------------------------------------------------------#  
      
      # Recover price coefficients
      b11 <- model_linear$coefficients["c11"]*model_linear$coefficients["c11"] * mean(data_boot$qq_cereals) / mean(data_boot$np_cereals) / mean(data_boot$PHI_cereals)
      b12 <- model_linear$coefficients["c11"]*model_linear$coefficients["c12"] * mean(data_boot$qq_cereals) / mean(data_boot$np_protein) / mean(data_boot$PHI_cereals)
      b13 <- model_linear$coefficients["c11"]*model_linear$coefficients["c13"] * mean(data_boot$qq_cereals) / mean(data_boot$np_oilseed) / mean(data_boot$PHI_cereals)
      b14 <- model_linear$coefficients["c11"]*model_linear$coefficients["c14"] * mean(data_boot$qq_cereals) / mean(data_boot$np_roots) / mean(data_boot$PHI_cereals)
      b15 <- model_linear$coefficients["c11"]*model_linear$coefficients["c15"] * mean(data_boot$qq_cereals) / mean(data_boot$np_corn) / mean(data_boot$PHI_cereals)
      b16 <- model_linear$coefficients["c11"]*model_linear$coefficients["c16"] * mean(data_boot$qq_cereals) / mean(data_boot$nw_fert) / mean(data_boot$PHI_cereals)
      
      b21 <- (model_linear$coefficients["c12"]*model_linear$coefficients["c11"]) * mean(data_boot$qq_protein) / mean(data_boot$np_cereals) / mean(data_boot$PHI_protein)
      b22 <- (model_linear$coefficients["c12"]*model_linear$coefficients["c12"]+model_linear$coefficients["c22"]*model_linear$coefficients["c22"]) * mean(data_boot$qq_protein) / mean(data_boot$np_protein) / mean(data_boot$PHI_protein)
      b23 <- (model_linear$coefficients["c12"]*model_linear$coefficients["c13"]+model_linear$coefficients["c22"]*model_linear$coefficients["c23"]) * mean(data_boot$qq_protein) / mean(data_boot$np_oilseed) / mean(data_boot$PHI_protein)
      b24 <- (model_linear$coefficients["c12"]*model_linear$coefficients["c14"]+model_linear$coefficients["c22"]*model_linear$coefficients["c24"]) * mean(data_boot$qq_protein) / mean(data_boot$np_roots) / mean(data_boot$PHI_protein)
      b25 <- (model_linear$coefficients["c12"]*model_linear$coefficients["c15"]+model_linear$coefficients["c22"]*model_linear$coefficients["c25"]) * mean(data_boot$qq_protein) / mean(data_boot$np_corn) / mean(data_boot$PHI_protein)
      b26 <- (model_linear$coefficients["c12"]*model_linear$coefficients["c16"]+model_linear$coefficients["c22"]*model_linear$coefficients["c26"]) * mean(data_boot$qq_protein) / mean(data_boot$nw_fert) / mean(data_boot$PHI_protein)
      
      b31 <- (model_linear$coefficients["c13"]*model_linear$coefficients["c11"]) * mean(data_boot$qq_oilseed) / mean(data_boot$np_cereals) / mean(data_boot$PHI_oilseed)
      b32 <- (model_linear$coefficients["c13"]*model_linear$coefficients["c12"]+model_linear$coefficients["c23"]*model_linear$coefficients["c22"]) * mean(data_boot$qq_oilseed) / mean(data_boot$np_protein) / mean(data_boot$PHI_oilseed)
      b33 <- (model_linear$coefficients["c13"]*model_linear$coefficients["c13"]+model_linear$coefficients["c23"]*model_linear$coefficients["c23"]+model_linear$coefficients["c33"]*model_linear$coefficients["c33"]) * mean(data_boot$qq_oilseed) / mean(data_boot$np_oilseed) / mean(data_boot$PHI_oilseed)
      b34 <- (model_linear$coefficients["c13"]*model_linear$coefficients["c14"]+model_linear$coefficients["c23"]*model_linear$coefficients["c24"]+model_linear$coefficients["c33"]*model_linear$coefficients["c34"]) * mean(data_boot$qq_oilseed) / mean(data_boot$np_roots) / mean(data_boot$PHI_oilseed)
      b35 <- (model_linear$coefficients["c13"]*model_linear$coefficients["c15"]+model_linear$coefficients["c23"]*model_linear$coefficients["c25"]+model_linear$coefficients["c33"]*model_linear$coefficients["c35"]) * mean(data_boot$qq_oilseed) / mean(data_boot$np_corn) / mean(data_boot$PHI_oilseed)
      b36 <- (model_linear$coefficients["c13"]*model_linear$coefficients["c16"]+model_linear$coefficients["c23"]*model_linear$coefficients["c26"]+model_linear$coefficients["c33"]*model_linear$coefficients["c36"]) * mean(data_boot$qq_oilseed) / mean(data_boot$nw_fert) / mean(data_boot$PHI_oilseed)
      
      b41 <- (model_linear$coefficients["c14"]*model_linear$coefficients["c11"]) * mean(data_boot$qq_roots) / mean(data_boot$np_cereals) / mean(data_boot$PHI_roots)
      b42 <- (model_linear$coefficients["c14"]*model_linear$coefficients["c12"]+model_linear$coefficients["c24"]*model_linear$coefficients["c22"]) * mean(data_boot$qq_roots) / mean(data_boot$np_protein) / mean(data_boot$PHI_roots)
      b43 <- (model_linear$coefficients["c14"]*model_linear$coefficients["c13"]+model_linear$coefficients["c24"]*model_linear$coefficients["c23"]+model_linear$coefficients["c34"]*model_linear$coefficients["c33"]) * mean(data_boot$qq_roots) / mean(data_boot$np_oilseed) / mean(data_boot$PHI_roots)
      b44 <- (model_linear$coefficients["c14"]*model_linear$coefficients["c14"]+model_linear$coefficients["c24"]*model_linear$coefficients["c24"]+model_linear$coefficients["c34"]*model_linear$coefficients["c34"]+model_linear$coefficients["c44"]*model_linear$coefficients["c44"]) * mean(data_boot$qq_roots) / mean(data_boot$np_roots) / mean(data_boot$PHI_roots)
      b45 <- (model_linear$coefficients["c14"]*model_linear$coefficients["c15"]+model_linear$coefficients["c24"]*model_linear$coefficients["c25"]+model_linear$coefficients["c34"]*model_linear$coefficients["c35"]+model_linear$coefficients["c44"]*model_linear$coefficients["c45"]) * mean(data_boot$qq_roots) / mean(data_boot$np_corn) / mean(data_boot$PHI_roots)
      b46 <- (model_linear$coefficients["c14"]*model_linear$coefficients["c16"]+model_linear$coefficients["c24"]*model_linear$coefficients["c26"]+model_linear$coefficients["c34"]*model_linear$coefficients["c36"]+model_linear$coefficients["c44"]*model_linear$coefficients["c46"]) * mean(data_boot$qq_roots) / mean(data_boot$nw_fert) / mean(data_boot$PHI_roots)
      
      b51 <- (model_linear$coefficients["c15"]*model_linear$coefficients["c11"]) * mean(data_boot$qq_corn) / mean(data_boot$np_cereals) / mean(data_boot$PHI_corn)
      b52 <- (model_linear$coefficients["c15"]*model_linear$coefficients["c12"]+model_linear$coefficients["c25"]*model_linear$coefficients["c22"]) * mean(data_boot$qq_corn) / mean(data_boot$np_protein) / mean(data_boot$PHI_corn)
      b53 <- (model_linear$coefficients["c15"]*model_linear$coefficients["c13"]+model_linear$coefficients["c25"]*model_linear$coefficients["c23"]+model_linear$coefficients["c35"]*model_linear$coefficients["c33"]) * mean(data_boot$qq_corn) / mean(data_boot$np_oilseed) / mean(data_boot$PHI_corn)
      b54 <- (model_linear$coefficients["c15"]*model_linear$coefficients["c14"]+model_linear$coefficients["c25"]*model_linear$coefficients["c24"]+model_linear$coefficients["c35"]*model_linear$coefficients["c34"]+model_linear$coefficients["c45"]*model_linear$coefficients["c44"]) * mean(data_boot$qq_corn) / mean(data_boot$np_roots) / mean(data_boot$PHI_corn)
      b55 <- (model_linear$coefficients["c15"]*model_linear$coefficients["c15"]+model_linear$coefficients["c25"]*model_linear$coefficients["c25"]+model_linear$coefficients["c35"]*model_linear$coefficients["c35"]+model_linear$coefficients["c45"]*model_linear$coefficients["c45"]) * mean(data_boot$qq_corn) / mean(data_boot$np_corn) / mean(data_boot$PHI_corn)
      b56 <- (model_linear$coefficients["c15"]*model_linear$coefficients["c16"]+model_linear$coefficients["c25"]*model_linear$coefficients["c26"]+model_linear$coefficients["c35"]*model_linear$coefficients["c36"]+model_linear$coefficients["c45"]*model_linear$coefficients["c46"]) * mean(data_boot$qq_corn) / mean(data_boot$nw_fert) / mean(data_boot$PHI_corn)
      
      b61 <- (model_linear$coefficients["c16"]*model_linear$coefficients["c11"]) * mean(data_boot$x_fert) / mean(data_boot$np_cereals)
      b62 <- (model_linear$coefficients["c16"]*model_linear$coefficients["c12"]+model_linear$coefficients["c26"]*model_linear$coefficients["c22"]) * mean(data_boot$x_fert) / mean(data_boot$np_protein)
      b63 <- (model_linear$coefficients["c16"]*model_linear$coefficients["c13"]+model_linear$coefficients["c26"]*model_linear$coefficients["c23"]+model_linear$coefficients["c36"]*model_linear$coefficients["c33"]) * mean(data_boot$x_fert) / mean(data_boot$np_oilseed)
      b64 <- (model_linear$coefficients["c16"]*model_linear$coefficients["c14"]+model_linear$coefficients["c26"]*model_linear$coefficients["c24"]+model_linear$coefficients["c36"]*model_linear$coefficients["c34"]+model_linear$coefficients["c46"]*model_linear$coefficients["c44"]) * mean(data_boot$x_fert) / mean(data_boot$np_roots)
      b65 <- (model_linear$coefficients["c16"]*model_linear$coefficients["c15"]+model_linear$coefficients["c26"]*model_linear$coefficients["c25"]+model_linear$coefficients["c36"]*model_linear$coefficients["c35"]+model_linear$coefficients["c46"]*model_linear$coefficients["c45"]) * mean(data_boot$x_fert) / mean(data_boot$np_corn)
      b66 <- (model_linear$coefficients["c16"]*model_linear$coefficients["c16"]+model_linear$coefficients["c26"]*model_linear$coefficients["c26"]+model_linear$coefficients["c36"]*model_linear$coefficients["c36"]+model_linear$coefficients["c46"]*model_linear$coefficients["c46"]) * mean(data_boot$x_fert) / mean(data_boot$nw_fert)
      
      # Recover weather coefficients
      w11 <- model_linear$coefficients["w11"] * mean(data_boot$qq_cereals) / mean(data_boot$gdd_obs) / mean(data_boot$PHI_cereals)
      w12 <- model_linear$coefficients["w12"] * mean(data_boot$qq_cereals) / mean(data_boot$prec_obs) / mean(data_boot$PHI_cereals)
      w13 <- model_linear$coefficients["w13"] * mean(data_boot$qq_cereals) / mean(data_boot$gddHigh_obs) / mean(data_boot$PHI_cereals)
      w14 <- model_linear$coefficients["w14"] * mean(data_boot$qq_cereals) / mean(data_boot$dd_obs) / mean(data_boot$PHI_cereals)
      w15 <- model_linear$coefficients["w15"] * mean(data_boot$qq_cereals) / mean(data_boot$gdd_1to3) / mean(data_boot$PHI_cereals)
      w16 <- model_linear$coefficients["w16"] * mean(data_boot$qq_cereals) / mean(data_boot$prec_1to3) / mean(data_boot$PHI_cereals)
      w17 <- model_linear$coefficients["w17"] * mean(data_boot$qq_cereals) / mean(data_boot$gddHigh_1to3) / mean(data_boot$PHI_cereals)
      w18 <- model_linear$coefficients["w18"] * mean(data_boot$qq_cereals) / mean(data_boot$dd_1to3) / mean(data_boot$PHI_cereals)
      w19 <- model_linear$coefficients["w19"] * mean(data_boot$qq_cereals) / mean(data_boot$gdd_4to10) / mean(data_boot$PHI_cereals)
      w110 <- model_linear$coefficients["w110"] * mean(data_boot$qq_cereals) / mean(data_boot$prec_4to10) / mean(data_boot$PHI_cereals)
      w111 <- model_linear$coefficients["w111"] * mean(data_boot$qq_cereals) / mean(data_boot$gddHigh_4to10) / mean(data_boot$PHI_cereals)
      w112 <- model_linear$coefficients["w112"] * mean(data_boot$qq_cereals) / mean(data_boot$dd_4to10) / mean(data_boot$PHI_cereals) #-2.80314
      
      w21 <- model_linear$coefficients["w21"] * mean(data_boot$qq_protein) / mean(data_boot$gdd_obs) / mean(data_boot$PHI_protein)
      w22 <- model_linear$coefficients["w22"] * mean(data_boot$qq_protein) / mean(data_boot$prec_obs) / mean(data_boot$PHI_protein)
      w23 <- model_linear$coefficients["w23"] * mean(data_boot$qq_protein) / mean(data_boot$gddHigh_obs) / mean(data_boot$PHI_protein)
      w24 <- model_linear$coefficients["w24"] * mean(data_boot$qq_protein) / mean(data_boot$dd_obs) / mean(data_boot$PHI_protein)
      w25 <- model_linear$coefficients["w25"] * mean(data_boot$qq_protein) / mean(data_boot$gdd_1to3) / mean(data_boot$PHI_protein)
      w26 <- model_linear$coefficients["w26"] * mean(data_boot$qq_protein) / mean(data_boot$prec_1to3) / mean(data_boot$PHI_protein)
      w27 <- model_linear$coefficients["w27"] * mean(data_boot$qq_protein) / mean(data_boot$gddHigh_1to3) / mean(data_boot$PHI_protein)
      w28 <- model_linear$coefficients["w28"] * mean(data_boot$qq_protein) / mean(data_boot$dd_1to3) / mean(data_boot$PHI_protein)
      w29 <- model_linear$coefficients["w29"] * mean(data_boot$qq_protein) / mean(data_boot$gdd_4to10) / mean(data_boot$PHI_protein)
      w210 <- model_linear$coefficients["w210"] * mean(data_boot$qq_protein) / mean(data_boot$prec_4to10) / mean(data_boot$PHI_protein)
      w211 <- model_linear$coefficients["w211"] * mean(data_boot$qq_protein) / mean(data_boot$gddHigh_4to10) / mean(data_boot$PHI_protein)
      w212 <- model_linear$coefficients["w212"] * mean(data_boot$qq_protein) / mean(data_boot$dd_4to10) / mean(data_boot$PHI_protein)
      
      w31 <- model_linear$coefficients["w31"] * mean(data_boot$qq_oilseed) / mean(data_boot$gdd_obs) / mean(data_boot$PHI_oilseed)
      w32 <- model_linear$coefficients["w32"] * mean(data_boot$qq_oilseed) / mean(data_boot$prec_obs) / mean(data_boot$PHI_oilseed)
      w33 <- model_linear$coefficients["w33"] * mean(data_boot$qq_oilseed) / mean(data_boot$gddHigh_obs) / mean(data_boot$PHI_oilseed)
      w34 <- model_linear$coefficients["w34"] * mean(data_boot$qq_oilseed) / mean(data_boot$dd_obs) / mean(data_boot$PHI_oilseed)
      w35 <- model_linear$coefficients["w35"] * mean(data_boot$qq_oilseed) / mean(data_boot$gdd_1to3) / mean(data_boot$PHI_oilseed)
      w36 <- model_linear$coefficients["w36"] * mean(data_boot$qq_oilseed) / mean(data_boot$prec_1to3) / mean(data_boot$PHI_oilseed)
      w37 <- model_linear$coefficients["w37"] * mean(data_boot$qq_oilseed) / mean(data_boot$gddHigh_1to3) / mean(data_boot$PHI_oilseed)
      w38 <- model_linear$coefficients["w38"] * mean(data_boot$qq_oilseed) / mean(data_boot$dd_1to3) / mean(data_boot$PHI_oilseed)
      w39 <- model_linear$coefficients["w39"] * mean(data_boot$qq_oilseed) / mean(data_boot$gdd_4to10) / mean(data_boot$PHI_oilseed)
      w310 <- model_linear$coefficients["w310"] * mean(data_boot$qq_oilseed) / mean(data_boot$prec_4to10) / mean(data_boot$PHI_oilseed)
      w311 <- model_linear$coefficients["w311"] * mean(data_boot$qq_oilseed) / mean(data_boot$gddHigh_4to10) / mean(data_boot$PHI_oilseed)
      w312 <- model_linear$coefficients["w312"] * mean(data_boot$qq_oilseed) / mean(data_boot$dd_4to10) / mean(data_boot$PHI_oilseed)
      
      w41 <- model_linear$coefficients["w41"] * mean(data_boot$qq_roots) / mean(data_boot$gdd_obs) / mean(data_boot$PHI_roots)
      w42 <- model_linear$coefficients["w42"] * mean(data_boot$qq_roots) / mean(data_boot$prec_obs) / mean(data_boot$PHI_roots)
      w43 <- model_linear$coefficients["w43"] * mean(data_boot$qq_roots) / mean(data_boot$gddHigh_obs) / mean(data_boot$PHI_roots)
      w44 <- model_linear$coefficients["w44"] * mean(data_boot$qq_roots) / mean(data_boot$dd_obs) / mean(data_boot$PHI_roots)
      w45 <- model_linear$coefficients["w45"] * mean(data_boot$qq_roots) / mean(data_boot$gdd_1to3) / mean(data_boot$PHI_roots)
      w46 <- model_linear$coefficients["w46"] * mean(data_boot$qq_roots) / mean(data_boot$prec_1to3) / mean(data_boot$PHI_roots)
      w47 <- model_linear$coefficients["w47"] * mean(data_boot$qq_roots) / mean(data_boot$gddHigh_1to3) / mean(data_boot$PHI_roots)
      w48 <- model_linear$coefficients["w48"] * mean(data_boot$qq_roots) / mean(data_boot$dd_1to3) / mean(data_boot$PHI_roots)
      w49 <- model_linear$coefficients["w49"] * mean(data_boot$qq_roots) / mean(data_boot$gdd_4to10) / mean(data_boot$PHI_roots)
      w410 <- model_linear$coefficients["w410"] * mean(data_boot$qq_roots) / mean(data_boot$prec_4to10) / mean(data_boot$PHI_roots)
      w411 <- model_linear$coefficients["w411"] * mean(data_boot$qq_roots) / mean(data_boot$gddHigh_4to10) / mean(data_boot$PHI_roots)
      w412 <- model_linear$coefficients["w412"] * mean(data_boot$qq_roots) / mean(data_boot$dd_4to10) / mean(data_boot$PHI_roots)
      
      w51 <- model_linear$coefficients["w51"] * mean(data_boot$qq_corn) / mean(data_boot$gdd_obs) / mean(data_boot$PHI_corn)
      w52 <- model_linear$coefficients["w52"] * mean(data_boot$qq_corn) / mean(data_boot$prec_obs) / mean(data_boot$PHI_corn)
      w53 <- model_linear$coefficients["w53"] * mean(data_boot$qq_corn) / mean(data_boot$gddHigh_obs) / mean(data_boot$PHI_corn)
      w54 <- model_linear$coefficients["w54"] * mean(data_boot$qq_corn) / mean(data_boot$dd_obs) / mean(data_boot$PHI_corn)
      w55 <- model_linear$coefficients["w55"] * mean(data_boot$qq_corn) / mean(data_boot$gdd_1to3) / mean(data_boot$PHI_corn)
      w56 <- model_linear$coefficients["w56"] * mean(data_boot$qq_corn) / mean(data_boot$prec_1to3) / mean(data_boot$PHI_corn)
      w57 <- model_linear$coefficients["w57"] * mean(data_boot$qq_corn) / mean(data_boot$gddHigh_1to3) / mean(data_boot$PHI_corn)
      w58 <- model_linear$coefficients["w58"] * mean(data_boot$qq_corn) / mean(data_boot$dd_1to3) / mean(data_boot$PHI_corn)
      w59 <- model_linear$coefficients["w59"] * mean(data_boot$qq_corn) / mean(data_boot$gdd_4to10) / mean(data_boot$PHI_corn)
      w510 <- model_linear$coefficients["w510"] * mean(data_boot$qq_corn) / mean(data_boot$prec_4to10) / mean(data_boot$PHI_corn)
      w511 <- model_linear$coefficients["w511"] * mean(data_boot$qq_corn) / mean(data_boot$gddHigh_4to10) / mean(data_boot$PHI_corn)
      w512 <- model_linear$coefficients["w512"] * mean(data_boot$qq_corn) / mean(data_boot$dd_4to10) / mean(data_boot$PHI_corn)
      
      w61 <- model_linear$coefficients["w61"] * mean(data_boot$x_fert) / mean(data_boot$gdd_obs) 
      w62 <- model_linear$coefficients["w62"] * mean(data_boot$x_fert) / mean(data_boot$prec_obs) 
      w63 <- model_linear$coefficients["w63"] * mean(data_boot$x_fert) / mean(data_boot$gddHigh_obs) 
      w64 <- model_linear$coefficients["w64"] * mean(data_boot$x_fert) / mean(data_boot$dd_obs) 
      w65 <- model_linear$coefficients["w65"] * mean(data_boot$x_fert) / mean(data_boot$gdd_1to3) 
      w66 <- model_linear$coefficients["w66"] * mean(data_boot$x_fert) / mean(data_boot$prec_1to3) 
      w67 <- model_linear$coefficients["w67"] * mean(data_boot$x_fert) / mean(data_boot$gddHigh_1to3) 
      w68 <- model_linear$coefficients["w68"] * mean(data_boot$x_fert) / mean(data_boot$dd_1to3) 
      w69 <- model_linear$coefficients["w69"] * mean(data_boot$x_fert) / mean(data_boot$gdd_4to10) 
      w610 <- model_linear$coefficients["w610"] * mean(data_boot$x_fert) / mean(data_boot$prec_4to10) 
      w611 <- model_linear$coefficients["w611"] * mean(data_boot$x_fert) / mean(data_boot$gddHigh_4to10) 
      w612 <- model_linear$coefficients["w612"] * mean(data_boot$x_fert) / mean(data_boot$dd_4to10)
      
      # recover phi_* coefficients
      phi1 <- model_linear$coefficients["phi11"] * mean(data_boot$qq_cereals) / mean(data_boot$phi_cereals)
      phi2 <- model_linear$coefficients["phi22"] * mean(data_boot$qq_protein) / mean(data_boot$phi_protein)
      phi3 <- model_linear$coefficients["phi33"] * mean(data_boot$qq_oilseed) / mean(data_boot$phi_oilseed)
      phi4 <- model_linear$coefficients["phi44"] * mean(data_boot$qq_roots) / mean(data_boot$phi_roots)
      phi5 <- model_linear$coefficients["phi55"] * mean(data_boot$qq_corn) / mean(data_boot$phi_corn)
      
      # Create an object coef_linear with the same names as in the unrestricted model
      
        # Load structural coefficients from unrestricted model so I can replace the coefficients there
        load("rOutput/coef_structural_linear.Rda")
        
        coef_linear["QQcereals_I(PHI_cereals * np_cereals)"] <- b11
        coef_linear["QQcereals_I(PHI_cereals * np_protein)"] <- b12
        coef_linear["QQcereals_I(PHI_cereals * np_oilseed)"] <- b13
        coef_linear["QQcereals_I(PHI_cereals * np_roots)"] <- b14
        coef_linear["QQcereals_I(PHI_cereals * np_corn)"] <- b15
        coef_linear["QQcereals_I(PHI_cereals * nw_fert)"] <- b16
        coef_linear["QQcereals_I(PHI_cereals * gdd_obs)"] <- w11
        coef_linear["QQcereals_I(PHI_cereals * prec_obs)"] <- w12
        coef_linear["QQcereals_I(PHI_cereals * gddHigh_obs)"] <- w13
        coef_linear["QQcereals_I(PHI_cereals * dd_obs)"] <- w14
        coef_linear["QQcereals_I(PHI_cereals * gdd_1to3)"] <- w15
        coef_linear["QQcereals_I(PHI_cereals * prec_1to3)"] <- w16
        coef_linear["QQcereals_I(PHI_cereals * gddHigh_1to3)"] <- w17
        coef_linear["QQcereals_I(PHI_cereals * dd_1to3)"] <- w18
        coef_linear["QQcereals_I(PHI_cereals * gdd_4to10)"] <- w19
        coef_linear["QQcereals_I(PHI_cereals * prec_4to10)"] <- w110
        coef_linear["QQcereals_I(PHI_cereals * gddHigh_4to10)"] <- w111
        coef_linear["QQcereals_I(PHI_cereals * dd_4to10)"] <- w112
        coef_linear["QQcereals_phi_cereals"] <- phi1
        
        coef_linear["QQprotein_I(PHI_protein * np_cereals)"] <- b21
        coef_linear["QQprotein_I(PHI_protein * np_protein)"] <- b22
        coef_linear["QQprotein_I(PHI_protein * np_oilseed)"] <- b23
        coef_linear["QQprotein_I(PHI_protein * np_roots)"] <- b24
        coef_linear["QQprotein_I(PHI_protein * np_corn)"] <- b25
        coef_linear["QQprotein_I(PHI_protein * nw_fert)"] <- b26
        coef_linear["QQprotein_I(PHI_protein * gdd_obs)"] <- w21
        coef_linear["QQprotein_I(PHI_protein * prec_obs)"] <- w22
        coef_linear["QQprotein_I(PHI_protein * gddHigh_obs)"] <- w23
        coef_linear["QQprotein_I(PHI_protein * dd_obs)"] <- w24
        coef_linear["QQprotein_I(PHI_protein * gdd_1to3)"] <- w25
        coef_linear["QQprotein_I(PHI_protein * prec_1to3)"] <- w26
        coef_linear["QQprotein_I(PHI_protein * gddHigh_1to3)"] <- w27
        coef_linear["QQprotein_I(PHI_protein * dd_1to3)"] <- w28
        coef_linear["QQprotein_I(PHI_protein * gdd_4to10)"] <- w29
        coef_linear["QQprotein_I(PHI_protein * prec_4to10)"] <- w210
        coef_linear["QQprotein_I(PHI_protein * gddHigh_4to10)"] <- w211
        coef_linear["QQprotein_I(PHI_protein * dd_4to10)"] <- w212
        coef_linear["QQprotein_phi_protein"] <- phi2
        
        coef_linear["QQoilseed_I(PHI_oilseed * np_cereals)"] <- b31
        coef_linear["QQoilseed_I(PHI_oilseed * np_protein)"] <- b32
        coef_linear["QQoilseed_I(PHI_oilseed * np_oilseed)"] <- b33
        coef_linear["QQoilseed_I(PHI_oilseed * np_roots)"] <- b34
        coef_linear["QQoilseed_I(PHI_oilseed * np_corn)"] <- b35
        coef_linear["QQoilseed_I(PHI_oilseed * nw_fert)"] <- b36
        coef_linear["QQoilseed_I(PHI_oilseed * gdd_obs)"] <- w31
        coef_linear["QQoilseed_I(PHI_oilseed * prec_obs)"] <- w32
        coef_linear["QQoilseed_I(PHI_oilseed * gddHigh_obs)"] <- w33
        coef_linear["QQoilseed_I(PHI_oilseed * dd_obs)"] <- w34
        coef_linear["QQoilseed_I(PHI_oilseed * gdd_1to3)"] <- w35
        coef_linear["QQoilseed_I(PHI_oilseed * prec_1to3)"] <- w36
        coef_linear["QQoilseed_I(PHI_oilseed * gddHigh_1to3)"] <- w37
        coef_linear["QQoilseed_I(PHI_oilseed * dd_1to3)"] <- w38
        coef_linear["QQoilseed_I(PHI_oilseed * gdd_4to10)"] <- w39
        coef_linear["QQoilseed_I(PHI_oilseed * prec_4to10)"] <- w310
        coef_linear["QQoilseed_I(PHI_oilseed * gddHigh_4to10)"] <- w311
        coef_linear["QQoilseed_I(PHI_oilseed * dd_4to10)"] <- w312
        coef_linear["QQoilseed_phi_oilseed"] <- phi3
        
        coef_linear["QQroots_I(PHI_roots * np_cereals)"] <- b41
        coef_linear["QQroots_I(PHI_roots * np_protein)"] <- b42
        coef_linear["QQroots_I(PHI_roots * np_oilseed)"] <- b43
        coef_linear["QQroots_I(PHI_roots * np_roots)"] <- b44
        coef_linear["QQroots_I(PHI_roots * np_corn)"] <- b45
        coef_linear["QQroots_I(PHI_roots * nw_fert)"] <- b46
        coef_linear["QQroots_I(PHI_roots * gdd_obs)"] <- w41
        coef_linear["QQroots_I(PHI_roots * prec_obs)"] <- w42
        coef_linear["QQroots_I(PHI_roots * gddHigh_obs)"] <- w43
        coef_linear["QQroots_I(PHI_roots * dd_obs)"] <- w44
        coef_linear["QQroots_I(PHI_roots * gdd_1to3)"] <- w45
        coef_linear["QQroots_I(PHI_roots * prec_1to3)"] <- w46
        coef_linear["QQroots_I(PHI_roots * gddHigh_1to3)"] <- w47
        coef_linear["QQroots_I(PHI_roots * dd_1to3)"] <- w48
        coef_linear["QQroots_I(PHI_roots * gdd_4to10)"] <- w49
        coef_linear["QQroots_I(PHI_roots * prec_4to10)"] <- w410
        coef_linear["QQroots_I(PHI_roots * gddHigh_4to10)"] <- w411
        coef_linear["QQroots_I(PHI_roots * dd_4to10)"] <- w412
        coef_linear["QQroots_phi_roots"] <- phi4
        
        coef_linear["QQcorn_I(PHI_corn * np_cereals)"] <- b51
        coef_linear["QQcorn_I(PHI_corn * np_protein)"] <- b52
        coef_linear["QQcorn_I(PHI_corn * np_oilseed)"] <- b53
        coef_linear["QQcorn_I(PHI_corn * np_roots)"] <- b54
        coef_linear["QQcorn_I(PHI_corn * np_corn)"] <- b55
        coef_linear["QQcorn_I(PHI_corn * nw_fert)"] <- b56
        coef_linear["QQcorn_I(PHI_corn * gdd_obs)"] <- w51
        coef_linear["QQcorn_I(PHI_corn * prec_obs)"] <- w52
        coef_linear["QQcorn_I(PHI_corn * gddHigh_obs)"] <- w53
        coef_linear["QQcorn_I(PHI_corn * dd_obs)"] <- w54
        coef_linear["QQcorn_I(PHI_corn * gdd_1to3)"] <- w55
        coef_linear["QQcorn_I(PHI_corn * prec_1to3)"] <- w56
        coef_linear["QQcorn_I(PHI_corn * gddHigh_1to3)"] <- w57
        coef_linear["QQcorn_I(PHI_corn * dd_1to3)"] <- w58
        coef_linear["QQcorn_I(PHI_corn * gdd_4to10)"] <- w59
        coef_linear["QQcorn_I(PHI_corn * prec_4to10)"] <-w510
        coef_linear["QQcorn_I(PHI_corn * gddHigh_4to10)"] <-w511
        coef_linear["QQcorn_I(PHI_corn * dd_4to10)"] <-w512
        coef_linear["QQcorn_phi_corn"] <- phi5
        
        coef_linear["NXfert_np_cereals"] <- b61
        coef_linear["NXfert_np_protein"] <- b62
        coef_linear["NXfert_np_oilseed"] <- b63
        coef_linear["NXfert_np_roots"] <- b64
        coef_linear["NXfert_np_corn"] <- b65
        coef_linear["NXfert_nw_fert"] <- b66
        coef_linear["NXfert_gdd_obs"] <- w61
        coef_linear["NXfert_prec_obs"] <- w62
        coef_linear["NXfert_gddHigh_obs"] <- w63
        coef_linear["NXfert_dd_obs"] <- w64
        coef_linear["NXfert_gdd_1to3"] <- w65
        coef_linear["NXfert_prec_1to3"] <- w66
        coef_linear["NXfert_gddHigh_1to3"] <- w67
        coef_linear["NXfert_dd_1to3"] <- w68
        coef_linear["NXfert_gdd_4to10"] <- w69
        coef_linear["NXfert_prec_4to10"] <- w610
        coef_linear["NXfert_gddHigh_4to10"] <- w611
        coef_linear["NXfert_dd_4to10"] <- w612
      
      
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
      
      result <- c(model_linear$coefficients, #structural coefficients
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
      
    }, error=function(err) {rep(NA,366)} ) # insert NA's if there is an optimization error
  }
  
}


#-------------------------------------------#
#### Obtain and store the actual results ####
#-------------------------------------------#

# run one draw with the original data
system.time(boot_actresults <- my.boot(data=df_farm, nrep=4, cluster=df_farm$key, nCores=1, actual=TRUE))

# save actual result
save(boot_actresults, file="rOutput/bootresults_restr_actual.Rda")
boot_actresults <- NULL

#------------------------------------------------#
#### Run bootstrap: Cluster at the farm-level ####
#------------------------------------------------#

# find number of cores
nCores <- parallel::detectCores() -1 

# run the bootstrap
set.seed(1234)
system.time(boot_results <- my.boot(data=df_farm, nrep=8, cluster=df_farm$key, nCores=nCores, actual=FALSE))

# save
save(boot_results, file="rOutput/bootresults_restr.Rda")
boot_results <- NULL



#----------------------------------#
#### 2) Presentation of results ####
#----------------------------------#

# Load bootstrap results
load("rOutput/bootresults_restr_actual.Rda")
load("rOutput/bootresults_restr.Rda")

# Turn matrix into data frame
boot_results <- as.data.frame(boot_results)
summary(boot_results$sel_qcer_gdd_obs)
table(!is.na(boot_results$sel_qcer_gdd_obs))

# --------------------------------- #
# Create list with original results #
# --------------------------------- #

names(boot_results)
boot_summary <- data.frame("Variable"=c(names(boot_results)),
                           "Original"=boot_actresults) 

# Add confidence intervals
CIs <- t(sapply(boot_results,function(i) quantile(i,c(0.025, 0.975), na.rm=TRUE)))
CIs <- as.data.frame(CIs)
colnames(CIs) <- c("L95","U95")

# Add CIs to boot_summary
boot_summary <- cbind(boot_summary, CIs)

# Organize
qqnxStructural <- boot_summary[c(1:269),]

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


#----------------------------------------------------#
#### Table: Parameter estimates for probit models ####
#----------------------------------------------------#

# (same as structural model) 

#------------------------------------------------------#
#### Table: Regression results structural equations ####
#------------------------------------------------------#

#Take original tables without Mundlak terms
tab_structural <- qqnxStructural %>% 
  filter(!grepl("_fm",Variable) &
           !startsWith(x = qqnxStructural$Variable, prefix = "m")) %>% 
  select("Variable","Coeff"="Original", "L95", "U95")

#Put into nice order
tab_structural$order <- seq(1:131)
tab_structural <- tab_structural[c(2:7, #cholesky
                                   9:13,
                                   15:18,
                                   20:22,
                                   1, #cereals
                                   25:27,
                                   43,
                                   49,
                                   55:66,
                                   127,
                                   8, #protein
                                   28:30,
                                   44,
                                   50,
                                   67:78,
                                   128,
                                   14, #oilseed
                                   31:33,
                                   45,
                                   51,
                                   79:90,
                                   129,
                                   19, #roots
                                   34:36,
                                   46,
                                   52,
                                   91:102,
                                   130,
                                   23, #corn
                                   37:39,
                                   47,
                                   53,
                                   103:114,
                                   131,
                                   24, #fertilizer
                                   40:42,
                                   48,
                                   54,
                                   115:126),]

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
tab_structural$Variable <- c("c11","c12","c13","c14","c15","c16",
                             "c22", "c23", "c24", "c25", "c26",
                             "c33", "c34", "c35", "c36",
                             "c44", "c45", "c46", 
                             rep(c("Intercept",
                                   "K Land",
                                   "K Labor",
                                   "K Capital",
                                   "Trend",
                                   "Trend^2",
                                   "GDD",
                                   "Prec",
                                   "GddHigh",
                                   "DryDays",
                                   "GDD1to3",
                                   "Prec1to3",
                                   "GddHigh1to3",
                                   "DryDays1to3",
                                   "GDD4to10",
                                   "Prec4to10",
                                   "GddHigh4to10",
                                   "DryDays4to10",
                                   "Phi"),5),
                             c("Intercept",
                               "K Land",
                               "K Labor",
                               "K Capital",
                               "Trend",
                               "Trend^2",
                               "GDD",
                               "Prec",
                               "GddHigh",
                               "DryDays",
                               "GDD1to3",
                               "Prec1to3",
                               "GddHigh1to3",
                               "DryDays1to3",
                               "GDD4to10",
                               "Prec4to10",
                               "GddHigh4to10",
                               "DryDays4to10"))


# Insert LHS indicators
tab_structural$order <- NULL
tab_structural <- berryFunctions::insertRows(tab_structural, 1 , new = c("Cholesky terms","",""))
tab_structural <- berryFunctions::insertRows(tab_structural, 20, new = c("Cereals supply","",""))
tab_structural <- berryFunctions::insertRows(tab_structural, 40, new = c("Protein crops supply","",""))
tab_structural <- berryFunctions::insertRows(tab_structural, 60, new = c("Oilseeds supply","",""))
tab_structural <- berryFunctions::insertRows(tab_structural, 80, new = c("Root crops supply","",""))
tab_structural <- berryFunctions::insertRows(tab_structural, 100, new = c("Corn supply","",""))
tab_structural <- berryFunctions::insertRows(tab_structural, 120, new = c("Fertilizer demand","",""))

# Print
write_xlsx(tab_structural,"Tables/Table_S6.5.1.xlsx")

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
write_xlsx(tab_PriceElast,"Tables/Table_S6.5.2.xlsx")

#---------------------------------------------------#
#### Table: Weather effects on output and inputs ####
#---------------------------------------------------#
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
write_xlsx(tab_WeatherElast,"Tables/Table_S6.5.3.xlsx")
