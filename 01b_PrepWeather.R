#------------------------------------------#
# This program prepares the weather data:  #
# 1) Lag structures                        #
# 2) Weather simulations                   #
#                                          #
#------------------------------------------#

#-----------------------------#
#### 0) Some preliminaries ####
#-----------------------------#

# open libraries
library(dplyr)
library(readxl) 
library(writexl)
library(stringr)
library(plm) 


# load weather data
df_weather <- read_excel("Data/weatherMarAug.xlsx", col_names = TRUE)



# Add a zero to the nuts5 names if it has only 7 digits
df_weather$nuts5 <- ifelse(str_count(df_weather$nuts5)==7,paste0("0", df_weather$nuts5),df_weather$nuts5)

# make panel data
df_weather <- pdata.frame(df_weather,index=c("nuts5","year"))  

#-------------------------#
#### 1) Lag structures ####
#-------------------------#

#Create dataframes in which lags can be stored
df_weather1to3 <- df_weather
df_weather4to10 <- df_weather
df_weather1to5 <- df_weather
df_weather6to10 <- df_weather
df_weather6to20 <- df_weather
df_weather1 <- df_weather
df_weather2to10 <- df_weather
df_weather1to10 <- df_weather
df_weather1to30 <- df_weather

# Lag1to3 for all variables

    # calculate and save in new data frame
    for(i in 3:6) {       # for-loop over columns that contain weather vars
      df_weather1to3[ , i] <- (   plm::lag(df_weather[ , i],1) +
                                  plm::lag(df_weather[ , i],2) +
                                  plm::lag(df_weather[ , i],3) ) / 3
    }
    
    #rename variables (to indicate which lags they represent)
    df_weather1to3 <- dplyr::rename(df_weather1to3, c("gdd_1to3"="gdd",
                                                      "gddHigh_1to3" = "gddHigh",
                                                      "prec_1to3" = "prec",
                                                      "dd_1to3" = "dd"))
    

# Lag4to10 for all variables
    
    # calculate and save in new data frame
    for(i in 3:6) {       # for-loop over columns
      df_weather4to10[ , i] <- (    plm::lag(df_weather[ , i],4) +
                                    plm::lag(df_weather[ , i],5) +
                                    plm::lag(df_weather[ , i],6) +
                                    plm::lag(df_weather[ , i],7) +
                                    plm::lag(df_weather[ , i],8) +
                                    plm::lag(df_weather[ , i],9) +
                                    plm::lag(df_weather[ , i],10) ) / 7
    }
    
    #rename variables
    df_weather4to10 <- dplyr::rename(df_weather4to10, c("gdd_4to10"="gdd",
                                                        "gddHigh_4to10" = "gddHigh",
                                                        "prec_4to10" = "prec",
                                                        "dd_4to10" = "dd"))
 
# Lag1to5 for all variables

    # calculate and save in new data frame
    for(i in 3:6) {       # for-loop over columns that contain weather vars
      df_weather1to5[ , i] <- (   plm::lag(df_weather[ , i],1) +
                                  plm::lag(df_weather[ , i],2) +
                                  plm::lag(df_weather[ , i],3) +
                                  plm::lag(df_weather[ , i],4) +
                                  plm::lag(df_weather[ , i],5) ) / 5
    }
    
    #rename variables (to indicate which lags they represent)
    df_weather1to5 <- dplyr::rename(df_weather1to5, c("gdd_1to5"="gdd",
                                                      "gddHigh_1to5" = "gddHigh",
                                                      "prec_1to5" = "prec",
                                                      "dd_1to5" = "dd"))

# Lag6to10 for all variables
    
    # calculate and save in new data frame
    for(i in 3:6) {       # for-loop over columns
      df_weather6to10[ , i] <- (    plm::lag(df_weather[ , i],6) +
                                    plm::lag(df_weather[ , i],7) +
                                    plm::lag(df_weather[ , i],8) +
                                    plm::lag(df_weather[ , i],9) +
                                    plm::lag(df_weather[ , i],10) ) / 5
    }
    
    #rename variables
    df_weather6to10 <- dplyr::rename(df_weather6to10, c("gdd_6to10"="gdd",
                                                        "gddHigh_6to10" = "gddHigh",
                                                        "prec_6to10" = "prec",
                                                        "dd_6to10" = "dd"))

# Lag6to20 for all variables
    
    # calculate and save in new data frame
    for(i in 3:6) {       # for-loop over columns
      df_weather6to20[ , i] <- (    plm::lag(df_weather[ , i],6) +
                                    plm::lag(df_weather[ , i],7) +
                                    plm::lag(df_weather[ , i],8) +
                                    plm::lag(df_weather[ , i],9) +
                                    plm::lag(df_weather[ , i],10) +
                                    plm::lag(df_weather[ , i],11) +
                                    plm::lag(df_weather[ , i],12) +
                                    plm::lag(df_weather[ , i],13) +
                                    plm::lag(df_weather[ , i],14) +
                                    plm::lag(df_weather[ , i],15) +
                                    plm::lag(df_weather[ , i],16) +
                                    plm::lag(df_weather[ , i],17) +
                                    plm::lag(df_weather[ , i],18) +
                                    plm::lag(df_weather[ , i],19) +
                                    plm::lag(df_weather[ , i],20) ) / 15
    }
    
    #rename variables
    df_weather6to20 <- dplyr::rename(df_weather6to20, c("gdd_6to20"="gdd",
                                                        "gddHigh_6to20" = "gddHigh",
                                                        "prec_6to20" = "prec",
                                                        "dd_6to20" = "dd"))

# Lag1 for all variables
    
    # calculate and save in new data frame
    for(i in 3:6) {       # for-loop over columns
      df_weather1[ , i] <-    plm::lag(df_weather[ , i],1) 
    }
    
    #rename variables
    df_weather1 <- dplyr::rename(df_weather1, c("gdd_1"="gdd",
                                                    "gddHigh_1" = "gddHigh",
                                                    "prec_1" = "prec",
                                                    "dd_1" = "dd"))

# Lag2to10 for all variables
    
    # calculate and save in new data frame
    for(i in 3:6) {       # for-loop over columns
      df_weather2to10[ , i] <- (    plm::lag(df_weather[ , i],2) +
                                      plm::lag(df_weather[ , i],3) +
                                      plm::lag(df_weather[ , i],4) +
                                      plm::lag(df_weather[ , i],5) +
                                      plm::lag(df_weather[ , i],6) +
                                      plm::lag(df_weather[ , i],7) +
                                      plm::lag(df_weather[ , i],8)  +
                                      plm::lag(df_weather[ , i],9)  +
                                      plm::lag(df_weather[ , i],10)  ) / 9
    }
    
    #rename variables
    df_weather2to10 <- dplyr::rename(df_weather2to10, c("gdd_2to10"="gdd",
                                                        "gddHigh_2to10" = "gddHigh",
                                                        "prec_2to10" = "prec",
                                                        "dd_2to10" = "dd"))
   
# Lag1to10 for all variables
    
    # calculate and save in new data frame
    for(i in 3:6) {       # for-loop over columns
      df_weather1to10[ , i] <- (    plm::lag(df_weather[ , i],1) +
                                      plm::lag(df_weather[ , i],2) +
                                      plm::lag(df_weather[ , i],3) +
                                      plm::lag(df_weather[ , i],4) +
                                      plm::lag(df_weather[ , i],5) +
                                      plm::lag(df_weather[ , i],6) +
                                      plm::lag(df_weather[ , i],7) +
                                      plm::lag(df_weather[ , i],8)  +
                                      plm::lag(df_weather[ , i],9)  +
                                      plm::lag(df_weather[ , i],10)  ) / 10
    }
    
    #rename variables
    df_weather1to10 <- dplyr::rename(df_weather1to10, c("gdd_1to10"="gdd",
                                                        "gddHigh_1to10" = "gddHigh",
                                                        "prec_1to10" = "prec",
                                                        "dd_1to10" = "dd"))

# Lag1to30 for all variables
    
    # calculate and save in new data frame
    for(i in 3:6) {       # for-loop over columns
      df_weather1to30[ , i] <- (    plm::lag(df_weather[ , i],1) +
                                      plm::lag(df_weather[ , i],2) +
                                      plm::lag(df_weather[ , i],3) +
                                      plm::lag(df_weather[ , i],4) +
                                      plm::lag(df_weather[ , i],5) +
                                      plm::lag(df_weather[ , i],6) +
                                      plm::lag(df_weather[ , i],7) +
                                      plm::lag(df_weather[ , i],8) +
                                      plm::lag(df_weather[ , i],9) +
                                      plm::lag(df_weather[ , i],10) +
                                      plm::lag(df_weather[ , i],11) +
                                      plm::lag(df_weather[ , i],12) +
                                      plm::lag(df_weather[ , i],13) +
                                      plm::lag(df_weather[ , i],14) +
                                      plm::lag(df_weather[ , i],15) +
                                      plm::lag(df_weather[ , i],16) +
                                      plm::lag(df_weather[ , i],17) +
                                      plm::lag(df_weather[ , i],18) +
                                      plm::lag(df_weather[ , i],19) +
                                      plm::lag(df_weather[ , i],20) +
                                      plm::lag(df_weather[ , i],21) +
                                      plm::lag(df_weather[ , i],22) + 
                                      plm::lag(df_weather[ , i],23) +
                                      plm::lag(df_weather[ , i],24) +
                                      plm::lag(df_weather[ , i],25) +
                                      plm::lag(df_weather[ , i],26) +
                                      plm::lag(df_weather[ , i],27) +
                                      plm::lag(df_weather[ , i],28) +
                                      plm::lag(df_weather[ , i],29) +
                                      plm::lag(df_weather[ , i],30)  ) / 30
    }
    
    #rename variables
    df_weather1to30 <- dplyr::rename(df_weather1to30, c("gdd_1to30"="gdd",
                                                        "gddHigh_1to30" = "gddHigh",
                                                        "prec_1to30" = "prec",
                                                        "dd_1to30" = "dd"))
    
#merge all weather dataframes
df_weather <- right_join(df_weather, df_weather1to3)    
df_weather <- right_join(df_weather, df_weather4to10)  
df_weather <- right_join(df_weather, df_weather1to5)
df_weather <- right_join(df_weather, df_weather6to10) 
df_weather <- right_join(df_weather, df_weather6to20)
df_weather <- right_join(df_weather, df_weather1)
df_weather <- right_join(df_weather, df_weather2to10)
df_weather <- right_join(df_weather, df_weather1to10)    
df_weather <- right_join(df_weather, df_weather1to30) 

#convert back to non-panel dataframe
df_weather <- as.data.frame(df_weather)

#make all variables except nuts5 numeric
cols.num <- names(df_weather)[-1]
df_weather[cols.num] <- sapply(df_weather[cols.num],as.character)
df_weather[cols.num] <- sapply(df_weather[cols.num],as.numeric)
df_weather$nuts5 <- as.character(df_weather$nuts5)

#------------------------------#
#### 2) Weather simulations ####
#------------------------------#

#define weather vars
vars <- c("gdd", "gddHigh","prec","dd")

#Weather averages by year
weather_by_year <- df_weather %>% 
  group_by(year) %>% 
  summarise_at(vars,
               list(mean=mean))

# Average values during sample period 
df_weather_lta <-  df_weather[df_weather$year>=1996 & df_weather$year < 2020,]
df_weather_lta <- lapply(df_weather_lta[,vars], 
                         function(var)
                           mean(var))

#Generate list with 2018-values
df_weather_2018 <-  df_weather[df_weather$year==2018,]
df_weather_2018 <- lapply(df_weather_2018[,vars], 
                          function(var)
                            mean(var))

#--------------------------------------#
##### 2a) Simulate: 1to3 and 4to10 #####
#--------------------------------------#

sim2018 <- data.frame(matrix(ncol = 9, nrow = 12))
colnames(sim2018) <- c("t","gdd_1to3","gdd_4to10","gddHigh_1to3","gddHigh_4to10","prec_1to3","prec_4to10","dd_1to3","dd_4to10")
sim2018$t <- seq(-1, 10)


sim2018$gdd_1to3 <- c( df_weather_lta$gdd, 
                      df_weather_lta$gdd,
                      (1/3) * df_weather_2018$gdd + (2/3) * df_weather_lta$gdd, 
                      (1/3) * df_weather_2018$gdd + (2/3) * df_weather_lta$gdd, 
                      (1/3) * df_weather_2018$gdd + (2/3) * df_weather_lta$gdd, 
                      df_weather_lta$gdd, 
                      df_weather_lta$gdd, 
                      df_weather_lta$gdd, 
                      df_weather_lta$gdd, 
                      df_weather_lta$gdd, 
                      df_weather_lta$gdd, 
                      df_weather_lta$gdd )

sim2018$gdd_4to10 <- c(df_weather_lta$gdd, 
                      df_weather_lta$gdd, 
                      df_weather_lta$gdd, 
                      df_weather_lta$gdd, 
                      df_weather_lta$gdd, 
                      (1/7)*df_weather_2018$gdd + (6/7)*df_weather_lta$gdd, 
                      (1/7)*df_weather_2018$gdd + (6/7)*df_weather_lta$gdd, 
                      (1/7)*df_weather_2018$gdd + (6/7)*df_weather_lta$gdd, 
                      (1/7)*df_weather_2018$gdd + (6/7)*df_weather_lta$gdd, 
                      (1/7)*df_weather_2018$gdd + (6/7)*df_weather_lta$gdd, 
                      (1/7)*df_weather_2018$gdd + (6/7)*df_weather_lta$gdd, 
                      (1/7)*df_weather_2018$gdd + (6/7)*df_weather_lta$gdd )

sim2018$gddHigh_1to3 <- c(df_weather_lta$gddHigh, 
                         df_weather_lta$gddHigh, 
                         (1/3) * df_weather_2018$gddHigh + (2/3) * df_weather_lta$gddHigh, 
                         (1/3) * df_weather_2018$gddHigh + (2/3) * df_weather_lta$gddHigh, 
                         (1/3) * df_weather_2018$gddHigh + (2/3) * df_weather_lta$gddHigh, 
                         df_weather_lta$gddHigh, 
                         df_weather_lta$gddHigh, 
                         df_weather_lta$gddHigh, 
                         df_weather_lta$gddHigh, 
                         df_weather_lta$gddHigh, 
                         df_weather_lta$gddHigh, 
                         df_weather_lta$gddHigh )

sim2018$gddHigh_4to10 <- c(df_weather_lta$gddHigh, 
                          df_weather_lta$gddHigh, 
                          df_weather_lta$gddHigh, 
                          df_weather_lta$gddHigh, 
                          df_weather_lta$gddHigh, 
                          (1/7)*df_weather_2018$gddHigh + (6/7)*df_weather_lta$gddHigh, 
                          (1/7)*df_weather_2018$gddHigh + (6/7)*df_weather_lta$gddHigh, 
                          (1/7)*df_weather_2018$gddHigh + (6/7)*df_weather_lta$gddHigh, 
                          (1/7)*df_weather_2018$gddHigh + (6/7)*df_weather_lta$gddHigh, 
                          (1/7)*df_weather_2018$gddHigh + (6/7)*df_weather_lta$gddHigh, 
                          (1/7)*df_weather_2018$gddHigh + (6/7)*df_weather_lta$gddHigh, 
                          (1/7)*df_weather_2018$gddHigh + (6/7)*df_weather_lta$gddHigh )

sim2018$prec_1to3 <- c( df_weather_lta$prec, 
                       df_weather_lta$prec, 
                       (1/3) * df_weather_2018$prec + (2/3) * df_weather_lta$prec, 
                       (1/3) * df_weather_2018$prec + (2/3) * df_weather_lta$prec, 
                       (1/3) * df_weather_2018$prec + (2/3) * df_weather_lta$prec, 
                       df_weather_lta$prec, 
                       df_weather_lta$prec, 
                       df_weather_lta$prec, 
                       df_weather_lta$prec, 
                       df_weather_lta$prec, 
                       df_weather_lta$prec, 
                       df_weather_lta$prec )

sim2018$prec_4to10 <- c( df_weather_lta$prec, 
                        df_weather_lta$prec, 
                        df_weather_lta$prec, 
                        df_weather_lta$prec, 
                        df_weather_lta$prec, 
                        (1/7)*df_weather_2018$prec + (6/7)*df_weather_lta$prec, 
                        (1/7)*df_weather_2018$prec + (6/7)*df_weather_lta$prec, 
                        (1/7)*df_weather_2018$prec + (6/7)*df_weather_lta$prec, 
                        (1/7)*df_weather_2018$prec + (6/7)*df_weather_lta$prec, 
                        (1/7)*df_weather_2018$prec + (6/7)*df_weather_lta$prec, 
                        (1/7)*df_weather_2018$prec + (6/7)*df_weather_lta$prec, 
                        (1/7)*df_weather_2018$prec + (6/7)*df_weather_lta$prec )

sim2018$dd_1to3 <- c( df_weather_lta$dd, 
                     df_weather_lta$dd, 
                     (1/3) * df_weather_2018$dd + (2/3) * df_weather_lta$dd, 
                     (1/3) * df_weather_2018$dd + (2/3) * df_weather_lta$dd, 
                     (1/3) * df_weather_2018$dd + (2/3) * df_weather_lta$dd, 
                     df_weather_lta$dd, 
                     df_weather_lta$dd, 
                     df_weather_lta$dd, 
                     df_weather_lta$dd, 
                     df_weather_lta$dd, 
                     df_weather_lta$dd, 
                     df_weather_lta$dd )

sim2018$dd_4to10 <- c( df_weather_lta$dd, 
                      df_weather_lta$dd, 
                      df_weather_lta$dd, 
                      df_weather_lta$dd, 
                      df_weather_lta$dd, 
                      (1/7)*df_weather_2018$dd + (6/7)*df_weather_lta$dd, 
                      (1/7)*df_weather_2018$dd + (6/7)*df_weather_lta$dd, 
                      (1/7)*df_weather_2018$dd + (6/7)*df_weather_lta$dd, 
                      (1/7)*df_weather_2018$dd + (6/7)*df_weather_lta$dd, 
                      (1/7)*df_weather_2018$dd + (6/7)*df_weather_lta$dd, 
                      (1/7)*df_weather_2018$dd + (6/7)*df_weather_lta$dd, 
                      (1/7)*df_weather_2018$dd + (6/7)*df_weather_lta$dd )

# Make table to be presented in paper

  tab_sim2018 <- sim2018
  
  # Add observed weather to simulated weather
  tab_sim2018$gdd_obs <- c(df_weather_lta$gdd, df_weather_2018$gdd, rep(df_weather_lta$gdd,10))
  tab_sim2018$prec_obs <- c(df_weather_lta$prec, df_weather_2018$prec, rep(df_weather_lta$prec,10))
  tab_sim2018$gddHigh_obs <- c(df_weather_lta$gddHigh, df_weather_2018$gddHigh, rep(df_weather_lta$gddHigh,10))
  tab_sim2018$dd_obs <- c(df_weather_lta$dd, df_weather_2018$dd, rep(df_weather_lta$dd,10))
  
  # Create table to be printed
  tab_sim2018 <- data.frame(gdd_obs = tab_sim2018$gdd_obs,
                            prec_obs = tab_sim2018$prec_obs,
                            gddHigh_obs = tab_sim2018$gddHigh_obs,
                            dd_obs = tab_sim2018$dd_obs,
                            gdd_1to3 = sim2018$gdd_1to3,
                            prec_1to3 = sim2018$prec_1to3,
                            gddHigh_1to3 = sim2018$gddHigh_1to3,
                            dd_1to3 = sim2018$dd_1to3,
                            gdd_4to10 = sim2018$gdd_4to10,
                            prec_4to10 = sim2018$prec_4to10,
                            gddHigh_4to10 = sim2018$gddHigh_4to10,
                            dd_4to10 = sim2018$dd_4to10)

  # Print
  write_xlsx(tab_sim2018,path="Tables/Table_S2.1.xlsx")
  
#--------------------------------------#
##### 2b) Simulate: 1to5 and 6to10 #####
#--------------------------------------#

sim2018_1to5_6to10 <- data.frame(matrix(ncol = 9, nrow = 12))
colnames(sim2018_1to5_6to10) <- c("t","gdd_1to5","gdd_6to10","gddHigh_1to5","gddHigh_6to10","prec_1to5","prec_6to10","dd_1to5","dd_6to10")
sim2018_1to5_6to10$t <- seq(-1, 10)

sim2018_1to5_6to10$gdd_1to5 <- c(df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                (1/5) * df_weather_2018$gdd + (4/5) * df_weather_lta$gdd, 
                                (1/5) * df_weather_2018$gdd + (4/5) * df_weather_lta$gdd, 
                                (1/5) * df_weather_2018$gdd + (4/5) * df_weather_lta$gdd,
                                (1/5) * df_weather_2018$gdd + (4/5) * df_weather_lta$gdd,
                                (1/5) * df_weather_2018$gdd + (4/5) * df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd)

sim2018_1to5_6to10$gdd_6to10 <- c(df_weather_lta$gdd, 
                                 df_weather_lta$gdd, 
                                 df_weather_lta$gdd, 
                                 df_weather_lta$gdd, 
                                 df_weather_lta$gdd, 
                                 df_weather_lta$gdd, 
                                 df_weather_lta$gdd, 
                                 (1/5)*df_weather_2018$gdd + (4/5)*df_weather_lta$gdd, 
                                 (1/5)*df_weather_2018$gdd + (4/5)*df_weather_lta$gdd, 
                                 (1/5)*df_weather_2018$gdd + (4/5)*df_weather_lta$gdd, 
                                 (1/5)*df_weather_2018$gdd + (4/5)*df_weather_lta$gdd, 
                                 (1/5)*df_weather_2018$gdd + (4/5)*df_weather_lta$gdd ) 

sim2018_1to5_6to10$gddHigh_1to5 <- c(df_weather_lta$gddHigh,
                                    df_weather_lta$gddHigh,
                                    (1/5) * df_weather_2018$gddHigh + (4/5) * df_weather_lta$gddHigh, 
                                    (1/5) * df_weather_2018$gddHigh + (4/5) * df_weather_lta$gddHigh, 
                                    (1/5) * df_weather_2018$gddHigh + (4/5) * df_weather_lta$gddHigh,
                                    (1/5) * df_weather_2018$gddHigh + (4/5) * df_weather_lta$gddHigh,
                                    (1/5) * df_weather_2018$gddHigh + (4/5) * df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh)

sim2018_1to5_6to10$gddHigh_6to10 <- c( df_weather_lta$gddHigh,
                                     df_weather_lta$gddHigh,
                                     df_weather_lta$gddHigh, 
                                     df_weather_lta$gddHigh, 
                                     df_weather_lta$gddHigh, 
                                     df_weather_lta$gddHigh, 
                                     df_weather_lta$gddHigh, 
                                     (1/5)*df_weather_2018$gddHigh + (4/5)*df_weather_lta$gddHigh, 
                                     (1/5)*df_weather_2018$gddHigh + (4/5)*df_weather_lta$gddHigh, 
                                     (1/5)*df_weather_2018$gddHigh + (4/5)*df_weather_lta$gddHigh, 
                                     (1/5)*df_weather_2018$gddHigh + (4/5)*df_weather_lta$gddHigh, 
                                     (1/5)*df_weather_2018$gddHigh + (4/5)*df_weather_lta$gddHigh) 

sim2018_1to5_6to10$prec_1to5 <- c( df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  (1/5) * df_weather_2018$prec + (4/5) * df_weather_lta$prec, 
                                  (1/5) * df_weather_2018$prec + (4/5) * df_weather_lta$prec, 
                                  (1/5) * df_weather_2018$prec + (4/5) * df_weather_lta$prec,
                                  (1/5) * df_weather_2018$prec + (4/5) * df_weather_lta$prec,
                                  (1/5) * df_weather_2018$prec + (4/5) * df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec)

sim2018_1to5_6to10$prec_6to10 <- c( df_weather_lta$prec, 
                                   df_weather_lta$prec, 
                                   df_weather_lta$prec, 
                                   df_weather_lta$prec, 
                                   df_weather_lta$prec, 
                                   df_weather_lta$prec, 
                                   df_weather_lta$prec, 
                                   (1/5)*df_weather_2018$prec + (4/5)*df_weather_lta$prec, 
                                   (1/5)*df_weather_2018$prec + (4/5)*df_weather_lta$prec, 
                                   (1/5)*df_weather_2018$prec + (4/5)*df_weather_lta$prec, 
                                   (1/5)*df_weather_2018$prec + (4/5)*df_weather_lta$prec, 
                                   (1/5)*df_weather_2018$prec + (4/5)*df_weather_lta$prec) 

sim2018_1to5_6to10$dd_1to5 <- c( df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                (1/5) * df_weather_2018$dd + (4/5) * df_weather_lta$dd, 
                                (1/5) * df_weather_2018$dd + (4/5) * df_weather_lta$dd, 
                                (1/5) * df_weather_2018$dd + (4/5) * df_weather_lta$dd,
                                (1/5) * df_weather_2018$dd + (4/5) * df_weather_lta$dd,
                                (1/5) * df_weather_2018$dd + (4/5) * df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd)

sim2018_1to5_6to10$dd_6to10 <- c( df_weather_lta$dd, 
                                 df_weather_lta$dd, 
                                 df_weather_lta$dd, 
                                 df_weather_lta$dd, 
                                 df_weather_lta$dd, 
                                 df_weather_lta$dd, 
                                 df_weather_lta$dd, 
                                 (1/5)*df_weather_2018$dd + (4/5)*df_weather_lta$dd, 
                                 (1/5)*df_weather_2018$dd + (4/5)*df_weather_lta$dd, 
                                 (1/5)*df_weather_2018$dd + (4/5)*df_weather_lta$dd, 
                                 (1/5)*df_weather_2018$dd + (4/5)*df_weather_lta$dd, 
                                 (1/5)*df_weather_2018$dd + (4/5)*df_weather_lta$dd) 


#--------------------------------------#
##### 2c) Simulate: 1to5 and 6to20 #####
#--------------------------------------#

sim2018_1to5_6to20 <- data.frame(matrix(ncol = 9, nrow = 22))
colnames(sim2018_1to5_6to20) <- c("t","gdd_1to5","gdd_6to20","gddHigh_1to5","gddHigh_6to20","prec_1to5","prec_6to20","dd_1to5","dd_6to20")
sim2018_1to5_6to20$t <- seq(-1, 20)

sim2018_1to5_6to20$gdd_1to5 <- c(df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                (1/5) * df_weather_2018$gdd + (4/5) * df_weather_lta$gdd, 
                                (1/5) * df_weather_2018$gdd + (4/5) * df_weather_lta$gdd, 
                                (1/5) * df_weather_2018$gdd + (4/5) * df_weather_lta$gdd,
                              	(1/5) * df_weather_2018$gdd + (4/5) * df_weather_lta$gdd,
                              	(1/5) * df_weather_2018$gdd + (4/5) * df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd, 
                                df_weather_lta$gdd)

sim2018_1to5_6to20$gdd_6to20 <- c(df_weather_lta$gdd, 
                                 df_weather_lta$gdd, 
                                 df_weather_lta$gdd, 
                                 df_weather_lta$gdd, 
                                 df_weather_lta$gdd, 
                                 df_weather_lta$gdd, 
                                 df_weather_lta$gdd, 
                                 (1/15)*df_weather_2018$gdd + (14/15)*df_weather_lta$gdd, 
                                 (1/15)*df_weather_2018$gdd + (14/15)*df_weather_lta$gdd, 
                                 (1/15)*df_weather_2018$gdd + (14/15)*df_weather_lta$gdd, 
                                 (1/15)*df_weather_2018$gdd + (14/15)*df_weather_lta$gdd, 
                                 (1/15)*df_weather_2018$gdd + (14/15)*df_weather_lta$gdd, 
                                 (1/15)*df_weather_2018$gdd + (14/15)*df_weather_lta$gdd, 
                                 (1/15)*df_weather_2018$gdd + (14/15)*df_weather_lta$gdd, 
                                 (1/15)*df_weather_2018$gdd + (14/15)*df_weather_lta$gdd, 
                                 (1/15)*df_weather_2018$gdd + (14/15)*df_weather_lta$gdd, 
                                 (1/15)*df_weather_2018$gdd + (14/15)*df_weather_lta$gdd, 
                                 (1/15)*df_weather_2018$gdd + (14/15)*df_weather_lta$gdd, 
                                 (1/15)*df_weather_2018$gdd + (14/15)*df_weather_lta$gdd, 
                                 (1/15)*df_weather_2018$gdd + (14/15)*df_weather_lta$gdd, 
                                 (1/15)*df_weather_2018$gdd + (14/15)*df_weather_lta$gdd, 
                                 (1/15)*df_weather_2018$gdd + (14/15)*df_weather_lta$gdd ) 

sim2018_1to5_6to20$gddHigh_1to5 <- c( df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    (1/5) * df_weather_2018$gddHigh + (4/5) * df_weather_lta$gddHigh, 
                                    (1/5) * df_weather_2018$gddHigh + (4/5) * df_weather_lta$gddHigh, 
                                    (1/5) * df_weather_2018$gddHigh + (4/5) * df_weather_lta$gddHigh,
                            			  (1/5) * df_weather_2018$gddHigh + (4/5) * df_weather_lta$gddHigh,
                            			  (1/5) * df_weather_2018$gddHigh + (4/5) * df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh, 
                                    df_weather_lta$gddHigh)

sim2018_1to5_6to20$gddHigh_6to20 <- c(df_weather_lta$gddHigh, 
                                     df_weather_lta$gddHigh, 
                                     df_weather_lta$gddHigh, 
                                     df_weather_lta$gddHigh, 
                                     df_weather_lta$gddHigh, 
                                     df_weather_lta$gddHigh, 
                                     df_weather_lta$gddHigh, 
                                     (1/15)*df_weather_2018$gddHigh + (14/15)*df_weather_lta$gddHigh, 
                                     (1/15)*df_weather_2018$gddHigh + (14/15)*df_weather_lta$gddHigh, 
                                     (1/15)*df_weather_2018$gddHigh + (14/15)*df_weather_lta$gddHigh, 
                                     (1/15)*df_weather_2018$gddHigh + (14/15)*df_weather_lta$gddHigh, 
                                     (1/15)*df_weather_2018$gddHigh + (14/15)*df_weather_lta$gddHigh, 
                                     (1/15)*df_weather_2018$gddHigh + (14/15)*df_weather_lta$gddHigh, 
                                     (1/15)*df_weather_2018$gddHigh + (14/15)*df_weather_lta$gddHigh, 
                                     (1/15)*df_weather_2018$gddHigh + (14/15)*df_weather_lta$gddHigh, 
                                     (1/15)*df_weather_2018$gddHigh + (14/15)*df_weather_lta$gddHigh, 
                                     (1/15)*df_weather_2018$gddHigh + (14/15)*df_weather_lta$gddHigh, 
                                     (1/15)*df_weather_2018$gddHigh + (14/15)*df_weather_lta$gddHigh, 
                                     (1/15)*df_weather_2018$gddHigh + (14/15)*df_weather_lta$gddHigh, 
                                     (1/15)*df_weather_2018$gddHigh + (14/15)*df_weather_lta$gddHigh, 
                                     (1/15)*df_weather_2018$gddHigh + (14/15)*df_weather_lta$gddHigh, 
                                     (1/15)*df_weather_2018$gddHigh + (14/15)*df_weather_lta$gddHigh ) 

sim2018_1to5_6to20$prec_1to5 <- c( df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  (1/5) * df_weather_2018$prec + (4/5) * df_weather_lta$prec, 
                                  (1/5) * df_weather_2018$prec + (4/5) * df_weather_lta$prec, 
                                  (1/5) * df_weather_2018$prec + (4/5) * df_weather_lta$prec,
                          			  (1/5) * df_weather_2018$prec + (4/5) * df_weather_lta$prec,
                          			  (1/5) * df_weather_2018$prec + (4/5) * df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec, 
                                  df_weather_lta$prec)

sim2018_1to5_6to20$prec_6to20 <- c( df_weather_lta$prec, 
                                   df_weather_lta$prec, 
                                   df_weather_lta$prec, 
                                   df_weather_lta$prec, 
                                   df_weather_lta$prec, 
                                   df_weather_lta$prec, 
                                   df_weather_lta$prec, 
                                   (1/15)*df_weather_2018$prec + (14/15)*df_weather_lta$prec, 
                                   (1/15)*df_weather_2018$prec + (14/15)*df_weather_lta$prec, 
                                   (1/15)*df_weather_2018$prec + (14/15)*df_weather_lta$prec, 
                                   (1/15)*df_weather_2018$prec + (14/15)*df_weather_lta$prec, 
                                   (1/15)*df_weather_2018$prec + (14/15)*df_weather_lta$prec, 
                                   (1/15)*df_weather_2018$prec + (14/15)*df_weather_lta$prec, 
                                   (1/15)*df_weather_2018$prec + (14/15)*df_weather_lta$prec, 
                                   (1/15)*df_weather_2018$prec + (14/15)*df_weather_lta$prec, 
                                   (1/15)*df_weather_2018$prec + (14/15)*df_weather_lta$prec, 
                                   (1/15)*df_weather_2018$prec + (14/15)*df_weather_lta$prec, 
                                   (1/15)*df_weather_2018$prec + (14/15)*df_weather_lta$prec, 
                                   (1/15)*df_weather_2018$prec + (14/15)*df_weather_lta$prec, 
                                   (1/15)*df_weather_2018$prec + (14/15)*df_weather_lta$prec, 
                                   (1/15)*df_weather_2018$prec + (14/15)*df_weather_lta$prec, 
                                   (1/15)*df_weather_2018$prec + (14/15)*df_weather_lta$prec ) 

sim2018_1to5_6to20$dd_1to5 <- c( df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                (1/5) * df_weather_2018$dd + (4/5) * df_weather_lta$dd, 
                                (1/5) * df_weather_2018$dd + (4/5) * df_weather_lta$dd, 
                                (1/5) * df_weather_2018$dd + (4/5) * df_weather_lta$dd,
                        			  (1/5) * df_weather_2018$dd + (4/5) * df_weather_lta$dd,
                        			  (1/5) * df_weather_2018$dd + (4/5) * df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd, 
                                df_weather_lta$dd)

sim2018_1to5_6to20$dd_6to20 <- c( df_weather_lta$dd, 
                                 df_weather_lta$dd, 
                                 df_weather_lta$dd, 
                                 df_weather_lta$dd, 
                                 df_weather_lta$dd, 
                                 df_weather_lta$dd, 
                                 df_weather_lta$dd, 
                                 (1/15)*df_weather_2018$dd + (14/15)*df_weather_lta$dd, 
                                 (1/15)*df_weather_2018$dd + (14/15)*df_weather_lta$dd, 
                                 (1/15)*df_weather_2018$dd + (14/15)*df_weather_lta$dd, 
                                 (1/15)*df_weather_2018$dd + (14/15)*df_weather_lta$dd, 
                                 (1/15)*df_weather_2018$dd + (14/15)*df_weather_lta$dd, 
                                 (1/15)*df_weather_2018$dd + (14/15)*df_weather_lta$dd, 
                                 (1/15)*df_weather_2018$dd + (14/15)*df_weather_lta$dd, 
                                 (1/15)*df_weather_2018$dd + (14/15)*df_weather_lta$dd, 
                                 (1/15)*df_weather_2018$dd + (14/15)*df_weather_lta$dd, 
                                 (1/15)*df_weather_2018$dd + (14/15)*df_weather_lta$dd, 
                                 (1/15)*df_weather_2018$dd + (14/15)*df_weather_lta$dd, 
                                 (1/15)*df_weather_2018$dd + (14/15)*df_weather_lta$dd, 
                                 (1/15)*df_weather_2018$dd + (14/15)*df_weather_lta$dd, 
                                 (1/15)*df_weather_2018$dd + (14/15)*df_weather_lta$dd, 
                                 (1/15)*df_weather_2018$dd + (14/15)*df_weather_lta$dd ) 

#-----------------------------------#
##### 2e) Simulate: 1 and 2to10 #####
#-----------------------------------#

sim2018_1_2to10 <- data.frame(matrix(ncol = 9, nrow = 12))
colnames(sim2018_1_2to10) <- c("t","gdd_1","gdd_2to10","gddHigh_1","gddHigh_2to10","prec_1","prec_2to10","dd_1","dd_2to10")
sim2018_1_2to10$t <- seq(-1, 10)

sim2018_1_2to10$gdd_1 <- c(df_weather_lta$gdd, 
                          df_weather_lta$gdd, 
                          df_weather_2018$gdd, 
                          df_weather_lta$gdd, 
                          df_weather_lta$gdd, 
                          df_weather_lta$gdd, 
                          df_weather_lta$gdd, 
                          df_weather_lta$gdd, 
                          df_weather_lta$gdd, 
                          df_weather_lta$gdd, 
                          df_weather_lta$gdd, 
                          df_weather_lta$gdd)

sim2018_1_2to10$gdd_2to10 <- c(df_weather_lta$gdd, 
                              df_weather_lta$gdd, 
                              df_weather_lta$gdd, 
                              (1/9)*df_weather_2018$gdd + (8/9)*df_weather_lta$gdd, 
                              (1/9)*df_weather_2018$gdd + (8/9)*df_weather_lta$gdd, 
                              (1/9)*df_weather_2018$gdd + (8/9)*df_weather_lta$gdd, 
                              (1/9)*df_weather_2018$gdd + (8/9)*df_weather_lta$gdd, 
                              (1/9)*df_weather_2018$gdd + (8/9)*df_weather_lta$gdd, 
                              (1/9)*df_weather_2018$gdd + (8/9)*df_weather_lta$gdd, 
                              (1/9)*df_weather_2018$gdd + (8/9)*df_weather_lta$gdd, 
                              (1/9)*df_weather_2018$gdd + (8/9)*df_weather_lta$gdd, 
                              (1/9)*df_weather_2018$gdd + (8/9)*df_weather_lta$gdd ) 

sim2018_1_2to10$gddHigh_1 <- c(df_weather_lta$gddHigh, 
                              df_weather_lta$gddHigh, 
                              df_weather_2018$gddHigh, 
                              df_weather_lta$gddHigh,  
                              df_weather_lta$gddHigh, 
                              df_weather_lta$gddHigh, 
                              df_weather_lta$gddHigh, 
                              df_weather_lta$gddHigh, 
                              df_weather_lta$gddHigh, 
                              df_weather_lta$gddHigh, 
                              df_weather_lta$gddHigh, 
                              df_weather_lta$gddHigh)

sim2018_1_2to10$gddHigh_2to10 <- c(df_weather_lta$gddHigh, 
                                  df_weather_lta$gddHigh, 
                                  df_weather_lta$gddHigh, 
                                  (1/9)*df_weather_2018$gddHigh + (8/9)*df_weather_lta$gddHigh, 
                                  (1/9)*df_weather_2018$gddHigh + (8/9)*df_weather_lta$gddHigh, 
                                  (1/9)*df_weather_2018$gddHigh + (8/9)*df_weather_lta$gddHigh, 
                                  (1/9)*df_weather_2018$gddHigh + (8/9)*df_weather_lta$gddHigh, 
                                  (1/9)*df_weather_2018$gddHigh + (8/9)*df_weather_lta$gddHigh, 
                                  (1/9)*df_weather_2018$gddHigh + (8/9)*df_weather_lta$gddHigh, 
                                  (1/9)*df_weather_2018$gddHigh + (8/9)*df_weather_lta$gddHigh, 
                                  (1/9)*df_weather_2018$gddHigh + (8/9)*df_weather_lta$gddHigh, 
                                  (1/9)*df_weather_2018$gddHigh + (8/9)*df_weather_lta$gddHigh) 

sim2018_1_2to10$prec_1 <- c( df_weather_lta$prec, 
                            df_weather_lta$prec, 
                            df_weather_2018$prec, 
                            df_weather_lta$prec, 
                            df_weather_lta$prec, 
                            df_weather_lta$prec, 
                            df_weather_lta$prec, 
                            df_weather_lta$prec, 
                            df_weather_lta$prec, 
                            df_weather_lta$prec, 
                            df_weather_lta$prec, 
                            df_weather_lta$prec)

sim2018_1_2to10$prec_2to10 <- c( df_weather_lta$prec, 
                                df_weather_lta$prec, 
                                df_weather_lta$prec, 
                                (1/9)*df_weather_2018$prec + (8/9)*df_weather_lta$prec,
                                (1/9)*df_weather_2018$prec + (8/9)*df_weather_lta$prec,
                                (1/9)*df_weather_2018$prec + (8/9)*df_weather_lta$prec,
                                (1/9)*df_weather_2018$prec + (8/9)*df_weather_lta$prec,
                                (1/9)*df_weather_2018$prec + (8/9)*df_weather_lta$prec, 
                                (1/9)*df_weather_2018$prec + (8/9)*df_weather_lta$prec, 
                                (1/9)*df_weather_2018$prec + (8/9)*df_weather_lta$prec, 
                                (1/9)*df_weather_2018$prec + (8/9)*df_weather_lta$prec, 
                                (1/9)*df_weather_2018$prec + (8/9)*df_weather_lta$prec) 

sim2018_1_2to10$dd_1 <- c( df_weather_lta$dd, 
                          df_weather_lta$dd, 
                          df_weather_2018$dd, 
                          df_weather_lta$dd, 
                          df_weather_lta$dd, 
                          df_weather_lta$dd, 
                          df_weather_lta$dd, 
                          df_weather_lta$dd, 
                          df_weather_lta$dd, 
                          df_weather_lta$dd, 
                          df_weather_lta$dd, 
                          df_weather_lta$dd)

sim2018_1_2to10$dd_2to10 <- c( df_weather_lta$dd, 
                              df_weather_lta$dd, 
                              df_weather_lta$dd, 
                              (1/9)*df_weather_2018$dd + (8/9)*df_weather_lta$dd, 
                              (1/9)*df_weather_2018$dd + (8/9)*df_weather_lta$dd, 
                              (1/9)*df_weather_2018$dd + (8/9)*df_weather_lta$dd, 
                              (1/9)*df_weather_2018$dd + (8/9)*df_weather_lta$dd, 
                              (1/9)*df_weather_2018$dd + (8/9)*df_weather_lta$dd, 
                              (1/9)*df_weather_2018$dd + (8/9)*df_weather_lta$dd, 
                              (1/9)*df_weather_2018$dd + (8/9)*df_weather_lta$dd, 
                              (1/9)*df_weather_2018$dd + (8/9)*df_weather_lta$dd, 
                              (1/9)*df_weather_2018$dd + (8/9)*df_weather_lta$dd) 

#----------------------------------#
#                                  #
# End of weather data preparation  # 
#                                  #
#----------------------------------#

# save data frame with weather vars (realized and lags, to be merged with farm data)
save(df_weather,file="rOutput/df_weather.Rda")

# save data frame with LTA weather
save(df_weather_lta,file="rOutput/df_weather_lta.Rda")

# save data frame with 2018 weather 
save(df_weather_2018,file="rOutput/df_weather_2018.Rda")

# save data frame with simulated 2018 drought shock
save(sim2018,file="rOutput/sim2018.Rda")

# save simulated 2018 with alternative lag structures (for robustness checks)
save(sim2018_1to5_6to10,file="rOutput/sim2018_1to5_6to10.Rda")
save(sim2018_1to5_6to20,file="rOutput/sim2018_1to5_6to20.Rda")
save(sim2018_1_2to10,file="rOutput/sim2018_1_2to10.Rda")
