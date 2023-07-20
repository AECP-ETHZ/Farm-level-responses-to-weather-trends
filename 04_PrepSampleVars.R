#-----------------------------------------------------------#
#                                                           #
# This program prepares the sample and all variables used   #
# in the estimation of the profit system.                   #
#                                                           #
#   0)  Some preliminaries                                  # 
#   1)  Define weather variables                            #
#   2)  Define price variables                              #
#   3)  Define input quantities                             #
#   4)  Prepare land shares and lags thereof                #
#   5)  Create dummies for crop production                  #
#   6)  Construct sample                                    #
#   7)  Keep only complete observations for the estimation  #
#   8)  Calculate descriptive statistics                    #
#   9)  Calculate interaction terms to assess heterogeneity #
#   10) Calculate farm means for Mundlak device             #
#                                                           #
#-----------------------------------------------------------#

#-----------------------------#
#### 0) Some preliminaries ####
#-----------------------------#

# Open packages
library(dplyr)
library(stargazer)
library(GJRM) 
library(sampleSelection) 
library(systemfit) 

# Load data
load("rOutput/farm_outp_inp.Rda")

#-----------------------------------#
#### 1) Define weather variables ####
#-----------------------------------#

#observed weather
df_farm$gdd_obs <- df_farm$gdd
df_farm$gdd2_obs <- df_farm$gdd*df_farm$gdd
df_farm$prec_obs <- df_farm$prec
df_farm$prec2_obs <- df_farm$prec*df_farm$prec
df_farm$gddHigh_obs <- df_farm$gddHigh
df_farm$dd_obs <- df_farm$dd

#non-linear terms for past weather

  #1to3
  df_farm$gdd2_1to3 <- df_farm$gdd_1to3*df_farm$gdd_1to3
  df_farm$prec2_1to3 <- df_farm$prec_1to3*df_farm$prec_1to3
  
  #4to10
  df_farm$gdd2_4to10 <- df_farm$gdd4to10*df_farm$gdd4to10
  df_farm$prec2_4to10 <- df_farm$prec_4to10 * df_farm$prec_4to10
  
  #1to5
  df_farm$gdd2_1to5 <- df_farm$gdd1to5*df_farm$gdd1to5
  df_farm$prec2_1to5 <- df_farm$prec_1to5*df_farm$prec_1to5
  
  #6to10
  df_farm$gdd2_6to10 <- df_farm$gdd6to10*df_farm$gdd6to10
  df_farm$prec2_6to10 <- df_farm$prec_6to10 * df_farm$prec_6to10
  
  #6to20
  df_farm$gdd2_6to20 <- df_farm$gdd6to20*df_farm$gdd6to20
  df_farm$prec2_6to20 <- df_farm$prec_6to20 * df_farm$prec_6to20
  
  #1
  df_farm$gdd2_1 <- df_farm$gdd1*df_farm$gdd1
  df_farm$prec2_1 <- df_farm$prec_1 * df_farm$prec_1
  
  #2to5
  df_farm$gdd2_2to5 <- df_farm$gdd2to5*df_farm$gdd2to5
  df_farm$prec2_2to5 <- df_farm$prec_2to5*df_farm$prec_2to5
  
  #2to10
  df_farm$gdd2_2to10 <- df_farm$gdd2to10*df_farm$gdd2to10
  df_farm$prec2_2to10 <- df_farm$prec_2to10*df_farm$prec_2to10
  
  #1to10
  df_farm$gdd2_1to10 <- df_farm$gdd1to10 * df_farm$gdd1to10
  df_farm$prec2_1to10 <- df_farm$prec_1to10 * df_farm$prec_1to10
  
  #1to30 
  df_farm$gdd2_1to30 <- df_farm$gdd1to30 * df_farm$gdd1to30
  df_farm$prec2_1to30 <- df_farm$prec_1to30 * df_farm$prec_1to30

#Interaction terms
df_farm$gdd_obs_1to3 <- df_farm$gdd * df_farm$gdd_1to3
df_farm$gdd2_obs_1to3 <- df_farm$gdd*df_farm$gdd * df_farm$gdd_1to3
df_farm$prec_obs_1to3 <- df_farm$prec * df_farm$prec_1to3
df_farm$prec2_obs_1to3 <- df_farm$prec*df_farm$prec * df_farm$prec_1to3
df_farm$gddHigh_obs_1to3 <- df_farm$gddHigh * df_farm$gddHigh_1to3
df_farm$dd_obs_1to3 <- df_farm$dd * df_farm$dd_1to3

df_farm$gdd_obs_4to10 <- df_farm$gdd * df_farm$gdd_4to10
df_farm$gdd2_obs_4to10 <- df_farm$gdd*df_farm$gdd * df_farm$gdd_4to10
df_farm$prec_obs_4to10 <- df_farm$prec * df_farm$prec_4to10
df_farm$prec2_obs_4to10 <- df_farm$prec*df_farm$prec * df_farm$prec_4to10
df_farm$gddHigh_obs_4to10 <- df_farm$gddHigh * df_farm$gddHigh_4to10
df_farm$dd_obs_4to10 <- df_farm$dd * df_farm$dd_4to10

df_farm$gdd_obs_1to5 <- df_farm$gdd * df_farm$gdd_1to5
df_farm$gdd2_obs_1to5 <- df_farm$gdd*df_farm$gdd * df_farm$gdd_1to5
df_farm$prec_obs_1to5 <- df_farm$prec * df_farm$prec_1to5
df_farm$prec2_obs_1to5 <- df_farm$prec*df_farm$prec * df_farm$prec_1to5
df_farm$gddHigh_obs_1to5 <- df_farm$gddHigh * df_farm$gddHigh_1to5
df_farm$dd_obs_1to5 <- df_farm$dd * df_farm$dd_1to5

df_farm$gdd_obs_6to10 <- df_farm$gdd * df_farm$gdd_6to10
df_farm$gdd2_obs_6to10 <- df_farm$gdd*df_farm$gdd * df_farm$gdd_6to10
df_farm$prec_obs_6to10 <- df_farm$prec * df_farm$prec_6to10
df_farm$prec2_obs_6to10 <- df_farm$prec*df_farm$prec * df_farm$prec_6to10
df_farm$gddHigh_obs_6to10 <- df_farm$gddHigh * df_farm$gddHigh_6to10
df_farm$dd_obs_6to10 <- df_farm$dd * df_farm$dd_6to10

df_farm$gdd_obs_6to20 <- df_farm$gdd * df_farm$gdd_6to20
df_farm$gdd2_obs_6to20 <- df_farm$gdd*df_farm$gdd * df_farm$gdd_6to20
df_farm$prec_obs_6to20 <- df_farm$prec * df_farm$prec_6to20
df_farm$prec2_obs_6to20 <- df_farm$prec*df_farm$prec * df_farm$prec_6to20
df_farm$gddHigh_obs_6to20 <- df_farm$gddHigh * df_farm$gddHigh_6to20
df_farm$dd_obs_6to20 <- df_farm$dd * df_farm$dd_6to20

df_farm$gdd_obs_1 <- df_farm$gdd * df_farm$gdd_1
df_farm$gdd2_obs_1 <- df_farm$gdd*df_farm$gdd * df_farm$gdd_1
df_farm$prec_obs_1 <- df_farm$prec * df_farm$prec_1
df_farm$prec2_obs_1 <- df_farm$prec*df_farm$prec * df_farm$prec_1
df_farm$gddHigh_obs_1 <- df_farm$gddHigh * df_farm$gddHigh_1
df_farm$dd_obs_1 <- df_farm$dd * df_farm$dd__1

df_farm$gdd_obs_2to5 <- df_farm$gdd * df_farm$gdd_2to5
df_farm$gdd2_obs_2to5 <- df_farm$gdd*df_farm$gdd * df_farm$gdd_2to5
df_farm$prec_obs_2to5 <- df_farm$prec * df_farm$prec_2to5
df_farm$prec2_obs_2to5 <- df_farm$prec*df_farm$prec * df_farm$prec_2to5
df_farm$gddHigh_obs_2to5 <- df_farm$gddHigh * df_farm$gddHigh_2to5
df_farm$dd_obs_2to5 <- df_farm$dd * df_farm$dd_2to5

df_farm$gdd_obs_2to10 <- df_farm$gdd * df_farm$gdd_2to10
df_farm$gdd2_obs_2to10 <- df_farm$gdd*df_farm$gdd * df_farm$gdd_2to10
df_farm$prec_obs_2to10 <- df_farm$prec * df_farm$prec_2to10
df_farm$prec2_obs_2to10 <- df_farm$prec*df_farm$prec * df_farm$prec_2to10
df_farm$gddHigh_obs_2to10 <- df_farm$gddHigh * df_farm$gddHigh_2to10
df_farm$dd_obs_2to10 <- df_farm$dd * df_farm$dd_2to10

#---------------------------------#
#### 2) Define price variables ####
#---------------------------------#

# assign regional prices to farms
df_farm$p_cereals <- ifelse(is.na(df_farm$n2lp_cereals), 
                            ifelse(is.na(df_farm$n1lp_cereals), 
                                   df_farm$n0lp_cereals, df_farm$n1lp_cereals), 
                            df_farm$n2lp_cereals)

df_farm$p_oilseed <- ifelse(is.na(df_farm$n2lp_oilseed), 
                            ifelse(is.na(df_farm$n1lp_oilseed), 
                                   df_farm$n0lp_oilseed, df_farm$n1lp_oilseed), 
                            df_farm$n2lp_oilseed)

df_farm$p_roots <- ifelse(is.na(df_farm$n2lp_root), 
                          ifelse(is.na(df_farm$n1lp_root), 
                                 df_farm$n0lp_root, df_farm$n1lp_root), 
                          df_farm$n2lp_root)

df_farm$p_protein <- ifelse(is.na(df_farm$n2lp_protein), 
                            ifelse(is.na(df_farm$n1lp_protein), 
                                   df_farm$n0lp_protein, df_farm$n1lp_protein), 
                            df_farm$n2lp_protein)

df_farm$p_corn <- ifelse(is.na(df_farm$n2lp_corn), 
                            ifelse(is.na(df_farm$n1lp_corn), 
                                   df_farm$n0lp_corn, df_farm$n1lp_corn), 
                            df_farm$n2lp_corn)

df_farm$w_fert <- df_farm$n0w_fert
df_farm$w_otherinp <- df_farm$n2w_tq_otherinp

# normalize prices by numeraire
df_farm$np_cereals <- df_farm$p_cereals / df_farm$w_otherinp
df_farm$np_oilseed <- df_farm$p_oilseed / df_farm$w_otherinp
df_farm$np_roots <- df_farm$p_roots / df_farm$w_otherinp
df_farm$np_protein <- df_farm$p_protein / df_farm$w_otherinp
df_farm$np_corn <- df_farm$p_corn / df_farm$w_otherinp
df_farm$nw_fert <- df_farm$w_fert / df_farm$w_otherinp

#----------------------------------#
#### 3) Define input quantities ####
#----------------------------------#

# Fertilizer (negative value needed for estimation)
df_farm$x_fert <- -df_farm$z2539s02 / df_farm$w_fert # Note: z2539s02 is a negative number
df_farm$nx_fert <- -df_farm$x_fert

# Other variable inputs (except fertilizer)
df_farm$x_otherinp <- (df_farm$c_cropspec + df_farm$c_animspec + 
                         df_farm$c_intermediates) / df_farm$w_otherinp


# Fixed inputs: land, labor, capital
df_farm$k_land <- df_farm$x_land
df_farm$k_labor <- df_farm$x_labor
df_farm$k_capital <- df_farm$x_capital

#-----------------------------------------------#
#### 4) Prepare land shares and lags thereof ####
#-----------------------------------------------#

#land shares
df_farm$aa_crops <- df_farm$aa_cereals + df_farm$aa_protein + df_farm$aa_oilseed + df_farm$aa_roots + df_farm$aa_corn
df_farm$sh_cereals <- df_farm$aa_cereals / df_farm$aa_crops
df_farm$sh_protein <- df_farm$aa_protein / df_farm$aa_crops
df_farm$sh_oilseed <- df_farm$aa_oilseed / df_farm$aa_crops
df_farm$sh_roots <- df_farm$aa_roots / df_farm$aa_crops
df_farm$sh_corn <- df_farm$aa_corn / df_farm$aa_crops

# lagged land shares
myvars <- c("key","year","sh_cereals", "sh_protein", "sh_oilseed", "sh_roots", "sh_corn")
lags <- df_farm[myvars]

lags <- plm::pdata.frame(lags,index=c("key","year")) 
lags$lsh_cereals <- plm::lag(lags$sh_cereals)
lags$lsh_oilseed <- plm::lag(lags$sh_oilseed)
lags$lsh_protein <- plm::lag(lags$sh_protein)
lags$lsh_roots <- plm::lag(lags$sh_roots)
lags$lsh_corn    <- plm::lag(lags$sh_corn)

lags[] <- lapply(lags, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})

df_farm <- left_join(df_farm, 
                     lags)
remove(lags)

#---------------------------------------------#
#### 5) Create dummies for crop production ####
#---------------------------------------------#

# Dummy variables for production / no production
df_farm$iCereals <- ifelse(df_farm$qq_cereals >0,1,0)
df_farm$iProtein <- ifelse(df_farm$qq_protein >0,1,0)
df_farm$iOilseed <- ifelse(df_farm$qq_oilseed >0,1,0)
df_farm$iRoots <- ifelse(df_farm$qq_roots >0,1,0)
df_farm$iCorn <- ifelse(df_farm$qq_corn >0,1,0)

#---------------------------#
#### 6) Construct sample ####
#---------------------------#

# Keep only specialized crop farms

  # Create indicator whether farm specialized crop farm
  df_farm$cropspec <- ifelse(df_farm$z2099s05>0 & # crops
                             (df_farm$z2199s05==0 | is.na(df_farm$z2199s05) ) & # livestock
                             (df_farm$z2209s05==0 | is.na(df_farm$z2209s05) ) & # fruits
                             (df_farm$z2269s05==0 | is.na(df_farm$z2269s05) ) & # horticulture
                             (df_farm$z2299s05==0 | is.na(df_farm$z2299s05) ),1,0) # viniculture

  # Keep only farms that are specialized crop farms in each year  
  df_farm <- df_farm %>% 
      group_by(key) %>% 
      mutate(min=min(cropspec)) %>% 
      filter(min==1) %>% 
      dplyr::select(-min)

  table(df_farm$year)

# Select years from 1996 - 2019 (1995 excluded because of lagged vars)
df_farm <- subset(df_farm, year > 1995) 

# Create the corresponding time trend and its squared term
df_farm$trend <- df_farm$year - 1996
df_farm$trend2 <- df_farm$trend * df_farm$trend

# Consider only farms that produce considered crops in each year
df_farm <- df_farm %>% 
  group_by(key) %>% 
  mutate(min=min(aa_crops)) %>% 
  filter(min>0) %>% 
  dplyr::select(-min)

# Consider only farms that appear at least three times
df_farm <- df_farm %>% 
  group_by(key) %>%
  mutate(count=n()) %>% 
  filter(count>2) 

#-------------------------------------------------------------#
#### 7) Keep only complete observations for the estimation ####
#-------------------------------------------------------------#

# Summarize variables that are required for the analysis
cols <- c("key", "year", "trend", "trend2",
          "qq_cereals", "qq_oilseed", "qq_roots", "qq_protein", "qq_corn", "nx_fert", "x_fert", "x_otherinp",
          "np_cereals", "np_oilseed", "np_roots", "np_protein", "np_corn", "nw_fert", 
          "k_land", "k_labor", "k_capital", 
          "aa_crops", "aa_cereals", "aa_oilseed", "aa_roots", "aa_protein", "aa_corn",
          "iCereals", "iOilseed", "iRoots", "iProtein", "iCorn",
          "sh_cereals", "sh_protein", "sh_oilseed", "sh_roots", "sh_corn",
          "lsh_cereals", "lsh_protein", "lsh_oilseed", "lsh_roots", "lsh_corn",
          "p_cereals", "p_oilseed", "p_roots", "p_protein", "p_corn", "w_fert", "w_otherinp",
          "gdd_obs", "prec_obs", "gddHigh_obs", "dd_obs",
          "gdd_1to3", "prec_1to3", "gddHigh_1to3", "dd_1to3",
          "gdd_4to10",  "prec_4to10", "gddHigh_4to10", "dd_4to10",
          "gdd_1to10",  "prec_1to10", "gddHigh_1to10", "dd_1to10",
          "gdd_1to30",  "prec_1to30", "gddHigh_1to30", "dd_1to30")
summary(df_farm[cols])

# Replace missing values by 0: capital (n=54)
df_farm$k_capital[is.na(df_farm$k_capital)] <- 0 

# Keep only required variables in the data
df_farm <- df_farm %>% 
  dplyr::select(key, year, trend, trend2, 
                qq_cereals, qq_oilseed, qq_roots, qq_protein, qq_corn, nx_fert, x_fert, x_otherinp,
                np_cereals, np_oilseed, np_roots, np_protein, np_corn, nw_fert, 
                k_land, k_labor, k_capital, 
                aa_crops, aa_cereals, aa_oilseed, aa_roots, aa_protein, aa_corn,
                sh_cereals, sh_protein, sh_oilseed, sh_roots, sh_corn,
                lsh_cereals, lsh_protein, lsh_oilseed, lsh_roots, lsh_corn,
                iCereals, iOilseed, iRoots, iProtein, iCorn,
                p_cereals, p_oilseed, p_roots, p_protein, p_corn, w_fert, w_otherinp,
                gdd_obs, gdd2_obs, prec_obs, prec2_obs, gddHigh_obs, dd_obs,
                gdd_1to3, gdd2_1to3, prec_1to3, prec2_1to3, gddHigh_1to3, dd_1to3,
                gdd_4to10, gdd2_4to10,  prec_4to10, prec2_4to10, gddHigh_4to10, dd_4to10,
                gdd_1to5, gdd2_1to5, prec_1to5, prec2_1to5, gddHigh_1to5, dd_1to5,
                gdd_6to10, gdd2_6to10,  prec_6to10, prec2_6to10, gddHigh_6to10, dd_6to10,
                gdd_6to20, gdd2_6to20,  prec_6to20, prec2_6to20, gddHigh_6to20, dd_6to20,
                gdd_1, gdd2_1, prec_1, prec2_1, gddHigh_1, dd_1,
                gdd_2to5, gdd2_2to5,  prec_2to5, prec2_2to5, gddHigh_2to5, dd_2to5,
                gdd_2to10, gdd2_2to10, prec_2to10, prec2_2to10, gddHigh_2to10, dd_2to10,
                gdd_1to10, gdd2_1to10, prec_1to10, prec2_1to10, gddHigh_1to10, dd_1to10,
                gdd_1to30, gdd2_1to30, prec_1to30, prec2_1to30,gddHigh_1to30, dd_1to30,
                gdd_obs_1to3, gdd2_obs_1to3, prec_obs_1to3, prec2_obs_1to3, gddHigh_obs_1to3, dd_obs_1to3, 
                gdd_obs_4to10, gdd2_obs_4to10, prec_obs_4to10, prec2_obs_4to10, gddHigh_obs_4to10, dd_obs_4to10,
                gdd_obs_1to5, gdd2_obs_1to5, prec_obs_1to5, prec2_obs_1to5, gddHigh_obs_1to5, dd_obs_1to5, 
                gdd_obs_6to10, gdd2_obs_6to10, prec_obs_6to10, prec2_obs_6to10, gddHigh_obs_6to10, dd_obs_6to10,
                gdd_obs_6to20, gdd2_obs_6to20, prec_obs_6to20, prec2_obs_6to20, gddHigh_obs_6to20, dd_obs_6to20,
                gdd_obs_1, gdd2_obs_1, prec_obs_1, prec2_obs_1, gddHigh_obs_1, dd_obs_1, 
                gdd_obs_2to10, gdd2_obs_2to10, prec_obs_2to10, prec2_obs_2to10, gddHigh_obs_2to10, dd_obs_2to10, 
                nuts1, nuts2, nuts3, AGS_total)

df_farm$complete <- complete.cases(df_farm)
df_farm$complete <- as.integer(as.logical(df_farm$complete))
summary(df_farm$complete)

# drop the remaining obs with missing values (they arise from the lagged shares)
df_farm <- df_farm %>% 
  filter(complete==1)

# There are now 24 farms with only 1 obs; This is because they are not observed in
# all years (i.e., data gaps) and hence have missing lag variables although they 
# are in the data >2 times. They need to be dropped: 
df_farm <- df_farm %>% 
  group_by(key) %>%
  mutate(count=n()) 

df_farm <- df_farm %>% 
  filter(count>1) 

#-------------------------------------------#
#### 8) Calculate descriptive statistics ####
#-------------------------------------------#

#-------------------------------#
# Describe nature of panel data #
#-------------------------------#

# Number of distinct farms
length(unique(df_farm$key))

# First year of sample appearance
df_firstyear <- df_farm %>% 
  group_by(key) %>% 
  filter(year == min(year)) %>% 
  select(key,
         first_appearenance=year)
table(df_firstyear$first_appearenance)

# Last year of sample appearance
df_lastyear <- df_farm %>% 
  group_by(key) %>% 
  filter(year == max(year)) %>% 
  select(key,
         last_appearenance=year)
table(df_lastyear$last_appearenance)

#------------------------------#
# Within- and between variance #
#------------------------------#

# Preparation
list_farms <- as.data.frame(df_farm) %>% dplyr::select(key) %>% distinct() # list of farms
n_farms <- list_farms %>% count() %>% as.numeric() # number of farms
list_years <- as.data.frame(df_farm) %>% dplyr::select(year) %>% arrange(year) %>% distinct() # list of years
n_years <- list_years %>% count() %>% as.numeric() # number of years

# Write function
fun.variability <- function(x) {
  
  str.x_between <- vector(,length=nrow(list_farms))
  for (i in 1:n_farms) { 
    cur_farm <- list_farms[i,1] 
    str.x_between[i] <- ( mean(x[df_farm$key==cur_farm]) - mean(x) )^2
  }
  x_between <- sum(str.x_between) / (n_farms - 1)
  
  mat.x_within = matrix(, nrow = nrow(list_farms), ncol = nrow(list_years))
  for (i in 1:n_farms) {
    for (t in 1:n_years) {
      cur_farm <- list_farms[i,1] 
      cur_year    <- list_years[t,1] 
      if (identical(x[df_farm$key==cur_farm & df_farm$year==cur_year], numeric(0)) ) {
        mat.x_within [i,t] <- NA
      }
      else {
        mat.x_within[i,t] <- ( x[df_farm$key==cur_farm & df_farm$year==cur_year] - mean(x[df_farm$key==cur_farm]) )^2
      }
    }
  } 
    x_within <- sum(mat.x_within, na.rm=TRUE) / (n_farms*(length(unique(df_farm$year))-1))
    
    result <- list()
    result$x_mean <- mean(x)
    result$x_sd <- sd(x)
    result$x_min <- min(x)
    result$x_max <- max(x)
    result$x_p10 <- quantile(x, probs=0.10)
    result$x_p90 <- quantile(x, probs=0.90)
    result$x_between <- x_between
    result$x_within <- x_within
    result$x_mean_1996 <- mean(x[df_farm$year==1996])
    result$x_sd_1996 <- sd(x[df_farm$year==1996])
    result$x_min_1996 <- min(x[df_farm$year==1996])
    result$x_max_1996 <- max(x[df_farm$year==1996])
    result$x_p10_1996 <- quantile(x[df_farm$year==1996], probs=0.10)
    result$x_p90_1996 <- quantile(x[df_farm$year==1996], probs=0.90)
    result$x_mean_2008  <- mean(x[df_farm$year==2008])
    result$x_sd_2008 <- sd(x[df_farm$year==2008])
    result$x_min_2008 <- min(x[df_farm$year==2008])
    result$x_max_2008 <- max(x[df_farm$year==2008])
    result$x_p10_2008 <- quantile(x[df_farm$year==2008], probs=0.10)
    result$x_p90_2008 <- quantile(x[df_farm$year==2008], probs=0.90)
    result$x_mean_2019  <- mean(x[df_farm$year==2019])
    result$x_sd_2019 <- sd(x[df_farm$year==2019])
    result$x_min_2019 <- min(x[df_farm$year==2019])
    result$x_max_2019 <- max(x[df_farm$year==2019])
    result$x_p10_2019 <- quantile(x[df_farm$year==2019], probs=0.10)
    result$x_p90_2019 <- quantile(x[df_farm$year==2019], probs=0.90)
    return(result)
  }
  
# Calculate stats using the new function
stats_qq_cereals <- fun.variability(df_farm$qq_cereals)
stats_qq_oilseed <- fun.variability(df_farm$qq_oilseed)
stats_qq_roots <- fun.variability(df_farm$qq_roots)
stats_qq_protein <- fun.variability(df_farm$qq_protein)
stats_qq_corn <- fun.variability(df_farm$qq_corn)
stats_x_fert <- fun.variability(df_farm$x_fert)
stats_x_otherinp <- fun.variability(df_farm$x_otherinp)
stats_p_cereals <- fun.variability(df_farm$p_cereals)
stats_p_oilseed <- fun.variability(df_farm$p_oilseed)
stats_p_roots <- fun.variability(df_farm$p_roots)
stats_p_protein <- fun.variability(df_farm$p_protein)
stats_p_corn <- fun.variability(df_farm$p_corn)
stats_w_fert <- fun.variability(df_farm$w_fert)
stats_w_otherinp <- fun.variability(df_farm$w_otherinp)
stats_iCereals <- fun.variability(df_farm$iCereals)
stats_iOilseed <- fun.variability(df_farm$iOilseed)
stats_iRoots <- fun.variability(df_farm$iRoots)
stats_iProtein <- fun.variability(df_farm$iProtein)
stats_iCorn <- fun.variability(df_farm$iCorn)
stats_sh_cereals <- fun.variability(df_farm$sh_cereals)
stats_sh_oilseed <- fun.variability(df_farm$sh_oilseed)
stats_sh_roots <- fun.variability(df_farm$sh_roots)
stats_sh_protein <- fun.variability(df_farm$sh_protein)
stats_sh_corn <- fun.variability(df_farm$sh_corn)
stats_k_land <- fun.variability(df_farm$k_land)
stats_k_labor <- fun.variability(df_farm$k_labor)
stats_k_capital <- fun.variability(df_farm$k_capital)
stats_gdd_obs <- fun.variability(df_farm$gdd_obs)
stats_prec_obs <- fun.variability(df_farm$prec_obs)
stats_gddHigh_obs <- fun.variability(df_farm$gddHigh_obs)
stats_dd_obs <- fun.variability(df_farm$dd_obs)
stats_gdd_1to3 <- fun.variability(df_farm$gdd_1to3)
stats_prec_1to3 <- fun.variability(df_farm$prec_1to3)
stats_gddHigh_1to3 <- fun.variability(df_farm$gddHigh_1to3)
stats_dd_1to3 <- fun.variability(df_farm$dd_1to3)
stats_gdd_4to10 <- fun.variability(df_farm$gdd_4to10)
stats_prec_4to10 <- fun.variability(df_farm$prec_4to10)
stats_gddHigh_4to10 <- fun.variability(df_farm$gddHigh_4to10)
stats_dd_4to10 <- fun.variability(df_farm$dd_4to10)

# Save stats in a list
stats <- list(stats_qq_cereals,
              stats_qq_oilseed,
              stats_qq_roots,
              stats_qq_protein,
              stats_qq_corn,
              stats_x_fert,
              stats_x_otherinp,
              stats_p_cereals,
              stats_p_oilseed,
              stats_p_roots,
              stats_p_protein,
              stats_p_corn,
              stats_w_fert,
              stats_w_otherinp,
              stats_iCereals,
              stats_iOilseed,
              stats_iRoots,
              stats_iProtein,
              stats_iCorn,
              stats_sh_cereals,
              stats_sh_oilseed,
              stats_sh_roots,
              stats_sh_protein,
              stats_sh_corn,
              stats_k_land,
              stats_k_labor,
              stats_k_capital,
              stats_gdd_obs,
              stats_prec_obs,
              stats_gddHigh_obs,
              stats_dd_obs,
              stats_gdd_1to3,
              stats_prec_1to3,
              stats_gddHigh_1to3,
              stats_dd_1to3,
              stats_gdd_4to10,
              stats_prec_4to10,
              stats_gddHigh_4to10,
              stats_dd_4to10)

names(stats) <- c("stats_qq_cereals",
                  "stats_qq_oilseed",
                  "stats_qq_roots",
                  "stats_qq_protein",
                  "stats_qq_corn",
                  "stats_x_fert",
                  "stats_x_otherinp",
                  "stats_p_cereals",
                  "stats_p_oilseed",
                  "stats_p_roots",
                  "stats_p_protein",
                  "stats_p_corn",
                  "stats_w_fert",
                  "stats_w_otherinp",
                  "stats_iCereals",
                  "stats_iOilseed",
                  "stats_iRoots",
                  "stats_iProtein",
                  "stats_iCorn",
                  "stats_sh_cereals",
                  "stats_sh_oilseed",
                  "stats_sh_roots",
                  "stats_sh_protein",
                  "stats_sh_corn",
                  "stats_k_land",
                  "stats_k_labor",
                  "stats_k_capital",
                  "stats_gdd_obs",
                  "stats_prec_obs",
                  "stats_gddHigh_obs",
                  "stats_dd_obs",
                  "stats_gdd_1to3",
                  "stats_prec_1to3",
                  "stats_gddHigh_1to3",
                  "stats_dd_1to3",
                  "stats_gdd_4to10",
                  "stats_prec_4to10",
                  "stats_gddHigh_4to10",
                  "stats_dd_4to10")

# ----------------------- #
# Table for the main text #
# ----------------------- #

Tab_stats_main <- data.frame(
  "Variable"=c("Cereals quantity", "Protein crops quantity", "Oilseeds quantity", "Root crops quantity", "Corn quantity", "Fertilizer quantity", "Other input quantity",
               "Cereals price",  "Protein crops price", "Oilseeds price", "Root crops price", "Corn price", "Fertilizer price", "Other input price",
               "Cereals > 0", "Protein > 0", "Oilseeds > 0", "Root crops > 0", "Corn >0",
               "Share cereals", "Share protein", "Share oilseeds", "Share root crops", "Share corn",
               "Land", "Labor", "Capital", 
               "Gdd_obs", "Prec_obs", "GddHigh_obs",  "DD_obs", 
               "Gdd_1to3", "Prec_1to3", "GddHigh_1to3",  "DD_1to3",
               "Gdd_4to10", "Prec_4to10", "GddHigh_4to10",  "DD_4to10"),
  "Mean" = c(stats$stats_qq_cereals$x_mean,
             stats$stats_qq_protein$x_mean,
             stats$stats_qq_oilseed$x_mean,
             stats$stats_qq_roots$x_mean,
             stats$stats_qq_corn$x_mean,
             stats$stats_x_fert$x_mean,
             stats$stats_x_otherinp$x_mean,
             stats$stats_p_cereals$x_mean,
             stats$stats_p_protein$x_mean,
             stats$stats_p_oilseed$x_mean,
             stats$stats_p_roots$x_mean,
             stats$stats_p_corn$x_mean,
             stats$stats_w_fert$x_mean,
             stats$stats_w_otherinp$x_mean,
             stats$stats_iCereals$x_mean,
             stats$stats_iProtein$x_mean,
             stats$stats_iOilseed$x_mean,
             stats$stats_iRoots$x_mean,
             stats$stats_iCorn$x_mean,
             stats$stats_sh_cereals$x_mean,
             stats$stats_sh_protein$x_mean,
             stats$stats_sh_oilseed$x_mean,
             stats$stats_sh_roots$x_mean,
             stats$stats_sh_corn$x_mean,
             stats$stats_k_land$x_mean,
             stats$stats_k_labor$x_mean,
             stats$stats_k_capital$x_mean,
             stats$stats_gdd_obs$x_mean,
             stats$stats_prec_obs$x_mean,
             stats$stats_gddHigh_obs$x_mean,
             stats$stats_dd_obs$x_mean,
             stats$stats_gdd_1to3$x_mean,
             stats$stats_prec_1to3$x_mean,
             stats$stats_gddHigh_1to3$x_mean,
             stats$stats_dd_1to3$x_mean,
             stats$stats_gdd_4to10$x_mean,
             stats$stats_prec_4to10$x_mean,
             stats$stats_gddHigh_4to10$x_mean,
             stats$stats_dd_4to10$x_mean),
  "Std. dev." = c(stats$stats_qq_cereals$x_sd,
                  stats$stats_qq_protein$x_sd,
                  stats$stats_qq_oilseed$x_sd,
                  stats$stats_qq_roots$x_sd,
                  stats$stats_qq_corn$x_sd,
                  stats$stats_x_fert$x_sd,
                  stats$stats_x_otherinp$x_sd,
                  stats$stats_p_cereals$x_sd,
                  stats$stats_p_protein$x_sd,
                  stats$stats_p_oilseed$x_sd,
                  stats$stats_p_roots$x_sd,
                  stats$stats_p_corn$x_sd,
                  stats$stats_w_fert$x_sd,
                  stats$stats_w_otherinp$x_sd,
                  stats$stats_iCereals$x_sd,
                  stats$stats_iProtein$x_sd,
                  stats$stats_iOilseed$x_sd,
                  stats$stats_iRoots$x_sd,
                  stats$stats_iCorn$x_sd,
                  stats$stats_sh_cereals$x_sd,
                  stats$stats_sh_protein$x_sd,
                  stats$stats_sh_oilseed$x_sd,
                  stats$stats_sh_roots$x_sd,
                  stats$stats_sh_corn$x_sd,
                  stats$stats_k_land$x_sd,
                  stats$stats_k_labor$x_sd,
                  stats$stats_k_capital$x_sd,
                  stats$stats_gdd_obs$x_sd,
                  stats$stats_prec_obs$x_sd,
                  stats$stats_gddHigh_obs$x_sd,
                  stats$stats_dd_obs$x_sd,
                  stats$stats_gdd_1to3$x_sd,
                  stats$stats_prec_1to3$x_sd,
                  stats$stats_gddHigh_1to3$x_sd,
                  stats$stats_dd_1to3$x_sd,
                  stats$stats_gdd_4to10$x_sd,
                  stats$stats_prec_4to10$x_sd,
                  stats$stats_gddHigh_4to10$x_sd,
                  stats$stats_dd_4to10$x_sd),
  "Within-std. dev." = sqrt(c(stats$stats_qq_cereals$x_within,
                              stats$stats_qq_protein$x_within,
                              stats$stats_qq_oilseed$x_within,
                              stats$stats_qq_roots$x_within,
                              stats$stats_qq_corn$x_within,
                              stats$stats_x_fert$x_within,
                              stats$stats_x_otherinp$x_within,
                              stats$stats_p_cereals$x_within,
                              stats$stats_p_protein$x_within,
                              stats$stats_p_oilseed$x_within,
                              stats$stats_p_roots$x_within,
                              stats$stats_p_corn$x_within,
                              stats$stats_w_fert$x_within,
                              stats$stats_w_otherinp$x_within,
                              stats$stats_iCereals$x_within,
                              stats$stats_iProtein$x_within,
                              stats$stats_iOilseed$x_within,
                              stats$stats_iRoots$x_within,
                              stats$stats_iCorn$x_within,
                              stats$stats_sh_cereals$x_within,
                              stats$stats_sh_protein$x_within,
                              stats$stats_sh_oilseed$x_within,
                              stats$stats_sh_roots$x_within,
                              stats$stats_sh_corn$x_within,
                              stats$stats_k_land$x_within,
                              stats$stats_k_labor$x_within,
                              stats$stats_k_capital$x_within,
                              stats$stats_gdd_obs$x_within,
                              stats$stats_prec_obs$x_within,
                              stats$stats_gddHigh_obs$x_within,
                              stats$stats_dd_obs$x_within,
                              stats$stats_gdd_1to3$x_within,
                              stats$stats_prec_1to3$x_within,
                              stats$stats_gddHigh_1to3$x_within,
                              stats$stats_dd_1to3$x_within,
                              stats$stats_gdd_4to10$x_within,
                              stats$stats_prec_4to10$x_within,
                              stats$stats_gddHigh_4to10$x_within,
                              stats$stats_dd_4to10$x_within)),
  "Between-std. devt." = sqrt(c(stats$stats_qq_cereals$x_between,
                                stats$stats_qq_protein$x_between,
                                stats$stats_qq_oilseed$x_between,
                                stats$stats_qq_roots$x_between,
                                stats$stats_qq_corn$x_between,
                                stats$stats_x_fert$x_between,
                                stats$stats_x_otherinp$x_between,
                                stats$stats_p_cereals$x_between,
                                stats$stats_p_protein$x_between,
                                stats$stats_p_oilseed$x_between,
                                stats$stats_p_roots$x_between,
                                stats$stats_p_corn$x_between,
                                stats$stats_w_fert$x_between,
                                stats$stats_w_otherinp$x_between,
                                stats$stats_iCereals$x_between,
                                stats$stats_iProtein$x_between,
                                stats$stats_iOilseed$x_between,
                                stats$stats_iRoots$x_between,
                                stats$stats_iCorn$x_between,
                                stats$stats_sh_cereals$x_between,
                                stats$stats_sh_protein$x_between,
                                stats$stats_sh_oilseed$x_between,
                                stats$stats_sh_roots$x_between,
                                stats$stats_sh_corn$x_between,
                                stats$stats_k_land$x_between,
                                stats$stats_k_labor$x_between,
                                stats$stats_k_capital$x_between,
                                stats$stats_gdd_obs$x_between,
                                stats$stats_prec_obs$x_between,
                                stats$stats_gddHigh_obs$x_between,
                                stats$stats_dd_obs$x_between,
                                stats$stats_gdd_1to3$x_between,
                                stats$stats_prec_1to3$x_between,
                                stats$stats_gddHigh_1to3$x_between,
                                stats$stats_dd_1to3$x_between,
                                stats$stats_gdd_4to10$x_between,
                                stats$stats_prec_4to10$x_between,
                                stats$stats_gddHigh_4to10$x_between,
                                stats$stats_dd_4to10$x_between))
)

#Numbers smaller than 0.0005 should be in scientific format; others in 0.001 etc.
Tab_stats_main <- as.data.frame(cbind(Tab_stats_main[,1], sapply(Tab_stats_main[,2:5], function(x) ifelse(abs(x)<0.0005&is.numeric(x),format(x,scientific=TRUE,digits=2),format(round(x, 3), nsmall = 3, scientific=FALSE)) )) )
Tab_stats_main

#Write to excel
write_xlsx(Tab_stats_main,"Tables/Table_1.xlsx")

# ---------------------- #
# Table for the appendix #
# ---------------------- #

Tab_stats_app <- data.frame(
  "Variable" = c("Cereals quantity", "", "Protein crops quantity", "", "Oilseeds quantity", "", "Root crops quantity", "", "Corn quantity", "", "Fertilizer quantity", "", "Other input quantity", "",
                 "Cereals price",  "", "Protein crops price", "", "Oilseeds price", "", "Root crops price", "", "Corn price", "", "Fertilizer price", "", "Other input price", "",
                 "Cereals > 0", "", "Protein > 0", "", "Oilseeds > 0", "", "Root crops > 0", "", "Corn >0", "",
                 "Share cereals", "", "Share protein", "", "Share oilseeds", "", "Share root crops", "", "Share corn", "",
                 "Land", "", "Labor", "", "Capital", "", 
                 "Gdd_obs", "", "Prec_obs", "", "GddHigh_obs", "",  "DD_obs", "", 
                 "Gdd_1to3", "", "Prec_1to3", "", "GddHigh_1to3", "",  "DD_1to3", "", 
                 "Gdd_4to10", "", "Prec_4to10", "", "GddHigh_4to10", "",  "DD_4to10", ""),
  "Average" = c(stats$stats_qq_cereals$x_mean, stats$stats_qq_cereals$x_sd,
                stats$stats_qq_protein$x_mean, stats$stats_qq_protein$x_sd,
                stats$stats_qq_oilseed$x_mean, stats$stats_qq_oilseed$x_sd,
                stats$stats_qq_roots$x_mean, stats$stats_qq_roots$x_sd,
                stats$stats_qq_corn$x_mean, stats$stats_qq_corn$x_sd,
                stats$stats_x_fert$x_mean, stats$stats_x_fert$x_sd,
                stats$stats_x_otherinp$x_mean, stats$stats_x_otherinp$x_sd,
                stats$stats_p_cereals$x_mean, stats$stats_p_cereals$x_sd,
                stats$stats_p_protein$x_mean, stats$stats_p_protein$x_sd,
                stats$stats_p_oilseed$x_mean, stats$stats_p_oilseed$x_sd,
                stats$stats_p_roots$x_mean, stats$stats_p_roots$x_sd,
                stats$stats_p_corn$x_mean, stats$stats_p_corn$x_sd,
                stats$stats_w_fert$x_mean, stats$stats_w_fert$x_sd,
                stats$stats_w_otherinp$x_mean, stats$stats_w_otherinp$x_sd,
                stats$stats_iCereals$x_mean, stats$stats_iCereals$x_sd,
                stats$stats_iProtein$x_mean, stats$stats_iProtein$x_sd,
                stats$stats_iOilseed$x_mean, stats$stats_iOilseed$x_sd,
                stats$stats_iRoots$x_mean, stats$stats_iRoots$x_sd,
                stats$stats_iCorn$x_mean, stats$stats_iCorn$x_sd,
                stats$stats_sh_cereals$x_mean, stats$stats_sh_cereals$x_sd,
                stats$stats_sh_protein$x_mean, stats$stats_sh_protein$x_sd,
                stats$stats_sh_oilseed$x_mean, stats$stats_sh_oilseed$x_sd,
                stats$stats_sh_roots$x_mean, stats$stats_sh_roots$x_sd,
                stats$stats_sh_corn$x_mean, stats$stats_sh_corn$x_sd,
                stats$stats_k_land$x_mean, stats$stats_k_land$x_sd,
                stats$stats_k_labor$x_mean, stats$stats_k_labor$x_sd,
                stats$stats_k_capital$x_mean, stats$stats_k_capital$x_sd,
                stats$stats_gdd_obs$x_mean, stats$stats_gdd_obs$x_sd,
                stats$stats_prec_obs$x_mean, stats$stats_prec_obs$x_sd,
                stats$stats_gddHigh_obs$x_mean, stats$stats_gddHigh_obs$x_sd,
                stats$stats_dd_obs$x_mean, stats$stats_dd_obs$x_sd,
                stats$stats_gdd_1to3$x_mean, stats$stats_gdd_1to3$x_sd,
                stats$stats_prec_1to3$x_mean, stats$stats_prec_1to3$x_sd,
                stats$stats_gddHigh_1to3$x_mean, stats$stats_gddHigh_1to3$x_sd,
                stats$stats_dd_1to3$x_mean, stats$stats_dd_1to3$x_sd,
                stats$stats_gdd_4to10$x_mean, stats$stats_gdd_4to10$x_sd,
                stats$stats_prec_4to10$x_mean, stats$stats_prec_4to10$x_sd,
                stats$stats_gddHigh_4to10$x_mean, stats$stats_gddHigh_4to10$x_sd,
                stats$stats_dd_4to10$x_mean, stats$stats_dd_4to10$x_sd),
  "Year 1996" = c(stats$stats_qq_cereals$x_mean_1996, stats$stats_qq_cereals$x_sd_1996,
                  stats$stats_qq_protein$x_mean_1996, stats$stats_qq_protein$x_sd_1996,
                  stats$stats_qq_oilseed$x_mean_1996, stats$stats_qq_oilseed$x_sd_1996,
                  stats$stats_qq_roots$x_mean_1996, stats$stats_qq_roots$x_sd_1996,
                  stats$stats_qq_corn$x_mean_1996, stats$stats_qq_corn$x_sd_1996,
                  stats$stats_x_fert$x_mean_1996, stats$stats_x_fert$x_sd_1996,
                  stats$stats_x_otherinp$x_mean_1996, stats$stats_x_otherinp$x_sd_1996,
                  stats$stats_p_cereals$x_mean_1996, stats$stats_p_cereals$x_sd_1996,
                  stats$stats_p_protein$x_mean_1996, stats$stats_p_protein$x_sd_1996,
                  stats$stats_p_oilseed$x_mean_1996, stats$stats_p_oilseed$x_sd_1996,
                  stats$stats_p_roots$x_mean_1996, stats$stats_p_roots$x_sd_1996,
                  stats$stats_p_corn$x_mean_1996, stats$stats_p_corn$x_sd_1996,
                  stats$stats_w_fert$x_mean_1996, stats$stats_w_fert$x_sd_1996,
                  stats$stats_w_otherinp$x_mean_1996, stats$stats_w_otherinp$x_sd_1996,
                  stats$stats_iCereals$x_mean_1996, stats$stats_iCereals$x_sd_1996,
                  stats$stats_iProtein$x_mean_1996, stats$stats_iProtein$x_sd_1996,
                  stats$stats_iOilseed$x_mean_1996, stats$stats_iOilseed$x_sd_1996,
                  stats$stats_iRoots$x_mean_1996, stats$stats_iRoots$x_sd_1996,
                  stats$stats_iCorn$x_mean_1996, stats$stats_iCorn$x_sd_1996,
                  stats$stats_sh_cereals$x_mean_1996, stats$stats_sh_cereals$x_sd_1996,
                  stats$stats_sh_protein$x_mean_1996, stats$stats_sh_protein$x_sd_1996,
                  stats$stats_sh_oilseed$x_mean_1996, stats$stats_sh_oilseed$x_sd_1996,
                  stats$stats_sh_roots$x_mean_1996, stats$stats_sh_roots$x_sd_1996,
                  stats$stats_sh_corn$x_mean_1996, stats$stats_sh_corn$x_sd_1996,
                  stats$stats_k_land$x_mean_1996, stats$stats_k_land$x_sd_1996,
                  stats$stats_k_labor$x_mean_1996, stats$stats_k_labor$x_sd_1996,
                  stats$stats_k_capital$x_mean_1996, stats$stats_k_capital$x_sd_1996,
                  stats$stats_gdd_obs$x_mean_1996, stats$stats_gdd_obs$x_sd_1996,
                  stats$stats_prec_obs$x_mean_1996, stats$stats_prec_obs$x_sd_1996,
                  stats$stats_gddHigh_obs$x_mean_1996, stats$stats_gddHigh_obs$x_sd_1996,
                  stats$stats_dd_obs$x_mean_1996, stats$stats_dd_obs$x_sd_1996,
                  stats$stats_gdd_1to3$x_mean_1996, stats$stats_gdd_1to3$x_sd_1996,
                  stats$stats_prec_1to3$x_mean_1996, stats$stats_prec_1to3$x_sd_1996,
                  stats$stats_gddHigh_1to3$x_mean_1996, stats$stats_gddHigh_1to3$x_sd_1996,
                  stats$stats_dd_1to3$x_mean_1996, stats$stats_dd_1to3$x_sd_1996,
                  stats$stats_gdd_4to10$x_mean_1996, stats$stats_gdd_4to10$x_sd_1996,
                  stats$stats_prec_4to10$x_mean_1996, stats$stats_prec_4to10$x_sd_1996,
                  stats$stats_gddHigh_4to10$x_mean_1996, stats$stats_gddHigh_4to10$x_sd_1996,
                  stats$stats_dd_4to10$x_mean_1996, stats$stats_dd_4to10$x_sd_1996),
  "Year 2008" = c(stats$stats_qq_cereals$x_mean_2008, stats$stats_qq_cereals$x_sd_2008,
                  stats$stats_qq_protein$x_mean_2008, stats$stats_qq_protein$x_sd_2008,
                  stats$stats_qq_oilseed$x_mean_2008, stats$stats_qq_oilseed$x_sd_2008,
                  stats$stats_qq_roots$x_mean_2008, stats$stats_qq_roots$x_sd_2008,
                  stats$stats_qq_corn$x_mean_2008, stats$stats_qq_corn$x_sd_2008,
                  stats$stats_x_fert$x_mean_2008, stats$stats_x_fert$x_sd_2008,
                  stats$stats_x_otherinp$x_mean_2008, stats$stats_x_otherinp$x_sd_2008,
                  stats$stats_p_cereals$x_mean_2008, stats$stats_p_cereals$x_sd_2008,
                  stats$stats_p_protein$x_mean_2008, stats$stats_p_protein$x_sd_2008,
                  stats$stats_p_oilseed$x_mean_2008, stats$stats_p_oilseed$x_sd_2008,
                  stats$stats_p_roots$x_mean_2008, stats$stats_p_roots$x_sd_2008,
                  stats$stats_p_corn$x_mean_2008, stats$stats_p_corn$x_sd_2008,
                  stats$stats_w_fert$x_mean_2008, stats$stats_w_fert$x_sd_2008,
                  stats$stats_w_otherinp$x_mean_2008, stats$stats_w_otherinp$x_sd_2008,
                  stats$stats_iCereals$x_mean_2008, stats$stats_iCereals$x_sd_2008,
                  stats$stats_iProtein$x_mean_2008, stats$stats_iProtein$x_sd_2008,
                  stats$stats_iOilseed$x_mean_2008, stats$stats_iOilseed$x_sd_2008,
                  stats$stats_iRoots$x_mean_2008, stats$stats_iRoots$x_sd_2008,
                  stats$stats_iCorn$x_mean_2008, stats$stats_iCorn$x_sd_2008,
                  stats$stats_sh_cereals$x_mean_2008, stats$stats_sh_cereals$x_sd_2008,
                  stats$stats_sh_protein$x_mean_2008, stats$stats_sh_protein$x_sd_2008,
                  stats$stats_sh_oilseed$x_mean_2008, stats$stats_sh_oilseed$x_sd_2008,
                  stats$stats_sh_roots$x_mean_2008, stats$stats_sh_roots$x_sd_2008,
                  stats$stats_sh_corn$x_mean_2008, stats$stats_sh_corn$x_sd_2008,
                  stats$stats_k_land$x_mean_2008, stats$stats_k_land$x_sd_2008,
                  stats$stats_k_labor$x_mean_2008, stats$stats_k_labor$x_sd_2008,
                  stats$stats_k_capital$x_mean_2008, stats$stats_k_capital$x_sd_2008,
                  stats$stats_gdd_obs$x_mean_2008, stats$stats_gdd_obs$x_sd_2008,
                  stats$stats_prec_obs$x_mean_2008, stats$stats_prec_obs$x_sd_2008,
                  stats$stats_gddHigh_obs$x_mean_2008, stats$stats_gddHigh_obs$x_sd_2008,
                  stats$stats_dd_obs$x_mean_2008, stats$stats_dd_obs$x_sd_2008,
                  stats$stats_gdd_1to3$x_mean_2008, stats$stats_gdd_1to3$x_sd_2008,
                  stats$stats_prec_1to3$x_mean_2008, stats$stats_prec_1to3$x_sd_2008,
                  stats$stats_gddHigh_1to3$x_mean_2008, stats$stats_gddHigh_1to3$x_sd_2008,
                  stats$stats_dd_1to3$x_mean_2008, stats$stats_dd_1to3$x_sd_2008,
                  stats$stats_gdd_4to10$x_mean_2008, stats$stats_gdd_4to10$x_sd_2008,
                  stats$stats_prec_4to10$x_mean_2008, stats$stats_prec_4to10$x_sd_2008,
                  stats$stats_gddHigh_4to10$x_mean_2008, stats$stats_gddHigh_4to10$x_sd_2008,
                  stats$stats_dd_4to10$x_mean_2008, stats$stats_dd_4to10$x_sd_2008),
  "Year 2019" = c(stats$stats_qq_cereals$x_mean_2019, stats$stats_qq_cereals$x_sd_2019,
                  stats$stats_qq_protein$x_mean_2019, stats$stats_qq_protein$x_sd_2019,
                  stats$stats_qq_oilseed$x_mean_2019, stats$stats_qq_oilseed$x_sd_2019,
                  stats$stats_qq_roots$x_mean_2019, stats$stats_qq_roots$x_sd_2019,
                  stats$stats_qq_corn$x_mean_2019, stats$stats_qq_corn$x_sd_2019,
                  stats$stats_x_fert$x_mean_2019, stats$stats_x_fert$x_sd_2019,
                  stats$stats_x_otherinp$x_mean_2019, stats$stats_x_otherinp$x_sd_2019,
                  stats$stats_p_cereals$x_mean_2019, stats$stats_p_cereals$x_sd_2019,
                  stats$stats_p_protein$x_mean_2019, stats$stats_p_protein$x_sd_2019,
                  stats$stats_p_oilseed$x_mean_2019, stats$stats_p_oilseed$x_sd_2019,
                  stats$stats_p_roots$x_mean_2019, stats$stats_p_roots$x_sd_2019,
                  stats$stats_p_corn$x_mean_2019, stats$stats_p_corn$x_sd_2019,
                  stats$stats_w_fert$x_mean_2019, stats$stats_w_fert$x_sd_2019,
                  stats$stats_w_otherinp$x_mean_2019, stats$stats_w_otherinp$x_sd_2019,
                  stats$stats_iCereals$x_mean_2019, stats$stats_iCereals$x_sd_2019,
                  stats$stats_iProtein$x_mean_2019, stats$stats_iProtein$x_sd_2019,
                  stats$stats_iOilseed$x_mean_2019, stats$stats_iOilseed$x_sd_2019,
                  stats$stats_iRoots$x_mean_2019, stats$stats_iRoots$x_sd_2019,
                  stats$stats_iCorn$x_mean_2019, stats$stats_iCorn$x_sd_2019,
                  stats$stats_sh_cereals$x_mean_2019, stats$stats_sh_cereals$x_sd_2019,
                  stats$stats_sh_protein$x_mean_2019, stats$stats_sh_protein$x_sd_2019,
                  stats$stats_sh_oilseed$x_mean_2019, stats$stats_sh_oilseed$x_sd_2019,
                  stats$stats_sh_roots$x_mean_2019, stats$stats_sh_roots$x_sd_2019,
                  stats$stats_sh_corn$x_mean_2019, stats$stats_sh_corn$x_sd_2019,
                  stats$stats_k_land$x_mean_2019, stats$stats_k_land$x_sd_2019,
                  stats$stats_k_labor$x_mean_2019, stats$stats_k_labor$x_sd_2019,
                  stats$stats_k_capital$x_mean_2019, stats$stats_k_capital$x_sd_2019,
                  stats$stats_gdd_obs$x_mean_2019, stats$stats_gdd_obs$x_sd_2019,
                  stats$stats_prec_obs$x_mean_2019, stats$stats_prec_obs$x_sd_2019,
                  stats$stats_gddHigh_obs$x_mean_2019, stats$stats_gddHigh_obs$x_sd_2019,
                  stats$stats_dd_obs$x_mean_2019, stats$stats_dd_obs$x_sd_2019,
                  stats$stats_gdd_1to3$x_mean_2019, stats$stats_gdd_1to3$x_sd_2019,
                  stats$stats_prec_1to3$x_mean_2019, stats$stats_prec_1to3$x_sd_2019,
                  stats$stats_gddHigh_1to3$x_mean_2019, stats$stats_gddHigh_1to3$x_sd_2019,
                  stats$stats_dd_1to3$x_mean_2019, stats$stats_dd_1to3$x_sd_2019,
                  stats$stats_gdd_4to10$x_mean_2019, stats$stats_gdd_4to10$x_sd_2019,
                  stats$stats_prec_4to10$x_mean_2019, stats$stats_prec_4to10$x_sd_2019,
                  stats$stats_gddHigh_4to10$x_mean_2019, stats$stats_gddHigh_4to10$x_sd_2019,
                  stats$stats_dd_4to10$x_mean_2019, stats$stats_dd_4to10$x_sd_2019)
)

#Numbers smaller than 0.0005 should be in scientific format; others in 0.001 etc.
Tab_stats_app <- as.data.frame(cbind(Tab_stats_app[,1], sapply(Tab_stats_app[,2:5], function(x) ifelse(abs(x)<0.0005,format(x,scientific=TRUE,digits=2),format(round(x, 3), nsmall = 3, scientific=FALSE)) )) )

# Put parentheses around standard deviations
seq(from = 2, to = nrow(Tab_stats_app), by = 2)

for (i in seq(from = 2, to = nrow(Tab_stats_app), by = 2) ) {
  Tab_stats_app[i,2:5] <- paste0("(",Tab_stats_app[i,2:5],")") 
}

# remove white space
Tab_stats_app$Average <- gsub(" ", "", Tab_stats_app$Average, fixed = TRUE)
Tab_stats_app$Year.1996 <- gsub(" ", "", Tab_stats_app$Year.1996, fixed = TRUE)
Tab_stats_app$Year.2008 <- gsub(" ", "", Tab_stats_app$Year.2008, fixed = TRUE)
Tab_stats_app$Year.2019 <- gsub(" ", "", Tab_stats_app$Year.2019, fixed = TRUE)

# Write to excel
write_xlsx(Tab_stats_app,"Tables/Table_S3.1.xlsx")

#---------------------------------------------------------------#
#### 9) Calculate interaction terms to assess heterogeneity ####
#---------------------------------------------------------------#

df_farm$gdd_obs_land <- df_farm$gdd_obs * df_farm$k_land
df_farm$prec_obs_land <- df_farm$prec_obs * df_farm$k_land
df_farm$gddHigh_obs_land <- df_farm$gddHigh_obs * df_farm$k_land
df_farm$dd_obs_land <- df_farm$dd_obs * df_farm$k_land

df_farm$gdd_1to3_land <- df_farm$gdd_1to3 * df_farm$k_land
df_farm$prec_1to3_land <- df_farm$prec_1to3 * df_farm$k_land
df_farm$gddHigh_1to3_land <- df_farm$gddHigh_1to3 * df_farm$k_land
df_farm$dd_1to3_land <- df_farm$dd_1to3 * df_farm$k_land

df_farm$gdd_4to10_land <- df_farm$gdd_4to10 * df_farm$k_land
df_farm$prec_4to10_land <- df_farm$prec_4to10 * df_farm$k_land
df_farm$gddHigh_4to10_land <- df_farm$gddHigh_4to10 * df_farm$k_land
df_farm$dd_4to10_land <- df_farm$dd_4to10 * df_farm$k_land

df_farm$gdd_1to5_land <- df_farm$gdd_1to5 * df_farm$k_land
df_farm$prec_1to5_land <- df_farm$prec_1to5 * df_farm$k_land
df_farm$gddHigh_1to5_land <- df_farm$gddHigh_1to5 * df_farm$k_land
df_farm$dd_1to5_land <- df_farm$dd_1to5 * df_farm$k_land

df_farm$gdd_6to10_land <- df_farm$gdd_6to10 * df_farm$k_land
df_farm$prec_6to10_land <- df_farm$prec_6to10 * df_farm$k_land
df_farm$gddHigh_6to10_land <- df_farm$gddHigh_6to10 * df_farm$k_land
df_farm$dd_6to10_land <- df_farm$dd_6to10 * df_farm$k_land

df_farm$gdd_6to20_land <- df_farm$gdd_6to20 * df_farm$k_land
df_farm$prec_6to20_land <- df_farm$prec_6to20 * df_farm$k_land
df_farm$gddHigh_6to20_land <- df_farm$gddHigh_6to20 * df_farm$k_land
df_farm$dd_6to20_land <- df_farm$dd_6to20 * df_farm$k_land

df_farm$gdd_1_land <- df_farm$gdd_1 * df_farm$k_land
df_farm$prec_1_land <- df_farm$prec_1 * df_farm$k_land
df_farm$gddHigh_1_land <- df_farm$gddHigh_1 * df_farm$k_land
df_farm$dd_1_land <- df_farm$dd_1 * df_farm$k_land

df_farm$gdd_2to5_land <- df_farm$gdd_2to5 * df_farm$k_land
df_farm$prec_2to5_land <- df_farm$prec_2to5 * df_farm$k_land
df_farm$gddHigh_2to5_land <- df_farm$gddHigh_2to5 * df_farm$k_land
df_farm$dd_2to5_land <- df_farm$dd_2to5 * df_farm$k_land

df_farm$gdd_2to10_land <- df_farm$gdd_2to10 * df_farm$k_land
df_farm$prec_2to10_land <- df_farm$prec_2to10 * df_farm$k_land
df_farm$gddHigh_2to10_land <- df_farm$gddHigh_2to10 * df_farm$k_land
df_farm$dd_2to10_land <- df_farm$dd_2to10 * df_farm$k_land

df_farm$gdd_1to10_land <- df_farm$gdd_1to10 * df_farm$k_land
df_farm$prec_1to10_land <- df_farm$prec_1to10 * df_farm$k_land
df_farm$gddHigh_1to10_land <- df_farm$gddHigh_1to10 * df_farm$k_land
df_farm$dd_1to10_land <- df_farm$dd_1to10 * df_farm$k_land

df_farm$gdd_1to30_land <- df_farm$gdd_1to30 * df_farm$k_land
df_farm$prec_1to30_land <- df_farm$prec_1to30 * df_farm$k_land
df_farm$gddHigh_1to30_land <- df_farm$gddHigh_1to30 * df_farm$k_land
df_farm$dd_1to30_land <- df_farm$dd_1to30 * df_farm$k_land

#---------------------------------------------------#
#### 10) Calculate farm means for Mundlak device ####
#---------------------------------------------------#

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

df_farm <- df_farm %>%
  group_by(key) %>%
  summarise_at(c(all_of(allindepvars)), list(fm = mean)) %>% 
  left_join(df_farm)

#-----------------------------------------------------#
#                                                     #
# End of sample construction and variable preparation # 
#                                                     #
#-----------------------------------------------------#

save(df_farm, file="rOutput/farm_ready.Rda")


