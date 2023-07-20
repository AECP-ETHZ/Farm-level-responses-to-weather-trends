#------------------------------------------------------------------#
# This program prepares the farm data:                             #
# 0) Some preliminaries                                            #
# 1) Merge weather and price data with farm data                   #
# 2) Rename area variables                                         #
# 3) Farm-level quantities, prices, and revenues                   #
# 4) Regional prices                                               #
# 5) Merge regional output price data to farm-level data           #
# 6) Quantities of fixed inputs                                    #
# 7) Prices of variable inputs                                     #
# 8) Merge regional price data to farm-level data                  #
#                                                                  #
# The outcome is saved in "rOutput/farm_outp_inp.Rda" and will be  #
# loaded in the next program.                                      #
#                                                                  #
#------------------------------------------------------------------#

#-----------------------------#
#### 0) Some preliminaries ####
#-----------------------------#

# open libraries
library(dplyr)
library(readxl) 
library(writexl) 

# load data
df_farm <- read_excel("Data/tbn_FINAL_all.xlsx", col_names = TRUE)
df_farm_raw <- df_farm

# define nuts levels
df_farm$nuts1 <- substr(df_farm$AGS_total, 1, 2)
df_farm$nuts2 <- substr(df_farm$AGS_total, 1, 3)
df_farm$nuts3 <- substr(df_farm$AGS_total, 1, 5)
df_farm$nuts5 <- df_farm$AGS_total 

# rename the year variable
df_farm <- dplyr::rename(df_farm, c("year"="jahr"))

# sort data by key and year
df_farm <- df_farm %>% 
  arrange(key, year)
  
#------------------------------------------------------#
#### 1) Merge weather and price data with farm data ####
#------------------------------------------------------#

#--------------------------#
##### 1a) Weather data #####
#--------------------------#
  
  # read weather data
  load("rOutput/df_weather.Rda")
  
  # make sure merging vars have same class
  class(df_weather$nuts5)
  class(df_farm$nuts5)
  df_weather$nuts5 <- as.character(df_weather$nuts5)
  
  # merge, then remove weather data frame
  df_farm <- left_join(df_farm, df_weather)
  remove(df_weather)
  
  # Drop obs which cannot be matched (<1%)
  df_farm <- df_farm %>% 
    filter(!is.na(gdd_832))

#------------------------#
##### 1b) Price data #####
#------------------------#
  
  # read price data
  load("rOutput/wi_nuts0.Rda")
  
  # merge
  df_farm <- left_join(df_farm, wi_nuts0, by="year") 

#-----------------------------#
#### 2) Rename area variables #
#-----------------------------#

  # Cereals (without corn)
  df_farm <- dplyr::rename(df_farm, 
                           c("aa_wwheat"="z4001s02",
                             "aa_swheat"="z4002s02",
                             "aa_rye"="z4004s02",
                             "aa_wbarley"="z4005s02",
                             "aa_sbarley"="z4006s02",
                             "aa_oat"="z4007s02",
                             "aa_scer"="z4008s02",
                             "aa_trit"="z4012s02",
                             "aa_othcer"="z4017s02",
                             "aa_encer" = "z4031s02"))
  
  # Protein crops 
  df_farm <- dplyr::rename(df_farm, 
                           c("aa_beans"="z4020s02",
                             "aa_peas"="z4021s02")) 
  
  # Oilseed crops
  df_farm <- dplyr::rename(df_farm, 
                           c("aa_wcanola"="z4024s02",
                             "aa_scanola"="z4025s02",
                             "aa_othoil"="z4028s02",
                             "aa_fibre"="z4029s02",
                             "aa_flax"="z4030s02",
                             "aa_enoil" = "z4034s02"))
  
  # Root crops
  df_farm <- dplyr::rename(df_farm, 
                           c("aa_potat"="z4039s02",
                             "aa_beets"="z4040s02",
                             "aa_enbeets"="z4035s02")) 
  
  # Corn
  df_farm <- dplyr::rename(df_farm, 
                           c("aa_corn"="z4010s02",
                             "aa_silmaize"="z4070s02", 
                             "aa_enmaize" = "z4032s02")) 
  
  
  # replace missing values by zeroes
  aa_allcrops <- c("aa_wwheat", "aa_swheat","aa_rye","aa_wbarley","aa_sbarley","aa_oat","aa_scer","aa_trit","aa_othcer","aa_encer",
                   "aa_beans","aa_peas",
                   "aa_wcanola","aa_scanola", "aa_othoil", "aa_fibre","aa_flax","aa_enoil",
                   "aa_potat","aa_beets","aa_enbeets",
                   "aa_corn", "aa_silmaize","aa_enmaize")
  df_farm <- mutate_at(df_farm, c(aa_allcrops), ~replace(., is.na(.), 0))
  

#------------------------------------------------------#
#### 3) Farm-level quantities, prices, and revenues ####
#------------------------------------------------------#
  
#----------------------------------------------#
##### 3a) Calculate farm-level quantities  #####
#----------------------------------------------#
    
  # Cereals
  df_farm$qq_wwheat <- df_farm$aa_wwheat * df_farm$z4001s03
  df_farm$qq_swheat <- df_farm$aa_swheat * df_farm$z4002s03
  df_farm$qq_rye <- df_farm$aa_rye * df_farm$z4004s03
  df_farm$qq_wbarley <- df_farm$aa_wbarley * df_farm$z4005s03
  df_farm$qq_sbarley <- df_farm$aa_sbarley * df_farm$z4006s03
  df_farm$qq_oat <- df_farm$aa_oat * df_farm$z4007s03
  df_farm$qq_scer <- df_farm$aa_scer * df_farm$z4008s03
  df_farm$qq_trit <- df_farm$aa_trit * df_farm$z4012s03
  df_farm$qq_othcer <- df_farm$aa_othcer * df_farm$z4017s03
  df_farm$qq_encer <- df_farm$aa_encer * df_farm$z4031s03
  
  # Protein crops 
  df_farm$qq_beans <- df_farm$aa_beans * df_farm$z4020s03
  df_farm$qq_peas <- df_farm$aa_peas * df_farm$z4021s03
  
  # Oilseed crops
  df_farm$qq_wcanola <- df_farm$aa_wcanola * df_farm$z4024s03
  df_farm$qq_scanola <- df_farm$aa_scanola * df_farm$z4025s03
  df_farm$qq_othoil <- df_farm$aa_othoil * df_farm$z4028s03
  df_farm$qq_fibre <- df_farm$aa_fibre * df_farm$z4029s03
  df_farm$qq_flax <- df_farm$aa_flax * df_farm$z4030s03
  df_farm$qq_enoil <- df_farm$aa_enoil * df_farm$z4034s03
  
  # Root crops
  df_farm$qq_potat <- df_farm$aa_potat * df_farm$z4039s03
  df_farm$qq_beets <- df_farm$aa_beets * df_farm$z4040s03
  df_farm$qq_enbeets <- df_farm$aa_enbeets * df_farm$z4035s03
  
  # Corn
  df_farm$qq_corn <- df_farm$aa_corn * df_farm$z4010s03
  df_farm$qq_silmaize <- df_farm$z5670s04 # reliable 
  df_farm$qq_enmaize <- df_farm$z5632s04 # reliable

  # Calculate area and quantity for total maize (energy and silage)
    
    # Area
    df_farm$aa_silmaize[is.na(df_farm$aa_silmaize)] <- 0 # replace missing values with zero
    df_farm$aa_enmaize[is.na(df_farm$aa_enmaize)] <- 0 # replace missing values with zero
    df_farm$aa_maize = df_farm$aa_silmaize + df_farm$aa_enmaize
    
    # Quantity
    df_farm$qq_silmaize[is.na(df_farm$qq_silmaize)] <- 0 # replace missing values with zero
    df_farm$qq_enmaize[is.na(df_farm$qq_enmaize)] <- 0 # replace missing values with zero
    df_farm$qq_maize = df_farm$qq_silmaize + df_farm$qq_enmaize
  
#-----------------------------------------#
##### 3b) Calculate farm-level prices #####
#-----------------------------------------#
  
  # Cereals
  df_farm <- dplyr::rename(df_farm, 
                           c("p_wwheat"="z4001s04",
                             "p_swheat"="z4002s04",
                             "p_rye"="z4004s04",
                             "p_wbarley"="z4005s04",
                             "p_sbarley"="z4006s04",
                             "p_oat"="z4007s04",
                             "p_scer"="z4008s04",
                             "p_trit"="z4012s04",
                             "p_othcer"="z4017s04",
                             "p_encer" = "z4031s04"))
  
  # Protein crops 
  df_farm <- dplyr::rename(df_farm, 
                           c("p_beans"="z4020s04",
                             "p_peas"="z4021s04"))
  
  # Oilseed crops
  df_farm <- dplyr::rename(df_farm, 
                           c("p_wcanola"="z4024s04",
                             "p_scanola"="z4025s04",
                             "p_othoil"="z4028s04",
                             "p_fibre"="z4029s04",
                             "p_flax"="z4030s04",
                             "p_enoil" = "z4034s04"))
  
  # Root crops
  df_farm <- dplyr::rename(df_farm, 
                           c("p_potat"="z4039s04",
                             "p_beets"="z4040s04",
                             "p_enbeets"="z4035s04"))
  
  # Corn
  df_farm <- dplyr::rename(df_farm, 
                           "p_corn"="z4010s04")
  df_farm$p_maize <- (df_farm$z2032s05 + df_farm$z2070s05) / 
                        (df_farm$z5632s05 + df_farm$z5670s05) # Note: Prices not directly reported, but can be calculated from earnings and sold quantities)
  
  # In case earnings are indicated but no sold quantities, p must be set to missing
  df_farm$p_maize[is.infinite(df_farm$p_maize)] <- NA 


#-------------------------------------------#
##### 3c) Calculate farm-level revenues #####
#-------------------------------------------#
  
  # Cereals
  df_farm$rev_wwheat <- df_farm$qq_wwheat * df_farm$p_wwheat
  df_farm$rev_swheat <- df_farm$qq_swheat * df_farm$p_swheat
  df_farm$rev_rye <- df_farm$qq_rye * df_farm$p_rye
  df_farm$rev_wbarley <- df_farm$qq_wbarley * df_farm$p_wbarley
  df_farm$rev_sbarley <- df_farm$qq_sbarley * df_farm$p_sbarley
  df_farm$rev_oat <- df_farm$qq_oat * df_farm$p_oat
  df_farm$rev_scer <- df_farm$qq_scer * df_farm$p_scer
  df_farm$rev_trit <- df_farm$qq_trit * df_farm$p_trit
  df_farm$rev_othcer <- df_farm$qq_othcer * df_farm$p_othcer
  df_farm$rev_encer <- df_farm$qq_encer * df_farm$p_encer
  
  # Protein crops 
  df_farm$rev_beans <- df_farm$qq_beans * df_farm$p_beans
  df_farm$rev_peas <- df_farm$qq_peas * df_farm$p_peas

  # Oilseed crops
  df_farm$rev_wcanola <- df_farm$qq_wcanola * df_farm$p_wcanola
  df_farm$rev_scanola <- df_farm$qq_scanola * df_farm$p_scanola
  df_farm$rev_othoil <- df_farm$qq_othoil * df_farm$p_othoil
  df_farm$rev_fibre <- df_farm$qq_fibre * df_farm$p_fibre
  df_farm$rev_flax <- df_farm$qq_flax * df_farm$p_flax
  df_farm$rev_enoil <- df_farm$qq_enoil * df_farm$p_enoil
  
  # Root crops
  df_farm$rev_potat <- df_farm$qq_potat * df_farm$p_potat
  df_farm$rev_beets <- df_farm$qq_beets * df_farm$p_beets
  df_farm$rev_enbeets <- df_farm$qq_enbeets * df_farm$p_enbeets
  
  # Corn
  df_farm$rev_corn <- df_farm$qq_corn * df_farm$p_corn
  df_farm$rev_maize <- df_farm$qq_maize * df_farm$p_maize

#--------------------------#
#### 4) Regional prices ####
#--------------------------#
  
  # Note: We calculate global averages. This means we take the sum of crop-revenues 
  #       in a specific region and divide this by the sum of crop-quantities in the 
  #       same region. These are used later to compute the prices for the crop
  #       categories; We use global average prices instead of farm-level prices
  #       because farm-level prices are missing if a certain farm does not produce 
  #       this crop. 

  # Merge individual crops into crop categories. 
  cernames <- c("wwheat","swheat","rye","wbarley","sbarley","oat",
                "scer","trit","othcer","encer")
  protnames <- c("beans", "peas") 
  oilnames <- c("wcanola", "scanola", "othoil", "fibre", "flax", "enoil")
  rootnames <- c("potat", "beets", "enbeets") 
  cornnames <- c("corn", "maize") #corn is grain; maize is for energy and silage 
  
#------------------------------------------------------------#
##### 4a) Calculate regional prices for individual crops #####
#------------------------------------------------------------#
  
  # define variable names
  cropnames <- c("wwheat","swheat","rye","wbarley","sbarley","oat", "scer","trit","othcer","encer",
                 "corn", 
                 "beans", "peas",
                 "wcanola","scanola", "othoil","fibre","flax","enoil",
                 "potat","beets","enbeets",
                 "silmaize", "enmaize", "maize","clover","silcer","otherfodd", "fodroots")
  
  qq_cropnames <- paste0("qq_", cropnames, sep="") #qq means quantity
  aa_cropnames <- paste0("aa_", cropnames, sep="") #aa means area
  rev_cropnames <- paste0("rev_", cropnames, sep="") #rev means revenue
  p_cropnames <- paste0("p_", cropnames, sep="") #p means prices
  
  # replace NAs for all qq, aa, and rev in main data set with zeros
  df_farm[qq_cropnames][is.na(df_farm[qq_cropnames])] <- 0
  df_farm[aa_cropnames][is.na(df_farm[aa_cropnames])] <- 0
  df_farm[rev_cropnames][is.na(df_farm[rev_cropnames])] <- 0
  
  #---------------------------------#
  # nuts2-prices for all crop items #
  #---------------------------------#
  
  p_allcrops_nuts2 <- df_farm %>%
    group_by(nuts2,year) %>% # group by nuts2 and year to obtain the global averages
    summarise_at(c(all_of(rev_cropnames),all_of(qq_cropnames)), 
                 list(sum), na.rm = TRUE) %>% 
    mutate(p_wwheat = rev_wwheat / qq_wwheat,
           p_swheat = rev_swheat / qq_swheat,
           p_rye = rev_rye / qq_rye,
           p_wbarley = rev_wbarley / qq_wbarley,
           p_sbarley = rev_sbarley / qq_sbarley,
           p_oat = rev_oat / qq_oat,
           p_scer = rev_scer / qq_scer,
           p_trit = rev_trit / qq_trit,
           p_othcer = rev_othcer / qq_othcer,
           p_encer = rev_encer / qq_encer,
           p_beans = rev_beans / qq_beans, 
           p_peas = rev_peas / qq_peas, 
           p_wcanola = rev_wcanola / qq_wcanola,
           p_scanola = rev_scanola / qq_scanola, 
           p_othoil = rev_othoil / qq_othoil, 
           p_fibre = rev_fibre / qq_fibre, 
           p_flax = rev_flax / qq_flax, 
           p_enoil = rev_enoil / qq_enoil,
           p_potat = rev_potat / qq_potat, 
           p_beets = rev_beets / qq_beets,
           p_enbeets = rev_enbeets / qq_enbeets,
           p_corn = rev_corn / qq_corn,
           p_maize = rev_maize / qq_maize) %>% 
    select(nuts2,year,all_of(p_cropnames))
  
  # In rare cases, calculation at the farm-level yields NaN, Inf, or 0. The reasons are:
    # 1) Calculating 0/0 --> NaN (fine to be missing as no price available --> change to NA)
    # 2) Calculating 1/0 --> Inf (should be missing --> change to NA) 
    # 3) Calculating 0/1 --> 0   (should be missing --> change to NA)
  
  # Replace NaN, Inf by NA (cases 1+2)
  p_allcrops_nuts2 <- do.call(data.frame,                      
                              lapply(p_allcrops_nuts2,
                                     function(x) replace(x, is.infinite(x) | 
                                                           is.nan(x), NA))) 
  # Replace 0 by NA (case 3)
  p_allcrops_nuts2[p_allcrops_nuts2 == 0] <- NA
  
  #---------------------------------#
  # nuts1-prices for all crop items #
  #---------------------------------#
  
  p_allcrops_nuts1 <- df_farm %>%
    group_by(nuts1,year) %>% # group by nuts1 and year to obtain the global averages
    summarise_at(c(all_of(rev_cropnames),all_of(qq_cropnames)), 
                 list(sum), na.rm = TRUE) %>% 
    mutate(p_wwheat = rev_wwheat / qq_wwheat,
           p_swheat = rev_swheat / qq_swheat,
           p_rye = rev_rye / qq_rye,
           p_wbarley = rev_wbarley / qq_wbarley,
           p_sbarley = rev_sbarley / qq_sbarley,
           p_oat = rev_oat / qq_oat,
           p_scer = rev_scer / qq_scer,
           p_trit = rev_trit / qq_trit,
           p_othcer = rev_othcer / qq_othcer,
           p_encer = rev_encer / qq_encer,
           p_beans = rev_beans / qq_beans, 
           p_peas = rev_peas / qq_peas, 
           p_wcanola = rev_wcanola / qq_wcanola,
           p_scanola = rev_scanola / qq_scanola, 
           p_othoil = rev_othoil / qq_othoil, 
           p_fibre = rev_fibre / qq_fibre, 
           p_flax = rev_flax / qq_flax, 
           p_enoil = rev_enoil / qq_enoil,
           p_potat = rev_potat / qq_potat, 
           p_beets = rev_beets / qq_beets,
           p_enbeets = rev_enbeets / qq_enbeets,
           p_corn = rev_corn / qq_corn,
           p_maize = rev_maize / qq_maize) %>% 
    select(nuts1,year,all_of(p_cropnames))
  
  # In rare cases, calculation at the farm-level yields NaN, Inf, or 0. The reasons are:
    # 1) Calculating 0/0 --> NaN (fine to be missing as no price available --> change to NA)
    # 2) Calculating 1/0 --> Inf (should be missing --> change to NA) 
    # 3) Calculating 0/1 --> 0   (should be missing --> change to NA)
    
  # Replace NaN, Inf by NA (cases 1+2)
  p_allcrops_nuts1 <- do.call(data.frame,                      
                              lapply(p_allcrops_nuts1,
                                     function(x) replace(x, is.infinite(x) | 
                                                           is.nan(x), NA))) 
  # Replace 0 by NA (case 3)
  p_allcrops_nuts1[p_allcrops_nuts1 == 0] <- NA
  
  #---------------------------------#
  # nuts0-prices for all crop items #
  #---------------------------------#
  
  p_allcrops_nuts0 <- df_farm %>%
    group_by(year) %>% # group by year to obtain the global averages
    summarise_at(c(all_of(rev_cropnames),all_of(qq_cropnames)), 
                 list(sum), na.rm = TRUE) %>% 
    mutate(p_wwheat = rev_wwheat / qq_wwheat,
           p_swheat = rev_swheat / qq_swheat,
           p_rye = rev_rye / qq_rye,
           p_wbarley = rev_wbarley / qq_wbarley,
           p_sbarley = rev_sbarley / qq_sbarley,
           p_oat = rev_oat / qq_oat,
           p_scer = rev_scer / qq_scer,
           p_trit = rev_trit / qq_trit,
           p_othcer = rev_othcer / qq_othcer,
           p_encer = rev_encer / qq_encer,
           p_beans = rev_beans / qq_beans, 
           p_peas = rev_peas / qq_peas, 
           p_wcanola = rev_wcanola / qq_wcanola,
           p_scanola = rev_scanola / qq_scanola, 
           p_othoil = rev_othoil / qq_othoil, 
           p_fibre = rev_fibre / qq_fibre, 
           p_flax = rev_flax / qq_flax, 
           p_enoil = rev_enoil / qq_enoil,
           p_potat = rev_potat / qq_potat, 
           p_beets = rev_beets / qq_beets,
           p_enbeets = rev_enbeets / qq_enbeets,
           p_corn = rev_corn / qq_corn,
           p_maize = rev_maize / qq_maize) %>% 
    select(year,all_of(p_cropnames))
  
  
  # In rare cases, calculation at the farm-level yields NaN, Inf, or 0. The reasons are:
    # 1) Calculating 0/0 --> NaN (fine to be missing as no price available --> change to NA)
    # 2) Calculating 1/0 --> Inf (should be missing --> change to NA) 
    # 3) Calculating 0/1 --> 0   (should be missing --> change to NA)
  
  # Replace NaN, Inf by NA (cases 1+2)
  p_allcrops_nuts0 <- do.call(data.frame,                      
                              lapply(p_allcrops_nuts0,
                                     function(x) replace(x, is.infinite(x) | 
                                                           is.nan(x), NA))) 
  # Replace 0 by NA (case 3)
  p_allcrops_nuts0[p_allcrops_nuts0 == 0] <- NA
  

#---------------------------------------------------------------------------#
##### 4b) Calculate farm-level areas and quantities for crop categories #####
#---------------------------------------------------------------------------#
  
  #-------------------------#
  ###### 4b.1) Cereals ######
  #-------------------------#
    
    # Create vectors with cereal names
    qq_cernames <- paste0("qq_", cernames, sep="")
    aa_cernames <- paste0("aa_", cernames, sep="")
    rev_cernames <- paste0("rev_", cernames, sep="")
    p_cernames <- paste0("p_", cernames, sep="")
    
    df_farm <- df_farm %>%
      rowwise() %>% 
      mutate(aa_cereals = sum(c_across(all_of(aa_cernames)))) 
    
    df_farm <- df_farm %>%
      rowwise() %>% 
      mutate(qq_cereals = sum(c_across(all_of(qq_cernames)))) 
    
  #-------------------------#
  ###### 4b.2) Protein ######
  #-------------------------#    
          
    # Create vectors with protein names
    qq_protnames <- paste0("qq_", protnames, sep="")
    aa_protnames <- paste0("aa_", protnames, sep="")
    rev_protnames <- paste0("rev_", protnames, sep="")
    p_protnames <- paste0("p_", protnames, sep="")
    
    df_farm <- df_farm %>%
      rowwise() %>% 
      mutate(aa_protein = sum(c_across(all_of(aa_protnames)))) 
    
    df_farm <- df_farm %>%
      rowwise() %>% 
      mutate(qq_protein = sum(c_across(all_of(qq_protnames)))) 
    
  #-------------------------#
  ###### 4b.3) Oilseed ######
  #-------------------------#  
          
    qq_oilnames <- paste0("qq_", oilnames, sep="")
    aa_oilnames <- paste0("aa_", oilnames, sep="")
    rev_oilnames <- paste0("rev_", oilnames, sep="")
    p_oilnames <- paste0("p_", oilnames, sep="")
    
    df_farm <- df_farm %>%
      rowwise() %>% 
      mutate(aa_oilseed = sum(c_across(all_of(aa_oilnames)))) 
    
    df_farm <- df_farm %>%
      rowwise() %>% 
      mutate(qq_oilseed = sum(c_across(all_of(qq_oilnames)))) 
    
  #----------------------------#
  ###### 4b.4) Root crops ######
  #----------------------------#  
    
    # Create vectors with root crops names
    qq_rootnames <- paste0("qq_", rootnames, sep="")
    aa_rootnames <- paste0("aa_", rootnames, sep="")
    rev_rootnames <- paste0("rev_", rootnames, sep="")
    p_rootnames <- paste0("p_", rootnames, sep="")
    
    df_farm <- df_farm %>%
      rowwise() %>% 
      mutate(aa_roots = sum(c_across(all_of(aa_rootnames)))) 
    
    df_farm <- df_farm %>%
      rowwise() %>% 
      mutate(qq_roots = sum(c_across(all_of(qq_rootnames)))) 
    
  #-----------------------------------------#
  ###### 4b.5) Corn (grain and silage) ######
  #-----------------------------------------# 
          
    # Create vectors with corn crops names
    qq_cornnames <- paste0("qq_", cornnames, sep="")
    aa_cornnames <- paste0("aa_", cornnames, sep="")
    rev_cornnames <- paste0("rev_", cornnames, sep="")
    p_cornnames <- paste0("p_", cornnames, sep="")
    
    df_farm <- df_farm %>%
      rowwise() %>% 
      mutate(aa_corn = sum(c_across(all_of(aa_cornnames)))) 
    
    df_farm <- df_farm %>%
      rowwise() %>% 
      mutate(qq_corn = sum(c_across(all_of(qq_cornnames)))) 
          
#-----------------------------------------------------------#
##### 4c) Calculate regional prices for crop categories #####
#-----------------------------------------------------------#

#-------------------------#
###### 4c.1) Cereals ######
#-------------------------#
          
  #-----------------------------------------#
  # Cereals (w/o corn) price at nuts2-level #
  #-----------------------------------------#

    # nuts2-land shares:
    
      #1) sum of individual cereal items, and total cereal area
      aa_cereals_nuts2 <- df_farm %>%
        group_by(nuts2,year) %>%
        summarise_at(all_of(aa_cernames), list(sum), na.rm = TRUE) %>% 
        rowwise() %>% 
        mutate(SumCerealArea = sum(c_across(all_of(aa_cernames)))) 
      
      #2) Share of individual cereal items over total cereal area
      aa_cereals_nuts2[ ,aa_cernames] <- lapply(aa_cereals_nuts2[ ,aa_cernames],
                                                function(x) 
                                                  x / aa_cereals_nuts2$SumCerealArea)  
    
    # nuts2-surface-weighted price for cereals
    p_cereals_nuts2 <- p_allcrops_nuts2[,c("year", "nuts2", p_cernames)] #subset cer
    p_cereals_nuts2 <- left_join(p_cereals_nuts2, 
                                 aa_cereals_nuts2) # merge with area shares
    p_cereals_nuts2[is.na(p_cereals_nuts2)] <- 0 # replace missings
    attach(p_cereals_nuts2)
    p_cereals_nuts2$n2p_cereals <- 
      p_wwheat * aa_wwheat + # note: aa_ are shares
      p_swheat * aa_swheat + 
      p_rye * aa_rye + 
      p_wbarley * aa_wbarley + 
      p_sbarley * aa_sbarley + 
      p_oat * aa_oat + 
      p_scer * aa_scer + 
      p_trit * aa_trit + 
      p_othcer * aa_othcer +
      p_encer * aa_encer
    p_cereals_nuts2 <- select(p_cereals_nuts2,nuts2,year,n2p_cereals)
    detach(p_cereals_nuts2)

  #-----------------------------------------#
  # Cereals (w/o corn) price at nuts1-level #
  #-----------------------------------------#

    # nuts1-land shares:
    
      #1) sum of individual cereal items, and total cereal area
      aa_cereals_nuts1 <- df_farm %>%
        group_by(nuts1,year) %>%
        summarise_at(all_of(aa_cernames), list(sum), na.rm = TRUE) %>% 
        rowwise() %>% 
        mutate(SumCerealArea = sum(c_across(all_of(aa_cernames)))) 
      
      #2) Share of individual cereal items over total cereal area
      aa_cereals_nuts1[ ,aa_cernames] <- lapply(aa_cereals_nuts1[ ,aa_cernames],
                                                function(x) 
                                                  x / aa_cereals_nuts1$SumCerealArea)  
      
    # nuts1-surface-weighted price for cereals
    p_cereals_nuts1 <- p_allcrops_nuts1[,c("year", "nuts1", p_cernames)] #subset cer
    p_cereals_nuts1 <- left_join(p_cereals_nuts1, 
                                 aa_cereals_nuts1) # merge with area shares
    p_cereals_nuts1[is.na(p_cereals_nuts1)] <- 0 # replace missings
    attach(p_cereals_nuts1)
    p_cereals_nuts1$n1p_cereals <- 
      p_wwheat * aa_wwheat + # note: aa_ are shares
      p_swheat * aa_swheat + 
      p_rye * aa_rye + 
      p_wbarley * aa_wbarley + 
      p_sbarley * aa_sbarley + 
      p_oat * aa_oat + 
      p_scer * aa_scer + 
      p_trit * aa_trit + 
      p_othcer * aa_othcer +
      p_encer * aa_encer
    p_cereals_nuts1 <- select(p_cereals_nuts1,nuts1,year,n1p_cereals)
    detach(p_cereals_nuts1)

  #-----------------------------------------#
  # Cereals (w/o corn) price at nuts0-level #
  #-----------------------------------------#

    # nuts0-land shares:
    
      #1) sum of individual cereal items, and total cereal area
      aa_cereals_nuts0 <- df_farm %>%
        group_by(year) %>%
        summarise_at(all_of(aa_cernames), list(sum), na.rm = TRUE) %>% 
        rowwise() %>% 
        mutate(SumCerealArea = sum(c_across(all_of(aa_cernames)))) 
      
      #2) Share of individual cereal items over total cereal area
      aa_cereals_nuts0[ ,aa_cernames] <- lapply(aa_cereals_nuts0[ ,aa_cernames],
                                                function(x) 
                                                  x / aa_cereals_nuts0$SumCerealArea)  
      
      # nuts0-surface-weighted price for cereals
      p_cereals_nuts0 <- p_allcrops_nuts0[,c("year", p_cernames)] #subset cer
      p_cereals_nuts0 <- left_join(p_cereals_nuts0, 
                                   aa_cereals_nuts0) # merge with area shares
      p_cereals_nuts0[is.na(p_cereals_nuts0)] <- 0 # replace missings
      attach(p_cereals_nuts0)
      p_cereals_nuts0$n0p_cereals <- 
        p_wwheat * aa_wwheat + # note: aa_ are shares
        p_swheat * aa_swheat + 
        p_rye * aa_rye + 
        p_wbarley * aa_wbarley + 
        p_sbarley * aa_sbarley + 
        p_oat * aa_oat + 
        p_scer * aa_scer + 
        p_trit * aa_trit + 
        p_othcer * aa_othcer +
        p_encer * aa_encer
      p_cereals_nuts0 <- select(p_cereals_nuts0,year,n0p_cereals)
      detach(p_cereals_nuts0)

      
#-------------------------#
###### 4c.2) Protein ######
#-------------------------#
      
  #------------------------------#
  # Protein price at nuts2-level #
  #------------------------------#

    # nuts2-land shares:
    
      #1) sum of individual protein items, and total protein area
      aa_protein_nuts2 <- df_farm %>%
        group_by(nuts2,year) %>%
        summarise_at(all_of(aa_protnames), list(sum), na.rm = TRUE) %>% 
        rowwise() %>% 
        mutate(SumProteinArea = sum(c_across(all_of(aa_protnames)))) 
      
      #2) Share of individual protein items over total protein area
      aa_protein_nuts2[ ,aa_protnames] <- lapply(aa_protein_nuts2[ ,aa_protnames],
                                                 function(x) 
                                                   x / aa_protein_nuts2$SumProteinArea)  
    
    # nuts2-surface-weighted price for protein crops
    p_protein_nuts2 <- p_allcrops_nuts2[,c("year", "nuts2", p_protnames)] #subset prot
    p_protein_nuts2 <- left_join(p_protein_nuts2, 
                                 aa_protein_nuts2) # merge with area shares
    p_protein_nuts2[is.na(p_protein_nuts2)] <- 0 # replace missings
    attach(p_protein_nuts2)
    p_protein_nuts2$n2p_protein <- 
      p_beans * aa_beans + # note: aa_ are shares
      p_peas * aa_peas
    p_protein_nuts2 <- select(p_protein_nuts2, c("nuts2","year",n2p_protein))
    detach(p_protein_nuts2)

  
  #------------------------------#
  # Protein price at nuts1-level #
  #------------------------------#

    # nuts1-land shares:
    
      #1) sum of individual protein items, and total protein area
      aa_protein_nuts1 <- df_farm %>%
        group_by(nuts1,year) %>%
        summarise_at(all_of(aa_protnames), list(sum), na.rm = TRUE) %>% 
        rowwise() %>% 
        mutate(SumProteinArea = sum(c_across(all_of(aa_protnames)))) 
      
      #2) Share of individual protein items over total protein area
      aa_protein_nuts1[ ,aa_protnames] <- lapply(aa_protein_nuts1[ ,aa_protnames],
                                                 function(x) 
                                                   x / aa_protein_nuts1$SumProteinArea)  
    
    # nuts1-surface-weighted price for protein crops
    p_protein_nuts1 <- p_allcrops_nuts1[,c("year", "nuts1", p_protnames)] #subset prot
    p_protein_nuts1 <- left_join(p_protein_nuts1, 
                                 aa_protein_nuts1) # merge with area shares
    p_protein_nuts1[is.na(p_protein_nuts1)] <- 0 # replace missings
    attach(p_protein_nuts1)
    p_protein_nuts1$n1p_protein <- 
      p_beans * aa_beans + # note: aa_ are shares
      p_peas * aa_peas
    p_protein_nuts1 <- select(p_protein_nuts1, c("nuts1","year",n1p_protein))
    detach(p_protein_nuts1)


#------------------------------#
# Protein price at nuts0-level #
#------------------------------#

    # nuts0-land shares:
    
      #1) sum of individual protein items, and total protein area
      aa_protein_nuts0 <- df_farm %>%
        group_by(year) %>%
        summarise_at(all_of(aa_protnames), list(sum), na.rm = TRUE) %>% 
        rowwise() %>% 
        mutate(SumProteinArea = sum(c_across(all_of(aa_protnames)))) 
      
      #2) Share of individual protein items over total protein area
      aa_protein_nuts0[ ,aa_protnames] <- lapply(aa_protein_nuts0[ ,aa_protnames],
                                                 function(x) 
                                                   x / aa_protein_nuts0$SumProteinArea)  
      
    # nuts0-surface-weighted price for protein
    p_protein_nuts0 <- p_allcrops_nuts0[,c("year", p_protnames)] #subset prot
    p_protein_nuts0 <- left_join(p_protein_nuts0, 
                                 aa_protein_nuts0) # merge with area shares
    p_protein_nuts0[is.na(p_protein_nuts0)] <- 0 # replace missings
    attach(p_protein_nuts0)
    p_protein_nuts0$n0p_protein <- 
      p_beans * aa_beans + # note: aa_ are shares
      p_peas * aa_peas
    p_protein_nuts0 <- select(p_protein_nuts0,year,n0p_protein)
    detach(p_protein_nuts0)


#-------------------------#
###### 4c.3) Oilseed ######
#-------------------------#

  #------------------------------#
  # Oilseed price at nuts2-level #
  #------------------------------#

    # nuts2-land shares:
    
      #1) sum of individual oilseed items, and total oilseed area
      aa_oilseed_nuts2 <- df_farm %>%
        group_by(nuts2,year) %>%
        summarise_at(all_of(aa_oilnames), list(sum), na.rm = TRUE) %>% 
        rowwise() %>% 
        mutate(SumOilseedArea = sum(c_across(all_of(aa_oilnames)))) 
      
      #2) Share of individual oilseed items over total oilseed area
      aa_oilseed_nuts2[ ,aa_oilnames] <- lapply(aa_oilseed_nuts2[ ,aa_oilnames],
                                                function(x) 
                                                  x / aa_oilseed_nuts2$SumOilseedArea)  
    
    # nuts2-surface-weighted price for oilseed
    p_oilseed_nuts2 <- p_allcrops_nuts2[,c("year", "nuts2", p_oilnames)] #subset oil
    p_oilseed_nuts2 <- left_join(p_oilseed_nuts2, 
                                 aa_oilseed_nuts2) # merge with area shares
    p_oilseed_nuts2[is.na(p_oilseed_nuts2)] <- 0 # replace missings
    attach(p_oilseed_nuts2)
    p_oilseed_nuts2$n2p_oilseed <- 
      p_wcanola * aa_wcanola + # note: aa_ are shares
      p_scanola * aa_scanola + 
      p_othoil * aa_othoil + 
      p_fibre * aa_fibre + 
      p_flax * aa_flax +
      p_enoil * aa_enoil
    p_oilseed_nuts2 <- select(p_oilseed_nuts2, c("nuts2","year",n2p_oilseed))
    detach(p_oilseed_nuts2)

  #------------------------------#
  # Oilseed price at nuts1-level #
  #------------------------------#
    
    # nuts1-land shares:
    
      #1) sum of individual oilseed items, and total oilseed area
      aa_oilseed_nuts1 <- df_farm %>%
        group_by(nuts1,year) %>%
        summarise_at(all_of(aa_oilnames), list(sum), na.rm = TRUE) %>% 
        rowwise() %>% 
        mutate(SumOilseedArea = sum(c_across(all_of(aa_oilnames)))) 
      
      #2) Share of individual oilseed items over total oilseed area
      aa_oilseed_nuts1[ ,aa_oilnames] <- lapply(aa_oilseed_nuts1[ ,aa_oilnames],
                                                function(x) 
                                                  x / aa_oilseed_nuts1$SumOilseedArea)  
    
    # nuts1-surface-weighted price for oilseed
    p_oilseed_nuts1 <- p_allcrops_nuts1[,c("year", "nuts1", p_oilnames)] #subset oil
    p_oilseed_nuts1 <- left_join(p_oilseed_nuts1, 
                                 aa_oilseed_nuts1) # merge with area shares
    p_oilseed_nuts1[is.na(p_oilseed_nuts1)] <- 0 # replace missings
    attach(p_oilseed_nuts1)
    p_oilseed_nuts1$n1p_oilseed <- 
      p_wcanola * aa_wcanola + # note: aa_ are shares
      p_scanola * aa_scanola + 
      p_othoil * aa_othoil + 
      p_fibre * aa_fibre + 
      p_flax * aa_flax +
      p_enoil * aa_enoil
    p_oilseed_nuts1 <- select(p_oilseed_nuts1, c("nuts1","year",n1p_oilseed))
    detach(p_oilseed_nuts1)
    

  #------------------------------#
  # Oilseed price at nuts0-level #
  #------------------------------#

    # nuts0-land shares:
    
      #1) sum of individual oilseed items, and total oilseed area
      aa_oilseed_nuts0 <- df_farm %>%
        group_by(year) %>%
        summarise_at(all_of(aa_oilnames), list(sum), na.rm = TRUE) %>% 
        rowwise() %>% 
        mutate(SumOilseedArea = sum(c_across(all_of(aa_oilnames)))) 
      
      #2) Share of individual oilseed items over total oilseed area
      aa_oilseed_nuts0[ ,aa_oilnames] <- lapply(aa_oilseed_nuts0[ ,aa_oilnames],
                                                function(x) 
                                                  x / aa_oilseed_nuts0$SumOilseedArea)  
    
    # nuts0-surface-weighted price for oilseed
    p_oilseed_nuts0 <- p_allcrops_nuts0[,c("year", p_oilnames)] #subset oil
    p_oilseed_nuts0 <- left_join(p_oilseed_nuts0, 
                                 aa_oilseed_nuts0) # merge with area shares
    p_oilseed_nuts0[is.na(p_oilseed_nuts0)] <- 0 # replace missings
    attach(p_oilseed_nuts0)
    p_oilseed_nuts0$n0p_oilseed <- 
      p_wcanola * aa_wcanola + # note: aa_ are shares
      p_scanola * aa_scanola + 
      p_othoil * aa_othoil + 
      p_fibre * aa_fibre + 
      p_flax * aa_flax +
      p_enoil * aa_enoil
    p_oilseed_nuts0 <- select(p_oilseed_nuts0,year,n0p_oilseed)
    detach(p_oilseed_nuts0)
    
#----------------------------#
###### 4c.4) Root crops ######
#----------------------------#

  #---------------------------#
  # Root price at nuts2-level #
  #---------------------------#
    
    # nuts2-land shares:
    
      #1) sum of individual root items, and total root area
      aa_root_nuts2 <- df_farm %>%
        group_by(nuts2,year) %>%
        summarise_at(all_of(aa_rootnames), list(sum), na.rm = TRUE) %>% 
        rowwise() %>% 
        mutate(SumRootArea = sum(c_across(all_of(aa_rootnames)))) 
      
      #2) Share of individual root items over total root area
      aa_root_nuts2[ ,aa_rootnames] <- lapply(aa_root_nuts2[ ,aa_rootnames],
                                              function(x) 
                                                x / aa_root_nuts2$SumRootArea)  
    
    # nuts2-surface-weighted price for roots
    p_root_nuts2 <- p_allcrops_nuts2[,c("year", "nuts2", p_rootnames)] #subset root
    p_root_nuts2 <- left_join(p_root_nuts2, 
                              aa_root_nuts2) # merge with area shares
    p_root_nuts2[is.na(p_root_nuts2)] <- 0 # replace missings
    attach(p_root_nuts2)
    p_root_nuts2$n2p_root <- 
      p_potat * aa_potat + # note: aa_ are shares
      p_beets * aa_beets + 
      p_enbeets * aa_enbeets
    p_root_nuts2 <- select(p_root_nuts2, c("nuts2","year",n2p_root))
    detach(p_root_nuts2)

  #---------------------------#
  # Root price at nuts1-level #
  #---------------------------#
    
    # nuts1-land shares:
    
      #1) sum of individual root items, and total root area
      aa_root_nuts1 <- df_farm %>%
        group_by(nuts1,year) %>%
        summarise_at(all_of(aa_rootnames), list(sum), na.rm = TRUE) %>% 
        rowwise() %>% 
        mutate(SumRootArea = sum(c_across(all_of(aa_rootnames)))) 
      
      #2) Share of individual root items over total root area
      aa_root_nuts1[ ,aa_rootnames] <- lapply(aa_root_nuts1[ ,aa_rootnames],
                                              function(x) 
                                                x / aa_root_nuts1$SumRootArea)  
      
    # nuts1-surface-weighted price for cereals
    p_root_nuts1 <- p_allcrops_nuts1[,c("year", "nuts1", p_rootnames)] #subset root
    p_root_nuts1 <- left_join(p_root_nuts1, 
                              aa_root_nuts1) # merge with area shares
    p_root_nuts1[is.na(p_root_nuts1)] <- 0 # replace missings
    attach(p_root_nuts1)
    p_root_nuts1$n1p_root <- 
      p_potat * aa_potat + # note: aa_ are shares
      p_beets * aa_beets + 
      p_enbeets * aa_enbeets
    p_root_nuts1 <- select(p_root_nuts1, c("nuts1","year",n1p_root))
    detach(p_root_nuts1)


  #---------------------------#
  # Root price at nuts0-level #
  #---------------------------#
    
    # nuts0-land shares:
    
      #1) sum of individual root items, and total root area
      aa_root_nuts0 <- df_farm %>%
        group_by(year) %>%
        summarise_at(all_of(aa_rootnames), list(sum), na.rm = TRUE) %>% 
        rowwise() %>% 
        mutate(SumRootArea = sum(c_across(all_of(aa_rootnames)))) 
      
      #2) Share of individual root items over total root area
      aa_root_nuts0[ ,aa_rootnames] <- lapply(aa_root_nuts0[ ,aa_rootnames],
                                              function(x) 
                                                x / aa_root_nuts0$SumRootArea)  
      
    # nuts0-surface-weighted price for root
    p_root_nuts0 <- p_allcrops_nuts0[,c("year", p_rootnames)] #subset root
    p_root_nuts0 <- left_join(p_root_nuts0, 
                              aa_root_nuts0) # merge with area shares
    p_root_nuts0[is.na(p_root_nuts0)] <- 0 # replace missings
    attach(p_root_nuts0)
    p_root_nuts0$n0p_root <- 
      p_potat * aa_potat + # note: aa_ are shares
      p_beets * aa_beets + 
      p_enbeets * aa_enbeets
    p_root_nuts0 <- select(p_root_nuts0,year,n0p_root)
    detach(p_root_nuts0)
    

#-----------------------------------------#
###### 4c.5) Corn (grain and silage) ######
#-----------------------------------------#
    
  #---------------------------#
  # Corn price at nuts2-level #
  #---------------------------#
    
    # nuts2-land shares:
    
    #1) sum of individual corn items, and total corn area
    aa_corn_nuts2 <- df_farm %>%
      group_by(nuts2,year) %>%
      summarise_at(all_of(aa_cornnames), list(sum), na.rm = TRUE) %>% 
      rowwise() %>% 
      mutate(SumCornArea = sum(c_across(all_of(aa_cornnames)))) 
    
    #2) Share of individual corn items over total corn area
    aa_corn_nuts2[ ,aa_cornnames] <- lapply(aa_corn_nuts2[ ,aa_cornnames],
                                            function(x) 
                                              x / aa_corn_nuts2$SumCornArea)  
    
    # nuts2-surface-weighted price for corn
    p_corn_nuts2 <- p_allcrops_nuts2[,c("year", "nuts2", p_cornnames)] #subset corn
    p_corn_nuts2 <- left_join(p_corn_nuts2, 
                              aa_corn_nuts2) # merge with area shares
    p_corn_nuts2[is.na(p_corn_nuts2)] <- 0 # replace missings
    attach(p_corn_nuts2)
    p_corn_nuts2$n2p_corn <- 
      p_maize * aa_maize + # note: aa_ are shares
      p_corn * aa_corn
    p_corn_nuts2 <- select(p_corn_nuts2, c("nuts2","year",n2p_corn))
    detach(p_corn_nuts2)
    
  #---------------------------#
  # Corn price at nuts1-level #
  #---------------------------#
    
    # nuts1-land shares:
    
    #1) sum of individual corn items, and total corn area
    aa_corn_nuts1 <- df_farm %>%
      group_by(nuts1,year) %>%
      summarise_at(all_of(aa_cornnames), list(sum), na.rm = TRUE) %>% 
      rowwise() %>% 
      mutate(SumCornArea = sum(c_across(all_of(aa_cornnames)))) 
    
    #2) Share of individual corn items over total corn area
    aa_corn_nuts1[ ,aa_cornnames] <- lapply(aa_corn_nuts1[ ,aa_cornnames],
                                            function(x) 
                                              x / aa_corn_nuts1$SumCornArea)  
    
    # nuts1-surface-weighted price for corns
    p_corn_nuts1 <- p_allcrops_nuts1[,c("year", "nuts1", p_cornnames)] #subset corn
    p_corn_nuts1 <- left_join(p_corn_nuts1, 
                              aa_corn_nuts1) # merge with area shares
    p_corn_nuts1[is.na(p_corn_nuts1)] <- 0 # replace missings
    attach(p_corn_nuts1)
    p_corn_nuts1$n1p_corn <- 
      p_maize * aa_maize + # note: aa_ are shares
      p_corn * aa_corn
    p_corn_nuts1 <- select(p_corn_nuts1, c("nuts1","year",n1p_corn))
    detach(p_corn_nuts1)
    
    
  #---------------------------#
  # Corn price at nuts0-level #
  #---------------------------#
    
    # nuts0-land shares:
    
    #1) sum of individual corn items, and total corn area
    aa_corn_nuts0 <- df_farm %>%
      group_by(year) %>%
      summarise_at(all_of(aa_cornnames), list(sum), na.rm = TRUE) %>% 
      rowwise() %>% 
      mutate(SumCornArea = sum(c_across(all_of(aa_cornnames)))) 
    
    #2) Share of individual corn items over total corn area
    aa_corn_nuts0[ ,aa_cornnames] <- lapply(aa_corn_nuts0[ ,aa_cornnames],
                                            function(x) 
                                              x / aa_corn_nuts0$SumCornArea)  
    
    # nuts0-surface-weighted price for corns
    p_corn_nuts0 <- p_allcrops_nuts0[,c("year", p_cornnames)] #subset corns
    p_corn_nuts0 <- left_join(p_corn_nuts0, 
                              aa_corn_nuts0) # merge with area shares
    p_corn_nuts0[is.na(p_corn_nuts0)] <- 0 # replace missings
    attach(p_corn_nuts0)
    p_corn_nuts0$n0p_corn <- 
      p_maize * aa_maize + # note: aa_ are shares
      p_corn * aa_corn
    p_corn_nuts0 <- select(p_corn_nuts0,year,n0p_corn)
    detach(p_corn_nuts0)       

#--------------------------------------------------------#
#### 5) Merge regional price data to farm-level data ####
#--------------------------------------------------------#
    
#------------------------------------------#
##### 5a) Join all regional price data #####
#-----------------------------------------#
  
  nuts0_p <- p_cereals_nuts0 %>% 
    left_join(p_protein_nuts0) %>% 
    left_join(p_oilseed_nuts0) %>% 
    left_join(p_root_nuts0) %>% 
    left_join(p_corn_nuts0)
  
  nuts1_p <- p_cereals_nuts1 %>% 
    left_join(p_protein_nuts1) %>% 
    left_join(p_oilseed_nuts1) %>% 
    left_join(p_root_nuts1) %>% 
    left_join(p_corn_nuts1)
  
  nuts2_p <- p_cereals_nuts2 %>% 
    left_join(p_protein_nuts2) %>% 
    left_join(p_oilseed_nuts2) %>% 
    left_join(p_root_nuts2) %>% 
    left_join(p_corn_nuts2)

  #replace zero prices with missings
  nuts0_p[nuts0_p == 0] <- NA
  nuts1_p[nuts1_p == 0] <- NA
  nuts2_p[nuts2_p == 0] <- NA
  
#----------------------------------#
##### 5b) Create lagged prices #####
#----------------------------------#
  
  #nuts0
  nuts0_p$key <- 1 #(auxiliary var; does not work without it; will be deleted)
  nuts0_p <- plm::pdata.frame(nuts0_p,index=c("key","year")) 
  nuts0_p$n0lp_cereals <- plm::lag(nuts0_p$n0p_cereals)
  nuts0_p$n0lp_oilseed <- plm::lag(nuts0_p$n0p_oilseed)
  nuts0_p$n0lp_protein <- plm::lag(nuts0_p$n0p_protein)
  nuts0_p$n0lp_root    <- plm::lag(nuts0_p$n0p_root)
  nuts0_p$n0lp_corn    <- plm::lag(nuts0_p$n0p_corn)
  nuts0_p$key <- NULL
  nuts0_p <- data.frame(nuts0_p)
  
  nuts0_p[] <- lapply(nuts0_p, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
  })
  sapply(nuts0_p, class)
  
  #nuts1
  nuts1_p <- plm::pdata.frame(nuts1_p,index=c("nuts1","year")) 
  nuts1_p$n1lp_cereals <- plm::lag(nuts1_p$n1p_cereals)
  nuts1_p$n1lp_oilseed <- plm::lag(nuts1_p$n1p_oilseed)
  nuts1_p$n1lp_protein <- plm::lag(nuts1_p$n1p_protein)
  nuts1_p$n1lp_root    <- plm::lag(nuts1_p$n1p_root)
  nuts1_p$n1lp_corn    <- plm::lag(nuts1_p$n1p_corn)
  nuts1_p <- data.frame(nuts1_p)
  
  nuts1_p[] <- lapply(nuts1_p, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
  })
  sapply(nuts1_p, class)
  
  #nuts2
  nuts2_p <- plm::pdata.frame(nuts2_p,index=c("nuts2","year")) 
  nuts2_p$n2lp_cereals <- plm::lag(nuts2_p$n2p_cereals)
  nuts2_p$n2lp_oilseed <- plm::lag(nuts2_p$n2p_oilseed)
  nuts2_p$n2lp_protein <- plm::lag(nuts2_p$n2p_protein)
  nuts2_p$n2lp_root <- plm::lag(nuts2_p$n2p_root)
  nuts2_p$n2lp_corn    <- plm::lag(nuts2_p$n2p_corn)
  nuts2_p <- data.frame(nuts2_p)
  
  nuts2_p[] <- lapply(nuts2_p, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
  })
  sapply(nuts2_p, class)


#----------------------------------------------------------#
##### 5c) Merge regional price data to farm-level data #####
#----------------------------------------------------------#
  
df_farm$nuts1 <- as.numeric(df_farm$nuts1)
df_farm$nuts2 <- as.numeric(df_farm$nuts2)

df_farm <- left_join(df_farm, 
                     nuts0_p)

df_farm <- left_join(df_farm, 
                     nuts1_p)

df_farm <- left_join(df_farm, 
                     nuts2_p)


#-------------------------------------#
#### 6) Quantities of fixed inputs ####
#-------------------------------------#

#---------------------------#
##### 6a) Land quantity #####
#---------------------------#

df_farm$x_land <- df_farm$z6119s07

#----------------------------#
##### 6b) Labor quantity #####
#----------------------------#

df_farm$x_labor <- df_farm$z7089s03 + df_farm$z7098s03 #fam. + hired labor
summary(df_farm$x_labor)

#------------------------------#
##### 6c) Capital quantity #####
#------------------------------#

# Note: Proxied by depreciation, and deflated by investment costs
df_farm$x_capital <- -df_farm$z2809s05 / as.numeric(df_farm$wi_invest)

#------------------------------------#
#### 7) Prices of variable inputs ####
#------------------------------------#

#------------------------------------------------------#
#### 7a) Price of material inputs except fertilizer ####
#------------------------------------------------------#

# Set missings values for material inputs expenditures to zero
df_farm$z2509s02[is.na(df_farm$z2509s02)] <- 0 #seed
df_farm$z2559s02[is.na(df_farm$z2559s02)] <- 0 #Plant protection
df_farm$z2570s02[is.na(df_farm$z2570s02)] <- 0 #Substrates
df_farm$z2591s02[is.na(df_farm$z2591s02)] <- 0 #Foil
df_farm$z2597s02[is.na(df_farm$z2597s02)] <- 0 #Other material for plant production
df_farm$z2598s02[is.na(df_farm$z2598s02)] <- 0 #Services for plant production
df_farm$z2705s02[is.na(df_farm$z2705s02)] <- 0 #Concentrates for cattle
df_farm$z2709s02[is.na(df_farm$z2709s02)] <- 0 #Concentrates for pigs
df_farm$z2712s02[is.na(df_farm$z2712s02)] <- 0 #Feed for poultry
df_farm$z2714s02[is.na(df_farm$z2714s02)] <- 0 #Feed for sheep and goats
df_farm$z2718s02[is.na(df_farm$z2718s02)] <- 0 #Feed for other animals
df_farm$z2706s02[is.na(df_farm$z2706s02)] <- 0 #other feed for cattle
df_farm$z2715s02[is.na(df_farm$z2715s02)] <- 0 #other feed for sheep and goats
df_farm$z2720s02[is.na(df_farm$z2720s02)] <- 0 #semen
df_farm$z2721s02[is.na(df_farm$z2721s02)] <- 0 #vet
df_farm$z2722s02[is.na(df_farm$z2722s02)] <- 0 #desinfection
df_farm$z2726s02[is.na(df_farm$z2726s02)] <- 0 #other material for livestock production
df_farm$z2728s02[is.na(df_farm$z2728s02)] <- 0 #services for livestock production
df_farm$z2770s02[is.na(df_farm$z2770s02)] <- 0 #heating
df_farm$z2771s02[is.na(df_farm$z2771s02)] <- 0 #electricity
df_farm$z2773s02[is.na(df_farm$z2773s02)] <- 0 #fuel
df_farm$z2774s02[is.na(df_farm$z2774s02)] <- 0 #lubricants
df_farm$z2772s02[is.na(df_farm$z2772s02)] <- 0 #water
df_farm$z2782s02[is.na(df_farm$z2782s02)] <- 0 #services  

# Costs for crop specific inputs (except fertilizer)
df_farm$c_cropspec <- -(df_farm$z2509s02+
                          df_farm$z2559s02+df_farm$z2570s02+df_farm$z2591s02+
                          df_farm$z2597s02+df_farm$z2598s02)

# Costs for animal-specific inputs
df_farm$c_animspec <- -(df_farm$z2705s02+df_farm$z2709s02+df_farm$z2712s02+
                          df_farm$z2714s02+df_farm$z2718s02+df_farm$z2706s02+
                          df_farm$z2715s02+df_farm$z2720s02+df_farm$z2721s02+
                          df_farm$z2722s02+df_farm$z2726s02+df_farm$z2728s02)

# Costs for intermediate inputs
df_farm$c_intermediates <- -(df_farm$z2770s02+df_farm$z2771s02+df_farm$z2773s02+
                               df_farm$z2774s02+df_farm$z2772s02+df_farm$z2782s02)

#load destatis prices  
load("rOutput/wi_nuts0.Rda")

#Note: Prices are already matched at the farm-level, but it is easier
#      to compute the regional price indexes at the regional level first

# Tornquist price index at the nuts2-level

# select all material input-variables except fertilizer
vec_otherinp  <- c( "z2509s02", #Seed
                    "z2559s02", #Plant protection
                    "z2570s02", #Substrates
                    "z2591s02", #Foil
                    "z2597s02", #Other material for plant production
                    "z2598s02", #Services for plant production
                    "z2705s02", #Concentrates for cattle
                    "z2709s02", #Concentrates for pigs
                    "z2712s02", #Feed for poultry
                    "z2714s02", #Feed for sheep and goats
                    "z2718s02", #Feed for other animals
                    "z2706s02", #other feed for cattle
                    "z2715s02", #other feed for sheep and goats
                    "z2720s02", #semen
                    "z2721s02", #vet
                    "z2722s02", #desinfection
                    "z2726s02", #other material for livestock production
                    "z2728s02", #services for livestock production
                    "z2770s02", #heating
                    "z2771s02", #electricity
                    "z2773s02", #fuel
                    "z2774s02", #lubricants
                    "z2772s02", #water
                    "z2782s02") #services

# compute yearly expenditures at regional level (incl. total: c_otherinp) 
n2w_tq_otherinp <- df_farm %>%
  group_by(nuts2,year) %>% 
  summarise_at(vec_otherinp, funs(sum), na.rm = TRUE) %>% 
  rowwise() %>% 
  mutate(c_otherinp = sum(c_across(vec_otherinp))) %>% 
  left_join(wi_nuts0) # adds the country-level price indexes

# add the newly calculated "c_otherinp" to the variable list
vec_otherinp <- c(vec_otherinp,"c_otherinp")

# add columns with the base year values (here: the first year)
n2w_tq_otherinp <- n2w_tq_otherinp %>% 
  group_by(nuts2) %>% 
  summarise_at(all_of(vec_otherinp), funs(t0=first)) %>% 
  left_join(n2w_tq_otherinp)

# Apply the Tornquist formula
n2w_tq_otherinp$n2w_tq_otherinp <- exp( 0.5 * 
                                          (   ( n2w_tq_otherinp$z2509s02 / n2w_tq_otherinp$c_otherinp +
                                                  n2w_tq_otherinp$z2509s02_t0 / n2w_tq_otherinp$c_otherinp_t0 ) 
                                              * log(n2w_tq_otherinp$wi_seed) + # seed
                                                ( n2w_tq_otherinp$z2559s02 / n2w_tq_otherinp$c_otherinp +
                                                    n2w_tq_otherinp$z2559s02_t0 / n2w_tq_otherinp$c_otherinp_t0 ) 
                                              * log(n2w_tq_otherinp$wi_pest) + # plant protection
                                                ( (n2w_tq_otherinp$z2570s02 + n2w_tq_otherinp$z2591s02 + 
                                                     n2w_tq_otherinp$z2597s02 + n2w_tq_otherinp$z2598s02) / n2w_tq_otherinp$c_otherinp +
                                                    (n2w_tq_otherinp$z2570s02_t0 + n2w_tq_otherinp$z2591s02_t0 + 
                                                       n2w_tq_otherinp$z2597s02_t0 + n2w_tq_otherinp$z2598s02_t0) / n2w_tq_otherinp$c_otherinp_t0 ) 
                                              * log(n2w_tq_otherinp$wi_inputs) + # other crop-specific inputs
                                                ( (n2w_tq_otherinp$z2705s02 + n2w_tq_otherinp$z2709s02 +
                                                     n2w_tq_otherinp$z2712s02 + n2w_tq_otherinp$z2714s02 +
                                                     n2w_tq_otherinp$z2718s02 + n2w_tq_otherinp$z2706s02 +
                                                     n2w_tq_otherinp$z2715s02) / n2w_tq_otherinp$c_otherinp +
                                                    (n2w_tq_otherinp$z2705s02_t0 + n2w_tq_otherinp$z2709s02_t0 +
                                                       n2w_tq_otherinp$z2712s02_t0 + n2w_tq_otherinp$z2714s02_t0 +
                                                       n2w_tq_otherinp$z2718s02_t0 + n2w_tq_otherinp$z2706s02_t0 +
                                                       n2w_tq_otherinp$z2715s02_t0) / n2w_tq_otherinp$c_otherinp_t0 ) 
                                              * log(n2w_tq_otherinp$wi_feed)  + # feed
                                                ( (n2w_tq_otherinp$z2720s02 + n2w_tq_otherinp$z2721s02) / n2w_tq_otherinp$c_otherinp +
                                                    (n2w_tq_otherinp$z2720s02_t0 + n2w_tq_otherinp$z2721s02_t0) / n2w_tq_otherinp$c_otherinp_t0 ) 
                                              * log(n2w_tq_otherinp$wi_vet) + # veterinary cost
                                                ( (n2w_tq_otherinp$z2722s02 + n2w_tq_otherinp$z2726s02 + 
                                                     n2w_tq_otherinp$z2728s02) / n2w_tq_otherinp$c_otherinp +
                                                    (n2w_tq_otherinp$z2722s02_t0 + n2w_tq_otherinp$z2726s02_t0 + 
                                                       n2w_tq_otherinp$z2728s02_t0) / n2w_tq_otherinp$c_otherinp_t0 ) 
                                              * log(n2w_tq_otherinp$wi_inputs) + # other anim-specific inputs
                                                ( n2w_tq_otherinp$z2770s02 / n2w_tq_otherinp$c_otherinp +
                                                    n2w_tq_otherinp$z2770s02_t0 / n2w_tq_otherinp$c_otherinp_t0 ) 
                                              * log(n2w_tq_otherinp$wi_heating) + # heating
                                                ( n2w_tq_otherinp$z2771s02 / n2w_tq_otherinp$c_otherinp +
                                                    n2w_tq_otherinp$z2771s02_t0 / n2w_tq_otherinp$c_otherinp_t0 ) 
                                              * log(n2w_tq_otherinp$wi_electricity)  + # electricity
                                                ( n2w_tq_otherinp$z2773s02 / n2w_tq_otherinp$c_otherinp +
                                                    n2w_tq_otherinp$z2773s02_t0 / n2w_tq_otherinp$c_otherinp_t0 ) 
                                              * log(n2w_tq_otherinp$wi_fuel)  + # fuel
                                                ( n2w_tq_otherinp$z2774s02 / n2w_tq_otherinp$c_otherinp +
                                                    n2w_tq_otherinp$z2774s02_t0 / n2w_tq_otherinp$c_otherinp_t0 ) 
                                              * log(n2w_tq_otherinp$wi_lubricants)  + # lubricants
                                                ( (n2w_tq_otherinp$z2772s02 + n2w_tq_otherinp$z2782s02) / n2w_tq_otherinp$c_otherinp +
                                                    (n2w_tq_otherinp$z2772s02_t0 + n2w_tq_otherinp$z2782s02_t0) / n2w_tq_otherinp$c_otherinp_t0 ) 
                                              * log(n2w_tq_otherinp$wi_mat_serv) ) ) # water and services

#--------------------------------------#
##### 7b) EUR-price for fertilizer #####
#--------------------------------------#

# load fertilizer data
library(readxl) #for opening excel files
n0w_fert <- read_excel("Data/Fertilizer prices and quantities_1992-2020.xlsx", 
                       col_names = TRUE)

# Calculate the share-weighted price for 1 ton of pure nutrient 
n0w_fert$q_fertilizer <- n0w_fert$q_stickstoff + n0w_fert$q_phosphat + 
  n0w_fert$q_kali + n0w_fert$q_kalk
n0w_fert$sh_stickstoff <- n0w_fert$q_stickstoff / n0w_fert$q_fertilizer
n0w_fert$sh_phosphat <- n0w_fert$q_phosphat / n0w_fert$q_fertilizer
n0w_fert$sh_kali <- n0w_fert$q_kali / n0w_fert$q_fertilizer
n0w_fert$sh_kalk <-  n0w_fert$q_kalk / n0w_fert$q_fertilizer

n0w_fert$n0w_fert =  n0w_fert$w_stickstoff * n0w_fert$sh_stickstoff +
  n0w_fert$w_phosphat * n0w_fert$sh_phosphat +
  n0w_fert$w_kali * n0w_fert$sh_kali +
  n0w_fert$w_kalk * n0w_fert$sh_kalk

# Set missing values for fertilizer expenditures to zero
df_farm$z2539s02[is.na(df_farm$z2539s02)] <- 0 

#-------------------------------------------------------#
#### 8) Merge regional price data to farm-level data ####
#-------------------------------------------------------#

# select as pre-step
n2w_tq_otherinp_select  <- n2w_tq_otherinp %>%  select(nuts2, year, n2w_tq_otherinp)
n0w_fert_select         <- n0w_fert %>%  select(year, n0w_fert)

#summarize input price indexes in one data frame
n2w <- n2w_tq_materials_select %>%  
  left_join(n0w_fert_select) %>%
  dplyr::select(nuts2, year,
                n2w_tq_otherinp,n0w_fert)

#append this price data to the farm-level data
df_farm <- left_join(df_farm, 
                     n2w)

#-------------------------------#
#                               #
# End of farm data preparation  # 
#                               #
#-------------------------------#

# save data frame
save(df_farm,file="rOutput/farm_outp_inp.Rda")




