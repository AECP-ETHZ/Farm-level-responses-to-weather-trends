#---------------------------------------------------#
#                                                   #
# This program organizes the German-level input and # 
# output price indexes from Destatis / Genesis      #
#                                                   #
#---------------------------------------------------#

# Open packages
library(dplyr)
library(readxl) 
library(cellranger) 
library(reshape2)

# load data, melt, rename
wi_raw <- read_excel("Data/destatis_inputs_61221-0002.xlsx", sheet = "61221-0002", range = "A6:BU40", col_names = TRUE)
wi_raw[wi_raw=="-"]<-NA
wi <- melt(wi_raw, id.vars="...1")
wi <- dplyr::rename(wi, c("input"="...1", "year"="variable", "wi"="value"))

#reshape
wi_nuts0 <- reshape(wi, idvar = "year", timevar = "input", direction = "wide")

#rename columns
colnames(wi_nuts0)
wi_nuts0 <- dplyr::rename(wi_nuts0, 
  c("wi_inputs"="wi.Landwirtschaftliche Betriebsmittel insgesamt",
    "wi_mat_serv"="wi.Waren und Dienstleist. des lfd. landw. Verbrauchs",
    "wi_seed"="wi.Saat- und Pflanzgut",
    "wi_energy"="wi.Energie und Schmierstoffe",
    "wi_heating"="wi.Heizstoffe",
    "wi_fuel"="wi.Treibstoffe",
    "wi_electricity"="wi.Elektrischer Strom",
    "wi_lubricants"="wi.Schmierstoffe",
    "wi_fert"="wi.Düngemittel",                                      
    "wi_pest"="wi.Pflanzenschutzmittel",
    "wi_fung"="wi.Fungizide",                                        
    "wi_insect"="wi.Insektizide",
    "wi_herb"="wi.Herbizide",                                        
    "wi_feed"="wi.Futtermittel",
    "wi_singlefeed"="wi.Einzelfuttermittel",                               
    "wi_wheatprod"="wi.Getreide und Mühlennachprodukte",
    "wi_oilkake"="wi.Ölkuchen und -schrot",                             
    "wi_compoundfeed"="wi.Mischfuttermittel", 
    "wi_compf_cattle"="wi.Mischfuttermittel für Rinder",                     
    "wi_compf_hogs"="wi.Mischfuttermittel für Schweine", 
    "wi_compf_poultry"="wi.Mischfuttermittel für Geflügel",                   
    "wi_vet"="wi.Veterinärleistungen", 
    "wi_maint_mach_mat"="wi.Instandhaltung von Maschinen und Material",        
    "wi_maint_build"="wi.Instandhaltung von Bauten", 
    "wi_othermatserv"="wi.Sonstige Waren und Dienstleistungen",              
    "wi_invest"="wi.Waren und Dienstleist. landwirt. Investitionen", 
    "wi_material"="wi.Material",                                         
    "wi_machinery"="wi.Maschinen und sonstige Ausrüstungsgüter",
    "wi_cropsmachinery"="wi.Maschinen und Geräte für Kulturen",                
    "wi_harvesters"="wi.Maschinen und Geräte für die Erntebergung",
    "wi_vehicles"="wi.Fahrzeuge",                                        
    "wi_tractors"="wi.Zugmaschinen",
    "wi_othervehicles"="wi.Sonstige Fahrzeuge",                               
    "wi_buildings"="wi.Bauten"))


#-------------------------------------------#
# Translate marketing years to actual years #
#-------------------------------------------#

# Note: The variable "year" in the farm data refers to the second part of the marketing year
#       (i.e., 2018 means marketing year 2017/18)

#take digits 1 and 2 as well as digits 6 and 7 from the original variable
wi_nuts0$year <- paste0(substr(wi_nuts0$year, 1, 2), 
                        substr(wi_nuts0$year, 6, 7))

#Correct the 1900 (from 1999/00)
wi_nuts0$year[wi_nuts0$year==1900] <- 2000

# Make data numeric
wi_nuts0 <- as.data.frame(lapply(wi_nuts0,as.numeric))

# save price data
save(wi_nuts0,file="rOutput/wi_nuts0.Rda")


															