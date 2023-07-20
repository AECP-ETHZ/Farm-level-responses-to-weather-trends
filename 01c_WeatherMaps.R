#------------------------------------------#
#                                          #
# This program creates the weather plots.  #
#                                          #
#------------------------------------------#

# load packages
library(dplyr)
library(sf)
library(ggplot2)
library(patchwork)

# load weather data
load("rOutput/df_weather.rda" )

# load shapefiles
shapefile_gem2 <- st_read("Shapefiles/gemeinde2007_bkg.shp") #municipalities
shapefile_bl <-st_read("Shapefiles/bdl2018_1.shp") #states (for borders)

#--------------------------------------------#
#### Prepare shape files and weather data ####
#--------------------------------------------#

#Merge weather data with shapefiles
shapefile_gem_weather <- shapefile_gem2 %>% left_join(df_weather,by=c("AGS"="nuts5"))

# Set limits for illustration of weather events
  
  #gdd
  shapefile_gem_weather$gdd_plot <- cut(shapefile_gem_weather$gdd, breaks=c(0, 600, 700, 800 , 900, 1000, 1100, Inf),
                                                labels = c('<= 600 ', '>600-700','>700-800','>800-900','>900-1000', '> 1000-1100', '>1100'), include.lowest = TRUE)
  shapefile_gem_weather$gdd_1to3_plot <- cut(shapefile_gem_weather$gdd_1to3, breaks=c(0, 600, 700, 800 , 900, 1000, 1100, Inf),
                                        labels = c('<= 600 ', '>600-700','>700-800','>800-900','>900-1000', '> 1000-1100', '>1100'), include.lowest = TRUE)
  shapefile_gem_weather$gdd_4to10_plot <- cut(shapefile_gem_weather$gdd_4to10, breaks=c(0, 600, 700, 800 , 900, 1000, 1100, Inf),
                                        labels = c('<= 600 ', '>600-700','>700-800','>800-900','>900-1000', '> 1000-1100', '>1100'), include.lowest = TRUE)
  
  #prec
  shapefile_gem_weather$prec_plot <- cut(shapefile_gem_weather$prec, breaks=c(0,200, 250,300,350, 400, 450, Inf),
                                                   labels = c('<= 200 ', '>200-250','>250-300','>300-350', '>350-400', '>400-450', '>450'), include.lowest = TRUE)
  shapefile_gem_weather$prec_1to3_plot <- cut(shapefile_gem_weather$prec_1to3, breaks=c(0,200, 250,300,350, 400, 450, Inf),
                                         labels = c('<= 200 ', '>200-250','>250-300','>300-350', '>350-400', '>400-450', '>450'), include.lowest = TRUE)
  shapefile_gem_weather$prec_4to10_plot <- cut(shapefile_gem_weather$prec_4to10, breaks=c(0,200, 250,300,350, 400, 450, Inf),
                                         labels = c('<= 200 ', '>200-250','>250-300','>300-350', '>350-400', '>400-450', '>450'), include.lowest = TRUE)
  
  
  #gddHigh
  shapefile_gem_weather$gddHigh_plot <- cut(shapefile_gem_weather$gddHigh, breaks=c(0,1,3,5,7,9, Inf),
                                                      labels = c('<= 1', '>1-3','>3-5','>5-7', '>7-9', '>9'), include.lowest = TRUE)
  shapefile_gem_weather$gddHigh_1to3_plot <- cut(shapefile_gem_weather$gddHigh_1to3, breaks=c(0,1,3,5,7,9, Inf),
                                            labels = c('<= 1', '>1-3','>3-5','>5-7', '>7-9', '>9'), include.lowest = TRUE)
  shapefile_gem_weather$gddHigh_4to10_plot <- cut(shapefile_gem_weather$gddHigh_4to10, breaks=c(0,1,3,5,7,9, Inf),
                                            labels = c('<= 1', '>1-3','>3-5','>5-7', '>7-9', '>9'), include.lowest = TRUE)
  
  #dd
  shapefile_gem_weather$dd_plot <- cut(shapefile_gem_weather$dd, breaks=c(0, 60,80, 100 ,120,140, Inf),
                                                  labels = c( '<=60', '>60-80','>80-100','>100-120',  '>120-140', '>140'), include.lowest = TRUE)
  shapefile_gem_weather$dd_1to3_plot <- cut(shapefile_gem_weather$dd_1to3, breaks=c(0, 60,80, 100 ,120,140, Inf),
                                       labels = c( '<=60', '>60-80','>80-100','>100-120',  '>120-140', '>140'), include.lowest = TRUE)
  shapefile_gem_weather$dd_4to10_plot <- cut(shapefile_gem_weather$dd_4to10, breaks=c(0, 60,80, 100 ,120,140, Inf),
                                       labels = c( '<=60', '>60-80','>80-100','>100-120',  '>120-140', '>140'), include.lowest = TRUE)

# Select years for illustration 
shapefile_gem_weather_1996 <- shapefile_gem_weather[which(shapefile_gem_weather$year == '1996'),]
shapefile_gem_weather_2008 <- shapefile_gem_weather[which(shapefile_gem_weather$year == '2008'),]
shapefile_gem_weather_2019 <- shapefile_gem_weather[which(shapefile_gem_weather$year == '2019'),]


# --------------------------------------------- #
#### Figure 1: Growing degree days variables ####
# --------------------------------------------- #

#-------------------------#
# Row 1: Realized weather #
#-------------------------#

# 1996
plot_gdd_1996 <- ggplot() +
  geom_sf(aes(fill = gdd_plot), color = 'transparent', data = shapefile_gem_weather_1996) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlGn", name = 'Growing Degree Days (10-30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


# 2008
plot_gdd_2008 <- ggplot() +
  geom_sf(aes(fill = gdd_plot), color = 'transparent', data = shapefile_gem_weather_2008) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlGn", name = 'Growing Degree Days (10-30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


#2019
plot_gdd_2019 <- ggplot() +
  geom_sf(aes(fill = gdd_plot), color = 'transparent', data = shapefile_gem_weather_2019) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlGn", name = 'Growing Degree Days (10-30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')

#---------------------#
# Row 2: Weather 1to3 #
#---------------------#

# 1996
plot_gdd1to3_1996 <- ggplot() +
  geom_sf(aes(fill = gdd_1to3_plot), color = 'transparent', data = shapefile_gem_weather_1996) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlGn", name = 'Growing Degree Days (10-30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


# 2008
plot_gdd1to3_2008 <- ggplot() +
  geom_sf(aes(fill = gdd_1to3_plot), color = 'transparent', data = shapefile_gem_weather_2008) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlGn", name = 'Growing Degree Days (10-30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


#2019
plot_gdd1to3_2019 <- ggplot() +
  geom_sf(aes(fill = gdd_1to3_plot), color = 'transparent', data = shapefile_gem_weather_2019) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlGn", name = 'Growing Degree Days (10-30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')

#----------------------#
# Row 3: Weather 4to10 #
#----------------------#

# 1996
plot_gdd4to10_1996 <- ggplot() +
  geom_sf(aes(fill = gdd_4to10_plot), color = 'transparent', data = shapefile_gem_weather_1996) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlGn", name = 'Growing Degree Days (10-30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


# 2008
plot_gdd4to10_2008 <- ggplot() +
  geom_sf(aes(fill = gdd_4to10_plot), color = 'transparent', data = shapefile_gem_weather_2008) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlGn", name = 'Growing Degree Days (10-30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


#2019
plot_gdd4to10_2019 <- ggplot() +
  geom_sf(aes(fill = gdd_4to10_plot), color = 'transparent', data = shapefile_gem_weather_2019) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlGn", name = 'Growing Degree Days (10-30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')



#----------------------------------#
# Create one figure with all plots #
#----------------------------------#

#create row and column labels
row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="GDD", angle = 90, size = 5) + theme_void()
row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="GDD1to3", angle = 90, size = 5) + theme_void() 
row3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="GDD4to10", angle = 90, size = 5) + theme_void()  

col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="1996", size = 5) + theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="2008", size = 5) + theme_void() 
col3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="2019", size = 5) + theme_void() 

#define plot list
plotlist <- list(a = row1, b = row2, c = row3, d = col1, e = col2, f = col3, 
                 g= plot_gdd_1996, h=plot_gdd_2008, i=plot_gdd_2019, 
                 j= plot_gdd1to3_1996, k=plot_gdd1to3_2008, l=plot_gdd1to3_2019,
                 m= plot_gdd4to10_1996, n=plot_gdd4to10_2008, o=plot_gdd4to10_2019)

#define plot arrangement (repeat letters for scaling)
layoutplot <- "
    #ddddddeeeeeeffffff
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    "

#Create and save plot
fig_gdd <- patchwork::wrap_plots(plotlist, guides = 'collect', design = layoutplot) & theme(legend.position = 'bottom')

ggsave("Figures/Figure_1.png", fig_gdd, device="png", width = 8, height = 8, units = "in", dpi=1200)
ggsave("Figures/Figure_1.eps", fig_gdd, device="eps", width = 8, height = 8, units = "in", dpi=1200)

# --------------------------------------- #
#### Figure 2: Precipitation variables ####
# --------------------------------------- #

#-------------------------#
# Row 1: Realized weather #
#-------------------------#

# 1996
plot_prec_1996 <- ggplot() +
  geom_sf(aes(fill = prec_plot), color = 'transparent', data = shapefile_gem_weather_1996) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "Blues", name = 'Precipitation (mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


# 2008
plot_prec_2008 <- ggplot() +
  geom_sf(aes(fill = prec_plot), color = 'transparent', data = shapefile_gem_weather_2008) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "Blues", name = 'Precipitation (mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


#2019
plot_prec_2019 <- ggplot() +
  geom_sf(aes(fill = prec_plot), color = 'transparent', data = shapefile_gem_weather_2019) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "Blues", name = 'Precipitation (mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')

#---------------------#
# Row 2: Weather 1to3 #
#---------------------#

# 1996
plot_prec1to3_1996 <- ggplot() +
  geom_sf(aes(fill = prec_1to3_plot), color = 'transparent', data = shapefile_gem_weather_1996) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "Blues", name = 'Precipitation (mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


# 2008
plot_prec1to3_2008 <- ggplot() +
  geom_sf(aes(fill = prec_1to3_plot), color = 'transparent', data = shapefile_gem_weather_2008) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "Blues", name = 'Precipitation (mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


#2019
plot_prec1to3_2019 <- ggplot() +
  geom_sf(aes(fill = prec_1to3_plot), color = 'transparent', data = shapefile_gem_weather_2019) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "Blues", name = 'Precipitation (mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')

#----------------------#
# Row 3: Weather 4to10 #
#----------------------#

# 1996
plot_prec4to10_1996 <- ggplot() +
  geom_sf(aes(fill = prec_4to10_plot), color = 'transparent', data = shapefile_gem_weather_1996) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "Blues", name = 'Precipitation (mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


# 2008
plot_prec4to10_2008 <- ggplot() +
  geom_sf(aes(fill = prec_4to10_plot), color = 'transparent', data = shapefile_gem_weather_2008) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "Blues", name = 'Precipitation (mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


#2019
plot_prec4to10_2019 <- ggplot() +
  geom_sf(aes(fill = prec_4to10_plot), color = 'transparent', data = shapefile_gem_weather_2019) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "Blues", name = 'Precipitation (mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')



#----------------------------------#
# Create one figure with all plots #
#----------------------------------#

#create row and column labels
row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="prec", angle = 90, size = 5) + theme_void()
row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="prec1to3", angle = 90, size = 5) + theme_void() 
row3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="prec4to10", angle = 90, size = 5) + theme_void()  

col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="1996", size = 5) + theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="2008", size = 5) + theme_void() 
col3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="2019", size = 5) + theme_void() 

#define plot list
plotlist <- list(a = row1, b = row2, c = row3, d = col1, e = col2, f = col3, 
                 g= plot_prec_1996, h=plot_prec_2008, i=plot_prec_2019, 
                 j= plot_prec1to3_1996, k=plot_prec1to3_2008, l=plot_prec1to3_2019,
                 m= plot_prec4to10_1996, n=plot_prec4to10_2008, o=plot_prec4to10_2019)

#define plot arrangement (repeat letters for scaling)
layoutplot <- "
    #ddddddeeeeeeffffff
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    "

#Create and save plot
fig_prec <- patchwork::wrap_plots(plotlist, guides = 'collect', design = layoutplot) & theme(legend.position = 'bottom')

ggsave("Figures/Figure_S3.1.png", fig_prec, device="png", width = 8, height = 8, units = "in", dpi=1200)
ggsave("Figures/Figure_S3.1.eps", fig_prec, device="eps", width = 8, height = 8, units = "in", dpi=1200)


# -------------------------------------------------- #
#### Figure 3: High Growing Degree Days Variables ####
# -------------------------------------------------- #

#-------------------------#
# Row 1: Realized weather #
#-------------------------#

# 1996
plot_gddHigh_1996 <- ggplot() +
  geom_sf(aes(fill = gddHigh_plot), color = 'transparent', data = shapefile_gem_weather_1996) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrRd", name = 'High Growing Degree Days (>30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


# 2008
plot_gddHigh_2008 <- ggplot() +
  geom_sf(aes(fill = gddHigh_plot), color = 'transparent', data = shapefile_gem_weather_2008) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrRd", name = 'High Growing Degree Days (>30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


#2019
plot_gddHigh_2019 <- ggplot() +
  geom_sf(aes(fill = gddHigh_plot), color = 'transparent', data = shapefile_gem_weather_2019) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrRd", name = 'High Growing Degree Days (>30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')

#---------------------#
# Row 2: Weather 1to3 #
#---------------------#

# 1996
plot_gddHigh1to3_1996 <- ggplot() +
  geom_sf(aes(fill = gddHigh_1to3_plot), color = 'transparent', data = shapefile_gem_weather_1996) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrRd", name = 'High Growing Degree Days (>30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


# 2008
plot_gddHigh1to3_2008 <- ggplot() +
  geom_sf(aes(fill = gddHigh_1to3_plot), color = 'transparent', data = shapefile_gem_weather_2008) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrRd", name = 'High Growing Degree Days (>30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


#2019
plot_gddHigh1to3_2019 <- ggplot() +
  geom_sf(aes(fill = gddHigh_1to3_plot), color = 'transparent', data = shapefile_gem_weather_2019) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrRd", name = 'High Growing Degree Days (>30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')

#----------------------#
# Row 3: Weather 4to10 #
#----------------------#

# 1996
plot_gddHigh4to10_1996 <- ggplot() +
  geom_sf(aes(fill = gddHigh_4to10_plot), color = 'transparent', data = shapefile_gem_weather_1996) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrRd", name = 'High Growing Degree Days (>30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


# 2008
plot_gddHigh4to10_2008 <- ggplot() +
  geom_sf(aes(fill = gddHigh_4to10_plot), color = 'transparent', data = shapefile_gem_weather_2008) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrRd", name = 'High Growing Degree Days (>30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


#2019
plot_gddHigh4to10_2019 <- ggplot() +
  geom_sf(aes(fill = gddHigh_4to10_plot), color = 'transparent', data = shapefile_gem_weather_2019) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrRd", name = 'High Growing Degree Days (>30°C)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')



#----------------------------------#
# Create one figure with all plots #
#----------------------------------#

#create row and column labels
row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="gddHigh", angle = 90, size = 5) + theme_void()
row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="gddHigh1to3", angle = 90, size = 5) + theme_void() 
row3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="gddHigh4to10", angle = 90, size = 5) + theme_void()  

col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="1996", size = 5) + theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="2008", size = 5) + theme_void() 
col3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="2019", size = 5) + theme_void() 

#define plot list
plotlist <- list(a = row1, b = row2, c = row3, d = col1, e = col2, f = col3, 
                 g= plot_gddHigh_1996, h=plot_gddHigh_2008, i=plot_gddHigh_2019, 
                 j= plot_gddHigh1to3_1996, k=plot_gddHigh1to3_2008, l=plot_gddHigh1to3_2019,
                 m= plot_gddHigh4to10_1996, n=plot_gddHigh4to10_2008, o=plot_gddHigh4to10_2019)

#define plot arrangement (repeat letters for scaling)
layoutplot <- "
    #ddddddeeeeeeffffff
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    "

#Create and save plot
fig_gddHigh <- patchwork::wrap_plots(plotlist, guides = 'collect', design = layoutplot) & theme(legend.position = 'bottom')

ggsave("Figures/Figure_S3.2.png", fig_gddHigh, device="png", width = 8, height = 8, units = "in", dpi=1200)
ggsave("Figures/Figure_S3.2.eps", fig_gddHigh, device="eps", width = 8, height = 8, units = "in", dpi=1200)

# ---------------------------------- #
#### Figure 4: Dry Days Variables ####
# ---------------------------------- #

#-------------------------#
# Row 1: Realized weather #
#-------------------------#


# 1996
plot_dd_1996 <- ggplot() +
  geom_sf(aes(fill = dd_plot), color = 'transparent', data = shapefile_gem_weather_1996) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrBr", name = 'Dry days (precipitation<10 mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


# 2008
plot_dd_2008 <- ggplot() +
  geom_sf(aes(fill = dd_plot), color = 'transparent', data = shapefile_gem_weather_2008) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrBr", name = 'Dry days (precipitation<10 mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


#2019
plot_dd_2019 <- ggplot() +
  geom_sf(aes(fill = dd_plot), color = 'transparent', data = shapefile_gem_weather_2019) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrBr", name = 'Dry days (precipitation<10 mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')

#---------------------#
# Row 2: Weather 1to3 #
#---------------------#

# 1996
plot_dd1to3_1996 <- ggplot() +
  geom_sf(aes(fill = dd_1to3_plot), color = 'transparent', data = shapefile_gem_weather_1996) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrBr", name = 'Dry days (precipitation<10 mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


# 2008
plot_dd1to3_2008 <- ggplot() +
  geom_sf(aes(fill = dd_1to3_plot), color = 'transparent', data = shapefile_gem_weather_2008) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrBr", name = 'Dry days (precipitation<10 mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


#2019
plot_dd1to3_2019 <- ggplot() +
  geom_sf(aes(fill = dd_1to3_plot), color = 'transparent', data = shapefile_gem_weather_2019) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrBr", name = 'Dry days (precipitation<10 mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')

#----------------------#
# Row 3: Weather 4to10 #
#----------------------#

# 1996
plot_dd4to10_1996 <- ggplot() +
  geom_sf(aes(fill = dd_4to10_plot), color = 'transparent', data = shapefile_gem_weather_1996) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrBr", name = 'Dry days (precipitation<10 mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


# 2008
plot_dd4to10_2008 <- ggplot() +
  geom_sf(aes(fill = dd_4to10_plot), color = 'transparent', data = shapefile_gem_weather_2008) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrBr", name = 'Dry days (precipitation<10 mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')


#2019
plot_dd4to10_2019 <- ggplot() +
  geom_sf(aes(fill = dd_4to10_plot), color = 'transparent', data = shapefile_gem_weather_2019) +
  geom_sf(fill = 'transparent', color = 'black', size = 0.7, data = shapefile_bl) +
  scale_fill_brewer(palette = "YlOrBr", name = 'Dry days (precipitation<10 mm)', drop = FALSE,
                    guide = guide_legend(
                      direction = 'horizontal',
                      title.position = 'top',
                      title.hjust = .5,
                      label.hjust = .3,
                      label.position = 'right',
                      keywidth = 1,
                      keyheight = .3,
                      nrow =1
                    )) +
  labs(title = expression('')) +
  theme_void() +
  theme(plot.title = element_text(size= 12, face = 'bold', hjust = 0.5), 
        plot.subtitle = element_text(size= 12, face = 'bold', hjust = 0.5),
        legend.title = element_text (size = 12),
        legend.text = element_text (size = 12),
        legend.position = 'bottom')



#----------------------------------#
# Create one figure with all plots #
#----------------------------------#

#create row and column labels
row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="dd", angle = 90, size = 5) + theme_void()
row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="dd1to3", angle = 90, size = 5) + theme_void() 
row3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="dd4to10", angle = 90, size = 5) + theme_void()  

col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="1996", size = 5) + theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="2008", size = 5) + theme_void() 
col3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="2019", size = 5) + theme_void() 

#define plot list
plotlist <- list(a = row1, b = row2, c = row3, d = col1, e = col2, f = col3, 
                 g= plot_dd_1996, h=plot_dd_2008, i=plot_dd_2019, 
                 j= plot_dd1to3_1996, k=plot_dd1to3_2008, l=plot_dd1to3_2019,
                 m= plot_dd4to10_1996, n=plot_dd4to10_2008, o=plot_dd4to10_2019)

#define plot arrangement (repeat letters for scaling)
layoutplot <- "
    #ddddddeeeeeeffffff
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    agggggghhhhhhiiiiii
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    bjjjjjjkkkkkkllllll
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    cmmmmmmnnnnnnoooooo
    "

#Create and save plot
fig_dd <- patchwork::wrap_plots(plotlist, guides = 'collect', design = layoutplot) & theme(legend.position = 'bottom')

ggsave("Figures/Figure_S3.3.png", fig_dd, device="png", width = 8, height = 8, units = "in", dpi=1200)
ggsave("Figures/Figure_S3.3.eps", fig_dd, device="eps", width = 8, height = 8, units = "in", dpi=1200)

