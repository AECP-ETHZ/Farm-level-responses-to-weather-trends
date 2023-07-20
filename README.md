# Farm-level-responses-to-weather-trends
This code replicates the figures and tables presented in the article “Farm-level responses to weather trends: A structural model”, published in the American Journal of Agricultural Economics: https://onlinelibrary.wiley.com/doi/10.1111/ajae.12421

The files are as follows:

-	Programs in 01a_OperWeather.txt operationalize the weather variables (in SAS). 
-	Programs in 01b_PrepWeather.R and 01c_WeatherMaps.R prepare and map the weather data used for the analysis. 
-	Programs in 02_PrepPrices.R prepare the outputs and input price indices.
-	Programs in 03_PrepFADN.R prepare the farm-level data.
-	Programs in 04_PrepSampleVars construct the sample and the variables used in the analysis.
-	Programs in 05a_EstSim_main estimates the main model and simulates farm-level responses to a drought shock; Programs in 05b_Bootstrap_main run the bootstrap and create the results figure and tables for the main model. 
-	Programs in 06a & 06b – 11a & 11b and 13a & 13b do the same for alternative models (heterogeneity with respect to farm size, different lag structures in the weather variables, nonlinearity in weather, curvature restriction, and a model without accounting for non-random crop selection).
-	Programs in 12_Est_OutOfSample.R estimate the model without the final year of the data and assesses out-of-sample prediction.


More detailed instructions for replication, including data availability statements, can be found in the file "README.pdf".
