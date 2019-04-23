
1. Introduction

A carbon stock changes model derived from land use changes was created to quantify uncertainty.

The model is developed to quantify uncertainty on the projections of LUC-related GHG emissions in Brazil towards 2030, given an expected increase in the global biofuel demand and distinct scenarios of LUC mitigation measures. The GHG emissions are evaluated based on carbon stock changes in Soil Organic Carbon (SOC) and Biomass Carbon (BC) between the initial (2012) and the final (2030) state of the land use system.

The model is set to run for five distinct mitigation measures scenarios, including the reference scenario: high agricultural productivity; change to the second generation of ethanol (sugar cane); change to the second generation of ethanol (eucalyptus); stricting conservation policies; and all measures together. Each scenario is evaluated twice: with and withouth a global increase in ethanol demand. This allows to evaluate the impact of additional ethanol demand in green house gas emissions derived from land use changes.

2. Specification

The model was built in Python 2.7 and uses the following libraries: gdal, osr, glob, os, numpy, pandas, scipy.stats and matplotlib. Also, it makes use of the PCRaster version 4.1.0 to run the land use change model (PLUC).

3. Description

The model main file is carbMod.py. Also, it includes the file carbParams.py that represents the parameters used in the main file. Therefore, to run the script, one should setup the former firstly. 

The uncertainty introduced in the model has two origins: the land use maps derived from the stochastic runs of the PLUC model (which includes a particle filter); and the input data based on Tier 1 IPCC values (IPCC Guidelines for National Greenhouse Gas Inventories, 2006) that is used to compute the carbon stocks in land, given the soil and climate characteristics of a given land use type (see climate.map and soil.map)

To quantify uncertainty, the model makes use of the Monte Carlo Simulation technique. Hence, the process of computing carbon stock changes is repeated in a loop based on the amount of Monte Carlo simulations defined by the user (hereto we use 10.000). For each sample, random values of SOC and Biomass Carbon Stocks are created by the model (based on IPCC input data) and then allocated in the stochastic land use maps, given the climate and soil maps.  

The model is built to run stochastically (i.e., accounting for uncertainty) but it can also run deterministically (no uncertanty). Also, to proceed with sensitivity analysis, the model is divided in three components in which each component can be 'isolated' to run either stochastically or deterministic. The model components are:

- SOC component: is the model component that computes SOC stocks. If stochastic, the model generate random values in each iteration . If deterministic, specific values are used in all iteration (no generation of random values).

- BC component: is the model component that computes BC stocks. It has the same behavior of "SOC" component: random values if stochastic and specific values if deterministic.

- LUC component (Land Use Change): is the model component that accounts for uncertainty or not in the land use maps. If stochastic, the stochastic maps of the PLUC model are used, based on a particle filter. If deterministic, then a unique map from the PLUC model is used (represented by the land use with the highest weight in the particle filter).

There are two types of computation of carbon stock changes and related GHG emissions: Overall and Cell-Based. Overall Carbon Stocks/Emissions are related to the total stocks/emissions in the whole study area (i.e., the stocks and emissions in each cell of the output maps are sum). This is used to compute the outputs that are not related to maps: boxplot, plot of sensitivity analysis, emissions w/ uncertainty quantification and txt files filles with statistics. Cell-based carbon stocks are related to the total stocks/emissions per cell. This is used to produce the maps of cell-based mean and standard deviation of carbon stocks/emissions. 

4. Inputs

4.1 Maps

- climate.map (deterministic) - represent the climate regions that are in line with the IPCC guidelines (2006). The data is derived from Hijmans et al. (2005) and NASA/NGA (2001).

- soil.map (deterministic) - based on the IPCC default soil classes derived from the Harmonized World Soil Data Base (Batjes, 2010)

- ini1.map (deterministic) - initial land use map representing the initial state of the system. It is derived from PLUC.

- Land use maps (stochastic) - derived from stochastic runs of PLUC (Verstegen, 2012)

All maps have the same extent, cell size (2500 x 2500), number of rows and columns (854, 885). The spatial reference is the WGS 1984 geographic coordinate system, with the Albers equal-area conic projection. To run the model with no error/bias in the computation, those properties should be shared between all maps.  

4.2 Text/csv files

- bcs_uncertainty.txt - it has codes with statistics of BC stocks used by the model to generate random values for BCS.
- socF_uncertainty.txt - it has codes with statistics of SOC factor (land use factor * land management factor * land inputs factor) used by the model to generate random values for SOCf  
- socR_uncertainty.txt - It has codes with statistics of SOC stock reference values  used by the model to generate random values for SOCr.  

- bcs_detSensAnalysis.txt - codes w/ deterministic values for BCS used to run sensitivity analysis.
- bcs_deterministic.txt - codes w/ IPCC's deterministic values for BCS used to run the model 'full-deterministic'. 
- socF_deterministic.txt - codes w/ deterministic values for SOCf to run the model either full-deterministic (IPCC) or for sensitivity analysis
- socR_deterministic.txt - codes w/ deterministic values for SOCr to run the model either full-deterministic (IPCC) or for sensitivity analysis

- climate_reclass.txt - file used to reclassify the original climate.map in order to produce the correct codes for the "climsoil" raster (climate + soil combination)
- magnet_output.txt - it has the total ethanol production up to 2030 for each scenario  
- legend_scenarios.txt - description of each scenario 
- particle_mapping.csv - particle mapping derived from PLUC. It is used to set which land use map should run for each iteration based on the weight of the particle filter.

5. Additional files:

- tables_stochastic.xlsx - file of which the uncertainty .txt files are created
- tables_deterministic.xlsx - file of which the deterministic .txt files are created
- input_codes.xlsx - description of each code given in the txt files

6. Python files

To run the model, the carbParams.py must be configured first. Then, the carbMod.py must be configured in the "SETUP" section.

6.1 carbParams.py

Configures the output paths, nr of monte carlo samples, the path for each scenario with the land use change maps, the configuration of model components, a reference raster to produce map outputs, the initial land use map, the particle filter mapping, etc. The 'main' functions of this script are:

- getMontecarloSamples(), where the user define the amount of samples for the Monte Carlo simulation

- configModelComponents(), where the user can define which component should run either stochastically or deterministically (mainly for sensitivity analysis). Important: if you want a full deterministic run based on the uncertainty files used for stochastic runs, you can't define here all components = 0. Therefore, you have to set the getMonteCarloSamples() = 1

- getMappingFromFile(), where the particle filter derived from PLUC is used

- getStochasticLUmaps(), where each LUC map is assigned to run a certain amount of times, based on the particle filter.

- getDeterministicLUmaps(), where the LUC map with the highest weight of the particle filter is set to be used when running in deterministic mode.

- LUtypesToMask(), where the codes representing LU types that are assumed to have no carbon stocks are defined.

6.2 carbMod.py

The carbMod is divided in Parts (explicitly shown in the script). The short description below shows how the model works when dealing with a  single scenario but keep in mind that this is repeated for each scenario. The "setup" part the carbMod Parts and it is used to basically turn on/turn off them, as follow:

> 1. Generating random values: If turned on, the random files are created based on a random seed and saved in folder "txt_files".

Inputs: bcs_uncertainty.txt, socF_uncertainty.txt, socR_uncertainty.txt
Outputs: txt files for bcs, socF and socR, according to the nr of Monte Carlo samples.

> 2. Processing overall carbon stocks: If turned on, computes the overall carbon stocks in the study area and save the data that is produced. If turned off, then the data is not produced but is loaded (this allows to not run the all the model again if necessary).

New inputs: climate.map, soil.map, init1.map, land use change map(s). If SOC_comp and/or BC_comp = 1 i.e. stochastic: it uses the txt files produced in part 1. If SOC_comp and/or BC_comp = 0 i.e. deterministic: it uses the deterministic text files (see item 4 - Inputs). 

Outputs: if the variable saveArrays4ov_stock = 1, saves a numpy array (npy) with the Overall carbon stocks for SOC/BC/SOC+BC(TC), for each Monte Carlo sample (used in Parts 3, 4 and 5). If the variable saveArrays4cb_stock == 1, saves for each Monte Carlo sample a compressed numpy array (npz) representing the cell-based carbon stocks for SOC and BC. This is afterwards used in Part 6.

In this Part, the difference between carbon stocks in scenarios with ethanol (scEth) and with additional ethanol (scAddEth) are also computed. Also, the function IPCCdet_getDiffStockPerPairScenarios() is run, which represents the processed descripted in this Part but considering a full deterministic run based on IPCC default values.

> 3. Getting overall statistics & saving:
Saves statistics of absolute carbon stocks and with regards to the difference in stocks between scenarios of scEth and scAddEth. It also prints the mean and uncertainty range (in percentage).

New inputs: none
Output: csv files with statistics

> 4. Computing overall emissions from ethanol production:
Converts the difference in carbon stocks between pair scenarios (scEth Vs scAddEth) in the final unit of measurement (gram CO2-eq Mj-1 EtOH)

New inputs: magnet_output.txt
Outputs: csv file with the overall emissions from ethanol production (that are also printed).

> 5. Outputs: represents the obtainment of the 'main' outputs: boxplot, plot of sensitivity analysis and data for statistical test. Important: for sensitivity analysis, the model must run four times: for each time,
the modelRunType variable must be set to 'stcSOC'/'stcBC'/'stcLUC'/'stcAll' (based on the carbParams.configModelComponents() function. The sensitivity analysis will only work if the variable  saveArrays4ov_stock = 1 before proceeding with each of the four runs

New inputs: none.
Outputs: already mentioned. 

> 6. Processing cell-based carbon stocks: Computes the cell-based mean & standard deviation regarding to absolute stocks in the study area. To do that, it loads the numpy compressed files that were saved in Part 2. 

Also it computes the cell-based mean & standard deviation difference between paired scenarios (scEth vc scAddEth) regarding to GHG emissions in the study area (this is not yet working in the script).

New inputs: none.
Output: arrays and maps of cell-based mean & standard deviation for absolute stocks and GHG emissions due to increse in ethanol production.