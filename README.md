
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

- Land use component: is the model component that accounts for uncertainty or not in the land use maps. If stochastic, the stochastic maps of the PLUC model are used, based on a particle filter. If deterministic, then a unique map from the PLUC model is used (represented by the land use with the highest weight in the particle filter).

4. Inputs

Maps:

climate.map (deterministic) - represent the climate regions that are in line with the IPCC guidelines (2006). The data is derived from Hijmans et al. (2005) and NASA/NGA (2001).

soil.map (deterministic) - based on the IPCC default soil classes derived from the Harmonized World Soil Data Base (Batjes, 2010)
land use maps (stochastic) - derived from stochastic runs of PLUC (Verstegen, 2012)

All maps has the same extent, cell size (2500 x 2500), number of rows (854) and columns (XXX). The used spatial reference is the WGS 1984 geographic coordinate system, with the Albers equal-area conic projection. To run the model with no error/bias in the computation, those properties should be shared between all maps.  


