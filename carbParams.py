""" 
Parameters to run "carbMod" model
"""

import os
import glob
import numpy as np

################################################################################

#-----------------------------     PARAMETERS     -----------------------------#

def setOutputPath():
    """Create necessary folders to run the script and return the folder in which
    the results will be stored"""
    results = os.path.join('carb_results')
    txt_files = os.path.join('txt_files')
    arr_cellBasedStocks = os.path.join(results, 'arr_cellBased')
    arr_totalStocks = os.path.join(results, 'arr_totalStocks')
    maps = os.path.join(results, 'maps')
    statistical_test = os.path.join('statistical_test')
    paths = [results, txt_files, arr_cellBasedStocks, 
             arr_totalStocks, maps, statistical_test]
    for path in paths:
        if not os.path.exists(path):
            os.makedirs(path)
    return results


def getScenariosPaths(path, wildcard):
    """Return a dictionary of all the scenarios (keys) with the paths to access 
    LUC maps (values). The dictio includes the initial state (sc0 = year 2012)"""
    scenariosList = sorted(glob.glob(os.path.join(path, wildcard)))
    scenariosDict = {
        int((scPath.split('sc'))[-1]): scPath for scPath in scenariosList}
    return scenariosDict


def getMonteCarloSamples():
    """Return number of Monte Carlo samples. If nrMCsamples == 1, the script
    will run deterministically. If nrMCsamples > 1, then is stochastic."""
    nrMCsamples = 10000
    return nrMCsamples


def configModelComponents():
    """Return the model components setup for sensitivity analysis: it returns 
    which component will run deterministically (0) or stochastically (1). """
    SOC_component = 1
    BC_component = 1
    LUC_component = 1
    lst = [SOC_component, BC_component, LUC_component]
    if lst == [1, 0, 0]:
        runType = 'stcSOC'
        description = 'SOC stochastic (for Sensitivity analysis)'
    elif lst == [0, 1, 0]:
        runType = 'stcBC'
        description = 'BC stochastic (for Sensitivity analysis)'
    elif lst == [0, 0, 1]:
        runType = 'stcLUC'
        description = 'LUC stochastic (for Sensitivity analysis)'
    elif lst == [1, 1, 1]:
        runType = 'stcAll'
        description = 'Full stochastic'
    elif lst == [0, 0, 0]:
        if getMonteCarloSamples() != 1:
            raise Exception("For a full deterministic run based on the\
            uncertainty  files used for 'full stochastic runs', you have to set\
            carbParams.getMonteCarloSamples() = 1") 
        runType = 'detAll'
        description = 'Full deterministic'        
    else:
        #runType = 'stcAll'
        #description = 'tst'
        raise Exception('Two components set as stochastic and one set as\
        deterministic. For Sensitivity Analysis, set just one component to run\
        stochastically.')
    print('Model components setup: SOC: {}, BC: {}, LUC: {}\t ==>\t{}\n'.format(
        SOC_component, BC_component, LUC_component, description))
    return runType, SOC_component, BC_component, LUC_component

def getInitialLandUseMap():
    initLUmap = os.path.join('PLUC', 'projIniMaps', 'ini1.map')
    return initLUmap

def referenceRaster():
    """Set the reference raster in which its geoparameters will be used to 
    create/save new rasters/numpy arrays."""
    refRaster = getInitialLandUseMap()
    return refRaster

def getMappingFromFile(filename, col1, col2):
    """Return mapping of new vs original particle number. Parameters: col1 == 
    column with the ID of LU maps; col2 = column with the weight of LU maps"""
    # read from csv file that is result of particle filter
    path = os.path.join(filename)
    # Original MCruns + multiply factor
    MCruns_CarbonMod = 10000
    MCruns_Orig = 5000
    multiplyFactor = MCruns_CarbonMod / MCruns_Orig
    # Particle filter mapping
    data = np.genfromtxt(path, dtype=int, delimiter='\t', names=True)
    # make dictionary of LUmap_code : LUmap_weight
    mapping = dict(zip(data[col1], data[col2] * multiplyFactor))
    return mapping


def getStochasticLUmaps(scenario_path, pfMapping, wildcard):
    """Assuming LU maps are the same for each SC, it returns a dictionary with 
    the total MC samples per PLUC PartFilterDict => e.g: if 03 PLUC runs (key), 
    then the Nr of runs == dict.values (assuming the model will run 10000x)"""
    LUmapsList = glob.glob(os.path.join(scenario_path, wildcard))
    LUmapsDict = {
        int((v.split('2030_'))[-1].split(".")[0]): v for v in LUmapsList}
    mappingSet = set(pfMapping)
    LUmapsDictSet = set(LUmapsDict)
    LUmapsToRun = {}
    nrRuns_acc = 0
    nrRuns_mcRange = 0
    # Preparing mapping to run based on the getMonteCarloSamples(), 
    ## In case getMonteCarloSamples() is less than the total runs of the first,
    ## do the following..
    if getMonteCarloSamples() < pfMapping[1]:
        LUmapsToRun[1] = [LUmapsDict[1], getMonteCarloSamples()]
    # Now creating a dictionary with the maps the will be used in the model and
    # related total runs
    else:
        for nrLU in sorted(mappingSet.intersection(LUmapsDictSet)):
            nrRuns = pfMapping[nrLU]
            LUmap = LUmapsDict[nrLU]
            #print nrLU, nrRuns, LUmap
            nrRuns_acc += nrRuns
            if nrRuns_acc in range(getMonteCarloSamples()):
                nrRuns_mcRange += nrRuns
                LUmapsToRun[nrLU] = [LUmap, nrRuns]
            # In case the getMonteCarloSamples() does not have the same value
            # than a particle filter. Next line adds a last map and subtracts
            # the Nr of accumulated Nr of runs
            if nrRuns_mcRange == nrRuns_acc:
                LUmapsToRun[nrLU + 1] = [LUmapsDict[nrLU + 1],
                                         getMonteCarloSamples() - nrRuns_acc]
    return LUmapsToRun


def getDeterministicLUmap(scenario_path, pfMapping):
    highestWeight = max(pfMapping.keys(), key=(lambda k: pfMapping[k]))
    LUmap = os.path.join(scenario_path, '2030_{}.map'.format(highestWeight))
    # print "Particle Filter: Highest weight between LU maps == code {} (weight\
    # == {}, considering 10000 runs)".format(highestWeight, pfMapping[highestWeight])
    return LUmap


def getConversionUnit():
    """Return conversion unit to get the the area of a cell. The initial results
    are given in tonne C/ha. Hence, the conversion unit is set to obtain the 
    carbon stocks in total hectares"""
    toHectares = 2500
    return toHectares


def getEmissionFactorConversion():
    """Return the factor used to convert Carbon stocks in Emission Factor (EF)."""
    EF = 44.0 / 12.0  # necessary to convert Carbon Stocks to Emission Factors
    return EF


def getScenariosType():
    initial_state = 0
    scenarios_no_additional_Eth = [1, 4, 7, 10, 13, 16]
    scenarios_additional_Eth = [3, 6, 9, 12, 15, 18]
    return initial_state, scenarios_no_additional_Eth, scenarios_additional_Eth


def figureColors():
    """1st => SOC; 2nd => BC"""
    colors = ['#d8b365', '#5ab4ac']
    return colors


def getScenariosNames():
    """Return list of scenarios  with name."""
    scenariosNames = ['No measures', 'High prod.', '2nd gen. SC',
                      '2nd gen. EU', 'Cons. policies', 'All measures']
    return scenariosNames


def setNrLUtypes():
    """Set the number of land use types for the case study. Zero is not LU type"""
    nrLUtypes = 11
    return nrLUtypes

def LUtypesToMask():
    """Set the land use types to use for masking numpy arrays including "0" 
    which represents NoData cells"""
    LuTypes = [0, 1, 2, 10]
    return LuTypes
