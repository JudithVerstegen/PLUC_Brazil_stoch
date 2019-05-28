""" 
Carbon emissions model derived from land use changes in Brazil
"""

import time
print("Starting... {}".format(time.asctime(time.localtime(time.time()))))

import gdal
import osr
import glob
import os
import numpy as np
import pandas as pd
import carbParams
import scipy.stats as ss
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from os.path import join as opj
from six import iteritems, itervalues


################################################################################

#--------------------------------     SETUP     -------------------------------#
np.set_printoptions(suppress=True)

# Creating outputs' path (all results are saved in resultsDir variable)
resultsDir = carbParams.setOutputPath()

# Setting nr of Monte Carlo Samples 
nrMCsamples = carbParams.getMonteCarloSamples()

# Setting model run type and components (used for turn on/turn off the stochastic/
# deterministic runs of each model components)
modelRunType, SOC_comp, BC_comp, LUC_comp = carbParams.configModelComponents()

# Model configuration (you must turn on/turn off, depending of which part you 
# want to run 
generateRandomValues = 0
getOverallCarbonStock = 1
getCellBasedCarbonStock = 0

# Defining if you want to save arrays. For overall carbon stocks: used to save the 
# arrays w/ the sum of CStocks per MCr. If 1, the arrays will be saved regardless 
# the 'modelRunType' setup. For cell based stocks: the model is set to just save
# if 'modelRunType' = 'stcAll' (full stochastic), to use in cell-based computations)
saveArrays4ov_stock = 0
saveArrays4cb_stock = 0

# Setup to plot/show/save figures...
# OBS: For sensitivity analysis, you must run the model four times: For each time,
# the modelRunType variable must be set to 'stcSOC'/'stcBC'/'stcLUC'/'stcAll'. 
# (see carbParams.configModelComponents() for better understanding). IMPORTANT: 
# the sensitivity analysis is only possible if you set saveArrays4ov_stock to 1 
# in each of the runs!
plotSensAnalysis = 0
plotBoxplot = 1
showPlots = 1
savePlots = 1

# Dictionary corresponding to the path for scenarios with LU maps,
scenariosDict = carbParams.getScenariosPaths(opj(os.path.join('PLUC', \
                                                'results_stoch')), 'sc*')
# Adding the deterministic initial LU map to scenariosDict (k = 0, v = LU path)
scenariosDict[0] = carbParams.getInitialLandUseMap()

print ('Nr of scenarios (including initial): {}.\nMC simulation: {} samples per\
 scenario.\n'.format(len(scenariosDict), nrMCsamples))

# Functions related to random values
def generateRandomSOCref(inpTxt, outPath):
    """For each Monte Carlo sample, create a text file filled with random values
    for 'SOC reference' w/ a factor of 1000. Arguments: inpTxt = Input file with
    uncertainty range; outPath = path to save the text file."""
    mcRange = np.arange(carbParams.getMonteCarloSamples())
    for i in mcRange:
        sample = i + 1
        outTxt = open(opj(outPath, 'socR_{}.txt'.format(str(sample))), 'w')
        with open(inpTxt, 'r') as inp:
            inpHeader = inp.readline()
            outHeader = outTxt.write('id\tvalue\n')
            for row in inp:
                rowValue = row.split()
                code_id = rowValue[0]
                SOCr_mean = float(rowValue[1])
                SOCr_std = float(rowValue[2])
                SOCr = np.random.normal(SOCr_mean, SOCr_std)
                while SOCr <= 0:
                    print (
                        'Neg/zero value computed for SOCr. Recomputing...')
                    SOCr = np.random.normal(SOCr_mean, SOCr_std)
                outTxt.write('{}\t{}\n'.format(code_id, str(int(SOCr * 1000))))
        outTxt.close()


def generateRandomSOCfactor(inpTxt, outPath):
    """For each Monte Carlo sample, create a text file filled with random values
    for 'SOC factors' w/ a factor of 1000. Arguments: inpTxt = Input file with 
    uncertainty range; outPath = path to save the text file."""
    mcRange = np.arange(carbParams.getMonteCarloSamples())
    for i in mcRange:
        sample = i + 1
        outTxt = open(opj(outPath, 'socF_{}.txt'.format(str(sample))), 'w')
        with open(inpTxt, 'r') as inp:
            inpHeader = inp.readline()
            outHeader = outTxt.write('id\tvalue\n')
            for row in inp:
                rowValue = row.split()
                code_id = rowValue[0]
                fLU_mean, fLU_std = float(rowValue[1]), float(rowValue[2])
                fMG_mean, fMG_std = float(rowValue[5]), float(rowValue[6])
                fI_mean, fI_std = float(rowValue[9]), float(rowValue[10])
                                # Setting SOCf = 0 for LU types w/out carbon stocks
                if (fLU_mean == 0 and fMG_mean == 0 and fI_mean == 0):
                    SOCf = 0
                # Setting SOCf = 1 for LU types w/out IPCC uncertainty range:
                elif (fLU_mean == 1 and fMG_mean == 1 and fI_mean == 1):
                    SOCf = 1
                # Computing SOCf for LU types w/ carbon stocks
                else:
                    # Getting random values for SOC factors (LU = land use;
                    # MG = land management; I = land inputs)
                    fLU = np.random.normal(fLU_mean, fLU_std)
                    fMG = np.random.normal(fMG_mean, fMG_std)
                    fI = np.random.normal(fI_mean, fI_std)
                    while (fLU <= 0 or fMG <= 0 or fI <= 0):
                        fLU = np.random.normal(fLU_mean, fLU_std)
                        fMG = np.random.normal(fMG_mean, fMG_std)
                        fI = np.random.normal(fI_mean, fI_std)
                        print (
                            'Neg.value computed for SOC factor. Recomputing...')
                    SOCf = fLU * fMG * fI
                outTxt.write('{}\t{}\n'.format(code_id, str(int(SOCf * 1000))))
        outTxt.close()


def generateRandomBCstock(inpTxt, outPath):
    """For each Monte Carlo sample, create a text file filled with random values
    for 'biomass carbon stocks' (bcs) w/ a factor of 1000. Arguments: inpTxt = 
    Input file with uncertainty range;outPath = path to save the text file."""
    mcRange = np.arange(carbParams.getMonteCarloSamples())
    for i in mcRange:
        sample = i + 1
        outTxt = open(opj(outPath, 'bcs_{}.txt'.format(str(sample))), 'w')
        with open(inpTxt, 'r') as inp:
            inpHeader = inp.readline()
            outHeader = outTxt.write('id\tvalue\n')
            for row in inp:
                rowValue = row.split()
                code_id = rowValue[0]
                agb_mean, agb_std = float(rowValue[1]), float(rowValue[2])
                r2s_mean, r2s_std = float(rowValue[5]), float(rowValue[6])
                # Set zero as min because the r2s PDFs have neg. values
                r2s_min, r2s_max = 0, float(rowValue[8])
                cf_mean, cf_std = float(rowValue[9]), float(rowValue[10])
                BCS_mean, BCS_std = float(rowValue[13]), float(rowValue[14])
                # Getting random values for above-ground biomass (AGB)
                agb = np.random.normal(agb_mean, agb_std)
                while agb < 0:
                    print ('Neg.value computed for AGB. Recomputing...')
                    agb = np.random.normal(agb_mean, agb_std)
                # Getting random values for root-to-shoot (r2s)
                r2s = np.random.normal(r2s_mean, r2s_std)
                while r2s < 0:
                    print (
                        'Neg.value computed for R2S. Recomputing by truncation ...')
                    if r2s_std != 0:  # scipy.stats does not allow std = 0,
                        # therefore I use truncation:
                        a, b = (r2s_min - r2s_mean) / r2s_std, (r2s_max - r2s_mean) / r2s_std  
                        r2s = ss.truncnorm.rvs(a, b, r2s_mean, r2s_std)
                    else:
                        r2s = np.random.normal(r2s_mean, r2s_std)
                # Getting random values for carbon fraction (CF)
                cf = np.random.normal(cf_mean, cf_std)
                while cf < 0:
                    print ('Neg. value computed for CF. Recomputing...')
                    cf = np.random.normal(cf_mean, cf_std)
                # Computing below-ground biomass (BGB)
                bgb = agb * r2s
                # Computing biomass carbon stocks  
                BCS = (agb + bgb) * cf
                # Now, if clause used to compute BCS for 'crops' and 'abandoned 
                # land' (IPCC doesn't give AGB, R2S and CF for those LU types)
                if BCS == 0:
                    BCS = (np.random.normal(BCS_mean, BCS_std))
                    while BCS < 0:
                        print (
                            'Neg. value identified for BCS. Getting positive value ...')
                        BCS = np.random.normal(BCS_mean, BCS_std)
                outTxt.write('{}\t{}\n'.format(code_id, str(int(BCS * 1000))))
        outTxt.close()

# Functions used to get overall total carbon stocks (not cell-based)
def checkNegRandomValues(txtFilesFolder):
    """Check the occurrence of negative values in all the txt files generated by
    the random values' functions"""
    files2check = sorted(glob.iglob(opj(txtFilesFolder, '*')))
    for txtFile in files2check:
        df = pd.read_csv(txtFile, sep='\t', header=0)
        if df[df.value < 0].empty == False:
            print (txtFile.split('txt_files/')[1], '\n', df[df.value < 0])
            raise Exception('Negative value found in file {}'.format(
                txtFile.split('txt_files/')[1]))


def rasterToArray(inpRaster):
    """Return a numpy array based on a raster file. The 'No data' of the raster 
    are converted to zero"""
    inpRaster = gdal.Open(inpRaster)
    unqBand = inpRaster.GetRasterBand(1)
    noData = unqBand.GetNoDataValue()
    npArr = inpRaster.ReadAsArray()
    npArr[npArr == noData] = 0  # Converting noData values to zero
    if npArr[npArr < 0].any() == True:  # Making sure no neg value is in array
        npArr[npArr < 0] = 0
    return npArr


def reclassArray(npRaster, txtFile):
    """ Return a reclassified array based on a text file representing codes / 
    random values ​for SOCF, SOCr or BCS. The codes given by the text files must 
    encompass all the codes given by the array. Also, it removes the factor of 
    1000 of the text files. Arguments: npRaster = numpy array representing codes
    of a raster file; txtFile = text file with codes and random values used for
    reclassification."""
    # 'Code-value' reference matrix: convert the .txt into array
    cvrm = np.loadtxt(txtFile, dtype=np.int32, skiprows=1)  
    code, value = cvrm.T  # Transposed array
    # Lookup table with all the codes of npRaster. '+1' encompasses all codes
    lut = np.arange(npRaster.max() + 1)
    # Assign a value for the lookup table based on the associated code
    lut[code] = value
    # Reclassify npRaster filled w/ codes according to lookup table's values
    reclassifArr = (lut[npRaster]).astype(np.float32) / 1000
    return reclassifArr

def getClimSoil(climate_raster, soil_raster):
    """ Return a numpy array representing the combination of climate and soil 
    raster (which is afterwards used to obtain soil organic carbon stocks. Args: 
    climate_raster = climate raster; soil_raster = soil raster"""
    # 1st step: reclassifying climate_raster array
    txtFile = 'climate_reclass.txt'
    # 'Code-value' reference matrix: convert the .txt into array
    cvrm = np.loadtxt(txtFile, dtype=np.int32, skiprows=1)  
    code, value = cvrm.T  # Transposed array
    # Lookup table with all the codes of npRaster. '+1' encompasses all codes
    lut = np.arange(climate_raster.max() + 1)
    # Assign a value for the lookup table based on the associated code
    lut[code] = value
    # Reclassify npRaster filled w/ codes according to lookup table's values
    climate_reclassif = (lut[climate_raster]).astype(np.int16)
    # Creating mask where soil == 0 (meaning that the cell represents water)
    mask = np.ma.masked_where(soil_raster == 0, soil_raster)
    # 2nd step: getting climate + soil combination
    climsoil = climate_reclassif + soil_raster
    # Replacing the output cells with the mask 
    climsoil[mask.mask] = 0 
    return climsoil


def maskNoCarbonCells(LUmap_array):
    """LUmask: masked array of a given land use map in which the cells that are 
    assumed to have no carbon (i.e. lu types with no carbon stocks) are masked. 
    Argument: LUmap_array = array of a given land use map"""
    noCarbonValues = carbParams.LUtypesToMask()
    mask = np.logical_or.reduce(
        [LUmap_array == value for value in noCarbonValues])
    maskedArray = np.ma.masked_where(mask == True, LUmap_array)
    return maskedArray


def getLUmap(sc, scenario_path, LUmapStoch=""):
    """Returns the Land Use map(s) to be used in the carbon stock computation. 
    The output depend if the 'LUC_comp' variable is set to run deterministically
    or stochastically. Args: sc = scenario; scenario_path = scenario folder; 
    LUmapStoch = only used when the LUC component is set to run stochastically."""
    if sc == 0:
        LUmap = carbParams.getInitialLandUseMap()
    if sc > 0:
        if nrMCsamples == 1:
            LUmap = carbParams.getDeterministicLUmap(
                    scenario_path, particleFilter_mapping)
        if nrMCsamples > 1:
            if LUC_comp == 1:  # stochastic
                LUmap = LUmapStoch
            if LUC_comp == 0:  # deterministic
                LUmap = carbParams.getDeterministicLUmap(
                    scenario_path, particleFilter_mapping)
    return LUmap

def getLUmapMCrange(nrLUmap, mcRuns_LUmap, mcRun_accumulated):
    """Return the Monte Carlo samples Range for each LUC map. It is used when 
    LUC_comp == 1 (stochastic)). Arguments: LUmapsDict = dictionary with the LUC
    maps to be used and related MC runs based on Particle Filter Mapping 
    (see carbParams.getStochasticLUmap function); nrLUmap = code of the LU map; 
    mcRuns_LUmap: total mcRuns of a given land use map;  mcRun_accumulated = 
    the current accumulated MCrun of the for loop (since this function is run 
    in a for loop)"""
    # Getting the exact Monte Carlo samples range in which a given land use map 
    # will be used 
    mcRunToStart = mcRun_accumulated + 1
    mcRun_accumulated += mcRuns_LUmap
    mcRunToStop = mcRun_accumulated + 1
    mcRange = np.arange(mcRunToStart, mcRunToStop)
    return mcRange

def getTotalSOCstock(climLu, climSoil, LUmask, SOC_comp, sample=None):
    """Returns the total SOC stock for a given land use map (soc == array with 
    total stock per cell, tSOC == sum of array). The outputs depend whether the 
    'SOC_comp' variable is set to run deterministically or stochastically. 
    Arguments: climLU = array of climate + lu maps combination (for SOC factor);
    climSoil = array of climate + soil maps combination (for SOC reference); 
    LUmask = it is the output of the maskNoCarbonCells() function; SOC_comp = 
    setup of the SOC component (see carbParams.configModelComponents)"""
    # Setting the path of txt files (depending if stochastic or not)
    if SOC_comp == 1:  # stochastic
        txtFile_socF = opj(txtFilesPath, 'socF_{}.txt'.format(str(sample)))
        txtFile_socR = opj(txtFilesPath, 'socR_{}.txt'.format(str(sample)))
    if SOC_comp == 0 or nrMCsamples == 1: # deterministic for sensitivity analysis
        txtFile_socF = opj('socF_deterministic.txt')
        txtFile_socR = opj('socR_deterministic.txt')
    if SOC_comp == 'IPCC':  # IPCC deterministic - this is only used inside the 
        #IPCCdet_getOverallCarbonStock() function (now files are == sensAnalysis
        # but we might differentiate afterwards)
        txtFile_socF = opj('socF_deterministic.txt')
        txtFile_socR = opj('socR_deterministic.txt')
    # Allocating random values in arrays (SOC factors, then SOC reference values)
    socF = reclassArray(climLu, txtFile_socF)
    socR = reclassArray(climSoil, txtFile_socR)
    soc = socF * socR
    # Ensuring that cells with no stocks remain with no stocks
    soc[LUmask.mask] = 0
    # Summing total stocks in the array
    tSOC = np.sum(soc)
    return soc, tSOC


def getTotalBCstock(climLu, LUmask, BC_comp, sample=None):
    """Returns the total biomass carbon stocks for a given land use map 
    (bc == array with total stock per cell, tBC == sum of array). 
    The outputs depend whether the 'BC_opjcomp' variable is set to run 
    deterministically or stochastically. Arguments: climLU = array of climate + 
    lu maps combination (for SOC factor); LUmask: it is the output of the 
    maskNoCarbonCells() function. BC_comp = setup of the BC component 
    (see carbParams.configModelComponents())"""
    # Setting the path of txt files (depending if stochastic or not)
    if BC_comp == 1:  # stochastic
        txtFile_bcs = opj(txtFilesPath, 'bcs_{}.txt'.format(str(sample)))
    if BC_comp == 0 or nrMCsamples == 1:  # deterministic for sensitivity analysis
        txtFile_bcs = opj('bcs_detSensAnalysis.txt')
    if BC_comp == 'IPCC':  # IPCC deterministic - this is only used inside the 
        #IPCCdet_getOverallCarbonStock() function
        txtFile_bcs = opj('bcs_deterministic.txt')
    # Allocating random values in arrays
    bc = (reclassArray(climLu, txtFile_bcs))
    # Ensuring sing that cells with no stocks remain with no stocks
    bc[LUmask.mask] = 0
    # Summing total stocks in the array
    tBC = np.sum(bc)
    return bc, tBC


def getTotalStock(tSOC, tBC):
    """Return the sum of SOC and biomass carbon stocks for a given land use map.
    The output is a single number representing the computed stocks for the study
    area. Arguments: tSOC = sum of SOC array (see getTotalSOCstock function); 
    tBC = sum of BC array (see getTotalBCstock function"""
    # Computing total stocks in the study area
    tC = tSOC + tBC
    return tC

# Functions used to compute the difference of carbon stocks between scenarios
# (considering overall total carbon stocks, not cell-based)
def getDifferenceInStocksPerPairScenario(dictio_totStockPerMCr):
    """Return a dictionary representing the difference of carbon stocks (t C/ha)
    between the scenario runs with and without ethanol increase (key == scenario
    , value == array w/ difference). Arguments: dictio_totStockPerMCr = dictio
    with all scenarios (keys) and respective arrays (values) in which each element
    is a Monte Carlo sample representing the total carbon stocks in the study area"""
    # Differentiating the scenarios with and without additional ethanol production
    scCode_initial, scCode_Eth, scCode_addEth = carbParams.getScenariosType()
    # Dividing input array in two dictionaries: 'Eth' and 'addEth' scenarios
    d_scEth = {k: np.asarray(
        v) for k, v in sorted(iteritems(dictio_totStockPerMCr)) if k in scCode_Eth}
    d_scAddEth = {k: np.asarray(
        v) for k, v in sorted(iteritems(dictio_totStockPerMCr)) if k in scCode_addEth}
    # Computing the subtraction: 'addEth' minus 'Eth' scenarios
    diff_totStockPerScPair = {
        k: d_scAddEth[k] - d_scEth.get(k - 2) for k in d_scAddEth.keys()}
    return diff_totStockPerScPair


def IPCCdet_getDiffStockPerPairScenarios():
    """Compute the carbon stocks deterministically, based on IPCC mean values
    and LU maps used by Floor"""
    # Setting SOC and BC comps = 'IPCC' (necessary to get the correct txt files 
    # in getTotalSOCstock()/getTotalBCstock() functions)
    SOC_comp = 'IPCC'
    BC_comp = 'IPCC'
    # Creating dict to fill each scenario w/ total carbon stocks
    d_tSOC_det = {sc: [] for sc in scenariosDict.keys()}
    d_tBC_det = {sc: [] for sc in scenariosDict.keys()}
    d_tC_det = {sc: [] for sc in scenariosDict.keys()}
    # Getting LUmaps used in Floor's paper 
    LUmaps_Path = opj('PLUC', 'results_det')
    LUscenarios = glob.glob(opj(LUmaps_Path , 'sc*'))
    luMapsDict_det = {int(lu_path.split('sc')[-1]):opj(lu_path, 'landUse0025.map') 
                      for lu_path in sorted(LUscenarios)} # the int here makes the
    # dictionary correctly sorted by the scenario number. If run in windows, 
    # you might have to change the slash punctuation
    luMapsDict_det[0] = carbParams.getInitialLandUseMap()
    for sc, LUmapPath in sorted(iteritems(luMapsDict_det)):
        LUmap = rasterToArray(LUmapPath)
        clim_lu = (climate + LUmap)
        mask = maskNoCarbonCells(LUmap)
        soc, tSOC = getTotalSOCstock(clim_lu, clim_soil, mask, SOC_comp)
        bc, tBC = getTotalBCstock(clim_lu, mask, BC_comp)
        tC = getTotalStock(tSOC, tBC)
        # Adding results to dictionaries
        d_tSOC_det[sc].append(tSOC)
        d_tBC_det[sc].append(tBC)
        d_tC_det[sc].append(tC)
    # Difference of total stock per MCr between paired scenarios (Eth Vs addEth).
    # Output: dictio key is the code of SC w/additional Eth
    diff_tSOCeth_det = getDifferenceInStocksPerPairScenario(d_tSOC_det)
    diff_tBCeth_det = getDifferenceInStocksPerPairScenario(d_tBC_det)
    diff_tCeth_det = getDifferenceInStocksPerPairScenario(d_tC_det)
    return diff_tSOCeth_det, diff_tBCeth_det, diff_tCeth_det

# Function to get overal carbon stock stats (not cell-based)
def getStats_OverallCarbonStock(statsType, outFile=""):
    """Compute carbon stock statistics for SOC, BC and TC (with no conversion 
    unit i.e.in tonnes C/ha). Arguments: statsType: the user has to type 
    "absolute" OR diff_ethanol" OR "diff_ethanol_det". The 1st == stats for 
    absolute stocks. The 2nd == stats of difference between scenarios with 
    ethanol and with additional ethanol. The 3rd == stats to compute difference 
    between scenarios and initial state; outFile: file to save the statistics"""
    if statsType == 'absolute':
        stockTypesList = ['SOC', 'BC', 'TC']
        stockDictiosList = [d_tSOCperMCr, d_tBCperMCr, d_tCperMCr]
    elif statsType == 'diff_ethanol':
        stockTypesList = ['d_SOC', 'd_BC', 'd_TC']
        stockDictiosList = [diff_tSOCeth, diff_tBCeth, diff_tCeth]
    elif statsType == 'diff_ethanol_det':
        stockTypesList = ['d_SOC', 'd_BC', 'd_TC']
        stockDictiosList = [diff_tSOCeth_det, diff_tBCeth_det, diff_tCeth_det]
    else:
        raise Exception(
            "The 'statsType' argument is wrong. The user must set it as\
            'absolute' OR 'diff_ethanol' OR 'diff_ethanol_det'")
    # Creating empty dictionaries to add statistics later
    statsToSave = {k: [] for k in stockDictiosList[0].keys()}
    meanDict = {k: [] for k in stockDictiosList[0].keys()}
    # Computing statistics
    for i, (stockType, dictio) in enumerate(zip(stockTypesList, stockDictiosList)):
        print ("\n{}: carbon stocks + uncertainty for  ''{}'' (tons C/ha):".
            format(statsType, stockType))
        for sc, arr in sorted(iteritems(dictio)): # arr = array w/ overall stocks per MCrun
            mean = np.mean(arr)
            median = np.median(arr)
            std = np.std(arr)
            var = np.var(arr)
            percLow = np.percentile(arr, 2.5)
            percHigh = np.percentile(arr, 97.5)
            CIwidth = (percHigh - percLow)
            uncLow = ((mean - percLow) / mean) * 100
            uncHigh = ((percHigh - mean) / mean) * 100
            toSave = {'{}mean'.format(stockType): round(mean, 2), 
                      '{}median'.format(stockType): round(median, 2), 
                      '{}std'.format(stockType): round(std, 2), 
                      '{}var'.format(stockType): round(var, 2), 
                      '{}percLow'.format(stockType): round(percLow, 2), 
                      '{}percHigh'.format(stockType): round(percHigh, 2), 
                      '{}uncLow'.format(stockType): round(uncLow, 2), 
                      '{}uncHigh'.format(stockType): round(uncHigh, 2)}
            meanToReturn = {'{}mean'.format(stockType): round(mean, 2)}
            if i == 0:
                statsToSave[sc] = toSave
                meanDict[sc] = meanToReturn
            else:
                for k, v in sorted(iteritems(toSave)):
                    statsToSave[sc][k] = v
                for k, v in sorted(iteritems(meanToReturn)):
                    meanDict[sc][k] = v
            print ('SC{} -\t|\tMean: {:.3f} ({:.2f}% , {:.2f}%)'.
                   format(sc, mean, uncLow, uncHigh))
    # Converting stats to dataframe then saving to a file
    if statsType != 'diff_ethanol_det':  # no need to save deterministic stats
        df = pd.DataFrame.from_dict(
            statsToSave, orient='columns', dtype=np.float64)
        df.to_csv(opj(resultsDir, outFile), sep=';',
                  decimal='.', header=True, index=True)
    return meanDict

# Functions to compute results in GHG emissions
def getEthanolDiffProduction(filename):
    """Return a dictionary with the difference in the production of ethanol 
    between scenarios with ethanol production (eth) and with additional ethanol 
    production (addEth), according to MAGNET output. Argument: filename = the 
    file which has the MAGNET ethanol production per scenario"""
    # Reading file with magnet ethanol production for all scenarios (in mln liters)
    data = pd.read_csv(filename, sep=",", header=0, usecols=[
                       'scenario', 'ethProd_mln_l'], index_col='scenario',
                        skiprows=[1], dtype=np.int32)
    # Computing the difference between 'eth' and 'addEth'
    ethDiff = data.diff(axis=0)
    # selecting only results of "addEth" minus "Eth"
    ethDiff = ethDiff[ethDiff.ethProd_mln_l > 0]
    # Converting results from 'mln l' to 'Gj'
    ethanol_conversion = 23.4 * 1000  # mln_l to Gj
    ethDiff['dEth_Gj'] = ethDiff['ethProd_mln_l'] * ethanol_conversion
    # Converting results to dictionary
    ethDiffDict = dict(ethDiff['dEth_Gj'].round(0))
    #print ("Additional eth. production (million litres):\n", ethDiff['dEth_Gj'])
    return ethDiffDict


def getEmissionsDueToEthanol(dictio_diffSOC, dictio_diffBC, dictio_diffTC, 
                             ethDiff, outFile):
    """Return dictionaries of SOC, BC and TC arrays representing the GHG emissions
    per MC run due to the increase in ethanol production between the scenarios 
    with eth (eth) and with additional eth (addEth) production. The final unit 
    of measurement is == gram CO2-eq Mj-1 EtOH. Arguments: dictio_diffSOC, 
    dictio_diffBC and dictio_diffTC: dictionaries where k = scenario and v = array
    w/ difference in carbon stocks per MCrun; ethDiff: output of 
    getEthanolDiffProduction(); outFile: name to save in csv the results with 
    the mean value of the arrays"""
    # Creating empty dictionaries to append the outputs
    SOC_ghgDueToEth = {k: "" for k in dictio_diffSOC.keys()}
    BC_ghgDueToEth = {k: "" for k in dictio_diffBC.keys()}
    TC_ghgDueToEth = {k: "" for k in dictio_diffTC.keys()}
    # Creating empty dictio to append the mean of BC, SOC and TC, per scenario
    ghgDueToEth_mean = {k: [] for k in dictio_diffSOC.keys()}
    # Start processing
    stockTypesList = ['SOC', 'BC', 'TC']
    stockDictiosList = dictio_diffSOC, dictio_diffBC, dictio_diffTC
    for stockType, dictio in zip(stockTypesList, stockDictiosList):
        for (k1, array_carbonStock), (k2, deltaEth) in zip(sorted(iteritems(dictio)), sorted(iteritems(ethDiff))):
            if k1 == k2:  # k1 and 2 are both scenario codes
                # Converting values to gram CO2-eq Mj-1 EtOH (final measurement unit)
                # a) converting cStocks (tonne/ha) in total tonnes
                array_hectares = array_carbonStock * carbParams.getConversionUnit()
                # b) converting cStocks to emission Factor
                array_EF = array_hectares * carbParams.getEmissionFactorConversion()
                # c) # dividing by amortization period ( = 20)
                array_20yrs = array_EF / 20  
                # d) dividing by the diff in ethanol production
                array_Eth = array_20yrs / deltaEth  
                # e) getting final measurement unity
                array_FinalUnit = array_Eth * -1 * 1000
                mean = np.mean(array_FinalUnit)
                meanToDict = {'{}mean'.format(stockType): round(mean, 2)}
                if stockType == 'SOC':
                    SOC_ghgDueToEth[k1] = array_FinalUnit
                    ghgDueToEth_mean[k1] = meanToDict
                if stockType == 'BC':
                    BC_ghgDueToEth[k1] = array_FinalUnit
                if stockType == 'TC':
                    TC_ghgDueToEth[k1] = array_FinalUnit
                for k, v in sorted(iteritems(meanToDict)):
                    ghgDueToEth_mean[k1][k] = v
    # Saving means to csv file
    data = pd.DataFrame.from_dict(ghgDueToEth_mean, orient='index')
    data = data.reindex(sorted(data.columns), axis=1)
    data = data.round(2)
    data.to_csv(opj(resultsDir, outFile), sep=';',
                decimal='.', header=True, index=True)
    print ('\n'),('{}:\tEmissions due to increase in eth. prod. (CO2-eq Mj-1 EtOH):'.format(outFile)),('\n'), (data)
    return SOC_ghgDueToEth, BC_ghgDueToEth, TC_ghgDueToEth


def saveData4statisticalTest(outFile):
    """Save SOC, BC and TC arrays of GHG emissions (CO2-eq Mj-1 EtOH) for 
    statistical test in R. Arg: outFile = file to save the data"""
    # Creating lists to add outputs
    soc4statsTest, bc4statsTest, tc4statsTest = [], [], []
    data4StatsTestList = [soc4statsTest, bc4statsTest, tc4statsTest]
    # Start processing
    stockTypesList = ['SOC', 'BC', 'TC']
    stockDictiosList = [socGHGe, bcGHGe, tcGHGe]
    i = 0
    for dictio in stockDictiosList:
        for array in list(itervalues(dictio)):
            data4StatsTestList[i].append(array)
        i += 1
    for stockType, data4StatsTest in zip(stockTypesList, data4StatsTestList):
        df = pd.DataFrame(data4StatsTest)
        df.to_csv(opj('statistical_test', '{}_{}.csv'.format(outFile, stockType)),
                  sep=';', decimal='.', header=False, index=False)

# Functions for plotting boxplot
def setupBoxplot(boxplot, color):
    """Customize the boxplot built by plotBoxplot() function"""
    if (len(color) == 1) or type(color) is str:
        plt.setp(boxplot['boxes'], facecolor=color,
             edgecolor='#585858', linewidth=0.4, alpha=0.8)
    else:
        for patch, col in zip(boxplot['boxes'], color):
            patch.set_facecolor(col)
        plt.setp(boxplot['boxes'], edgecolor='#585858', linewidth=0.4, \
                 alpha=0.8)
    plt.setp(boxplot['whiskers'], color='#848484',
             linewidth=0.6, linestyle="--")
    plt.setp(boxplot['caps'], color='#848484', linewidth=1)
    plt.setp(boxplot['medians'], linewidth=1, color='#848484', linestyle="-")
    plt.setp(boxplot['means'], linestyle="-", color='k',
             linewidth=1.6, marker="+", markersize=6)
    


def getBoxplot(socGHGe, bcGHGe, socGHGe_det, bcGHGe_det): # OBS: I  couldn't 
    #figure out how to use a marker instead of a line in the "usermedians" 
    # argument of boxplot, that's why in the plot legend and boxplot don't match
    # for "deterministic" in the plot 
    """Create a boxplot of SOC/BC GHG emissions due to ethanol increase between 
    scenarios with eth ('eth')and with additional eth ('addEth')."""
    if plotBoxplot == 1:
        plt.show()
        fig, ax = plt.subplots(figsize=(8, 5))
        colors = carbParams.figureColors()
        # Setting input data for boxplots (both stoch and det results  
        socGHGe = list(itervalues(socGHGe))
        bcGHGe = list(itervalues(bcGHGe))
        socGHGe_det = [v[0] for v in list(itervalues(socGHGe_det))]
        bcGHGe_det = [v[0] for v in list(itervalues(bcGHGe_det))]
        # building boxplots
        socBoxes = plt.boxplot(socGHGe, positions=np.array(range(len(socGHGe)))
                               * 2 + 0.4, sym='', #usermedians=socGHGe_det, 
                               meanline=True, showmeans=False, patch_artist=True, 
                               widths=0.6)
        bcBoxes = plt.boxplot(bcGHGe, positions=np.array(range(len(bcGHGe))) 
                              * 2 - 0.4, sym='', #usermedians=bcGHGe_det, 
                              meanline=True, showmeans=False, patch_artist=True, 
                              widths=0.6)
        # Customizing boxplots
        setupBoxplot(bcBoxes, colors[1])
        setupBoxplot(socBoxes, colors[0])
        ymin = int(np.max(np.stack((socGHGe, bcGHGe))) * -1)
        ymax = int(np.max(np.stack((socGHGe, bcGHGe))))
        plt.yticks(np.arange(0, ymax + 1, step=10))
        plt.xticks(np.arange(0, len(socGHGe) * 2, 2), 
                   carbParams.getScenariosNames(), fontsize=8.5)
        plt.xlim(-1, len(socGHGe) * 2 - 1)
        plt.ylim(0, ymax + 1)
        # Drawing bars and plots to use in legend
        plt.bar(1, [0], color=colors[0], label='SOC', alpha=0.80)
        plt.bar(1, [0], color=colors[1], label='Biomass', alpha=0.80)
        #plt.plot([], color='k', linewidth=2, marker="+",
        #         markersize=6, label='Mean - stochastic')
        ax.plot(np.array(range(len(socGHGe)))* 2 + 0.4, socGHGe_det,
                 marker="^", markerfacecolor=colors[0],
                 markeredgecolor='#848484', linestyle='None',
                 markersize=8, color="w", label="Deterministic", zorder=99)
        ax.plot(np.array(range(len(bcGHGe)))* 2 - 0.4, bcGHGe_det,
                 marker="^", markerfacecolor=colors[1],
                 markeredgecolor='#848484', linestyle='None',
                 markersize=8, color="w", zorder=100)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], fontsize=10, loc='upper right')
        #plt.legend(fontsize=7, loc='upper right', frameon=True)
        ax.grid(which="minor", axis='x', color='#848484', linewidth=0.15)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        plt.ylabel('GHG emissions per component\n' +
                   r'($gram$ $CO_2$-$eq$/$MJ$$_E$$_t$$_O$$_H$)')
        if savePlots == 1:
            plt.savefig(opj(resultsDir, 'plot_boxplot'), dpi=700)
        if showPlots == 1:
            plt.show()

def getBoxplotThreshold(tcGHGe, tcGHGe_det):  
    """Create a boxplot of total GHG emissions due to ethanol increase between 
    scenarios with eth ('eth') and with additional eth ('addEth'). Also, includes
    the five thresholds related to Directive EU 2018/2001 """
    if plotBoxplot == 1:
        plt.show()
        # Setting input data for boxplots (both stoch and det results  
        tcGHGe = list(itervalues(tcGHGe))
        tcGHGe_det = [v[0] for v in list(itervalues(tcGHGe_det))]

        
        # The five thresholds based on (EU) 2018/2001 (as percentage)
        percent = np.array([0.8])#[0.5, 0.6, 0.65, 0.7, 0.8]
        threshold_colors = ['k']#['#c7e9b4', '#7fcdbb', '#41b6c4', '#1d91c0', '#0c2c84']        
        labels = ['threshold RED (80%)']#['50 % threshold', '60 % threshold', 
                  #'65 % threshold', '70 % threshold', '80 % threshold']
        #  Assumed Life cycle GHG emissions of eth. prod. from sugar cane in BRA
        SC_emissions = 20 # gram CO2-eq/MJ
        # Assumed default GHG emissions from gasoline 
        GSLN_emissions = 94 # gram CO2-eq/MJ
        threshold = (GSLN_emissions*(1-percent))
        # Create figure
        fig, ax = plt.subplots(figsize=(8, 5))
        # New: color based on above/below threshold
        ##colors = carbParams.figureColors()
        upper = np.percentile(np.array(tcGHGe), 97.5, axis=1)
        lower = np.percentile(np.array(tcGHGe), 2.5, axis=1)
        colors = np.where(upper < threshold, 'green', '0')
        colors = np.where(lower > threshold, 'red', colors)
        colors = np.where(colors == '0', 'orange', colors)
        # Box plots
        tcBoxes = ax.boxplot(tcGHGe, positions=np.array(range(len(tcGHGe))) 
                              * 2, sym='', #usermedians=tcGHGe_det, 
                              meanline=True, showmeans=False, patch_artist=True, 
                              widths=0.6)

        # Customizing boxplots, using previously defined colors
        setupBoxplot(tcBoxes, colors)
        ymin = int(np.max(tcGHGe)) * -1
        ymax = int(np.max(tcGHGe))
        plt.yticks(np.arange(0, ymax + 1, step=10))
        plt.xticks(np.arange(0, len(tcGHGe) * 2, 2), 
                   carbParams.getScenariosNames(), fontsize=8.5)
        plt.xlim(-1, len(tcGHGe) * 2 - 1)
        plt.ylim(0, ymax + 1)
        # Drawing bars and plots to use in legend
        #plt.bar(1, [0], color=colors[-1], label='', alpha=0.80)
        #ax.plot([], color='k', linewidth=2, marker="+",
        #         markersize=6, label='Mean - stochastic')
        ax.plot(np.array(range(len(tcGHGe)))* 2, tcGHGe_det, marker="^", \
                markerfacecolor='w', markeredgecolor='#848484', linestyle='None', \
                markersize=8, label="Deterministic", zorder=99)
        
        # Adding thresholds in plot 
        for threshold, color, label in zip(threshold, threshold_colors, labels):
            ax.axhline(y=threshold, linewidth = 0.9, 
                        linestyle = '--', color = color, label = label,
                       zorder=100)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], fontsize=10, loc='upper right')
        #plt.legend(fontsize=7, loc='upper right', frameon=True)
        plt.grid(which="minor", axis='x', color='#848484', linewidth=0.15)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        plt.ylabel('total GHG emissions\n' + r'($gram$ $CO_2$-$eq$/$MJ$$_E$$_t$$_O$$_H$)')
        if savePlots == 1:
            plt.savefig(opj(resultsDir, 'plot_boxplot_threshold'), dpi=700)
        if showPlots == 1:
            plt.show()


# Function to run sensitivity analysis
def getSensitivityAnalysis():
    """Run the sensitivity analysis based on three model components, then plot 
    it (stcSOC', 'stcBC', 'stcLUC')"""
    if plotSensAnalysis == 1:
        # Loading total CS per MCrun, getting variance and appending to dictios
        components = ['stcSOC', 'stcBC', 'stcLUC', 'stcAll']
        # Creating empty dictionary to append variances per scenario, per component
        varDict = {k: "" for k in components}
        for comp in components:
            try:
                # Loading array per SC w/ all MC runs results 
                #(representing carbon stocks in tons C/ha)
                path = opj(resultsDir, 'arr_totalStocks')
                totalStocks = {sc: np.load(
                    opj(path, "{}_sc{}_tC.npy".format(comp, sc))) for sc in scenariosDict.keys()}
            except:
                raise Exception("'{}' files not found in {}. Before run sens.\
                analysis, you have to consider running the model with the '{}'\
                component set to run stochastically. To save arrays of that run\
                the 'saveArrays' variable must be = 1".format(comp, path, comp))
            # Getting the array of the difference between scenarios with eth 
            # ('eth') and with additional eth ('addEth')
            diffStocks = getDifferenceInStocksPerPairScenario(totalStocks)
            # Getting diffStocks variance per scenario, for each component
            variancePerSc = {k: {} for k in diffStocks.keys()}
            for sc in variancePerSc.keys():
                var = np.var(diffStocks[sc]).astype(np.float64)
                variancePerSc[sc] = var
            # Adding VariancePerSc to varDict
            varDict[comp] = variancePerSc
        varData = pd.DataFrame.from_dict(
            varDict, orient='columns', dtype=np.float64)
        # Sensitivity analysis: obtaining variance fractions per component based
        #on the total variance of "stcAll" (i.e, full stochastic run of the model)
        varFract = (varData.div(varData['stcAll'], axis=0)) * 100
        # Getting the total contribution of all components in overall variance
        varFract['tot'] = varFract['stcBC'] + \
            varFract['stcSOC'] + varFract['stcLUC']
        # Getting the model interactions contribution in overall variance
        varFract['diff'] = varFract['stcAll'] - varFract['tot']
        # Plotting resultsplot_sensAnalysis
        fig, ax = plt.subplots(figsize=(10, 6))
        colors = ['#5ab4ac', '#d8b365', '#f6e8c3', '#610B0B']
        bars = np.arange(6)
        barWidth = 0.75
        labels = ['Reference\nscenario', 'High\nproductivity.', 
                  r'2$^n$$^d$ generation' + '\nsugar cane',
                  r'2$^n$$^d$ generation' + '\neucalyptus',
                  'Conservation\npolicies', 'All\nmeasures']
        # Plotting bars
        barBC = plt.bar(bars, varFract.stcBC, color=colors[0], width=barWidth)
        barSOC = plt.bar(bars, varFract.stcSOC,
                         bottom=varFract.stcBC, color=colors[1], width=barWidth)
        barLUC = plt.bar(bars, varFract.stcLUC, bottom=varFract.stcBC +
                         varFract.stcSOC, color=colors[2], width=barWidth)
        barOther = plt.bar(bars, varFract['diff'], bottom=varFract['tot'],
                           color=colors[3], width=barWidth)
        # Plot setup
        plt.xticks(bars, labels, fontsize=10)
        plt.yticks(np.arange(0, 101, step=50), fontsize=9.5)
        plt.ylabel("Component contribution in the variance of GHGE from ethanol" +
                   "production/\nTotal variance of GHGE from ethanol production (%)*",
                   fontsize=9)
        plt.legend((barBC[0], barSOC[0], barLUC[0], barOther[0]), 
                   ('Biomass', 'SOC', 'LUC', 'Model interactions'), 
                   bbox_to_anchor=(0.78, -0.12), fontsize=10, ncol=4)
        plt.grid(which="major", axis="y", color='0.7', linewidth=0.3)
        plt.subplots_adjust(bottom=0.2)
        if savePlots == 1:
            plt.savefig(opj(resultsDir, 'plot_sensAnalysis'), dpi=700)
        if showPlots == 1:
            plt.show()


# Function to transform array to map (tiff extension)
def arrayToMap(inpArray, outRaster):
    """Convert an array to a raster map and save to disk (tiff format). 
    Arguments: inpArray = array to be converted; outRaster = string with the 
    filename representing the raster map"""
    # Getting geoparameters of the original raster
    refRaster = gdal.Open(carbParams.referenceRaster())
    geotransform = refRaster.GetGeoTransform()
    originX = geotransform[0]  # top left-x
    originY = geotransform[3]  # top left-y
    pixelWidth = geotransform[1]  # w-e pixel resolution
    pixelHeight = geotransform[5]  # n-s pixel resolution
    rows = refRaster.RasterYSize
    cols = refRaster.RasterXSize
    # Generating a masked array representing cells with null value 
    # (further used when converting array to raster)
    # reference-array for masking null value
    nvArr = (refRaster.ReadAsArray()).astype(np.int32)
    # Converts to -9999 any value that doesn't represent a LU type code
    nvArr[(nvArr < 0) | (nvArr > carbParams.setNrLUtypes())] = -9999
    nvMsk = np.ma.masked_where(nvArr != -9999, nvArr)  # null value mask
    # Updating inpArray with the null value of the MskArr
    inpArray[~nvMsk.mask] = -9999
    # Creating new raster, then adding geoparameters
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(
        opj(resultsDir, "maps", outRaster + ".tiff"), cols, rows, 1, gdal.GDT_Float64)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(inpArray)
    outband.SetNoDataValue(-9999)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(refRaster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

#-------------------     1. GENERATING RANDOM VALUES     ----------------------#

txtFilesPath = opj('txt_files')
if generateRandomValues == 1:
    # Setting random seed
    randomSeed = np.random.seed(seed=62871)
    print ('Generating random values per MC...')
    generateRandomSOCref('socR_uncertainty.txt', txtFilesPath)
    generateRandomSOCfactor('socF_uncertainty.txt', txtFilesPath)
    generateRandomBCstock('bcs_uncertainty.txt', txtFilesPath)
    print ('\nChecking occurrence of negative values...')
    #checkNegRandomValues(txtFilesPath)

#-----------------     2. PROCESSING OVERALL CARBON STOCKS     ----------------#

# Input rasters
climate = rasterToArray('climate.map')
soil = rasterToArray('soil.map')

# Getting climate + soil codes' combination
clim_soil = getClimSoil(climate, soil)

# Setting particle filter mapping
particleFilter_mapping = carbParams.getMappingFromFile(
    'particle_mapping.csv', 'New_ID', 'Nr_particles_weight')

if getOverallCarbonStock == 1:
    # Creating dict w/ lists for each scenario >>> they will contain the sum of 
    # array elements of tSOC, tBC & tC, for each MC run
    d_tSOCperMCr = {sc: [] for sc in scenariosDict.keys()}
    d_tBCperMCr = {sc: [] for sc in scenariosDict.keys()}
    d_tCperMCr = {sc: [] for sc in scenariosDict.keys()}
    
    # Running the model stochastically according to 'modelRunType' variable
    for sc, scPath in sorted(iteritems(scenariosDict)):
        print ('\nSC{} - computing overall carbon stocks (t C/ha)\t...\tStart: {}'.
               format(sc, time.asctime(time.localtime(time.time()))))
        if LUC_comp == 0 or sc == 0:
            LUmapPath = getLUmap(sc, scPath)
            LUmap = rasterToArray(LUmapPath)
            mask = maskNoCarbonCells(LUmap)
            print ('\n\tLUmap (deterministic):\tsc{}/{}\n'.
                   format(sc,LUmapPath.split("/")[-1]))
            clim_lu = (climate + LUmap)
            for sample in range(1, nrMCsamples + 1):
                # Computing total stocks in the study area
                soc, tSOC = getTotalSOCstock(clim_lu, clim_soil, mask, SOC_comp,
                                             int(sample))
                bc, tBC = getTotalBCstock(clim_lu, mask, BC_comp,
                                          int(sample))
                tC = getTotalStock(tSOC, tBC)
                # Appending totals to list
                d_tSOCperMCr[sc].append(tSOC)
                d_tBCperMCr[sc].append(tBC)
                d_tCperMCr[sc].append(tC)
                if saveArrays4cb_stock == 1 and modelRunType == 'stcAll':
                    np.savez_compressed(opj(resultsDir, 'arr_cellBased', 'cb_sc{}mc{}'.
                            format(sc, sample)), soc=soc, bc=bc)                
                print ('\tSC{} MC{}\t|\t{}\t|\ttSOC: {}\ttBCS: {}\ttCS: {} (t C/ha)'.
                       format(sc, sample, modelRunType, tSOC, tBC, tC))
        if LUC_comp == 1 and sc > 0:
            # setting iterator to count nr of mcRuns (necessary to change the 
            # LU map according to particle filter
            mcRun_acc = 0
            LUmapsDict = carbParams.getStochasticLUmaps(
                scPath, particleFilter_mapping, '2030_*')
            for descr in sorted(LUmapsDict.items()):
                mapCode, mapPath, mapMCruns = descr[0], descr[1][0], descr[1][1]
                LUmapPath = getLUmap(sc, scPath, mapPath)
                LUmap = rasterToArray(LUmapPath)
                clim_lu = (climate + LUmap)
                mask = maskNoCarbonCells(LUmap)
                # the MC range used to run the current LU map
                mcRange_curr = getLUmapMCrange(mapCode, mapMCruns, mcRun_acc)
                print(mcRange_curr)
                mcRun_acc += mapMCruns
                print ('\n\tLUmap nr {} (stoch): {}\t|\tsamples: {} (runs {}-{})\n'.
                       format(mapCode, LUmapPath.split("ts/")[-1], mapMCruns, 
                              min(mcRange_curr), max(mcRange_curr)))
                for sample in mcRange_curr:
                    # Computing total stocks in the study area
                    soc, tSOC = getTotalSOCstock(
                        clim_lu, clim_soil, mask, SOC_comp, int(sample))
                    bc, tBC = getTotalBCstock(clim_lu, mask, BC_comp, int(sample))
                    tC = getTotalStock(tSOC, tBC)
                    # Appending totals to list
                    d_tSOCperMCr[sc].append(tSOC)
                    d_tBCperMCr[sc].append(tBC)
                    d_tCperMCr[sc].append(tC)
                    if saveArrays4cb_stock == 1 and modelRunType == 'stcAll':
                        np.savez_compressed(opj(resultsDir, 'arr_cellBased', 'cb_sc{}mc{}'.
                                format(sc, sample)), soc=soc, bc=bc)
                    print ('\tSC{} MC{}\t|\t{}\t|\ttSOC: {}\ttBCS: {}\ttCS: {} (ton C/ha)'.
                           format(sc, sample, modelRunType, tSOC, tBC, tC))
        if saveArrays4ov_stock == 1:
            # Saving list w/ total stock per MCr as array
            np.save(opj(resultsDir, 'arr_totalStocks', "{}_sc{}_tSOC".format(
                modelRunType, sc)), np.asarray(d_tSOCperMCr[sc]))
            np.save(opj(resultsDir, 'arr_totalStocks', "{}_sc{}_tBC".format(
                modelRunType, sc)), np.asarray(d_tBCperMCr[sc]))
            np.save(opj(resultsDir, 'arr_totalStocks', "{}_sc{}_tC".format(
                modelRunType, sc)), np.asarray(d_tCperMCr[sc]))

if getOverallCarbonStock == 0:
    # Opening numpy arrays with total stocks of each MC sample
    d_tSOCperMCr = {sc: np.load(opj(resultsDir, 'arr_totalStocks', '{}_sc{}_tSOC.npy'.
                    format(modelRunType, sc))) for sc in scenariosDict.keys()}
    d_tBCperMCr = {sc: np.load(opj(resultsDir, 'arr_totalStocks', '{}_sc{}_tBC.npy'.
                    format(modelRunType, sc))) for sc in scenariosDict.keys()}
    d_tCperMCr = {sc: np.load(opj(resultsDir, 'arr_totalStocks', '{}_sc{}_tC.npy'.
                    format(modelRunType, sc))) for sc in scenariosDict.keys()}

# Difference of total stock per MCr between paired scenarios (Eth Vs addEth). 
# Output: dictio key is the code of SC w/additional Eth
diff_tSOCeth = getDifferenceInStocksPerPairScenario(d_tSOCperMCr)
diff_tBCeth = getDifferenceInStocksPerPairScenario(d_tBCperMCr)
diff_tCeth = getDifferenceInStocksPerPairScenario(d_tCperMCr)

# Getting the difference but now running FULL deterministic(based on IPCC means)
diff_tSOCeth_det, diff_tBCeth_det, diff_tCeth_det = IPCCdet_getDiffStockPerPairScenarios()

#---------------     3. GETTING OVERALL STATISTICS & SAVING     ---------------#

scStats = getStats_OverallCarbonStock(
    'absolute', outFile = '{}_absolute_Stats.csv'.format(modelRunType))
scDiffStats = getStats_OverallCarbonStock(
    'diff_ethanol', outFile = "{}_dEth_Stats.csv".format(modelRunType))
scDiffStats_det = getStats_OverallCarbonStock('diff_ethanol_det')

#-------     4. COMPUTING OVERALL EMISSIONS FROM ETHANOL PRODUCTION     -------#

# Getting emissions due to ethanol production (for SOC, BC and TC) - STOCHASTIC
ethDiff = getEthanolDiffProduction('magnet_output.txt')
socGHGe, bcGHGe, tcGHGe = getEmissionsDueToEthanol(diff_tSOCeth, diff_tBCeth,
                        diff_tCeth, ethDiff, '{}_GHGe.csv'.format(modelRunType))


# Now getting deterministic emissions (IPCC)  
socGHGe_det, bcGHGe_det, tcGHGe_det = getEmissionsDueToEthanol(
    diff_tSOCeth_det, diff_tBCeth_det, diff_tCeth_det, ethDiff, 'determ_GHGe.csv')

#-----   5. OUTPUTS (BOXPLOT, SENS. ANALYSIS, DATA FOR STATISTICAL TEST)  -----#

getBoxplot(socGHGe, bcGHGe, socGHGe_det, bcGHGe_det)
getBoxplotThreshold(tcGHGe, tcGHGe_det)
getSensitivityAnalysis()
if modelRunType == "stcAll":
    saveData4statisticalTest('data4StatsTest')

#---------------     6. PROCESSING CELL-BASED CARBON STOCKS     ---------------#

# Functions used to get cell-based carbon stocks with uncertainty
def createArr4cellBasedStocks(refRaster):
    """Create an array filled with zero, necessary to compute mean and std of 
    cell-based carbon stocks. The shape size is based on the number of columns 
    and rows of a raster file."""
    refRaster = gdal.Open(refRaster)
    rows = refRaster.RasterYSize
    cols = refRaster.RasterXSize
    zerosArr = np.zeros([rows, cols], dtype=np.float64)
    return zerosArr

def getCellBasedMeanStocks(csType, sc, saveToFile):
    """Compute the Cell-based mean of a given scenario. Arguments: csType = it 
    must be 'soc' OR 'bc' OR 'tc'; sc: scenario number; saveToFile: string with  
    the name that array and map will be saved"""
    #if not csType == 'soc' or csType == 'bc' or csType == 'tc':
        #raise Exception("csType argument set to '{}'. It must be 'soc', 'bc' or\
         #'tc'.".format(csType))
    npzFiles = glob.glob(opj(resultsDir, 'arr_cellBased', 'cb_sc{}mc*'.
                             format(sc)))
    npzFilesDict = {
        int((v.split('mc'))[-1].split('.')[0]): v for v in sorted(npzFiles)}
    # Computing total cell-based mean and per LU type
    cs_accum = createArr4cellBasedStocks(carbParams.referenceRaster())
    if sc == 0:
        # Creating Mask
        LUmapPath = carbParams.getInitialLandUseMap()
        LUmap = rasterToArray(LUmapPath)
        mask = maskNoCarbonCells(LUmap)
        # Starting computation
        for sample, fl in sorted(iteritems(npzFilesDict)):
            npData = np.load(fl)
            if csType == 'soc' or csType == 'bc':
                cs = npData[csType]
            if csType == 'tc':
                soc = npData['soc']
                bc = npData['bc']
                cs = soc + bc
            cs[mask.mask] = 0
            cs_accum = cs + cs_accum
            print ('\tmean - file {} ...'.format(npzFilesDict[sample].split("_")[-1]))
            if np.isnan(cs_accum).any() == True:
                print ('SC{} - NaN in cs_accum when summing {} with file {}'.
                       format(sc, csType, npzFilesDict[sample]))                
    if sc > 0:
        # set iterator of mcRuns (used to change the LUmap based on part. filter)
        mcRun_acc = 0
        LUmapsDict = carbParams.getStochasticLUmaps(
            scPath, particleFilter_mapping, '2030_*')
        for descr in sorted(LUmapsDict.items()):
            mapCode, mapPath, mapMCruns = descr[0], descr[1][0], descr[1][1]
            # Creating LU mask
            LUmapPath = getLUmap(sc, scPath, mapPath)
            LUmap = rasterToArray(LUmapPath)
            mask = maskNoCarbonCells(LUmap)
            # Getting the MC range to run the current LU map
            mcRange_curr = getLUmapMCrange(mapCode, mapMCruns, mcRun_acc)
            mcRun_acc += mapMCruns
            #print ('runs {}-{}'.format(min(mcRange_curr), max(mcRange_curr)))
            # Starting computation
            for sample in mcRange_curr:
                npData = np.load(npzFilesDict[sample])
                if csType == 'soc' or csType == 'bc':
                    cs = npData[csType]
                if csType == 'tc':
                    soc = npData['soc']
                    bc = npData['bc']
                    cs = soc + bc
                cs[mask.mask] = 0
                cs_accum = cs + cs_accum                    
                print ('\tmean - {}, {}, {}, file {} ...'.format(
                    csType, sample, LUmapPath, npzFilesDict[sample].split("_")[-1]))
                if np.isnan(cs_accum).any() == True:
                    print ('SC{} - NaN in cs_accum when summing {} with file {}'.
                           format(sc, csType, npzFilesDict[sample]))                    
    csMean_cell = (cs_accum / len(npzFiles))  # OK! (It has 0 values)
    # Saving array 
    np.save(opj(resultsDir, 'maps', saveToFile), csMean_cell)
    # Converting array to map and saving map
    arrayToMap(csMean_cell, saveToFile)
    return csMean_cell

def getCellBasedStdStocks(csType, cellBased_mean, sc, saveToFile):
    """Compute the Cell-based standard deviation of a given scenario. Arguments:
    csType = it must be 'soc' OR 'bc' OR 'tc'; cellBased_mean = cell-based mean
    array; sc =  scenario number; saveToFile =  string with the name that array 
    and map will be saved"""
    #if not csType == 'soc' or csType == 'bc' or csType == 'tc':
        #raise Exception("csType argument set to '{}'. It must be 'soc', 'bc' or\
         #'tc'.".format(csType))
    npzFiles = glob.glob(opj(resultsDir, 'arr_cellBased', 'cb_sc{}mc*'.
                             format(sc)))
    npzFilesDict = {
        int((v.split('mc'))[-1].split('.')[0]): v for v in sorted(npzFiles)}
    # Computing total cell-based mean and per LU type
    cs_sqAccum = createArr4cellBasedStocks(carbParams.referenceRaster())
    if sc == 0:
        # Creating mask 
        LUmapPath = carbParams.getInitialLandUseMap()
        LUmap = rasterToArray(LUmapPath)
        mask = maskNoCarbonCells(LUmap)
        # Starting computation
        for sample, fl in sorted(iteritems(npzFilesDict)):
            npData = np.load(fl)
            if csType == 'soc' or csType == 'bc':
                cs = npData[csType]
            if csType == 'tc':
                soc = npData['soc']
                bc = npData['bc']
                cs = soc + bc
            # Using mask to remove zero values when computing squares
            csMask = np.ma.masked_values(cs,0)
            csSq = (csMask - cellBased_mean) ** 2
            cs_sqAccum = csSq + cs_sqAccum            
            print ('\tstd - file {} ...'.format(npzFilesDict[sample].split("_")[-1]))
            if np.isnan(cs_sqAccum).any() == True:
                print ('SC{} - NaN in cs_accum when summing {} with file {}'.
                       format(sc, csType, npzFilesDict[sample]))                
    if sc > 0:
        # set iterator of mcRuns (used to change the LUmap based on part. filter)
        mcRun_acc = 0
        LUmapsDict = carbParams.getStochasticLUmaps(
            scPath, particleFilter_mapping, '2030_*')
        for descr in sorted(LUmapsDict.items()):
            mapCode, mapPath, mapMCruns = descr[0], descr[1][0], descr[1][1]
            # Creating mask 
            LUmapPath = getLUmap(sc, scPath, mapPath)
            LUmap = rasterToArray(LUmapPath)
            mask = maskNoCarbonCells(LUmap)
            # Getting the MC range to run the current LU map
            mcRange_curr = getLUmapMCrange(mapCode, mapMCruns, mcRun_acc)
            mcRun_acc += mapMCruns
            #print ('runs {}-{}'.format(min(mcRange_curr), max(mcRange_curr)))
            # Starting computation 
            for sample in mcRange_curr:
                npData = np.load(npzFilesDict[sample])
                if csType == 'soc' or csType == 'bc':
                    cs = npData[csType]
                if csType == 'tc':
                    soc = npData['soc']
                    bc = npData['bc']
                    cs = soc + bc
                # Using mask to remove zero values when computing squares
                csMask = np.ma.masked_values(cs,0) 
                csSq = (csMask - cellBased_mean) ** 2
                cs_sqAccum = csSq + cs_sqAccum                  
                print ('\tstd - {}, {}, {}, file {} ...'.format(
                    csType, sample, LUmapPath, npzFilesDict[sample].split("_")[-1]))
                if np.isnan(cs_sqAccum).any() == True:
                    print ('SC{} - NaN in cs_accum when summing {} with file {}'.
                           format(sc, csType, npzFilesDict[sample]))                    
    csStd_cell = np.sqrt(cs_sqAccum/len(npzFiles))
    csStd_cell = np.asarray(csStd_cell)
    # Saving array 
    np.save(opj(resultsDir, 'maps', saveToFile), csStd_cell)
    # Converting array to map and saving map
    arrayToMap(csStd_cell, saveToFile)
    return csStd_cell

def getDiffCellBasedMeanStocks(csType, scEth, scAddEth, saveToFile):
    """ """
    npzFilesEth = glob.glob(opj(resultsDir, 'arr_cellBased', 'cb_sc{}mc*'.
                             format(scEth)))
    
    npzFilesAddEth = glob.glob(opj(resultsDir, 'arr_cellBased', 'cb_sc{}mc*'.
                             format(scAddEth)))    
    npzDictEth = {
        int((v.split('mc'))[-1].split('.')[0]): v for v in sorted(npzFilesEth)}
    npzDictAddEth = {
        int((v.split('mc'))[-1].split('.')[0]): v for v in sorted(npzFilesAddEth)}    
    # To confirm if arrays are sharing the same MCruun, use:
    #for (k1, v1), (k2,v2) in zip(npzDictEth.iteritems(), npzDictAddEth.iteritems()):
        #print (k1, k2,"-----", v1,v2)     
    cs_accum = createArr4cellBasedStocks(carbParams.referenceRaster())
    # set iterator of mcRuns (used to change the LUmap based on part. filter)
    mcRun_acc = 0
    LUmapsDict = carbParams.getStochasticLUmaps(scAddEthPath, particleFilter_mapping, '2030_*')
    for descr in sorted(LUmapsDict.items()):
        mapCode, mapPath, mapMCruns = descr[0], descr[1][0], descr[1][1]
        # Creating mask of LU map regarding scenario w/ additional ethanol
        LUmapPath = getLUmap(scAddEth, scAddEthPath, mapPath)
        LUmap = rasterToArray(LUmapPath)
        mask = maskNoCarbonCells(LUmap)
        # the MC range used to run the current LU map
        mcRange_curr = getLUmapMCrange(mapCode, mapMCruns, mcRun_acc)
        mcRun_acc += mapMCruns
        #print ('runs {}-{}'.format(min(mcRange_curr), max(mcRange_curr)))
        #Starting computations
        for sample in mcRange_curr:
            npDataEth = np.load(npzDictEth[sample])
            npDataAddEth = np.load(npzDictAddEth[sample])
            if csType == 'soc' or csType == 'bc':
                csEth = npDataEth[csType]
                csAddEth = npDataAddEth[csType]
            if csType == 'tc':
                socEth = npDataEth['soc']
                bcEth = npDataEth['bc']
                csEth = socEth + bcEth
                socAddEth = npDataAddEth['soc']
                bcAddEth = npDataAddEth['bc']
                csAddEth = socAddEth + bcAddEth            
            csDiff = csAddEth - csEth
            csDiff[mask.mask] = 0
            cs_accum = csDiff + cs_accum                    
            print ('\tDiff mean of {} ({} vs {}) - {}, {}, files {} & {} ...'.format(
                csType, scAddEth, scEth, sample, LUmapPath, 
                npzDictAddEth[sample].split("_")[-1], npzDictEth[sample].split("_")[-1]))
            if np.isnan(cs_accum).any() == True:
                print ('SC{} - NaN in cs_accum when summing {} with file {}'.
                       format(sc, csType, npzDictAddEth[sample]))                    
    csDiffMean_cell = (cs_accum / len(npzFilesAddEth))  # OK! (It has 0 values)
    # Saving array 
    np.save(opj(resultsDir, 'maps', saveToFile), csDiffMean_cell)
    # Converting array to map and saving map
    arrayToMap(csDiffMean_cell, saveToFile)
    return csDiffMean_cell
    
def getDiffCellBasedStdStocks(csType, scEth, scAddEth, csDiffMean_cell, saveToFile):
    """ """
    npzFilesEth = glob.glob(opj(resultsDir, 'arr_cellBased', 'cb_sc{}mc*'.
                             format(scEth)))
    
    npzFilesAddEth = glob.glob(opj(resultsDir, 'arr_cellBased', 'cb_sc{}mc*'.
                             format(scAddEth)))    
    npzDictEth = {
        int((v.split('mc'))[-1].split('.')[0]): v for v in sorted(npzFilesEth)}
    npzDictAddEth = {
        int((v.split('mc'))[-1].split('.')[0]): v for v in sorted(npzFilesAddEth)}    
    # To confirm if arrays are sharing the same MCruun, use:
    #for (k1, v1), (k2,v2) in zip(npzDictEth.iteritems(), npzDictAddEth.iteritems()):
        #print (k1, k2,"-----", v1,v2)     
    cs_sqAccum = createArr4cellBasedStocks(carbParams.referenceRaster())
    # set iterator of mcRuns (used to change the LUmap based on part. filter)
    mcRun_acc = 0
    LUmapsDict = carbParams.getStochasticLUmaps(scAddEthPath, particleFilter_mapping, '2030_*')
    for descr in sorted(LUmapsDict.items()):
        mapCode, mapPath, mapMCruns = descr[0], descr[1][0], descr[1][1]
        # Creating mask of LU map regarding scenario w/ additional ethanol
        LUmapPath = getLUmap(scAddEth, scAddEthPath, mapPath)
        LUmap = rasterToArray(LUmapPath)
        mask = maskNoCarbonCells(LUmap)
        # the MC range used to run the current LU map
        mcRange_curr = getLUmapMCrange(mapCode, mapMCruns, mcRun_acc)
        mcRun_acc += mapMCruns
        #print ('runs {}-{}'.format(min(mcRange_curr), max(mcRange_curr)))
        #Starting computations
        for sample in mcRange_curr:
            npDataEth = np.load(npzDictEth[sample])
            npDataAddEth = np.load(npzDictAddEth[sample])
            if csType == 'soc' or csType == 'bc':
                csEth = npDataEth[csType]
                csAddEth = npDataAddEth[csType]
            if csType == 'tc':
                socEth = npDataEth['soc']
                bcEth = npDataEth['bc']
                csEth = socEth + bcEth
                socAddEth = npDataAddEth['soc']
                bcAddEth = npDataAddEth['bc']
                csAddEth = socAddEth + bcAddEth            
            csDiff = csAddEth - csEth
            # mask to remove zero values when computing squares
            maskDiff = (np.ma.masked_values(csDiff, 0))
            sqDiff = (maskDiff - csDiffMean_cell) ** 2
            cs_sqAccum = sqDiff + cs_sqAccum
            print ('\tDiff std of {} ({} vs {}) - {}, {}, files {} & {} ...'.format(
                csType, scAddEth, scEth, sample, LUmapPath, 
                npzDictAddEth[sample].split("_")[-1], npzDictEth[sample].split("_")[-1]))
            if np.isnan(cs_accum).any() == True:
                print ('SC{} - NaN in cs_accum when summing {} with file {}'.
                       format(sc, csType, npzDictAddEth[sample]))                    
    csDiffStd_cell = np.sqrt(cs_sqAccum / len(npzFilesAddEth))
    csDiffStd_cell = np.asarray(csDiffStd_cell)
    # Saving array 
    np.save(opj(resultsDir, 'maps', saveToFile), csDiffStd_cell)
    # Converting array to map and saving map
    arrayToMap(csDiffStd_cell, saveToFile)
    return csDiffStd_cell    
    
#------------     DATA PROCESSING (cell-based carbon stocks)     --------------#

if getCellBasedCarbonStock == 1:
    # OBS: to use all the 10.0000 saved cell-based arrays per scenario in the 
    # cell-based functions, the  carbParams.nrMCsamples must be set to 10000 
    # (I did not have time to make it more "flexible")
    LUC_comp = 1
    #for sc, scPath in sorted(iteritems(scenariosDict)):
        #print ('\nSC{} - computing cell-based carbon stocks (tonne C/ha)\t...\t\
        #Start time: {}'.format(sc, time.asctime(time.localtime(time.time()))))
        ## Getting mean & std mapsmaps
        #mSOC_cell = getCellBasedMeanStocks('soc', sc, 
                            #saveToFile = "sc{}_{}_mean".format(sc, 'SOC'))
        #stdSOC_cell = getCellBasedStdStocks('soc', mSOC_cell, sc, 
                                #saveToFile = "sc{}_{}_std".format(sc, 'SOC'))
        #mBC_cell = getCellBasedMeanStocks('bc', sc, 
                                #saveToFile = "sc{}_{}_mean".format(sc, 'BC'))
        #stdBC_cell = getCellBasedStdStocks('bc', mBC_cell, sc, 
                                #saveToFile = "sc{}_{}_std".format(sc, 'BC'))
        #mTC_cell = getCellBasedMeanStocks('tc', sc, 
                                #saveToFile = "sc{}_{}_mean".format(sc, 'TC'))
        #stdTC_cell = getCellBasedStdStocks('tc', mTC_cell, sc, 
                                #saveToFile = "sc{}_{}_std".format(sc, 'TC'))
        #print ('end', time.asctime(time.localtime(time.time())))
        
    # Cell-based diffs
    scInitial, scEth, scAddEth = carbParams.getScenariosType()
    scEthDict = {k: v for k, v in sorted(iteritems(scenariosDict)) if k in scEth}
    scAddEthDict = {k: v for k, v in sorted(iteritems(scenariosDict)) if k in scAddEth}
    for (scEth, scEthPath), (scAddEth, scAddEthPath) in zip(
        sorted(iteritems(scEthDict)), sorted(iteritems(scAddEthDict))):
        socDiffMean = getDiffCellBasedMeanStocks('soc', scEth, scAddEth, "sc{}_{}_diff_mean".format(scAddEth, 'SOC'))
        socDiffStd = getDiffCellBasedStdStocks('soc', scEth, scAddEth, csDiffMean_cell, "sc{}_{}_diff_std".format(scAddEth, 'SOC'))
        bcDiffMean = csgetDiffCellBasedMeanStocks('bc', scEth, scAddEth, "sc{}_{}_diff_mean".format(scAddEth, 'BC'))
        bcDiffStd = getDiffCellBasedStdStocks('bc', scEth, scAddEth, csDiffMean_cell, "sc{}_{}_diff_std".format(scAddEth, 'BC'))
        tcDiffMean = getDiffCellBasedMeanStocks('tc', scEth, scAddEth, "sc{}_{}_diff_mean".format(scAddEth, 'TC'))
        tcDiffStd = getDiffCellBasedStdStocks('tc', scEth, scAddEth, csDiffMean_cell, "sc{}_{}_diff_std".format(scAddEth, 'TC'))

print("Finished at: {}".format(time.asctime(time.localtime(time.time()))))
