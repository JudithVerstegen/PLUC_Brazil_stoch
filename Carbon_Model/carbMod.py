import time
start = time.time()

print 'Starting ...\n'

import gdal
import osr
import glob
import os
import numpy as np
import pandas as pd
import parameters
import pickle
import scipy.stats as ss
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
np.set_printoptions(suppress=True)

#-------------------     PARAMETERS     -------------------#

# Creating outputs' path 
resultsDir = parameters.setOutputPath()

# Dictionary corresponding to the path for scenarios with LUC maps
scenariosDict = parameters.getScenariosPaths(r'/home/rber/Work/PLUC_Brazil_stoch-master/Results','sc*')
scenariosDict[0] = parameters.getInitialLandUseMap() # adding to dictio the deterministic initial LUmap from PLUC (key: 0, valeu: LUmap path)

# Setting model run type and components 
modelRunType, SOC_comp, BC_comp, LUC_comp = parameters.configModelComponents()

nrMCsamples = parameters.getMonteCarloSamples()
generateRandomValues = 0 
getOverallCarbonStock = 0 # if = 0, then the model loads saved files in 'arr_totalStocks' folder
#getCellBasedStock = 0 --- Not ready
plotFigures = 1
showPlots = 1
savePlots = 0
saveArrays = 0

print 'Nr of scenarios (including initial): {}.\nMC simulation: {} samples per scenario.'.format(len(scenariosDict), nrMCsamples)

# Functions related to random values
def generateRandomSOCref(nrMCsamples, inpTxt, outPath):
    """For each Monte Carlo sample, create a text file filled with random values for 'SOC reference' w/ a factor of 1000. Arguments: inpTxt = Input file with uncertainty range"""
    mcRange = np.arange(nrMCsamples)
    for i in mcRange:
        sample = i+1
        outTxt = open(os.path.join(outPath, 'socR_{}.txt'.format(str(sample))), 'w')
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
                    print 'Neg/zero value identified for SOCr...Getting positive value ...'
                    SOCr = np.random.normal(SOCr_mean, SOCr_std)
                outTxt.write('{}\t{}\n'.format(code_id, str(int(SOCr * 1000))))
        outTxt.close()
                
def generateRandomSOCfactor(inpTxt, nrMCsamples, outPath):
    """For each Monte Carlo sample, create a text file filled with random values for 'SOC factors' w/ a factor of 1000. Arguments: inpTxt = Input file with uncertainty range"""
    mcRange = np.arange(nrMCsamples)
    for i in mcRange:
        sample = i+1
        outTxt = open(os.path.join(outPath, 'socF_{}.txt'.format(str(sample))), 'w')
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
                    fLU = np.random.normal(fLU_mean, fLU_std)
                    fMG = np.random.normal(fMG_mean, fMG_std)
                    fI = np.random.normal(fI_mean, fI_std)
                    while (fLU <= 0 or fMG <= 0 or fI <= 0):
                        fLU = np.random.normal(fLU_mean, fLU_std)
                        fMG = np.random.normal(fMG_mean, fMG_std)
                        fI = np.random.normal(fI_mean, fI_std)
                        print 'Neg. value identified for SOC factor(s). Getting positive value(s)...'
                    SOCf = fLU * fMG * fI
                outTxt.write('{}\t{}\n'.format(code_id, str(int(SOCf * 1000))))
        outTxt.close()

def generateRandomBCstock(inpTxt, nrMCsamples, outPath):
    """For each Monte Carlo sample, create a text file filled with random values for 'biomass carbon stocks' (bcs) w/ a factor of 1000. Arguments: inpTxt = Input file with uncertainty range"""
    mcRange = np.arange(nrMCsamples)
    for i in mcRange:
        sample = i + 1
        outTxt = open(os.path.join(outPath, 'bcs_{}.txt'.format(str(sample))), 'w')
        with open(inpTxt, 'r') as inp:
            inpHeader = inp.readline()
            outHeader = outTxt.write('id\tvalue\n')
            for row in inp:
                rowValue = row.split()
                code_id = rowValue[0]
                agb_mean, agb_std = float(rowValue[1]), float(rowValue[2])
                r2s_mean, r2s_std = float(rowValue[5]), float(rowValue[6])
                r2s_min, r2s_max = 0, float(rowValue[8])  # Set zero as min because the r2s PDFs have neg. values
                cf_mean, cf_std = float(rowValue[9]), float(rowValue[10])
                BCS_mean, BCS_std = float(rowValue[13]), float(rowValue[14])
                agb = np.random.normal(agb_mean, agb_std)
                while agb < 0:
                    print 'Neg. value identified for AGB. Getting positive value ...'
                    agb = np.random.normal(agb_mean, agb_std)
                r2s = np.random.normal(r2s_mean, r2s_std)
                while r2s < 0:
                    print 'Neg. value identified for r2s. Getting positive value by truncation ...'
                    if r2s_std != 0:  # scipy.stats does not allow computation of std = 0
                        a, b = (r2s_min - r2s_mean) / r2s_std, (r2s_max - r2s_mean) / r2s_std  # setting truncation
                        r2s = ss.truncnorm.rvs(a, b, r2s_mean, r2s_std)
                    else:
                        r2s = np.random.normal(r2s_mean, r2s_std)
                cf = np.random.normal(cf_mean, cf_std)
                while cf < 0:
                    print 'Neg. value identified for CF. Getting positive value ...'
                    cf = np.random.normal(cf_mean, cf_std)
                bgb = agb * r2s
                BCS = (agb + bgb) * cf
                if BCS == 0:  # if clause necessary to compute BCS for CROPs and ABANDONNED LAND (IPCC doesn't give AGB, R2S and CF for those LU types)
                    BCS = (np.random.normal(BCS_mean, BCS_std))
                    while BCS < 0:
                        print 'Neg. value identified for BCS. Getting positive value ...'
                        BCS = np.random.normal(BCS_mean, BCS_std)
                outTxt.write('{}\t{}\n'.format(code_id, str(int(BCS * 1000))))
        outTxt.close()

# Functions used to get overall total carbon stocks (not cell-based)
def checkNegRandomValues(txtFilesFolder):
    """Check the occurrence of negative values in all the txt files generated by the random values' functions"""  
    files2check = sorted(glob.iglob(os.path.join(txtFilesFolder, '*')))
    for txtFile in files2check:
        df = pd.read_csv(txtFile, sep ='\t', header=0)
        if df[df.value < 0].empty == False:
            print txtFile.split('txt_files/')[1], '\n', df[df.value < 0]
            raise Exception('Negative value found in file {}'.format(txtFile.split('txt_files/')[1]))    
def rasterToArray(inpRaster):
    """Return a numpy array based on a raster file. The 'No data' of the raster are converted to zero"""
    inpRaster = gdal.Open(inpRaster)
    unqBand = inpRaster.GetRasterBand(1)
    noData = unqBand.GetNoDataValue()
    npArr = inpRaster.ReadAsArray()
    npArr[npArr == noData] = 0 # Converting noData values to zero
    if npArr[npArr < 0].any() == True: # Making sure that no negative value is within the array
        npArr[npArr < 0] = 0  
    return npArr

def reclassArray(npRaster, txtFile):
    """ Return a reclassified array based on a text file representing codes / random values ​for SOCF, SOCr or BCS. The codes given by the text files must encompass all the codes given by the array. Also, it removes the factor of 1000 of the text files. 
    Arguments: npRaster = numpy array representing codes of a raster file; txtFile = text file with codes and random values used for reclassification."""
    cvrm = np.loadtxt(txtFile, dtype=np.int32, skiprows=1)  # 'Code-value' reference matrix: convert the .txt into array  
    code, value = cvrm.T  # Transposed array representing code and value  
    lut = np.arange(npRaster.max() + 1)  # Lookup table with all the codes of npRaster. '+1' is necessary to encompass all the codes
    lut[code] = value # Assign a value for the lookup table based on the associated code   
    reclassifArr = (lut[npRaster]).astype(np.float32) / 1000 # Reclassify the npRaster filled with codes according to the values of the lookup tables 
    return reclassifArr

def maskNoCarbonCells(LUmap_array):
    """LUmask: masked array of a given land use map in which the cells that are assumed to have no carbon (i.e. lu types with no carbon stocks) are masked. Argument: LUmap_array = array of a given land use map"""
    noCarbonValues = parameters.LUtypesToMask()
    mask = np.logical_or.reduce([LUmap_array== value for value in noCarbonValues])
    maskedArray = np.ma.masked_where(mask == True, LUmap_array)
    return maskedArray

def getLUmap(sc, scenario_path, LUmapStoch=""): ### <============== Function to get the LUmap depending on the scenario nr LUC component
    """Returns the Land Use map(s) to be used in the carbon stock computation. The output depend whether the 'LUC_comp' variable is set to run deterministically or stochastically. Arguments: sc = scenario; scenario_path = scenario folder; LUmapStoch = only used when the LUC component is set to run stochastically."""
    if sc == 0:
        LUmap = parameters.getInitialLandUseMap()
    if sc > 0:
        if LUC_comp == 1: # stochastic
            LUmap = LUmapStoch
        if LUC_comp == 0: # deterministic
            LUmap = parameters.getDeterministicLUmap(scenario_path, particleFilter_mapping)  
    return LUmap

def getLUmapMCrange(nrLUmap, mcRuns_LUmap, mcRun_accumulated):
    """Return the Monte Carlo samples Range for each LUC map. It is used when LUC_comp == 1 (stochastic)). Arguments: LUmapsDict = dictionary with the LUC maps to be used and related MC runs based on Particle Filter Mapping (see parameters.getStochasticLUmap function); nrLUmap = code of the LU map; mcRuns_LUmap: total mcRuns of a given land use map;  mcRun_accumulated = the current accumulated MCrun of the for loop (since this function is run in a for loop)"""
    # Getting the exact Monte Carlo samples range in which the given land use map will run  
    mcRunToStart = mcRun_accumulated+1
    mcRun_accumulated += mcRuns_LUmap
    mcRunToStop = mcRun_accumulated+1
    mcRange = np.arange(mcRunToStart, mcRunToStop)
    return mcRange 

def getNrCellsPerLUtype(arrayLU): ### <============== Function to get the Nr of Cells per LU type
    """Return a dictionary with the total amount of cells (dictio values) per LU type (dictio keys). 
    Arguments: arrayLU = array representing a land use map"""
    luUnq, cellCounts = np.unique(arrayLU, return_counts=True)
    CellPerLu = dict(zip(luUnq, cellCounts))
    #if 0 in AreaPerLu:  # Zero represents NoData, therefore is not LU type
        #del AreaPerLu[0]
    return CellPerLu
def getTotalSOCstock(climLu, climSoil, LUmask, SOC_comp):
    """Returns the total SOC stock for a given land use map (soc == array with total stock per cell, tSOC == sum of array). The outputs depend whether the 'SOC_comp' variable is set to run deterministically or stochastically. Arguments: climLU = array of climate + lu maps combination (for SOC factor); climSoil: array of climate + soil maps combination (for SOC reference); LUmask: it is the output of the maskNoCarbonCells() function; SOC_comp = setup of the SOC component (see parameters.configModelComponents())"""
    # Setting the path of txt files (depending if stochastic or not)
    if SOC_comp == 1: # stochastic
        txtFile_socF = os.path.join(txtFilesPath, 'socF_{}.txt'.format(str(sample)))
        txtFile_socR = os.path.join(txtFilesPath, 'socR_{}.txt'.format(str(sample)))
    if SOC_comp == 0: # deterministic for sensitivity analysis
        txtFile_socF = os.path.join('socF_deterministic.txt')
        txtFile_socR = os.path.join('socR_deterministic.txt')
    if SOC_comp == 'IPCC': # # IPCC deterministic - this is only used in the IPCCdet_getOverallCarbonStock() function (now files are == sensAnalysis but we might differentiate afterwards)   
        txtFile_socF = os.path.join('socF_deterministic.txt')
        txtFile_socR = os.path.join('socR_deterministic.txt')    
    #print txtFile_socF, txtFile_socR
    # Allocating random values in arrays (SOC factors, then SOC reference values)
    socF = reclassArray(climLu, txtFile_socF)
    socR = reclassArray(climSoil, txtFile_socR)
    soc = socF * socR
    # Ensuring that cells with no stocks remain with no stocks 
    soc[LUmask.mask] = 0 
    # Summing total stocks in the array 
    tSOC = np.sum(soc)
    return soc, tSOC
    
def getTotalBCstock(climLu, LUmask, BC_comp):
    """Returns the total biomass carbon stocks for a given land use map (bc == array with total stock per cell, tBC == sum of array). The outputs depend whether the 'BC_comp' variable is set to run deterministically or stochastically. Arguments: climLU = array of climate + lu maps combination (for SOC factor); LUmask: it is the output of the maskNoCarbonCells() function. BC_comp = setup of the BC component (see parameters.configModelComponents())"""
    # Setting the path of txt files (depending if stochastic or not)
    if BC_comp == 1: # stochastic
        txtFile_bcs = os.path.join(txtFilesPath, 'bcs_{}.txt'.format(str(sample)))        
    if BC_comp == 0: # deterministic for sensitivity analysis
        txtFile_bcs = os.path.join('bcs_sensAnalysis.txt')
    if BC_comp == 'IPCC': # IPCC deterministic - this is only used in the IPCCdet_getOverallCarbonStock() function
        txtFile_bcs = os.path.join('bcs_deterministic.txt')
    # Allocating random values in arrays
    bc = (reclassArray(climLu, txtFile_bcs))
    # Ensuring sing that cells with no stocks remain with no stocks 
    bc[LUmask.mask] = 0
    # Summing total stocks in the array 
    tBC = np.sum(bc)
    return bc, tBC
        
def getTotalStock(tSOC, tBC):
    """Return the sum of SOC and biomass carbon stocks for a given land use map. The output is a single number representing the computed stocks for the study area. Arguments: tSOC = sum of SOC array (see getTotalSOCstock function); tBC = sum of BC array (see getTotalBCstock function"""     
    # Computing total stocks in the study area 
    tC = tSOC + tBC
    return tC

# Functions used to compute the difference of carbon stocks between scenarios (considering overall total carbon stocks, not cell-based)
def getDifferenceInStocksPerPairScenario(dictio_tSOCperMCr, dictio_tBCperMCr, dictio_tCperMCr):
    """Return a dictionary representing the difference of carbon stocks (tonne C/ha) between the scenario runs with and without ethanol increase (key == scenario, value == array w/ difference). Arguments: dictio_totStockPerMCr = dictionary with all scenarios (keys) and respective arrays (values) in which each element is a Monte Carlo sample representing the total carbon stocks in the study area"""
    stockDictiosInput = [dictio_tSOCperMCr, dictio_tBCperMCr, dictio_tCperMCr]    
    stockDictiosOutput = []
    
    # Differentiating the scenarios with and without additional ethanol production 
    scCode_initial, scCode_Eth, scCode_addEth = parameters.getScenariosType()
    # Getting the difference 
    for i, dictio in enumerate(stockDictiosInput):
        # Dividing input array in two dictionaries: 'Eth' and 'addEth' scenarios
        d_scEth = {k:np.asarray(v) for k,v in dictio.iteritems() if k in scCode_Eth}
        d_scAddEth = {k:np.asarray(v) for k,v in dictio.iteritems() if k in scCode_addEth}
        # Computing the subtraction: 'addEth' minus 'Eth' scenarios
        diff_totStockPerScPair = {k: d_scAddEth[k] - d_scEth.get(k-2) for k in d_scAddEth.keys()}
        stockDictiosOutput.append(diff_totStockPerScPair)
    SOCdiff = stockDictiosOutput[0]
    BCdiff = stockDictiosOutput[1]
    tCdiff = stockDictiosOutput[2]    
    return SOCdiff, BCdiff, tCdiff

def IPCCdet_getDiffStockPerPairScenarios():
    # Setting particle filter
    particleFilter_mapping = parameters.getMappingFromFile('particle_mapping.csv', 'New_ID', 'Nr_particles_weight')
    """Compute the carbon stocks deterministically, based on IPCC mean values""" 
    # Setting SOC and BC comps = 'IPCC' (necessary to get the correct txt files in getTotalSOCstock()/getTotalBCstock() functions)
    SOC_comp = 'IPCC'
    BC_comp = 'IPCC'
    # Creating dict to fill each scenario w/ total carbon stocks
    d_tSOC_det = {sc: [] for sc in scenariosDict.keys()}
    d_tBC_det = {sc: [] for sc in scenariosDict.keys()}
    d_tC_det = {sc: [] for sc in scenariosDict.keys()}        
    for sc, scPath in scenariosDict.iteritems():
        if sc == 0:
            LUmapPath = parameters.getInitialLandUseMap()
        if sc > 0:
            LUmapPath = parameters.getDeterministicLUmap(scPath, particleFilter_mapping)
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
    # Difference of total stock per MCr between paired scenarios (Eth Vs addEth). Output: dictio key is the code of SC w/additional Eth     
    diff_tSOCeth_det , diff_tBCeth_det, diff_tCeth_det = getDifferenceInStocksPerPairScenario(d_tSOC_det, d_tBC_det, d_tC_det)
    return diff_tSOCeth_det, diff_tBCeth_det , diff_tCeth_det 

# Function to get overal carbon stock stats (not cell-based)
def getStats_OverallCarbonStock(statsType, outFile=""):
    """Compute carbon stock statistics for SOC, BC and TC (with no conversion unit i.e.in tonnes C/ha). Arguments: statsType: the user has to type "absolute" OR diff_ethanol" OR "diff_ethanol_det". The 1st == stats for absolute stocks. The 2nd == stats of difference between scenarios with ethanol and with additional ethanol. The 3rd == stats to compute difference between scenarios and initial state; outFile: filename to save the statistics"""  
    if statsType == 'absolute':
        stockTypesList= ['SOC', 'BC', 'TC']
        stockDictiosList = [d_tSOCperMCr, d_tBCperMCr, d_tCperMCr]
    elif statsType == 'diff_ethanol':
        stockTypesList= ['d_SOC', 'd_BC', 'd_TC']
        stockDictiosList = [diff_tSOCeth, diff_tBCeth, diff_tCeth]
    elif statsType == 'diff_ethanol_det':
        stockTypesList= ['d_SOC', 'd_BC', 'd_TC']
        stockDictiosList = [diff_tSOCeth_det, diff_tBCeth_det, diff_tCeth_det]
    else:
        raise Exception("The 'statsType' argument is wrong. The user must set it as 'absolute' OR 'diff_ethanol' OR 'diff_ethanol_det'")
    # Creating empty dictionaries to add statistics later 
    statsToSave = {k: [] for k in stockDictiosList[0].keys()}
    meanDict = {k: [] for k in stockDictiosList[0].keys()}
    # Computing statistics
    for i, (stockType, dictio) in enumerate(zip(stockTypesList, stockDictiosList)):
        print '\n{}: carbon stocks + uncertainty for < {} > (tons C/ha):'.format(statsType, stockType)
        for sc, arr in dictio.iteritems():  # arr = array containing overall totStocks of each MCrun 
            mean = np.mean(arr)
            median = np.median(arr)
            std = np.std(arr)
            var = np.var(arr)
            percLow = np.percentile(arr, 2.5)
            percHigh = np.percentile(arr, 97.5)
            CIwidth = (percHigh - percLow)
            uncLow = ((mean - percLow) / mean) * 100
            uncHigh = ((percHigh - mean) / mean) * 100
            toSave = {'{}mean'.format(stockType): round(mean,2), '{}median'.format(stockType):round(median,2), '{}std'.format(stockType): round(std,2),'{}var'.format(stockType): round(var,2), '{}percLow'.format(stockType): round(percLow,2), '{}percHigh'.format(stockType): round(percHigh,2),'{}uncLow'.format(stockType): round(uncLow,2), '{}uncHigh'.format(stockType):round(uncHigh, 2)}
            meanToReturn = {'{}mean'.format(stockType): round(mean,2)}    
            if i == 0:
                statsToSave[sc] = toSave
                meanDict[sc] = meanToReturn
            else:
                for k, v in toSave.iteritems():
                    statsToSave[sc][k] = v
                for k, v in meanToReturn.iteritems():
                    meanDict[sc][k] = v
            print 'SC{} -\t|\tMean: {:.3f} ({:.2f}% , {:.2f}%)'.format(sc, mean, uncLow, uncHigh)
    # Converting stats to dataframe then saving to a file
    if statsType != 'diff_ethanol_det': # no need to save deterministic stats
        df = pd.DataFrame.from_dict(statsToSave, orient='columns', dtype=np.float64)
        df.to_csv(os.path.join(resultsDir, outFile), sep=';', decimal='.', header=True, index=True)
    return meanDict

# Functions to compute results in GHG emissions
def getEthanolDiffProduction(filename):
    """Return a dictionary with the difference in the production of ethanol between scenarios with ethanol production (eth) and with additional ethanol production (addEth), according to MAGNET output. Argument: filename = the file which has the MAGNET ethanol production per scenario"""  
    # Reading file with magnet ethanol production for all scenarios (in million liters)
    data = pd.read_csv(filename,  sep=",", header=0, usecols=['scenario','ethProd_mln_l'], index_col='scenario', skiprows=[1], dtype = np.int32)
    # Computing the difference between 'eth' and 'addEth'
    ethDiff = data.diff(axis=0)
    ethDiff = ethDiff[ethDiff.ethProd_mln_l > 0] # selecting only results of "addEth" minus "Eth" 
    # Converting results from 'mln l' to 'Gj'  
    ethanol_conversion = 23.4 * 1000 # mln_l to Gj    
    ethDiff['dEth_Gj'] = ethDiff['ethProd_mln_l']*ethanol_conversion  
    # Converting results to dictionary 
    ethDiffDict = dict(ethDiff['dEth_Gj'].round(0))
    #print "Additional ethanol production (million litres):\n", ethDiff['dEth_Gj']
    return ethDiffDict  

def getEmissionsDueToEthanol(dictio_diffSOC, dictio_diffBC, dictio_diffTC, ethDiff, outFile):
    """Return dictionaries of SOC, BC and TC arrays representing the GHG emissions per MC run due to the increase in ethanol production between the scenarios with eth (eth) and with additional eth (addEth) production. The final unit of measurement is == gram CO2-eq Mj-1 EtOH. Arguments: dictio_diffSOC, dictio_diffBC and dictio_diffTC: dictionaries where k = scenario and v = array w/ difference in carbon stocks per MCrun; ethDiff: output of getEthanolDiffProduction(); outFile: name to save in csv the results with the mean value of the arrays"""
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
        for (k1, array_carbonStock), (k2, deltaEth) in zip(dictio.iteritems() , ethDiff.iteritems()):
            if k1 == k2: # k1 and 2 are both scenario codes
                # Converting values in array to gram CO2-eq Mj-1 EtOH (our final measurement unit)
                array_hectares = array_carbonStock*parameters.getConversionUnit() # converting cStocks (tonne/ha) in total tonnes 
                array_EF = array_hectares*parameters.getEmissionFactorConversion() # converting cStocks to emission Factor
                array_20yrs = array_EF/20 # dividing by amortization period
                array_Eth = array_20yrs/deltaEth # dividing by the diff in ethanol production   
                array_FinalUnit = array_Eth*-1*1000 
                mean = np.mean(array_FinalUnit) 
                meanToDict = {'{}mean'.format(stockType): round(mean,2)}    
                if stockType == 'SOC':
                    SOC_ghgDueToEth[k1] = array_FinalUnit
                    ghgDueToEth_mean[k1] = meanToDict
                if stockType == 'BC':
                    BC_ghgDueToEth[k1] = array_FinalUnit 
                if stockType == 'TC':
                    TC_ghgDueToEth[k1] = array_FinalUnit
                for k, v in meanToDict.iteritems():
                    ghgDueToEth_mean[k1][k] = v   
    # Saving means to csv file 
    data = pd.DataFrame.from_dict(ghgDueToEth_mean, orient='index') 
    data = data.reindex(sorted(data.columns), axis=1)
    data.to_csv(os.path.join(resultsDir, outFile), sep=';', decimal='.', header=True, index=True)
    print '\n{}:\tEmissions due to increase in ethanol production (CO2-eq Mj-1 EtOH):\n\n'.format(outFile), data 
    return SOC_ghgDueToEth, BC_ghgDueToEth, TC_ghgDueToEth

def saveData4statisticalTest(outFile):
    """Save SOC, BC and TC arrays of GHG emissions (CO2-eq Mj-1 EtOH) for statistical test in R"""
    # Creating lists to add outputs
    soc4statsTest, bc4statsTest, tc4statsTest = [], [], []
    data4StatsTestList = [soc4statsTest, bc4statsTest, tc4statsTest]
    # Start processing
    stockTypesList = ['SOC', 'BC', 'TC']
    stockDictiosList = [socGHGe, bcGHGe, tcGHGe] 
    i=0
    for dictio in stockDictiosList:
        for array in dictio.values():
            data4StatsTestList[i].append(array)
        i+=1
    for stockType, data4StatsTest in zip(stockTypesList, data4StatsTestList):
        df = pd.DataFrame(data4StatsTest)
        df.to_csv(os.path.join(resultsDir, '{}_{}.txt'.format(outFile, stockType)), sep=';', decimal='.', header=False, index=False)
# Functions to plot
def setupBoxplot(boxplot, color): 
    """Customize the boxplot built by  showBoxplot() function"""
    plt.setp(boxplot['boxes'], facecolor=color, edgecolor='#585858', linewidth=0.4, alpha=0.8)
    plt.setp(boxplot['whiskers'], color='#848484', linewidth=0.6, linestyle="--")
    plt.setp(boxplot['caps'], color='#848484', linewidth=1)
    plt.setp(boxplot['medians'], linewidth=1 ,color ='#B40404', linestyle="-")
    plt.setp(boxplot['means'], linestyle="-",color='k', linewidth=1.6, marker="+", markersize=6)

def createBoxplot(): # OBS: I still couldn't figure out how to use a marker instead of a line in the "usermedians" argument of boxplot, that's why legend and boxplot don't match for "deterministic"
    """Create a boxplot of GHG emissions due to ethanol increase between scenarios with eth ('eth')and with additional eth ('addEth')"""
    if plotFigures == 1:
        plt.show()     
        fig, ax = plt.subplots(figsize=(8, 5))
        colors = parameters.figureColors()
        socBoxes = plt.boxplot(socGHGe.values(), positions=np.array(xrange(len(socGHGe.values())))*2+0.4, sym='', usermedians=socGHGe_det.values(), meanline=True, showmeans=True, patch_artist=True, widths=0.6)
        bcBoxes = plt.boxplot(bcGHGe.values(), positions=np.array(xrange(len(bcGHGe.values())))*2-0.4, sym='', usermedians=bcGHGe_det.values(), meanline=True, showmeans=True, patch_artist=True, widths=0.6)
        setupBoxplot(bcBoxes, colors[1])
        setupBoxplot(socBoxes, colors[0])
        ymin, ymax = (int(np.max(np.stack((socGHGe.values(), bcGHGe.values())))*-1), int(np.max(np.stack((socGHGe.values(), bcGHGe.values())))))
        print ymin, ymax
        plt.yticks(np.arange(ymin, ymax+1, step=5))
        plt.xticks(np.arange(0, len(socGHGe.values())*2, 2), parameters.getScenariosNames(), fontsize=8.5)
        plt.xlim(-1, len(socGHGe.values())*2-1)
        # draw bar/plots to use as legend
        plt.bar(1,[0], color=colors[0], label='SOC', alpha=0.80)
        plt.bar(1,[0], color=colors[1], label='Biomass', alpha=0.80)
        plt.plot([], color='k', linewidth=2, marker="+", markersize=6, label='Mean - stochastic') 
        plt.plot([],[],  marker="^", markerfacecolor='w', markeredgecolor = "#c10202", markersize=8, color ="w",label="Mean - deterministic")
        plt.legend(fontsize=7, loc = 'upper adj', frameon=True)
        plt.grid(which="minor", axis='x',color='0.3', linewidth=0.15)
        plt.title('GHG emissions from SOC/BC resulting \nfrom an increased ethanol production\n',  color='#6E6E6E', fontsize=10, fontweight='bold', pad=2)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        plt.ylabel('GHG emissions\n'+r'$gram$ $CO_2$-$eq$/$Mj$$_E$$_t$$_O$$_H$')
        if savePlots == 1:
            plt.savefig(os.path.join(resultsDir, 'boxplot'), dpi=700)
        if showPlots == 1:
            plt.show()    

#-------------------     GENERATING RANDOM VALUES     -------------------#

txtFilesPath = os.path.join('txt_files')
if generateRandomValues == 1:
    # Setting random seed 
    randomSeed = np.random.seed(seed=62871) # Obs: I added random seed after the thesis. If you need the thesis txtFiles, see folder 'txt_files_thesis' 
    print 'Generating random values per MC...'
    generateRandomSOCref(nrMCsamples,'socR_uncertainty.txt', txtFilesPath)
    generateRandomSOCfactor(nrMCsamples,'socF_uncertainty.txt', txtFilesPath)
    generateRandomBCstock(nrMCsamples,'bcs_uncertainty.txt', txtFilesPath)
    print '\nChecking occurrence of negative values...'
    checkNegRandomValues(txtFilesPath)

#-------------------     DATA PROCESSING (total carbon stocks)     -------------------#

# Input rasters
climate = rasterToArray('climate')
clim_soil = rasterToArray('climsoil')   

if getOverallCarbonStock == 1:
    # Setting particle filter mapping
    particleFilter_mapping = parameters.getMappingFromFile('particle_mapping.csv', 'New_ID', 'Nr_particles_weight')
    # Creating dict w/ lists for each scenario >>> they will contain the sum of array elements of tSOC, tBC & tC, for each MC run
    d_tSOCperMCr = {sc: [] for sc in scenariosDict.keys()}
    d_tBCperMCr = {sc: [] for sc in scenariosDict.keys()}
    d_tCperMCr = {sc: [] for sc in scenariosDict.keys()}
    d_CellsPerLU = {sc:'' for sc in scenariosDict.keys()} ### <=============== Here is the dictionary which I compute nr of Cells per LU
    # Running stochastic 
    for sc, scPath in scenariosDict.iteritems():
        print '\nSC{} - computing carbon stocks (tonne C/ha)\t...\tStart time: {}'.format(sc, time.asctime(time.localtime(time.time())))
        if LUC_comp == 0 or sc == 0:
            LUmapPath = getLUmap(sc, scPath)
            LUmap = rasterToArray(LUmapPath)
            mask = maskNoCarbonCells(LUmap)
            print '\n\tLUmap (deterministic):\tsc{}/{}\n'.format(sc, LUmapPath.split("/")[-1])  
            d_CellsPerLU[sc] = getNrCellsPerLUtype(LUmap) #### <=============== # Here is where I call the function to compute Nr of Cells, according to the LU map that is got by getLUmap()
            print 'sc{}, path: {}...'.format(sc, LUmapPath),"dictio w/ cells per LU type:\n\t", d_CellsPerLU[sc] ### <=========== print of the results
            clim_lu = (climate + LUmap)  # np.int32 (0-511)
            for sample in range(1, nrMCsamples+1):
                # Computing total stocks in the study area
                soc, tSOC = getTotalSOCstock(clim_lu, clim_soil, mask, SOC_comp)
                bc, tBC = getTotalBCstock(clim_lu, mask, BC_comp)
                tC = getTotalStock(tSOC, tBC)
                # Appending totals to list
                d_tSOCperMCr[sc].append(tSOC)
                d_tBCperMCr[sc].append(tBC)
                d_tCperMCr[sc].append(tC)
                #print '\tSC{} MC{}\t|\t{}\t|\ttSOC: {}\ttBCS: {}\ttCS: {} (ton C/ha)'.format(sc, sample, modelRunType, tSOC, tBC, tC)
        if LUC_comp == 1 and sc > 0:
            mcRun_acc = 0 # iterator to count nr of mcRuns (necessary to change the LU map according to particle filter 
            LUmapsDict = parameters.getStochasticLUmaps(scPath, particleFilter_mapping, '2030_*')
            for mapDescr in sorted(LUmapsDict.items()):
                mapCode, mapPath, mapMCruns = mapDescr[0], mapDescr[1][0], mapDescr[1][1]
                LUmapPath = getLUmap(sc, scPath, mapPath)
                LUmap = rasterToArray(LUmapPath)
                clim_lu = (climate + LUmap)
                mask = maskNoCarbonCells(LUmap)
                mcRange_curr  = getLUmapMCrange(mapCode, mapMCruns, mcRun_acc) # the MC range used to run the current LU map
                mcRun_acc += mapMCruns  
                print '\n\tLUmap nr {} (stochastic): {}\t|\tsamples: {} (runs {}-{})\n'.format(mapCode, LUmapPath.split("ts/")[-1], mapMCruns, min(mcRange_curr), max(mcRange_curr))
                for sample in mcRange_curr:
                    # Computing total stocks in the study area
                    soc, tSOC = getTotalSOCstock(clim_lu, clim_soil, mask, SOC_comp)
                    bc, tBC = getTotalBCstock(clim_lu, mask, BC_comp)
                    tC = getTotalStock(tSOC, tBC)
                    # Appending totals to list
                    d_tSOCperMCr[sc].append(tSOC)
                    d_tBCperMCr[sc].append(tBC)
                    d_tCperMCr[sc].append(tC)
                    if saveArrays == 1 and modelRunType == 'Full stochastic':
                        'cbStocks_sc{}mc{}'.format(sc,i+1)
                        np.savez_compressed(os.path.join(resultsDir, 'arr_cellBased', 'cb_sc{}mc{}'.format(sc,i+1)), soc=soc, bc=bc) # next: add the map array thing mentioned by Judith to reduce disk storage
                    print '\tSC{} MC{}\t|\t{}\t|\ttSOC: {}\ttBCS: {}\ttCS: {} (ton C/ha)'.format(sc, sample, modelRunType, tSOC, tBC, tC)
        if saveArrays == 1:
            # Saving list w/ total stock per MCr as array
            np.save(os.path.join(resultsDir, 'arr_totalStocks', "{}_sc{}_tSOC".format(modelRunType, sc)), np.asarray(d_tSOCperMCr[sc]))
            np.save(os.path.join(resultsDir, 'arr_totalStocks', "{}_sc{}_tBC".format(modelRunType, sc)), np.asarray(d_tBCperMCr[sc]))
            np.save(os.path.join(resultsDir, 'arr_totalStocks', "{}_sc{}_tC".format(modelRunType, sc)), np.asarray(d_tCperMCr[sc]))
    
    ### Saving dictionary with cells per LU type (meaning that it will save the result of the LUC map with more weight (i.e., nr 43)
    if LUC_comp == 0:
        df = pd.DataFrame(d_CellsPerLU).T
        df.to_csv(os.path.join(resultsDir, "d_CellsPerLU_det.txt"), sep=',')    


if getOverallCarbonStock == 0:
    # Opening numpy arrays with total stocks of each MC sample
    d_tSOCperMCr = {sc:np.load(os.path.join(resultsDir, 'arr_totalStocks', '{}_sc{}_tSOC.npy'.format(modelRunType, sc))) for sc in scenariosDict.keys()} 
    d_tBCperMCr = {sc:np.load(os.path.join(resultsDir, 'arr_totalStocks', '{}_sc{}_tBC.npy'.format(modelRunType, sc))) for sc in scenariosDict.keys()} 
    d_tCperMCr = {sc:np.load(os.path.join(resultsDir, 'arr_totalStocks', '{}_sc{}_tC.npy'.format(modelRunType, sc))) for sc in scenariosDict.keys()} 


# Difference of total stock per MCr between paired scenarios (Eth Vs addEth). Output: dictio key is the code of SC w/additional Eth     
diff_tSOCeth, diff_tBCeth, diff_tCeth = getDifferenceInStocksPerPairScenario(d_tSOCperMCr, d_tBCperMCr, d_tCperMCr)

# Getting the difference but now running the model deterministically (based on IPCC mean values)
diff_tSOCeth_det , diff_tBCeth_det , diff_tCeth_det = IPCCdet_getDiffStockPerPairScenarios()

#-------------------     OBTAINING STATISTICS & SAVING     -------------------#

scStats = getStats_OverallCarbonStock('absolute', '{}_absolute_Stats.csv'.format(modelRunType))
scDiffStats = getStats_OverallCarbonStock('diff_ethanol', "{}_dEth_Stats.csv".format(modelRunType))
scDiffStats_det = getStats_OverallCarbonStock('diff_ethanol_det')

#-------------------     COMPUTING RESULTS IN GHG EMISSIONS BASED ON ETHANOL PRODUCTION     -------------------#

# Getting emissions due to ethanol production (for SOC, BC and TC) (1st line == stochastic, 2nd line == deterministic)
ethDiff = getEthanolDiffProduction('magnet_output.txt')
socGHGe, bcGHGe, tcGHGe = getEmissionsDueToEthanol(diff_tSOCeth, diff_tBCeth, diff_tCeth, ethDiff, '{}_GHGe.csv'.format(modelRunType))
socGHGe_det, bcGHGe_det, tcGHGe_det = getEmissionsDueToEthanol(diff_tSOCeth_det, diff_tBCeth_det, diff_tCeth_det, ethDiff, 'determ_GHGe.csv')
# Saving socGHGe, bcGHGe, tcGHGe variables for statistical test in R
if modelRunType == "stcAll":
    saveData4statisticalTest('data4StatsTest')

#-------------------     BOXPLOT     -------------------#

createBoxplot()

print ("\nDone.\n\n---\t{0:.2f} sec ({1:.2f} min)\t---".format(float(time.time() - start), float((time.time() - start)/60)))

