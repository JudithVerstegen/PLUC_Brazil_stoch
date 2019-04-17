import time
start = time.time()

print ('Starting ...\n')

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
getOverallCarbonStock = 0
getCellBasedCarbonStock = 1
plotFigures = 0
showPlots = 0
savePlots = 0
saveArrays = 0
runSensAnalysis = 0

print ('Nr of scenarios (including initial): {}.\nMC simulation: {} samples per scenario.'.format(len(scenariosDict), nrMCsamples))

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
                    print ('Neg/zero value identified for SOCr...Getting positive value ...')
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
                        print ('Neg. value identified for SOC factor(s). Getting positive value(s)...')
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
                    print ('Neg. value identified for AGB. Getting positive value ...')
                    agb = np.random.normal(agb_mean, agb_std)
                r2s = np.random.normal(r2s_mean, r2s_std)
                while r2s < 0:
                    print ('Neg. value identified for r2s. Getting positive value by truncation ...')
                    if r2s_std != 0:  # scipy.stats does not allow computation of std = 0
                        a, b = (r2s_min - r2s_mean) / r2s_std, (r2s_max - r2s_mean) / r2s_std  # setting truncation
                        r2s = ss.truncnorm.rvs(a, b, r2s_mean, r2s_std)
                    else:
                        r2s = np.random.normal(r2s_mean, r2s_std)
                cf = np.random.normal(cf_mean, cf_std)
                while cf < 0:
                    print ('Neg. value identified for CF. Getting positive value ...')
                    cf = np.random.normal(cf_mean, cf_std)
                bgb = agb * r2s
                BCS = (agb + bgb) * cf
                if BCS == 0:  # if clause necessary to compute BCS for CROPs and ABANDONNED LAND (IPCC doesn't give AGB, R2S and CF for those LU types)
                    BCS = (np.random.normal(BCS_mean, BCS_std))
                    while BCS < 0:
                        print ('Neg. value identified for BCS. Getting positive value ...')
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
            print (txtFile.split('txt_files/')[1], '\n', df[df.value < 0])
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
def getDifferenceInStocksPerPairScenario(dictio_totStockPerMCr):
    """Return a dictionary representing the difference of carbon stocks (tonne C/ha) between the scenario runs with and without ethanol increase (key == scenario, value == array w/ difference). Arguments: dictio_totStockPerMCr = dictionary with all scenarios (keys) and respective arrays (values) in which each element is a Monte Carlo sample representing the total carbon stocks in the study area"""
    # Differentiating the scenarios with and without additional ethanol production 
    scCode_initial, scCode_Eth, scCode_addEth = parameters.getScenariosType()
    # Dividing input array in two dictionaries: 'Eth' and 'addEth' scenarios
    d_scEth = {k:np.asarray(v) for k,v in dictio_totStockPerMCr.iteritems() if k in scCode_Eth}
    d_scAddEth = {k:np.asarray(v) for k,v in dictio_totStockPerMCr.iteritems() if k in scCode_addEth}
    # Computing the subtraction: 'addEth' minus 'Eth' scenarios
    diff_totStockPerScPair = {k: d_scAddEth[k] - d_scEth.get(k-2) for k in d_scAddEth.keys()}
    return diff_totStockPerScPair

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
            print (LUmapPath)
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
    diff_tSOCeth_det = getDifferenceInStocksPerPairScenario(d_tSOC_det)
    diff_tBCeth_det = getDifferenceInStocksPerPairScenario(d_tBC_det)
    diff_tCeth_det = getDifferenceInStocksPerPairScenario(d_tC_det)
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
        print ('\n{}: carbon stocks + uncertainty for < {} > (tons C/ha):'.format(statsType, stockType))
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
            print ('SC{} -\t|\tMean: {:.3f} ({:.2f}% , {:.2f}%)'.format(sc, mean, uncLow, uncHigh))
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
    #print ("Additional ethanol production (million litres):\n", ethDiff['dEth_Gj'])
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
    print ('\n{}:\tEmissions due to increase in ethanol production (CO2-eq Mj-1 EtOH):\n\n'.format(outFile), data) 
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
# Functions for boxplot
def setupBoxplot(boxplot, color): 
    """Customize the boxplot built by plotBoxplot() function"""
    plt.setp(boxplot['boxes'], facecolor=color, edgecolor='#585858', linewidth=0.4, alpha=0.8)
    plt.setp(boxplot['whiskers'], color='#848484', linewidth=0.6, linestyle="--")
    plt.setp(boxplot['caps'], color='#848484', linewidth=1)
    plt.setp(boxplot['medians'], linewidth=1 ,color ='#B40404', linestyle="-")
    plt.setp(boxplot['means'], linestyle="-",color='k', linewidth=1.6, marker="+", markersize=6)

def plotBoxplot(): # OBS: I still couldn't figure out how to use a marker instead of a line in the "usermedians" argument of boxplot, that's why legend and boxplot don't match for "deterministic"
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
            plt.savefig(os.path.join(resultsDir, 'plot_boxplot'), dpi=700)
        if showPlots == 1:
            plt.show()    

# Function to run sensitivity analysis 
def plotSensitivityAnalysis():
    """Run the sensitivity analysis based on three model components, then plot it (stcSOC', 'stcBC', 'stcLUC')"""
    if runSensAnalysis == 1:
        # Loading total CS per MCrun, getting variance and appending to dictios
        components = ['stcSOC', 'stcBC', 'stcLUC', 'stcAll']
        # Creating empty dictionary to append variances per scenario, per component
        varDict = {k:"" for k in components}
        for comp in components:
            try: 
                # Loading array per SC w/ all MC runs results (representing carbon stocks in tons C/ha)
                path = os.path.join(resultsDir, 'arr_totalStocks')
                totalStocks = {sc:np.load(os.path.join(path, "{}_sc{}_tC.npy".format(comp, sc))) for sc in scenariosDict.keys()}
            except:
                raise Exception("'{}' files not found in {}. Before the sensitivity analysis, you have to consider running the model with the '{}' component set to run stochastically and other components to run deterministically. To save the arrays of this run, the 'saveArrays' variable must be == 1".format(comp, path, comp))
            # Getting the array of the difference between scenarios with eth ('eth') and with additional eth ('addEth') 
            diffStocks = getDifferenceInStocksPerPairScenario(totalStocks)
            # Getting diffStocks variance per scenario, for each component
            variancePerSc = {k:{} for k in diffStocks.keys()}
            for sc in variancePerSc.keys():
                var = np.var(diffStocks[sc]).astype(np.float64)
                variancePerSc[sc] = var
            # Adding VariancePerSc to varDict
            varDict[comp] = variancePerSc
        varData = pd.DataFrame.from_dict(varDict, orient='columns', dtype=np.float64)
        # Sensitivity analysis: obtaining variance fractions per component based on the total variance of "stcAll" (i.e, full stochastic run of the model) 
        varFract = (varData.div(varData['stcAll'], axis=0))*100
        # Getting the total contribution of all components in overall variance    
        varFract['tot'] = varFract['stcBC'] + varFract['stcSOC'] + varFract['stcLUC']
        # Getting the model interactions contribution in overall variance 
        varFract['diff'] = varFract['stcAll'] - varFract['tot']
            # Plotting resultsplot_sensAnalysis
        fig, ax = plt.subplots(figsize=(10, 6))
        colors = ['#5ab4ac', '#d8b365', '#f6e8c3','#610B0B']  
        bars = np.arange(6)
        barWidth = 0.75
        labels = ['Reference\nscenario', 'High\nproductivity.', r'2$^n$$^d$ generation'+'\nsugar cane', r'2$^n$$^d$ generation'+'\neucalyptus', 'Conservation\npolicies', 'All\nmeasures'] 
        # Plotting bars  
        barBC = plt.bar(bars, varFract.stcBC, color=colors[0], width=barWidth)
        barSOC = plt.bar(bars, varFract.stcSOC, bottom=varFract.stcBC, color=colors[1], width=barWidth)
        barLUC = plt.bar(bars, varFract.stcLUC, bottom=varFract.stcBC+varFract.stcSOC, color=colors[2], width=barWidth)
        barOther = plt.bar(bars, varFract['diff'], bottom=varFract['tot'], color=colors[3], width=barWidth)
        # Plot setup
        plt.xticks(bars, labels, fontsize=10)
        plt.yticks(np.arange(0, 101, step=50), fontsize=9.5)
        plt.ylabel("Component contribution in the variance of GHGE from ethanol production/\nTotal variance of GHGE from ethanol production (%)*", fontsize=9)
        plt.legend((barBC[0], barSOC[0], barLUC[0], barOther[0]), ('Biomass', 'SOC', 'LUC', 'Model interactions'), bbox_to_anchor=(0.78,-0.12), fontsize=10, ncol=4)
        plt.grid(which="major", axis="y", color='0.7', linewidth=0.3)
        plt.subplots_adjust(bottom=0.2)
        if savePlots == 1:
            plt.savefig(os.path.join(resultsDir, 'plot_sensAnalysis'), dpi=700)
        if showPlots == 1:
            plt.show()
    





# Function to transform array to map (tiff extension)
def arrayToMap(inpArray, outRaster):
    """Convert an array to a raster map and save to disk (tiff format). Arguments: inpArray = array to be converted; outRaster = string with the filename representing the raster map"""
    # Getting geoparameters of the original raster
    refRaster = gdal.Open(parameters.referenceRaster)
    geotransform = refRaster.GetGeoTransform()
    originX = geotransform[0]  # top left-x
    originY = geotransform[3]  # top left-y
    pixelWidth = geotransform[1]  # w-e pixel resolution
    pixelHeight = geotransform[5]  # n-s pixel resolution
    rows = refRaster.RasterYSize
    cols = refRaster.RasterXSize
    # Generating a masked array representing cells with null value (further used when converting array to raster)
    nvArr = (refRaster.ReadAsArray()).astype(np.int32)  # reference-array for masking null value
    nvArr[(nvArr < 0) | (nvArr > nrLUtypes)] = -9999  # Converts to -9999 any value that doesn't represent a LU type code
    nvMsk = np.ma.masked_where(nvArr != -9999, nvArr)  # null value mask
    inpArray[~nvMsk.mask] = -9999  # Updating inpArray with the null value of the MskArr
    # Creating new raster, then adding geoparameters
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(os.path.join(resultsDir, "maps", outRaster + ".tiff"), cols, rows, 1, gdal.GDT_Float64)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(inpArray)
    outband.SetNoDataValue(-9999)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(refRaster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

#-------------------     GENERATING RANDOM VALUES     -------------------#

txtFilesPath = os.path.join('txt_files')
if generateRandomValues == 1:
    # Setting random seed 
    randomSeed = np.random.seed(seed=62871) # Obs: I added random seed after the thesis. If you need the thesis txtFiles, see folder 'txt_files_thesis' 
    print ('Generating random values per MC...')
    generateRandomSOCref(nrMCsamples,'socR_uncertainty.txt', txtFilesPath)
    generateRandomSOCfactor(nrMCsamples,'socF_uncertainty.txt', txtFilesPath)
    generateRandomBCstock(nrMCsamples,'bcs_uncertainty.txt', txtFilesPath)
    print ('\nChecking occurrence of negative values...')
    checkNegRandomValues(txtFilesPath)

#-------------------     DATA PROCESSING (total carbon stocks)     -------------------#

# Input rasters
climate = rasterToArray('climate')
clim_soil = rasterToArray('climsoil')   

# Setting particle filter mapping
particleFilter_mapping = parameters.getMappingFromFile('particle_mapping.csv', 'New_ID', 'Nr_particles_weight')

if getOverallCarbonStock == 1:
    # Creating dict w/ lists for each scenario >>> they will contain the sum of array elements of tSOC, tBC & tC, for each MC run
    d_tSOCperMCr = {sc: [] for sc in scenariosDict.keys()}
    d_tBCperMCr = {sc: [] for sc in scenariosDict.keys()}
    d_tCperMCr = {sc: [] for sc in scenariosDict.keys()}
    d_CellsPerLU = {sc:'' for sc in scenariosDict.keys()} ### <=============== Here is the dictionary which I compute nr of Cells per LU
    # Running stochastic 
    for sc, scPath in scenariosDict.iteritems():
        print ('\nSC{} - computing carbon stocks (tonne C/ha)\t...\tStart time: {}'.format(sc, time.asctime(time.localtime(time.time()))))
        if LUC_comp == 0 or sc == 0:
            LUmapPath = getLUmap(sc, scPath)
            LUmap = rasterToArray(LUmapPath)
            mask = maskNoCarbonCells(LUmap)
            print ('\n\tLUmap (deterministic):\tsc{}/{}\n'.format(sc, LUmapPath.split("/")[-1]))
            d_CellsPerLU[sc] = getNrCellsPerLUtype(LUmap) #### <=============== # Here is where I call the function to compute Nr of Cells, according to the LU map that is got by getLUmap()
            print ('sc{}, path: {}...'.format(sc, LUmapPath),"dictio w/ cells per LU type:\n\t", d_CellsPerLU[sc]) ### <=========== print of the results
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
                #print ('\tSC{} MC{}\t|\t{}\t|\ttSOC: {}\ttBCS: {}\ttCS: {} (ton C/ha)'.format(sc, sample, modelRunType, tSOC, tBC, tC))
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
                print ('\n\tLUmap nr {} (stochastic): {}\t|\tsamples: {} (runs {}-{})\n'.format(mapCode, LUmapPath.split("ts/")[-1], mapMCruns, min(mcRange_curr), max(mcRange_curr)))
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
                        np.savez_compressed(os.path.join(resultsDir, 'arr_cellBased', 'cb_sc{}mc{}'.format(sc,sample)), soc=soc, bc=bc) # next: add the map array thing mentioned by Judith to reduce disk storage
                    print ('\tSC{} MC{}\t|\t{}\t|\ttSOC: {}\ttBCS: {}\ttCS: {} (ton C/ha)'.format(sc, sample, modelRunType, tSOC, tBC, tC))
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
diff_tSOCeth = getDifferenceInStocksPerPairScenario(d_tSOCperMCr)
diff_tBCeth = getDifferenceInStocksPerPairScenario(d_tBCperMCr)
diff_tCeth = getDifferenceInStocksPerPairScenario(d_tCperMCr)

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

#-------------------     GETTING OUTPUTS (BOXPLOT, SENSITIVITY ANALYSIS, DATA FOR STATISTICAL TEST      -------------------#

# Boxplot
plotBoxplot()

# For sensAnalysis, the model must run four times (see parameters.configModelComponents()): "full stochastic", "stcSOC", "stcBC" and "stcAll".   
plotSensitivityAnalysis()

# Saving socGHGe, bcGHGe, tcGHGe variables for statistical test in R (only if the model setup runs full stochastically (i.e., 'stcAll')
if modelRunType == "stcAll":
    saveData4statisticalTest('data4StatsTest')

print ("\nDone.\n\n---\t{0:.2f} sec ({1:.2f} min)\t---".format(float(time.time() - start), float((time.time() - start)/60)))

#-------------------     CELL-BASED FUNCTIONS (NEED TO ORGANIZE)    -------------------#
# Functions used to get cell-based carbon stocks with uncertainty
def getCellMean(csType, sc, saveTo):
    """XXXXXX"""
    npzFiles = sorted(glob.glob(os.path.join(resultsDir, 'arr_cellBased', 'cbStocks_sc{}mc*'.format(sc))))
    npzFilesDict = {int((v.split('mc'))[-1].split('.')[0]): v for v in sorted(npzFiles)}
    # Computing total cell-based mean and per LU type 
    cs_accum = createArr4cellBasedStocks() 
    cs_sqAccum = createArr4cellBasedStocks()
    if sc == 0:
        LUmapPath = parameters.getInitialLandUseMap()
        LUmap = rasterToArray(LUmapPath)
        mask = maskNoCarbonCells(LUmap)
        for sample, fl in npzFilesDict.iteritems():
            npData = np.load(fl)
            cs = npData[csType]
            cs[mask.mask] = 0 
            cs_accum = cs + cs_accum
            #print ('Processing file {} ...'.format(npzFilesDict[sample].split("_")[-1]))
            if np.isnan(cs_accum).any() == True:
                print ('SC{} - NaN found in cs_accum when summing {} with file {}'.format(sc, csType, npzFilesDict[sample]))
    if sc > 0:
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
            #print ('runs {}-{}'.format(min(mcRange_curr), max(mcRange_curr)))
            for sample in mcRange_curr:
                npData = np.load(npzFilesDict[sample])
                cs = npData[csType]
                cs[mask.mask] = 0 
                cs_accum = cs + cs_accum
                #print ('Processing file {} ...'.format(npzFilesDict[sample].split("_")[-1]))
                if np.isnan(cs_accum).any() == True:
                    print ('SC{} - NaN found in cs_accum when summing {} with file {}'.format(sc, csType, npzFilesDict[sample]))
        csMean_cell = (cs_accum/len(npzFiles)) # OK! (It has 0 values)
        # Saving array and transforming to map
        np.save(os.path.join(results, "sc{}_{}_mean".format(sc,csType)), csMean_cell)
        arrayToMap(csMean_cell, "sc{}_{}_mean".format(sc,csType))
        return csMean_cell
      
        
        
        for mapDescr in sorted(LUmapsDict.items()):
            mapCode, mapPath, mapMCruns = mapDescr[0], mapDescr[1][0], mapDescr[1][1]
            LUmapPath = getLUmap(sc, scPath, mapPath)
            LUmap = rasterToArray(LUmapPath)
            clim_lu = (climate + LUmap)
            mask = maskNoCarbonCells(LUmap)
            mcRange_curr  = getLUmapMCrange(mapCode, mapMCruns, mcRun_acc) # the MC range used to run the current LU map
            mcRun_acc += mapMCruns  
            #print ('runs {}-{}'.format(min(mcRange_curr), max(mcRange_curr))
            for sample in mcRange_curr:
                npData = np.load(npzFilesDict[sample])
                cs = npData[csType]
                cs[mask.mask] = 0 
                cs_accum = cs + cs_accum
                # Computation Per LU
                print ('Processing file {} ...'.format(npzFilesDict[sample].split("_")[-1]))
                if np.isnan(cs_accum).any() == True:
                    print ('SC{} - NaN found in cs_accum when summing {} with file {}'.format(sc, csType, npzFilesDict[sample]))
        
        
                #for mcRun in npFiles:
                #npData = np.load(mcRun)
                #cs = npData[csType]
                #csMask = (np.ma.masked_values(cs,0)) # mask to remove zero values when doing the computation of squares
                #csSq = (csMask - csMean_cell)**2
                #cs_sqAccum = csSq + cs_sqAccum
                ## Computation Per LU
                #csPerLU = getStockPerLUtype(cs, lu)   
                #csPerLU_mask = (np.ma.masked_values(csPerLU,0)) # mask to remove zero values when doing the computation of squares
                #csPerLU_sq = (csPerLU_mask - csMean_cell_PerLU)**2
                #csPerLU_sqAccum = csPerLU_sq + csPerLU_sqAccum
            #csStd_cell = np.sqrt(cs_sqAccum/len(npFiles))
            #csStd_cell = np.asarray(csStd_cell)
            #csStd_cell[luMask.mask] = 0 
            #csStd_cell_PerLU = np.sqrt(csPerLU_sqAccum/len(npFiles))
            #csStd_cell_PerLU = np.asarray(csStd_cell_PerLU)
            #np.save(os.path.join(saveTo, "sc{}_m{}perCell".format(sc,csType)), csMean_cell.astype(np.float32))
            #np.save(os.path.join(saveTo, "sc{}_std{}perCell".format(sc,csType)), csStd_cell.astype(np.float32))
            ## Savin CS per LU
            #if os.path.exists(os.path.join(outDir, "cs_{}perLU_mean.txt".format(csType))) == True and os.path.exists(os.path.join(outDir, "cs_{}perLU_std.txt".format(csType))) == True:
                #with open((os.path.join(outDir, "cs_{}perLU_mean.txt".format(csType))),"a+") as mean:
                    #mean.write('sc{};'.format(sc))           
                    #np.savetxt(mean,csMean_cell_PerLU.T, delimiter=';')    
                #with open((os.path.join(outDir, "cs_{}perLU_std.txt".format(csType))),"a+") as std:
                    #std.write('sc{};'.format(sc))           
                    #np.savetxt(std,csStd_cell_PerLU.T, delimiter=';')      
            #else:
                #with open((os.path.join(outDir, "cs_{}perLU_mean.txt".format(csType))),"w+") as mean:
                    #mean.write('MEAN stock per LU\nsc{};'.format(sc))        
                    #np.savetxt(mean,csMean_cell_PerLU.T, delimiter=';')    
                #with open((os.path.join(outDir, "cs_{}perLU_std.txt".format(csType))),"w+") as std:
                    #std.write('STD stock per LU\nsc{};'.format(sc))      
                    #np.savetxt(std,csStd_cell_PerLU.T, delimiter=';')              
            
            
              
              
                

            #for sample in mcRange_curr:         
            #.return csMean_cell, csStd_cell


def getCellStdDev(csType, sc, cellMean, saveTo):
    """XXXXXX"""
    npzFiles = sorted(glob.glob(os.path.join(resultsDir, 'arr_cellBased', 'cbStocks_sc{}mc*'.format(sc))))
    npzFilesDict = {int((v.split('mc'))[-1].split('.')[0]): v for v in sorted(npzFiles)}
    # Computing total cell-based mean and per LU type 
    cs_accum = createArr4cellBasedStocks() 
    cs_sqAccum = createArr4cellBasedStocks()
    if sc == 0:
        LUmapPath = parameters.getInitialLandUseMap()
        LUmap = rasterToArray(LUmapPath)
        mask = maskNoCarbonCells(LUmap)
        for sample, fl in npzFilesDict.iteritems():
            npData = np.load(fl)
            cs = npData[csType]
            cs[mask.mask] = 0 
            cs_accum = cs + cs_accum
            print ('Processing file {} ...'.format(npzFilesDict[sample].split("_")[-1]))
            if np.isnan(cs_accum).any() == True:
                print ('SC{} - NaN found in cs_accum when summing {} with file {}'.format(sc, csType, npzFilesDict[sample]))
    if sc > 0:
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
            #print ('runs {}-{}'.format(min(mcRange_curr), max(mcRange_curr))
            for sample in mcRange_curr:
                npData = np.load(npzFilesDict[sample])
                cs = npData[csType]
                cs[mask.mask] = 0 
                cs_accum = cs + cs_accum
                #print ('Processing file {} ...'.format(npzFilesDict[sample].split("_")[-1])
                if np.isnan(cs_accum).any() == True:
                    print ('SC{} - NaN found in cs_accum when summing {} with file {}'.format(sc, csType, npzFilesDict[sample]))
        csMean_cell = (cs_accum/len(npzFiles)) # OK! (It has 0 values)
        # Saving array and transforming to map
        np.save(os.path.join(results, "sc{}_{}_mean".format(sc,csType)), csMean_cell)
        arrayToMap(csMean_cell, "sc{}_{}_mean".format(sc,csType))
        return csMean_cell




def createArr4cellBasedStocks():
    """ Create an array filled with zero, necessary to compute mean and std of cell-based carbon stocks. The shape size is based on the number of columns and rows of a raster file. Obs: use numpy.zeros and do not use numpy.empty (it does not work for the computations used in this script)."""
    refRaster = gdal.Open(parameters.referenceRaster())
    rows = refRaster.RasterYSize
    cols = refRaster.RasterXSize
    zerosArr = np.zeros([rows, cols], dtype=np.float64)
    return zerosArr




def diffEth_cellStocks_analysis(scEth, csType, luRaster, inPath, saveTo):
    "XXXXX"
    luEth = rasterToArray(luRaster)  
    luNoEth = rasterToArray(luScenarios[(luScenarios.index(luRaster))-1])
    luMask = np.ma.masked_where(((luEth <= 2) | (luEth == 10)), luEth)    
    ethFiles = [npFile for npFile in sorted(glob.iglob(os.path.join(inPath, "arrays_sc{}mc*".format(scEth))))]
    noEthFiles = [npFile for npFile in sorted(glob.iglob(os.path.join(inPath, "arrays_sc{}mc*".format(scEth-1))))]
    # Computing MEAN cell-based diff
    acc_diff = np.zeros([885, 854], dtype=np.float64) 
    acc_diffPerLU = np.zeros([12, 1], dtype=np.float64)
    for mcRun_eth, mcRun_noEth in zip(ethFiles, noEthFiles):
        csEth = np.load(mcRun_eth)[csType]
        csNoEth = np.load(mcRun_noEth)[csType]
        diff = csEth-csNoEth
        acc_diff = diff + acc_diff
        # Computation Per LU
        csEthPerLU = getStockPerLUtype(csEth, luEth)        
        csNoEthPerLU = getStockPerLUtype(csNoEth, luNoEth)                
        diffPerLU = (csEthPerLU - csNoEthPerLU).astype(np.float64)
        acc_diffPerLU = diffPerLU + acc_diffPerLU
        #print (mcRun_Eth.split("_")[4], 'mean', np.unique(acc_diff))
        print (scEth, '- mean -', mcRun_eth.split("_")[5], mcRun_noEth.split("_")[5], np.isnan(acc_diff).any())
    mDiff = acc_diff/len(ethFiles)
    mDiffperLU = acc_diffPerLU/len(ethFiles)
    # Computing STD cell-based diff
    acc_sqDiff = np.zeros([885, 854], dtype=np.float64) 
    acc_sqDiffPerLU = np.zeros([12, 1], dtype=np.float64)
    for mcRun_eth, mcRun_noEth in zip(ethFiles, noEthFiles):
        csEth = np.load(mcRun_eth)[csType]
        csNoEth = np.load(mcRun_noEth)[csType]
        diff = csEth-csNoEth        
        maskDiff = (np.ma.masked_values(diff,0)) #mask to remove zero values when doing the computation of squares
        sqDiff = (maskDiff - mDiff)**2
        acc_sqDiff = sqDiff + acc_sqDiff 
        #print (mcRun_eth.split("_")[4], np.unique(acc_sqDiff))
        print (scEth, '- std -', mcRun_eth.split("_")[5], mcRun_noEth.split("_")[5], np.isnan(acc_diff).any())
        # Computation Per LU
        csEthPerLU = getStockPerLUtype(csEth, luEth)        
        csNoEthPerLU = getStockPerLUtype(csNoEth, luNoEth)                
        diffPerLU = (csEthPerLU - csNoEthPerLU).astype(np.float64)
        #maskDiffperLU = (np.ma.masked_values(diffPerLU,0)) #mask to remove zero values when doing the computation of squares
        sqDiffPerLU = (diffPerLU - mDiffperLU)**2
        acc_sqDiffPerLU = sqDiffPerLU + acc_sqDiffPerLU    
    stdDiff = np.sqrt(acc_sqDiff/len(ethFiles)) 
    stdDiff = np.asarray(stdDiff)
    stdDiff[luMask.mask] = 0    
    stdDiffperLU = np.sqrt(acc_sqDiffPerLU/len(ethFiles))
    np.save(os.path.join(saveTo, "dsc{}Xsc{}_mTC_cell".format(scEth,scEth-1)), mDiff.astype(np.float32))
    np.save(os.path.join(saveTo, "dSc{}Xsc{}_stdTC_Cell".format(scEth,scEth-1)), stdDiff.astype(np.float32))
    #print ('SC{} vs SC{}\t|\t ... {} cell-based mean diffs are now saved'.format(csType, scEth, scEth-1))
    if os.path.exists(os.path.join(outDir, "csDiff_{}perLU_mean.txt".format(csType))) == True and os.path.exists(os.path.join(outDir, "csDiff_{}perLU_std.txt".format(csType))) == True:
        with open((os.path.join(outDir, "csDiff_{}perLU_mean.txt".format(csType))),"a+") as mean:
            mean.write('sc{}Xsc{};'.format(sc, sc-1))           
            np.savetxt(mean,mDiffperLU.T, delimiter=';')    
        with open((os.path.join(outDir, "csDiff_{}perLU_std.txt".format(csType))),"a+") as std:
            std.write('sc{}Xsc{};'.format(sc, sc-1))           
            np.savetxt(std,stdDiffperLU.T, delimiter=';')      
    else:
        with open((os.path.join(outDir, "csDiff_{}perLU_mean.txt".format(csType))),"w+") as mean:
            mean.write("scen;"+";".join([str(i) for i in  range(len(mDiffperLU))])+'\nsc{}Xsc{};'.format(sc, sc-1))
            np.savetxt(mean,mDiffperLU.T, delimiter=';')    
        with open((os.path.join(outDir, "csDiff_{}perLU_std.txt".format(csType))),"w+") as std:
            std.write("scen;"+";".join([str(i) for i in  range(len(stdDiffperLU))])+'\nsc{}Xsc{};'.format(sc, sc-1))
            np.savetxt(std,stdDiffperLU.T, delimiter=';')      
    #return mDiff, stdDiff

#-------------------     DATA PROCESSING (cell-based carbon stocks)     -------------------#

if getCellBasedCarbonStock == 1:
    LUC_comp = 1
    for sc, scPath in scenariosDict.iteritems():
        print ('\nSC{} - computing cell-based carbon stocks (tonne C/ha)\t...\tStart time: {}'.format(sc, time.asctime(time.localtime(time.time()))))
        mSOC_cell = cellStocks_analysis('soc', sc, saveTo=resultsDir)
        #mBC_cell = cellStocks_analysis('bc', rLU, sc, inPath=os.path.join(tempDir, 'npz_arrays'), saveTo=tempDir)    
        #mTC_cell = cellStocks_analysis('tc', rLU, sc, inPath=os.path.join(tempDir, 'npz_arrays'), saveTo=tempDir)
        #print ("SC{}\t|\t Cell-based processing took {:.2f} min\n".format(sc, float(time.time()-starttime)/60)))

time_cs_cell = time.time()
print ("Finished. Execution time: {:.2f} min (process) / {:.2f} min (total)".format(float((time_cs_cell - time_cs)/60), float((time_cs_cell - start)/60)))

#-------------------     COMPUTING DIFFERENCE BETWEEN SCENARIOS     -------------------#

# Cell-based diffs
for sc, rLU in enumerate(luScenarios):
    starttime = time.time()
    if sc != 0:
        #print ("SC{}\t|\tGetting cell-based diff of stocks Vs 2012, per MCr...\n".format(sc))
        #diffInit_cellStocks_analysis(sc, "SOC", rLU, inPath=os.path.join(tempDir, 'npz_arrays'), saveTo=outDir)
        #diffInit_cellStocks_analysis(sc, "BC", rLU, inPath=os.path.join(tempDir, 'npz_arrays'), saveTo=outDir)
        #diffInit_cellStocks_analysis(sc, "TC", rLU, inPath=os.path.join(tempDir, 'npz_arrays'), saveTo=outDir)
        if sc < 5 and sc % 2 == 0:
        #if sc % 2 == 0 and sc != 0: # to get scEth
            print ("SC{}\t|\tGetting cell-based diff of stocks (Eth-noEth), per MCr...\n".format(sc))
            #diffEth_cellStocks_analysis(sc, "SOC", rLU, inPath=os.path.join(tempDir, 'npz_arrays'), saveTo=outDir)
            #diffEth_cellStocks_analysis(sc, "BC", rLU, inPath=os.path.join(tempDir, 'npz_arrays'), saveTo=outDir)
            diffEth_cellStocks_analysis(sc, "TC", rLU, inPath=os.path.join(tempDir, 'npz_arrays'), saveTo=outDir)
        print ("SC{}\t|\t processing took {:.2f} min\n".format(sc, float(time.time()-starttime)/60))

# Barchart of delta GHG emissions per LU type (NOT READY because it depends on the cell-based results)

dataLU = ['csDiff_TCperLU_mean.txt', 'csDiff_TCperLU_std.txt']
dataLU_eth_plot = []
dataLU_sc_plot = []

for data in dataLU:
    data_ha = pd.read_csv(os.path.join(outDir, data),  sep=";", header=0, index_col=0)
    data_ha = data_ha.loc[:, (data_ha != 0).any(axis=0)]
    data_GHG_20yrs = (data_ha*2500*44/12.0)/20
    data_GHG_20yrs.columns = ['Natural forest', 'Rangeland', 'Planted forest', 'Crops', 'Grass and shrubs', 'Sugar cane', 'Planted pasture', 'Abandoned']
    data_GHG_20yrs = data_GHG_20yrs.reindex(['Natural forest', 'Planted forest', 'Grass and shrubs', 'Rangeland', 'Crops', 'Planted pasture', 'Sugar cane', 'Abandoned'], axis=1)    
    # To plot considering division by increase in ethanol production... 
    df_dMgt.set_index(data_ha.index, inplace=True)
    data_GHG_20yrs_ETH = (data_GHG_20yrs.div(df_dMgt['dEth_Gj'], axis=0))*-1*1000
    dataLU_eth_plot.append(data_GHG_20yrs_ETH.values.astype(np.float32))
    # To plot considering division by increase in sugar cane area... 
    df_dLU.set_index(data_ha.index, inplace=True)
    data_GHG_20yrs_SC = (data_GHG_20yrs.div(df_dLU['8'], axis=0))*-1
    dataLU_sc_plot.append(data_GHG_20yrs_SC.values.astype(np.float32))
   
fig2, ax = plt.subplots(figsize=(13, 7))
plt.title('GHG emissions per land use type resulting \nfrom an increased ethanol production\n',  color='#6E6E6E', fontsize=11, fontweight='bold', pad=2)
lu_type = ['Natural forest', 'Planted forest', 'Grass and shrubs', 'Rangeland', 'Crops', 'Planted pasture', 'Sugar cane', 'Abandoned']         
lu_colors = ['#385623', '#75b54a','#2efe2e', '#ed7d31', '#cfaa13', '#843c0b', '#ff3787', '#6e6e6e']
scenarios = ['Reference\nscenario', 'High\nproductivity', r'2$^n$$^d$ generation'+'\nSugar cane', r'2$^n$$^d$ generation'+'\nEucalyptus', 'Conservation\npolicies', 'All\nmeasures']
# If plotting total per ethanol production
mTC_lu = dataLU_eth_plot[0].T
stdTC_lu = dataLU_eth_plot[1].T
plt_name = "eth"
plt.ylabel('GHG emissions\n'+r'$gram$ $CO_2$-$eq$/ $Mj$$_E$$_t$$_O$$_H$', fontsize=11)
plt.yticks(np.arange(round(np.min(mTC_lu+stdTC_lu))-2, round(np.max(mTC_lu-stdTC_lu))+1, step=10))
#If plotting total per sugar cane area
mTC_lu = dataLU_sc_plot[0].T
stdTC_lu = dataLU_sc_plot[1].T
plt_name = "scane"
plt.ylabel('GHG emissions\n'+r'$gram$ $CO_2$-$eq$/ $sugar\,cane\,area\,(ha)$', fontsize=11)
plt.yticks(np.arange(round(np.min(mTC_lu+stdTC_lu))-2, round(np.max(mTC_lu-stdTC_lu))+1, step=5))
errorLine = dict(lw=0.5, capsize=1.5, ecolor = '#2E2E2E')
gap = 0.8/len(mTC_lu)
for (i, mean), std in zip((enumerate(mTC_lu)),stdTC_lu):
    lu_classes = np.arange(len(mean))
    plt.bar(0.15+lu_classes+i*gap, mean, width=gap, yerr=std, color=lu_colors[i % len(lu_colors)], edgecolor='#898989', linewidth=0.3, alpha=0.8, error_kw=errorLine)
    plt.legend(lu_type, fontsize=10, frameon=True) #loc='lower right'
plt.grid(which="major", color='0.5', linewidth=0.3)
ax.set_xticklabels('')
ax.set_xticks([0.5,1.5,2.5,3.5,4.5, 5.5],minor=True)
ax.set_xticklabels(scenarios, minor=True, fontsize=11)
if saveFigures == 1:
    plt.savefig(os.path.join(outDir, 'barchart_LU_{}'.format(plt_name)), dpi=700)
plt.show()
