"""Land use change model, designed for Brazil
Judith Verstegen, 2013-10-04

"""

def getNrTimesteps():
  """Return nr of time steps.

  e.g. 2005 to 2030 is 26 time steps."""

  timesteps = 25
  return timesteps

def getFirstTimeStep():
  """Return nr of the first time step.

  e.g. 8 when the first 7 have been used for calibration."""
  
  firstTimeStep = 8
  return firstTimeStep

def getNrSamples():
  """Return nr of Monte Carlo samples required.

  If Monte Carlo isn't required fill in 1; no statistics will be calculated."""

  samples = 2
  return samples

def getScenario():
  """Return the scenario number for reading the demand file."""
  
  scenario = 3
  return scenario

def getConversionUnit():
  """Return conversion unit for max yield unit to square meters.

  e.g. when max yield in ton/ha fill in 10000."""

  toMeters = 10000
  return toMeters

def getLandUseList():
  """Return list of landuse types in ORDER of 'who gets to choose first'."""
  landUseList = [9, 5, 8, 4, 6]
  return landUseList

def getForestNr():
  """Return class number of land use types considered forest.

  Abandoned land of this type will be called 'deforestation'
  and will grow back to forest after 10 years."""

  forest = 3
  return forest

def getRelatedTypeDict():
  """Return dictionary which type (key) is related to which others (items).

  e.g. relatedTypeDict[3] = [1, 2, 3, 7] means:
  land use type 3 is related to types 1, 2, 3 and 7.
  This is used in suitability factor 1 about neighbors
  of the same or a related type."""
  
  relatedTypeDict = {}
  relatedTypeDict[8] = [8]
  relatedTypeDict[6] = [6]
  relatedTypeDict[5] = [5]
  relatedTypeDict[9] = [9]
  relatedTypeDict[4] = [4]
  
  return relatedTypeDict

def getSuitFactorDict():
  """Return dictionary which type (key) has which suit factors (items).

  e.g. suitFactorDict[1] = [1, 2, 4, 5, 6, 9] means:
  land use type 1 uses suitability factors 1, 2, 4, 5, 6 and 9."""
  
  suitFactorDict = {}
  suitFactorDict[8] = [1, 5, 6, 9]
  suitFactorDict[6] = [1, 4, 5, 6, 7]
  suitFactorDict[5] = [1, 2, 7]
  suitFactorDict[9] = [1, 4, 5]
  suitFactorDict[4] = [1, 2, 5]
  return suitFactorDict

def getWeightDict():
  """Return dictionary how a type (key) weights (items) its suit factors.

  e.g. weightDict[1] = [0.3, 0.1, 0.2, 0.1, 0.2, 0.1] means:
  land use type 1 has suitability factor - weight:
  1 - 0.3
  2 - 0.1
  4 - 0.2
  5 - 0.1
  6 - 0.2
  9 - 0.1

  Note that the number and order of weights has to correspond to the
  suitbility factors in the previous method."""
  
  weightDict = {}
  ## A list with weights in the same order as the suit factors above
  weightDict[8] = [0.29,0.22,0.21,0.28]  
  weightDict[6] = [0.22,0.14,0.11,0.23,0.30]
  weightDict[5] = [0.29,0.34,0.37]
  weightDict[9] = [0.53,0.45,0.02]
  weightDict[4] = [0.46,0.35,0.19]
  return weightDict

def getVariableSuperDict():
  """Return nested dictionary for which type (key1) which factor (item1
  and key2) uses which parameters (items2; a list).

  e.g. variableDict1[2] = [-1, 10000, 1, 2] means:
  land use type 1 uses in suitability factor 2 the parameters:
  -1 for direction of the relation (decreasing)
  10000 for the maximum distance of influence
  1 for friction
  and relation type 'inversely proportional' (=2).

  An explanation of which parameters are required for which suitability
  factor is given in the manual of the model."""

  variableSuperDict = {}
  variableDict8 = {}
  variableDict8[1] = [25000, 0.5]
  variableDict8[5] = [1]
  variableDict8[6] = {6:1, 5:1, 9:1, 11:1}
  variableDict8[9] = [1]
  variableSuperDict[8] = variableDict8
  variableDict6 = {}
  variableDict6[1] = [25000, 0.5]
  variableDict6[4] = [1]
  variableDict6[5] = [1]
  variableDict6[6] = {5:1, 9:1, 11:1}
  variableDict6[7] = [1]
  variableSuperDict[6] = variableDict6
  variableDict5 = {}
  variableDict5[1] = [25000, 0.5]
  variableDict5[2] = [1]
  variableDict5[5] = [1]
  variableDict5[7] = [1]
  variableSuperDict[5] = variableDict5
  variableDict9 = {}
  variableDict9[1] = [25000, 0.5]
  variableDict9[4] = [1]
  variableDict9[5] = [1]
  variableSuperDict[9] = variableDict9
  variableDict4 = {}
  variableDict4[1] = [25000, 0.5]
  variableDict4[2] = [1]
  variableDict4[5] = [1]
  variableSuperDict[4] = variableDict4
  return variableSuperDict

def getNoGoLanduseTypes():
  """Return a list of land use type numbers that cannot be changed

  At the moment this only works for static land uses
  The noGo map is calculated once at the beginning of the run."""

  noGoLanduse = [1, 2]
  return noGoLanduse

def getPrivateNoGoSlopeDict():
  """Return dictionary of a type's (key) slope contrains (items; mapname).

  e.g. privateNoGoDict[1] = 0.16 means:
  land use type 1 can not allocate on locations that are have a slope
  of more than 16%."""
  
  privateNoGoDict = {}
  privateNoGoDict[8] = 0.12
  privateNoGoDict[6] = 0.16
  return privateNoGoDict

def getHubMapName(typeNr):
  """Return the names of the hubs map for this land use type (mandatory)."""
  hubMapNameDict = {}
  hubMapNameDict[6] = ['crop_facil', 'fert_facil']
  hubMapNameDict[9] = ['livestock_facil', 'slaugther_houses']
  needed = hubMapNameDict.get(typeNr)
  return needed


def getYieldMapName(typeNr):
  """Return the name of the yield map for this land use type (mandatory)."""
  yieldMapNameDict = {}
  yieldMapNameDict[8] = 'scYield'
  yieldMapNameDict[6] = 'cropYield'
  yieldMapNameDict[5] = 'growthPeriod'
  yieldMapNameDict[9] = 'pastureSuit'
  yieldMapNameDict[4] = 'pastureSuit'

  needed = yieldMapNameDict.get(typeNr)
  return needed
