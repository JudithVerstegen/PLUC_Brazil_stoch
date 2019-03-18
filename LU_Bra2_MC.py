"""Land use change model of Brazil
Judith Verstegen, 2019-03-18

"""

from pcraster import *
from pcraster.framework import *
import Parameters
import ParametersProjection
import pickle
import shutil
import math

#######################################

class LandUseType:
  def __init__(self, typeNr, environment, relatedTypeList, suitFactorList, \
               weightList, variableDict, noise, nullMask, \
               windowLengthRealization):
    """Create LandUseType object that represents a class on the land use map.

    Takes ten arguments:
    typeNr -- class nr of the land use type on the land use map
    environment -- global land use map that will evolve
    relatedTypeList -- list with land use type next to which growth is preferred
    suitFactorList -- list of suitability factors the type takes into account
    weightList -- list of relative weights for those factors
    variableDict -- dictionary in which inputs for factors are found
    noise -- very small random noise to ensure cells can't get same suitability
    nullMask -- map with value 0 for study area and No Data outside
    windowLengthRealization -- window length for neighborhood operation in m.
    
    """
    
    self.typeNr = typeNr
    self.environment = environment
    self.relatedTypeList = relatedTypeList
    self.suitFactorList = suitFactorList
    self.weightList = weightList
    self.variableDict = variableDict
    self.nullMask = nullMask

    self.noise = noise
    self.toMeters = Parameters.getConversionUnit()
##    self.stochDistance = Parameters.getStochDistance()
##    self.stochWindow = Parameters.getStochWindow()
    self.windowLengthRealization = windowLengthRealization
    # This new yieldmap approach is a problem for the stochastic mode
    yieldMapName = Parameters.getYieldMapName(typeNr)
    self.yieldMapSuit = scalar(readmap(yieldMapName))
    self.yieldMapSuit = self.yieldMapSuit / mapmaximum(self.yieldMapSuit)
    self.yieldFrac = self.nullMask + 1
    if self.typeNr == Parameters.getForestNr():
      self.forest = True
    else:
      self.forest = False

  def setEnvironment(self, environment):
    """Update the environment (land use map)."""
    self.environment = environment

  def createInitialMask(self, globalMapNoGo, privateMapsNoGo):
    """Combine the global no-go map with areas unsuitable for this land use."""
    self.mask = globalMapNoGo
    if privateMapsNoGo is not None:
      self.mask = pcror(self.mask, privateMapsNoGo)
##        report(self.mask, 'privMask')

  def normalizeMap(self, aMap):
    """Return a normalized version of the input map."""
    mapMax = mapmaximum(aMap)
    mapMin = mapminimum(aMap)
    diff = float(mapMax - mapMin)
    if diff < 0.000001:
      normalizedMap = (aMap - mapMin) / 0.000001
    else:
      normalizedMap = (aMap - mapMin) / diff
    return normalizedMap
  
  # 1
  def getNeighborSuitability(self):
    """Return suitability map based on nr of neighors with a related type."""
    booleanSelf = pcreq(self.environment, self.typeNr)
    for aType in self.relatedTypeList:
      booleanMap = pcreq(self.environment, aType)
      booleanSelf = pcror(booleanSelf, booleanMap)
    scalarSelf = scalar(booleanSelf)
    # Count nr of neighbors with 'true' in a window with length from parameters
    # and assign this value to the centre cell
    variableList = self.variableDict.get(1)
    windowLength = variableList[0]
    nrNeighborsSameLU = windowtotal(scalarSelf, windowLength) - scalarSelf
    # The nr of neighbors are turned into suitability values between 0 and 1
    maxNr = windowtotal(self.nullMask + 1, windowLength) - 1
    # NEW
    # f [0,1]
    f = variableList[1]
    neighborSuitability = -1*(nrNeighborsSameLU**2) + f * 2 * maxNr *\
                               nrNeighborsSameLU
    neighborSuitability = self.normalizeMap(neighborSuitability)
##    report(neighborSuitability, 'neighborSuit' + str(self.typeNr))
    return neighborSuitability

  # 2
  def getDistanceRoadSuitability(self, spreadMapRoads):
    """Return suitability map based on distance to roads."""
    variableList = self.variableDict.get(2)
    a = variableList[0]
    normalized = self.normalizeMap(spreadMapRoads)
    roadSuitability = 1 - (normalized ** a)  
##    report(roadSuitability, 'roadSuit' + str(self.typeNr))
    return roadSuitability

  # 3
  def getDistanceWaterSuitability(self, spreadMapWater):
    """Return suitability map based on distance to water."""
    variableList = self.variableDict.get(3)
    a = variableList[0]
    normalized = self.normalizeMap(spreadMapWater)
    waterSuitability = 1 - (normalized ** a)  
##    report(waterSuitability, 'waterSuit' + str(self.typeNr))
    return waterSuitability

  # 4
  def getDistanceHubsSuitability(self, friction):
    """Return suitability map based on distance to large cities."""
    hubMapNames = Parameters.getHubMapName(self.typeNr)
    hubs = boolean(self.nullMask)
    for aHubMap in hubMapNames:
      hubs = pcror(hubs, cover(readmap(aHubMap),boolean(self.nullMask)))
    distHubs = spread(hubs, 0, friction)
    # The usual way
    variableList = self.variableDict.get(4)
    a = variableList[0]
    normalized = self.normalizeMap(distHubs)
    hubSuitability = 1 - (normalized ** a)  
##    report(hubSuitability, 'hubSuit' + str(self.typeNr))
    return hubSuitability

  # 5
  def getYieldSuitability(self):
    """Return suitability map based on yield for crops or cattle."""
    variableList = self.variableDict.get(5)
    a = variableList[0]
    yieldSuitability = self.yieldMapSuit ** a
##    report(yieldSuitability, 'yieldSuit')
    return yieldSuitability

  # 6
  def getCurrentLandUseSuitability(self):
    """Return suitability map based on current land use type."""
    variableDict = self.variableDict.get(6) 
    current = self.nullMask
    for aKey in variableDict.keys():
      current = ifthenelse(pcreq(self.environment, aKey), \
                           variableDict.get(aKey), current)
    currentLandUseSuitability = self.normalizeMap(current)
##    report(currentLandUseSuitability, 'autoSuit' + str(self.typeNr))
    return currentLandUseSuitability

  # 7
  def getDoubleCropSuitability(self, growthPeriod):
    """Return suitability map based on slope."""
    normalized = self.normalizeMap(growthPeriod)
    variableList = self.variableDict.get(7)
    a = variableList[0]
    growthSuitability = 1 - (normalized ** a)
    return growthSuitability

  # 8
  def getRandomSuitability(self):
    """Return suitability that is completely random."""
    randomSuitability = uniform(boolean(self.yieldFrac))
    return randomSuitability

  # 9
  def getDistanceMillsSuitability(self, mills, friction):
    """Return suitability map based on distance to sugar cane mills."""
    distMills = spread(mills, 0, friction)
    variableList = self.variableDict.get(9)
    a = variableList[0]
    normalized = self.normalizeMap(distMills)
    millSuitability = 1 - (normalized ** a)  
##    report(millSuitability, 'millSuit' + str(self.typeNr))
    return millSuitability
  
  def createInitialSuitabilityMap(self, distRoads, relativeFriction, \
                                  growthPeriod):
    """Return the initial suitability map, i.e. for static factors.

    Given the maps:
    distRoads -- distances to roads
    relativeFriction -- map of road friction: low friction is high speed
    growthPeriod -- length of growing season, proxy double cropping potential
    
    Uses a list and two dictionaries created at construction of the object:
    factors -- the names (nrs) of the suitability factors (methods) needed
    parameters -- the input parameters for those factors
    weights -- the weights that belong to those factors (how they're combined).

    """

    self.weightInitialSuitabilityMap = 0
    self.initialSuitabilityMap = spatial(scalar(0))
    i = 0
    # For every number in the suitability factor list
    # that belongs to a STATIC factor
    # the corresponding function is called providing the necessary parameters
    # and the partial suitability map is added to the total
    # taking into account its relative importance (weight)
    for aFactor in self.suitFactorList:
      if aFactor == 2:
        self.initialSuitabilityMap += self.weightList[i] * \
                                 self.getDistanceRoadSuitability(distRoads)
        self.weightInitialSuitabilityMap += self.weightList[i]
      elif aFactor == 3:
        self.initialSuitabilityMap += self.weightList[i] * \
                                 self.getDistanceWaterSuitability()
        self.weightInitialSuitabilityMap += self.weightList[i]
      elif aFactor == 4:
        self.initialSuitabilityMap += self.weightList[i] * \
                                 self.getDistanceHubsSuitability(relativeFriction)
        self.weightInitialSuitabilityMap += self.weightList[i]
      elif aFactor == 5:
        self.initialSuitabilityMap += self.weightList[i] * \
                                      self.getYieldSuitability()
        self.weightInitialSuitabilityMap += self.weightList[i]
      elif aFactor == 7:
        self.initialSuitabilityMap += self.weightList[i] * \
                                 self.getDoubleCropSuitability(growthPeriod)
        self.weightInitialSuitabilityMap += self.weightList[i]
      elif aFactor in (1, 6, 8, 9):
        ## Dynamic factors are captured in the total suitability map
        pass
      else:
        print('ERROR: unknown suitability factor for landuse', self.typeNr)
      i += 1
    print('weight of initial factors of', self.typeNr, \
          'is', self.weightInitialSuitabilityMap)
    self.initialSuitabilityMap += self.noise
##    report(self.initialSuitabilityMap, 'iniSuit' + str(self.typeNr))

  def getTotalSuitabilityMap(self, mills, relativeFriction):
    """Return the total suitability map for the land use type.

    Uses a lists and two dictionaries:
    factors -- the names (nrs) of the suitability factors (methods) needed
    parameters -- the input parameters for those factors
    weights -- the weights that belong to those factors (how they're combined).

    """

    suitabilityMap = spatial(scalar(0))
    i = 0
    # For every number in the suitability factor list
    # that belongs to a DYNAMIC factor
    # the corresponding function is called providing the necessary parameters
    # and the partial suitability map is added to the total
    # taking into account its relative importance (weight)
    for aFactor in self.suitFactorList:
      if aFactor == 1:
        suitabilityMap += self.weightList[i] * self.getNeighborSuitability()
      elif aFactor == 6:
        suitabilityMap += self.weightList[i] * \
                          self.getCurrentLandUseSuitability()
      elif aFactor == 8:
        suitabilityMap += self.weightList[i] * \
                          self.getRandomSuitability()
      elif aFactor == 9:
        self.initialSuitabilityMap += self.weightList[i] * \
                          self.getDistanceMillsSuitability(mills, \
                                                           relativeFriction)
        self.weightInitialSuitabilityMap += self.weightList[i]
      elif aFactor in (2, 3, 4, 5, 7, 9):
        # Static factors already captured in the initial suitability map
        pass
      else:
        print('ERROR: unknown suitability factor for landuse', self.typeNr)
      i += 1
    suitabilityMap += self.weightInitialSuitabilityMap * \
                      self.initialSuitabilityMap
    self.totalSuitabilityMap = self.normalizeMap(suitabilityMap)
    return self.totalSuitabilityMap

  def setMaxYield(self, maxYield):
    """Set the maximum yield in this time step using the input from the tss."""
    convertedMaxYield = (maxYield / self.toMeters) * cellarea()
    ownMaxYield = ifthen(self.environment == self.typeNr, convertedMaxYield)
    ## maximum yield PER CELL
    self.maxYield = float(mapmaximum(ownMaxYield))
    self.yieldMap = self.yieldFrac * self.maxYield
    
  def updateYield(self, env):
    """Calculate total yield generated by cells occupied by this land use."""
    ## Current cells taken by this land use type
    self.currentYield = ifthen(env == self.typeNr, self.yieldMap)
##    report(self.currentYield, 'currentYield' + str(self.typeNr))
    self.totalYield = float(maptotal(self.currentYield))

  def allocate(self, demand, tempEnvironment, immutables):
    """ Assess total yield, compare with demand and add or remove difference."""
    self.setEnvironment(tempEnvironment)
    self.updateYield(tempEnvironment)
    ownDemand = ifthen(self.environment == self.typeNr, demand)
    self.demand = float(mapmaximum(ownDemand))
    if self.demand < 0:
      self.demand = 0.0
    print('\nland use type', self.typeNr)
    print('demand is:', self.demand)
    if self.forest:
      print('forest,', self.typeNr,'so remove')
      self.removeForest()
    else:
      print('total yield is:', self.totalYield)
      if ((self.totalYield > self.demand) and \
         ((self.totalYield + self.maxYield) < self.demand)) or \
         ((self.totalYield < self.demand) and \
         ((self.totalYield + self.maxYield) > self.demand)):
        print('do nothing')
      elif self.totalYield > self.demand:
        print('remove')
        self.remove()
      elif self.totalYield < self.demand:
        print('add')
        self.add(immutables)
      else:
        print('problem')
    newImmutables = ifthenelse(self.environment == self.typeNr, boolean(1),\
                               immutables)
    return self.environment, newImmutables
    
  def add(self, immutables):
    """Add cells of this land use type until demand is fullfilled."""
    ## Remove cells from immutables (already changed)
    totalSuitabilityMap = ifthen(pcrnot(immutables), \
                                      self.totalSuitabilityMap)
    ## Remove cells already occupied by this land use
    totalSuitabilityMap = ifthen(self.environment != self.typeNr, \
                                      totalSuitabilityMap)
    ## Determine maximum suitability and allocate new cells there
    mapMax = mapmaximum(totalSuitabilityMap)
    print('start mapMax =', float(mapMax))
    ordered = order(totalSuitabilityMap)
    maxIndex = int(mapmaximum(ordered))
    diff = float(self.demand - self.totalYield)
    x = int(maxIndex - diff / self.maxYield)
    xPrev = maxIndex
    i = 0
    tempEnv = self.environment
    while diff > 0 and xPrev > x:
      print('cells to add', int(maxIndex - x))
      if x < 0:
        print('No space left for land use', self.typeNr)
        break
      else:
        ## The key: cells with maximum suitability are turned into THIS type
        tempEnvironment = ifthen(ordered > x, nominal(self.typeNr))
        tempEnv = cover(tempEnvironment, self.environment)

        ## Check the yield of the land use type now that more land is occupied
        self.updateYield(tempEnv)
        i += 1
        xPrev = x
        ## Number of cells to be allocated
        diff = float(self.demand - self.totalYield)
        x -= int(diff / self.maxYield)
    self.setEnvironment(tempEnv)
    print('iterations', i, 'end yield is', self.totalYield)


  def remove(self):
    """Remove cells of this land use type until demand is fullfilled."""
    ## Only cells already occupied by this land use can be removed
    totalSuitabilityMap = ifthen(self.environment == self.typeNr, \
                                      self.totalSuitabilityMap)
    ordered = order(totalSuitabilityMap)
    mapMin = mapminimum(totalSuitabilityMap)
    print('start mapMin =', float(mapMin))
    diff = float(self.totalYield - self.demand)
    # changed maxYield * 0.8 to maxYield, because yield is always 1
    x = int(diff / (self.maxYield * 0.8))
    xPrev = 0
    i = 0
    tempEnv = self.environment
    while diff > 0 and xPrev < x and i < 100:
      print('cells to remove', x)
      ## The key: cells with minimum suitability are turned into 'abandoned'
      tempEnvironment = ifthen(ordered < x, nominal(11))
      tempEnv = cover(tempEnvironment, self.environment)
      
      ## Check the yield of the land use type now that less land is occupied
      self.updateYield(tempEnv)
      i += 1
      xPrev = x
      diff = float(self.totalYield - self.demand)
      if math.fmod(i, 40) == 0:
        print('NOT getting there...')
        ## Number of cells to be allocated
        x = 2 * (x + int(diff / self.maxYield))      
      else:
        ## Number of cells to be allocated
        x += int(diff / self.maxYield)
    self.setEnvironment(tempEnv)
    print('iterations', i, 'end yield is', self.totalYield)
##    report(self.environment, 'newEnv' + str(self.typeNr))

  def removeForest(self):
    """Remove area of forest indicated in time series."""
    if self.demand < 0.01:
      print('nothing to remove')
    else:
      ## Only cells already occupied by this land use can be removed
      self.totalSuitabilityMap = ifthen(self.environment == self.typeNr, \
                                        self.totalSuitabilityMap)
      ordered = order(self.totalSuitabilityMap)
      mapMin = mapminimum(self.totalSuitabilityMap)
      removedBiomass = self.nullMask
      diff = 1
      tempEnv = self.environment
      print('start mapMin =', float(mapMin))
      x = int(self.demand / self.maxYield * 0.8)
      xPrev = 0
      i = 0
      while diff > 0 and xPrev < x and i < 100:
        print('cells to remove', x)
        ## The key: cells with minimum suitability are turned into 'abandoned'
        tempEnvironment = ifthen(ordered < x, nominal(98))
        tempEnv = cover(tempEnvironment, self.environment)
        removed = ifthen(tempEnvironment == 98, nominal(self.typeNr))
        ## Check the yield of the land use type now that less land is occupied
        self.updateYield(removed)
        i += 1
        xPrev = x
        diff = float(self.demand - self.totalYield)
        if math.fmod(i, 40) == 0:
          print('NOT getting there...')
          ## Number of cells to be allocated
          x = 2 * (x + int(diff / self.maxYield))      
        else:
          ## Number of cells to be allocated
          x += int(diff / self.maxYield)
      self.setEnvironment(tempEnv)
      print('iterations', i, 'removed biomass is', self.totalYield)

#######################################

class LandUse:
  def __init__(self, types, environment, nullMask, macroRegions, states):
    """Construct a land use object with a nr of types and an environment."""
    self.types = types
    self.nrOfTypes = len(types)
    print('\nnr of dynamic land use types is:', self.nrOfTypes)
    self.environment = environment
    ## Map with 0 in study area and No Data outside, used for cover() functions
    self.nullMask = nullMask
    self.toMeters = Parameters.getConversionUnit()
    self.yearsDeforestated = nullMask
    self.forest = Parameters.getForestNr()
    self.states = states
    self.macroRegions = macroRegions
    self.macroRegionList = range(1, int(mapmaximum(ordinal(\
      self.macroRegions)))+1)
    self.formerEnvironment = None

  def setEnvironment(self, environment):
    """Update environment of the 'overall' class and separate land use types."""
    self.environment = environment
    for aType in self.landUseTypes:
      aType.setEnvironment(self.environment)
    
  def createLandUseTypeObjects(self, relatedTypeDict, suitabilityDict, \
                               weightDict, variableSuperDict, noise):
    """Generate an object for every dynamic land use type.

    Make objects with:
    typeNr -- class nr in land use map
    environment -- global land use map
    relatedTypes -- list with land use types next to which growth is preferred
    suitFactors -- list with nrs of the needed suitability factors
    weights -- list with relative weights for those factors
    variables -- dictionary with inputs for those factors
    noise -- small random noise that determines order when same suitability

    """
    ## List with the land use type OBJECTS
    self.landUseTypes = []
    self.weightDict = weightDict
    self.variableSuperDict = variableSuperDict
    self.relatedTypeDict = relatedTypeDict
    self.suitabilityDict = suitabilityDict
    self.noise = noise
    
    windowLengthRealization = float(mapnormal())
    
    for aType in self.types:
      ## Get the list that states witch types the current types relates to
      relatedTypeList = relatedTypeDict.get(aType)
      ## Get the right list of suitability factors out of the dictionary
      suitabilityList = suitabilityDict.get(aType)
      ## Get the weights and variables out of the weight dictionary
      weightList = self.weightDict.get(aType)
      variableDict = self.variableSuperDict.get(aType)
      ## Parameter list is notincluded yet
      self.landUseTypes.append(LandUseType(aType, self.environment, \
                                           relatedTypeList, suitabilityList, \
                                           weightList, variableDict, noise, \
                                           self.nullMask, \
                                           windowLengthRealization))
      
  def determineNoGoAreas(self, noGoMap, noGoLanduseList, privateNoGoSlopeDict,\
                         dem):
    """Create global no-go map, pass it to the types that add own no-go areas."""
    self.slopeMap = slope(dem)
    self.excluded = noGoMap
    privateNoGoAreas = None
    ## Check the list with immutable land uses
    if noGoLanduseList is not None:
      for aNumber in noGoLanduseList:
        booleanNoGo = pcreq(self.environment, aNumber)
        self.excluded = pcror(self.excluded, booleanNoGo)
##    report(self.excluded, 'excluded')
    i = 0
    for aType in self.types:
      ## Get land use type specific no-go areas based on slope from dictionary
      ## If not present the variable privateNoGoAreas is 'None'
      aSlope = privateNoGoSlopeDict.get(aType)
      if aSlope is not None:
        privateNoGoAreas = pcrgt(self.slopeMap, aSlope)
      self.landUseTypes[i].createInitialMask(self.excluded, privateNoGoAreas)
      i += 1

  def determineDistanceToRoads(self, booleanMapRoads):
    """Create map with distance to roads, given a boolean map with roads."""
    self.distRoads = spread(booleanMapRoads, 0, 1)
    report(self.distRoads, 'distRoads.map')

  def determineSpeedRoads(self, nominalMapRoads):
    """Create map with relative speed on raods, using boolean map with roads."""
    # By using the part below one can make a map of relative time to
    # reach a hub, giving roads a lower friction
    speed = cover(lookupscalar('speed.txt', nominalMapRoads), \
                  self.nullMask + 5)
    self.relativeFriction = 1.0/speed
    report(self.relativeFriction, 'relativeFriction.map')

  def loadDistanceMaps(self):
    """load the distance maps, when they cannot be kept in memory (fork)"""
##    print os.getcwd()
    self.distRoads = readmap('distRoads')
    self.relativeFriction = readmap('relativeFriction')
  
  def calculateStaticSuitabilityMaps(self, growthPeriod):
    """Get the part of the suitability maps that remains the same."""
    for aType in self.landUseTypes:
      ## Check whether the type has static suitability factors
      ## Those have to be calculated only once (in initial)
      aType.createInitialSuitabilityMap(self.distRoads, \
                                        self.relativeFriction, \
                                        growthPeriod)

  def calculateSuitabilityMaps(self, mills):      
    """Get the total suitability maps (static plus dynamic part)."""
    suitMaps = []
    for aType in self.landUseTypes:
      suitabilityMap = aType.getTotalSuitabilityMap(mills, self.relativeFriction)
      suitMaps.append(suitabilityMap)

  def allocate(self, maxYield, demand):
    """Allocate as much of a land use type as indicated in the demand tss."""
    self.formerEnvironment = self.environment
    for aRegion in self.macroRegionList:
      print('\nREGION', aRegion)
      region = self.macroRegions == aRegion
      immutables = ifthen(region, self.excluded)
      tempEnvironment = ifthen(region, self.environment)
      for aType in self.landUseTypes:
        aType.setMaxYield(maxYield)
        tempEnvironment, immutables = aType.allocate(demand, tempEnvironment, \
                                                     immutables)
      tempEnvironment = cover(tempEnvironment, self.environment)
      self.setEnvironment(tempEnvironment)    

  def growForest(self):
    """Regrow forest at deforestated areas after 10 years."""
    ## Get all cells that are abandoned in the timestep
    deforestated = pcreq(self.environment, 98)
    ## Update map that counts the years a cell is deforestated
    increment = ifthen(deforestated, self.yearsDeforestated + 1)
    self.yearsDeforestated = cover(increment, self.yearsDeforestated)
    ## Regrow forest after 9 years of abandonement, so it's available
    ## again after 10 years
    regrown = ifthen(self.yearsDeforestated == 9, nominal(self.forest))
    reset = ifthen(regrown == nominal(self.forest), scalar(0))
    self.yearsDeforestated = cover(reset, self.yearsDeforestated)
    ## Update environment
    filledEnvironment = cover(regrown, self.environment)
    self.setEnvironment(filledEnvironment)
    
  def getEnvironment(self):
    """Return the current land use map."""
    return self.environment


######################################

class LandUseChangeModel(DynamicModel, MonteCarloModel,\
                         ParticleFilterModel):
  def __init__(self):
    DynamicModel.__init__(self)
    MonteCarloModel.__init__(self)
    ParticleFilterModel.__init__(self)
    setclone('landuse')

  def premcloop(self):
    # mapping of current MC sample to calibration particle
    self.mapping = ParametersProjection.getMappingFromFile('particle_mapping.csv',\
                                                        'New_ID', 'ID_MC5000_2012')
    # Load or construct all input maps
    self.initialEnvironment = self.readmap('landuse')
    self.nullMask = self.readmap('nullMask')
    self.roadsNom = cover(self.readmap('roadsNom'), nominal(self.nullMask))
    self.roads = self.roadsNom != 0
    noGoMap = readmap('noGo')
    yieldFrac = self.readmap('oneMask')
    self.dem = self.readmap('dem')
    self.regions = self.readmap('provinces')
    self.macroRegions = self.readmap('macro')
    self.noGoMap = cover(noGoMap, boolean(self.nullMask))
    self.growthPeriod = cover(self.readmap('growthPeriod'), self.nullMask)
    self.scenario = ParametersProjection.getScenario()
    print('Scenario is', self.scenario)

    # Input values from Parameters file
    # 1. List of landuse types in order of 'who gets to choose first'
    ##self.landUseList = Parameters.getLandUseList()
    self.relatedTypeDict = Parameters.getRelatedTypeDict()
    self.suitFactorDict = Parameters.getSuitFactorDict()
    ##self.weightDict = Parameters.getWeightDict()
    self.variableSuperDictionary = Parameters.getVariableSuperDict()
    self.noGoLanduseList = Parameters.getNoGoLanduseTypes()
    # NEW: SCENARIOS >12 DO NOT ALLOW DEFORESTATION
    if self.scenario > 12:
      forest_nr = Parameters.getForestNr()
      self.noGoLanduseList.append(forest_nr)
    self.privateNoGoSlopeDict = Parameters.getPrivateNoGoSlopeDict()

    # Time step from which the projection starts (before is calibration)
    self.startTime = Parameters.getFirstTimeStep()

    # Uniform map of very small numbers, used to avoid equal suitabilities
    self.noise = uniform(1)/10000


  def initial(self):
    print('-----------------------\n\n')
    # Fixate random seed so that scenario and counterfactual have same
    # demandStoch
    setrandomseed(self.currentSampleNumber()+28)

    # FROM CALIBRATION PERIOD!
    calibration_sample_nr = self.mapping.get(self.currentSampleNumber())
    print(calibration_sample_nr)
    self.landUseList = ParametersProjection.getSampleInput(\
      calibration_sample_nr, 'LUTypes')
    self.weightDict = ParametersProjection.getSampleInput(\
      calibration_sample_nr, 'weights', self.landUseList)
    ##print(self.weightDict)

    
    # Land use map from calibration period (only one, no uncertainty)
    fileName = os.path.join('projIniMaps','ini1')
    self.environment = readmap(fileName)
    ##self.report(self.environment, 'initLanduse')
    self.landUse = LandUse(self.landUseList, self.environment, self.nullMask,\
                           self.macroRegions, self.regions)

    # Create an object for every landuse type in the list
    self.landUse.createLandUseTypeObjects(self.relatedTypeDict, \
                                          self.suitFactorDict, \
                                          self.weightDict, \
                                          self.variableSuperDictionary, \
                                          self.noise)

    # Static suitability factors
    self.landUse.determineNoGoAreas(self.noGoMap, self.noGoLanduseList, \
                                    self.privateNoGoSlopeDict, self.dem)
    self.landUse.determineDistanceToRoads(self.roads)
    self.landUse.determineSpeedRoads(self.roadsNom)
    self.landUse.calculateStaticSuitabilityMaps(self.growthPeriod)
          
    # Draw random numbers between zero and one
    # To determine demand, not used when run deterministically
    self.demandStoch = 0.5


  def dynamic(self):
    timeStep = self.currentTimeStep()
    print('\ntime step', timeStep)

    # start where the calibration has stopped
    # but also save the map of the last time step before calibration
    # to see the transition
    if timeStep == (self.startTime - 1):
      self.report(self.environment, 'landUse')
      os.system('legend --clone landuse.map -f \"legendLU.txt\" %s ' \
                %generateNameST('landUse', self.currentSampleNumber(),timeStep))
    
    if timeStep >= self.startTime:

      # Load the dynamic attributes (mills)
      mills = cover(self.readDeterministic('mills'), boolean(self.nullMask))

      # Get max yield and demand per land use type
      # Yield not taken into account for Brazil, so 1
      maxYield = spatial(scalar(1))

      # demand per macro region!
      # from file with MAGNET results
      fileName = 'demand_projection' + str(self.scenario) + '.txt'
      demand = lookupscalar(fileName, timeStep, self.macroRegions,\
                            self.environment)

      # Suibility maps are calculated
      self.landUse.calculateSuitabilityMaps(mills)

      # Allocate new land use using demands of current time step
      self.landUse.allocate(maxYield, demand)
      self.landUse.growForest()
      self.environment = self.landUse.getEnvironment()

      # Only save 2020 and 2030 to disk
      if timeStep in [15,25]:
        self.report(self.environment, 'landUse')
        os.system('legend --clone landuse.map -f \"legendLU.txt\" %s ' \
                  %generateNameST('landUse', self.currentSampleNumber(),timeStep))

    
  def postmcloop(self):
    print('\nrunning postmcloop...')
    if int(self.nrSamples()) > 1:
      # Stochastic variables for which mean, var and percentiles are needed
      print('...calculating statistics...')
      sampleNumbers = self.sampleNumbers()
      timeSteps = range(1, nrOfTimeSteps + 1)
      variables = []

      # Probability per cell of having a certain land use type
      print('...calculating average land use fractions')
      command = "python LUaverage.py"
      os.system(command)

      # Copy all 2030 maps to the result folder for further processing
      print('...copying land use 2030 maps')
      command = "python copy2030maps.py"
      os.system(command)

      # Also copy the input scripts, to remember the settings of this run
      for aScript in ['LU_Bra2_MC', 'Parameters', 'ParametersProjection']:
        src = str(aScript + '.py')
        dst = os.path.join("results", aScript + '.py')
        shutil.copy2(src, dst)

    print('\n...done')


    
nrOfTimeSteps = Parameters.getNrTimesteps()
nrOfSamples = ParametersProjection.getNrSamplesFromFile('particle_mapping.csv')
myModel = LandUseChangeModel()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
mcModel = MonteCarloFramework(dynamicModel, nrOfSamples)
# Forking does not work on Windows
##mcModel.setForkSamples(True,16)
mcModel.run()
