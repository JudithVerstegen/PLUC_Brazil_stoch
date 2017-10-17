"""Test transition matrix to create for Magnet
Judith Verstegen, 2014-01-22

"""

from pcraster import *
from pcraster.framework import *
import csv
import Parameters

class CheckModel(DynamicModel):
   def __init__(self):
      DynamicModel.__init__(self)
      setclone('nullMask')
      
   def initial(self):
      self.nullMask = self.readmap('nullMask')
      self.startTime = Parameters.getFirstTimeStep()
      
   def dynamic(self):
      timestep = self.currentTimeStep()
      if timestep >= self.startTime:
         for aType in LUtypes:
            subDict = {}
            summap = self.nullMask
            for aSample in samples:
               path = os.path.join(str(aSample),'landUse')
               landuse = self.readmap(path)
               thisLU = scalar(landuse == aType)
               summap = summap + thisLU
            averageMap = summap / float(len(samples))
            path2 = os.path.join('results', str(aType) + 'L')
            self.report(averageMap, path2)
      
      
samples = range(1,Parameters.getNrSamples() + 1)                                 
LUtypes = range(1,12)
nrOfTimeSteps = Parameters.getNrTimesteps()
myModel = CheckModel()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
dynamicModel.run()
