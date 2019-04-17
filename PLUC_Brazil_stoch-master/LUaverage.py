"""Test transition matrix to create for Magnet
Judith Verstegen, 2014-01-22

"""

from pcraster import *
from pcraster.framework import *
import csv
import Parameters
import ParametersProjection

class CheckModel(DynamicModel):
   def __init__(self):
      DynamicModel.__init__(self)
      setclone('nullMask')
      
   def initial(self):
      self.nullMask = self.readmap('nullMask')
      self.startTime = Parameters.getFirstTimeStep()
      self.mapping = ParametersProjection.getMappingFromFile('particle_mapping.csv',\
                                                      'New_ID', 'Nr_particles_weight')
      
   def dynamic(self):
      timestep = self.currentTimeStep()
      if timestep in [15,25]:
         for aType in LUtypes:
            subDict = {}
            summap = self.nullMask
            for aSample in samples:
               path = os.path.join(str(aSample),'landUse')
               landuse = self.readmap(path)
               thisLU = scalar(landuse == aType)
               weight = self.mapping.get(aSample)
               summap = summap + (thisLU * float(weight))
            total = sum(list(self.mapping.values())[0:nr_samples])
            ##print(total)
            averageMap = summap / float(total)
            path2 = os.path.join('results', str(aType) + 'L')
            self.report(averageMap, path2)
      
nr_samples = ParametersProjection.getNrSamplesFromFile('particle_mapping.csv')      
samples = range(1,nr_samples + 1)                                 
LUtypes = range(1,12)
nrOfTimeSteps = Parameters.getNrTimesteps()
myModel = CheckModel()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
dynamicModel.run()
