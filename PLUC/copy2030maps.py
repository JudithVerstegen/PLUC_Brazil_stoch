import os
import numpy as np
import Parameters
import ParametersProjection
import sys
import shutil
from pcraster import *

nrOfSamples= ParametersProjection.getNrSamplesFromFile('particle_mapping.csv')    
sampleNumbers=range(1,nrOfSamples+1,1)

for aSample in sampleNumbers:
  src = os.path.join(str(aSample), 'landUse0.025')
  dst = os.path.join('results', 'sc' + str(ParametersProjection.getScenario()),'2030_' + str(aSample) + '.map') # Renan: I added the "sc" folder here
  shutil.copy2(src, dst)
  # also 2020, can be removed if required
  src = os.path.join(str(aSample), 'landUse0.015')
  dst = os.path.join('results', 'sc' + str(ParametersProjection.getScenario()), '2020_' + str(aSample) + '.map') # Renan: I added the "sc" folder here
  shutil.copy2(src, dst)
