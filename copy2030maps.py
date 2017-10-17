import os
import numpy as np
import Parameters
import sys
import shutil
from pcraster import *

nrOfSamples=Parameters.getNrSamples()
nrOfTimesteps=Parameters.getNrTimesteps()

sampleNumbers=range(1,nrOfSamples+1,1)
timeSteps=range(1,nrOfTimesteps+1,1)
print timeSteps

for aSample in sampleNumbers:
  src = os.path.join(str(aSample), 'landUse0.025')
  dst = os.path.join('results', '2030_' + str(aSample) + '.map')
  shutil.copy2(src, dst)
  # also for 2015, 2020, and 2025 for the optimization
  src = os.path.join(str(aSample), 'landUse0.020')
  dst = os.path.join('results', '2025_' + str(aSample) + '.map')
  shutil.copy2(src, dst)
  src = os.path.join(str(aSample), 'landUse0.015')
  dst = os.path.join('results', '2020_' + str(aSample) + '.map')
  shutil.copy2(src, dst)
  src = os.path.join(str(aSample), 'landUse0.010')
  dst = os.path.join('results', '2015_' + str(aSample) + '.map')
  shutil.copy2(src, dst)