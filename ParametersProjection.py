"""Land use change model, designed for Brazil
Judith Verstegen, 2013-10-04

"""
import numpy as np
import os
from pcraster import *
setclone('landuse')

def getScenario():
  scenario = 3
  return scenario

def getFirstTimeStep():
  firstTimeStep = 8
  return firstTimeStep

def getLandUseDict():
  """Return list of landuses  with name."""
  landUseDict = {8: 'sc', 6: 'cr', 5: 'pf', 9: 'pp', 4: 'ra'}
  return landUseDict

def getSampleInput(sampleNr, filename, landUseList = None):

  if landUseList is not None:
    inputs = {}
    for aName in landUseList:
      fileName = os.path.join('inputsFromCalibration', \
                              filename + str(aName) + '.npy')
      array = np.load(fileName)
##      print(array.shape)
##      print array
      # take from array:
      # timestep (final of calibration, sampleNr, x (y is nothing)
      parameter = array[-1, sampleNr - 1, :]
      parameter = np.squeeze(parameter)
##      print parameter
      inputs[aName] = parameter.tolist()
  else:
    fileName = os.path.join('inputsFromCalibration', \
                            filename + '.npy')
    array = np.load(fileName)
##    print(array.shape)
##    print(array)
    # take from array:
    # timestep (final of calibration, sampleNr, x (y is nothing)
    parameter = array[-1, sampleNr - 1, :]
    parameter = np.squeeze(parameter)
##    print parameter
    inputs = parameter.tolist()
      
  return inputs

def getInitialLandUseMap(sampleNr):
    fileName = os.path.join('projIniMaps','ini' + str(sampleNr))
    theMap = readmap(fileName)
    return theMap
    
# test
##getInitialLandUseMap(1)
