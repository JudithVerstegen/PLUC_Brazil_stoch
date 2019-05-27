"""Land use change model, designed for Brazil
Judith Verstegen, 2013-10-04

"""
import numpy as np
import os
from pcraster import *
setclone('landuse')

def getScenario():
  scenario = 15
  return scenario

def getFirstTimeStep():
  firstTimeStep = 8
  return firstTimeStep

def getLandUseDict():
  """Return list of landuses  with name."""
  landUseDict = {8: 'sc', 6: 'cr', 5: 'pf', 9: 'pp', 4: 'ra'}
  return landUseDict

def getNrSamplesFromFile(filename):
  """Return number of samples to run."""
  # read from csv file that is result of particle filter
  path = os.path.join('inputsFromCalibration',filename)
  data = np.genfromtxt(path, dtype=int, delimiter='\t', names=True)
  # count number of unique particles
  unique = data.shape[0]
  return unique

def getMappingFromFile(filename, col1, col2):
  """Return mapping of new vs original particle number."""
  # read from csv file that is result of particle filter
  path = os.path.join('inputsFromCalibration',filename)
  data = np.genfromtxt(path, dtype=int, delimiter='\t', names=True)
  # make dictionary of new id : old id
  mapping = dict(zip(data[col1], data[col2]))
  ##print(mapping)
  return mapping

def getSampleInput(sampleNr, filename, landUseList = None):

  if landUseList is not None:
    inputs = {}
    for aName in landUseList:
      fileName = os.path.join('inputsFromCalibration', \
                              filename + str(aName) + '.npy')
      array = np.load(fileName)
##      print(array.shape)
##      print(array)
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
##mapping = getMappingFromFile('particle_mapping.csv','New_ID', 'ID_MC5000_2012')
##calibration_sample_nr = mapping.get(43)
##print(calibration_sample_nr)

##j = 4650
##print(j)
##print(getNrSamplesFromFile('particle_mapping.csv'))
##landUseList = getSampleInput(j, 'LUTypes')
##print(landUseList)
##weights = getSampleInput(j, 'weights', landUseList)
##for i in weights.keys():
##  print(i, weights[i])

