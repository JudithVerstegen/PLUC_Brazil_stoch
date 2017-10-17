"""Land use map
Making a norminal map from maps with fractions per land use type
The total number of cells allocated should be the same as the
total of the fractions (demand)
Judith Verstegen, 2013-09-03

"""
import random
from pcraster import *
#from PCRaster.Framework import *
from pcraster.numpy_operations import *
import numpy as np
import os

def getInitialLandUseMap(sd, uncertain = True):
  if uncertain == False:
    landusemap = readmap('landuse')
  else:
  ##setclone('nullMask')
    nullMask = readmap('nullMask')
    provinces = readmap('provinces')

    # Landuses and the number they should have in the nominal map
    # Numbers are chosen for the color they will get in Aguila
    landuses = {1:'urban', 2:'water', 3:'naturalForest', 4:'naturalPast', \
                5:'plantedForest', 6: 'crop', 7:'grass', 8:'sc', 9:'plantedPast', \
                10:'bare'}

    # To be able to allocate in a specific order we need a list, as
    # a dictionary has no order.
    order = [1,2,10,8,6,5,9,4,3,7]
    provinceList = [11,12,13,14,15,16,17,21,22,23,24,25,26,27,28,29,31,\
                 32,33,35,41,42,43,50,51,52,53]


    total = nullMask
    landusemap = nominal(nullMask + 11)
    for aKey in order:
      fracmap = readmap(os.path.join('initialMap', landuses[aKey] \
                                     + 'Frac5'))
      if aKey in [1,2,10]:
        noise = uniform(boolean(nullMask + 1)) / 10000
        fracmap += noise
        # TEST
##        fracmap = fracmap + normal(1) * sd * fracmap
##        fracmap = ifthenelse(pcrle(fracmap, 0), scalar(0), fracmap)
##        fracmap = ifthenelse(pcrge(fracmap, 1), scalar(1), fracmap)
      else:
        fracmap = fracmap + normal(1) * sd * fracmap
        fracmap = ifthenelse(pcrle(fracmap, 0), scalar(0), fracmap)
        fracmap = ifthenelse(pcrge(fracmap, 1), scalar(1), fracmap)
      # The total that should be allocated
      # Now the total of the fractions, but can also be a manual input
      maptotalfrac = float(maptotal(fracmap))
##      print landuses[aKey], (maptotalfrac * 25)
    ##  # Cells still available
      fracmap1 = ifthen(total < 1, fracmap)
      mask = ifthen(total < 1, nullMask)
      zvalues = mask
      
      for aProv in provinceList:
        # Calculate the demand as a fraction of the total region area
        thisProvince = provinces == aProv
        thisMask = ifthen(thisProvince, mask)
        provFrac = ifthen(thisProvince, fracmap)
    ##    provFrac = provFrac + noise
        provtotalfrac = float(maptotal(provFrac))
        if float(maptotal(thisMask + 1)) > 0:
          fraction = (provtotalfrac / float(maptotal(thisMask + 1))) * 100
        else:
          fraction = 100
    ##    print 'fraction', fraction
        # Calculate the percentage above which a cell should become cropland
        # pcr2numpy(map,mv)
        fracArray = pcr2numpy(ifthen(thisMask == 0, provFrac), -1)
        fracArray = fracArray[fracArray != -1]
        if fraction == 0.0:
          zvalue == 1
    ##      print 'here'
        elif fraction < 100:
          zvalue = float(np.percentile(fracArray, 100 - fraction))
        else:
          zvalue = 0.001
    ##      print 'here1'
    ##    print 'z', zvalue
        zvalues = ifthenelse(thisProvince, zvalues + zvalue, zvalues)
      # New method with the z-values
      booleanmap = ifthenelse(fracmap1 >= zvalues , boolean(1), boolean(0))
    ##  report(booleanmap, landuses[aKey] + '.map')
##      print float(maptotal(scalar(booleanmap)) * 25)
      total += cover(scalar(booleanmap), nullMask)
      landusemap = cover(ifthen(booleanmap, nominal(aKey)), landusemap)
                         

##    report(total, 'total')
    count = ifthen(total == 0, scalar(1))
##    print 'number of cells with no landuse is', float(maptotal(count))
    landusemap = ifthen(pcrnot(pcreq(landusemap, 11)), landusemap)
    majority = windowmajority(landusemap, 15000)
    landusemap = ifthenelse(total == 0, majority, landusemap)
    # Something with nominal value 11
  return landusemap

