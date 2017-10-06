'''
This file describes the basic moments, calculated using
the delta-function method identified in norris-etal-JPCM-2009.
'''

import sys
import numpy as np
#from pycraters.statistics.GENERIC import *


class base_moments(object):


  # initialization
  def __init__(self):
    pass


  def collect(self, idata, edata, rdata, params):

    lstats = dict()
    species = len(params.target)

    # storage
    lstats["iidM0"]   = 0
    lstats["iidM1"]   = np.zeros((3))
    lstats["iidM2"]   = np.zeros((3, 3))

    lstats["evdM0"]   = np.zeros((species+1))
    lstats["evdM1"]   = np.zeros((species+1, 3))
    lstats["evdM2"]   = np.zeros((species+1, 3, 3))

    lstats["m0e"]     = np.zeros((species+1))
    lstats["m1e"]     = np.zeros((species+1, 3))
    lstats["m2e"]     = np.zeros((species+1, 3, 3))

    lstats["rvdM0"]   = np.zeros((species+1, 3))
    lstats["rvdM1"]   = np.zeros((species+1, 3))
    lstats["rvdM2"]   = np.zeros((species+1, 3, 3))

    lstats["ridM0"]   = np.zeros((species+1, 3))
    lstats["ridM1"]   = np.zeros((species+1, 3))
    lstats["ridM2"]   = np.zeros((species+1, 3, 3))

    lstats["rddM0"]   = np.zeros((species+1, 3))
    lstats["rddM1"]   = np.zeros((species+1, 3))
    lstats["rddM2"]   = np.zeros((species+1, 3, 3))

    lstats["m0r"]     = np.zeros((species+1))
    lstats["m1r"]     = np.zeros((species+1, 3))
    lstats["m2r"]     = np.zeros((species+1, 3, 3))


    # how do implanted ions affect the statistics?
    if idata != None and idata != []:
      sid = 0
      finals = np.array( [entry[1] for entry in idata if entry[0] == sid] )
      if len(finals) > 0:
        lstats["iidM0"] = len(finals)
        lstats["iidM1"] = np.sum(finals, axis=0)
        lstats["iidM2"] = np.einsum("ab,ac->bc", finals, finals)


    if edata != None and edata != []:
      for sid in xrange(1, species+1):
        starts = np.array( [entry[1] for entry in edata if entry[0] == sid] )
        if len(starts) == 0:
          continue

        lstats["evdM0"][sid] = len(starts)
        lstats["evdM1"][sid] = np.sum(starts, axis=0)
        lstats["evdM2"][sid] = np.einsum("ab,ac->bc", starts, starts)


    if rdata != None and rdata != []:
      for sid in xrange(1, species+1):
        starts = np.array( [entry[1] for entry in rdata if entry[0] == sid] )
        finals = np.array( [entry[2] for entry in rdata if entry[0] == sid] )
        displs = finals - starts
        if len(displs) == 0:
          continue

        lstats["rvdM0"][sid] = len(starts) 
        lstats["rvdM1"][sid] = np.sum(starts, axis=0)
        lstats["rvdM2"][sid] = np.einsum("ab,ac->bc", starts, starts)

        lstats["ridM0"][sid] = len(finals)
        lstats["ridM1"][sid] = np.sum(finals, axis=0)
        lstats["ridM2"][sid] = np.einsum("ab,ac->bc", finals, finals)

        lstats["rddM0"][sid] = len(displs)
        lstats["rddM1"][sid] = np.sum(displs, axis=0)
        lstats["rddM2"][sid] = np.einsum("ab,ac->bc", displs, displs)
     

    # "old-style" moments
    lstats["m0e"] = - lstats["evdM0"]
    lstats["m1e"] = - lstats["evdM1"]
    lstats["m2e"] = - lstats["evdM2"]
    lstats["m0r"] = - lstats["rvdM0"] + lstats["ridM0"]
    lstats["m1r"] = - lstats["rvdM1"] + lstats["ridM1"]
    lstats["m2r"] = - lstats["rvdM2"] + lstats["ridM2"]



    # now update moment values (total of all species in position '0')
    for part1 in ["evd", "rvd", "rid", "rdd"]:
      for part2 in ["M0", "M1", "M2"]:
        key = "%s%s" % (part1, part2)
        lstats[key][0] = np.sum(lstats[key][1:], axis=0)

    for part1 in ["0", "1", "2"]:
      for part2 in ["e", "r"]:
        key = "m%s%s" % (part1, part2)
        lstats[key][0] = np.sum(lstats[key][1:], axis=0)

    return lstats


