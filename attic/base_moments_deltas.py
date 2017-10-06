'''
This file describes the basic moments, calculated using
the delta-function method identified in norris-etal-JPCM-2009.
'''

import sys
import numpy as np
from pycraters.statistics.GENERIC_Statistics import *


class base_moments_deltas(GENERIC_Statistics):


  # initialization
  def __init__(self, species=1, bdata=None):

    self.species = species
    self.bdata = bdata
    pass



  def collect(self, idata, edata, rdata):

    #print "  collecting statistics ...",
    #sys.stdout.flush()

    lstats = dict()
    species = self.species
    bdata   = self.bdata

    # storage
    lstats["iidM0"]   = 0
    lstats["iidM1"]   = np.zeros((3))
    lstats["iidM2"]   = np.zeros((3, 3))

    lstats["evdM0"]   = np.zeros((species+1))
    lstats["evdM1"]   = np.zeros((species+1, 3))
    lstats["evdM2"]   = np.zeros((species+1, 3, 3))

    lstats["rvdM0"]   = np.zeros((species+1, 3))
    lstats["rvdM1"]   = np.zeros((species+1, 3))
    lstats["rvdM2"]   = np.zeros((species+1, 3, 3))

    lstats["ridM0"]   = np.zeros((species+1, 3))
    lstats["ridM1"]   = np.zeros((species+1, 3))
    lstats["ridM2"]   = np.zeros((species+1, 3, 3))

    lstats["rddM0"]   = np.zeros((species+1, 3))
    lstats["rddM1"]   = np.zeros((species+1, 3))
    lstats["rddM2"]   = np.zeros((species+1, 3, 3))

    # ---

    if bdata != None:

      xbins = bdata.bincount[0]
      ybins = bdata.bincount[1]
      zbins = bdata.bincount[2]

      lstats["actXY"]   = np.zeros((xbins, ybins))
      lstats["actXZ"]   = np.zeros((xbins, zbins))
      lstats["actYZ"]   = np.zeros((ybins, zbins))

      lstats["iidXY"]   = np.zeros((xbins, ybins))
      lstats["iidXZ"]   = np.zeros((xbins, zbins))
      lstats["iidYZ"]   = np.zeros((ybins, zbins))

      lstats["evdXY"]   = np.zeros((xbins, ybins))
      lstats["evdXZ"]   = np.zeros((xbins, zbins))
      lstats["evdYZ"]   = np.zeros((ybins, zbins))

      lstats["rvdXY"]   = np.zeros((xbins, ybins))
      lstats["rvdXZ"]   = np.zeros((xbins, zbins))
      lstats["rvdYZ"]   = np.zeros((ybins, zbins))

      lstats["ridXY"]   = np.zeros((xbins, ybins))
      lstats["ridXZ"]   = np.zeros((xbins, zbins))
      lstats["ridYZ"]   = np.zeros((ybins, zbins))

      lstats["rddXYx"]   = np.zeros((xbins, ybins))
      lstats["rddXYy"]   = np.zeros((xbins, ybins))
      lstats["rddXZx"]   = np.zeros((xbins, zbins))
      lstats["rddXZz"]   = np.zeros((xbins, zbins))
      lstats["rddYZy"]   = np.zeros((ybins, zbins))
      lstats["rddYZz"]   = np.zeros((ybins, zbins))


    # how do implanted ions affect the statistics?
    if idata != None:
      for entry in idata:

        # extract stored data
        sid  = entry[0]	# species
        pos0 = entry[1]	# final position

        # update moment values
        lstats["iidM0"] += 1
        lstats["iidM1"][:] += pos0
        lstats["iidM2"][:, :] += np.outer(pos0, pos0)

        if bdata==None: continue

        # 2D projections of distributions
        binx = floor((pos0[0]-geom.xMin) / bdata.binsize[0]) ;
        biny = floor((pos0[1]-geom.yMin) / bdata.binsize[1]) ;
        binz = floor((pos0[2]-geom.zMin) / bdata.binsize[2]) ;

        if (binx >=0 and binx < xbins and biny >=0 and biny < ybins and binz >=0 and binz < zbins):

          lstats["iidXY"][sid, binx, biny] += 1
          lstats["iidXZ"][sid, binx, binz] += 1
          lstats["iidYZ"][sid, biny, binz] += 1


    # how do eroded atoms affect the statistics?
    if edata != None:
      for entry in edata:

        # extract stored data
        sid  = entry[0]	# species ID
        pos0 = entry[1]	# initial position of sputtered atom

        # update moment values (per species)
        lstats["evdM0"][sid] += 1
        lstats["evdM1"][sid,:] += pos0
        lstats["evdM2"][sid,:,:] += np.outer(pos0, pos0)
  
        if bdata==None: continue

        # 2D projections of distributions
        binx = floor((pos0[0]-geom.xMin) / bdata.binsize[0]) ;
        biny = floor((pos0[1]-geom.yMin) / bdata.binsize[1]) ;
        binz = floor((pos0[2]-geom.zMin) / bdata.binsize[2]) ;

        if (binx >=0 and binx < xbins and biny >=0 and biny < ybins and binz >=0 and binz < zbins):

          lstats["actXY"][sid, binx, biny] += 1
          lstats["actXZ"][sid, binx, binz] += 1
          lstats["actYZ"][sid, biny, binz] += 1

          lstats["evdXY"][sid, binx, biny] += 1
          lstats["evdXZ"][sid, binx, binz] += 1
          lstats["evdYZ"][sid, biny, binz] += 1


    # how do redistributed atoms affect the statistics?
    if rdata != None:
      for entry in rdata:

        # extract stored data
        sid  = entry[0]	# species ID
        pos0 = entry[1]	# initial position
        pos1 = entry[2]	# final position
        dpos = pos1-pos0	# displacement

        # update moment values (per species)
        lstats["rvdM0"][sid] += 1
        lstats["rvdM1"][sid, :] += pos0
        lstats["rvdM2"][sid, :, :] += np.outer(pos0, pos0)

        lstats["ridM0"][sid] += 1
        lstats["ridM1"][sid, :] += pos1
        lstats["ridM2"][sid, :, :] += np.outer(pos1, pos1)

        lstats["rddM0"][sid] += 1
        lstats["rddM1"][sid, :] += dpos
        lstats["rddM2"][sid, :, :] += np.outer(dpos, dpos)

        if bdata==None: continue

        # 2D projections of distributions
        binx = floor((pos0[0]-geom.xMin) / bdata.binsize[0]) ;
        biny = floor((pos0[1]-geom.yMin) / bdata.binsize[1]) ;
        binz = floor((pos0[2]-geom.zMin) / bdata.binsize[2]) ;

        if (binx >=0 and binx < xbins and biny >=0 and biny < ybins and binz >=0 and binz < zbins):

          lstats["actXY"][sid, binx, biny] += 1
          lstats["actXZ"][sid, binx, binz] += 1
          lstats["actYZ"][sid, biny, binz] += 1

          lstats["rvdXY"][sid, binx, biny] += 1
          lstats["rvdXZ"][sid, binx, binz] += 1
          lstats["rvdYZ"][sid, biny, binz] += 1

          lstats["ridXY"][sid, binx, biny] += 1
          lstats["ridXZ"][sid, binx, binz] += 1
          lstats["ridYZ"][sid, biny, binz] += 1

          lstats["rddXYx"][sid, binx, biny] += dpos[0]
          lstats["rddXYy"][sid, binx, biny] += dpos[1]
          lstats["rddXZx"][sid, binx, binz] += dpos[0]
          lstats["rddXZz"][sid, binx, binz] += dpos[2]
          lstats["rddYZy"][sid, biny, binz] += dpos[1]
          lstats["rddYZz"][sid, biny, binz] += dpos[2]


    # now update moment values (total of all species in position '0')
    for part1 in ["evd", "rvd", "rid", "rdd"]:
      for part2 in ["M0", "M1", "M2"]:
        key = "%s%s" % (part1, part2)
        lstats[key][0] = np.sum(lstats[key][1:], axis=0)

    if bdata != None:

      # and update 2D projected distributions (also total of all species in position '0')
      for part1 in ["act", "evd", "rvd", "rid"]:
        for part2 in ["XY", "XZ", "YZ"]:
          key = "%s%s" % (part1, part2)
          lstats[key][0] = np.sum(lstats[key][1:], axis=0)

      # and update 2D projected distributions (also total of all species in position '0')
      for part1 in ["rdd"]:
        for part2 in ["XYx", "XYy", "XZx", "XZz", "YZy", "YZz"]:
          key = "%s%s" % (part1, part2)
          lstats[key][0] = np.sum(lstats[key][1:], axis=0)



    # now provide some x, y, z data reflecting the axes we just created





    #print "done."
    #sys.stdout.flush()

    return lstats





