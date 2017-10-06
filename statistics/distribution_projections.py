'''
This file describes the basic moments, calculated using
the delta-function method identified in norris-etal-JPCM-2009.
'''

import sys
import numpy as np
#from pycraters.statistics.GENERIC import *



class distribution_projections(object):


  # initialization
  def __init__(self, boxdims, bincounts):

    self.xdim  =   float(boxdims[0])
    self.ydim  =   float(boxdims[1])
    self.zdim  =   float(boxdims[2])

    self.xbins =   bincounts[0]
    self.ybins =   bincounts[1]
    self.zbins =   bincounts[2]

    self.dx    =   self.xdim / self.xbins
    self.dy    =   self.ydim / self.ybins
    self.dz    =   self.zdim / self.zbins

    self.xMin  = - boxdims[0] / 2.0
    self.yMin  = - boxdims[1] / 2.0
    self.zMin  = - boxdims[2] / 2.0

    self.xMax  =   boxdims[0] / 2.0
    self.yMax  =   boxdims[1] / 2.0
    self.zMax  =   boxdims[2] / 2.0

    xi = np.arange( self.xMin + self.dx/2.0, self.xMax, self.dx )
    yi = np.arange( self.yMin + self.dy/2.0, self.yMax, self.dy )
    zi = np.arange( self.zMin + self.dz/2.0, self.zMax, self.dz )

    self.XYx, self.XYy = np.meshgrid(xi, yi, indexing='ij')
    self.YZy, self.YZz = np.meshgrid(yi, zi, indexing='ij')
    self.XZx, self.XZz = np.meshgrid(xi, zi, indexing='ij')

    




  def collect(self, idata, edata, rdata, params):

    #print "  collecting statistics ...",
    #sys.stdout.flush()

    lstats = dict()
    species = len(params.target)

    # co-ordinate meshes
    lstats["XYx"] = self.XYx
    lstats["XYy"] = self.XYy
    lstats["YZy"] = self.YZy
    lstats["YZz"] = self.YZz
    lstats["XZx"] = self.XZx
    lstats["XZz"] = self.XZz

    # allocate storage for the other arrays
    lstats["iidXY"]   = np.zeros((self.xbins, self.ybins))
    lstats["iidXZ"]   = np.zeros((self.xbins, self.zbins))
    lstats["iidYZ"]   = np.zeros((self.ybins, self.zbins))

    lstats["evdXY"]   = np.zeros((species+1, self.xbins, self.ybins))
    lstats["evdXZ"]   = np.zeros((species+1, self.xbins, self.zbins))
    lstats["evdYZ"]   = np.zeros((species+1, self.ybins, self.zbins))

    lstats["rvdXY"]   = np.zeros((species+1, self.xbins, self.ybins))
    lstats["rvdXZ"]   = np.zeros((species+1, self.xbins, self.zbins))
    lstats["rvdYZ"]   = np.zeros((species+1, self.ybins, self.zbins))

    lstats["ridXY"]   = np.zeros((species+1, self.xbins, self.ybins))
    lstats["ridXZ"]   = np.zeros((species+1, self.xbins, self.zbins))
    lstats["ridYZ"]   = np.zeros((species+1, self.ybins, self.zbins))

    lstats["rddXYx"]   = np.zeros((species+1, self.xbins, self.ybins))
    lstats["rddXYy"]   = np.zeros((species+1, self.xbins, self.ybins))
    lstats["rddXZx"]   = np.zeros((species+1, self.xbins, self.zbins))
    lstats["rddXZz"]   = np.zeros((species+1, self.xbins, self.zbins))
    lstats["rddYZy"]   = np.zeros((species+1, self.ybins, self.zbins))
    lstats["rddYZz"]   = np.zeros((species+1, self.ybins, self.zbins))


    # how do implanted ions affect the statistics?
    if idata != None:
      for entry in idata:

        # extract stored data
        sid  = entry[0]	# species
        pos0 = entry[1]	# final position

        # 2D projections of distributions
        binx = np.floor( (pos0[0]-self.xMin) / self.dx)
        biny = np.floor( (pos0[1]-self.yMin) / self.dy)
        binz = np.floor( (pos0[2]-self.zMin) / self.dz)

        if (binx >=0 and binx < self.xbins and biny >=0 and biny < self.ybins and binz >=0 and binz < self.zbins):

          lstats["iidXY"][binx, biny] += 1
          lstats["iidXZ"][binx, binz] += 1
          lstats["iidYZ"][biny, binz] += 1


    # how do eroded atoms affect the statistics?
    if edata != None:
      for entry in edata:

        # extract stored data
        sid  = entry[0]	# species ID
        pos0 = entry[1]	# initial position of sputtered atom

        # 2D projections of distributions
        binx = np.floor((pos0[0]-self.xMin) / self.dx)
        biny = np.floor((pos0[1]-self.yMin) / self.dy)
        binz = np.floor((pos0[2]-self.zMin) / self.dz)

        if (binx >=0 and binx < self.xbins and biny >=0 and biny < self.ybins and binz >=0 and binz < self.zbins):

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

        # 2D projections of distributions
        binx = np.floor((pos0[0]-self.xMin) / self.dx)
        biny = np.floor((pos0[1]-self.yMin) / self.dy)
        binz = np.floor((pos0[2]-self.zMin) / self.dz)

        if (binx >=0 and binx < self.xbins and biny >=0 and biny < self.ybins and binz >=0 and binz < self.zbins):

          lstats["rvdXY"][sid, binx, biny] += 1
          lstats["rvdXZ"][sid, binx, binz] += 1
          lstats["rvdYZ"][sid, biny, binz] += 1

          lstats["rddXYx"][sid, binx, biny] += dpos[0]
          lstats["rddXYy"][sid, binx, biny] += dpos[1]
          lstats["rddXZx"][sid, binx, binz] += dpos[0]
          lstats["rddXZz"][sid, binx, binz] += dpos[2]
          lstats["rddYZy"][sid, biny, binz] += dpos[1]
          lstats["rddYZz"][sid, biny, binz] += dpos[2]

        # 2D projections of distributions
        binx = np.floor((pos1[0]-self.xMin) / self.dx)
        biny = np.floor((pos1[1]-self.yMin) / self.dy)
        binz = np.floor((pos1[2]-self.zMin) / self.dz)

        if (binx >=0 and binx < self.xbins and biny >=0 and biny < self.ybins and binz >=0 and binz < self.zbins):

          lstats["ridXY"][sid, binx, biny] += 1
          lstats["ridXZ"][sid, binx, binz] += 1
          lstats["ridYZ"][sid, biny, binz] += 1


      # and update 2D projected distributions (also total of all species in position '0')
      for part1 in ["iid", "evd", "rvd", "rid"]:
        for part2 in ["XY", "XZ", "YZ"]:
          key = "%s%s" % (part1, part2)
          lstats[key][0] = np.sum(lstats[key][1:], axis=0)

      # and update 2D projected distributions (also total of all species in position '0')
      for part1 in ["rdd"]:
        for part2 in ["XYx", "XYy", "XZx", "XZz", "YZy", "YZz"]:
          key = "%s%s" % (part1, part2)
          lstats[key][0] = np.sum(lstats[key][1:], axis=0)


    return lstats





