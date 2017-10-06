'''
Updated on 2014-03-25 to include execline in the construction of SDTRIMSP_wrapper()
'''

# Standard Imports
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Import tools which form libCraters
from libCraters import *
from wrappers.libSDTrimSP import *
from libPlotTools import *

# location of data
data_path = './'

# instantiate the solver
SDTRIM_location = sys.argv[1]
foo = SDTrimSP_Wrapper(SDTRIM_location)

# Base parameter setttings
params = ParameterBlob()
params.target    = [["Si", 1.0]]
params.beam 	 = None		# will be iterated over
params.energy	 = None		# will be iterated over
params.angle	 = 0.0
params.edispl	 = 5.0
params.impacts   = 1000


# things to iterate over
beams    = ['Ne']
energies = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]

# okay, here we go
for bb in beams:

  print "simulations for %s\n" % (bb)
  params.beam = bb

  # do the simulations
  for ee in energies:
    print "%d eV ..." % (ee), 
    params.energy = ee
    foo.go(params)
    print "done"

  # get relevant data
  params.energy = None
  m0e_avg = np.array( array_range(data_path, params, "energy", energies, "m0e_avg") )
  m0e_std = np.array( array_range(data_path, params, "energy", energies, "m0e_std") ) / np.sqrt(params.impacts) * 1.97
  data = [a[0] for a in m0e_avg]
  errs = [b[0] for b in m0e_std]

  # plot the results
  plt.errorbar(np.log10(energies), data, yerr=errs, label=bb)

# finish and save the figure
plt.legend(loc=1)
plt.savefig('yields-vs-energy.svg')
plt.show()


