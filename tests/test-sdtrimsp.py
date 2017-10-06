# basic system and numerical packages
import sys
import numpy as np
import matplotlib.pyplot as plt

# libcraters basics
import pycraters.wrappers.SDTRIMSP as wrap
import pycraters.IO as io
import pycraters.helpers as help

# build wrapper and params
exec_location = sys.argv[1]
wrapper = wrap.SDTRIMSP_Wrapper(exec_location)
params  = wrap.SDTRIMSP_Parameters()

#basic parameter setup
params.beam = "Ar"
params.energy = 1000
params.angle = None
params.impacts = 1000
params.target = [["Si", 1.0]]

# do the simulations
angles = np.linspace(0,85,18)
for aa in angles:
  print "running angle %d" % (aa)
  params.angle = aa
  wrapper.go(params)

# plot graph of M0
m0e_avg = io.array_range( './', params, 'angle', angles, 'm0e_avg' )
m0e_std = io.array_range( './', params, 'angle', angles, 'm0e_std' ) 
yval = [ a[0] for a in m0e_avg ]
yerr = [ b[0]/np.sqrt(params.impacts) * 1.97 for b in m0e_std ]
plt.figure(1)
plt.errorbar(angles, yval, yerr=yerr, label='m0e')
plt.show()
plt.savefig('%s.svg' %(params.fname()))
