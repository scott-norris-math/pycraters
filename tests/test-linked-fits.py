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


# analyze data
finedeg = np.linspace(0, 90, 91)
thefits = help.linked_PDE_coefficients(foo, params, angles, finedeg)

# plot data
theplot = help.plot_single_flat_angle_dependence_summary(thefits[1])
plt.savefig('%s.svg' %(params.fname()))
plt.show()

