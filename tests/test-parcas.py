# basic system and numerical packages
import sys
import numpy as np
import matplotlib.pyplot as plt

# pycraters basics
import pycraters as pc
import pycraters.IO as io
import pycraters.helpers as help

# build wrapper and params
exec_location = None
wrapper = pc.create_wrapper("PARCAS", exec_location)
#wrapper.set_parameter("endfile_exp", "00?.xyz")
params  = wrapper.create_parameters()

#basic parameter setup
params.target = [["Si", 1.0]]
params.beam = "Ar"
params.energy = 250
params.impacts = 300
params.angle = None

# do the simulations
angles = np.linspace(0,80,9)
for aa in angles:
  print "running angle %d" % (aa)
  params.angle = aa
  wrapper.go(params)

finedeg = np.linspace(0,90,91)
params.angle = None
fits   = help.linked_PDE_coefficients_1D(wrapper, params, angles, finedeg)
myplot = help.plot_single_flat_angle_dependence_summary(fits[0])
plt.show()
exit()


# plot graph of M0
m0e_avg = io.array_range( './', params, 'angle', angles, 'evdM0_avg' )
m0e_std = io.array_range( './', params, 'angle', angles, 'evdM0_std' ) 
yval = [ a[0] for a in m0e_avg ]
yerr = [ b[0]/np.sqrt(params.impacts) * 1.97 for b in m0e_std ]
plt.figure(1)
plt.errorbar(angles, yval, yerr=yerr, label='m0e')
plt.show()
plt.savefig('%s.svg' %(params.fname()))
