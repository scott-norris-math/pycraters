# basic system and numerical packages
import sys
import numpy as np
import matplotlib.pyplot as plt

# libcraters basics
import pycraters.wrappers.TRI3DST as wrap
import pycraters.IO as io
import pycraters.helpers as help

# build wrapper and params
exec_location = sys.argv[1]
wrapper = wrap.TRI3DST_Wrapper(exec_location)
params  = wrap.TRI3DST_Parameters()

# get some more arguments
target_species = sys.argv[2]
ion_species = sys.argv[3]


#basic parameter setup
params.target = [[target_species, 1.0]]
params.beam = ion_species
params.energy = None
params.angle = None
params.impacts = 1000

# build the wrapper





# do the simulations
angles = np.linspace(0,85,18)
powertops = np.linspace(0, 20, 21)
energies = 100*10.**(powertops/10.)

tlistlist = []
finedeg = np.linspace(0, 90, 91)

for ee in energies:
  print "running energy %4f" % (ee)
  for aa in angles:
    params.energy = ee
    params.angle = aa
    wrapper.go(params)

  fits  = help.linked_PDE_coefficients_1D(foo, params, angles, finedeg)

  rvals = fits[1]
  plt.close()
  plt.figure(figsize=(6.5, 4.333))
  plt.plot(rvals["finedeg"], rvals["sx_coeffs"], 'c-', linewidth=3, label="SX")
  plt.plot(rvals["finedeg"], rvals["sy_coeffs"], 'm-', linewidth=3, label="SY")
  plt.plot([0, 90], [0, 0], 'k--', linewidth=1)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{X,Y} \left( \theta \right)$')
  plt.title('%5d eV' % (np.round(ee)))
  plt.legend(loc=3)
  plt.xlim(0,90)
  yval = ee * 3. / 100.
  plt.ylim(-yval, yval)

  filename = params.fname(['target', 'beam', 'energy'])
  plt.savefig('%s-coeffs.svg' % (filename))
  plt.savefig('%s-coeffs.pdf' % (filename))
  plt.savefig('%s-coeffs.png' % (filename))

  tdata = help.find_pattern_transitions(ee, angles, finedeg, fits[1]["sx_coeffs"], fits[1]["sy_coeffs"])
  tlistlist.append(tdata)


# plot the phase diagram
help.plot_energy_angle_phase_diagram(tlistlist)
plt.title(r'%s $\to$ %s' % (ion_species, target_species))
plt.savefig('%s%s-phase-diagram.svg' % (ion_species, target_species))
plt.show()

