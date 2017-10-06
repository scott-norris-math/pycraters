# basic numerical packages
import sys
import numpy as np
import matplotlib.pyplot as plt

# libcraters basics
import pycraters.wrappers.TRI3DST as wrap
import pycraters.IO as io
import pycraters.helpers as help

import pylab
def __init__():
  params = {'axes.labelsize': 10,
            'text.fontsize': 8,
            'legend.fontsize': 6,
            'legend.labelspacing':0.25,
            'xtick.labelsize': 7,
            'ytick.labelsize': 7}#,
  #         'figure.figsize': fig_size}
  pylab.rcParams.update(params)





# Get executable location and create wrapper
exec_location = sys.argv[1]
foo = wrap.TRI3DST_Wrapper(exec_location)
params = wrap.TRI3DST_Parameters()

# set parameter values
params.target  = [["Si", 1.0]]
params.beam    = "Ar"
params.energy  = 1000
params.impacts = 10000
angles   = np.linspace(0, 85, 18)

finedeg = np.linspace(0,90,91)

dks = [ .004 ]

for dk in dks:

  thefits = help.linked_PDE_coefficients_2D(foo, params, angles, finedeg, dk = dk)
  theplot = help.plot_single_curved_angle_dependence_summary(thefits[1])
  plt.savefig('curvature-dependence-dk-%f.svg' % (dk))
  plt.close()

  rvals = thefits[1]

  # plot the data and the fit
  plt.figure(1, figsize=(12,4))

  plt.subplot(131)
  plt.errorbar(rvals["angles"], rvals["m0e_avg"],    yerr=rvals["m0e_err"],    fmt='bs', label='flat')
  plt.errorbar(rvals["angles"], rvals["m0ek1p_avg"], yerr=rvals["m0ek1p_err"], fmt='gs', label='k11=%0.3f'%(dk))
  plt.errorbar(rvals["angles"], rvals["m0ek1m_avg"], yerr=rvals["m0ek1m_err"], fmt='rs', label='k11=%0.3f'%(-dk))
  plt.plot(rvals["finedeg"], rvals["m0e_vals"], 'b-', linewidth=2)
  plt.plot(rvals["finedeg"], rvals["m0ek1p_vals"], 'g-', linewidth=2)
  plt.plot(rvals["finedeg"], rvals["m0ek1m_vals"], 'r-', linewidth=2)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$M^{(0)}_{\mathsf{eros.}} \left( \theta \right)$')
  plt.legend(loc=3)
  plt.title('Sputter Yield')

  plt.subplot(132)
  plt.errorbar(rvals["angles"], rvals["m0e_avg"],    yerr=rvals["m0e_err"],    fmt='bs', label='flat')
  plt.errorbar(rvals["angles"], rvals["m0ek2p_avg"], yerr=rvals["m0ek2p_err"], fmt='gs', label='k22=%0.3f'%(dk))
  plt.errorbar(rvals["angles"], rvals["m0ek2m_avg"], yerr=rvals["m0ek2m_err"], fmt='rs', label='k22=%0.3f'%(-dk))
  plt.plot(rvals["finedeg"], rvals["m0e_vals"], 'b-', linewidth=2)
  plt.plot(rvals["finedeg"], rvals["m0ek2p_vals"], 'g-', linewidth=2)
  plt.plot(rvals["finedeg"], rvals["m0ek2m_vals"], 'r-', linewidth=2)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$M^{(0)}_{\mathsf{eros.}} \left( \theta \right)$')
  plt.legend(loc=3)
  plt.title('Sputter Yield')


  plt.subplot(133)
  plt.errorbar(rvals["angles"], rvals["m1e_avg"], yerr=rvals["m1e_err"], fmt='rs', label='eros.')
  plt.errorbar(rvals["angles"], rvals["m1r_avg"], yerr=rvals["m1r_err"], fmt='bs', label='redist.')
  plt.plot(rvals["finedeg"], rvals["m1e_vals"], 'r-', linewidth=2)
  plt.plot(rvals["finedeg"], rvals["m1r_vals"], 'b-', linewidth=2)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$M^{(1)} \left( \theta \right)$')
  plt.legend(loc=3)
  plt.title('First Moments')

  plt.tight_layout()
  plt.savefig('curvature-dependence-moments-dk-%f.svg' % (dk))


  # plot the data and the fit
  plt.figure(2, figsize=(8,4))

  plt.subplot(121)
  plt.plot(rvals["finedeg"], rvals["sxe_coeffs"], 'r-', linewidth=2, label=r"eros.")
  plt.plot(rvals["finedeg"], rvals["sxr_coeffs"], 'b-', linewidth=2, label=r"redist.")
  plt.plot(rvals["finedeg"], rvals["sxc_coeffs"], 'g-', linewidth=2, label=r"curv.")
  plt.plot(rvals["finedeg"], rvals["sx_coeffs"], 'k--', linewidth=2, label=r"total")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{X} \left( \theta \right)$')
  plt.legend(loc=3)
  plt.title('SX components')
  plt.ylim(-100,100)

  plt.subplot(122)
  plt.plot(rvals["finedeg"], rvals["sye_coeffs"], 'r-', linewidth=2, label=r"eros.")
  plt.plot(rvals["finedeg"], rvals["syr_coeffs"], 'b-', linewidth=2, label=r"redist.")
  plt.plot(rvals["finedeg"], rvals["syc_coeffs"], 'g-', linewidth=2, label=r"curv.")
  plt.plot(rvals["finedeg"], rvals["sy_coeffs"], 'k--', linewidth=2, label=r"total")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{Y} \left( \theta \right)$')
  plt.legend(loc=3)
  plt.title('SY components')
  plt.ylim(-100,100)

  plt.tight_layout()
  plt.savefig('curvature-dependence-coefficients-dk-%f.svg' % (dk))


  # plot the data and the fit
  plt.figure(3, figsize=(8,4))
  
  plt.subplot(121)
  plt.plot(rvals["finedeg"], rvals["sx_coeffs"], 'c-', linewidth=2, label="SX")
  plt.plot(rvals["finedeg"], rvals["sy_coeffs"], 'm-', linewidth=2, label="SY")
  plt.plot([0, 90], [0, 0], 'k--', linewidth=1)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{X,Y} \left( \theta \right)$')
  plt.legend(loc=3)
  plt.title('SX and SY (with curvature)')
  plt.ylim(-60,60)

  plt.subplot(122)
  plt.plot(rvals["finedeg"], rvals["sxe_coeffs"] + rvals["sxr_coeffs"], 'c-', linewidth=2, label="SX")
  plt.plot(rvals["finedeg"], rvals["sye_coeffs"] + rvals["syr_coeffs"], 'm-', linewidth=2, label="SY")
  plt.plot([0, 90], [0, 0], 'k--', linewidth=1)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{X,Y} \left( \theta \right)$')
  plt.legend(loc=3)
  plt.title('SX and SY (no curvature)')
  plt.ylim(-60,60)

  plt.tight_layout()
  plt.savefig('curvature-dependence-comparisons-dk-%f.svg' % (dk))


plt.show()

