# Basic Imports
import sys
import os
import copy


# pylab
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Parameters, minimize

# PyCraters imports
import pycraters.wrappers.TRI3DST as wrap
import pycraters.IO as io
import pycraters.helpers as help
import pycraters.fits.polynomial as pol








# -----------------------------------------
#         simulation section
# -----------------------------------------



# get the executable location and create wrapper
path = './'
execline = sys.argv[1]
wrapper = wrap.TRI3DST_Wrapper(execline)


# define surface binding energy models for each compound
# Heats of sublimation -- from ???
dHs_Ga = 2.82
dHs_Sb = 2.74
dHs_As = 3.13
dHs_P  = 3.28

# Heats of formation -- from Lev I Berger, "Semiconductor Materials" CRC Press 1997
dHf_GaSb = 0.46
dHf_GaAs = 0.87
dHf_GaP  = 0.75

AB_GaSb = 0.5*(dHs_Ga + dHs_Sb) + dHf_GaSb
AB_GaAs = 0.5*(dHs_Ga + dHs_As) + dHf_GaAs
AB_GaP  = 0.5*(dHs_Ga + dHs_P) + dHf_GaP

SBVgasb = np.array([[0, 0, 0],[0, dHs_Ga, AB_GaSb],[0, AB_GaSb, dHs_Sb]])
SBVgaas = np.array([[0, 0, 0],[0, dHs_Ga, AB_GaAs],[0, AB_GaAs, dHs_As]])
SBVgap  = np.array([[0, 0, 0],[0, dHs_Ga, AB_GaP ],[0, AB_GaP , dHs_P ]])



# do the simulations
angles = np.linspace(0,80,9)
finedeg = np.linspace(0, 90, 91)
powertops = np.linspace(0, 20, 21)
#energies = 100*10.**(powertops/10.)

energies = [1000]


concentrations = np.linspace(0.05, 0.95, 19)

Aconcs = concentrations
Bconcs = 1-concentrations
fineconcs = np.linspace(0.0, 1.0, 51)
Afconcs = fineconcs
Bfconcs = 1-fineconcs


gp3 = "Ga"
gp5 = "Sb"
compound = "GaSb"
bb  = "Ne"
for ee in energies:

  YAalist = []
  YBalist = []
  CAalist = []
  CBalist = []

  YAslist = []
  YBslist = []
  CAslist = []
  CBslist = []


  for cc in concentrations:

    print "running simulations for %s%f %s%f" % (gp3, cc, gp5, 1-cc)

    params  = wrap.TRI3DST_Parameters()
    params.target  = [[gp3, cc],[gp5, 1.0-cc]]
    params.beam    = bb
    params.energy  = ee
    params.angle   = None
    params.impacts = 1000

    results = help.linked_PDE_coefficients_1D(wrapper, params, angles, finedeg)

    #YAa0 = results[1]['m0e_vals'][0]
    #YBa0 = results[2]['m0e_vals'][0]
    #CAa0 = results[1]['sx_coeffs'][0]
    #CBa0 = results[2]['sx_coeffs'][0]

    #YAs0 = None
    #YBs0 = None
    #CAs0 = None
    #CBs0 = None

    YAa0 = results[1]['m0e_fit']["c0"].value
    YBa0 = results[2]['m0e_fit']["c0"].value 
    CAa0 = results[1]['m1e_fit']["c1"].value / 2.0 + results[1]['m1r_fit']["c1"].value
    CBa0 = results[2]['m1e_fit']["c1"].value / 2.0 + results[2]['m1r_fit']["c1"].value

    YAs0 = results[1]['m0e_fit']["c0"].stderr
    YBs0 = results[2]['m0e_fit']["c0"].stderr 
    CAs0 = results[1]['m1e_fit']["c1"].stderr / 2.0  + results[1]['m1r_fit']["c1"].stderr
    CBs0 = results[2]['m1e_fit']["c1"].stderr / 2.0 + results[2]['m1r_fit']["c1"].stderr

    YAalist.append(0.1 * YAa0)
    YBalist.append(0.1 * YBa0)
    CAalist.append(0.1 * CAa0)
    CBalist.append(0.1 * CBa0)

    YAslist.append(0.1 * YAs0)
    YBslist.append(0.1 * YBs0)
    CAslist.append(0.1 * CAs0)
    CBslist.append(0.1 * CBs0)


  YAvals = np.array(YAalist)
  YBvals = np.array(YBalist)
  CAvals = np.array(CAalist)
  CBvals = np.array(CBalist)

  YAstds = np.array(YAslist)
  YBstds = np.array(YBslist)
  CAstds = np.array(CAslist)
  CBstds = np.array(CBslist)

  #YAstds = None
  #YBstds = None
  #CAstds = None
  #CBstds = None

  fparams = Parameters()
  fparams.add("c1", value=0.0)
  fparams.add("c2", value=0.0)
  #fparams.add("c3", value=0.0)

  fitYA = minimize(pol.polynomial_residual, copy.deepcopy(fparams), args=(Aconcs, YAvals, YAstds))
  fitYB = minimize(pol.polynomial_residual, copy.deepcopy(fparams), args=(Bconcs, YBvals, YBstds))
  fitCA = minimize(pol.polynomial_residual, copy.deepcopy(fparams), args=(Aconcs, CAvals, CAstds))
  fitCB = minimize(pol.polynomial_residual, copy.deepcopy(fparams), args=(Bconcs, CBvals, CBstds))

  fittedYA = pol.polynomial(fitYA.params, Afconcs)	# note: this is a function of *cA*
  fittedYB = pol.polynomial(fitYB.params, Bfconcs)	# note: this is a function of *cB*
  fittedCA = pol.polynomial(fitCA.params, Afconcs)	# note: this is a function of *cA*
  fittedCB = pol.polynomial(fitCB.params, Bfconcs)	# note: this is a function of *cB*

  fittedYAp = pol.polynomial_D1(fitYA.params, Afconcs)	# note: this is a function of *cA*
  fittedYBp = pol.polynomial_D1(fitYB.params, Bfconcs)	# note: this is a function of *cB*

  cbulk = [0.5, 0.5]
  film_height = 5.0

  coeffA  = - fittedYAp + fittedYBp
  coeffC  =   fittedCA  + fittedCB
  coeffAp = ( - cbulk[1]*fittedYAp - cbulk[0]*fittedYBp ) / film_height
  coeffCp = (   cbulk[1]*fittedCA  - cbulk[0]*fittedCB  ) / film_height
  Ceff    = coeffC - coeffA * coeffCp / coeffAp







  f1  = plt.figure(1, figsize=(8,8))
  ax0 = f1.add_subplot(211)
  plt.errorbar(Aconcs, -YAvals, yerr=YAstds, fmt='go', linewidth=2, label="A")
  #plt.plot(Aconcs, -YAvals, 'gs', linewidth=2, label="A")
  plt.plot(Afconcs, -fittedYA, 'g-', linewidth=2)
  plt.errorbar(Aconcs, -YBvals, yerr=YBstds, fmt='bo', linewidth=2, label="B")
  #plt.plot(Aconcs, -YBvals, 'bs', linewidth=2, label="B")
  plt.plot(Afconcs, -fittedYB, 'b-', linewidth=2)
  plt.plot(Afconcs, -fittedYA - fittedYB, 'k--', linewidth=2, label="A+B")
  plt.legend(loc='lower center', prop={'size':8})
  plt.title("YA vs YB:  %s --> %s" % (bb, compound))
  plt.xlabel("%s concentration" % (gp3))
  plt.ylabel("value [atom / ion]")


  ax0 = f1.add_subplot(212)
  plt.errorbar(Aconcs, CAvals, yerr=CAstds, fmt='go', linewidth=2, label="A")
  #plt.plot(Aconcs, CAvals, 'gs', linewidth=2, label="A")
  plt.plot(Afconcs, fittedCA, 'g-', linewidth=2)
  plt.errorbar(Aconcs, CBvals, yerr=CBstds, fmt='bo', linewidth=2, label="B")
  #plt.plot(Aconcs, CBvals, 'bs', linewidth=2, label="B")
  plt.plot(Afconcs, fittedCB, 'b-', linewidth=2)
  plt.plot(Afconcs, fittedCA + fittedCB, 'k--', linewidth=2, label="A+B")
  plt.legend(loc='lower center', prop={'size':8})
  plt.title("CA vs CB:  %s --> %s" % (bb, compound))
  plt.xlabel("%s concentration" % (gp3))
  plt.ylabel("value [atom * nm / ion]")


  plt.tight_layout()
  plt.savefig("YAvsYB-%seV.png" % (ee))
  plt.close()



  f1  = plt.figure(1, figsize=(8,6))
  ax1 = f1.add_subplot(221)
  plt.plot(Afconcs, coeffA, 'b-', label="A", linewidth=2)
  plt.legend(loc="best", prop={'size':8})

  ax1 = f1.add_subplot(222)
  plt.plot(Afconcs, coeffAp, 'b-', label="Ap", linewidth=2)
  plt.legend(loc="best", prop={'size':8})

  ax1 = f1.add_subplot(223)
  plt.plot(Afconcs, coeffC, 'b-', label="C", linewidth=2)
  plt.legend(loc="best", prop={'size':8})

  ax1 = f1.add_subplot(224)
  plt.plot(Afconcs, coeffCp, 'b-', label="Cp", linewidth=2)
  plt.legend(loc="best", prop={'size':8})

  plt.tight_layout()
  plt.savefig("4coeffs-%seV.png" % (ee))
  plt.close()


  f1  = plt.figure(1, figsize=(6,4))
  ax1 = f1.add_subplot(111)
  plt.plot(Afconcs, coeffC, 'bo-', label="C", linewidth=2)
  plt.plot(Afconcs, Ceff, 'ro-', label="Ceff", linewidth=2)
  plt.legend(loc=1)
  plt.title("C vs Ceff:  %s --> %s" % (bb, compound))
  plt.xlabel("%s concentration" % (gp3))
  plt.ylabel("value [atom * nm / ion]")
  plt.tight_layout()
  plt.savefig("cvsceff-%seV.png" % (ee))
  plt.close()
    








