import os
import sys
import math
import shelve
import numpy as np
import matplotlib.pyplot as plt

import pycraters.IO as io		# allows data extraction from files.
import fittools				# external library for data fitting
import pycraters.fits.poly_yamamura_1D as pyam
import lmfit as lmf
import pycraters.fits.NestedLmfit as fitter




import pylab
def __init__():
  params = {'axes.labelsize': 10,
            'text.fontsize': 10,
            'legend.fontsize': 7,
            'legend.labelspacing':0.33,
            'xtick.labelsize': 7,
            'ytick.labelsize': 7}#,
  #         'figure.figsize': fig_size}
  pylab.rcParams.update(params)




def ensure_data(params, geom, execline, path):
  targetfile = "%s%s.moms" % (path, params.fname())
  if (os.path.isfile(targetfile) == False):
    print "No existing data found under %s. Proceeding with simulation" % (targetfile)
    sdtrimsp_run(execline, params, geom)
    sdtrimsp_get_statistics(params, geom)
  else:
    print "%s found." % (targetfile)






# -------------------------------------------------
#      find_equilibrium_concentration()
# -------------------------------------------------

def find_equilibrium_concentration(params, geom, execline, targets, path='./', funcparams=[0]):

  '''
  For a given set of parameters, this function simulates sputtering over a range of concentrations,
  and identifies that concentration at which the sputter rate is exactly equal to the supply rate;
  i.e., the steady concentration.  Identifying this concentration is necessary for estimation of
  parameters that appear in two-component systems.
  '''



  if (len(funcparams) == 3):
    alpha = funcparams[0]
    beta = funcparams[1]
    gamma = funcparams[2]
  for tt in targets:
    params.target = tt
    if (params.funcname and len(funcparams)==3):
      gm = gamma_max(tt[1][1], alpha, beta)
      if gm <= 1.0:
        gamma *= gm
      params.cfunc = lambda depth : phistar_list(depth, tt[1][1], alpha, beta, gamma)
    ensure_data(params, geom, execline, path)

  m0e_avgs = io.array_range(path, params, "target", targets, "m0e_avg")

  n=0 #find the intersection point, where the index [n] refers to the point after intersection, and at phi[n], a is less than b
  if (m0e_avgs[n][0] > m0e_avgs[n][1]):
    a = 0
    b = 1
  else:
    a = 1
    b = 0
  while (m0e_avgs[n][a] > m0e_avgs[n][b]):
    n += 1
  #currently the above logic is only necessary to set up the while loop
  #targets has form [ [["Ga", phi],["Sb", 1-phi]], etc ]
  #for a certain range of phi. We can access these values directly through
  #targets[n] --> [["Ga", phi], ["Sb", 1-phi]]
  #targets[n][0][1]

  phi0 = targets[n-1][0][1]
  phi1 = targets[n][0][1]

  equilibrium_concentration = (phi1*(m0e_avgs[n-1][a] - m0e_avgs[n-1][b]) - phi0*(m0e_avgs[n][a] - m0e_avgs[n][b])) / (m0e_avgs[n-1][a] + m0e_avgs[n][b] - m0e_avgs[n-1][b] - m0e_avgs[n][a])

  return equilibrium_concentration





#def calculate_energy_angle_phase_diagram(wrapper, params, energies, angles, fitmethod, guess, finedeg=None, path='./', curvature_effects_dk=None):
#  for ee in energies:
#    for aa in angles:
#      params.energy = ee
#      params.angle  = aa
#      wrapper.go(params)
#  for ee in energies:






# -------------------------------------------------
#      find_pattern_transitions()
# -------------------------------------------------

def find_pattern_transitions(energy, angles, finedeg, sx_values, sy_values):

  # allocate some storage
  transitions = dict()
  transitions["energy"] = energy
  transitions["angles"] = angles

  # identify the pattern *near* theta == 0
  oldangle  = finedeg[1]
  oldvector = np.array([0, sx_values[1], sy_values[1]])
  oldpattern   = np.argmin(oldvector)

  # sweep through angles
  for kk in range(2, len(finedeg)-1):

    # identify the pattern at each new angle
    newangle  = finedeg[kk]
    newvector = np.array([0, sx_values[kk], sy_values[kk]])
    newpattern   = np.argmin(newvector)

    # if pattern changes ...
    if (newpattern != oldpattern):

      # data from old pattern
      x0 = oldangle
      y0 = oldvector[oldpattern]
      z0 = oldvector[newpattern]

      # data from new pattern
      x1 = newangle
      y1 = newvector[oldpattern]
      z1 = newvector[newpattern]

      # calculate transition angle (linear approx.)
      dxstar = - (z0-y0)*(x1-x0) / ((z1-y1)-(z0-y0))
      if dxstar < 0:  print "WARNING: dxstar was found to be %f" % (dxstar)
      tangle = oldangle + dxstar
      ttype = "%s%s" % (oldpattern, newpattern)

      if ttype not in transitions:  transitions[ttype] = tangle


    # cycle pattern information
    oldangle = newangle
    oldvector = newvector
    oldpattern = newpattern

  return transitions





# -------------------------------------------------
#      plot_energy_angle_phase_diagram()
# -------------------------------------------------

def plot_energy_angle_phase_diagram(tlist):


  thefigure = plt.figure(figsize=(6.5, 4.333))
  for t in tlist:
    anglesim = t["angles"]
    energysim = t["energy"]*np.ones(len(anglesim)) 
    plt.semilogy(anglesim, energysim, 'ko', markersize=0.5)

  # set up a list of colors
  colors = dict()
  colors['01'] = 'b'
  colors['12'] = 'r'

  # set up a dictionary of descriptors
  descriptions = dict()
  descriptions['01'] = r'smooth $\to$ parallel'
  descriptions['12'] = r'parallel $\to$ perp.'
  
  # sweep through all possible transitions
  for ii in range(3):
    for jj in range(3):

      # string describing transition type
      ttype = "%s%s" % (ii,jj)

      # get all energies for which this type is present
      energies = [ t["energy"] for t in tlist if ttype in t ]

      # get the corresponding angles at which this type occurs
      angles   = [ t[ttype] for t in tlist if ttype in t ]

      # plot (WEAK IMPLEMENTATION -- does not distinguish between multiple types.

      cstring = colors[ttype] if ttype in colors else ''
      dstring = descriptions[ttype] if ttype in descriptions else ''

      plt.semilogy(angles, energies, '%ss-' % (cstring), linewidth=3, markersize=6)

  plt.xlim(0,90)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'energy $E$')
  #plt.legend(loc=4)
  plt.tight_layout()

  return thefigure







# --------------------------------------------------
#    fit_moments_1D()
# --------------------------------------------------


def fit_moments_1D(params, angles, fitmethod=None, fitguess=None, path='./'):

  '''
  This function takes a set of parametric dependencies and a *range of angles,*
  and fits the zeroth, first erosive, and first redistributive moments over 
  those angles to a Yamamura-modified polynomial.

  The function also stores the fits to a file at the end of the function, and
  checks for the existence of that file at the beginning of the function.  This
  allows scripts to be re-run during testing without re-calculating the fits.
  However, it may be that this function is not the best place for this 
  functionality.

  In principle different fit methods can be provided, but this is not currently
  implemented.
  '''



  # see if fit values already exist
  fname = params.fname(['target', 'beam', 'energy'])
  if os.path.isfile('%s.momfits' % (fname)):
    f = shelve.open('%s.momfits' % (fname))
    results_list = list(f['results_list'])
    f.close()
    return results_list


  # determine number of species, and number of impacts
  species = len(params.target)
  impacts = params.impacts
  finedeg = np.linspace(0,90,91)

  # configure storage for pure or mutli-species results
  extra_species = 0
  if species > 1:  extra_species = 1
  total_species = species + extra_species

  # load moment data
  m0e_avg = -np.array( io.array_range( './', params, 'angle', angles, 'evdM0_avg' ))
  m0e_std =  np.array( io.array_range( './', params, 'angle', angles, 'evdM0_std' ))
  m1e_avg = -np.array( io.array_range( './', params, 'angle', angles, 'evdM1_avg' ))
  m1e_std =  np.array( io.array_range( './', params, 'angle', angles, 'evdM1_std' ))
  m1r_avg = np.array(io.array_range( './', params, 'angle', angles, 'rddM1_avg' ))
  m1r_std = np.array(io.array_range( './', params, 'angle', angles, 'rddM1_std' ))

  # create even fit parameter
  fp_even = lmf.Parameters()
  fp_even.add("p" , value=0.25, min=0.0)
  fp_even.add("c0", value=0.00)
  fp_even.add("c2", value=0.00)
  fp_even.add("c4", value=0.00)

  # create odd fit parameter
  fp_odd  = lmf.Parameters()
  fp_odd.add("p" , value=0.25, min=0.0)
  fp_odd.add("c1", value=0.00)
  fp_odd.add("c3", value=0.00)
  fp_odd.add("c5", value=0.00)

  # greate a global fit parameter
  fp_global  = lmf.Parameters()
  fp_global.add("p", value=0.25, min=0.0)
  

  # allocate storage
  results_list = []

  # iterate through species
  for kk in range(total_species):

    # storage for returned values
    rvals = dict()
    rvals["angles"] = angles
    rvals["finedeg"] = finedeg

    # storage for nodelist
    nodelist = []

    # M0e data
    m0e_vals = np.array( [ a[kk] for a in m0e_avg ] )
    m0e_errs = np.array( [ b[kk]/np.sqrt(params.impacts) * 1.97 for b in m0e_std ] )
    rvals["m0e_dats"] = m0e_vals
    rvals["m0e_errs"] = m0e_errs
    m0e_node = fitter.FitNode( copy.deepcopy(fp_even), pyam.poly_yamamura_1D, angles, m0e_vals, m0e_errs )
    nodelist.append(m0e_node)

    # M1e data
    m1e_vals = np.array( [ a[kk][0] for a in m1e_avg ] )
    m1e_errs = np.array( [ b[kk][0]/np.sqrt(params.impacts) * 1.97 for b in m1e_std ] )
    rvals["m1e_dats"] = m1e_vals
    rvals["m1e_errs"] = m1e_errs
    m1e_node = fitter.FitNode( copy.deepcopy(fp_odd), pyam.poly_yamamura_1D, angles, m1e_vals, m1e_errs )
    nodelist.append(m1e_node)

    # M1r data
    m1r_vals = np.array( [ a[kk][0] for a in m1r_avg ] )
    m1r_errs = np.array( [ b[kk][0]/np.sqrt(params.impacts) * 1.97 for b in m1r_std ] )
    rvals["m1r_dats"] = m1r_vals
    rvals["m1r_errs"] = m1r_errs
    m1r_node = fitter.FitNode( copy.deepcopy(fp_odd), pyam.poly_yamamura_1D, angles, m1r_vals, m1r_errs )
    nodelist.append(m1r_node)

    # SuperNode and Fitting
    snode   = fitter.FitSuperNode(nodelist, fp_global)
    result  = snode.fit()
    fparams = result.params

    rvals['m0e_fit'] = m0e_node.result.params
    rvals['m1e_fit'] = m1e_node.result.params
    rvals['m1r_fit'] = m1r_node.result.params

    results_list.append(rvals)

  # populate a shelf
  fname = params.fname(['target', 'beam', 'energy'])
  f = shelve.open('%s.momfits' % (fname))
  f['results_list'] = results_list
  f.close()
  return results_list








# --------------------------------------------------
#    linked_PDE_coefficients_1D()
# --------------------------------------------------



def linked_PDE_coefficients_1D(wrapper, params, angles, finedeg, fitmethod=None, fitguess=None, path=None):

  '''
  This function calculates *flat-target* estimates of the values of the coefficients
  through second order in the linearized equation of motion.  It provides the total
  values of these coefficients, and also the contribution from each component (if
  the target is a multi-component target).  
  
  TODO:  Provide the *total* value within the zeroth entry field.
  TODO:  Modify SX and SY to take into account estimate of curvature correction.
  TODO:  Also calculate the coefficient of the first-order hx term.
  '''

  # change to the path
  if path:
    oldpath = os.getcwd()
    os.chdir( path )

  # run simulations if needed
  for aa in angles:
    params.angle = aa
    wrapper.go(params)

  # functions we will use for fitting
  fitfuncD0 = pyam.poly_yamamura
  fitfuncD1 = pyam.poly_yamamura_derivative

  # get fitted values
  results_list = fit_moments_1D(params, angles)

  # calculate PDE coefficients in 1D
  for rvals in results_list:

    rvals["finedeg"] = finedeg
    finerad = finedeg * np.pi / 180.

    # record fitted values of the functions, themselves
    rvals["m0e_vals"] = fitfuncD0(rvals["m0e_fit"], finedeg, evenopts)
    rvals["m1e_vals"] = fitfuncD0(rvals["m1e_fit"], finedeg,  oddopts)
    rvals["m1r_vals"] = fitfuncD0(rvals["m1r_fit"], finedeg,  oddopts)
    rvals["m1_vals"]  = rvals["m1e_vals"] + rvals["m1r_vals"]

    # record fitted values of the derivatives of M1
    rvals["m0ep_vals"] = fitfuncD1(rvals["m0e_fit"], finedeg, evenopts)
    rvals["m1ep_vals"] = fitfuncD1(rvals["m1e_fit"], finedeg, oddopts )
    rvals["m1rp_vals"] = fitfuncD1(rvals["m1r_fit"], finedeg, oddopts )
    rvals["m1p_vals"]  = rvals["m1ep_vals"] + rvals["m1rp_vals"]

    # construct SX, using the approximation for explicit curvature dependence
    rvals["sxe_coeffs"] =   np.cos(finerad) * rvals["m1ep_vals"] - np.sin(finerad) * rvals["m1e_vals"]
    rvals["sxr_coeffs"] =   np.cos(finerad) * rvals["m1rp_vals"] - np.sin(finerad) * rvals["m1r_vals"]
    rvals["sxc_coeffs"] = - np.cos(finerad) * rvals["m1ep_vals"] / 2.0
    rvals["sx_coeffs"]  = rvals["sxe_coeffs"] + rvals["sxr_coeffs"] + rvals["sxc_coeffs"]

    # construct SY, using the approximation for explicit curvature dependence
    rvals["sye_coeffs"] =   np.cos(finerad)**2 / np.sin(finerad) * rvals["m1e_vals"] 
    rvals["syr_coeffs"] =   np.cos(finerad)**2 / np.sin(finerad) * rvals["m1r_vals"] 
    rvals["syc_coeffs"] = - np.cos(finerad)**2 / np.sin(finerad) * rvals["m1e_vals"] / 2.0
    rvals["sy_coeffs"]  = rvals["sye_coeffs"] + rvals["syr_coeffs"] + rvals["syc_coeffs"]

  # change to the path
  if path:
    os.chdir( oldpath )

  # construct the sum dictionary
  return results_list
  




# --------------------------------------------------
#    binary_PDE_coefficients_1D()
# --------------------------------------------------



def binary_PDE_coefficients_1D(wrapper, params, angles, finedeg, cbulk=None, film_height=None, fitmethod=None, fitguess=None, path=None, dc=None):

  '''
  This function calculates *flat-target* estimates of the values of the coefficients
  through second order in the linearized equation of motion.  It provides the total
  values of these coefficients, and also the contribution from each component (if
  the target is a multi-component target).  *FINALLY, IT USES THE CONTRIBUTIONS TO
  THE HEIGHT-FIELD EQUATION OF MOTION TO GENERATE COEFFICIENTS FOR THE COUPLED PDE
  EQUATIONS OF MOTION*
  
  TODO:  Also calculate the coefficient of the first-order hx term.
  TODO:  Document the units of the returned quantities
  TODO:  identify whether we need to provide additional physical properties to this function
  '''

  # make sure bulk values of the concentrations are provided
  if (cbulk == None):
    print "no bulk concentration value provided in binary_PDE_coefficients_1D!  exiting."
    quit()

  # make film height is provided
  if (film_height == None):
    print "no film height provided in binary_PDE_coefficients_1D!  exiting."
    quit()

  # get the values of the coefficients at the requested concentration
  #print "fitting center point"; sys.stdout.flush()
  results = linked_PDE_coefficients_1D(wrapper, params, angles, finedeg)
  SXeA = results[1]['sxe_coeffs']
  SXeB = results[2]['sxe_coeffs']
  SXrA = results[1]['sxr_coeffs']
  SXrB = results[2]['sxr_coeffs']
  SXcA = results[1]['sxc_coeffs']
  SXcB = results[2]['sxc_coeffs']
  SXA  = results[1]['sx_coeffs']
  SXB  = results[2]['sx_coeffs']

  # identify the difference type
  difftype = "center"
  diffdenom = 2*dc

  Aconc = params.target[0][1]
  Bconc = params.target[1][1]

  YA0 = results[1]['m0e_vals']
  YB0 = results[2]['m0e_vals']

  if Aconc - dc < 0.0:
    # if the left neighbor has negative Aconc, then we are at the left edge
    # and must use a rightward finite difference
    difftype = "right"
    diffdenom = dc
    YAm = results[1]['m0e_vals']
    YBp = results[2]['m0e_vals']
  else:
    # otherwise, we can get the values of the coefficients
    # at a smaller concentration of species A
    params.target[0][1] = Aconc - dc
    params.target[1][1] = Bconc + dc
    resultsM = linked_PDE_coefficients_1D(wrapper, params, angles, finedeg)
    YAm = resultsM[1]['m0e_vals']
    YBp = resultsM[2]['m0e_vals']
    params.target[0][1] = Aconc
    params.target[1][1] = Bconc



  if Aconc + dc > 1.0:
    # if the right neighbor has Aconc > 1.0, then we are at the right edge
    # and must use a leftward finite difference
    difftype = "left"
    diffdenom = dc
    YAp = results[1]['m0e_vals']
    YBm = results[2]['m0e_vals']
  else:
    # otherwise, we can get the values of the coefficients
    # at a larger concentration of species A
    params.target[0][1] = Aconc + dc
    params.target[1][1] = Bconc - dc
    resultsP = linked_PDE_coefficients_1D(wrapper, params, angles, finedeg)
    YAp = resultsP[1]['m0e_vals']
    YBm = resultsP[2]['m0e_vals']
    params.target[0][1] = Aconc
    params.target[1][1] = Bconc


  print "Aconc = %4f,  Bconc = %4f.  YAp = %4f, YAm = %4f, YBp = %4f, YBm = %4f." % (Aconc, Bconc, YAp[0], YAm[0], YBp[0], YBm[0])

  # now calculate coefficients at theta = 0.
  YAprime = (YAp[0] - YAm[0]) / (diffdenom)
  YBprime = (YBp[0] - YBm[0]) / (diffdenom)

  # extract values of A,C,Ap,Cp at theta=0
  retvals = dict()
  retvals["YA"] = YA0[0]
  retvals["YB"] = YB0[0]
  retvals["A"]  = - YAprime + YBprime
  retvals["C"]  = ( SXA[0] + SXB[0] )
  retvals["Ap"] = ( - cbulk[1]*YAprime - cbulk[0]*YBprime ) / film_height
  retvals["Cp"] = (   cbulk[1]*SXA[0]  - cbulk[0]*SXB[0]  ) / film_height
  retvals["Ceff"] = retvals["C"] - retvals["A"] * retvals["Cp"] / retvals["Ap"]
  return retvals
  




# --------------------------------------------------
#    plot_angle_dependence_summary()
# --------------------------------------------------


def plot_single_flat_angle_dependence_summary(fitted_values):

  # rename the argument
  rvals = fitted_values

  # plot the data and the fit
  theplot = plt.figure(figsize=(12,9))

  plt.subplot(231)
  plt.errorbar(rvals["angles"], rvals["m0e_dats"], yerr=rvals["m0e_errs"], fmt='gs', label='data')
  plt.plot(rvals["finedeg"], rvals["m0e_vals"], 'g-', linewidth=2, label="fit")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$M^{(0)}_{\mathsf{eros.}} \left( \theta \right)$')
  plt.legend(loc=3)
  plt.title('Sputter Yield')

  plt.subplot(232)
  plt.errorbar(rvals["angles"], rvals["m1e_dats"], yerr=rvals["m1e_errs"], fmt='rs', label='data')
  plt.plot(rvals["finedeg"], rvals["m1e_vals"], 'r-', linewidth=2, label="fit")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$M^{(1)}_{\mathsf{eros.}} \left( \theta \right)$')
  plt.legend(loc=3)
  plt.title('1st erosive moment')

  plt.subplot(233)
  plt.errorbar(rvals["angles"], rvals["m1r_dats"], yerr=rvals["m1r_errs"], fmt='bs', label='data')
  plt.plot(rvals["finedeg"], rvals["m1r_vals"], 'b-', linewidth=2, label="fit")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$M^{(1)}_{\mathsf{redist.}} \left( \theta \right)$')
  plt.legend(loc=2)
  plt.title('1st redist. moment')

  plt.subplot(234)
  plt.plot(rvals["finedeg"], rvals["sxe_coeffs"], 'r-', linewidth=2, label=r"$S_{X,\mathsf{eros.}}$")
  plt.plot(rvals["finedeg"], rvals["sxr_coeffs"], 'b-', linewidth=2, label=r"$S_{X,\mathsf{redist.}}$")
  plt.plot(rvals["finedeg"], rvals["sxc_coeffs"], 'g-', linewidth=2, label=r"$S_{X,\mathsf{curv.}}$")
  plt.plot(rvals["finedeg"], rvals["sx_coeffs"], 'k--', linewidth=2, label=r"$S_{X,\mathsf{total}}$")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{X} \left( \theta \right)$')
  plt.legend(loc=1)
  plt.title('SX components')

  plt.subplot(235)
  plt.plot(rvals["finedeg"], rvals["sye_coeffs"], 'r-', linewidth=2, label=r"$S_{Y,\mathsf{eros.}}$")
  plt.plot(rvals["finedeg"], rvals["syr_coeffs"], 'b-', linewidth=2, label=r"$S_{Y,\mathsf{redist.}}$")
  plt.plot(rvals["finedeg"], rvals["syc_coeffs"], 'g-', linewidth=2, label=r"$S_{Y,\mathsf{curv.}}$")
  plt.plot(rvals["finedeg"], rvals["sy_coeffs"], 'k--', linewidth=2, label=r"$S_{Y,\mathsf{total}}$")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{Y} \left( \theta \right)$')
  plt.legend(loc=1)
  plt.title('SY components')

  plt.subplot(236)
  plt.plot(rvals["finedeg"], rvals["sx_coeffs"], 'c-', linewidth=2, label="SX")
  plt.plot(rvals["finedeg"], rvals["sy_coeffs"], 'm-', linewidth=2, label="SY")
  plt.plot([0, 90], [0, 0], 'k--', linewidth=1)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{X,Y} \left( \theta \right)$')
  plt.legend(loc=3)
  plt.title('SX and SY')

  plt.tight_layout()
  return theplot




# --------------------------------------------------
#    plot_binary_flat_angle_dependence_summary()
# --------------------------------------------------


def plot_binary_flat_angle_dependence_summary(fitted_values):

  # rename the argument
  rvals1 = fitted_values[1]
  rvals2 = fitted_values[2]

  # plot the data and the fit
  theplot = plt.figure(figsize=(20,15))

  plt.subplot(331)
  plt.errorbar(rvals1["angles"],  rvals1["m0e_dats"],  yerr=rvals1["m0e_errs"], fmt='gs', label='data')
  plt.plot    (rvals1["finedeg"], rvals1["m0e_vals"], 'g-', linewidth=2, label="fit")
  plt.errorbar(rvals2["angles"],  rvals2["m0e_dats"],  yerr=rvals2["m0e_errs"], fmt='bs', label='data')
  plt.plot    (rvals2["finedeg"], rvals2["m0e_vals"], 'b-', linewidth=2, label="fit")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$M^{(0)}_{\mathsf{eros.}} \left( \theta \right)$')
  plt.legend(loc=2)
  plt.title('Sputter Yield')

  plt.subplot(332)
  plt.errorbar(rvals1["angles"],  rvals1["m1e_dats"], yerr=rvals1["m1e_errs"], fmt='rs', label='data')
  plt.plot    (rvals1["finedeg"], rvals1["m1e_vals"], 'r-', linewidth=2, label="fit")
  plt.errorbar(rvals2["angles"],  rvals2["m1e_dats"], yerr=rvals2["m1e_errs"], fmt='rs', label='data')
  plt.plot    (rvals2["finedeg"], rvals2["m1e_vals"], 'r-', linewidth=2, label="fit")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$M^{(1)}_{\mathsf{eros.}} \left( \theta \right)$')
  plt.legend(loc=2)
  plt.title('1st erosive moment')

  plt.subplot(333)
  plt.errorbar(rvals1["angles"],  rvals1["m1r_dats"], yerr=rvals1["m1r_errs"], fmt='bs', label='data')
  plt.plot    (rvals1["finedeg"], rvals1["m1r_vals"], 'b-', linewidth=2, label="fit")
  plt.errorbar(rvals2["angles"],  rvals2["m1r_dats"], yerr=rvals2["m1r_errs"], fmt='bs', label='data')
  plt.plot    (rvals2["finedeg"], rvals2["m1r_vals"], 'b-', linewidth=2, label="fit")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$M^{(1)}_{\mathsf{redist.}} \left( \theta \right)$')
  plt.legend(loc=2)
  plt.title('1st redist. moment')

  plt.subplot(334)
  plt.plot(rvals1["finedeg"], rvals1["sxe_coeffs"], 'm-', linewidth=2, label=r"$S_{X,\mathsf{eros.}}$")
  plt.plot(rvals1["finedeg"], rvals1["sxr_coeffs"], 'c-', linewidth=2, label=r"$S_{X,\mathsf{redist.}}$")
  plt.plot(rvals1["finedeg"], rvals1["sx_coeffs"],  'k--', linewidth=2, label=r"$S_{X,\mathsf{total}}$")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{X} \left( \theta \right)$')
  plt.legend(loc=2)
  plt.title('component 1 SX')

  plt.subplot(335)
  plt.plot(rvals2["finedeg"], rvals2["sxe_coeffs"], 'm-', linewidth=2, label=r"$S_{Y,\mathsf{eros.}}$")
  plt.plot(rvals2["finedeg"], rvals2["sxr_coeffs"], 'c-', linewidth=2, label=r"$S_{Y,\mathsf{redist.}}$")
  plt.plot(rvals2["finedeg"], rvals2["sx_coeffs"],  'k--', linewidth=2, label=r"$S_{Y,\mathsf{total}}$")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{Y} \left( \theta \right)$')
  plt.legend(loc=2)
  plt.title('component 2 SX')

  plt.subplot(336)
  plt.plot(rvals1["finedeg"], rvals1["sxe_coeffs"] + rvals2["sxe_coeffs"], 'm-', linewidth=2, label="SX")
  plt.plot(rvals1["finedeg"], rvals1["sxr_coeffs"] + rvals2["sxr_coeffs"], 'c-', linewidth=2, label="SY")
  plt.plot(rvals1["finedeg"], rvals1["sx_coeffs"]  + rvals2["sx_coeffs"], 'k--', linewidth=2, label="SY")
  plt.plot([0, 90], [0, 0], 'k--', linewidth=1)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{X,Y} \left( \theta \right)$')
  plt.legend(loc=2)
  plt.title('total SX')

  plt.subplot(337)
  plt.plot(rvals1["finedeg"], rvals1["sye_coeffs"], 'm-', linewidth=2, label=r"$S_{X,\mathsf{eros.}}$")
  plt.plot(rvals1["finedeg"], rvals1["syr_coeffs"], 'c-', linewidth=2, label=r"$S_{X,\mathsf{redist.}}$")
  plt.plot(rvals1["finedeg"], rvals1["sy_coeffs"],  'k--', linewidth=2, label=r"$S_{X,\mathsf{total}}$")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{X} \left( \theta \right)$')
  plt.legend(loc=2)
  plt.title('component 1 SY')

  plt.subplot(338)
  plt.plot(rvals2["finedeg"], rvals2["sye_coeffs"], 'm-', linewidth=2, label=r"$S_{Y,\mathsf{eros.}}$")
  plt.plot(rvals2["finedeg"], rvals2["syr_coeffs"], 'c-', linewidth=2, label=r"$S_{Y,\mathsf{redist.}}$")
  plt.plot(rvals2["finedeg"], rvals2["sy_coeffs"],  'k--', linewidth=2, label=r"$S_{Y,\mathsf{total}}$")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{Y} \left( \theta \right)$')
  plt.legend(loc=2)
  plt.title('component 2 SY')

  plt.subplot(339)
  plt.plot(rvals1["finedeg"], rvals1["sye_coeffs"] + rvals2["sye_coeffs"], 'm-', linewidth=2, label="SX")
  plt.plot(rvals1["finedeg"], rvals1["syr_coeffs"] + rvals2["syr_coeffs"], 'c-', linewidth=2, label="SY")
  plt.plot(rvals1["finedeg"], rvals1["sy_coeffs"]  + rvals2["sy_coeffs"], 'k--', linewidth=2, label="SY")
  plt.plot([0, 90], [0, 0], 'k--', linewidth=1)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{X,Y} \left( \theta \right)$')
  plt.legend(loc=2)
  plt.title('total SY')

  return theplot














# --------------------------------------------------
#    linked_PDE_coefficients_2D()
# --------------------------------------------------


def linked_PDE_coefficients_2D(wrapper, params, angles, finedeg, fitmethod=None, fitguess=None, path=None, dk=None):

  '''
  This function calculates *curved-target* estimates of the values of the coefficients
  through second order in the linearized equation of motion.  It provides the total
  values of these coefficients, and also the contribution from each component (if
  the target is a multi-component target).  
  
  TODO:  FIGURE OUT EXACT METHOD USED, AND HOW THIS DIFFERS FROM LATER FUNCTION
  TODO:  Provide the *total* value within the zeroth entry field.
  TODO:  Also calculate the coefficient of the first-order hx term.
  '''

  oldpath = None
  if path != None:
    oldpath = os.getcwd()
    os.chdir(path)

  # run simulations if needed
  params.k11 = 0 ; params.k22 = 0
  for aa in angles:
    params.angle = aa
    wrapper.go(params)
  params.angle = None

  if dk:

    params.k11 = dk; params.k22 = 0
    for aa in angles:
      params.angle = aa
      wrapper.go(params)

    params.k11 = -dk; params.k22 = 0
    for aa in angles:
      params.angle = aa
      wrapper.go(params)

    params.k11 = 0; params.k22 = dk
    for aa in angles:
      params.angle = aa
      wrapper.go(params)

    params.k11 = 0; params.k22 = -dk
    for aa in angles:
      params.angle = aa
      wrapper.go(params)

    params.k11 = 0; params.k22 = 0



  # determine number of species, and number of impacts
  species = len(params.target)
  impacts = params.impacts

  # allocate storage
  results_list = []
  totals = dict()
  results_list.append(totals)
  
  for kk in range(species):
  
    # storage for returned values
    rvals = dict()
    rvals["angles"] = angles
    rvals["finedeg"] = finedeg

    # functions we will use for fitting
    fitfuncD0 = pyam.poly_yamamura
    fitfuncD1 = pyam.poly_yamamura_derivative

    # setup the FitMethod
    evenopts = {"symmetry":"even", "coeff_list":['A', 'B', 'C']}
    oddopts  = {"symmetry":"odd" , "coeff_list":['A', 'B', 'C']}
    fitnames = ["A", "B", "C"]
    fitguess = [1.0, 0.0, 0.0]
    evenfit = fittools.FitMethod( pyam.poly_yamamura, fitnames, evenopts, pyam.poly_yamamura_cerror )
    oddfit  = fittools.FitMethod( pyam.poly_yamamura, fitnames, oddopts, pyam.poly_yamamura_cerror )

    # storage for nodelist
    nodelist = []

    # M0e data
    m0e_avg = io.array_range( './', params, 'angle', angles, 'm0e_avg' )
    m0e_std = io.array_range( './', params, 'angle', angles, 'm0e_std' ) 
    m0e_vals = np.array( [ a[kk] for a in m0e_avg ] )
    m0e_errs = np.array( [ b[kk]/np.sqrt(params.impacts) * 1.97 for b in m0e_std ] )
    m0e_node = evenfit.CreateFitNode( list(fitguess), angles, m0e_vals, m0e_errs )
    rvals["m0e_avg"] = m0e_vals
    rvals["m0e_err"] = m0e_errs
    nodelist.append(m0e_node)

    # M1e data
    m1e_avg = io.array_range( './', params, 'angle', angles, 'm1e_avg' )
    m1e_std = io.array_range( './', params, 'angle', angles, 'm1e_std' ) 
    m1e_vals = np.array( [ a[kk][0] for a in m1e_avg ] )
    m1e_errs = np.array( [ b[kk][0]/np.sqrt(params.impacts) * 1.97 for b in m1e_std ] )
    m1e_node = oddfit.CreateFitNode( list(fitguess), angles, m1e_vals, m1e_errs )
    rvals["m1e_avg"] = m1e_vals
    rvals["m1e_err"] = m1e_errs
    nodelist.append(m1e_node)

    # M1r data
    m1r_avg = io.array_range( './', params, 'angle', angles, 'm1r_avg' )
    m1r_std = io.array_range( './', params, 'angle', angles, 'm1r_std' ) 
    m1r_vals = np.array( [ a[kk][0] for a in m1r_avg ] )
    m1r_errs = np.array( [ b[kk][0]/np.sqrt(params.impacts) * 1.97 for b in m1r_std ] )
    m1r_node = oddfit.CreateFitNode( list(fitguess), angles, m1r_vals, m1r_errs )
    rvals["m1r_avg"] = m1r_vals
    rvals["m1r_err"] = m1r_errs
    nodelist.append(m1r_node)

    
    if (dk):

      rvals["dk"] = dk

      # K11+
      params.k11 = dk; params.k22 = 0
      m0ek1p_avg = io.array_range( './', params, 'angle', angles, 'm0e_avg' )
      m0ek1p_std = io.array_range( './', params, 'angle', angles, 'm0e_std' ) 
      m0ek1p_vals = np.array( [ a[kk] for a in m0ek1p_avg ] )
      m0ek1p_errs = np.array( [ b[kk]/np.sqrt(params.impacts) * 1.97 for b in m0ek1p_std ] )
      m0ek1p_node = evenfit.CreateFitNode( list(fitguess), angles, m0ek1p_vals, m0ek1p_errs )
      rvals["m0ek1p_avg"] = m0ek1p_vals
      rvals["m0ek1p_err"] = m0ek1p_errs
      nodelist.append(m0ek1p_node)

      # K11-
      params.k11 = -dk; params.k22 = 0
      m0ek1m_avg = io.array_range( './', params, 'angle', angles, 'm0e_avg' )
      m0ek1m_std = io.array_range( './', params, 'angle', angles, 'm0e_std' ) 
      m0ek1m_vals = np.array( [ a[kk] for a in m0ek1m_avg ] )
      m0ek1m_errs = np.array( [ b[kk]/np.sqrt(params.impacts) * 1.97 for b in m0ek1m_std ] )
      m0ek1m_node = evenfit.CreateFitNode( list(fitguess), angles, m0ek1m_vals, m0ek1m_errs )
      rvals["m0ek1m_avg"] = m0ek1m_vals
      rvals["m0ek1m_err"] = m0ek1m_errs
      nodelist.append(m0ek1m_node)

      # K22+
      params.k11 = 0; params.k22 = dk
      m0ek2p_avg = io.array_range( './', params, 'angle', angles, 'm0e_avg' )
      m0ek2p_std = io.array_range( './', params, 'angle', angles, 'm0e_std' ) 
      m0ek2p_vals = np.array( [ a[kk] for a in m0ek2p_avg ] )
      m0ek2p_errs = np.array( [ b[kk]/np.sqrt(params.impacts) * 1.97 for b in m0ek2p_std ] )
      m0ek2p_node = evenfit.CreateFitNode( list(fitguess), angles, m0ek2p_vals, m0ek2p_errs )
      rvals["m0ek2p_avg"] = m0ek2p_vals
      rvals["m0ek2p_err"] = m0ek2p_errs
      nodelist.append(m0ek2p_node)

      # K22-
      params.k11 = 0; params.k22 = -dk
      m0ek2m_avg = io.array_range( './', params, 'angle', angles, 'm0e_avg' )
      m0ek2m_std = io.array_range( './', params, 'angle', angles, 'm0e_std' ) 
      m0ek2m_vals = np.array( [ a[kk] for a in m0ek2m_avg ] )
      m0ek2m_errs = np.array( [ b[kk]/np.sqrt(params.impacts) * 1.97 for b in m0ek2m_std ] )
      m0ek2m_node = evenfit.CreateFitNode( list(fitguess), angles, m0ek2m_vals, m0ek2m_errs )
      rvals["m0ek2m_avg"] = m0ek2m_vals
      rvals["m0ek2m_err"] = m0ek2m_errs
      nodelist.append(m0ek2m_node)

      params.k11 = 0; params.k22 = 0



    # SuperNode and Fitting
    superfit = fittools.FitSuperNode("p", [0.5], nodelist)
    fits, err = superfit.fit_data()
    m0e_fit = fits[0]
    m1e_fit = fits[1]
    m1r_fit = fits[2]

    #finedeg = np.linspace(0,90,91)
    finerad = finedeg * np.pi / 180.

    # record fitted values of the functions, themselves
    rvals["m0e_vals"] = fitfuncD0(m0e_fit, finedeg, evenopts)

    rvals["m1e_vals"] = fitfuncD0(m1e_fit, finedeg,  oddopts)
    rvals["m1r_vals"] = fitfuncD0(m1r_fit, finedeg,  oddopts)
    rvals["m1_vals"]  = rvals["m1e_vals"] + rvals["m1r_vals"]

    # record fitted values of the derivatives of M1
    rvals["m1ep_vals"] = fitfuncD1( m1e_fit, finedeg, oddopts )
    rvals["m1rp_vals"] = fitfuncD1( m1r_fit, finedeg, oddopts )
    rvals["m1p_vals"]  = rvals["m1ep_vals"] + rvals["m1rp_vals"]

    # construct SX
    rvals["sxe_coeffs"] = np.cos(finerad) * rvals["m1ep_vals"] - np.sin(finerad) * rvals["m1e_vals"]
    rvals["sxr_coeffs"] = np.cos(finerad) * rvals["m1rp_vals"] - np.sin(finerad) * rvals["m1r_vals"]
    rvals["sxc_coeffs_approx"] = - np.cos(finerad) * rvals["m1ep_vals"] / 2.0

    rvals["sx_coeffs"] = rvals["sxe_coeffs"] + rvals["sxr_coeffs"]

    # construct SY
    rvals["sye_coeffs"] = np.cos(finerad)**2 / np.sin(finerad) * rvals["m1e_vals"]
    rvals["syr_coeffs"] = np.cos(finerad)**2 / np.sin(finerad) * rvals["m1r_vals"]
    rvals["syc_coeffs_approx"] = - np.cos(finerad)**2 / np.sin(finerad) * rvals["m1e_vals"] / 2.0
    rvals["sy_coeffs"] = rvals["sye_coeffs"] + rvals["syr_coeffs"]

    if dk:

      k1p_fit = fits[3]
      k1m_fit = fits[4]
      k2p_fit = fits[5]
      k2m_fit = fits[6]

      print m0e_fit
      print m1e_fit
      print m1r_fit
      print k1p_fit
      print k1m_fit
      print k2p_fit
      print k2m_fit

      # record fitted values of the functions, themselves
      rvals["m0ek1p_vals"] = fitfuncD0(k1p_fit, finedeg, evenopts)
      rvals["m0ek1m_vals"] = fitfuncD0(k1m_fit, finedeg, evenopts)
      rvals["m0ek2p_vals"] = fitfuncD0(k2p_fit, finedeg, evenopts)
      rvals["m0ek2m_vals"] = fitfuncD0(k2m_fit, finedeg, evenopts)

      # construct corrections to SX, SY
      rvals["sxc_coeffs"] = - np.cos(finerad) * ( rvals["m0ek1p_vals"] - rvals["m0ek1m_vals"] ) / (2.0*dk)
      rvals["syc_coeffs"] = - np.cos(finerad) * ( rvals["m0ek2p_vals"] - rvals["m0ek2m_vals"] ) / (2.0*dk)
      rvals["sx_coeffs"] += rvals["sxc_coeffs"]
      rvals["sy_coeffs"] += rvals["syc_coeffs"]

    results_list.append(rvals)


  if oldpath != None:
    os.chdir(oldpath)

  # construct the sum dictionary
  return results_list
  








def plot_single_curved_angle_dependence_summary(fitted_values):

  # rename the argument
  rvals = fitted_values
  dk = rvals["dk"]


  # plot the data and the fit
  theplot = plt.figure(figsize=(16,9))

  plt.subplot(231)
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

  plt.subplot(232)
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


  plt.subplot(233)
  plt.errorbar(rvals["angles"], rvals["m1e_avg"], yerr=rvals["m1e_err"], fmt='rs', label='eros.')
  plt.errorbar(rvals["angles"], rvals["m1r_avg"], yerr=rvals["m1r_err"], fmt='bs', label='redist.')
  plt.plot(rvals["finedeg"], rvals["m1e_vals"], 'r-', linewidth=2)
  plt.plot(rvals["finedeg"], rvals["m1r_vals"], 'b-', linewidth=2)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$M^{(1)} \left( \theta \right)$')
  plt.legend(loc=3)
  plt.title('First Moments')


  plt.subplot(234)
  plt.plot(rvals["finedeg"], rvals["sxe_coeffs"], 'r-', linewidth=2, label=r"$S_{X,\mathsf{eros.}}$")
  plt.plot(rvals["finedeg"], rvals["sxr_coeffs"], 'b-', linewidth=2, label=r"$S_{X,\mathsf{redist.}}$")
  plt.plot(rvals["finedeg"], rvals["sxc_coeffs"], 'g-', linewidth=2, label=r"$S_{X,\mathsf{curv.}}$")
  plt.plot(rvals["finedeg"], rvals["sx_coeffs"], 'k--', linewidth=2, label=r"$S_{X,\mathsf{total}}$")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{X} \left( \theta \right)$')
  plt.legend(loc=2)
  plt.title('SX components')

  plt.subplot(235)
  plt.plot(rvals["finedeg"], rvals["sye_coeffs"], 'r-', linewidth=2, label=r"$S_{Y,\mathsf{eros.}}$")
  plt.plot(rvals["finedeg"], rvals["syr_coeffs"], 'b-', linewidth=2, label=r"$S_{Y,\mathsf{redist.}}$")
  plt.plot(rvals["finedeg"], rvals["syc_coeffs"], 'g-', linewidth=2, label=r"$S_{X,\mathsf{curv.}}$")
  plt.plot(rvals["finedeg"], rvals["sy_coeffs"], 'k--', linewidth=2, label=r"$S_{Y,\mathsf{total}}$")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{Y} \left( \theta \right)$')
  plt.legend(loc=2)
  plt.title('SY components')

  plt.subplot(236)
  plt.plot(rvals["finedeg"], rvals["sx_coeffs"], 'c-', linewidth=2, label="SX")
  plt.plot(rvals["finedeg"], rvals["sy_coeffs"], 'm-', linewidth=2, label="SY")
  plt.plot([0, 90], [0, 0], 'k--', linewidth=1)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{X,Y} \left( \theta \right)$')
  plt.legend(loc=3)
  plt.title('SX and SY')

  plt.tight_layout()
  return theplot





























# --------------------------------------------------
#    calculate_PDE_coefficients()
# --------------------------------------------------

def calculate_PDE_coefficients(wrapper, params, angles, fitmethod, guess, finedeg=None, path='./', curvature_effects_dk=None):

  '''
  This function performs most of the angle-dependent fitting and differenting needed to 
  estimate values of coefficients in PDE terms.  The user must supply 

  - a BCA_wrapper for running the simulations.
  - a params object, in which everything is specified except the angles
  - a list of angles to use in fitting
  - a FitMethod object
  - a guess at paramters

  OPTIONAL
  - a fine mesh on which to return fits and their derivatives (defaults to linspace(0,90,91))
  - path to the data files (defaults to './')
  '''

  # set up fine degree mesh and fine radian mesh
  if finedeg==None:
    finedeg = np.linspace(0,90,91) ;
  finerad = finedeg * np.pi / 180.0

  # identify number of species
  species = len(params.target)

  # if the files are not already present, then run the simulations
  for aa in angles:
    params.angle = aa;  wrapper.go(params)
    if (curvature_effects_dk):

      dk = curvature_effects_dk
      params.k22 = 0.0
      params.k11 =  dk;  wrapper.go(params);
      params.k11 = -dk;  wrapper.go(params);
      params.k11 = 0.0
      params.k22 =  dk;  wrapper.go(params);
      params.k22 = -dk;  wrapper.go(params);
      params.k22 = 0.0    


  # array for number of impacts
  impacts = np.array( io.array_range(path, params, "angle", angles, "impacts") )


  # now, a really big loop, to generate fits on the moments and coefficients for all species
  species_list = []
  for kk in range(species):

    #m1e_dats = np.array( m1e_avg[:,kk,1] )
    #m1e_errs = np.array( m1e_err[:,kk,1] )
    #m1r_dats = np.array( m1r_avg[:,kk,1] )
    #m1r_errs = np.array( m1r_err[:,kk,1] )


    rvals = dict()

    # extract the data for m1e (values, uncertainty), and try to fit it with the poly_yamamura() function
    fitmethod.params["symmetry"] = "even"
    m0e_avg = np.array( io.array_range(path, params, "angle", angles, "m0e_avg") )
    m0e_std = np.array( io.array_range(path, params, "angle", angles, "m0e_std") )
    m0e_err = np.array( [a / np.sqrt(b) * 1.97 for a,b in zip(m0e_std,impacts)] )
    m0e_fit = fitmethod.fit_data( guess, angles, np.array( m0e_avg[:,kk] ),   stderrs=np.array( m0e_err[:,kk] ) )
    rvals["m0e_avg"]  = m0e_avg[:,kk]
    rvals["m0e_err"]  = m0e_err[:,kk]
    rvals["m0e_vals"] = fitmethod.eval_function( m0e_fit, finedeg )


    # extract the data for m1e (values, uncertainty), and try to fit it with the poly_yamamura() function
    fitmethod.params["symmetry"] = "odd"
    m1e_avg = np.array( io.array_range(path, params, "angle", angles, "m1e_avg") )
    m1e_std = np.array( io.array_range(path, params, "angle", angles, "m1e_std") )
    m1e_err = np.array( [a / np.sqrt(b) * 1.97 for a,b in zip(m1e_std,impacts)] )
    m1e_fit = fitmethod.fit_data( guess, angles, np.array( m1e_avg[:,kk,1] ), stderrs=np.array( m1e_err[:,kk,1] ) )
    rvals["m1e_avg"]  = m1e_avg[:,kk,1]
    rvals["m1e_err"]  = m1e_err[:,kk,1]
    rvals["m1e_vals"] = fitmethod.eval_function( m1e_fit, finedeg )

    # extract the data for m1r (values, uncertainty), and try to fit it with the poly_yamamura() function
    m1r_avg = np.array( io.array_range(path, params, "angle", angles, "m1r_avg") )
    m1r_std = np.array( io.array_range(path, params, "angle", angles, "m1r_std") )
    m1r_err = np.array( [a / np.sqrt(b) * 1.97 for a,b in zip(m1r_std,impacts)] )
    m1r_fit = fitmethod.fit_data( guess, angles, np.array( m1r_avg[:,kk,1] ), stderrs=np.array( m1r_err[:,kk,1] ) )
    rvals["m1r_avg"]  = m1r_avg[:,kk,1]
    rvals["m1r_err"]  = m1r_err[:,kk,1]
    rvals["m1r_vals"] = fitmethod.eval_function( m1r_fit, finedeg )

    rvals["m1_vals"]  = rvals["m1e_vals"] + rvals["m1r_vals"]

    rvals["m1ep_vals"] = fitmethod.eval_derivative( m1e_fit, finedeg )
    rvals["m1rp_vals"] = fitmethod.eval_derivative( m1r_fit, finedeg )
    rvals["m1p_vals"]  = rvals["m1ep_vals"] + rvals["m1rp_vals"]

    rvals["sxe_coeffs"] = cos(finerad) * rvals["m1ep_vals"] - sin(finerad) * rvals["m1e_vals"]
    rvals["sxr_coeffs"] = cos(finerad) * rvals["m1rp_vals"] - sin(finerad) * rvals["m1r_vals"]
    rvals["sx_coeffs"] = rvals["sxe_coeffs"] + rvals["sxr_coeffs"]

    rvals["sye_coeffs"] = cos(finerad)**2 / sin(finerad) * rvals["m1e_vals"]
    rvals["syr_coeffs"] = cos(finerad)**2 / sin(finerad) * rvals["m1r_vals"]
    rvals["sy_coeffs"] = rvals["sye_coeffs"] + rvals["syr_coeffs"]



    if (curvature_effects_dk):
      dk = curvature_effects_dk
      fitmethod.params["symmetry"] = "even"

      params.k22 = 0.0
      params.k11 =  dk;
      m0e_avg  = np.array( io.array_range(path, params, "angle", angles, "m0e_avg") )
      m0e_std  = np.array( io.array_range(path, params, "angle", angles, "m0e_std") )
      m0e_err  = np.array( [a / np.sqrt(b) * 1.97 for a,b in zip(m0e_std,impacts)] )
      m0e_fit  = fitmethod.fit_data( guess, angles, np.array( m0e_avg[:,kk] ),   stderrs=np.array( m0e_err[:,kk] ) )
      k11p     = fitmethod.eval_function( m0e_fit, finedeg )
      rvals["m0ek1p_avg"] = m0e_avg[:,kk] ;
      rvals["m0ek1p_err"] = m0e_err[:,kk] ;
      rvals["m0ek1p_vals"] = k11p
        
      params.k11 = -dk
      m0e_avg  = np.array( io.array_range(path, params, "angle", angles, "m0e_avg") )
      m0e_std  = np.array( io.array_range(path, params, "angle", angles, "m0e_std") )
      m0e_err  = np.array( [a / np.sqrt(b) * 1.97 for a,b in zip(m0e_std,impacts)] )
      m0e_fit  = fitmethod.fit_data( guess, angles, np.array( m0e_avg[:,kk] ),   stderrs=np.array( m0e_err[:,kk] ) )
      k11m     = fitmethod.eval_function( m0e_fit, finedeg )
      rvals["m0ek1m_avg"] = m0e_avg[:,kk]
      rvals["m0ek1m_err"] = m0e_err[:,kk]
      rvals["m0ek1m_vals"] = k11m
      params.k11 = 0.0

      params.k11 = 0.0
      params.k22 =  dk
      m0e_avg  = np.array( io.array_range(path, params, "angle", angles, "m0e_avg") )
      m0e_std  = np.array( io.array_range(path, params, "angle", angles, "m0e_std") )
      m0e_err  = np.array( [a / np.sqrt(b) * 1.97 for a,b in zip(m0e_std,impacts)] )
      m0e_fit  = fitmethod.fit_data( guess, angles, np.array( m0e_avg[:,kk] ),   stderrs=np.array( m0e_err[:,kk] ) )
      k22p     = fitmethod.eval_function( m0e_fit, finedeg )
      rvals["m0ek2p_avg"] = m0e_avg[:,kk]
      rvals["m0ek2p_err"] = m0e_err[:,kk]
      rvals["m0ek2p_vals"] = k22p

      params.k22 = -dk
      m0e_avg  = np.array( io.array_range(path, params, "angle", angles, "m0e_avg") )
      m0e_std  = np.array( io.array_range(path, params, "angle", angles, "m0e_std") )
      m0e_err  = np.array( [a / np.sqrt(b) * 1.97 for a,b in zip(m0e_std,impacts)] )
      m0e_fit  = fitmethod.fit_data( guess, angles, np.array( m0e_avg[:,kk] ),   stderrs=np.array( m0e_err[:,kk] ) )
      k22m     = fitmethod.eval_function( m0e_fit, finedeg )
      rvals["m0ek2m_avg"] = m0e_avg[:,kk]
      rvals["m0ek2m_err"] = m0e_err[:,kk]
      rvals["m0ek2m_vals"] = k22m
      params.k22 = 0.0
      fitmethod.params["symmetry"] = "odd"


      rvals["sxe_curvcorr"] = (k11p - k11m) / (2*dk) * np.cos(finerad)
      rvals["sye_curvcorr"] = (k22p - k22m) / (2*dk) * np.cos(finerad)
      #rvals["sxe_curvcorr"] = (k11p - rvals["m0e_vals"]) / (dk) * np.cos(finerad)
      #rvals["sye_curvcorr"] = (k22p - rvals["m0e_vals"]) / (dk) * np.cos(finerad)
      rvals["sx_coeffs"] += rvals["sxe_curvcorr"]
      rvals["sy_coeffs"] += rvals["sye_curvcorr"]

    species_list.append(rvals)

  return species_list


