import os
import sys
import copy
import uncertainties.umath  # math  ### ilinov
import shelve
import numpy as np
import uncertainties.unumpy as unp  ### ilinov
import matplotlib.pyplot as plt

import IO as io		# allows data extraction from files.
import lmfit as lmf
import fits.NestedLmfit as fitter
import fits.poly_yamamura as pyam
import fits.polynomial as poly






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
  m0e_avg = -np.array(io.array_range( './', params, 'angle', angles, 'evdM0_avg' ))
  m0e_std =  np.array(io.array_range( './', params, 'angle', angles, 'evdM0_std' ))
  m1e_avg = -np.array(io.array_range( './', params, 'angle', angles, 'evdM1_avg' ))
  m1e_std =  np.array(io.array_range( './', params, 'angle', angles, 'evdM1_std' ))
  m1r_avg =  np.array(io.array_range( './', params, 'angle', angles, 'rddM1_avg' ))
  m1r_std =  np.array(io.array_range( './', params, 'angle', angles, 'rddM1_std' ))

  # create even fit parameter
  fp_even = lmf.Parameters()
  fp_even.add("p" , value=0.25, min=0.0)
  fp_even.add("c0", value=0.00)
  fp_even.add("c2", value=0.00)
  fp_even.add("c4", value=0.00)
  #fp_even.add("c6", value=0.00)


  # create odd fit parameter
  fp_odd  = lmf.Parameters()
  fp_odd.add("p" , value=0.25, min=0.0)
  fp_odd.add("c1", value=0.00)
  fp_odd.add("c3", value=0.00)
  fp_odd.add("c5", value=0.00)
  #fp_odd.add("c7", value=0.00)

  # greate a global fit parameter
  fp_global  = None
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
    m0e_node = fitter.FitNode( copy.deepcopy(fp_even), pyam.poly_yamamura, angles, m0e_vals, m0e_errs )
    nodelist.append(m0e_node)

    # M1e data
    m1e_vals = np.array( [ a[kk][0] for a in m1e_avg ] )
    m1e_errs = np.array( [ b[kk][0]/np.sqrt(params.impacts) * 1.97 for b in m1e_std ] )
    rvals["m1e_dats"] = m1e_vals
    rvals["m1e_errs"] = m1e_errs
    m1e_node = fitter.FitNode( copy.deepcopy(fp_odd), pyam.poly_yamamura, angles, m1e_vals, m1e_errs )
    nodelist.append(m1e_node)

    # M1r data
    m1r_vals = np.array( [ a[kk][0] for a in m1r_avg ] )
    m1r_errs = np.array( [ b[kk][0]/np.sqrt(params.impacts) * 1.97 for b in m1r_std ] )
    rvals["m1r_dats"] = m1r_vals
    rvals["m1r_errs"] = m1r_errs
    m1r_node = fitter.FitNode( copy.deepcopy(fp_odd), pyam.poly_yamamura, angles, m1r_vals, m1r_errs )
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







def fit_moments_2D(params, angles, curvatures, fitmethod=None, fitguess=None, path='./'):

  '''
  '''

  finedeg = np.linspace(0,90,91)
  curvatures = np.array(curvatures)

  # see if fit values already exist
  fname = params.fname(['target', 'beam', 'energy'])
  if os.path.isfile('%s.momfits' % (fname)):
    f = shelve.open('%s.momfits' % (fname))
    results_list = list(f['results_list'])
    f.close()
    return results_list

  # determine number of species, configure storage for pure or mutli-species results
  species = len(params.target)
  extra_species = 0
  if species > 1:  extra_species = 1
  total_species = species + extra_species


  # ----- FLAT DATA -----
  # load moment data
  m0e_avg = -np.array( io.array_range( './', params, 'angle', angles, 'evdM0_avg' ))
  m0e_std =  np.array( io.array_range( './', params, 'angle', angles, 'evdM0_std' ))
  m1e_avg = -np.array( io.array_range( './', params, 'angle', angles, 'evdM1_avg' ))
  m1e_std =  np.array( io.array_range( './', params, 'angle', angles, 'evdM1_std' ))
  m1r_avg =  np.array(io.array_range( './', params, 'angle', angles, 'rddM1_avg' ))
  m1r_std =  np.array(io.array_range( './', params, 'angle', angles, 'rddM1_std' ))

  # create even fit parameter
  fp_even = lmf.Parameters()
  fp_even.add("p" , value=0.25, min=0.0)
  fp_even.add("c0", value=0.00)
  fp_even.add("c2", value=0.00)
  fp_even.add("c4", value=0.00)
  fp_even.add("c6", value=0.00)

  # create odd fit parameter
  fp_odd  = lmf.Parameters()
  fp_odd.add("p" , value=0.25, min=0.0)
  fp_odd.add("c1", value=0.00)
  fp_odd.add("c3", value=0.00)
  fp_odd.add("c5", value=0.00)
  fp_odd.add("c7", value=0.00)

  # greate a global fit parameter
  fp_global  = None
  fp_global  = lmf.Parameters()
  fp_global.add("p", value=0.25, min=0.0)
  

  # allocate storage
  results_list = []



  for sid in xrange(total_species):

    # storage for returned values
    rvals = dict()
    rvals["angles"] = angles
    rvals["finedeg"] = finedeg

    # storage
    dk11_vals = []
    dk11_errs = []
    dk22_vals = []
    dk22_errs = []

    params.angle = None
    for aa in angles:

      params.angle = aa

      # ----- K11 ------
      # load moment data
      params.k22 = 0.0
      m0ek11_avg = np.array( io.array_range( './', params, 'k11', curvatures, 'evdM0_avg' ))
      m0ek11_std = np.array( io.array_range( './', params, 'k11', curvatures, 'evdM0_std' ))
      m0ek11_vals = np.array( [ a[sid] for a in m0ek11_avg ] )
      m0ek11_errs = np.array( [ b[sid]/np.sqrt(params.impacts) * 1.97 for b in m0ek11_std ] )

      c0_guess = np.mean(m0ek11_avg[sid])
      c1_guess = (m0ek11_avg[sid][-1]-m0ek11_avg[sid][0]) / (curvatures[-1]-curvatures[0])

      # create fit parameter
      fp1 = lmf.Parameters()
      fp1.add("c0", value=c0_guess)
      fp1.add("c1", value=c1_guess)
      fp1.add("c2", value=0.00)

      # fit data to parabola
      node1 = fitter.FitNode( fp1, poly.polynomial, curvatures, m0ek11_vals, m0ek11_errs )
      result1 = node1.fit()
      fitparams1 = node1.result.params

      # save the fitted values
      dk11_vals.append(fitparams1['c1'].value)
      dk11_errs.append(fitparams1['c1'].stderr)


      # ----- K22 ------
      # load moment data
      params.k11 = 0.0
      m0ek22_avg = np.array( io.array_range( './', params, 'k22', curvatures, 'evdM0_avg' ))
      m0ek22_std = np.array( io.array_range( './', params, 'k22', curvatures, 'evdM0_std' ))
      m0ek22_vals = np.array( [ a[sid] for a in m0ek22_avg ] )
      m0ek22_errs = np.array( [ b[sid]/np.sqrt(params.impacts) * 1.97 for b in m0ek22_std ] )

      c0_guess = np.mean(m0ek22_avg[sid])
      c1_guess = (m0ek22_avg[sid][-1]-m0ek22_avg[sid][0]) / (curvatures[-1]-curvatures[0])

      # create fit parameter
      fp2 = lmf.Parameters()
      fp2.add("c0", value=c0_guess)
      fp2.add("c1", value=c1_guess)
      fp2.add("c2", value=0.00)

      # fit data to parabola
      node2 = fitter.FitNode( fp2, poly.polynomial, curvatures, m0ek22_vals, m0ek22_errs )
      result2 = node2.fit()
      fitparams2 = node2.result.params

      # save the fitted values
      dk22_vals.append(fitparams2['c1'].value)
      dk22_errs.append(fitparams2['c1'].stderr)

      # diagnostic plot
      if (sid == 0):
        finek = np.linspace(curvatures[0], curvatures[-1], 51)

        # plot the data and the fit
        theplot = plt.figure(figsize=(10,3))
        plt.rcParams.update({'axes.titlesize': 12})
        plt.rcParams.update({'axes.labelsize': 9})
        plt.rcParams.update({'xtick.labelsize': 8})
        plt.rcParams.update({'ytick.labelsize': 8})
        textxloc = 0.475
        textyloc = 0.875

        # plot dependence of yield on curvature in X-direction
        plt.subplot(121)
        plt.errorbar(curvatures, m0ek11_vals, yerr=m0ek11_errs, fmt='bs', label="sims")
        plt.plot(finek, poly.polynomial( fitparams1, finek), 'r', linewidth=2, label="fit")
        plt.xlim(finek[0]*1.2, finek[-1]*1.2)
        ylim1 = plt.ylim()
        plt.xlabel(r'K11 curvature [A$^{-1}$]')
        plt.ylabel(r'$M^{\left(0\right)}_{\mathrm{eros.}} \left( K_{11} \right)$')
        plt.legend(loc='best',prop={'size':9})
        plt.text(textxloc, textyloc, "(b)", fontweight="bold", transform=plt.gca().transAxes)

        # plot dependence of yield on curvature in Y-direction
        plt.subplot(122)
        plt.errorbar(curvatures, m0ek22_vals, yerr=m0ek22_errs, fmt='bs', label="sims")
        plt.plot(finek, poly.polynomial( fitparams2, finek), 'r', linewidth=2, label="fit")
        plt.xlim(finek[0]*1.2, finek[-1]*1.2)
        ylim2 = plt.ylim()
        plt.xlabel(r'K22 curvature [A$^{-1}$]')
        plt.ylabel(r'$M^{\left(0\right)}_{\mathrm{eros.}} \left( K_{22} \right)$')
        plt.legend(loc='best',prop={'size':9})
        plt.text(textxloc, textyloc, "(c)", fontweight="bold", transform=plt.gca().transAxes)

        # get the plots to have the same limits on y (for easier comparison)
        global_ylim = [ min(ylim1[0], ylim2[0]), max(ylim1[1], ylim2[1]) ]
        plt.subplot(121) ; plt.ylim(global_ylim)
        plt.subplot(122) ; plt.ylim(global_ylim)

        plt.tight_layout()
        plt.savefig('curvature-fit-angle=%02d.png' % (aa))
        plt.close()





    # ----- FITTING -----

    # storage for nodelist
    nodelist = []

    # M0e data
    m0e_vals = np.array( [ a[sid] for a in m0e_avg ] )
    m0e_errs = np.array( [ b[sid]/np.sqrt(params.impacts) * 1.97 for b in m0e_std ] )
    rvals["m0e_dats"] = m0e_vals
    rvals["m0e_errs"] = m0e_errs
    m0e_node = fitter.FitNode( copy.deepcopy(fp_even), pyam.poly_yamamura, angles, m0e_vals, m0e_errs )
    nodelist.append(m0e_node)

    # M1e data
    m1e_vals = np.array( [ a[sid][0] for a in m1e_avg ] )
    m1e_errs = np.array( [ b[sid][0]/np.sqrt(params.impacts) * 1.97 for b in m1e_std ] )
    rvals["m1e_dats"] = m1e_vals
    rvals["m1e_errs"] = m1e_errs
    m1e_node = fitter.FitNode( copy.deepcopy(fp_odd), pyam.poly_yamamura, angles, m1e_vals, m1e_errs )
    nodelist.append(m1e_node)

    # M1r data
    m1r_vals = np.array( [ a[sid][0] for a in m1r_avg ] )
    m1r_errs = np.array( [ b[sid][0]/np.sqrt(params.impacts) * 1.97 for b in m1r_std ] )
    rvals["m1r_dats"] = m1r_vals
    rvals["m1r_errs"] = m1r_errs
    m1r_node = fitter.FitNode( copy.deepcopy(fp_odd), pyam.poly_yamamura, angles, m1r_vals, m1r_errs )
    nodelist.append(m1r_node)

    # dM0e/dK11 data
    dk11_vals = np.array(dk11_vals)
    dk11_errs = np.array(dk11_errs)
    rvals["dk11_dats"] = dk11_vals
    rvals["dk11_errs"] = dk11_errs
    dk11_node = fitter.FitNode( copy.deepcopy(fp_even), pyam.poly_yamamura, angles, dk11_vals, dk11_errs )
    nodelist.append(dk11_node)

    # dM0e/dK22 data
    dk22_vals = np.array(dk22_vals)
    dk22_errs = np.array(dk22_errs)
    rvals["dk22_dats"] = dk22_vals
    rvals["dk22_errs"] = dk22_errs
    dk22_node = fitter.FitNode( copy.deepcopy(fp_even), pyam.poly_yamamura, angles, dk22_vals, dk22_errs )
    nodelist.append(dk22_node)

    # SuperNode and Fitting
    snode   = fitter.FitSuperNode(nodelist, copy.deepcopy(fp_global))
    result  = snode.fit()
    fparams = result.params

    rvals['m0e_fit'] = m0e_node.result.params
    rvals['m1e_fit'] = m1e_node.result.params
    rvals['m1r_fit'] = m1r_node.result.params
    rvals['dk11_fit'] = dk11_node.result.params
    rvals['dk22_fit'] = dk22_node.result.params
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
  fitfuncD1 = pyam.poly_yamamura_D1

  # get fitted values
  results_list = fit_moments_1D(params, angles)

  # calculate PDE coefficients in 1D
  for rvals in results_list:

    rvals["finedeg"] = finedeg
    finerad = finedeg * np.pi / 180.

    ### ilinov. Now the arrays _vals and _coeffs will contain 1 std_err uncertanties
    # record fitted values of the functions, themselves
    rvals["m0e_vals"] = fitfuncD0(rvals["m0e_fit"], finedeg)
    rvals["m1e_vals"] = fitfuncD0(rvals["m1e_fit"], finedeg)
    rvals["m1r_vals"] = fitfuncD0(rvals["m1r_fit"], finedeg)
    rvals["m1_vals"]  = rvals["m1e_vals"] + rvals["m1r_vals"]

    # record fitted values of the derivatives of M1
    rvals["m0ep_vals"] = fitfuncD1(rvals["m0e_fit"], finedeg)
    rvals["m1ep_vals"] = fitfuncD1(rvals["m1e_fit"], finedeg)
    rvals["m1rp_vals"] = fitfuncD1(rvals["m1r_fit"], finedeg)
    rvals["m1p_vals"]  = rvals["m1ep_vals"] + rvals["m1rp_vals"]

    # construct SX, using the approximation for explicit curvature dependence
    rvals["sxe_coeffs"] =   np.cos(finerad) * rvals["m1ep_vals"] - np.sin(finerad) * rvals["m1e_vals"]
    rvals["sxr_coeffs"] =   np.cos(finerad) * rvals["m1rp_vals"] - np.sin(finerad) * rvals["m1r_vals"]
    rvals["sxc_coeffs"] = - np.cos(finerad) * rvals["m1ep_vals"] / 2.0
    rvals["sx_coeffs"]  = rvals["sxe_coeffs"] + rvals["sxr_coeffs"] + rvals["sxc_coeffs"]

    # construct SY, using the approximation for explicit curvature dependence
    rvals["sye_coeffs"]  =   copy.deepcopy( rvals["sxe_coeffs"] )
    rvals["syr_coeffs"]  =   copy.deepcopy( rvals["sxr_coeffs"] )
    rvals["syc_coeffs"]  =   copy.deepcopy( rvals["sxc_coeffs"] )
    rvals["sye_coeffs"][1:] =   np.cos(finerad[1:])**2 / np.sin(finerad[1:]) * rvals["m1e_vals"][1:] 
    rvals["syr_coeffs"][1:] =   np.cos(finerad[1:])**2 / np.sin(finerad[1:]) * rvals["m1r_vals"][1:] 
    rvals["syc_coeffs"][1:] = - np.cos(finerad[1:])**2 / np.sin(finerad[1:]) * rvals["m1e_vals"][1:] / 2.0
    rvals["sy_coeffs"]  = rvals["sye_coeffs"] + rvals["syr_coeffs"] + rvals["syc_coeffs"]
    

  # change to the path
  if path:
    os.chdir( oldpath )

  # construct the sum dictionary
  return results_list
  




def linked_PDE_coefficients_2D(wrapper, params, angles, curvatures, fitmethod=None, fitguess=None, path=None):

  '''
  This function calculates *flat-target* estimates of the values of the coefficients
  through second order in the linearized equation of motion.  It provides the total
  values of these coefficients, and also the contribution from each component (if
  the target is a multi-component target).  
  
  TODO:  Provide the *total* value within the zeroth entry field.
  TODO:  Modify SX and SY to take into account estimate of curvature correction.
  TODO:  Also calculate the coefficient of the first-order hx term.
  '''

  # functions we will use for fitting
  fitfuncD0 = pyam.poly_yamamura
  fitfuncD1 = pyam.poly_yamamura_D1

  # get fitted values
  results_list = fit_moments_2D(params, angles, curvatures)

  # calculate PDE coefficients in 1D
  for rvals in results_list:

    finedeg = np.linspace(0,90,91)
    rvals["finedeg"] = finedeg
    finerad = finedeg * np.pi / 180.

    # record fitted values of the functions, themselves
    rvals["m0e_vals"] = fitfuncD0(rvals["m0e_fit"], finedeg)
    rvals["m1e_vals"] = fitfuncD0(rvals["m1e_fit"], finedeg)
    rvals["m1r_vals"] = fitfuncD0(rvals["m1r_fit"], finedeg)
    rvals["m1_vals"]  = rvals["m1e_vals"] + rvals["m1r_vals"]
    rvals["dk11_vals"] = fitfuncD0(rvals["dk11_fit"], finedeg)
    rvals["dk22_vals"] = fitfuncD0(rvals["dk22_fit"], finedeg)

    # record fitted values of the derivatives of M1
    rvals["m1ep_vals"] = fitfuncD1(rvals["m1e_fit"], finedeg)
    rvals["m1rp_vals"] = fitfuncD1(rvals["m1r_fit"], finedeg)
    rvals["m1p_vals"]  = rvals["m1ep_vals"] + rvals["m1rp_vals"]

    # construct SX, using the explicit curvature dependence
    rvals["sxe_coeffs"] =  np.cos(finerad) * rvals["m1ep_vals"] - np.sin(finerad) * rvals["m1e_vals"]
    rvals["sxr_coeffs"] =  np.cos(finerad) * rvals["m1rp_vals"] - np.sin(finerad) * rvals["m1r_vals"]
    rvals["sxc_coeffs"] =  np.cos(finerad) * rvals["dk11_vals"]
    rvals["sxc_coeffs_approx"] = - np.cos(finerad) * rvals["m1ep_vals"] / 2.0
    rvals["sx_coeffs"]  = rvals["sxe_coeffs"] + rvals["sxr_coeffs"] + rvals["sxc_coeffs"]

    # construct SY, using the explicit curvature dependence
    rvals["syc_coeffs"] =  np.cos(finerad) * rvals["dk22_vals"]

    rvals["sye_coeffs"]  =   copy.deepcopy( rvals["sxe_coeffs"] )
    rvals["syr_coeffs"]  =   copy.deepcopy( rvals["sxr_coeffs"] )
    rvals["syc_coeffs_approx"]  =   copy.deepcopy( rvals["sxc_coeffs_approx"] )

    rvals["sye_coeffs"][1:] =   np.cos(finerad[1:])**2 / np.sin(finerad[1:]) * rvals["m1e_vals"][1:] 
    rvals["syr_coeffs"][1:] =   np.cos(finerad[1:])**2 / np.sin(finerad[1:]) * rvals["m1r_vals"][1:] 
    rvals["syc_coeffs_approx"][1:] = - np.cos(finerad[1:])**2 / np.sin(finerad[1:]) * rvals["m1e_vals"][1:] / 2.0

    rvals["sy_coeffs"]  = rvals["sye_coeffs"] + rvals["syr_coeffs"] + rvals["syc_coeffs"]


  # change to the path
  if path:
    os.chdir( oldpath )

  # construct the sum dictionary
  return results_list
  







# --------------------------------------------------
#    plot_angle_dependence_summary()
# --------------------------------------------------


def plot_single_flat_angle_dependence_summary(fitted_values):

  rvals = fitted_values

  # plot the data and the fit
  ### ilinov. All arrays with _vals and _coeffs contain uncertanties information. the initial errors originate
  ### from the moments fitting with the Yamamura polynomial -> we should get 1 standard error confidence bounds 
  theplot = plt.figure(figsize=(10,5.5))
  plt.rcParams.update({'axes.titlesize': 12})
  plt.rcParams.update({'axes.labelsize': 9})
  plt.rcParams.update({'xtick.labelsize': 8})
  plt.rcParams.update({'ytick.labelsize': 8})
  textxloc = 0.9
  textyloc = 0.875

  plt.subplot(221)
  plt.errorbar(rvals["angles"],  -rvals["m0e_dats"], yerr=rvals["m0e_errs"], fmt='gs', label='data')
  plt.plot    (rvals["finedeg"], -unp.nominal_values(rvals["m0e_vals"]), 'g-', linewidth=2, label="fit")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$M^{(0)}_{\mathsf{eros.}} \left( \theta \right)$   [atom / ion]')
  plt.legend(loc=2,prop={'size':9})
  plt.title('Sputter Yield')
  plt.text(textxloc, textyloc, "(a)", fontweight="bold", transform=plt.gca().transAxes)

  plt.subplot(222)
  plt.errorbar(rvals["angles"], rvals["m1e_dats"] / 10.0, yerr=rvals["m1e_errs"] / 10.0, fmt='rs', label='eros.')
  plt.plot    (rvals["finedeg"], unp.nominal_values(rvals["m1e_vals"]) / 10.0, 'r-', linewidth=2)
  plt.errorbar(rvals["angles"], rvals["m1r_dats"] / 10.0, yerr=rvals["m1r_errs"] / 10.0, fmt='bs', label='redist.')
  plt.plot    (rvals["finedeg"], unp.nominal_values(rvals["m1r_vals"]) / 10.0, 'b-', linewidth=2)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$M^{(1)} \left( \theta \right)$   [atom * nm / ion]')
  plt.legend(loc=3,prop={'size':9})
  plt.title('First Moments')
  plt.text(textxloc, textyloc, "(b)", fontweight="bold", transform=plt.gca().transAxes)

  plt.subplot(223)
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["sxe_coeffs"]) / 10.0, 'r-', linewidth=2, label="eros.")
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["sxr_coeffs"]) / 10.0, 'b-', linewidth=2, label="redist.")
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["sxc_coeffs"]) / 10.0, 'g-', linewidth=2, label="curv.*")
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["sx_coeffs"])  / 10.0,  'k--', linewidth=2, label="total")
  ### ilinov: add 1 std err confidence bounds
  xx = rvals["sxe_coeffs"]
  plt.fill_between(rvals["finedeg"], (unp.nominal_values(xx)-unp.std_devs(xx))/10.0, (unp.nominal_values(xx)+unp.std_devs(xx))/10.0, facecolor='r', alpha=0.3)
  xx = rvals["sxr_coeffs"]
  plt.fill_between(rvals["finedeg"], (unp.nominal_values(xx)-unp.std_devs(xx))/10.0, (unp.nominal_values(xx)+unp.std_devs(xx))/10.0, facecolor='b', alpha=0.3)
  xx = rvals["sxc_coeffs"]
  plt.fill_between(rvals["finedeg"], (unp.nominal_values(xx)-unp.std_devs(xx))/10.0, (unp.nominal_values(xx)+unp.std_devs(xx))/10.0, facecolor='g', alpha=0.3)
  xx = rvals["sx_coeffs"]
  plt.fill_between(rvals["finedeg"], (unp.nominal_values(xx)-unp.std_devs(xx))/10.0, (unp.nominal_values(xx)+unp.std_devs(xx))/10.0, facecolor='k', alpha=0.3)
  ###  
  ylimits2 = np.array(plt.ylim())
  ylimits2[0] = -ylimits2[1]
  plt.ylim(ylimits2)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{X} = C_{11} \left( \theta \right)$   [atom * nm / ion]')
  plt.legend(loc=3,prop={'size':9}, ncol=2)
  plt.title('Components of $S_{X}=C_{11}$')
  ymin, ymax = plt.ylim()
  plt.text(textxloc, textyloc, "(c)", fontweight="bold", transform=plt.gca().transAxes)

  plt.subplot(224)
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["sye_coeffs"]) / 10.0, 'r-', linewidth=2, label="eros.")
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["syr_coeffs"]) / 10.0, 'b-', linewidth=2, label="redist.")
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["syc_coeffs"]) / 10.0, 'g-', linewidth=2, label="curv.*")
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["sy_coeffs"])  / 10.0,  'k--', linewidth=2, label="total")
  ### ilinov: add 1 std err confidence bounds
  xx = rvals["sye_coeffs"]
  plt.fill_between(rvals["finedeg"], (unp.nominal_values(xx)-unp.std_devs(xx))/10.0, (unp.nominal_values(xx)+unp.std_devs(xx))/10.0, facecolor='r', alpha=0.3)
  xx = rvals["syr_coeffs"]
  plt.fill_between(rvals["finedeg"], (unp.nominal_values(xx)-unp.std_devs(xx))/10.0, (unp.nominal_values(xx)+unp.std_devs(xx))/10.0, facecolor='b', alpha=0.3)
  xx = rvals["syc_coeffs"]
  plt.fill_between(rvals["finedeg"], (unp.nominal_values(xx)-unp.std_devs(xx))/10.0, (unp.nominal_values(xx)+unp.std_devs(xx))/10.0, facecolor='g', alpha=0.3)
  xx = rvals["sy_coeffs"]
  plt.fill_between(rvals["finedeg"], (unp.nominal_values(xx)-unp.std_devs(xx))/10.0, (unp.nominal_values(xx)+unp.std_devs(xx))/10.0, facecolor='k', alpha=0.3)
  ###
  plt.ylim(ylimits2)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{Y} = C_{22} \left( \theta \right)$   [atom * nm / ion]')
  plt.legend(loc=3,prop={'size':9},ncol=2)
  plt.title(r'Components of $S_{Y}=C_{22}$')
  plt.ylim(ymin, ymax)
  plt.text(textxloc, textyloc, "(d)", fontweight="bold", transform=plt.gca().transAxes)

  #plt.subplot(313)
  #plt.plot(rvals["finedeg"], rvals["sx_coeffs"], 'c-', linewidth=2, label="SX")
  #plt.plot(rvals["finedeg"], rvals["sy_coeffs"], 'm-', linewidth=2, label="SY")
  #plt.plot([0, 90], [0, 0], 'k--', linewidth=1)
  #plt.xlabel(r'angle $\theta$')
  #plt.ylabel(r'$S_{X,Y} \left( \theta \right)$')
  #plt.legend(loc=3,prop={'size':9})
  #plt.title(r'$S_{X}$ vs $S_{Y}$')
  plt.tight_layout()

  return theplot







def plot_single_curved_angle_dependence_summary(fitted_values):

  # rename the argument
  rvals = fitted_values

  # plot the data and the fit
  theplot = plt.figure(figsize=(10,5.5))
  plt.rcParams.update({'axes.titlesize': 10})
  plt.rcParams.update({'axes.labelsize': 9})
  plt.rcParams.update({'xtick.labelsize': 8})
  plt.rcParams.update({'ytick.labelsize': 8})
  textxloc = 0.9
  textyloc = 0.875


  plt.subplot(221)
  plt.errorbar(rvals["angles"], rvals["m1e_dats"] / 10.0, yerr=rvals["m1e_errs"] / 10.0, fmt='rs', label='eros.')
  plt.errorbar(rvals["angles"], rvals["m1r_dats"] / 10.0, yerr=rvals["m1r_errs"] / 10.0, fmt='bs', label='redist.')
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["m1e_vals"]) / 10.0, 'r-', linewidth=2)
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["m1r_vals"]) / 10.0, 'b-', linewidth=2)
  plt.xlim((0,90))
  ylimits1 = plt.ylim()
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$M^{(1)} \left( \theta \right)$  [atom nm / ion]')
  plt.legend(loc=3,prop={'size':9})
  plt.title('Erosive and Redistributive First Moments')
  plt.text(textxloc, textyloc, "(a)", fontweight="bold", transform=plt.gca().transAxes)


  plt.subplot(222)
  #plt.errorbar(rvals["angles"], rvals["m0e_dats"],    yerr=rvals["m0e_errs"],    fmt='bs', label='flat')
  #plt.errorbar(rvals["angles"], rvals["m0ek1p_avg"], yerr=rvals["m0ek1p_err"], fmt='gs', label='k11=%0.3f'%(dk))
  #plt.errorbar(rvals["angles"], rvals["m0ek1m_avg"], yerr=rvals["m0ek1m_err"], fmt='rs', label='k11=%0.3f'%(-dk))
  plt.errorbar(rvals["angles"], rvals["dk11_dats"] / 10.0, yerr=rvals['dk11_errs'] / 10.0, fmt='gs', label='K11')
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["dk11_vals"]) / 10.0, 'g-', linewidth=2)
  #plt.plot(rvals["finedeg"], rvals["m0e_vals"], 'b-', linewidth=2)
  #plt.plot(rvals["finedeg"], rvals["m0ek1p_vals"], 'g-', linewidth=2)
  #plt.plot(rvals["finedeg"], rvals["m0ek1m_vals"], 'r-', linewidth=2)


  #plt.errorbar(rvals["angles"], rvals["m0e_dats"],    yerr=rvals["m0e_errs"],    fmt='bs', label='flat')
  #plt.errorbar(rvals["angles"], rvals["m0ek2p_avg"], yerr=rvals["m0ek2p_err"], fmt='gs', label='k22=%0.3f'%(dk))
  #plt.errorbar(rvals["angles"], rvals["m0ek2m_avg"], yerr=rvals["m0ek2m_err"], fmt='rs', label='k22=%0.3f'%(-dk))
  plt.errorbar(rvals["angles"], rvals["dk22_dats"] / 10.0, yerr=rvals['dk22_errs'] / 10.0, fmt='ys', label='K22')
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["dk22_vals"]) / 10.0, 'y-', linewidth=2)
  #plt.plot(rvals["finedeg"], rvals["m0e_vals"], 'b-', linewidth=2)
  #plt.plot(rvals["finedeg"], rvals["m0ek2p_vals"], 'g-', linewidth=2)
  #plt.plot(rvals["finedeg"], rvals["m0ek2m_vals"], 'r-', linewidth=2)

  plt.ylim(ylimits1)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$\partial M^{(0)}_{\mathsf{eros.}} / \partial K_{ii} \left( \theta \right)$ [atom nm / ion]')
  plt.legend(loc=3,prop={'size':9})
  plt.title('Curvature-Derivatives of Zeroth Moment')
  plt.text(textxloc, textyloc, "(b)", fontweight="bold", transform=plt.gca().transAxes)



  plt.subplot(223)
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["sxe_coeffs"]) / 10.0, 'r-', linewidth=2, label="from $M^{(1)}_{eros.}$")
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["sxr_coeffs"]) / 10.0, 'b-', linewidth=2, label="from $M^{(1)}_{redist.}$")
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["sxc_coeffs"]) / 10.0, 'g-', linewidth=2, label="from $M^{(0)}_{eros.}$")
  #plt.plot(rvals["finedeg"], rvals["sxc_coeffs_approx"], 'y-', linewidth=2, label="curv.**")
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["sx_coeffs"]) / 10.0, 'k--', linewidth=2, label="total")
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{X}=C_{11} \left( \theta \right)$\, [nm$^{4}$ / ion]')
  ylimits2 = np.array(plt.ylim())
  ylimits2[0] = -ylimits2[1]
  plt.ylim(ylimits2)
  plt.legend(loc=3,prop={'size':8},ncol=2)
  plt.title(r'Components of $S_X = C_{11}$ ' )
  plt.text(textxloc, textyloc, "(c)", fontweight="bold", transform=plt.gca().transAxes)

  plt.subplot(224)
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["sye_coeffs"]) / 10.0, 'r-', linewidth=2, label="from $M^{(1)}_{eros.}$")
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["syr_coeffs"]) / 10.0, 'b-', linewidth=2, label="from $M^{(1)}_{redist.}$")
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["syc_coeffs"]) / 10.0, 'y-', linewidth=2, label="from $M^{(0)}_{eros.}$")
  #plt.plot(rvals["finedeg"], rvals["syc_coeffs_approx"], 'y-', linewidth=2, label="curv.**")
  plt.plot(rvals["finedeg"], unp.nominal_values(rvals["sy_coeffs"]) / 10.0, 'k--', linewidth=2, label="total")
  plt.ylim(ylimits2)
  plt.xlabel(r'angle $\theta$')
  plt.ylabel(r'$S_{Y}=C_{22} \left( \theta \right)$\, [nm$^{4}$ / ion]')
  plt.legend(loc=3,prop={'size':8},ncol=2)
  plt.title('Components of $S_Y = C_{22}$')
  plt.text(textxloc, textyloc, "(d)", fontweight="bold", transform=plt.gca().transAxes)


  plt.tight_layout()
  return theplot







