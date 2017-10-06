import numpy as np


# ------------------------------------------
#    Polynomial-Yamamura
# ------------------------------------------


def poly_yamamura(params, angles, curvatures, options):

  # this function defines a generic curve with Yamamura-like decay.
  # useful for fitting moments as a function of angle.

  radt      = angles * np.pi / 180.	# convert angles to radians
  vals      = np.zeros(np.shape(angles))	# storage for evaluated function
  symmetry  = options["symmetry"]	# a parameter dictating even or odd symmetry
  clist     = options["coeff_list"]	# list of strings, giving names of coefficients
  pp        = params["p"]		# Yamamura power

  # an offset for the value of each power  
  power_adjustment = 1 if symmetry == "odd" else 0

  # loop through power, creating the function
  for kk,cc in enumerate(clist):
    coeff = params[cc]
    power = 2*kk + power_adjustment
    vals += coeff * radt**power * (np.exp(-1./np.cos(radt)) / np.cos(radt))**pp

  return vals











def poly_yamamura_derivative(params, angles, options):
  # this function returns the angle-derivative of a generic curve with Yamamura-like decay.
  # useful for calculating coefficients that depend on the angle derivative of moments.

  radt      = angles * np.pi / 180.	# convert angles to radians
  vals      = np.zeros(np.shape(angles))	# storage for evaluated function
  symmetry  = options["symmetry"]	# a parameter dictating even or odd symmetry
  clist     = options["coeff_list"]	# list of strings, giving names of coefficients
  pp        = params["p"]		# Yamamura power


  # an offset for the value of each power  
  power_adjustment = 1 if symmetry == "odd" else 0

  # loop through power, creating the function
  for kk,cc in enumerate(clist):
    coeff = params[cc]
    power = 2*kk + power_adjustment
    vals += coeff * (np.exp(-1./np.cos(radt)) / np.cos(radt))**pp * (				\
            power * radt**(power-1) + pp*radt**power*np.tan(radt)*(1 - 1./np.cos(radt))	\
            )

  return vals



def poly_yamamura_cerror(params):
  # this function defines a penalty for Yamamura coefficients that are negative
  if params['p'] >= 0:  return 0
  if params['p'] <  0:  return 10**9*params['p']**2



