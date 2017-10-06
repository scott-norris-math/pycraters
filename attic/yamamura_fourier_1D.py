import numpy as np


# ------------------------------------------
#    Polynomial-Yamamura
# ------------------------------------------


def yamamura_fourier(params, angles, options):
  # this function defines a generic curve with Yamamura-like decay.
  # useful for fitting moments as a function of angle.

  radt      = angles * np.pi / 180.	# convert angles to radians
  vals      = np.zeros(np.shape(angles))	# storage for evaluated function
  symmetry  = options["symmetry"]	# a parameter dictating even or odd symmetry
  clist     = options["coeff_list"]	# list of strings, giving names of coefficients
  pp        = params["p"]		# Yamamura power


  # loop through terms, creating the function
  for kk,cc in enumerate(clist):
    coeff = params[cc]

    if symmetry == "odd":
      nn = 2*kk+1
      vals += coeff * np.sin(nn*radt) * (np.exp(-1./np.cos(radt)) / np.cos(radt))**pp

    if symmetry == "even":
      nn = 2*kk
      vals += coeff * np.cos(nn*radt) * (np.exp(-1./np.cos(radt)) / np.cos(radt))**pp

  return vals



def yamamura_fourier_derivative(params, angles, options):
  # this function returns the angle-derivative of a generic curve with Yamamura-like decay.
  # useful for calculating coefficients that depend on the angle derivative of moments.

  radt      = angles * np.pi / 180.	# convert angles to radians
  vals      = np.zeros(np.shape(angles))	# storage for evaluated function
  symmetry  = options["symmetry"]	# a parameter dictating even or odd symmetry
  clist     = options["coeff_list"]	# list of strings, giving names of coefficients
  pp        = params["p"]		# Yamamura power


  # loop through terms, creating the function
  for kk,cc in enumerate(clist):
    coeff = params[cc]

    if symmetry == "odd":  
      nn = 2*kk+1
      vals += coeff * (np.exp(-1./np.cos(radt)) / np.cos(radt))**pp * (			\
            nn*np.cos(nn*radt) + np.sin(nn*radt)*pp*np.tan(radt)*(1 - 1./np.cos(radt))	\
            )

    if symmetry == "even":
      nn = 2*kk
      vals += coeff * (np.exp(-1./np.cos(radt)) / np.cos(radt))**pp * (			\
            -nn*np.sin(nn*radt) + np.cos(nn*radt)*pp*np.tan(radt)*(1 - 1./np.cos(radt))	\
            )

  return vals



def poly_yamamura_cerror(params):
  # this function defines a penalty for Yamamura coefficients that are negative
  if params['p'] >= 0:  return 0
  if params['p'] <  0:  return 10**9*params['p']**2



