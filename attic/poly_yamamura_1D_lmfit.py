import numpy as np


# ------------------------------------------
#    Polynomial-Yamamura
# ------------------------------------------


def poly_yamamura(LMparams, angles):
  # this function defines a generic angular function with Yamamura-like decay.

  radt      = angles * np.pi / 180.	# convert angles to radians
  yp        = LMparams["p"]		# Yamamura power

  # construct the function
  model  = np.zeros(np.shape(angles))
  for pp in range(0,10):
    pstring = "c%1d" % (pp)
    if pstring in LMparams:
      coeff   = LMparams[pstring].value
      model  += coeff * radt**pp * (np.exp(-1./np.cos(radt)) / np.cos(radt))**yp

  return model



def poly_yamamura_D1(LMparams, angles):
  # this function returns the angle-derivative of a generic curve with Yamamura-like decay.

  radt      = angles * np.pi / 180.	# convert angles to radians
  yp        = LMparams["p"]		# Yamamura power

  model     = np.zeros(np.shape(angles))
  for pp in range(0,10):
    pstring = "c%1d" % (pp)
    if pstring in LMparams:
      coeff   = LMparams[pstring].value
      model  += coeff * (np.exp(-1./np.cos(radt)) / np.cos(radt))**yp * (		\
            pp * radt**(pp-1) + yp*(radt**pp)*np.tan(radt)*(1 - 1./np.cos(radt))	\
            )

  return model




