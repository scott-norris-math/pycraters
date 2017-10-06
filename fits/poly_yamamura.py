import numpy as np
from uncertainties import ufloat
import uncertainties.umath  
import uncertainties.unumpy as unp  # numpy  ### ilinov


# ------------------------------------------
#    Polynomial-Yamamura
# ------------------------------------------


def poly_yamamura(LMparams, angles):
  # this function defines a generic angular function with Yamamura-like decay.

  yp        = LMparams["p"].value		# Yamamura power
  radt      = angles * np.pi / 180.		# convert angles to radians
  Yterm     = np.exp(1 - 1./np.cos(radt))**yp	# Yamamura decay term

  # construct the function
  model  = np.zeros(np.shape(angles))
  #model = unp.uarray(model, model)
  for pp in range(0,10):
    pstring = "c%1d" % (pp)
    if pstring in LMparams:
      #err = 0.0 if LMparams[pstring].stderr == None else LMparams[pstring].stderr == None 
      #print LMparams[pstring].value, "+/-", LMparams[pstring].stderr
      #coeff = ufloat(LMparams[pstring].value, err)
      coeff = LMparams[pstring].value
      model  += coeff * radt**pp * Yterm

  return model



def poly_yamamura_D1(LMparams, angles):
  # this function returns the angle-derivative of a generic curve with Yamamura-like decay.

  yp        = LMparams["p"].value			# Yamamura power
  radt      = angles * np.pi / 180.	# convert angles to radians
  Yterm     = np.exp(1 - 1./np.cos(radt))**yp	# Yamamura decay term
  YPterm    = Yterm * -yp * np.sin(radt) / np.cos(radt)**2 

  model  = np.zeros(len(angles))
  model = unp.uarray(model, model)
#  print np.cos(radt)**2
  for pp in range(0,10):
    pstring = "c%1d" % (pp)
    if pstring in LMparams:
      err = 0.0 if LMparams[pstring].stderr == None else LMparams[pstring].stderr == None 
      coeff = ufloat(LMparams[pstring].value, err)
            
      if (pp >  0):  model += coeff * pp * radt**(pp-1) * Yterm
      if (pp >= 0):  model += coeff * radt**pp * YPterm 

  return model




