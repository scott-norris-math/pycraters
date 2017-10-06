import numpy as np


def zero_polynomial(params, xi):

  terms  = len(params.keys())
  coeffs = np.empty(terms)
  model  = np.zeros(len(xi))
  for pp in range(1,terms+1):
    pstring = "c%1d" % (pp)
    coeff   = params[pstring].value
    model  += coeff * xi**pp

  return model


def zero_polynomial_D1(params, xi):

  terms  = len(params.keys())
  coeffs = np.empty(terms)
  model  = np.zeros(len(xi))
  for pp in range(1,terms+1):
    pstring = "c%1d" % (pp)
    coeff   = params[pstring].value
    model  += coeff * pp * xi**(pp-1)

  return model


 

def zero_polynomial_residual(params, xi, data, errs=None):

  if errs != None:
    return (zero_polynomial(params, xi) - data) / errs
  else:
    return (zero_polynomial(params, xi) - data)





  
