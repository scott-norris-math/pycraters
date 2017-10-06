import numpy as np


def polynomial(LMparams, xi):

  model  = np.zeros(len(xi))
  for pp in range(0,10):
    pstring = "c%1d" % (pp)
    if pstring in LMparams:
      coeff   = LMparams[pstring].value
      model  += coeff * xi**pp

  return model


def polynomial_D1(LMparams, xi):

  model  = np.zeros(len(xi))
  for pp in range(1,10):
    pstring = "c%1d" % (pp)
    if pstring in LMparams:
      coeff   = LMparams[pstring].value
      model  += coeff * pp * xi**(pp-1)

  return model



def polynomial_residual(params, xi, data, errs=None):

  if errs != None:
    return (polynomial(params, xi) - data) / errs
  else:
    return (polynomial(params, xi) - data)





