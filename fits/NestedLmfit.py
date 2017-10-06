import copy
import numpy as np
from lmfit import Parameter, Parameters, minimize


class FitNode(object):
  '''
  This class stores
  - a single lmfit.Parameters() object
  - a single supplied fitting function 
  - a list of indepdnent variable values
  - a list of dependent variable values
  - a list of uncertanty values (optional)
  and then optimizes the parameters using lmfit.minimize().

  Optionally, the standard errors of the data can be provided to guide the fit. 
  '''

  def __init__(self, params, ffunc, xvals, yvals, uvals=None):

    self.params = params
    self.ffunc  = ffunc
    self.xvals  = xvals
    self.yvals  = yvals
    self.uvals  = uvals
    self.result = None


  def evaluate(self, params=None, xvals=None):
    # evaluate the fitting function for a given set of independent variables (useful for plotting)

    if params == None:  params = self.result.params
    if xvals  == None:  xvals  = self.xvals
    return self.ffunc(params, xvals)


  def residual(self, params=None, xvals=None, yvals=None, uvals=None):
    # simple residual function with optional uncertainty weighting

    if params == None:  params = self.params
    if xvals  == None:  xvals  = self.xvals
    if yvals  == None:  yvals  = self.yvals
    if uvals  == None:  uvals  = self.uvals

    if uvals  == None:
      return (self.ffunc(params, xvals) - yvals)
    else:
      avgerror = np.mean(uvals[np.nonzero(uvals)])
      uvals[np.where(uvals==0)] = avgerror
      return (self.ffunc(params, xvals) - yvals) / uvals


  def fit(self):
    guess = copy.deepcopy(self.params)
    result = minimize(self.residual, guess, args=(self.xvals, self.yvals, self.uvals), method="leastsq", **{'xtol':1e-8})
    #print result.params
    self.result = result
    return result








class FitSuperNode(object):
  '''
  This class stores
  - a list of FitNodes (or other FitSuperNodes!)
  - a global Parameters() object that applies to all nodes
  and then optimizes the global parameters using lmfit.minimize().
  During this process, the local parameters in each FitNode are also optimized.
  '''

  def __init__(self, nodelist=None, params=None):

    self.nodelist = nodelist
    self.params   = params
    self.result   = None


  def residual(self, params=None):

    if params == None:  params = self.params
    #print params

    residual_list = []
    for node in self.nodelist:
      for key in params.keys():				# add global params to each local node
        node.params.add(key, value=params[key].value, vary=False)  
      result = node.fit()				# fit each node
      residual_list.append(np.array(result.residual))	# get local residual

    globalres = np.hstack(residual_list)
    return globalres



  def fit(self):

    if self.params == None:
      for node in self.nodelist: node.fit()
      return None

      
    guess = copy.deepcopy(self.params)
    result = minimize(self.residual, guess, method="leastsq", **{'xtol':1e-8})
    self.result = result

    for node in self.nodelist:
      for key in self.result.params.keys():
        node.result.params[key].value = self.result.params[key].value
        node.result.params[key].stderr = self.result.params[key].stderr
    return self.result

  
    




