'''
This file describes a generic instance of the Statistic object.
'''


class GENERIC_Statistics(object):

  # initialization of data structures
  def __init__(self):
    pass

  # how extract statistics from implanted, eroded, redistributed atoms
  def collect(idata, edata, rdata):
    '''
    This function recieves lists of impact data, erosion data, and
    redistribution data from a single impact as three separate lists,
    and uses these lists of atomic positions to calculate all numbers
    sought by the user.  

    This function must be able to accept (None, None, None) as arguments,
    in which case it returns values of the statistical objects with all
    entries set to be zero.  This is used as an initialization step by 
    the driver code in the wrapper.
    '''
    pass








