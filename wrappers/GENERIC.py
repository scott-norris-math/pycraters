'''
IDEA:

parameter blob definitions

I/O code for SDTrimSP (running simulations, reading output files, storing statistics)
I/O code for PARCAS (...., reading output files, storing statistics)
I/O code for ***

reading / writing moments to standard format
fitting and coefficient extraction code
plotting code


then use this code 
'''



import os
import sys
import shelve
import numpy as np
from subprocess import call


import pycraters.statistics.base_moments as bm





# ----------------------------------------------------
#       Generic Parameter Object Definition
# ----------------------------------------------------






class EnvironmentParameterError(Exception):
  def __init__(self, string):
    self.string = string
  def __str__(self):
    return self.string





class GeometryData(object):
  '''
  Storage for description of the target, and binning
  '''

  def __init__(self):
    self.range	= None
    self.window	= None
    self.bins	= None






class GENERIC_Parameters(object):
  '''
  Storage for various kinds of parameters
  '''

  # empty initialization routine
  def __init__(self):

    # text string for diagnostic messages
    self.description = "Generic BCA Environmental Parameters"

    # parameters related to the beam  (PRIMARY)
    self.beam	 = None  # should be a single text string
    self.energy	 = None  # should be a single float
    self.angle	 = None  # should be a single float
    self.impacts = None # should be an integer

    # parameters describing the target composition and energetics
    self.target	 = None	# should be a list of target/concentration pairs

    # BCA parameters describing the target energies
    self.SBE_mat = None	# a matrix of size nxn (where n is the number of species)
    self.SBE_vec = None	# a vector of size n (if matrix not supplied)
    self.BBE_vec = None	# a vector of size n
    self.ethresh = None	# a vector of size n
    self.ecutoff = None	# a vector of size n

    # parameters describing the target geometry
    self.k11     = 0.0	# curvature in x-direction
    self.k12	 = 0.0	# mixed second derivative (zero if irradiating in principle curvature direction)
    self.k22     = 0.0	# curvature in y-direction

    #  Dictionary for Additional Parameters (solver-specific)
    self.additional_parameters = dict()




    # deprecated functions
    # -----------------------------------
    # simulation domain geometry  (deprecated, for now)
    self.geometry = GeometryData() #adds a geometry object to help integrate with BCA Wrapper
    self.geometry.xwidth = 100.0
    self.geometry.ywidth = 100.0
    self.geometry.zwidth = 100.0
    self.geometry.xbins = 20
    self.geometry.ybins = 20
    self.geometry.zbins = 20

    # Deprecated (for now -- was used for z-dependent concentrations in SDTRIMSP -- should be moved to solver-specific parameters)
    self.depth  = None  	# should be a single float 
    self.dep_res = 1    	# should be an integer (depth resolution, default 1)
    self.cfunc  = None  	# should be a function of a float which returns a list of size (# of targets)
    self.funcname = None	# should be a string uniquely describing a cfunc






  def set_parameter(self, key, value, fstring=None):
    '''
    THIS FUNCTION SHOULD NOT NEED TO BE MODIFIED BY ANY DAUGHTER OBJECTS!
    '''
    if key in self.additional_parameters.keys():
      self.additional_parameters[key] = value
    else:
      raise EnvironmentParameterError("Environment Parameter %s not found in ParameterHolder of type: %s" % (key, self.description))

    '''
    Here we need a mechanism to specify how the specified additional parameter 
    should affect the name of the file containing the moments.  One might well
    want to sweep over a parameter that can only be specified in this function,
    and the resulting file names need to be unique.

    Probably, we should err on the side of filenames that are too long -- i.e.,
    any additional parameter specified by the user will appear in the filename,
    unless the user explicitly forbids it.
    '''

    return



  def fname(self, keys=None):
    '''
    Builds a unique filename from the parameter blob.
    Used for saving or loading files when doing loops
    over many different parameter values.
    '''

    if keys == None:  keys = ['target', 'beam', 'energy', 'angle', 'k11', 'k22']

    fstring = ''

    if 'target' in keys:
      for ts in self.target:
        spec = ts[0]
        conc = ts[1]
        cstr = '' if (conc == 1.0) else "%02d" % (round(conc*100))
        fstring += '%s%s-' % (spec, cstr)

    if 'beam' in keys:
      fstring += '%s-' % (self.beam)

    if 'energy' in keys:
      fstring += '%05deV-' % (int(self.energy))

    if 'angle' in keys:
      fstring += '%02ddeg-' % (int(self.angle))

    if 'k11' in keys:
      fstring += '%0.3fK11-' % (self.k11)

    if 'k22' in keys:
      fstring += '%0.3fK22' % (self.k22)


    #fstring = ''

    #for ts in self.target:
    #  spec = ts[0]
    #  conc = ts[1]
    #  cstr = '' if (conc == 1.0) else "%02d" % (round(conc*100))
    #  fstring += '%s%s-' % (spec, cstr)

    #fstring += '%s-' % (self.beam)
    #fstring += '%05deV-' % (int(self.energy))
    #fstring += '%02ddeg-' % (int(self.angle))
    #fstring += '%0.3fK11-' % (self.k11)
    #fstring += '%0.3fK22' % (self.k22)
    #fstring += '%02dedisp' % (int(self.edispl))

    '''
    NEEDED:  A way to include additional, non-default parameters in the filename
    (i.e., binding energies, displacement energies, etc..).  Even better, the user
    may eventually define parametric dependencies of certain functions themeselves,
    and there should be a way to record values of that parameter in the filename.

    UPDATE:  The above-described functionality should work together with set_parameter().
    '''
    return fstring    




  def load_default_energies():
    '''
    At some point, we should have this function load "default" values of the

      -- surface binding energy
      -- bulk binding energy
      -- displacement threshold energy
      -- electronic cutoff energy

    for each species in the target (and possibly also in the beam).  These should
    be extracted from some master database file in a human-readable format, that 
    can be included in the distribution and edited over time as better data become
    available.  The database should include citation information, etc...

    If the user wishes to change any of this information, then there should be an 
    interface to do so, and a way to inform the target program to update its own, 
    private list (if such a list is used internally).

    SQLite???  with Java-based SQuirreL?  with OpenOffice for forms?  or with
    Views / "instead of" triggers using raw SQL?  This could quickly become a 
    rabbit hole ...
    '''
    pass


  def __str__(self):
      return repr(self.value)



















# ----------------------------------------------------
#       Generic Simulation Wrapper Definition
# ----------------------------------------------------







class InputError(Exception):
  def __init__(self):
    pass
  def __str__(self):
    return "Location of executable file not provided!"



class SimulationParameterError(Exception):
  def __init__(self, string):
    self.string = string
  def __str__(self):
    return self.string




class GENERIC_Impact_Iterator():

  '''
  This method class needs to needs to read the output of the program,
  and return an object ITER that, upon each call to ITER.next_impact(),
  returns three lists: one each for implantation, erosive, and 
  redistributive data
  '''


  def next_impact():
    '''
    This method needs to identify the next set of impact information,
    and put it into a set of standard formats:


    implantation data:

          list[ N-list[ int, 3-array ] ]

    where N is the number of implanted particles, and for each particle we store, in order: 
    -- the ID number of the atom species
    -- and the final position of the particle.


    erosive data:

          list[ N-list[ int, 3-array ] ]

    where N is the number of sputtered particles, and for each particle we store, in order: 
    -- the ID number of the atom species
    -- and the initial position of the particle.


    redistributive data:

          list[ N-list[ int, 3-array, 3-array ] ]

    where N is the number of sputtered particles, and for each particle we store , in order:
    -- the ID number of the atom species
    -- the initial position of the particle, 
    -- and the final position, in that order.
    '''








class GENERIC_Wrapper(object):
  '''
  This is essentially an interface definition.
  Each solver will extend this class.
  '''



  # ------------------------------------------------------
  # Set things up so that the simulation can be run
  # ------------------------------------------------------
  def __init__(self):
    '''
    This needs to do everything needed to *allow* the simulation to be run.
    This may include 
    -- storing the location of the executable 
    -- storing the command needed to run the executable
    -- storing the location of any needed data files
    -- storing any instructions regarding capabilities that are specific to this solver
    '''

    self.wrapper_type = "Generic, Nonspecific BCA Wrapper"
    self.execline = ''

    self.simulation_parameters = dict()

    self.statistics_routines = list()
    self.add_statistics_routine(bm.base_moments())

    pass




  def add_statistics_routine(self, routine):
    '''
    THIS FUNCTION SHOULD NOT NEED TO BE MODIFIED BY ANY DAUGHTER OBJECTS!
    '''
    self.statistics_routines.append(routine)




  def set_parameter(self, key, value):
    '''
    THIS FUNCTION SHOULD NOT NEED TO BE MODIFIED BY ANY DAUGHTER OBJECTS!
    '''
    if key in self.simulation_parameters.keys():
      self.simulation_parameters[key] = value
    else:
      raise SimulationParameterError("Parameter %s not found in wrapper of type: %s" % (key, self.wrapper_type))
    return




  # ------------------------------------------------------
  #    single top-level running method
  # ------------------------------------------------------
  def go(self, params, save_raw_data=True):
    '''
    THIS FUNCTION SHOULD NOT NEED TO BE MODIFIED BY ANY DAUGHTER OBJECTS!
    '''

    if (self.check_for_data(params) == False):	# see if data already exists 
      print "  running simulations."
      sys.stdout.flush()
      self.run_simulation(params)		# actually runs the simulation
      print "  extracting statistics."
      sys.stdout.flush()
      self.extract_statistics(params)		# extracts the moments from raw data files
      if (save_raw_data == False):
        self.clean_up(params)			# deletes the raw data

    return 



  def check_for_data(self, params):
    '''
    THIS FUNCTION SHOULD NOT NEED TO BE MODIFIED BY ANY DAUGHTER OBJECTS!
    '''

    targetfile = "%s.moms" % (params.fname())
    return os.path.isfile(targetfile)





  # ------------------------------------------------------
  #    Actually run the simulation
  # ------------------------------------------------------


  def run_simulation(self, params):
    '''
    This method needs to do everything necessary to *actually run* the simulation.
    This may include
    -- generating an input file, from the supplied parameter file
    -- generating a command needed to run the simulation
    -- calling a shell to actually run the command

    This method should run the simulation in a **sandboxed** manner, meaning that
    multiple instances of this command should be able to run from the same directory.
    Otherwise, capabilities designed for clusters will not work!!!

    Some BCA solvers expect input files with fixed names, and so simply running
    multiple instances of this command will result in the input files being overwritten.
    You will have to be creative in getting around this (i.e., creating a uniquely-named
    subdirectory in which to actually run the script.

    This method should also delete any files output by the program that will not be needed
    to extract statistics (i.e., anything that does not conceivably need to be saved.)
    '''
    pass




  # ----------------------------------------------------------
  #    Get an iterator that will return data for each impact
  # ----------------------------------------------------------



  def get_impact_iterator():

    '''
    This method needs to read the output of the program,
    and return an object ITER that, upon each call to ITER.next(),
    returns three objects:


    and extract the implantation data into a standard format:

          list[ N-list[ 3-array ] ]

    where N is the number of implanted particles, and for each particle we store, in order: 
    -- and the final position of the particle.

    and extract the erosive data into a standard format:

          list[ N-list[ int, 3-array ] ]

    where N is the number of sputtered particles, and for each particle we store, in order: 
    -- the ID number of the atom species
    -- and the initial position of the particle.


    and extract the redistributive data into a standard format:

          list[ N-list[ int, 3-array, 3-array ] ]

    where N is the number of sputtered particles, and for each particle we store , in order:
    -- the ID number of the atom species
    -- the initial position of the particle, 
    -- and the final position, in that order.
  
    '''

    return GENERIC_Impact_Iterator()



  # Extract and save the statistics (a generic function that does not need to be overridden)
  def extract_statistics(self, params):
    '''
    THIS FUNCTION SHOULD NOT NEED TO BE MODIFIED BY ANY DAUGHTER OBJECTS!
    '''

    '''
    This file uses other, package-specific routines to put the data into a standard format, 
    and then extracts some standard statistical data in a ***package-independent way.***
    This method therefore *should not be overridden.*  Doing so will likely break the 
    ability of other users to read your moments files.

    Note that if additional, package-dependent statistics are desired in the moments file,
    routines to collect this information can be appended to the user_statistics_routines
    list during the package-specific __init__() method.
    '''

    # initialization of storage space, and zero base values of all statistics
    globalavgs = dict()
    globalstds = dict()
    for stat in self.statistics_routines:
      stats = stat.collect(None, None, None, params)
      for key,value in stats.items():
        globalavgs[key] = value
        globalstds[key] = value

    # now get the statistics and incrementally update the avg and std
    II = self.get_impact_iterator(params)
    impact_data = II.next_impact()
    impact_counter = 1 
    while (impact_data != None):
      idata   = impact_data[0]
      edata   = impact_data[1]
      rdata   = impact_data[2]

      for stat in self.statistics_routines:
        istats = stat.collect(idata, edata, rdata, params) ;
        for key,value in istats.items():
          newavg, newstd = incremental_update(globalavgs[key], globalstds[key], value, impact_counter)
          globalavgs[key] = newavg
          globalstds[key] = newstd

      impact_counter += 1
      impact_data = II.next_impact()


    # now store the statistics in a shelf file
    mfile = shelve.open("%s.moms" % (params.fname()))
    mfile["impacts"] = impact_counter - 1
    for key,value in globalavgs.items():
      mfile["%s_avg"%(key)] = value
    for key,value in globalstds.items():
      mfile["%s_std"%(key)] = value
    mfile.close()


    return






  # Clean things up after everything is done.
  def clean_up(self, params):
    '''
    By this point everything the user wants stored has been stored in *new* files generated by
    the method extract_statistics().  So this file simply needs to delete everything that was 
    originally output by the program.
    '''
    pass

  


def incremental_update(oldmean, oldstd, newvalue, newcount):

  if (newcount==1):
    newmean = newvalue
    newstd  = np.zeros(np.shape(newvalue))
  else:
    n = float(newcount)
    oldvar  = oldstd**2
    newmean = ((n-1)*oldmean + newvalue) / n
    newvar  = (n-2)/(n-1) * oldvar + (newvalue - oldmean)**2 / n
    newstd  = np.sqrt(newvar)

  return newmean, newstd




