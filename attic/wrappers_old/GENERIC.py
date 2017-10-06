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
    self.xwidth = None
    self.ywidth = None
    self.zwidth = None
    self.xbins  = None
    self.ybins  = None
    self.zbins  = None









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

    # parameters describing the target energies
    self.SBE_mat = None	# a matrix of size nxn (where n is the number of species)
    self.SBE_vec = None	# a vector of size n (if matrix not supplied)
    self.BBE_vec = None	# a vector of size n
    self.ethresh = None	# a vector of size n
    self.ecutoff = None	# a vector of size n

    # parameters describing the target geometry
    self.k11     = 0.0	# curvature in x-direction
    self.k12	 = 0.0	# mixed second derivative
    self.k22     = 0.0	# curvature in y-direction


    # -----------------------------------

    # simulation domain geometry  (deprecated, for now)
    self.geometry = GeometryData() #adds a geometry object to help integrate with BCA Wrapper
    self.geometry.xwidth = 100.0
    self.geometry.ywidth = 100.0
    self.geometry.zwidth = 100.0
    self.geometry.xbins = 20
    self.geometry.ybins = 20
    self.geometry.zbins = 20

    # Deprecated (for now)
    self.depth  = None  	# should be a single float 
    self.dep_res = 1    	# should be an integer (depth resolution, default 1)
    self.cfunc  = None  	# should be a function of a float which returns a list of size (# of targets)
    self.funcname = None	# should be a string uniquely describing a cfunc




    #  Dictionary for Additional Parameters (solver-specific)
    # ------------------------------------------------------------
    self.additional_parameters = dict()




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





class GENERIC_Wrapper(object):
  '''
  This is essentially an interface definition.
  Each solver will extend this class.
  '''


  # Set things up so that the simulation can be run
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
    self.user_statistics_routines = list()

    pass



  def add_statistics_routine(self, routine):
    '''
    THIS FUNCTION SHOULD NOT NEED TO BE MODIFIED BY ANY DAUGHTER OBJECTS!
    '''
    self.user_statistics_routines.append(routine)


  def set_parameter(self, key, value):
    '''
    THIS FUNCTION SHOULD NOT NEED TO BE MODIFIED BY ANY DAUGHTER OBJECTS!
    '''
    if key in self.simulation_parameters.keys():
      self.simulation_parameters[key] = value
    else:
      raise SimulationParameterError("Parameter %s not found in BCA_Wrapper of type: %s" % (key, self.wrapper_type))
    return


    



  # do everything at once (does not need to be modified!!!)
  def go(self, params, save_raw_data = False):
    '''
    THIS FUNCTION SHOULD NOT NEED TO BE MODIFIED BY ANY DAUGHTER OBJECTS!
    '''

    if (self.check_for_data(params) == False):	# see if data already exists 

      print "  running simulations."
      self.run_simulation(params)		# actually runs the simulation
      print "  extracting statistics."
      self.extract_moments(params)		# extracts the moments from raw data files
      if (save_raw_data == False):
        self.clean_up(params)			# cleans up output, leaving only raw data
        #self.delete_raw_output(params)		# deletes the raw data

    return 


  def check_for_data(self, params):
    '''
    THIS FUNCTION SHOULD NOT NEED TO BE MODIFIED BY ANY DAUGHTER OBJECTS!
    '''

    targetfile = "%s.moms" % (params.fname())
    return os.path.isfile(targetfile)




  # Actually run the simulation
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



  # Extract erosive data from the output, into a standard format
  def get_erosive_data(self, params):
    '''
    This method needs to read the output of the program, 
    and extract the erosive data into a standard format:

          list[ N-list[ int, int, 3-array ] ]

    where N is the number of sputtered particles, and for each particle we store, in order: 
    -- the impact in which it occurred, 
    -- the ID number of the atom species
    -- and the initial position of the particle.

    '''
    pass



  # Extract redistributive data from the output, into a standard format
  def get_redistributive_data(self, params):
    '''
    This method needs to read the output of the program, 
    and extract the redistributive data into a standard format:

          list[ N-list[ int, int, 3-array, 3-array ] ]

    where N is the number of sputtered particles, and for each particle we store , in order:
    -- the impact in which it occurred, 
    -- the ID number of the atom species
    -- the initial position of the particle, 
    -- and the final position, in that order.
    '''
    pass



  # Clean things up after everything is done.
  def clean_up(self, params):
    '''
    By this point everything the user wants stored has been stored in *new* files generated by
    the method extract_moments().  So this file simply needs to delete everything that was 
    originally output by the program.
    '''
    pass

  

  # Extract and save the moments (a generic function that does not need to be overridden)
  def extract_moments(self, params):
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

    # preparation
    # ---------------------------------

    mfilename = '%s.moms' % (params.fname())
    mfile = shelve.open(mfilename)

    # get number of target species and impacts
    species = len(params.target)
    impacts = params.impacts

    mfile["species"] = species
    mfile["impacts"] = impacts

    # configure storage for pure or mutli-species results
    extra_species = 0
    if species > 1:  extra_species = 1
    total_species = species + extra_species

    # storage for implantation interstitial distribution
    iidM0   = np.zeros((impacts))
    iidM1   = np.zeros((impacts, 3))
    iidM2   = np.zeros((impacts, 3, 3))

    # storage for erosive vacancy distribution
    evdM0   = np.zeros((impacts, total_species))
    evdM1   = np.zeros((impacts, total_species, 3))
    evdM2   = np.zeros((impacts, total_species, 3, 3))

    # storage for redistributive vacancy distribution
    rvdM0   = np.zeros((impacts, total_species, 3))
    rvdM1   = np.zeros((impacts, total_species, 3))
    rvdM2   = np.zeros((impacts, total_species, 3, 3))
  
    # storage for redistributive interstitial distribution
    ridM0   = np.zeros((impacts, total_species, 3))
    ridM1   = np.zeros((impacts, total_species, 3))
    ridM2   = np.zeros((impacts, total_species, 3, 3))

    # storage for redistributive displacement distribution
    rddM0   = np.zeros((impacts, total_species, 3))
    rddM1   = np.zeros((impacts, total_species, 3))
    rddM2   = np.zeros((impacts, total_species, 3, 3))




    # get implantation moments
    # ----------------------------------
    idata = self.get_implantation_data(params)
    for entry in idata:

      # extract stored data
      iid  = entry[0]	# impact ID
      pos0 = entry[1]	# initial position

      # update moment values
      iidM0[iid] += 1
      iidM1[iid, :] += pos0
      iidM2[iid, :, :] += np.outer(pos0, pos0)

    # now calculate averages and std. dev.
    mfile["iidM0_avg"] = np.mean(iidM0, axis=0)
    mfile["iidM0_std"] = np.std (iidM0, axis=0)
    mfile["iidM1_avg"] = np.mean(iidM1, axis=0)
    mfile["iidM1_std"] = np.std (iidM1, axis=0)
    mfile["iidM2_avg"] = np.mean(iidM2, axis=0)
    mfile["iidM2_std"] = np.std (iidM2, axis=0)



    # get scattering moments
    # ------------------------------------
    #sdata = self.get_scattering_data(params)
    #for entry in sdata:

    #  # extract stored data
    #  iid  = entry[0]
    #  pos0 = entry[1]

    #  # update moment values
    #  isdM0[iid] += 1
    #  isdM1[iid, :] += pos0
    #  isdM2[iid, :, :] += np.outer(pos0, pos0)

    ## now calculate averages and std. dev.
    #mfile["isdM0_avg"] = np.mean(iidM0, axis=0)
    #mfile["isdM0_std"] = np.std (iidM0, axis=0)
    #mfile["isdM1_avg"] = np.mean(iidM1, axis=0)
    #mfile["isdM1_std"] = np.std (iidM1, axis=0)
    #mfile["isdM2_avg"] = np.mean(iidM2, axis=0)
    #mfile["isdM2_std"] = np.std (iidM2, axis=0)
    


    # get erosive moments
    # ----------------------------------
    edata = self.get_erosive_data(params)
    for entry in edata:

      # extract stored data
      iid  = entry[0]	# impact ID
      sid  = entry[1]+1	# species ID
      pos0 = entry[2]	# initial position of sputtered atom

      # update total moment values
      evdM0[iid,0] += 1
      evdM1[iid,0,:] += pos0
      evdM2[iid,0,:,:] += np.outer(pos0, pos0)

      # update component-based moment values
      if extra_species > 0:
        evdM0[iid,sid] += 1
        evdM1[iid,sid,:] += pos0
        evdM2[iid,sid,:,:] += np.outer(pos0, pos0)

    # now calculate averages and std. dev.
    mfile["evdM0_avg"] = np.mean(evdM0, axis=0)
    mfile["evdM0_std"] = np.std (evdM0, axis=0)
    mfile["evdM1_avg"] = np.mean(evdM1, axis=0)
    mfile["evdM1_std"] = np.std (evdM1, axis=0)
    mfile["evdM2_avg"] = np.mean(evdM2, axis=0)
    mfile["evdM2_std"] = np.std (evdM2, axis=0)

    # "old-style" statistics
    mfile["m0e_avg"] = -np.mean(evdM0, axis=0)
    mfile["m0e_std"] = np.std (evdM0, axis=0)
    mfile["m1e_avg"] = -np.mean(evdM1, axis=0)
    mfile["m1e_std"] = np.std (evdM1, axis=0)
    mfile["m2e_avg"] = -np.mean(evdM2, axis=0)
    mfile["m2e_std"] = np.std (evdM2, axis=0)




    # get redistributive moments
    # ----------------------------------
    rdata = self.get_redistributive_data(params)
    for entry in rdata:

      # extract stored data
      iid  = entry[0]	# impact ID
      sid  = entry[1]+1	# species ID
      pos0 = entry[2]	# initial position
      pos1 = entry[3]	# final position
      dpos = pos1-pos0	# displacement

      # update total moment values
      rvdM0[iid, 0] += 1
      rvdM1[iid, 0, :] += pos0
      rvdM2[iid, 0, :, :] += np.outer(pos0, pos0)

      ridM0[iid, 0] += 1
      ridM1[iid, 0, :] += pos1
      ridM2[iid, 0, :, :] += np.outer(pos1, pos1)

      rddM0[iid, 0] += 1
      rddM1[iid, 0, :] += dpos
      rddM2[iid, 0, :, :] += np.outer(dpos, dpos)

      # update component-based moment values
      if extra_species > 0:

        rvdM0[iid, sid] += 1
        rvdM1[iid, sid, :] += pos0
        rvdM2[iid, sid, :, :] += np.outer(pos0, pos0)

        ridM0[iid, sid] += 1
        ridM1[iid, sid, :] += pos1
        ridM2[iid, sid, :, :] += np.outer(pos1, pos1)

        rddM0[iid, sid] += 1
        rddM1[iid, sid, :] += dpos
        rddM2[iid, sid, :, :] += np.outer(dpos, dpos)


    # now calculate averages and std. dev.
    mfile["rvdM0_avg"] = np.mean(rvdM0, axis=0)
    mfile["rvdM0_std"] = np.std (rvdM0, axis=0)
    mfile["rvdM1_avg"] = np.mean(rvdM1, axis=0)
    mfile["rvdM1_std"] = np.std (rvdM1, axis=0)
    mfile["rvdM2_avg"] = np.mean(rvdM2, axis=0)
    mfile["rvdM2_std"] = np.std (rvdM2, axis=0)

    mfile["ridM0_avg"] = np.mean(ridM0, axis=0)
    mfile["ridM0_std"] = np.std (ridM0, axis=0)
    mfile["ridM1_avg"] = np.mean(ridM1, axis=0)
    mfile["ridM1_std"] = np.std (ridM1, axis=0)
    mfile["ridM2_avg"] = np.mean(ridM2, axis=0)
    mfile["ridM2_std"] = np.std (ridM2, axis=0)

    mfile["rddM0_avg"] = np.mean(rddM0, axis=0)
    mfile["rddM0_std"] = np.std (rddM0, axis=0)
    mfile["rddM1_avg"] = np.mean(rddM1, axis=0)
    mfile["rddM1_std"] = np.std (rddM1, axis=0)
    mfile["rddM2_avg"] = np.mean(rddM2, axis=0)
    mfile["rddM2_std"] = np.std (rddM2, axis=0)

    # "old-style" statistics
    mfile["m0r_avg"] = np.mean(ridM0 - rvdM0, axis=0)
    mfile["m0r_std"] = np.std (ridM0 - rvdM0, axis=0)
    mfile["m1r_avg"] = np.mean(ridM1 - rvdM1, axis=0)
    mfile["m1r_std"] = np.std (ridM1 - rvdM1, axis=0)
    mfile["m2r_avg"] = np.mean(ridM2 - rvdM2, axis=0)
    mfile["m2r_std"] = np.std (ridM2 - rvdM2, axis=0)




    # if the user has supplied external statistics routines, then call them here
    if self.user_statistics_routines != []:
      for routine in self.user_statistics_routines:
        user_stats = routine(params, edata, rdata)
        for key,value in user_stats.items():
          mfile[key] = value


    mfile.close()








