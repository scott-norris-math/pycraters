
import os
import sys
import glob
import shelve
import numpy as np
from subprocess import call

#import pycraters.wrappers_new.PARCAS as wrap
from pycraters.wrappers.GENERIC import *
import pycraters.IO as io







class SDTRIMSP_Parameters(GENERIC_Parameters):

  def __init__(self):

    # call the initialization of the parent class
    super(SDTRIMSP_Parameters, self).__init__() 

    # text string for diagnostic messages
    self.description = "SDTRIMSP environmental parameters"

    # default values for simulation parameters
    self.additional_parameters["cfunc"] = None







# =============================
#    Impact Iterator
# =============================
class SDTRIMSP_Impact_Iterator(GENERIC_Impact_Iterator):


  # set up -- open files and read first lines
  # ---------------------------------------------
  def __init__(self, params=None):

    # filenames
    self.ifname = "%s-implant.out" % (params.fname())
    self.efname = "%s-eros.out" % (params.fname())
    self.rfname = "%s-redist.out" % (params.fname())

    lastline = " ende\n"

    # open the files
    self.ifile = open(self.ifname, 'r') ; 
    for _ in xrange(4):  self.ifile.readline()	# skip first 4 lines    
    self.efile = open(self.efname, 'r') ; 
    for _ in xrange(4):  self.efile.readline()	# skip first 4 lines
    self.rfile = open(self.rfname, 'r') ; 
    for _ in xrange(4):  self.rfile.readline()	# skip first 4 lines

    # line position holders
    self.iline  = self.ifile.readline() ; 
    if self.iline != lastline:  self.iline  = self.ifile.readline() ; # deal with SDTRIMSP inconsistency in 5th line
    self.eline  = self.efile.readline() ; 
    if self.eline != lastline:  self.eline  = self.efile.readline() ; # deal with SDTRIMSP inconsistency in 5th line
    self.rline  = self.rfile.readline() ; 
    if self.rline != lastline:  self.rline  = self.rfile.readline() ; # deal with SDTRIMSP inconsistency in 5th line

    # impact counter
    self.icounter = 0




  #   next impact
  # -----------------
  def next_impact(self):

    lastline = " ende\n"
    if self.iline == self.eline == self.rline == lastline:  return None


    # storage for implantation data
    # --------------------------------
    idata = []
    while self.iline != lastline:
      words = self.iline.split()
      sid = int(words[0])-1  		#species ID
      iid = round(float(words[2]))-1 	#impact ID
      if (iid != self.icounter):  break ;
      pos1 = np.array( [ float(words[8]), float(words[9]), -float(words[7]) ] )	# (transformed) initial co-ordinates
      idata.append([ sid, pos1 ]) 
      self.iline = self.ifile.readline()


    # storage for erosive data
    # --------------------------------
    edata = []
    while self.eline != lastline:
      words = self.eline.split()
      sid = int(words[0])-1  		#species ID
      iid = round(float(words[2]))-1 	#impact ID
      if (iid != self.icounter):  break ;
      pos0 = np.array( [ float(words[5]), float(words[6]), -float(words[4]) ] )	# (transformed) initial co-ordinates
      edata.append([ sid, pos0 ]) 
      self.eline = self.efile.readline()

    
    # storage for redistributive data
    # ---------------------------------
    rdata = []
    while self.rline != lastline:
      words = self.rline.split()
      sid = int(words[0])-1		# species ID
      iid = round(float(words[2]))-1	# impact ID
      if (iid != self.icounter):  break ;
      pos0 = np.array( [ float(words[5]), float(words[6]), -float(words[4]) ] )	# (transformed) initial co-ordinates
      pos1 = np.array( [ float(words[8]), float(words[9]), -float(words[7]) ] )	# (transformed) final co-ordinates
      rdata.append([ sid, pos0, pos1 ]) # store transformed co-ordiantes
      self.rline = self.rfile.readline()
 

    #   increment the run counter and return
    self.icounter += 1
    return [idata, edata, rdata]













#========================
#    SDTrimSP Wrapper
#========================

class SDTRIMSP_Wrapper(GENERIC_Wrapper):


  # Set things up so that the simulation can be run
  def __init__(self, execline):
    '''
    This needs to do everything needed to *allow* the simulation to be run.
    This may include 
    -- storing the location of the executable 
    -- storing the command needed to run the executable
    -- storing the location of any needed data files
    -- storing any instructions regarding capabilities that are specific to this solver
    '''

    # call the initialization of the parent class
    super(SDTRIMSP_Wrapper, self).__init__() 

    self.wrapper_type = "SDTrimSP wrapper"
    self.execline = execline

    # default values for simulation parameters
    self.simulation_parameters["ipot"] = 1		# corresponds to the Kr-C potential
    self.simulation_parameters["iintegral"] = 0		# corresponds to the MAGIC integration method



  # Actually run the simulation
  def run_simulation(self, params):
    '''
    This method needs to do everything necessary to *actually run* the simulation.
    This may include
    -- generating an in put file, from the supplied parameter file
    -- generating a command needed to run the simulation
    -- calling a shell to actually run the command
    '''

    # get the base path from the execline
    path = ''
    dtree = self.execline.split('/')
    for dd in dtree[:-3]: path += (dd + '/')

    # change the current working directory
    subfolder = params.fname()
    call("mkdir %s" % (subfolder), shell=True)
    os.chdir( './%s' % (subfolder))

    # write the input file(s)
    sdtrimsp_write_tri_inp(path, params, self.simulation_parameters)
  
    # run the trimsp command
    call("%s > %s.log" % (self.execline, params.fname()), shell=True)
  
    # clean up the results: rename only the files we want, delete the rest
    call("mv partic_stop_p.dat %s-implant.out" %(params.fname()), shell=True)
    call("mv partic_back_r.dat %s-eros.out"   %(params.fname()), shell=True)
    call("mv partic_stop_r.dat %s-redist.out" %(params.fname()), shell=True)
    call("mv partic_tran_r.dat %s-trans.out" %(params.fname()), shell=True)
    call("mv tri.inp %s-tri.inp" %(params.fname()), shell=True)
    if os.path.isfile('layer.inp'):  call("mv layer.inp %s-layer.inp" %(params.fname()), shell=True)
    call("rm *.dat", shell=True)

    # move key results back into parent, and delete the subdirectory
    call("mv *.out *.inp *.log ../", shell = True)
    os.chdir('../')
    call("rmdir %s" % (subfolder), shell=True)

    return
 


  def get_impact_iterator(self, params):

    return SDTRIMSP_Impact_Iterator(params)




  # Clean things up after everything is done.
  def clean_up(self, params):
    '''
    By this point everything the user wants stored has been stored in *new* files generated by
    the method extract_moments().  So this file simply needs to delete everything that was 
    originally output by the program.
    '''
    call("rm *.out", shell=True)












# ----------------------------------------------------
#       Generic Parameter Object Definition
# ----------------------------------------------------


class Geometry(object):
  def __init__(self):
    self.y = 400
    self.z = 400
    self.x = 400
    self.ybound = 800
    self.zbound = 800
    self.xbound = 800

  def set_attr(self, mid, bound):
    self.y = mid[0]
    self.z = mid[1]
    self.x = mid[2]
    self.ybound = bound[0]
    self.zbound = bound[1]
    self.xbound = bound[2]






class EnvironmentParameterError(Exception):
  def __init__(self, string):
    self.string = string
  def __str__(self):
    return self.string
















# ============================================
#    SDTrimSP-specific Run Functions
# ============================================



def sdtrimsp_write_tri_inp(path, params, sdtrim_parameters):
  '''
  Writes the files tri.inp for a given parameter blob.
  '''
  #os.chdir(subfolder)
  f = open('tri.inp', "w")
  
  f.write("auto-generated file\n")
  f.write("&TRI_INP\n\n")


  # specifying locations
  # -------------------------------------------------
  f.write("tableinp = '%stables'\n" % (path))


  # text for describing integration methods
  # -------------------------------------------------
  f.write("text='---basic operation mode parameters---'\n")
  f.write('idrel = 1\n')						# simulation mode (0 - dynamic, 1 - static)
  f.write('ipot      = %s\n' % sdtrim_parameters["ipot"])		# interaction potential
  f.write('iintegral = %s\n' % sdtrim_parameters["iintegral"])		# integration method
  f.write('\n\n')



  # text for specifying elements involved in the simulations
  # -------------------------------------------------
  f.write("text='---elements used---'\n")
  f.write("ncp = %d\n" % (len(params.target)+1))

  stext = 'symbol = "%s"' % (params.beam)
  for ts in params.target:
    stext += ', "%s"' % (ts[0])
  stext += '\n'
  f.write(stext)
  f.write('\n\n')


  # text for describing the beam and its properties
  # -------------------------------------------------
  f.write("text='---beam properties---'\n")

  f.write('case_e0 = 0\n')
  etext = 'e0 = %f' % (params.energy)
  for ts in params.target:
    etext += ', 0.00'
  etext += '\n'
  f.write(etext)

  f.write('case_alpha=0\n')
  atext = 'alpha0 = %f' % (params.angle)
  for ts in params.target:
    atext += ', 0.00'
  atext += '\n'
  f.write(atext)

  ctext = 'qubeam = 1.00'
  for ts in params.target:
    ctext += ', 0.00'
  ctext += '\n'
  f.write(ctext)
  f.write('\n\n')


  # text for describing the target thickness
  # -------------------------------------------------
  f.write("text='---target properties---'\n")
  targ_string = 'ttarget =  %f \n' % (params.energy)  # ASSUME that one angstrom per eV is thick enough
  f.write(targ_string)
  depth_string = 'nqx = %d \n' % (params.energy / 10.0) # one layer per nanometer (recommended minimum)
  f.write(depth_string)	         		



  # text for describing target composition
  # -------------------------------------------------
  cfunc = params.additional_parameters["cfunc"]
  if cfunc != None:
    # if a cfunc is provided, then also write a layer.inp function
    f.write('iq0 = -1\n')
    sdtrimsp_write_layer_inp(params)

  else:
    f.write('iq0 = 0\n') 
    ctext = 'qu = 0.00'				# list of concentrations (zero for beam species)
    for ts in params.target:
      ctext += ', %f' % (ts[1])
    ctext += '\n'
    f.write(ctext)

  ctext = 'qumax = 0.00'			# list of max concentrations (only in dynamic mode)
  for ts in params.target:
    ctext += ', 1.00'
  ctext += '\n'
  f.write(ctext)
  f.write('\n\n')


  # text for describing program control
  f.write("text='---impacts and control options---'\n")
  f.write('flc = %d\n' % (params.impacts))	# 'fluence' (set to nh in static mode, possibly means something in dynamic)
  f.write('nh = %d\n' % (params.impacts))	# number of 'histories' (i.e., impacts)
  f.write('nr_pproj= 1\n')			# applies mostly in dynamic mode, also possibly something about parallel (???)
  f.write('idout=1\n')				# how many impacts between outputs


  # output all displacements (no displacement energy)

  f.write("text='---BCA energies---'\n")
  f.write('isbv  = 1\n')			# flag to describe surface binding model
  f.write('sfin  = 0\n')			# description of inelastic energy losses (0 - no loss)

  edispl = params.ethresh
  if edispl == None:
    f.write('irc0  = -1\n')
    f.write('lpart_r_ed = .false.\n')
    etext = 'e_displ = 0.0'
    for ts in params.target:
      etext += ', 0.00'
    etext += '\n'
    f.write(etext)

  else:
    f.write('irc0  = 1\n')
    f.write('lpart_r_ed = .true.\n')
    etext = 'e_displ = 0.0'
    for ee in edispl[1:]:
      etext += ', %f' % ee
    etext += '\n'
    f.write(etext)





  # text for describing output
  f.write("text='---output parameter---'\n")
  f.write('lparticle_r= .true.\n')
  f.write('lparticle_p= .true.\n')
  f.write('ltraj_r= .false.\n')
  f.write('ltraj_p= .false.\n')
  f.write('numb_hist= 100000\n')
  f.write('ioutput_hist = 10000000, 10000000, 10000000, 10000000, 10000000, 10000000\n')
  f.write('ioutput_part= 10000000, 10000000, 10000000, 10000000, 10000000, 10000000\n')
  f.write('lmatrices = .true.\n')
  f.write('/\n\n\n\n\n')

  f.close()



def sdtrimsp_write_layer_inp(params):
  # profile_func needs to return an iterable list of concentrations
  #os.chdir(subfolder)
  print "Writing layer.inp"
  
  f = open("layer.inp", "w")

  f.write("Below me stand:\n") #headers
  f.write("#Layers, Thickness(Angstrom), Composition (mole fraction)\n")
  
  step_depth = 10.0					# We assume layers are 1nm thick (recommended minimum)
  current_depth = 5.0 					# We start at the midpoint of the top layer
  max_depth = params.energy				# We ASSUME that one angstron per eV is deep enough
  cfunc = params.additional_parameters["cfunc"]

  while (current_depth + step_depth) < max_depth: 	#as long as the next step isn't too big
    param_string = '1 '+ str(step_depth) 		#make a string for 1 layer, X angstroms thick
    for a in cfunc(params, current_depth):
      param_string += ' %1.8f' % (a) 			#the relative concentrations of elements.
    param_string += '\n'
    f.write(param_string) 
    current_depth += step_depth

  #create last layer, if necessary. 
  if (max_depth != current_depth):
    param_string = '1 ' + str(max_depth - current_depth) #make a string for last layer
    for a in cfunc(params, current_depth):
      param_string += ' ' + str(a) 			#the relative concentrations of elements.
    param_string += '\n'
    f.write(param_string)

  f.write("0 0 0 0\n")
  f.close()


