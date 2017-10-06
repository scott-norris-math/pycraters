
import os
import sys
import glob
import shelve
import numpy as np
from subprocess import call

#import pycraters.wrappers_new.PARCAS as wrap
from pycraters.wrappers.GENERIC import *
import pycraters.IO as io






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














# ----------------------------------------------------
#       Generic Simulation Wrapper Definition
# ----------------------------------------------------




class TRI3DST_Parameters(GENERIC_Parameters):

  def __init__(self):

    # call the initialization of the parent class
    super(TRI3DST_Parameters, self).__init__() 

    # text string for diagnostic messages
    self.description = "TRI3DST environmental parameters"

    # default values for simulation parameters
    self.additional_parameters["shape"] = 5			# corresponds to parabolic contour
    self.additional_parameters["azimuthal_angle"] = 0.0
    self.additional_parameters["film_thickness"] = 0.0
    self.additional_parameters["film_composition_vector"] = None
    self.additional_parameters["SBV_matrix"] = None

    #self.geom = Geometry()
    #self.geom.set_attr([400, 400, 400], [800, 800, 800])



  def runid(self):
    return "0123456789" #TODO: update this





# =============================
#    Impact Iterator
# =============================
class TRI3DST_Impact_Iterator(GENERIC_Impact_Iterator):


  # set up -- open files and read first lines
  def __init__(self, params=None):

    # filenames
    self.ifname = "%s-implant.out" % (params.fname())
    self.efname = "%s-eros.out" % (params.fname())
    self.rfname = "%s-redist.out" % (params.fname())

    # open the files -- do *NOT* read all the lines
    self.ifile = open(self.ifname, 'r')
    self.efile = open(self.efname, 'r')
    self.rfile = open(self.rfname, 'r')

    # line position holders
    self.iline  = self.ifile.readline()
    self.eline  = self.efile.readline()
    self.rline  = self.rfile.readline()

    # impact counter
    self.icounter = 0



  #   next impact
  # -----------------
  def next_impact(self):

    if self.iline == self.eline == self.rline == "":  return None

    # get implantation data
    idata = []
    while self.iline != "":
      words = self.iline.split()
      iid  = int(words[1])-1	# impact ID
      sid  = int(words[0])-1	# species ID
      if (iid != self.icounter):  break ;
      pos0 = np.array( [ -float(words[-2]), -float(words[-1]), -float(words[-3]) ] )	# (transformed) initial co-ordinates as of 3/15/14
      idata.append([ sid, pos0 ]) 
      self.iline = self.ifile.readline()

    # get erosive data
    edata = []
    while self.eline != "":
      words = self.eline.split()
      iid  = int(words[1])-1 	# impact  ID
      sid  = int(words[0])-1	# species ID
      if (iid != self.icounter):  break ;
      pos0 = np.array( [ -float(words[-2]), -float(words[-1]), -float(words[-3]) ] )	# (transformed) initial co-ordinates as of 3/15/14
      edata.append([ sid, pos0 ]) 
      self.eline = self.efile.readline()

    # get redistributive data
    rdata = []
    while self.rline != "":
      words = self.rline.split()
      iid = int(words[1])-1 	# impact ID
      sid = int(words[0])-1	# species ID
      if (iid != self.icounter):  break ;
      pos0 = np.array( [ -float(words[-2]), -float(words[-1]), -float(words[-3]) ] )	# (transformed) initial position as of 3/15/14
      pos1 = np.array( [ -float(words[-5]), -float(words[-4]), -float(words[-6]) ] )	# (transformed) final position as of 3/15/14
      rdata.append([ sid, pos0, pos1 ])
      self.rline = self.rfile.readline()


    #   increment the run counter and return
    self.icounter += 1
    return [idata, edata, rdata]










#========================
#    Tri3DST Wrapper
#========================

class TRI3DST_Wrapper(GENERIC_Wrapper):

  #set things up so the simulation can run
  def __init__(self, execline):

    # call the initialization of the parent class
    super(TRI3DST_Wrapper, self).__init__() 

    self.wrapper_type = "TRI3DST wrapper"
    self.execline = execline



  def run_simulation(self, params):

    # get the base path from the execline
    path = ''
    dtree = self.execline.split('/')
    for dd in dtree[:-1]: path += (dd + '/')

    # change the current working directory
    subfolder = params.fname()
    call("mkdir %s" % (subfolder), shell=True)
    os.chdir( './%s' % (subfolder))

    # write the input file
    tri3dst_write_input(path, params, self.simulation_parameters)

    # frequent strings
    fname = params.fname()
    runid = params.runid()

    # run the command
    callstr = "%s < %s.in > %s.log" %(self.execline, fname, fname)
    call(callstr, shell=True)

    # clean up the results
    call("mv %s_rglst.dat %s-implant.out" % (runid, fname), shell=True)
    call("mv %s_splst.dat %s-eros.out"    % (runid, fname), shell=True)
    call("mv %s_rclst.dat %s-redist.out"  % (runid, fname), shell=True)
    call("rm *.dat", shell=True)

    # move key results back into parent, and delete the subdirectory
    call("mv *.in *.out *.log ../", shell = True)
    os.chdir('../')
    call("rmdir %s" % (subfolder), shell=True)
    return



  # Return an Impact Iterator
  def get_impact_iterator(self, params):
    return TRI3DST_Impact_Iterator(params)


  def create_parameters(self):
    return TRI3DST_Parameters()



  # Clean things up after everything is done.
  def clean_up(self, params):
    '''
    By this point everything the user wants stored has been stored in *new* files generated by
    the method extract_moments().  So this file simply needs to delete everything that was 
    originally output by the program.
    '''
    call("rm %s-*.out" % (params.fname()), shell=True)
  


















# -------------------------------------



def tri3dst_write_input (path, environment_parameters, tri3dst_parameters): #meant to take a parameter object

  #if (isinstance(environment_parameters.shape, basestring) == True):
  #  parameters.shape = shapedict[parameters.shape]
  params = environment_parameters
  extras = params.additional_parameters

  # build target based on ion energy
  geom   = Geometry()
  energy = params.energy
  geom.set_attr([energy, energy, energy], [2*energy, 2*energy, 2*energy])


  fname = params.fname() + '.in'
  f = open(fname, 'w')

  # THIS CODE IS FOR THE VERSION WOLFHARD E-MAILED TO ME ON 2014-07-25
  f.write("%s automatically generated file\n" % (params.runid()) )
  f.write('%d %d %d %.2f 12345\n' % (params.impacts, len(params.target)+1, extras["shape"], extras["film_thickness"]) )
  f.write('1 3 1 0 0. 0 0 0 3\n') #ask professor norris about these x-z atomic density averages
  f.write('%d %d %d -1 -1\n' % (geom.xbound, geom.ybound, geom.zbound)) #these are 400A xyz geometries and open (not periodic) boundaries
  f.write('%d %d %d %.8f %.8f\n' % (geom.x, geom.y, geom.z, params.k11, params.k22) ) 
  f.write('100 1 1 1 1 1 0 0 0 0 0 0  0 0 0\n') #output control. 
  f.write('-1 -1 1 5 0 0 0 0\n') #listed as diagnostic parameters


  # read standard values of various quantities for many elements
  table_location = path + 'table1'
  g = open(table_location, 'r')
  elements = g.readlines()
  g.close()


  # text describing the beam
  ap = io.extract_atomic_properties(params.beam, elements)
  if (params.ethresh != None):
    ap.displ_energy = params.ethresh[0]

  f.write("%s %s" % (ap.atomic_number, ap.atomic_mass))
  f.write(" %s %s %s" % (ap.bulk_bind_energy, ap.displ_energy, ap.cutoff_energy)) 
  f.write(" 0. 0. %s" % (ap.density_aa)) #fractional composition of bulk material and thin film, respectively (0 for beam) and atomic density.
  f.write(" 1.0\n") #electronic stopping correction factor

  f.write("%d. %d. %d." % (params.energy, params.angle, extras["azimuthal_angle"]) )
  f.write(" 1.00") 				# beam fraction

  # specify the surface binding energies; default case or with SBV provided
  sbe_row = ""
  if (extras["SBV_matrix"] == None):
    sbe_row = " %s" % (ap.surf_bind_energy)	# rows in surface binding energy matrix
    for kk in range(len(params.target)):
      sbe_row += " 0.00"
  else:
    my_row = extras["SBV_matrix"][0]
    for entry in my_row:
      sbe_row += " %8f" % (entry)
  f.write(sbe_row + "\n")


  f.write("-5 0 %d %d 25 0. 0. 0.\n" % (geom.y, geom.z) )



  for kk,element in enumerate(params.target):

    ap = io.extract_atomic_properties(element[0], elements)
    if (params.ethresh != None):
      ap.displ_energy = params.ethresh[kk+1]

    f.write("%s %s" % (ap.atomic_number, ap.atomic_mass))
    #f.write(" %s %s %s" % (ap.bulk_bind_energy, ap.displ_energy, ap.cutoff_energy)) 
    f.write(" %s %s %s" % (0.0, 2.0, 1.0)) 


    bulk_concentration = element[1]
    film_concentration = element[1] if extras["film_composition_vector"] == None else extras["film_composition_vector"][kk+1]
    f.write(" %s %s %s" % (bulk_concentration, film_concentration, ap.density_aa))   #fractional composition of bulk material and thin film, respectively (0 for beam) and atomic density.
    f.write(" 1.0\n") #electronic stopping correction factor

    f.write("%d. %d. %d." % (params.energy, params.angle, extras["azimuthal_angle"]) )
    f.write(" 0.00") 	# beam fraction

    # rows in surface binding energy matrix
    sbe_row = ""
    if (extras["SBV_matrix"] == None):
      sbe_row = " 0.0"
      for ii in range(len(params.target)):
        sbe_row += " %s" % (ap.surf_bind_energy)
    else:
      my_row = extras["SBV_matrix"][kk+1]
      for entry in my_row:
        sbe_row += " %8f" % (entry)
    f.write(sbe_row + "\n")


    f.write("-5 0 %d %d 25 0. 0. 0.\n" % (geom.y, geom.z) )

  f.close() 

























