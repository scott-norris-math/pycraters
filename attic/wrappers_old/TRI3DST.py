#!/usr/bin/python
# coding: utf-8
'''

This file is the library that extends libCraters for the BCA solver tri3Dst.

The main points are the extension of the BCA_Wrapper class, a new parameter blob, and tri3Dst-specific IO functions.

'''

from pycraters.wrappers_old.GENERIC import *
import pycraters.IO as io
from subprocess import call
import shelve
import numpy as np
import sys



# In the TRI3DST co-ordinate system, the target consists of the half-space x > 0.
# Positive values of the incidence angle send the from the *positive* y-direction,
# with a projected ion direction in the negative y-direction.  This leads to the 
# following co-ordinate transform
# x (Craters) = -y (TRI3DST)
# y (Craters) = -z (TRI3DST)
# z (Craters) = -x (TRI3DST)
TRI3DST_coordinate_transform = np.matrix([[0, -1, 0],[0, 0, -1],[-1, 0, 0]])








#========================
#    Tri3DST Library
#========================

#shapedict = {"Flat Plane" : 1, "Sphere" : 2, "Horizontal Cylinder" : 3, "Cube" : 4, "3D Parabolic Contour" : 5}

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


# ---------------------------------



#========================
#    Tri3DST Wrapper
#========================

class TRI3DST_Wrapper(GENERIC_Wrapper):

  def __init__(self, execline):
    #set things up so the simulation can run

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



  def get_implantation_data(self, params):
    #Function to extract output files of the name run-id_SPLST.DAT

    filename = "%s-implant.out" % (params.fname())
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    idata = []

    for linenum,line in enumerate(lines):
      # as of 3/15/14, the initial position of sputtered atoms is located 
      # within the last three columns.  
      
      words = line.split()
      sid  = int(words[0])-2	# species ID
      iid  = int(words[1])-1 	# impact  ID
      pos0 = np.array( [ -float(words[-2]), -float(words[-1]), -float(words[-3]) ] )	# (transformed) initial co-ordinates
      idata.append([ iid, pos0 ]) 
  
    return idata



  def get_erosive_data(self, params):
    #Function to extract output files of the name run-id_SPLST.DAT

    filename = "%s-eros.out" % (params.fname())
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    edata = []

    for linenum,line in enumerate(lines):
      # as of 3/15/14, the initial position of sputtered atoms is located 
      # within the last three columns.  
      
      words = line.split()
      sid  = int(words[0])-2	# species ID
      iid  = int(words[1])-1 	# impact  ID
      pos0 = np.array( [ -float(words[-2]), -float(words[-1]), -float(words[-3]) ] )	# (transformed) initial co-ordinates
      edata.append([ iid, sid, pos0 ]) 
  
    return edata




  def get_redistributive_data(self, params):
    #Function to extract from output files of the name run-id_RCLST.DAT

    runid = params.runid()
    filename = "%s-redist.out" % (params.fname())
    f = open(filename, 'r')
    lines = f.readlines()
    f.close

    rdata = []

    for linenum,line in enumerate(lines):
      # as of 3/15/14, the initial position of redistributed atoms is located 
      # within the last three columns, and the final position is located in the
      # three columns before that
      
      words = line.split()
      sid = int(words[0])-2	# species ID
      iid = int(words[1])-1 	# impact ID
      pos0 = np.array( [ -float(words[-2]), -float(words[-1]), -float(words[-3]) ] )	# (transformed) initial position
      pos1 = np.array( [ -float(words[-5]), -float(words[-4]), -float(words[-6]) ] )	# (transformed) final position
      rdata.append([ iid, sid, pos0, pos1 ])

    return rdata



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
  f.write('%d %d %d -1 -1\n' % (geom.xbound, geom.ybound, geom.zbound)) #these are 400Å xyz geometries and open (not periodic) boundaries
  f.write('%d %d %d %.5f %.5f\n' % (geom.x, geom.y, geom.z, params.k11, params.k22) ) 
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







