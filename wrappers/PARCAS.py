#
# this script analyzes a set of input/output files from MD simulation to generate
# average erosive and redistributive moments, as well as average displacement vector
# fields.
# 
# usage:  python process-md.py initial_pos.xyz impact_points.dat final_pos*.xyz
# argument  1  is the file containing the initial atomic positions
# argument  2  is the file containing the impact points and angles
# arguments 3  is a regular expression matching all files with final atomic positions; e.g. 0??.xyz 
#
# this script will save all of the averaged data into a file called 'averages.db',
# which can then be used with the averages from other angles to generate moments
# as a function of angle.
#
# Require modules: numpy, matplotlib
# 



import os
import sys
import glob
import shelve
import numpy as np
from subprocess import call

#import pycraters.wrappers_new.PARCAS as wrap
from pycraters.wrappers.GENERIC import *



# empty class for holding geometry
class Geom: 
    def __init__(self):
        self.atoms = 0 ;
        self.boxX  = 0 ;
        self.boxY  = 0 ;
        self.boxZ  = 0 ;

        self.xMin = 0 ;
        self.yMin = 0 ;
        self.zMin = 0 ;
        self.xMax = 0 ;
        self.yMax = 0 ;
        self.zMax = 0 ;

    def __eq__(self, otherG):
        val = True ;
        if (self.atoms != otherG.atoms): val = False ;
        if (abs(self.boxX - otherG.boxX) > .00001):  val = False ;
        if (abs(self.boxY - otherG.boxY) > .00001):  val = False ;
        if (abs(self.boxZ - otherG.boxZ) > .00001):  val = False ;
        return val ;

    def __ne__(self, otherG):
        if (self == otherG): return False ;
        return True ;





# empty class for holding bin structure
class BinData:
    def __init__(self, geom, binsize, rbinsize):

        # define some geometries
        xbins = int(round(geom.boxX / binsize[0])) ;
        ybins = int(round(geom.boxY / binsize[1])) ;
        zbins = int(round(geom.boxZ / binsize[2]) + 1) ;
        self.bincount = (xbins, ybins, zbins) ;

        xbinwidth = geom.boxX / self.bincount[0] ;
        ybinwidth = geom.boxY / self.bincount[1] ;
        zbinwidth = geom.boxZ / (self.bincount[2] - 1) ;
        self.binsize = (xbinwidth, ybinwidth, zbinwidth) ;
 
        self.rbins = round(min(geom.boxX, geom.boxY) / rbinsize) ;
        self.rbinwidth = min(geom.boxX, geom.boxY) / self.rbins ;

        # define a co-ordinate mesh corresponding to the bin centers
        mgrid = np.lib.index_tricks.nd_grid()

        RR = mgrid[0:self.rbins] ;
        RR = self.rbinwidth * (0.5 + RR) ;

        X1 = mgrid[0:xbins]
        Y1 = mgrid[0:ybins]
        Z1 = mgrid[0:zbins]
        X1 = geom.xMin + xbinwidth * (0.5 + X1) ;
        Y1 = geom.yMin + ybinwidth * (0.5 + Y1) ;
        Z1 = geom.zMin + zbinwidth * (0.5 + Z1) ;

        XX,YY = mgrid[0:xbins,0:ybins]
        XX = geom.xMin + xbinwidth * (0.5 + XX) ;
        YY = geom.yMin + ybinwidth * (0.5 + YY) ;

        xx,yy,zz = mgrid[0:xbins,0:ybins,0:zbins]
        xx = geom.xMin + xbinwidth * (0.5 + xx) ;
        yy = geom.yMin + ybinwidth * (0.5 + yy) ;
        zz = geom.zMin + zbinwidth * (0.5 + zz) ;

        self.x1 = X1 ;
        self.y1 = Y1 ;
        self.z1 = Z1 ;
        self.XX = XX ;
        self.YY = YY ;
        self.RR = RR ;
        self.xx = xx ;
        self.yy = yy ;
        self.zz = zz ;




# read geometry from data files
def read_geometry_data(fname):
    g = Geom() ;
    geodata = np.loadtxt(fname, usecols=(1,))

    #myfile   = open(fname, 'r') 
    #mytext   = myfile.readline() 
    #try: g.atoms = int(mytext) - 1 ;
    #except: print 'File', fname, 'improperly formatted.' ; exit() ;
    #mytext   = myfile.readline()
    #myarray  = mytext.split()
    g.atoms = int(geodata[0])
    g.xMin  = geodata[1] ;
    g.xMax  = geodata[2] ;
    g.yMin  = geodata[3] ;
    g.yMax  = geodata[4] ;
    g.zMin  = geodata[5] ;
    g.zMax  = geodata[6] ;

    g.boxX  = float(g.xMax-g.xMin)
    g.boxY  = float(g.yMax-g.yMin)
    g.boxZ  = float(g.zMax-g.zMin)

    return g



def my_ndgrid2(x1, x2):

    n1 = len(x1) ;
    n2 = len(x2) ;
    M1 = zeros((n1, n2)) ;
    M2 = zeros((n1, n2)) ;
    for ii in range(0,n1):
        for jj in range(0,n2):
            M1[ii,jj] = x1[ii] ;
            M2[ii,jj] = x2[jj] ;

    return M1,M2












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




class PARCAS_Parameters(GENERIC_Parameters):

  def __init__(self):
    # call the initialization of the parent class
    super(PARCAS_Parameters, self).__init__() 

    # now customize to fit the needs of PARCAS
    self.description = "PARCAS Parameters"

    # additional PARCAS parameters describing *impact* environment
    #self.additional_parameters[""] = ...  # PARCAS-specific parameter 1
    #self.additional_parameters[""] = ...  # PARCAS-specific parameter 1
    #self.additional_parameters[""] = ...  # PARCAS-specific parameter 1



class PARCAS_Impact_Iterator(GENERIC_Impact_Iterator):


  # ------------------------------------------------------
  # Set things up so that the simulation can be run
  # ------------------------------------------------------
  def __init__(self, parcas_wrapper=None, parcas_params=None):
    '''
    '''
    self.parcas_wrapper = parcas_wrapper
    self.parcas_params  = parcas_params

    print "Initializing: ",
    sys.stdout.flush()

    datadir = parcas_params.fname()
    geomfile = "%s/%s" % (datadir, parcas_wrapper.simulation_parameters["geomfile"])
    startfile = "%s/%s" % (datadir, parcas_wrapper.simulation_parameters["startfile"])
    impactfile = "%s/%s" % (datadir, parcas_wrapper.simulation_parameters["impactfile"])

    # information needed to access and read the ending files
    self.datadir = datadir
    os.chdir(datadir)
    self.endfiles = sorted(glob.glob(parcas_wrapper.simulation_parameters["endfile_exp"]))
    os.chdir("../")
    self.rowskips = parcas_wrapper.simulation_parameters["rowskips"]
    self.current_run = 0

    # read and compare geometry files; build bins and grids
    self.geodata  = np.loadtxt(geomfile, usecols=(1,))
    self.geom     = read_geometry_data(geomfile)
    self.bdata    = BinData(self.geom, [4,4,4], 2) ;
    self.cell_width = ((self.geom.xMax - self.geom.xMin) + (self.geom.yMax - self.geom.yMin)) / 2.0
    self.surf_height = self.geom.zMax 

    # load atomic species and generate unique list
    self.species_list = [aa[0] for aa in parcas_params.target]
    self.species_lookup = dict()
    species_count = 1
    for aa in self.species_list:
      self.species_lookup[aa] = species_count
      species_count += 1
    species_raw  = np.loadtxt(startfile, skiprows=self.rowskips, usecols=(0,), dtype='string')
    self.sp0 = [self.species_lookup[aa] for aa in species_raw[0:self.geom.atoms]]

    # load initial positions, correct for wrapping
    self.pos0    = np.loadtxt(startfile, skiprows=self.rowskips, usecols=[1,2,3]) 
    self.Ipos0   = self.pos0[-1,:]
    self.pos0    = self.pos0[0:self.geom.atoms,:]
    self.x0 = self.pos0[:,0] ;
    self.x0 += (self.x0 < self.geom.xMin)*self.geom.boxX ;
    self.x0 -= (self.x0 > self.geom.xMax)*self.geom.boxX ;
    self.y0 = self.pos0[:,1] ;
    self.y0 += (self.y0 < self.geom.yMin)*self.geom.boxY ;
    self.y0 -= (self.y0 > self.geom.yMax)*self.geom.boxY ;
    self.z0 = self.pos0[:,2] ;

    # load impact points 
    self.ipoints = np.loadtxt(impactfile, usecols=[1,2,3,4])
    print "done."
    sys.stdout.flush()

    # get background displacements
    self.dxP = None
    self.dyP = None
    self.dzP = None
    if (parcas_wrapper.simulation_parameters["remove_background"]):
      background_file = "%s/%s" % (self.datadir, self.parcas_wrapper.simulation_parameters["background"])
      if os.path.isfile(background_file):
        print 'Loading previously-gathered BG displacements.'
      else:
        print 'Collecting Background displacements. This may take some time.'
        self.process_background()

      f = shelve.open(background_file)
      self.dxP = f['dxP']
      self.dyP = f['dyP']
      self.dzP = f['dzP']
      f.close() ;

    return 



  #      process background (optional)
  # ------------------------------------
  def process_background(self):

    # get local instances of global variables
    geodata  = self.geodata
    geom     = self.geom
    bdata    = self.bdata
    ipoints  = self.ipoints
    rowskips = self.rowskips
    sp_count = len(self.species_list)

    # get local instances of loaded values
    x0 = self.x0
    y0 = self.y0
    z0 = self.z0
    sp0 = self.sp0

    # initialize some perm storage
    dxP  = np.zeros(geom.atoms) ;
    dyP  = np.zeros(geom.atoms) ;
    dzP  = np.zeros(geom.atoms) ;
    wtP  = np.zeros(geom.atoms) ;

    actP = np.zeros(np.append(sp_count+1, bdata.bincount)) ;
    vfxP = np.zeros(np.append(sp_count+1, bdata.bincount)) ;   
    vfyP = np.zeros(np.append(sp_count+1, bdata.bincount)) ;  
    vfzP = np.zeros(np.append(sp_count+1, bdata.bincount)) ; 

    # read final data to get background displacements
    for ff,endfile in enumerate(self.endfiles):

      realendfile = "%s/%s" % (self.datadir, endfile)
      print realendfile, ':  ',
      print 'reading ...',
      sys.stdout.flush()

      # read data from finalpos file  (CURRENTLY WASTEFUL!!!)
      ele2 = np.loadtxt(realendfile, skiprows=rowskips, usecols=[0], dtype='string')
      pos2 = np.loadtxt(realendfile, skiprows=rowskips, usecols=(1,2,3))
      pos2 = pos2[0:self.geom.atoms,:]

      print 'analyzing ...',
      sys.stdout.flush()

      # get impact info based on this filename
      parts  = endfile.split('.')
      pp = int(parts[0]) - 1;
      incidence = ipoints[pp,0] * np.pi / 180.0 ;
      azimuthal = ipoints[pp,1] * np.pi / 180.0 + np.pi ;
      xshift    = ipoints[pp,2] ;
      yshift    = ipoints[pp,3] ;

      # get shifted initial conditions and wrap
      x1 = x0 + xshift * geom.boxX ;
      x1 = x1 - geom.boxX * (x1 > geom.xMax)
      x1 = x1 + geom.boxX * (x1 < geom.xMin)
      y1 = y0 + yshift * geom.boxY ;
      y1 = y1 - geom.boxY * (y1 > geom.yMax)
      y1 = y1 + geom.boxY * (y1 < geom.yMin)
      z1 = z0

      # get final positions and displacements
      x2 = pos2[:,0] ;
      y2 = pos2[:,1] ;
      z2 = pos2[:,2] ;
      dx = x2 - x1 ;
      dy = y2 - y1 ;
      dz = z2 - z1 ;

      # try and correct for wrapping
      dx = dx - geom.boxX*(dx >  geom.boxX*.5) ;
      dx = dx + geom.boxX*(dx < -geom.boxX*.5) ;
      dy = dy - geom.boxY*(dy >  geom.boxY*.5) ;
      dy = dy + geom.boxY*(dy < -geom.boxY*.5) ;


      # identify atoms away from the impact point  
      # change to a sphere centered at the anticipated point of first collision
      #    x = ximpact - zmax/2*sin(theta)
      #    z = - zmax/2*cos(theta)
      dist        = np.sqrt(x1**2 + y1**2)    # + (z0 - geom.zMax)**2 ;
      away_impact = (dist > self.cell_width/4.0)

      # identify atoms away from boundary
      near_bdry   = np.zeros(geom.atoms) ;
      near_bdry  += (np.abs(x1) > geom.xMax - 10) ;
      near_bdry  += (np.abs(y1) > geom.yMax - 10) ;
      near_bdry  += (z1 < geom.zMin + 20) ;
      away_bdry   = (near_bdry == 0) ;

      # identify small displacements
      mag         = np.sqrt(dx**2 + dy**2 + dz**2)
      mag75       = np.percentile(mag, 75.0)
      small_disp  = (mag < mag75)

      # weight only those atoms far from the impact and not in cooling bdry
      weight     = away_impact * away_bdry * small_disp 

      # save weighted background effect
      wtP += weight ;
      dxP += (dx*weight) ;
      dyP += (dy*weight) ;
      dzP += (dz*weight) ;
      print 'done.', '\r'
      sys.stdout.flush()


    print 'Final analysis ...'
    # now analyze the data, atom by atom
    binx0 = np.floor((x0-geom.xMin) / bdata.binsize[0]) ;
    biny0 = np.floor((y0-geom.yMin) / bdata.binsize[1]) ;
    binz0 = np.floor((z0-geom.zMin) / bdata.binsize[2]) ;

    for aa in range(0, geom.atoms):

      if (wtP[aa] == 0): continue ;

      # divide weighted avg. displacement by total weight
      dxP[aa] /= wtP[aa]
      dyP[aa] /= wtP[aa]
      dzP[aa] /= wtP[aa]

      # add displacement to correct bin
      actP[sp0[aa], binx0[aa], biny0[aa], binz0[aa]] +=  1.0 ;
      vfxP[sp0[aa], binx0[aa], biny0[aa], binz0[aa]] +=  dxP[aa] ;
      vfyP[sp0[aa], binx0[aa], biny0[aa], binz0[aa]] +=  dyP[aa] ;
      vfzP[sp0[aa], binx0[aa], biny0[aa], binz0[aa]] +=  dzP[aa] ;

    print 'done.', '\r'
    sys.stdout.flush()


    # save the resulting data to the file
    print 'saving ...'
    background_file = "%s/%s" % (self.datadir, self.parcas_wrapper.simulation_parameters["background"])
    f = shelve.open(background_file)

    f['params'] = self.parcas_params
    f['geometry'] = [geom.atoms, geom.xMin, geom.xMax, geom.yMin, geom.yMax, geom.zMin, geom.zMax] ;

    f['x1'] = bdata.x1
    f['y1'] = bdata.y1
    f['z1'] = bdata.z1
    f['XX'] = bdata.XX
    f['YY'] = bdata.YY
    f['xx'] = bdata.xx
    f['yy'] = bdata.yy
    f['zz'] = bdata.zz

    f['dxP'] = dxP
    f['dyP'] = dyP
    f['dzP'] = dzP
    f['acX'] = actP
    f['vfx'] = vfxP
    f['vfy'] = vfyP
    f['vfz'] = vfzP

    f.close() ;
    print 'done.'







  #   process signal (next_impact)
  # ------------------------------------
  def next_impact(self):
    '''
    
    '''


    if (self.current_run >= len(self.endfiles)):
      return None


    # get local instances of global variables
    geodata  = self.geodata
    geom     = self.geom
    bdata    = self.bdata
    ipoints  = self.ipoints
    rowskips = self.rowskips
    sp_count = len(self.species_list)

    # get local instances of loaded values
    x0 = self.x0
    y0 = self.y0
    z0 = self.z0
    sp0 = self.sp0




    # --------------------------------------------
    # identify and read data from finalpos file
    # --------------------------------------------



    endfile = "%s/%s" % (self.datadir, self.endfiles[self.current_run])
    print '  reading final position file', endfile, '...',
    sys.stdout.flush()

    pos2 = np.loadtxt(endfile, skiprows=self.rowskips, usecols=[1,2,3])
    Ipos2 = pos2[-1,:]
    pos2 = pos2[0:self.geom.atoms,:]

    print 'done.'
    sys.stdout.flush()

    # --------------------------------------------
    #   coordinate transformation and wrapping
    # --------------------------------------------

    print '  coordinate transformation and wrapping  ...',
    sys.stdout.flush()

    # get impact info based on this filename
    parts  = self.endfiles[self.current_run].split('.')
    pp = int(parts[0]) - 1;
    incidence = ipoints[pp,0] * np.pi / 180 ;
    azimuthal = ipoints[pp,1] * np.pi / 180 + np.pi ;
    xshift    = ipoints[pp,2] ;
    yshift    = ipoints[pp,3] ;

    # get shifted initial conditions
    x1  = x0 + xshift * geom.boxX ;
    x1 -= geom.boxX * (x1 > geom.xMax)		# unnecessary?
    x1 += geom.boxX * (x1 < geom.xMin)		# unnecessary?
    y1  = y0 + yshift * geom.boxY ;
    y1 -= geom.boxY * (y1 > geom.yMax)		# unnecessary?
    y1 += geom.boxY * (y1 < geom.yMin)		# unnecessary?
    z1 =  z0

    # get rotated initial co-ordinates
    u1 =   x1 * np.cos(azimuthal) + y1 * np.sin(azimuthal) ;
    v1 = - x1 * np.sin(azimuthal) + y1 * np.cos(azimuthal) ;

    # get final co-ordinates
    x2  = pos2[:,0] ;
    x2 -= geom.boxX * (x2 > geom.xMax)		# unnecessary?
    x2 += geom.boxX * (x2 < geom.xMin)		# unnecessary?
    y2  = pos2[:,1] ;
    y2 -= geom.boxY * (y2 > geom.yMax)		# unnecessary?
    y2 += geom.boxY * (y2 < geom.yMin)		# unnecessary?
    z2  = pos2[:,2] ;

    # get displacements
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1

    # adjust for background
    if self.parcas_wrapper.simulation_parameters["remove_background"]:
      dx -= self.dxP
      dy -= self.dyP
      dz -= self.dzP

    # try to adjust for wrapping
    dx -= geom.boxX*(dx >  geom.boxX*.5) ;
    dx += geom.boxX*(dx < -geom.boxX*.5) ;
    dy -= geom.boxY*(dy >  geom.boxY*.5) ;
    dy += geom.boxY*(dy < -geom.boxY*.5) ;
    x2 = x1 + dx ;
    y2 = y1 + dy ;
    z2 = z1 + dz ;

    # get rotated final co-ordinates
    u2 =   x2 * np.cos(azimuthal) + y2 * np.sin(azimuthal) ;
    v2 = - x2 * np.sin(azimuthal) + y2 * np.cos(azimuthal) ;

    print 'done.'
    sys.stdout.flush()

    

    # --------------------------------------------
    #   identify implants, erosion, and redist.
    # --------------------------------------------

    print '  sorting  ...',
    sys.stdout.flush()

    Iu2 =   Ipos2[0] * np.cos(azimuthal) + Ipos2[1] * np.sin(azimuthal) ;
    Iv2 = - Ipos2[0] * np.sin(azimuthal) + Ipos2[1] * np.cos(azimuthal) ;
    Iz2 =   Ipos2[2]
    idata = [[  0, np.array([Iu2, Iv2, Iz2])  ]]

    eroded_atoms = np.where(z2 - self.surf_height >  5.0)
    elist = eroded_atoms[0]
    edata = [
            [ self.sp0[aa], np.array( [u1[aa], v1[aa], z1[aa]] ) ]
            for aa in elist
            ]

    redist_atoms = np.where(z2 - self.surf_height <= 5.0)
    rlist = redist_atoms[0]
    rdata = [ 
            [ self.sp0[aa], np.array([u1[aa], v1[aa], z1[aa]]), np.array([u2[aa], v2[aa], z2[aa]]) ] 
            for aa in rlist
            ]

    print 'done.'
    sys.stdout.flush()

    # --------------------------------------------
    #   increment the run counter and return
    # --------------------------------------------
    self.current_run += 1
    sys.stdout.flush()
    return [idata, edata, rdata]












class PARCAS_Wrapper(GENERIC_Wrapper):


  def __init__(self, exec_location):
    '''
    This needs to do everything needed to *allow* the simulation to be run.
    This may include 
    -- storing the location of the executable 
    -- storing the command needed to run the executable
    -- storing the location of any needed data files
    -- storing any instructions regarding capabilities that are specific to this solver
    '''

    # call the initialization of the parent class
    super(PARCAS_Wrapper, self).__init__() 


    # GLOBAL PARAMETERS COMMON TO ALL WRAPPERS
    # ------------------------------------------------------
    self.wrapper_type = "PARCAS MD Wrapper"
    self.execline = exec_location

    # In MD, we want to save raw data by default
    self.simulation_parameters["save_raw_data"] = True

    # whether to perform certain pre-analysis procedures
    self.simulation_parameters["remove_background"] = False
    self.simulation_parameters["correct_for_shear"] = False

    # set some parameters associated with directory structure
    self.simulation_parameters["geomfile"]    = "geometry.dat"
    self.simulation_parameters["impactfile"]  = "impact_points.dat"
    self.simulation_parameters["startfile"]   = "initial.xyz"
    self.simulation_parameters["endfile_exp"] = "???.xyz"
    self.simulation_parameters["rowskips"]    = 2

    # default location of the background file
    self.simulation_parameters["background"]  = "background.pkl"

    return



  # ------------------------------------------------------
  #    Actually run the simulation
  # ------------------------------------------------------


  def run_simulation(self, params):
    '''
    This section will need to be written by Andrey/Alvaro.

    The argument $params is here a user-supplied argument defined in
      PARCAS.py:  PARCAS_Parameters()
    which inherits some generic functionality from
      GENERIC.py:  GENERIC_Parameters()
    Essentially, it should contain all the information about the *simulation environment*
    -- anything that could change from run to run on a given user's machine.  The most
    basic elements are fields of params:  params.target, params.ion, params.energy, etc.
    and are defined in GENERIC.py:GENERIC_Parameters()  [you shouldn't need to change that file].
    Anything specific to PARCAS is contained in params.additional_parameters[""]; examples could
    include information on which pre-prepared target to use, which pre-calculated potential
    to use, and so on.

    There is another source of parameters that you might need, defined in 
      PARCAS.py:  PARCAS_Wrapper().simulation_parameters (a dictionary)
    which inherits some generic functionality from
      GENERIC.py:  GENERIC_Wrapper()
    This dictionary describes the *software environment* -- i.e., the name of the slurm queue,
    the location of the executable, the structure of the output directory, and so on.
    '''

    base_params = params
    extra_params = params.additional_parameters
    software_params = self.simulation_parameters

    '''
    This function will have to read the content of the three parameter objects defined above, and then

    (1)  Write any input files needed by PARCAS
    (2)  Run the commands to do simulation / submit jobs
    (3)  ensure that the result is in the directory structure assumed by PARCAS_Impact_Iterator()

    For instance, this could include

    -- creating a directory using the output of params.fname()
    -- copying in appropriate initial.xyz, geometry.dat, and impacts.dat
    -- creating an *external* script that runs the simulation and produces ???-final.xyz
    -- submitting n instances of that script to the slurm scheduler

    after slurm runs all of the jobs, the directory will be in the form expected by PARCAS_Impact_Iterator()
    '''


    pass




  # ----------------------------------------------------------
  #    Get an iterator that will return data for each impact
  # ----------------------------------------------------------



  def create_parameters(self):

    return PARCAS_Parameters()



  def get_impact_iterator(self, params):

    '''
    This method needs to read the output of the program,
    and return an object ITER that, upon each call to ITER.next(),
    returns three objects:  
    '''

    return PARCAS_Impact_Iterator(self, params)





  # Clean things up after everything is done.
  def clean_up(self, params):
    '''
    By this point everything the user wants stored has been stored in *new* files generated by
    the method extract_statistics().  So this file simply needs to delete everything that was 
    originally output by the program.
    '''
    pass

  








