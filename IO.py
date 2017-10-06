import os
import shelve


class MissingFileError(Exception):
  def __init__(self, filename):
    self.filename = filename
  def __str__(self):
    return "File not found: %s" % (filename) #repr(self.value)





def read_value(path, params, quantity):
  # Given a directory $path full of moment files, this function
  # takes a $params object, reads the associated data file, and returns
  # the associated $quantity.

    value = None

    # specify the target file
    targetfile = '%s%s.moms' % (path, params.fname())
    print "targetfile is %s" % (targetfile)
    if (os.path.isfile(targetfile) == False):
      print "File not found: ", targetfile, "\n\n"
      raise MissingFileError(targetfile)      

    # now try and extract the specified value
    try:
      f = shelve.open(targetfile)
      value = f[quantity]
      f.close()
    except KeyError as ee:
      print "Key Error"
      print "Working directory is %s." % (os.getcwd())
      print "Target file is %s" % (targetfile)
      print "Did not find key %s." % (ee)
      print "Available keys are", f.keys(), "\n\n"

    return value






def array_range(path, params, field, values, quantity):
  # Given a directory $path full of moment files, this function
  # takes a basic set of $params, sets $params.$field equal to each value
  # in the list $values, and retrieves the associated $quantity.  Returns
  # a list of the same size as $values.

  yvals = []

  for kk,val in enumerate(values):
    setattr(params, field, val)
    targetfile = '%s%s.moms' % (path, params.fname())
    if (os.path.isfile(targetfile) == False):
      print "File not found: ", targetfile, "\n\n"
      raise MissingFileError(targetfile)      
    try:
      f = shelve.open(targetfile)
      yvals.append(f[quantity])
      f.close()
    except KeyError as ee:
      print "Key Error"
      print "Working directory is %s." % (os.getcwd())
      print "Target file is %s" % (targetfile)
      print "Did not find key %s." % (ee)
      print "Available keys are", f.keys(), "\n\n"

  return yvals









class AtomicProperties(object):

  def __init__(self):
    pass


def extract_atomic_properties(atomstring, lines):

  #find the beam properties in table1
  linenum = 12
  while (True):
    entries = lines[linenum].split()
    species = entries[0] 

    if (species == atomstring):
      break
    elif (linenum < 121):
      linenum += 1
      continue
    else:
      print "\n\nFatal: Element %s not found.\n\n" % (atomstring)
      exit()


  ap = AtomicProperties()

  # Real Atomic Properties
  ap.atomic_number	= entries[1]
  ap.atomic_mass	= entries[2]
  ap.density_gcc	= entries[3]
  ap.density_aa		= entries[4]

  ## BCA Energy Parameters
  ap.surf_bind_energy	= entries[5]
  ap.displ_energy	= entries[6]
  ap.cutoff_energy	= entries[7]
  #ap.bulk_bind_energy	= entries[12]
  ap.bulk_bind_energy   = 0.0

  return ap






