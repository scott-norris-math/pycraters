import shelve
import sys
from subprocess import call

'''
This file serves a dual purpose: 

(1) it provides a very simple means of generating files that can be supplied
    to Condor for the purpose of running many simulations in batch mode.

(2) it provides a simple script called by the Condor scheduler that actually
    runs the simulations.

The code in this file is agnostic as to the simulation tool used.
'''


def run_condor_jobs(wrapper, params_list, email=None):

  # write the files
  write_condor_files(wrapper, params_list, email)

  # submit all jobs to the condor scheduler
  # -----------------------------------------------
  call("condor_submit condor.job" % ())






# This function writes a condor jobs script, as well as a shelf file
def write_condor_files(wrapper, params_list, email=None):



  # write a python script to be called by Condor
  # -----------------------------------------------
  fstring = \
"""
import shelve
import sys
import pycraters.wrappers.GENERIC

# get arguments from stdin, supplied by Condor
shelffile = sys.argv[1]
pindex = int(sys.argv[2])

# extract the wrapper and params from the shelf
f = shelve.open(shelffile)
wrapper = f['wrapper']
plist   = f['params_list']
params  = plist[pindex]
f.close()

# run the simulation
wrapper.go(params)
"""
  f = open('condor_exec.py', "w")
  f.write(fstring)
  f.close()


  # store relevant data in shelf
  # ------------------------------------------
  f = shelve.open('condor.shelf')
  f['wrapper'] = wrapper ;
  f['params_list'] = params_list
  f.close()


  # write condor jobfile
  # ------------------------------------------
  f = open('condor.job', "w")

  # general-purpose things
  f.write('universe     = vanilla\n') 
  f.write('getenv       = true\n') 
  if email:
    f.write('notification = always\n') 
    f.write('notify_user  = %s\n' % (email)) 
  else:
    f.write('notification = never\n')

  # log and err files
  f.write('log          = condor.log\n')
  f.write('error        = condor.err\n')

  # define the executable 
  f.write('executable   = %s\n\n' % (sys.executable)) 

  # add a list of jobs to be executed
  for kk,pp in enumerate(params_list):
    fprefix = pp.fname()
    f.write('log          = %s.clog\n' % (fprefix)) 
    f.write('error        = %s.err\n' % (fprefix)) 
    f.write('output       = %s.out\n' % (fprefix)) 
    f.write('arguments    = "condor_exec.py condor.shelf %d"\n' % (kk)) 
    f.write('queue\n\n') 

  # clean up and return
  f.close()



  return

