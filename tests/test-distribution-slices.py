# basic system and numerical packages
import sys
import numpy as np
import shelve
import matplotlib.pyplot as plt
import pycraters as pc

# build wrapper and params
exec_location = sys.argv[1]
wrapper = pc.create_wrapper("TRI3DST", exec_location)
params  = wrapper.create_parameters()

# add a custom statistic
import pycraters.statistics.distribution_projections as dps
dps_stat = dps.distribution_projections([400.0,399.4,205.6],[80,80,40])
wrapper.add_statistics_routine(dps_stat)

#basic parameter setup
params.beam = "Ar"
params.energy = 1000
params.angle = 60
params.impacts = 1000
params.target = [["Si", 1.0]]

# run the simulation
wrapper.go(params)

# extract the statistics
f = shelve.open( "%s.moms" % (params.fname()) )

XZx = f['XZx_avg']
XZz = f['XZz_avg']
rddXZx = f['rddXZx_avg'][0]
rddXZz = f['rddXZz_avg'][0]

mag = (rddXZx**2 + rddXZz**2)**(0.01)

plt.figure(1, figsize=(12,6))
plt.quiver(XZx/10.0, XZz/10.0, rddXZx/mag, rddXZz/mag, scale=50.0)
plt.plot([0],[0],'bo',markersize=10)
plt.plot([-20, 20], [0,0], 'k--')
plt.xlim(-20,20)
plt.ylim(-12,12)
plt.show()
