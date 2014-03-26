import struct
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import gzip
from array import *
import os
from pcraster.numpy import *
from pcraster import *
import grads, datetime
import matplotlib.pyplot as plt

def datetime2gradstime(date):

 #Convert datetime to grads time
 str = date.strftime('%HZ%d%b%Y')

 return str

def gradstime2datetime(str):

 #Convert grads time to datetime
 date = datetime.datetime.strptime(str,'%HZ%d%b%Y')

 return date

def Check_and_Make_Directory(dir):

 #Check if a directory exists and if not create it
 if os.path.exists(dir) == False:
  os.system("mkdir %s" % dir)
  
 return

setclone("../maps/africa.map")

#Open grads and read in the VIC data
grads_exe = '/opt/opengrads/grads'
ga = grads.GrADS(Bin=grads_exe,Window=False,Echo=False)

datafile = "/home/latent2/mpan/vic/scripts/output_grid_19480101.ctl"
ga("open %s" % datafile)

#Get variable info
qh = ga.query("file")
vars = qh.vars
vars_info = qh.var_titles

#print vars 
#print vars_info

#Determine the last date
ga("set time 00z01Jan2005")
fdate = gradstime2datetime(ga.exp(vars[0]).grid.time[0])
ga("set time 00z01Jan2000")
idate = gradstime2datetime(ga.exp(vars[0]).grid.time[0])

print "Date range for conversion:"
print idate
print fdate

ga("set lat -89.875 89.625")
ga("set lon -179.875 179.625")

nt = (fdate - idate).days +1
date = idate
dt = datetime.timedelta(days=1)

print "Aggregate to daily data and save to map files:"

dir = "../VIC"
Check_and_Make_Directory(dir)

for t in xrange(0,nt):
  timestamp = datetime2gradstime(date+t*dt)
  timestamp2 = datetime2gradstime(date+t*dt + 7*datetime.timedelta(hours=3))
  dir = "../VIC/pcrglobwb_VIC_%s" % ((date+t*dt).strftime('%Y'))
  Check_and_Make_Directory(dir)

  print "Processing timestep %d, %s to %s" %(t,timestamp,timestamp2)
  ga("set time %s" % timestamp)
  
  #Read and resample precip
  ga("dprec = maskout(sum(prec, time=%s, time=%s),prec)" % (timestamp, timestamp2))
  ga("drprec = re(dprec,720,linear,-179.875,0.5,360,linear,-89.875,0.5,ba)")
  precdata = ga.exp("drprec")
  np.ma.set_fill_value(precdata,-9.999e+02)
  prec = np.true_divide(np.ma.getdata(precdata),1000)

  #Read and resample runoff
  ga("drunoff = maskout(sum(runoff, time=%s, time=%s),runoff)" % (timestamp, timestamp2))
  ga("drrunoff = re(drunoff,720,linear,-179.875,0.5,360,linear,-89.875,0.5,ba)")
  runoffdata = ga.exp("drrunoff")
  np.ma.set_fill_value(runoffdata,-9.999e+02)
  runoff = np.true_divide(np.ma.getdata(runoffdata),1000)

  #Read and resample baseflow
  ga("dbaseflow = maskout(sum(baseflow, time=%s, time=%s),baseflow)" % (timestamp, timestamp2))
  ga("drbaseflow = re(dbaseflow,720,linear,-179.875,0.5,360,linear,-89.875,0.5,ba)")
  baseflowdata = ga.exp("drbaseflow")
  np.ma.set_fill_value(baseflowdata,-9.999e+02)
  baseflow = np.true_divide(np.ma.getdata(baseflowdata),1000)

  #Produce daily map files
  precmap = numpy2pcr(Scalar,np.flipud(precdata),-9.999e+02)
  report(precmap,"%s/qw000000.%s" % (dir,(date+t*dt).strftime('%j')))
  runoffmap = numpy2pcr(Scalar,np.flipud(runoffdata+baseflowdata),-9.999e+02)
  report(runoffmap,"%s/qloc0000.%s" % (dir,(date+t*dt).strftime('%j')))

