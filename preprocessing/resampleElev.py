import subprocess
import numpy as np
import os

import netCDF4 as nc
from pcraster.numpy import *
from pcraster import *

## Function to get mapAttributes of the cloneMap
def getMapAttribute(cloneMap):
  co = "mapattr -p "+ cloneMap
  cOut,err = subprocess.Popen(co, stdout=subprocess.PIPE,stderr=open('/dev/null'),shell=True).communicate()
  mapAttr = {'cellsize': float(cOut.split()[7]) ,\
    'rows'    : float(cOut.split()[3]) ,\
    'cols'    : float(cOut.split()[5]) ,\
    'xUL'     : float(cOut.split()[17]),\
    'yUL'     : float(cOut.split()[19])}
  return(mapAttr)

## Clonemap and NetCDF files for resampling
cloneMap = "DEM.map"
ncFile = "Elevation.nc"
outFolder = "../maps/25deg/"
levels = [0.0,0.01,0.05,0.10,0.20,0.30,0.40, 0.50,0.60,0.70,0.80,0.90,1.00]
mapsLevels = ["000","001","005","010","020","030","040","050","060","070","080","090","100"]

### Setclone and start loading data
setclone(cloneMap)
foo = subprocess.Popen("mkdir "+outFolder, stdout=subprocess.PIPE,stderr=open('/dev/null'),shell=True).communicate()

mapProps = getMapAttribute(cloneMap)

mask = scalar(cloneMap)
numpyMask = pcr2numpy(mask, -999)

f = nc.Dataset(ncFile)

lons = np.array(f.variables["Longitude"])
lats = np.array(f.variables["Latitude"])
elev = f.variables["Elev"]

### Create output arrays
outMaps = np.zeros((mapProps["rows"], mapProps["cols"], len(levels))) - 999.0

### Loop over global cells of mask to none-missing values
for x in range(int(mapProps["cols"])):
  xDown = mapProps["xUL"] + x * mapProps["cellsize"] - mapProps["cellsize"]/2.0
  xUp = mapProps["xUL"] + x * mapProps["cellsize"] + mapProps["cellsize"]/2.0
  for y in range(int(mapProps["rows"])):
    if numpyMask[y,x] != -999:
      yDown = mapProps["yUL"] - y * mapProps["cellsize"] - mapProps["cellsize"]/2.0
      yUp = mapProps["yUL"] - y * mapProps["cellsize"] + mapProps["cellsize"]/2.0
      print str(xDown) +" "+str(yUp)
      lonWindowSelect = np.nonzero((lons > xDown) & (lons < xUp))[0]
      latWindowSelect = np.nonzero((lats > yDown) & (lats < yUp))[0]
      vals = elev[latWindowSelect[0]:latWindowSelect[-1],lonWindowSelect[0]:lonWindowSelect[-1]]
      try:
	elevProfile = np.array(vals)
      except:
	elevProfile = np.zeros(len(lonWindowSelect)*len(latWindowSelect))
      for l in range(len(levels)):
	outMaps[y,x,l] = np.percentile(elevProfile, levels[l]*100)

### Write data      
dataZero = numpy2pcr(Scalar, outMaps[:,:,0], -999.0)
report(dataZero, outFolder+"elev0000."+mapsLevels[0])

for l in range(len(levels)):
  data = numpy2pcr(Scalar, outMaps[:,:,l], -999.0)
  report(data-dataZero, outFolder+"elev0000."+mapsLevels[l])

