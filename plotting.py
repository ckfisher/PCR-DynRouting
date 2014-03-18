import os
from osgeo import gdal, gdalconst, osr
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib import colors

def convertXY(xy_source, inproj, outproj):
    # function to convert coordinates

    shape = xy_source[0,:,:].shape
    size = xy_source[0,:,:].size

    # the ct object takes and returns pairs of x,y, not 2d grids
    # so the the grid needs to be reshaped (flattened) and back.
    ct = osr.CoordinateTransformation(inproj, outproj)
    xy_target = np.array(ct.TransformPoints(xy_source.reshape(2, size).T))

    xx = xy_target[:,0].reshape(shape)
    yy = xy_target[:,1].reshape(shape)

    return xx, yy

# Read the data and metadata
file = os.path.abspath("results/qc000000.181")
dataset = gdal.Open(file, gdalconst.GA_ReadOnly)
if dataset is None:
	print "Dataset not read\n"

indata = dataset.ReadAsArray()
band = dataset.GetRasterBand(1)
ndv = band.GetNoDataValue()
data = np.ma.masked_where(indata == ndv, indata)
gt = dataset.GetGeoTransform()
proj = dataset.GetProjection()

xres = gt[1]
yres = gt[5]

# Get the edge coordinates and shift to center coordinates
xmin = gt[0] + xres * 0.5
xmax = gt[0] + (xres * dataset.RasterXSize) - xres * 0.5
ymin = gt[3] + (yres * dataset.RasterYSize) + yres * 0.5
ymax = gt[3] - yres * 0.5

dataset = None

#Grid of xy coordinates
xy_source = np.mgrid[xmin:xmax+xres:xres, ymax+yres:ymin:yres]

fig = plt.figure(figsize=(12, 6))
m = Basemap(projection='cyl', lon_0=0, llcrnrlat=ymin, urcrnrlat=ymax, \
	    llcrnrlon=xmin, urcrnrlon=xmax, resolution='c')

inproj = osr.SpatialReference()
inproj.ImportFromProj4(m.proj4string)
outproj = osr.SpatialReference()
outproj.ImportFromProj4(m.proj4string)

xx, yy = convertXY(xy_source, inproj, outproj)

im1 = m.pcolormesh(xx, yy, data[:,:].T, cmap=plt.cm.jet_r,norm=colors.LogNorm())
plt.clim(0.001,100000)
cb = m.colorbar(im1,"right")
cb.set_label('Discharge (CMS)')
m.drawcoastlines()
plt.show()


