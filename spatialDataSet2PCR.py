import os, sys, tempfile, shutil
import subprocess
import pcraster as pcr

class getSpatialAttributes:
	#-retrieves attributes of a spatial data set

	def __init__(self,fileName):
		#-character strings
		whiteSpace= ['\t',' ','']
		coordinatePairs= ['Lower Left','Upper Right']
		#-get resolution and type
		co= 'gdalinfo %s' % fileName
		cOut,err = subprocess.Popen(co, stdout= subprocess.PIPE,\
			stderr= subprocess.PIPE,shell=True).communicate()  
		if err != '':
				sys.exit('information on type and resolution for %s could not be retrieved' % fileName)
		rawList= cOut.split('\n')

		for iCnt in xrange(len(rawList)-1):
			entry= rawList[iCnt]
			try:
				entry= entry.strip('\r')
			except:
				pass
			for coordinatePair in coordinatePairs:
				if coordinatePair in entry:
					txtList= []
					for separator in ['(',',',')']:
						charPos= entry.find(separator)
						txtStr= entry[:charPos].strip()
						txtStr= txtStr.replace(' ','_').lower()
						txtList.append(txtStr)
						entry= entry[charPos+1:]
					#set attribute
					setattr(self,txtList[0],[float(txtList[1]),float(txtList[2])])

##    #-get resolution and type
##    co= 'gdalinfo %s' % fileName
##    cOut,err = subprocess.Popen(co, stdout=subprocess.PIPE,stderr=open('/dev/null'),shell=True).communicate()
##    if err != None:
##        sys.exit('information on type and resolution could not be retrieved')
##    #print '* %s *' % cOut
##    infoList= cOut.replace('\n','\t').split()
##    iCnt= infoList.index('Pixel')
##    iCnt+= infoList[iCnt:].index('Size')
##    resolutions= infoList[iCnt+2].strip(' ()').split(',')
##    dx_o= abs(float(resolutions[0]))
##    dy_o= abs(float(resolutions[1]))
##    typeStr= 'None'
##    for iCnt in xrange(len(infoList)):
##      if 'Type' in infoList[iCnt]:
##        typeStr= infoList[iCnt]
##        typeStr= typeStr.replace('Type','')
##        typeStr= typeStr.strip(',=')

class spatialDataSet:
	#-handles the processing of spatial datasets in GDAL, stores data as numpy array in memory under the variable name specified

	def __init__(self,fileName,variable,typeStr,outputFormat,\
				x0,y0,x1,y1,dx,dy,resampleRatio= 1.,attribute= None):
		#-requires as input ...
		#-read in the source file name and decide on processing
		tempFileName= os.path.split(fileName)[1].lower()
		tempFileName= 'temp_%s' % os.path.splitext(tempFileName)[0]
		#-set output format
		outputFormat= 'PCRASTER_VALUESCALE=VS_%s' % outputFormat.upper()
		if '.shp' in fileName:
			#-shape file: create empty geotiff and rasterize prior to further processing
			shutil.copy('empty.tif','%s.tif' % tempFileName)
##      co= 'gdal_rasterize -a %s -ot %s -tr %f %f -te %f %f %f %f -l %s %s %s.tif -quiet' %\
##          (attribute,typeStr,dx,dy,x0,y0,x1,y1,os.path.split(fileName)[1].replace('.shp',''),fileName,tempFileName)
			co= 'gdal_rasterize -at  -a %s -l %s %s %s.tif -quiet ' %\
					(attribute,os.path.split(fileName)[1].replace('.shp',''),fileName,tempFileName)
			subprocess.check_call(co,stdout= subprocess.PIPE,stderr= subprocess.PIPE,shell= True)
			fileName= '%s.tif' % tempFileName
		#-process grid
		co= 'gdal_translate -ot %s -of PCRaster -mo %s %s %s.map -projwin %f %f %f %f -outsize %s %s -quiet' %\
				(typeStr,outputFormat,fileName,tempFileName,\
					x0,y1,x1,y0,'%.3f%s' % (100.*resampleRatio,'%'),'%.3f%s' % (100.*resampleRatio,'%'))
		#print co
		subprocess.check_call(co,stdout= subprocess.PIPE,stderr= subprocess.PIPE,shell= True)
		#-read resulting map
		setattr(self,variable,pcr.readmap('%s.map' % tempFileName))
		for suffix in ['map','tif']:
			try:
				os.remove('%s.%s' %(tempFileName,suffix))
				os.remove('%s.%s.aux.xml' %(tempFileName,suffix))
			except:
					pass

##    #-process grid
##    #-set output format
##    if 'FLOAT' in typeStr.upper():
##      outputFormat= 'PCRASTER_VALUESCALE=VS_SCALAR'
##    elif 'INT32' in typeStr.upper():
##      outputFormat= 'PCRASTER_VALUESCALE=VS_NOMINAL'
##    elif 'BYTE' in typeStr.upper():
##      outputFormat= 'PCRASTER_VALUESCALE=VS_BOOLEAN'
##    else:
##      sys.exit('grid format %s  is ambiguous' % typeStr)
