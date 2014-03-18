#-contains some basic PCRaster functions associated with PCR-GLOBWB
import  os, shutil, tarfile
import numpy as np
import pcraster as pcr
from pcraster.framework import generateNameT
from spatialDataSet2PCR import spatialDataSet

class readTSS:
  #-reads a TSS file and processes the data per time step accordingly

  def __init__(self,tssFileName):
    #-open text file for read and get contents
    txtFile= open(tssFileName)
    rawList= txtFile.read().split('\n')
    txtFile.close()
    #-determine whehther the file contains a header or not and return entries
    header= False
    iCnt= -1
    while len(rawList[iCnt]) == 0:
      iCnt-= 1   
    separator= ''
    if ',' in rawList[iCnt]:
      #-comma separated
      separator= ','
    try:
     nrEntries= int(rawList[1])
     header= True
    except:
      if separator != ',':
        nrEntries= len(rawList[1].split())
      else:
        nrEntries= len(rawList[1].split(separator))
    self.IDs= np.arange(1,nrEntries)
    if header == True:
      for iCnt in xrange(nrEntries-1):
        self.IDs[iCnt]= int(rawList[iCnt+3])
    #-read values
    headerRows= 0
    if header:
      headerRows= nrEntries+2
    if separator != ',':
      self.values= np.loadtxt(tssFileName, skiprows= headerRows)
    else:
      self.values= np.loadtxt(tssFileName, skiprows= headerRows,delimiter= separator)

  def assignID(self,ids,d,MV= -999.9):
    values= np.ones(ids.shape)*MV
    for ID in ids:
      if ID in self.IDs:
        values[ids == ID]= self.values[d-1,1:][self.IDs == ID]
    return values

def stackAverage(path,Root,StackStart,StackEnd):
	#calculates the average from a stack of maps, missing maps are skipped
	#Initialization
	MV= pcr.scalar(-999)
	NCount= pcr.scalar(0)
	SumStack= pcr.scalar(0)
	for StackNumber in range(StackStart,StackEnd):
		try:
			InMap= pcr.readmap(generateNameT(os.path.join(path,Root),StackNumber))
		except:
			InMap= MV;
		InMap= pcr.cover(InMap,MV)
		SumStack= SumStack+InMap
		NCount= NCount+pcr.ifthenelse(InMap <> MV,pcr.scalar(1),pcr.scalar(0))
	AvgStack= pcr.ifthenelse(NCount>0,SumStack/NCount,MV)
	return AvgStack

def createOutputDirs(DirList):
	#creates empty output directories
	for DirOut in DirList:
		if os.path.exists(DirOut):
			for Root, Dirs, Files in os.walk(DirOut):
				for File in Files:
					FileName= os.path.join(Root,File)
					try:
						os.remove(FileName)
					except:
						pass
		else:
			os.makedirs(DirOut)
	
def generateIniFiles(pathRes,YearDays):
	#this procedure does not take any arguments other than the input path and the number of days in
	#the year; it copies files reported at the last day of the year to the new input
	#files with specified extension
	SFX= 'ini'
	#setting last daily output to new initial conditions
	for FileNameIn in os.listdir(pathRes):
		if FileNameIn[-3:] == ('%03d' % YearDays):
			FileNameOut= os.path.join(pathRes,FileNameIn.replace('-','_'))
			FileNameOut= FileNameOut.replace('%03d' % YearDays ,SFX)
			FileNameIn= os.path.join(pathRes,FileNameIn)
			try:
				#copying input files
				shutil.copy(FileNameIn,FileNameOut)
			except:
				print 'file ', FileNameIn,' not found'
				pass

def copyIniFiles(pathIn,pathOut,IniExtList= [],SpecMapList= [],recur= 0):
	#copy files from in to out on the basis of the provided extension list; if empty all is copied
	#creating file inventory
	for Root, Dirs, Files in os.walk(pathIn):
		#Dstpath created if absent
		Dir= Root.find('/')+1
		if Dir > 0 and recur:
			DstPath= Root[Dir:len(Root)]
			DstPath= os.path.join(pathOut,DstPath)
		else:
			DstPath= pathOut
		#checking whether destination directory exist, if not it is created
		if not os.path.exists(DstPath):
			os.mkdir(DstPath)
		for File in Files:
			#check if File is a valid inital file
			if len(IniExtList) >0:
				CopyFile= False
			else:
				CopyFile=  True
			for ncnt in range(0,len(IniExtList)):
				Ext= File[-3:]
				if Ext == IniExtList[ncnt]:
					CopyFile=  True
			if CopyFile:
				SrcName = os.path.join(Root,File)
				DstName = os.path.join(DstPath,File)
				try:
					shutil.copy(SrcName,DstName)
					print '\tcopying %s => %s' %(SrcName,DstName)
				except:
					print '%s could not be copied' % SrcName
					pass
	for File in SpecMapList:
		SrcName = os.path.join(pathIn,File)
		DstName = os.path.join(pathOut,File)
		try:
			shutil.copy(SrcName,DstName)
			print '\tcopying %s => %s' %(SrcName,DstName)
		except:
			print '%s could not be copied' % SrcName
			pass

def extractDailyData(arcPath,arcFile,varPath,varList,leapYear):
	#-extracts daily data from tar file
	nDays= 365
	if leapYear: nDays= nDays+1
	arcList= []
	for var in varList:
		for day in range(1,nDays+1):
			arcList.append(generateNameT(var,day))
	archive= tarfile.open(os.path.join(arcPath,arcFile),'r:gz')
	for tarInfo in archive:
			if tarInfo.name in arcList:
				archive.extract(tarInfo,varPath)
	archive.close()

class pcrPartialRead(object):
	"""resamples and scales a map to the extent specified and reads the resulting map for the specified clone map"""
	
	def __init__(self,xMin,xMax,yMin,yMax,gridResolution,resampleRatio,clone):
		#sets the clone and coordinates
		self.xMin= xMin
		self.xMax= xMax
		self.yMin= yMin
		self.yMax= yMax
		self.gridResolution= gridResolution
		self.resampleRatio= resampleRatio
		self.clone= clone
		self.varTypes= {}
		self.varTypes['boolean']= 'BYTE' 
		self.varTypes['ldd']= 'BYTE'
		self.varTypes['nominal']= 'INT32'
		self.varTypes['nominal']= 'INT32'
		self.varTypes['ordinal']= 'INT32'
		self.varTypes['scalar']= 'FLOAT32'
		self.varTypes['directional']= 'FLOAT32'
		
	def get(self,fileName,valueScale= 'scalar'):
		#-set variable and varType
		varType= self.varTypes.get(valueScale.lower(),'FLOAT32')
		variable= 'valueField'
		field= getattr(spatialDataSet(fileName,variable,varType,valueScale,\
			self.xMin,self.yMin,self.xMax,self.yMax,self.gridResolution,self.gridResolution,self.resampleRatio),variable)
		#if '.map' in fileName:
			#pcr.report(pcr.ifthen(self.clone,field),'temp_%s' % os.path.split(fileName)[1])
		return pcr.ifthen(self.clone,field)
