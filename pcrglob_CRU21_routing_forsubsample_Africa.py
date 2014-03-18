#-Runs PCR-GLOBWB using the ERA-Interim Reanalysis over the period 1979-2010 inclusive
# using mean daily temperature and daily GPCP-corrected rainfall totals as well as monthly reference
# potential evapotranspiration based on Penman-Monteith
# NOTE: the model is run by a generic water balance model that takes daily input, including reference potential evapotranspiration
# 1) daily meto input is retrieved by symbolic links
#     monthly potential evapotranspiration is scaled on the basis of Hamon's temperature based method;
# 2) the model is initially run over the entire, available reanalysis  period  until the cell-based average
#     specific runoff from the land surface converges within 1% for the last year (note that runoff over the
#     water surface is instantaneaous and does not change from year-to-year); for thespin-up, the report interval
#     is set to yearly only, except for the last year of the reanalysis period; for the definite run, all years are reported
#     on a daily timestep and stored;
# 3) in the definite run, also the routing model is applied (standard parameterization); from the runs, monthly storages
#     are retrieved for comparsion with GRACE

##-Modules used
import  os, sys,shutil, math, calendar,tarfile,datetime, zlib, zipfile
import pcraster as pcr
import pcraster.framework as pcrm
from pcrBasicFunctions import *
from pcrglobwb_WaterBodies_forsubsample import *
from random import choice

#-functions
def estimateBankFulDischarge(Q,constant= 9.8874,exponent= 0.8363):
	return constant*Q**exponent

def extractWaterBodyParameterization(zipFileName,extractionPath,matchList):
	#-function without child, extracting matches from archive
	zipArchive= zipfile.ZipFile(zipFileName)
	matchedFiles= []
	for sourceFileName in zipArchive.namelist():
		for match in matchList:
			if match in sourceFileName:
				matchedFiles.append(sourceFileName)
	for sourceFileName in matchedFiles:
		targetFileName= os.path.join(extractionPath,os.path.split(sourceFileName)[1])
		tempFile= open(targetFileName,'wb')
		tempFile.write(zipArchive.read(sourceFileName))
		tempFile.close()
	#-close file
	zipArchive.close()
	
#-Main
#-Initialization
#-forcing dataset to be used
forcingDataSet= 'VIC'
domainStr= 'Africa'

#-date management: start and end date and corresponding time steps
# and duration (days) of a model step
duration= 1.
timeSec= 86400
startYear= 2000
endYear= 2005
yearRange= range(startYear,endYear)
lastYear= endYear-1
initYear= startYear-1
initYears= range(startYear,startYear+1)
julianDay= [1,32,60,91,121,152,182,213,244,274,305,335,366]
addDay= [0,0,1,1,1,1,1,1,1,1,1,1,1]
#-paths
#-results and storage
mapsDir= 'maps'
tempDir= 'temp'
pathIni= 'initialconditions'
pathRes= 'results'
pathResMonthly= 'resmonthly'
pathResZip= 'reszip'
pathRemote= ''

#-model outputs for import and archiving
resStackList= ['qloc','qw','qc','wst','vel','fldd','fldf','wata','watd','rest','mbe','mbr']
qlocStart= 0
iniOutVarStart= 2
addOutVarStart= 4
extFile= '.001'
# - archives with daily results stored as tar.gz and zip file with parameterization
dailyArchiveFileName= 'pcrglobwb_%s_%04d_30min.tar.gz'
reservoirParameterizationZipFileName= '../reservoirparameterization.zip'

##-file names used by the routing model, this is messy!
#-initialization: boolean settings
noReservoirs= False
testLocalWaterBalance= False
reportLocalWaterBalance= False
getSurfaceWaterAttributes= True
initializeRoutingModel= False
#-number of runs to initialize the model
nrRunsInitialization= 2
#-weight to update long-term average discharge and bankful values
updateQWeight= 0.2
#-maps: clone map and cell area
clone= pcr.readmap(os.path.join(mapsDir,'cloneAfrica.map'))
#-set sample domain
xMin= -25.5
xMax= 57.5	
yMin= -35.0
yMax= 37.5
sampleResolution= 0.5
resampleRatio= 1.
#-class to clip map during read
clippedRead= pcrPartialRead(xMin,xMax,yMin,yMax,sampleResolution,resampleRatio,clone)
cellAreaMap= os.path.join(mapsDir,'cellarea30.map')
cellArea= clippedRead.get(cellAreaMap,'scalar')
#-channel and floodplain characteristics
LDDMap= os.path.join(mapsDir,'lddsound_30min.map')	
LDD= pcr.lddrepair(clippedRead.get(LDDMap,'ldd'))	
channelGradient= clippedRead.get(os.path.join(mapsDir,'globalgradchannel.map'))
channelWidth= clippedRead.get(os.path.join(mapsDir,'channel_width.map'))
channelLength= clippedRead.get(os.path.join(mapsDir,'channel_length.map'))
channelDepth= clippedRead.get(os.path.join(mapsDir,'channel_depth.map'))
floodplainMask= pcr.spatial(pcr.boolean(1)) # NOTE: set to zero for static, to one for dynamic floodplains
channelManN= 0.04
floodplainManN= 0.10
#-flood plain parameterization
#-root of file name with maps of relative elvation above floodplain
# and associated fractions
relZFileName= 'elev%04d.map'
areaFractions=[0.0,0.01,0.05,0.10,0.20,0.30,0.40,\
  0.50,0.60,0.70,0.80,0.90,1.00]
# reduction parameter of smoothing interval and error threshold
reductionKK= 0.5
criterionKK= 40.
#-modelSignature
if pcr.cellvalue(pcr.mapmaximum(pcr.scalar(floodplainMask)),1)[0] == 1:
	modelSignature= forcingDataSet+'_dynamic-routing'
else:
	modelSignature= forcingDataSet+'_static-routing'
if noReservoirs:
  modelSignature= modelSignature+'_noreservoirs'
modelSignature+= '_%s' % domainStr
#-stacks of specific runoff and direct flux over water surface [m per unit area and per time step]
# for model input and root of maps of initial storage [m3], discharge [m3/s], flood depth [m] and
# flood fraction [m2/m2]
landSurfaceQFileName= os.path.join(pathRes,'qloc')
waterSurfaceQFileName= os.path.join(pathRes,'qw')
actualStorageFileName= os.path.join(pathRes,'wst')
QFileName= os.path.join(pathRes,'qc')
flowVelocityFileName= os.path.join(pathRes,'vel')
floodedDepthFileName= os.path.join(pathRes,'fldd')
floodedFractionFileName= os.path.join(pathRes,'fldf')
surfaceWaterAreaFileName= os.path.join(pathRes,'wata')
surfaceWaterDepthFileName= os.path.join(pathRes,'watd')
surfaceWaterResidenceTimeFileName= os.path.join(pathRes,'rest')
mbeFileName= os.path.join(pathRes,'mbe.map')
mbrFileName= os.path.join(pathRes,'mbr.map')
  
#-Start
#-echo to screen and file
currentDate= datetime.datetime.now()
logFileName= 'pcr-globwb_%s_%04d-%02d-%02d_30min.log' % (modelSignature,currentDate.year,currentDate.month,currentDate.day)
logFile= open(logFileName,'w')
wStr= 'Running routing for PCR-GLOBWB large-scale hydrological model forced with %s over %d-%d inclusive\n' %\
	(modelSignature.replace('_',' '),startYear,lastYear)
print wStr
logFile.write(wStr+'\n')
#-create output directory for initial conditions
try:
	os.makedirs(os.path.join(pathIni,modelSignature))
except:
	pass
#-create temporary output directory
try:
	createOutputDirs([tempDir])
except:
	pass
#-read in maps with water body parameterization for initialization
extractWaterBodyParameterization(reservoirParameterizationZipFileName,tempDir,[str(initYear)])
#-initialize the model if specified
if initializeRoutingModel:
	#-initialization
	pathIniSource= os.path.join(pathIni,modelSignature,'ini_%04d' % initYear)
	if not os.path.exists(pathIniSource):
		createOutputDirs([pathIniSource])
	#-long-term values of evapotranspiration and discharge from monthly files
	wStr= '* initialization:  extracting long-term evaporation and discharge for %d-%d' % (min(initYears),max(initYears))
	print wStr
	logFile.write(wStr+'\n')
	#-fraction water surface needed to initialize long-term discharge and open water evaporation
	fractionWater= clippedRead.get(os.path.join(tempDir,'fracwat_%d.map' % initYear))
	#-discharge in m3/s, runoff over water surface in m3/year
	averageQ= pcr.scalar(0)
	averageQWat= pcr.scalar(0)
	yearCnt= 0
	for year in initYears:
		yearCnt+= 1
		#-extract input from archive
		dailyArchiveInputFileName=  dailyArchiveFileName % (forcingDataSet,year)
		print '- processing archive %s' %  dailyArchiveInputFileName
		extractDailyData(os.path.join(pathRemote,forcingDataSet),dailyArchiveInputFileName,pathRes,\
			['q1mon','q2mon','q3mon','qwmon'],calendar.isleap(year))
		#-update long-term values
		for month in xrange(1,13):
			day= julianDay[month]
			if calendar.isleap(year):
				day+= addDay[month]
			day-= 1
			monthLength= calendar.monthrange(year,month)[1]
			#-local runoff and evaporation in m3/month
			qloc= pcr.scalar(0)
			qwat= pcr.scalar(0)
			for fileRoot in ['q1mon','q2mon','q3mon']:
				qloc+= (1.-fractionWater)*cellArea*clippedRead.get(generateNameT(os.path.join(pathRes,fileRoot),day))
			for fileRoot in ['qwmon']:
				qwat+= fractionWater*cellArea*clippedRead.get(generateNameT(os.path.join(pathRes,fileRoot),day))
			averageQ+= pcr.accuthresholdflux(LDD,qloc+pcr.max(0.,qwat),pcr.max(0.,-qwat))/(12*monthLength*24*3600)
			averageQWat+= qwat
	#-update maps for number of years and compute bankful discharge
	averageQ/= yearCnt
	averageQIni= averageQ
	averageQ= pcr.ifthenelse(clippedRead.get(os.path.join(tempDir,'waterbodyoutlet_%d.map' % initYear),'nominal') != 0,\
		averageQ,0)
	averageQWat/= yearCnt
	bankfulQ= estimateBankFulDischarge(averageQ)	
	#-report maps
	targetFileName= os.path.join(pathIniSource,'qavg_longterm.map')
	pcr.report(averageQ,targetFileName)
	targetFileName= os.path.join(pathIniSource,'qbank_longterm.map')
	pcr.report(bankfulQ,targetFileName)
	targetFileName= os.path.join(pathIniSource,'qwat_avg_longterm.map')
	pcr.report(averageQWat,targetFileName)
#-define initital values for waterbodies
#-initializing class of water bodies for the first year
if not initializeRoutingModel:
	pathIniSource= os.path.join(pathIni,modelSignature,'ini_%04d' % startYear)
	fractionWater= clippedRead.get(os.path.join(tempDir,'fracwat_%d.map' % initYear))
waterBodies= pcrglobWaterBodies(os.path.join(tempDir,'waterbodyid_%d.map' % initYear),\
	os.path.join(tempDir,'waterbodyoutlet_%d.map' % initYear),os.path.join(tempDir,'waterbodytype_%d.map' % initYear),\
	os.path.join(mapsDir,'channel_width.map'),os.path.join(pathIniSource,'qavg_longterm.map'),\
	os.path.join(pathIniSource,'qbank_longterm.map'),LDDMap,os.path.join(tempDir,'reservoirparameters_%d_30min_.tbl' % initYear),timeSec,clippedRead)
#-create initial files if the routing model is initialized
if initializeRoutingModel:
	for variable in resStackList[iniOutVarStart:addOutVarStart]:
		targetFileName= os.path.join(pathIniSource,pcrm.generateNameT(variable,0).replace('.000','.ini'))
		if 'qc' in variable:
			pcr.report(pcr.scalar(0),targetFileName)
		elif 'wst'in variable:
			waterBodiesType= pcr.nominal(waterBodies.returnMapValue(pcr.spatial(pcr.scalar(0.)),waterBodies.type))
			endorheicLakes= pcr.ifthenelse((pcr.areatotal(pcr.scalar(waterBodies.outlet != 0),waterBodies.distribution) == 0) & \
				(waterBodies.distribution != 0),waterBodies.distribution,0)
			#-storage of water bodies is assigned at outlets
			#-storage over endorheic lakes
			actualStorage= pcr.ifthen(endorheicLakes != 0,\
				pcr.upstream(LDD,pcr.ifthenelse((endorheicLakes == 0) & ((waterBodies.distribution == 0) | (waterBodies.outlet != 0)),\
				averageQ,0)))
			actualStorage= pcr.ifthen((waterBodies.location != 0) & (endorheicLakes != 0),\
				pcr.max(0.,pcr.areatotal(pcr.cover(actualStorage*timeSec*365.25+averageQWat,0),endorheicLakes)))
			#-storage over other water bodies
			actualStorage= pcr.cover(actualStorage,\
				pcr.ifthen((waterBodiesType == 2) & (waterBodies.location != 0),\
				waterBodies.returnMapValue(pcr.spatial(pcr.scalar(0.)),waterBodies.capacity*waterBodies.maxLimit)),\
				pcr.ifthen((waterBodiesType == 1) & (waterBodies.location != 0),\
				pcr.areatotal(fractionWater*cellArea,waterBodies.distribution)*(averageQ/(waterBodies.cLake*channelWidth))**(2./3.)))
			actualStorage= pcr.ifthen(waterBodies.distribution != 0,pcr.cover(actualStorage,0))
			#-storage over rivers
			actualStorage= pcr.cover(actualStorage,channelLength* (channelManN*averageQ*channelWidth**(2./3.)*channelGradient**-0.5 )**0.6)
			pcr.report(actualStorage,targetFileName)
#-initial conditions created, execute the followingl until exitCode becomes true
# set run counter and completion identifier
runCnt= 0
exitCode= False
while not exitCode:
	#-pre-processing: selection of year and in- and output locations depends on
	# whether the model is initialized or actually ran for a time series
	createOutputDirs([pathRes])
	if initializeRoutingModel:
		#-initialization: take random year from selection and set in and output of initial conditions
		# to year prior to start run
		wStr= 'initialization routing model'
		year= choice(initYears)
		pathIniSource=os.path.join(pathIni,modelSignature,'ini_%04d' % initYear)
		if (runCnt+1) < nrRunsInitialization:
			pathIniTarget= pathIniSource
		else:
			pathIniTarget= os.path.join(pathIni,modelSignature,'ini_%04d' % (startYear))
			createOutputDirs([pathIniTarget])
	else:
		wStr= 'actual run'
		#-run the model for actual sequence of years, output location incremented by one year
		year= yearRange.pop(0)
		pathIniSource=os.path.join(pathIni,modelSignature,'ini_%04d' % year)
		pathIniTarget= os.path.join(pathIni,modelSignature,'ini_%04d' % (year+1))
		createOutputDirs([pathIniTarget])
	#-initializing, calling and running the model:
	#-set start and end date and number of days
	startDate= datetime.datetime(year,1,1)
	endDate= datetime.datetime(year,12,31)
	yearDays= endDate.toordinal()-startDate.toordinal()+1
	##-extract input from archive
	dailyArchiveInputFileName=  dailyArchiveFileName % (forcingDataSet,year)
	print '* processing archive %s' %  dailyArchiveInputFileName
	#extractDailyData(os.path.join(pathRemote,forcingDataSet),dailyArchiveInputFileName,pathRes,\
	#	resStackList[qlocStart:iniOutVarStart],calendar.isleap(year))
	#-copy initial conditions
	copyIniFiles(pathIniSource,pathRes,['ini'])
	#-retrieve average long-term discharge and estimate bankful discharge
	averageQ=  pcr.readmap(os.path.join(pathIniSource,'qavg_longterm.map'))
	bankfulQ= estimateBankFulDischarge(averageQ)
	targetFileName= os.path.join(pathIniSource,'qbank_longterm.map')
	pcr.report(bankfulQ,targetFileName)
	#-call model
	execfile('pcr_rout_dynfloodplain_v3_forsubsample.py')
	routingscheme= pcrglobRoutingDynFloodPlains(startDate,endDate)
	dynRouting = pcrm.DynamicFramework(routingscheme,lastTimeStep= int((endDate.toordinal()-startDate.toordinal()+1)/duration),\
		firstTimestep= 1)
	dynRouting.run()
	#-run completed, update waterbodies for current year and post-process
	currentDate= datetime.datetime.now()
	wStr= '* run:%6d, year:%6d, end: %04d-%02d-%02d %02d:%02d %30s' %\
		(runCnt,year,currentDate.year,currentDate.month,currentDate.day,currentDate.hour,currentDate.minute,wStr)
	#-update waterbody parameters
	if not initializeRoutingModel:
		#-read in maps with water body parameterization for initialization and initialize water bodies for next year
		createOutputDirs([tempDir])
		extractWaterBodyParameterization(reservoirParameterizationZipFileName,tempDir,[str(year)])
	#-update long-term average discharge when initialization complete, resetting it to zero for changed outlets
	if runCnt > 0:
		if initializeRoutingModel:
			changedOutlet= (waterBodies.outlet != 0) ^ (clippedRead.get(os.path.join(tempDir,'waterbodyoutlet_%d.map' % initYear)) != 0)
		else:
			changedOutlet= (waterBodies.outlet != 0) ^ (clippedRead.get(os.path.join(tempDir,'waterbodyoutlet_%d.map' % year)) != 0)
		targetFileName= os.path.join(pathIniTarget,'updated_outlet.map')
		pcr.report(changedOutlet,targetFileName)	
		try:
			targetFileName= os.path.join(pathIniTarget,'qavg_longterm.map')
			Q= stackAverage(pathRes,resStackList[iniOutVarStart],1,yearDays)
			averageQ= (1.-updateQWeight)*averageQ+updateQWeight*Q
			averageQ= pcr.ifthenelse(changedOutlet,averageQIni,averageQ)
			pcr.report(averageQ,targetFileName)
			#-update wStr
			wStr+= 'total runoff [km3]: %.2e' % \
				pcr.cellvalue(pcr.maptotal(pcr.ifthen(LDD ==  5,Q*1.e-9*yearDays*timeSec)),1)[0]
		except:
			wStr+= '\n - long-term discharge has not been updated!'
	if not initializeRoutingModel:
		#-read in maps with water body parameterization for initialization and initialize water bodies for next year
		del waterBodies
		waterBodies= pcrglobWaterBodies(os.path.join(tempDir,'waterbodyid_%d.map' % year),\
			os.path.join(tempDir,'waterbodyoutlet_%d.map' % year),os.path.join(tempDir,'waterbodytype_%d.map' % year),\
			os.path.join(mapsDir,'channel_width.map'),os.path.join(pathIniSource,'qavg_longterm.map'),\
			os.path.join(pathIniSource,'qbank_longterm.map'),LDDMap,os.path.join(tempDir,'reservoirparameters_%d_30min_.tbl' % year),timeSec,clippedRead)
		fractionWater= clippedRead.get(os.path.join(tempDir,'fracwat_%d.map' % year))
		wStr+= '\n - waterbodies updated for %d' % (year+1)
	#-write log string to screen and file
	print wStr
	wStr+= '\n'
	logFile.write(wStr)
	logFile.flush()
	#-update initial conditions
	targetFileName=generateNameT(actualStorageFileName,0).replace('.000','.ini')
	shutil.copy(targetFileName,targetFileName.replace('.ini','org.ini'))
	actualStorage= pcr.readmap(targetFileName)
	actualStorage= pcr.ifthenelse(waterBodies.distribution != 0,\
		pcr.ifthenelse(waterBodies.location != 0, pcr.areatotal(actualStorage,waterBodies.distribution),0.),\
		actualStorage)
	pcr.report(actualStorage,targetFileName)
	targetFileName= generateNameT(QFileName,0).replace('.000','.ini')
	shutil.copy(targetFileName,targetFileName.replace('.ini','org.ini'))
	Q= pcr.readmap(targetFileName)
	Q= pcr.ifthenelse(waterBodies.distribution != 0,\
		pcr.ifthenelse(waterBodies.location != 0, Q,0),Q)
	pcr.report(Q,targetFileName)
	#-move initial velocity estimate
	sourceFileName= generateNameT(flowVelocityFileName,0).replace('.000','.ini')
	targetFileName= sourceFileName.replace(pathRes,pathIniSource)
	shutil.move(sourceFileName,targetFileName)
	#-copy initial conditions and archive output
	copyIniFiles(pathRes,pathIniTarget,['ini'])
	"""if not initializeRoutingModel:
		dailyArchiveOutputFileName=  os.path.join(pathResZip,\
			dailyArchiveFileName % (modelSignature,year))
		print '\t* archiving daily output to %s' % dailyArchiveOutputFileName
		dailyArchive= tarfile.open(dailyArchiveOutputFileName,'w:gz')
		#-add files
		for fileName in os.listdir(pathRes):
			sourceFileName= os.path.join(pathRes,fileName)
			if not '.ini' in fileName:
				for variable in resStackList[iniOutVarStart:]:
					if variable in fileName:
						dailyArchive.add(sourceFileName,fileName)
		for fileName in os.listdir(pathIniSource):
			sourceFileName= os.path.join(pathIniSource,fileName)
			dailyArchive.add(sourceFileName,sourceFileName)
		dailyArchive.close()
		try:
			print '\t* transferring output',
			shutil.move(os.path.join(pathResZip,dailyArchiveFileName % (modelSignature,year)),\
				os.path.join(pathRemote,forcingDataSet,dailyArchiveFileName % (modelSignature,year)))
			print
		except:
			print 'has not succeeded'
	"""
        #-increment runCnt and reset initialization
	# and decide to stop run when no actual years are left
	runCnt+= 1
	if runCnt >= nrRunsInitialization: initializeRoutingModel= False
	if len(yearRange) == 0: exitCode= True
#-close log file, empty output directories and terminate
#createOutputDirs([pathRes,tempDir])
#-close log file and moving files
wStr= 'run terminated at %s' % datetime.datetime.now()
print wStr
logFile.write(wStr+'\n')
logFile.close()
"""print 'moving output'
for fileName in os.listdir(pathResZip):
	if '.zip' in fileName or '.tar' in fileName:
		sourceFileName= os.path.join(pathResZip,fileName)
		targetFileName= os.path.join(pathRemote,forcingDataSet,fileName)
		shutil.move(sourceFileName,targetFileName)
shutil.copy(logFileName,os.path.join(pathRemote,forcingDataSet,logFileName))
print 'all processed!'
"""
