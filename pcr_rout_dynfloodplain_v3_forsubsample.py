#-PCR-GLOBWB Routing scheme: basic implementation for python version
# with the following characteristics:
# * total water storage [m3] of which volume in excess of channel storage is used as flood volume
#   and distributed over the surrounding area using smooth floodplain transitions based on Kavetski and Kuczera (2007)
# * flood wave propagation can be described by means of the kinematic wave or a simpler traveltime function
# * lakes and reservoirs are included by means of the waterbodies class

import os, datetime, calendar
import pcraster as pcr
import pcraster.framework as pcrm

#############################################################################
#-start of model script #######################################################
################################################################################

class pcrglobRoutingDynFloodPlains(pcrm.DynamicModel):

	def __init__(self,startDate,endDate,nrIterDefault= 12.):
		pcrm.DynamicModel.__init__(self)
		##############
		# * __init__ #
		##############
		#-echo to screen
		print 'PCR-GLOBWB dynamic floodplain model - version 3.0 June 2012'
		print '\tbeta version with smoothed floodplain elevations'
		print '\tincluding lakes and reservoirs'
		#-constants: duration of time step (day) in seconds 
		self.timeSec= duration*timeSec
		#-passing global variables to local ones
		#-cloone and cell area
		self.clone= clone
		self.cellArea= cellArea
		#-model settings
		self.nrIterDefault= nrIterDefault
		self.duration= duration
		self.currentDate= startDate
		self.endDate= endDate
		self.areaFractions= areaFractions
		self.relZFileName= os.path.join(mapsDir,relZFileName)
		self.reductionKK= reductionKK
		self.criterionKK= criterionKK
		#-number of entries in list and derived slopes and volumes
		self.nrEntries= len(self.areaFractions)
		self.relZ= [0.]*self.nrEntries
		self.floodVolume= [0.]*(self.nrEntries)    
		#-flood plain and channel characteristics
		self.channelGradient= pcr.max(1.e-7,channelGradient)
		self.LDD= LDD
		self.channelWidth= channelWidth
		self.channelLength= channelLength
		self.channelDepth= channelDepth
		self.channelManN= channelManN
		self.floodplainManN= floodplainManN
		self.floodplainMask= floodplainMask
		self.waterFractionMask= fractionWater
		#-waterBodies
		self.waterBodies= waterBodies
		self.waterBodies.actualArea= self.waterBodies.retrieveMapValue(pcr.areatotal(self.waterFractionMask*\
			self.cellArea,self.waterBodies.distribution))
		#self.reservoirDemandTSS= readTSS(reservoirDemandTSS)
		#-map names: initial maps of discharge and storage
		self.QIniMap= pcrm.generateNameT(QFileName,0).replace('.000','.ini')
		self.actualStorageIniMap= pcrm.generateNameT(actualStorageFileName,0).replace('.000','.ini')
		self.averageQ= averageQ
		self.bankfulQ= bankfulQ
		#-patch elevations: those that are part of sills are updated on the basis of the floodplain gradient
		# using local distances deltaX per increment upto z[N] and the sum over sills
		#-fill all lists including smoothing interval and slopes
		for iCnt in range(1,self.nrEntries):
			self.relZ[iCnt]= clippedRead.get(self.relZFileName %\
				(self.areaFractions[iCnt]*100))
		#-minimum slope of floodplain, being defined as the longest sill, first used to retrieve
		# longest cumulative distance 
		deltaX= [self.cellArea**0.5]*self.nrEntries
		deltaX[0]= 0.
		sumX= deltaX[:]
		minSlope= 0.
		for iCnt in range(self.nrEntries):
			if iCnt < self.nrEntries-1:
				deltaX[iCnt]= (self.areaFractions[iCnt+1]**0.5-self.areaFractions[iCnt]**0.5)*deltaX[iCnt]
			else:
				deltaX[iCnt]= (1.-self.areaFractions[iCnt-1]**0.5)*deltaX[iCnt]
			if iCnt > 0:
				sumX[iCnt]= pcr.ifthenelse(self.relZ[iCnt] == self.relZ[iCnt-1],sumX[iCnt-1]+deltaX[iCnt],0.)
				minSlope= pcr.ifthenelse(self.relZ[iCnt] == self.relZ[iCnt-1],\
					pcr.max(sumX[iCnt],minSlope),minSlope)
		minSlope= pcr.min(self.channelGradient,0.5*pcr.max(deltaX[1],minSlope)**-1.)
		#-add small increment to elevations to each sill except in the case of lakes
		for iCnt in range(self.nrEntries):
			self.relZ[iCnt]= self.relZ[iCnt]+sumX[iCnt]*pcr.ifthenelse(self.relZ[self.nrEntries-1] > 0.,\
				minSlope,0.)
		#-set slope and smoothing interval between dy= y(i+1)-y(i) and dx= x(i+1)-x(i)
		# on the basis of volume
		#-slope and smoothing interval
		self.kSlope=  [0.]*(self.nrEntries)
		self.mInterval= [0.]*(self.nrEntries)
		for iCnt in range(1,self.nrEntries):
			self.floodVolume[iCnt]= self.floodVolume[iCnt-1]+\
				0.5*(self.areaFractions[iCnt]+self.areaFractions[iCnt-1])*\
				(self.relZ[iCnt]-self.relZ[iCnt-1])*self.cellArea
			self.kSlope[iCnt-1]= (self.areaFractions[iCnt]-self.areaFractions[iCnt-1])/\
				pcr.max(0.001,self.floodVolume[iCnt]-self.floodVolume[iCnt-1])
		for iCnt in range(1,self.nrEntries):
			if iCnt < (self.nrEntries-1):
				self.mInterval[iCnt]= 0.5*self.reductionKK*pcr.min(self.floodVolume[iCnt+1]-self.floodVolume[iCnt],\
					self.floodVolume[iCnt]-self.floodVolume[iCnt-1])
			else:
				self.mInterval[iCnt]= 0.5*self.reductionKK*(self.floodVolume[iCnt]-self.floodVolume[iCnt-1])

	def integralLogisticFunction(self,x):
		#-returns a tupple of two values holding the integral of the logistic functions
		# of (x) and (-x)
		logInt=pcr.ln(pcr.exp(-x)+1)
		return logInt,x+logInt

	def returnFloodedFraction(self,volume):
		#-returns the flooded fraction given the flood volume and the associated water height
		# using a logistic smoother near intersections (K&K, 2007)
		#-find the match on the basis of the shortest distance to the available intersections or steps
		deltaXMin= self.floodVolume[self.nrEntries-1]
		y_i= pcr.scalar(1.)
		k= [pcr.scalar(0.)]*2
		mInt= pcr.scalar(0.)
		for iCnt in range(self.nrEntries-1,0,-1):
			#-find x_i for current volume and update match if applicable
			# also update slope and intercept
			deltaX= volume-self.floodVolume[iCnt]
			mask= pcr.abs(deltaX) < pcr.abs(deltaXMin)
			deltaXMin= pcr.ifthenelse(mask,deltaX,deltaXMin)
			y_i= pcr.ifthenelse(mask,self.areaFractions[iCnt],y_i)
			k[0]= pcr.ifthenelse(mask,self.kSlope[iCnt-1],k[0])
			k[1]= pcr.ifthenelse(mask,self.kSlope[iCnt],k[1])
			mInt= pcr.ifthenelse(mask,self.mInterval[iCnt],mInt)
		#-all values returned, process data: calculate scaled deltaX and smoothed function
		# on the basis of the integrated logistic functions PHI(x) and 1-PHI(x)
		deltaX= deltaXMin
		deltaXScaled= pcr.ifthenelse(deltaX < 0.,pcr.scalar(-1.),1.)*\
			pcr.min(criterionKK,pcr.abs(deltaX/pcr.max(1.,mInt)))
		logInt= self.integralLogisticFunction(deltaXScaled)
		#-compute fractional flooded area and flooded depth
		floodedFraction= pcr.ifthenelse(volume > 0.,\
			pcr.ifthenelse(pcr.abs(deltaXScaled) < criterionKK,\
			y_i-k[0]*mInt*logInt[0]+k[1]*mInt*logInt[1],\
			y_i+pcr.ifthenelse(deltaX < 0.,k[0],k[1])*deltaX),0.)
		floodedFraction= pcr.max(0.,pcr.min(1.,floodedFraction))
		floodDepth= pcr.ifthenelse(floodedFraction > 0.,volume/(floodedFraction*self.cellArea),0.)
		return floodedFraction, floodDepth

	def kinAlphaStatic(self,watStor):
		#-given the total water storage in the cell, returns the Q-A relationship
		# for the kinematic wave and required parameters using a static floodplain extent
		wetA= watStor/self.channelLength
		wetP= 2.*wetA/self.channelWidth+self.channelWidth
		alphaQ= (self.channelManN*wetP**(2./3.)*self.channelGradient**-0.5)**self.betaQ
		#-returning variable of interest: flooded fraction, cross-sectional area
		# and alphaQ
		return wetA,alphaQ

	def kinAlphaDynamic(self,watStor):
		#-given the total water storage in the cell, returns the Q-A relationship
		# for the kinematic wave and required parameters
		floodVol= pcr.max(0,watStor-self.channelStorageCapacity)
		floodFrac, floodZ= self.returnFloodedFraction(floodVol)
		#-wetted perimeter, cross-sectional area and
		# corresponding mannings' n
		wetA= watStor/self.channelLength
		#-wetted perimeter, alpha and composite manning's n
		wetPFld= pcr.max(0.,floodFrac*self.cellArea/self.channelLength-\
			self.channelWidth)+2.*floodZ
		wetPCh= self.channelWidth+\
			2.*pcr.min(self.channelDepth,watStor/self.channelArea)
		wetP= wetPFld+wetPCh
		manQ= (wetPCh/wetP*self.channelManN**1.5+\
			wetPFld/wetP*self.floodplainManN**1.5)**(2./3.)
		alphaQ= (manQ*wetP**(2./3.)*self.channelGradient**-0.5)**self.betaQ
		#-returning variables of interest: flooded fraction, cross-sectional area
		# and alphaQ
		return floodFrac,floodZ,wetA,alphaQ

	def kinAlphaComposite(self,watStor,mask):
		#-given the total water storage and the mask specifying the occurrence of
		# floodplain conditions, retrns the Q-A relationship for the kinematic
		# wave and the associated parameters
		floodplainStorage= pcr.ifthen(mask,watStor)
		floodFrac, floodZ, dynamicWetA, dynamicAlphaQ= self.kinAlphaDynamic(floodplainStorage)
		staticWetA, staticAlphaQ= self.kinAlphaStatic(watStor)
		floodFrac= pcr.ifthenelse(mask,floodFrac,0.)
		floodZ= pcr.ifthenelse(mask,floodZ,0.)
		wetA= pcr.ifthenelse(mask,dynamicWetA,staticWetA)
		alphaQ= pcr.ifthenelse(mask,dynamicAlphaQ,staticAlphaQ)
		return floodFrac,floodZ,wetA,alphaQ

	def initial(self):
		#####################
		# * initial section #
		#####################
		#-constants
		# betaQ [-]: constant of kinematic wave momentum equation
		self.betaQ= 0.6
		#-channel LDD
		self.channelLDD= pcr.ifthenelse(self.waterBodies.distribution != 0,\
			pcr.ldd(5),self.LDD)
		#-channel area and storage
		self.channelArea= self.channelWidth*self.channelLength
		self.channelStorageCapacity= pcr.ifthenelse(self.waterBodies.distribution == 0,\
			self.channelArea*self.channelDepth,pcr.scalar(0.))
		#-basin outlets
		self.basinOutlet= pcr.pit(self.LDD) != 0
		#-read initial conditions
		self.Q= clippedRead.get(self.QIniMap)
		self.actualStorage= clippedRead.get(self.actualStorageIniMap)
		self.actualStorage= pcr.ifthenelse(self.waterBodies.distribution != 0,\
			pcr.ifthenelse(self.waterBodies.location != 0,\
				pcr.areatotal(self.actualStorage,self.waterBodies.distribution),0),\
					self.actualStorage)   
		self.waterBodies.actualStorage= self.waterBodies.retrieveMapValue(self.actualStorage)
		#-update targets of average and bankful discharge
		self.waterBodies.averageQ= self.waterBodies.retrieveMapValue(self.averageQ)
		self.waterBodies.bankfulQ= self.waterBodies.retrieveMapValue(self.bankfulQ)
		#-return the parameters for the kinematic wave,
		# including alpha, wetted area, flood fraction, flood volume and depth
		# and the corresponding land area
		floodedFraction,floodedDepth,\
			self.wettedArea,self.alphaQ= self.kinAlphaComposite(self.actualStorage,self.floodplainMask)
		self.wettedArea= self.waterBodies.returnMapValue(self.wettedArea,\
			self.waterBodies.channelWidth+2.*self.waterBodies.updateWaterHeight())
		self.waterFraction= pcr.ifthenelse(self.waterBodies.distribution == 0,\
			pcr.max(self.waterFractionMask,floodedFraction),self.waterFractionMask)
		self.landFraction= pcr.max(0.,1.-self.waterFraction)
		#-update on velocity and check on Q - NOTE: does not work in case of reservoirs!
		self.flowVelocity= pcr.ifthenelse(self.wettedArea > 0,self.Q/self.wettedArea,0.)
		pcr.report(self.flowVelocity,pcrm.generateNameT(flowVelocityFileName,0).replace('.000','.ini'))
		#-setting initial values for specific runoff and surface water extraction
		self.landSurfaceQ= pcr.scalar(0.)
		self.potWaterSurfaceQ= pcr.scalar(0.)
		self.surfaceWaterExtraction= pcr.scalar(0.)
		#-budget check: setting initial values for cumulative discharge and 
		# net cumulative input, including initial storage [m3]   
		self.totalDischarge= pcr.scalar(0.)
		self.cumulativeDeltaStorage= pcr.catchmenttotal(self.actualStorage,self.LDD)

	def dynamic(self):
		#####################
		# * dynamic section #
		#####################
		#-evaluation of the current date: return current month and the time step used
		#-reading in fluxes over land and water area for current time step [m/d]
		# and read in reservoir demand and surface water extraction [m3]
		try:
			self.landSurfaceQ= clippedRead.get(pcrm.generateNameT(landSurfaceQFileName,self.currentTimeStep()))
		except:
			pass
		try:
			self.potWaterSurfaceQ= clippedRead.get(pcrm.generateNameT(waterSurfaceQFileName,self.currentTimeStep()))
		except:
			pass
		#-surface water extraction and reservoir demand currently set to zero, should
		# be computed automatically and updated to reservoirs
		self.potSurfaceWaterExtraction= pcr.spatial(pcr.scalar(0.))
		#self.waterBodies.demand=  #self.reservoirDemandTSS.assignID(self.waterBodies.ID,self.currentTimeStep(),0.)*self.timeSec
		#-initialization of cumulative values of actual water extractions
		self.actWaterSurfaceQ= pcr.spatial(pcr.scalar(0.))
		self.actSurfaceWaterExtraction= pcr.spatial(pcr.scalar(0.))    
		#-definition of sub-loop for routing scheme - explicit scheme has to satisfy Courant condition
		timeLimit= pcr.cellvalue(pcr.mapminimum((pcr.cover(pcr.ifthen(self.waterBodies.distribution == 0,\
			self.channelLength/self.flowVelocity),\
				self.timeSec/self.nrIterDefault)*self.timeSec/self.nrIterDefault)**0.5),1)[0]
		nrIter= int(self.timeSec/timeLimit)
		nrIter= min(nrIter,int(self.timeSec/300.))
		while float(self.timeSec/nrIter) % 1 <> 0:
			nrIter+= 1
		deltaTime= self.timeSec/nrIter
		#-sub-loop for current time step
		if self.currentDate.day == 1 or nrIter >= 24:
			print '\n*\tprocessing %s, currently using %d substeps of %d seconds\n' % \
				(self.currentDate.date(),nrIter,deltaTime)
		#-update discharge and storage
		for nrICur in range(nrIter):
			#-initializing discharge for the current sub-timestep and fill in values
			# for channels and at outlets of waterbodies
			# * channels *
			estQ= pcr.ifthenelse((self.actualStorage > 0.) & (self.waterBodies.distribution == 0) ,\
				(self.wettedArea/self.alphaQ)**(1./self.betaQ),0.)
			#estQ= pcr.ifthenelse((self.actualStorage > 0.) & (self.waterBodies.distribution == 0) ,\
				#0.5*(self.Q+(self.wettedArea/self.alphaQ)**(1./self.betaQ)),0.)
			#estQ= pcr.min(estQ,self.actualStorage/deltaTime)
			self.report(estQ,'results/qest')
			self.Q= pcr.spatial(pcr.scalar(0.))
			self.Q= pcr.ifthenelse(self.waterBodies.distribution == 0,\
				pcr.kinematic(self.channelLDD,estQ,0.,self.alphaQ,\
					self.betaQ,1,deltaTime,self.channelLength),self.Q)
			# * water bodies *
			self.waterBodies.dischargeUpdate()
			self.Q= self.waterBodies.returnMapValue(self.Q,self.waterBodies.actualQ)
			#-fluxes and resulting change in storage: first the local fluxes are evaluated
			# and aggregated over the water bodies where applicable; this includes the specific runoff [m/day/m2]
			# from input and the estimated extraction from surface water as volume per day [m3/day];
			# specific runoff from the land surface is always positive whereas the fluxes over the water surface
			# are potential, including discharge, and are adjusted to match the availabe storage; to this end,
			# surface water storage and fluxes over water bodies are totalized and assigned to the outlet;
			# discharge is updated in a separate step, after vertical fluxes are compared to the actual storage
			deltaActualStorage= ((self.landFraction*self.landSurfaceQ+\
				self.waterFraction*self.potWaterSurfaceQ)*self.cellArea-\
				self.potSurfaceWaterExtraction)*float(self.duration)/nrIter
			deltaActualStorage= pcr.ifthenelse(self.waterBodies.distribution != 0,\
				pcr.ifthenelse(self.waterBodies.location != 0,\
					pcr.areatotal(deltaActualStorage,self.waterBodies.distribution),0),\
						deltaActualStorage)   
			adjustmentRatio= pcr.ifthenelse(deltaActualStorage < 0.,\
				pcr.min(1.,-self.actualStorage/deltaActualStorage),1.)
			self.actWaterSurfaceQ+= adjustmentRatio*self.potWaterSurfaceQ
			self.actSurfaceWaterExtraction+= adjustmentRatio*self.actSurfaceWaterExtraction
			deltaActualStorage*= adjustmentRatio
			#-local water balance check
			if testLocalWaterBalance:
				differenceActualStorage= self.actualStorage
				differenceActualStorage+= deltaActualStorage
			#-overall water balance check: net input
			self.cumulativeDeltaStorage+= pcr.catchmenttotal(deltaActualStorage,self.LDD)
			#-update storage first with local changes, then balance discharge with storage and update storage
			# with lateral flow and return value to water bodies
			self.actualStorage+= deltaActualStorage
			self.actualStorage= pcr.max(0.,self.actualStorage)
			self.Q= pcr.min(self.Q,self.actualStorage/deltaTime)
			deltaActualStorage= (-self.Q+pcr.upstream(self.LDD,self.Q))*deltaTime
			deltaActualStorage= pcr.ifthenelse(self.waterBodies.distribution != 0,\
				pcr.ifthenelse(self.waterBodies.location != 0,\
					pcr.areatotal(deltaActualStorage,self.waterBodies.distribution),0),\
						deltaActualStorage)
			self.actualStorage+= deltaActualStorage
			self.actualStorage= pcr.max(0.,self.actualStorage)
			self.waterBodies.actualStorage= self.waterBodies.retrieveMapValue(self.actualStorage)
			#-flooded fraction returned
			floodedFraction,floodedDepth,\
					self.wettedArea,self.alphaQ= self.kinAlphaComposite(self.actualStorage,self.floodplainMask)
			self.wettedArea= self.waterBodies.returnMapValue(self.wettedArea,\
				self.waterBodies.channelWidth+2.*self.waterBodies.updateWaterHeight())
			self.waterFraction= pcr.ifthenelse(self.waterBodies.distribution == 0,\
				pcr.max(self.waterFractionMask,floodedFraction),self.waterFractionMask)
			self.landFraction= pcr.max(0.,1.-self.waterFraction)
			self.flowVelocity= pcr.ifthenelse(self.wettedArea > 0,self.Q/self.wettedArea,0.)
			#-local water balance check
			if testLocalWaterBalance:
				differenceActualStorage+= deltaActualStorage
				differenceActualStorage-= self.actualStorage
				totalDifference= pcr.cellvalue(pcr.maptotal(differenceActualStorage),1)[0]
				minimumDifference= pcr.cellvalue(pcr.mapminimum(differenceActualStorage),1)[0]
				maximumDifference= pcr.cellvalue(pcr.mapmaximum(differenceActualStorage),1)[0]
				if abs(totalDifference) > 1.e-3:
					print 'water balance error: total %e; min %e; max %e' %\
						(totalDifference,minimumDifference,maximumDifference)
					if  reportLocalWaterBalance:           
						pcr.report(differenceActualStorage,'mbe_%s.map' % self.currentDate.date())
			#-overall water balance check: updating cumulative discharge and total storage [m3]
			self.totalDischarge+= self.Q*deltaTime
			self.totalStorage= pcr.catchmenttotal(self.actualStorage,self.LDD)
		#-check on occurrence of last day and report mass balance
		if self.currentDate == self.endDate:
			#-report initial maps
			pcr.report(self.Q,self.QIniMap)
			pcr.report(self.actualStorage,self.actualStorageIniMap)
			#-return relative and absolute water balance error per cell and
			# as total at basin outlets
			self.totalDischarge= pcr.ifthen((self.waterBodies.distribution == 0) | \
				(self.waterBodies.location != 0),self.totalDischarge)
			self.cumulativeDeltaStorage= pcr.ifthen((self.waterBodies.distribution == 0) | \
				(self.waterBodies.location != 0),self.cumulativeDeltaStorage)
			massBalanceError= self.totalStorage+self.totalDischarge-\
				self.cumulativeDeltaStorage
			relMassBalanceError= 1.+pcr.ifthenelse(self.cumulativeDeltaStorage <> 0.,
				massBalanceError/self.cumulativeDeltaStorage,0.)
			totalMassBalanceError= pcr.cellvalue(pcr.maptotal(pcr.ifthen(self.basinOutlet,\
				massBalanceError)),1)[0]
			totalCumulativeDeltaStorage= pcr.cellvalue(pcr.maptotal(pcr.ifthen(self.basinOutlet,\
				self.cumulativeDeltaStorage)),1)[0]
			if totalCumulativeDeltaStorage > 0:
				totalRelativeMassBalanceError= 1.+totalMassBalanceError/totalCumulativeDeltaStorage
			else:
				totalRelativeMassBalanceError= 1.
			#-report maps and echo value
			pcr.report(massBalanceError,mbeFileName)
			pcr.report(relMassBalanceError,mbrFileName)
			print '\n*\ttotal global mass balance error [m3]: %8.3g' % totalMassBalanceError
			print '\n*\trelative global mass balance error [-]: %5.3f' % totalRelativeMassBalanceError     
			#-echo to screen: total mass balance error and completion of run
			print '\trun completed'
		#-end of day: return states and fluxes
		#-get surface water attributes?
		if getSurfaceWaterAttributes:
			#-compute the following secondary variables:
			# surface water area [m2]: area given dynamic surface water fraction
			# residence time [days]: volume over discharge, assigned -1 in case discharge is zero
			# surface water depth [m], weighed by channel and floodplain volume
			surfaceWaterArea= self.waterFraction*self.cellArea
			surfaceWaterArea= pcr.ifthenelse(self.waterBodies.distribution != 0,\
					pcr.ifthenelse(self.waterBodies.location != 0,\
						pcr.areatotal(surfaceWaterArea,self.waterBodies.distribution),0),\
							surfaceWaterArea)     
			surfaceWaterResidenceTime= pcr.ifthenelse(self.Q > 0.,self.actualStorage/(self.Q*self.timeSec),-1)
			surfaceWaterDepth= pcr.ifthenelse(self.actualStorage > 0.,\
				pcr.max(0.,self.actualStorage-self.channelStorageCapacity)**2/\
					(self.actualStorage*surfaceWaterArea),0.)
			surfaceWaterDepth+= pcr.ifthenelse(self.actualStorage > 0.,\
				pcr.min(self.channelStorageCapacity,self.actualStorage)**2/(self.waterFractionMask*\
				self.cellArea*self.actualStorage),0.)
			#-reports: values at outlet of lakes or reservoirs are assigned to their full extent
			self.report(pcr.ifthenelse(self.waterBodies.distribution != 0,\
				pcr.areamaximum(surfaceWaterArea,self.waterBodies.distribution),surfaceWaterArea),\
					surfaceWaterAreaFileName)   
			self.report(pcr.ifthenelse(self.waterBodies.distribution != 0,\
				pcr.areamaximum(surfaceWaterResidenceTime,self.waterBodies.distribution),surfaceWaterResidenceTime),\
					surfaceWaterResidenceTimeFileName)
			self.report(pcr.ifthenelse(self.waterBodies.distribution != 0,\
				pcr.areamaximum(surfaceWaterDepth,self.waterBodies.distribution),surfaceWaterDepth),\
					surfaceWaterDepthFileName)
		#-reports on standard output: values at outlet of lakes or reservoirs are assigned to their full extent
		self.report(pcr.ifthenelse(self.waterBodies.distribution != 0,
			pcr.areamaximum(self.flowVelocity,self.waterBodies.distribution),self.flowVelocity),flowVelocityFileName)
		self.report(pcr.ifthenelse(self.waterBodies.distribution != 0,
			pcr.areamaximum(self.Q,self.waterBodies.distribution),self.Q),QFileName)
		self.report(pcr.ifthenelse(self.waterBodies.distribution == 0,\
			floodedFraction,0.),floodedFractionFileName)
		self.report(pcr.ifthenelse(self.waterBodies.distribution == 0,\
			floodedDepth,0.),floodedDepthFileName)
		self.report(self.actualStorage,actualStorageFileName)
		#-update date for time step and report relevant daily output 
		self.currentDate= self.currentDate+datetime.timedelta(self.duration)

		############################################################################
		#-end of model script ######################################################
		############################################################################
