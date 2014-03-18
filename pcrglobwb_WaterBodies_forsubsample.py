################################################################################
#-water bodies: class definition ###############################################
################################################################################

import numpy as np
import pcraster as pcr
from pcraster.numpy import pcr2numpy, numpy2pcr

class pcrglobWaterBodies:
	#-class defining water bodies that interrupt channel flow
	# all attributes are contained by arrays except for the spatial field specifying the distribution:
	# 0) ID linked to 1) distribution: a spatial field of nominal IDs delineating the extent of the respective waterbodies
	# 2) Coordinates: row and column number
	# 3) Type: a nominal classification identifying types (1: lake, 2: reservoir, 3: wetland)
	# 4) Capacity (m3): represents the maximum -operational- capacity of reservoirs, the volume not contributing to
	#      discharge of lakes (death storage)- by default, this value is set to zero for lakes and wetlands
	# 5) Upper limit to operations (-): ratio to maximum storage, set for reservoirs on basis of purpose,
	#      set to unity by default for lakes and wetlands; note that minimum storage is set globally (e.g., 0.1 or 10%)
	#     - capacity and upper limit are read from the text file propertiesTable
	# 6) Parameter relating area to volume, based on actual characteristics, or, in absence thereof, a generic value
	# 7) Channel breadth (m)
	# 8) Average discharge (m3/s)
	# 9) Time-variant total demand (m3)
	#10) Time-variant total storage (m3)
	#11) Time-variant total area (m2)
	#12) Time-variant discharge (m3/s)
	
	def __init__(self,distributionMap,outletMap,typeMap,channelBreadthMap,averageQMap,bankfulQMap,LDDMap,parameterTBL,deltaTime,clippedRead):
		#-clippedRead
		self.clippedRead= clippedRead
		#-constants
		self.deltaTime= deltaTime
		self.MV= -999.9
		self.cLake= 1.7
		self.minLimit= 0.0
		#-spatial field of nominal IDs delineating the extent of the respective waterbodies
		self.distribution= self.clippedRead.get(distributionMap,'nominal')
		self.outlet= self.clippedRead.get(outletMap,'nominal')
		waterBodiesType= self.clippedRead.get(typeMap,'nominal')
		LDD= self.clippedRead.get(LDDMap,'ldd')
		endorheicLakes= pcr.ifthen((pcr.areatotal(pcr.scalar(self.outlet != 0),self.distribution) == 0) & (self.distribution != 0),\
			self.distribution)
		self.location= pcr.cover(pcr.ifthen(self.outlet != 0,self.outlet),\
			pcr.ifthenelse((self.distribution != 0) & (LDD == 5),self.distribution,0))
		pcr.report(self.location,'maps/waterbodies_reportlocations.map')
		#-extract ID and location as row and column number for outlets; these are subsequently used to extract type from other maps
		tempIDArray= pcr2numpy(self.location,self.MV)
		nrRows,nrCols= tempIDArray.shape
		iCnt= 0
		self.ID= np.ones((1))*self.MV
		self.coordinates= np.ones((2))*self.MV
		for row in xrange(nrRows):
			for col in xrange(nrCols):
				if tempIDArray[row,col] > 0:
					iCnt+= 1
					if iCnt > 1:
						self.ID= np.append(self.ID,tempIDArray[row,col])
						self.coordinates= np.vstack((self.coordinates,\
							np.array([row,col])))     
					else:
						self.ID= np.array([tempIDArray[row,col]])
						self.coordinates= np.array([row,col])
		#-process valid entries
		if self.ID[0] <> self.MV:
			#-reset nrRows to nr of rowwise entries in ID
			self.nrEntries= self.ID.shape[0]
			#-sort on ID
			indices= self.ID.argsort()
			self.ID= self.ID[indices]
			self.coordinates= self.coordinates[indices]
		else:
			self.nrEntries= None
			self.coordinates= np.array([])
		#-read from maps: type, channel breadth and average discharge
		self.type= self.retrieveMapValue(self.clippedRead.get(typeMap,'nominal'))
		self.channelWidth= self.retrieveMapValue(self.clippedRead.get(channelBreadthMap))
		self.averageQ= self.retrieveMapValue(self.clippedRead.get(averageQMap))
		self.bankfulQ= self.retrieveMapValue(self.clippedRead.get(bankfulQMap))
		self.endorheic= self.retrieveMapValue(pcr.cover(endorheicLakes,0))
		#-set to zero, to be updated in script
		self.demand= np.zeros((self.nrEntries))
		self.actualStorage= np.zeros((self.nrEntries))
		self.actualArea= np.zeros((self.nrEntries))
		self.actualQ= np.zeros((self.nrEntries))
		#-read from table
		self.capacity= np.ones((self.nrEntries))*self.MV
		self.maxLimit= np.ones((self.nrEntries))*self.MV
		self.avParameter= np.ones((self.nrEntries))*self.MV
		tempIDArray= np.loadtxt(parameterTBL)
		if self.nrEntries <> None:
			for iCnt in xrange(self.nrEntries):
				ID= self.ID[iCnt]
				mask= tempIDArray[:,0] == ID
				if np.any(mask):
					#-entry in table with reservoir properties
					#-set capacity, fractional maximum storage limit and area-volume parameter
					self.capacity[iCnt]= tempIDArray[:,1][mask]
					self.maxLimit[iCnt]= tempIDArray[:,2][mask]
					self.avParameter[iCnt]= tempIDArray[:,3][mask]
				else:
					#-lake or wetland: set capacity and upper limit to zero and the area-volume parameter
					# to default values
					self.capacity[iCnt]= 0.
					self.maxLimit[iCnt]= 0.
					if self.type[iCnt] == 1:
						#-lake
						self.avParameter[iCnt]= 210.5
					else:
						#-wetland
						self.avParameter[iCnt]= 1407.2

	def retrieveMapValue(self,pcrX):
		#-retrieves values from a map and returns an array conform the IDs stored in properties
		if self.nrEntries <> None:
			x= np.ones((self.nrEntries))*self.MV
			tmpIDArray= pcr2numpy(pcrX,self.MV)
			for iCnt in xrange(self.nrEntries):
				row,col= self.coordinates[iCnt,:]
				x[iCnt]= tmpIDArray[row,col]
			return x
		else:
			return np.zeros(self.ID.shape)

	def returnMapValue(self,pcrX,x):
		#-retrieves value from an array and update values in the map
		if self.nrEntries <> None:
			tempIDArray= pcr2numpy(pcrX,self.MV)
			for iCnt in xrange(self.nrEntries):
				row,col= self.coordinates[iCnt,:]
				tempIDArray[row,col]= x[iCnt]
			pcrX= numpy2pcr(pcr.Scalar,tempIDArray,self.MV)
		return pcrX

	def updateWaterHeight(self):
		#-computes water height and inserts this into the map
		if self.nrEntries <> None:
			waterHeight= np.zeros((self.nrEntries))
			mask= self.type != 2
			waterHeight[mask]= np.maximum(0.,self.actualStorage[mask]-self.capacity[mask])/\
				np.maximum(1.,self.actualArea[mask])
			mask= self.type == 2
			waterHeight[mask]= self.actualStorage[mask]/np.maximum(1.,self.actualArea[mask])
		else:
			waterHeight= np.zeros(self.ID.shape)
		return waterHeight
	
	def dischargeUpdate(self):
		#-updates the discharge taking the actual storage into account
		#-water height is a temporal variable, giving the level contributing to discharge
		waterHeight= self.updateWaterHeight()
		self.actualQ= np.zeros(self.ID.shape)
		if self.nrEntries <> None:
			#-lakes and reservoirs: discharge depends on water height: cLake*WH**1.5*B
			mask= (self.type != 2) & (self.endorheic == 0)
			self.actualQ[mask]= self.cLake*waterHeight[mask]**1.5*self.channelWidth[mask]
			#-reservoirs: release is a function of demand, storage and limits
			mask= self.type == 2
			#-Q as function of storage and average discharge
			self.actualQ[mask]= np.minimum(1.,np.maximum(0.,\
				self.actualStorage[mask]-self.minLimit*self.capacity[mask])/\
				((self.maxLimit[mask]-self.minLimit)*self.capacity[mask]))*self.averageQ[mask]
			##-Q as a function of demand and available storage
			#self.actualQ[mask]= np.maximum(self.actualQ[mask],\
				#np.minimum(self.demand[mask]*np.minimum(1.,self.actualStorage[mask]/\
				#(self.minLimit*self.capacity[mask])),self.actualStorage[mask])/self.deltaTime)
			#-update reservoir release to account for storage in excess of the upper limit
			# or above total storage capacity
			self.actualQ[mask]+= np.maximum(0.,self.actualStorage[mask]-self.capacity[mask])/self.deltaTime+\
				np.maximum(0,self.actualStorage[mask]-self.maxLimit[mask]*self.capacity[mask])/\
					((1.-self.maxLimit[mask])*self.capacity[mask])*\
						np.maximum(0.,self.bankfulQ[mask]-self.actualQ[mask])
			#-limit discharge to available storage
			self.actualQ= np.minimum(self.actualStorage/self.deltaTime,self.actualQ)

	def updateRoutingParameters(self,actualStorage,characteristicDistance):
		#-passes the actual storage available for routing and the corresponding travel distance
		# for water bodies to the corresponding map
		#-water height is a temporal variable, giving the level contributing to discharge
		#-shape
		nrRows= self.ID.shape[0]
		#-water height is a temporal variable, giving the level contributing to discharge
		waterHeight= self.updateWaterHeight()
		#-assigns the routable storage to the lake outlet and update the characteristic distance
		#-the fraction of storage that is routed is the ratio of the actual flow over available storage
		fractionRouted= np.ones((nrRows))
		fractionRouted= np.minimum(1.,self.actualQ*self.deltaTime/np.maximum(1.e-12,self.actualStorage))
		actualStorage=  self.returnMapValue(actualStorage,fractionRouted*self.actualStorage)
		#-update remaining storage 
		self.actualStorage*= np.maximum(0.,1.-fractionRouted)
		#-characteristic travel distance
		characteristicDistance= self.returnMapValue(characteristicDistance,\
			self.actualQ*self.deltaTime/np.maximum(1.e-3,waterHeight*self.channelWidth))
		#-return updated maps
		return actualStorage,characteristicDistance

	# /*end of class defining water bodies*/
