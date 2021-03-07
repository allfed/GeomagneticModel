import Params
import fits

import glob
import Model.TFsite
from Model.TFsite import TFsite
import Model.MTsite
from Model.MTsite import MTsite

import Model.GCmodel
from Model.GCmodel import GCmodel

import h5py
import spacepy.coordinates as coord
from spacepy.time import Ticktock
import matplotlib.pyplot as plt
import itertools
import numpy as np

class EarthModel:

	# initialize a site by providing a list of frequencies w to determine the transfer function.
	# if a path to a TF site EDI document is provided, the init function 
	def __init__(self):
		Params.importIfNotAlready()
		self.refratesperyear=[]
		self.refEfields=[]
		self.refwindowperiods=[]
		self.combinedslopes=[]
		self.combinedexponents=[]
		self.refconductivity=Params.refconductivity
		self.refmaglat=Params.refmaglat
		self.apparentCondMap=[]
		self.windowedEmaps=[]

	def loadApparentCond(self):
		self.apparentCondMap=np.load(Params.apparentcondloc)

	def initTFsites(self):
		tfsites=[]
		#set up default frequencies for TF sites
		f=[1.367187E-01,1.093750E-01,8.593753E-02,6.640627E-02,5.078124E-02,3.906250E-02 ,\
		3.027344E-02,2.343750E-02,1.855469E-02,1.464844E-02,1.171875E-02,9.765625E-03 ,\
		7.568361E-03,5.859374E-03,4.638671E-03,3.662109E-03,2.929688E-03,2.441406E-03 ,\
		1.892090E-03,1.464844E-03,1.159668E-03,9.155271E-04,7.324221E-04,6.103516E-04 ,\
		4.425049E-04,3.204346E-04,2.136230E-04,1.373291E-04,8.392331E-05,5.340577E-05]
		w = [x*2.0*np.pi for x in f] #rad s^-1 (angular frequency)
		rows=[0,0,0,0,0,0,\
		1,1,1,1,1,1,\
		2,2,2,2,2,2,\
		3,3,3,3,3,3,\
		4,4,4,4,4,4]
		cols=[0,1,2,3,4,5,\
		0,1,2,3,4,5,\
		0,1,2,3,4,5,\
		0,1,2,3,4,5,\
		0,1,2,3,4,5]
		for i in range(0,len(Params.tfsitenames)):
			folder=Params.tfsitesdir
			filename=Params.tfsitenames[i]
			tfsites=tfsites + [TFsite(folder+filename+'.edi',f,rows,cols)]
		return tfsites


	def initMTsites(self):
		mtsites=[]
		if(not Params.mtsitenames): #if mt sites were not supplied as params, just use all existing mtsites in the mtsites directory as default data
			loaddirectory=Params.mtRepeatRatesDir
			allfiles=glob.glob(loaddirectory+'*.np[yz]')
			for i in range(0,len(allfiles)):
				mtsites = mtsites + [MTsite('',i)]
		else:
			for i in range(0,len(Params.mtsitenames)):
				folder=Params.mtsitesdir
				filename=Params.mtsitenames[i]
				mtsites=mtsites + [MTsite(folder+filename+'.netcdf',i)]

		return mtsites

	def initGCmodel(self):
		return GCmodel()

	#the purpose of this function is to process the MT sites while calculating and saving E fields at all durations
	def processMTsites(self,mtsites,tfsites):
		for i in range(0,len(tfsites)):
			mtsite=mtsites[i]
			tfsite=tfsites[i]
			mtsite.importSite()
			for i in range(0,mtsite.nchunks):
				mtsite.createChunk(i)
				mtsite.cleanChunkBfields(i)
				print('')
				mtsite.calcChunkEfields(tfsite,i)
				print('')
				mtsite.saveChunkEfields(i)
			# mtsite.calcEratesPerYear()
			# mtste.saveEratesPerYear()

	def calcandplotEratesPerYear(self,mtsites):
		for i in range(0,len(mtsites)):
			mtsite=mtsites[i]
			mtsite.importSite()
			mtsite.createChunks()
			mtsite.loadEfields()
			mtsite.calcEratesPerYear()
			mtsite.fitEratesPerYear()
			mtsite.plotEratesPerYear()



	def calcAndPlotMTEfields(self,tfsites,mtsites):
		for i in range(0,len(mtsites)):
			mtsite=mtsites[i]
			tfsite=tfsites[i]
			mtsite.importSite()
			mtsite.createChunks()
			chunkindex=1
			mtsite.cleanChunkBfields(chunkindex)
			mtsite.calcChunkEfields(tfsite,chunkindex)
			#for love, 2019 vaq55,FRD plot
			# mtsite.plotEfields(chunkindex,732000-5*60-20+12*60,732000-5*60-20+84*60)
			
			#for love, 2018 vaq58,FRD plot
			mtsite.plotEfields(chunkindex,732000-5*60-20+12*60,732000-5*60-20+108*60)

	def loadPreviousMTfits(self,mtsites):
		for i in range(0,len(mtsites)):
			mtsites[i].loadWindowedRates()
			mtsites[i].fitEratesPerYear()

	# see https://stackoverflow.com/questions/7948450/conversion-from-geographic-to-geomagnetic-coordinates
	def geotomag(self,lat,lon):
		#call with altitude at sea level (6371.0 km) and lat/lon in degrees 
		cvals = coord.Coords([1, np.float(lat), np.float(lon)], 'GEO', 'sph',['Re','deg','deg'])
		#set time epoch for coordinates:
		cvals.ticks=Ticktock(['2021-02-01T12:00:00'], 'ISO')
		#return the magnetic coords in the same units as the geographic:
		return cvals.convert('MAG','sph')

	# magic formula to get the divisor for the field to go from estimate at 72 to at any magnetic latitude 
	def getMagLatDivisor(self,ml):
		divisor = (0.211 \
		+ 0.000681*ml \
		+ -0.000148*(ml**2) \
		+ -0.00000226*(ml**3) \
		+ -0.0000000347*(ml**4) \
		+ 0.00000000087*(ml**5) \
		+ 0.00000000011*(ml**6) \
		+ -3.49E-14*(ml**7) \
		+ -2.54E-14*(ml**8) \
		+ -7.53E-18*(ml**9) \
		+ 1.61E-18*(ml**10))/1.66675
		return divisor

	# find E fields at all locations on conductivity map for each duration
	# Correct for geomagnetic latitude and apparent conductivity
	def calcGlobalEfields(self, ratePerYears):
		self.windowedEmaps=[]
		for r in ratePerYears:
			windowedEmapsatrate=[]
			#save a map of e fields for each duration at this rate per Year
			for i in range(0,len(self.refwindowperiods)):
				windowperiod=self.refwindowperiods[i]
				print(r)
				print(self.combinedslopes[i])
				print(self.combinedexponents[i])
				refField=fits.powerlawxfromy(r,self.combinedslopes[i],self.combinedexponents[i])
				print('refField level rest of map is proportional to this field (V/km):'+str(refField))
				col=0
				longInterval=int(np.floor(Params.longituderes/0.25))
				latInterval=int(np.floor(Params.latituderes/0.25))
				with h5py.File(Params.globalcondloc, 'r') as f:
					latitudes = f.get('LATGRID')[0]
					longitudes = f.get('LONGRID')[:,0]
					lenlat=len(latitudes)
					lenlong=len(longitudes)
					EArr = np.zeros((int(lenlong/longInterval), int(lenlat/latInterval)))

					minE=10^10
					maxE=-10^10
					for longkey in range(0,lenlong-1,longInterval):
						row=0
						
						for latkey in range(0,lenlat-1,latInterval):
							latitude = f['LATGRID'][longkey,latkey]
							longitude=f['LONGRID'][longkey,latkey]
							c=self.apparentCondMap[row,col]

							magcoords = self.geotomag(latitude, longitude)
							maglat=magcoords.data[0][1]
							maglong=magcoords.data[0][2]

							# equivalentML=getEquivalentMagLat(maglat,ratePerYear)
							mld = self.getMagLatDivisor(maglat)

							E = refField*np.sqrt(self.refconductivity)/np.sqrt(abs(c))*mld

							if maxE < E:
								maxE = E

							if minE > E:
								minE = E

							EArr[row][col]=E

							row = row+1	
						
						col = col+1
				plt.figure(1)
				plt.imshow(np.flipud(EArr[8:,:]), cmap='hot', interpolation='nearest')
				plt.title('E field levels, '+str(r)+' per year, '+str(windowperiod)+'s')
				cbar = plt.colorbar()
				cbar.set_label('E field (V/km)')
				plt.savefig(Params.globalEfieldPlots+'Results'+str(r)+'perYearWindow'+str(windowperiod)+'s.png')
				if(i==0):
					plt.show()
				windowedEmapsatrate = windowedEmapsatrate + [windowperiod,EArr]
			self.windowedEmaps = self.windowedEmaps + [r,windowedEmapsatrate]

		np.save(Params.globalEfieldData,self.windowedEmaps)

	#adjusts ratesperyear to reference maglat and conductivity and saves the result as properties of this EarthModel instance
	def calcReferenceRateperyear(self,tfsites,mtsites):
		self.refratesperyear=[]
		self.refEfields=[]
		self.refwindowperiods=[]

		allEatwindow=[]
		allratesatwindow=[]
		allwindows=[]
		for i in range(0,len(mtsites)):
			tfsite=tfsites[i]
			mtsite=mtsites[i]

			for j in range(0,int(np.floor(len(mtsite.windowedRates))/3)):
				durationindex = j*3
				efieldindex = j*3+1
				ratesindex = j*3+2

				windowperiod = mtsite.windowedRates[durationindex]
				Efields = mtsite.windowedRates[efieldindex]
				rates = mtsite.windowedRates[ratesindex]
				windowbucket=-1
				apparentcond=tfsite.getApparentcCloseToWindowPeriod(windowperiod)
				adjustedtosite=self.adjustEfieldsToRefCond(apparentcond,Efields,windowperiod)
				adjustedtomaglat=np.real(np.array(adjustedtosite)/self.getMagLatDivisor(mtsite.maglat))

				for i in range(0,len(allwindows)):
					if(windowperiod==allwindows[i]):
						windowbucket=i
						allratesatwindow = allratesatwindow + rates
						allEatwindow = allEatwindow+adjustedtomaglat
						break
				if(windowbucket==-1): 
					allratesatwindow = allratesatwindow+[rates]
					allwindows = allwindows+[windowperiod]
					allEatwindow = allEatwindow+[adjustedtomaglat]

		self.refwindowperiods=allwindows
		self.refratesperyear=allratesatwindow
		self.refEfields=allEatwindow

	#adjust mtsite fits to reference conductivity.
	#this assumes they all have
	def adjustEfieldsToRefCond(self,apparentcond,Efields,windowperiod):
		#find the closest TFsite frequency match to the given window period.
		#we want the apparent conductivity at the average frequency of the EMP, which can be approximated by 1/windowperiod

		return np.sqrt(apparentcond)/np.sqrt(self.refconductivity)*np.array(Efields)

	def plotCombinedRates(self):
		plt.figure()
		plt.clf()
		lines=['o', 'v', '^', '<', '>', 's', '8', 'p']*100
		plt.loglog()
		ax = plt.gca()
		self.combinedslopes=[]
		self.combinedexponents=[]
		for i in range(0,len(self.refwindowperiods)):
			windowperiod = self.refwindowperiods[i]
			rates  = self.refratesperyear[i]
			Efields= self.refEfields[i]
			[slope,exponent]=fits.fitPower(Efields,rates)
			color = next(ax._get_lines.prop_cycler)['color']
			plt.plot(Efields,rates, linestyle='', markeredgecolor='none', marker=lines[i], color=color)
			plt.plot(Efields,fits.powerlaw(Efields,slope,exponent), linestyle='-', color = color, lw=1,label = "Powerfit, field averaged over "+str(windowperiod)+" seconds")
			self.combinedslopes=self.combinedslopes+[slope]
			self.combinedexponents=self.combinedexponents+[exponent]
			print('plottingcombinedlopes')
			print(self.combinedslopes)
			print('plottingcombinedexps')
			print(self.combinedexponents)
			# plt.plot(Efields,rates, lw=1,label = "Field averaged over "+str(windowperiod)+" seconds")

		plt.legend()
		plt.title('Rate geoelectric field is above threshold')
		plt.xlabel('Geoelectric Field (V/km)')
		plt.ylabel('Average rate per year distinct contiguous sample average is above E field (counts/year)')
	
		plt.show()
