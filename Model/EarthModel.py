import Params
import fits

import matplotlib.colors as colors
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
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import contextily as ctx
class EarthModel:

	# initialize a site by providing a list of frequencies w to determine the transfer function.
	# if a path to a TF site EDI document is provided, the init function 
	def __init__(self):
		Params.importIfNotAlready()
		self.refcumsum=[]
		self.refEfields=[]
		self.combinedpowerfits=[]
		self.combinedlogfits=[]
		self.refconductivity=Params.refconductivity
		self.refmaglat=Params.refmaglat
		self.apparentCondMap=[]
		self.windowedEmaps=[]
		self.allwindowperiods=[]
		self.averagedRatios=[]
		self.GClatitudes=[]
		self.GClongitudes=[]

	def loadApparentCond(self):
		self.apparentCondMap=np.load(Params.apparentcondloc)
		self.f=1/120

	def initTFsites(self):
		tfsites=[]
		#set up default frequencies for TF sites
		# f=[1.367187E-01,1.093750E-01,8.593753E-02,6.640627E-02,5.078124E-02,3.906250E-02 ,\
		# 3.027344E-02,2.343750E-02,1.855469E-02,1.464844E-02,1.171875E-02,9.765625E-03 ,\
		# 7.568361E-03,5.859374E-03,4.638671E-03,3.662109E-03,2.929688E-03,2.441406E-03 ,\
		# 1.892090E-03,1.464844E-03,1.159668E-03,9.155271E-04,7.324221E-04,6.103516E-04 ,\
		# 4.425049E-04,3.204346E-04,2.136230E-04,1.373291E-04,8.392331E-05,5.340577E-05]
		# w = [x*2.0*np.pi for x in f] #rad s^-1 (angular frequency)
		# rows=[0,0,0,0,0,0,\
		# 1,1,1,1,1,1,\
		# 2,2,2,2,2,2,\
		# 3,3,3,3,3,3,\
		# 4,4,4,4,4,4]
		# cols=[0,1,2,3,4,5,\
		# 0,1,2,3,4,5,\
		# 0,1,2,3,4,5,\
		# 0,1,2,3,4,5,\
		# 0,1,2,3,4,5]
		for i in range(0,len(Params.tfsitenames)):
			tfsite=TFsite()
			tfsite.initWithID(i)
			tfsites=tfsites + [tfsite]
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
				if(Params.useMTsite[i]):
					mtsites=mtsites + [MTsite(folder+filename+'.netcdf',i)]
				else:
					mtsites=mtsites + [[]]
		return mtsites

	def initGCmodel(self):
		return GCmodel()

	#the purpose of this function is to process the MT sites while calculating and saving E fields at all durations
	def calcAndSaveEfields(self,mtsites,tfsites):
		for i in range(0,len(tfsites)):
			if(not Params.useMTsite[i]):
				continue
			mtsite=mtsites[i]
			tfsite=tfsites[i]
			mtsite.importSite()
			mtsite.createChunks()
			for j in range(0,mtsite.nchunks):
				mtsite.cleanChunkBfields(j)
				print('')
				mtsite.calcChunkEfields(tfsite,j)
				print('')
				print('saving chunk '+str(j))
				mtsite.saveChunkEfields(j)

	def calcAndSaveRecurrence(self,mtsites,tfsites):
		for i in range(0,len(tfsites)):
			if(not Params.useMTsite[i]):
				continue
			mtsite=mtsites[i]
			tfsite=tfsites[i]
			mtsite.importSite()
			mtsite.createChunks()
			mtsite.loadEfields()
			mtsite.calcEratesPerYear(False)
			mtsite.saveEratesPerYear()

	#the purpose of this function is to process the MT sites while calculating and saving E fields at all durations
	def processMTsites(self,mtsites,tfsites):
		for i in range(0,len(tfsites)):
			if(not Params.useMTsite[i]):
				continue
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

			mtsite.calcEratesPerYear(False)
			mtsite.saveEratesPerYear()

	def peakEvsDuration(self,mtsites,plot):
		if(plot):
			plt.figure()
		durationEratios=[]
		self.allwindowperiods=[]
		self.cumulativeyears=[]
		for k in range(0,len(mtsites)):
			if(not Params.useMTsite[k]):
				continue
			mtsite=mtsites[k]
			mtsite.importSite()
			mtsite.createChunks()
			mtsite.loadWindowedCounts()
			mtsite.fitEratesPerYear()
			averageOccurrence=mtsite.calcPeakEvsDuration(plot)

			#add this duration data to the existing windowedCounts
			for i in range(0,int(np.floor(len(mtsite.windowedCounts))/5)):
				durationindex = i*5
				efieldindex = i*5+1
				countsindex = i*5+2
				cumsumindex = i*5+3
				netyearsindex = i*5+4

				duration = mtsite.windowedCounts[durationindex]
				Efields = mtsite.windowedCounts[efieldindex]
				counts = mtsite.windowedCounts[countsindex]
				cumsum = mtsite.windowedCounts[cumsumindex]
				netyears = mtsite.windowedCounts[netyearsindex]
				rpy=cumsum/mtsite.cumulativeyears

				matchid=-1
				for j in range(0,len(self.allwindowperiods)):
					wp=self.allwindowperiods[j]
					if(wp==duration):
						matchid=j
						durationEratios[j].append(averageOccurrence[i])
				if(matchid==-1):
					self.allwindowperiods.append(duration)
					durationEratios.append([averageOccurrence[i]])
					self.cumulativeyears.append(netyears)

	
		#average the occurrence with all the other windows
		self.averagedRatios=[]
		for ratios in durationEratios:
			self.averagedRatios.append(np.mean(ratios))

		if(plot):
			# plt.show()
			# plt.figure()
			plt.plot(self.allwindowperiods,self.averagedRatios,lw=7,label = " all ratios averaged")
			plt.legend()
			plt.title('Once per decade ratio of durations')
			plt.xlabel('Duration (s)')
			plt.ylabel('Ratio of peak field to field at 60s')

			plt.show()

	def calcandplotEratesPerYear(self,mtsites,fittype):
		plt.figure()
		for i in range(0,len(mtsites)):
			if( not Params.useMTsite[i]):
				continue
			mtsite=mtsites[i]
			mtsite.importSite()
			mtsite.createChunks()

			#recalcs
			# mtsite.loadEfields()
			# mtsite.calcEratesPerYear()
			# mtsite.saveEratesPerYear()
			# continue

			#reload
			mtsite.loadWindowedCounts()
			
			# fitx=np.array(mtsite.windowedCounts[1])

			# # plt.figure()
			# # plt.loglog(fitx,mtsite.windowedCounts[2])
			# # plt.show()
			# # mtsite.loadWindowedCounts()
			# # print(np.array(mtsite.windowedCounts[1]))
			# # plt.figure()

			# # offsets=range(-100,100,100)
			# # widths=range(1,10,5)
			# # for offset in offsets:
			# # 	for width in widths:
			# # 		fity=fits.importedlognormal(np.array(mtsite.windowedCounts[1]),-57+offset,(18.4*width)**2)
			# # 		plt.loglog(fitx,fity,lw=1,label = "fit imp offset:"+str(offset)+" width"+str(width))
			# # 		fity=1.1*fits.lognormal(np.array(mtsite.windowedCounts[1]),-57+offset,(18.4*width)**2)
			
			# params=fits.fitLognormalwithloc(np.array(mtsite.windowedCounts[1]),np.array(mtsite.windowedCounts[2]))

			# # start=fits.lognormal(np.array(mtsite.windowedCounts[1]),-11,4.6)
			
			# upsilon=params[0]
			# epsilon=params[1]
			# loc=params[2]
			# fity=fits.locimportedlognormal(np.array(mtsite.windowedCounts[1]),upsilon,np.abs(epsilon),loc)

			# # plt.loglog(fitx,start,'-')
			# plt.loglog(fitx,fity,lw=1,label = "fit upsilon:"+str(upsilon)+" epsilon"+str(epsilon)+'loc'+str(loc))
			# # plt.loglog(fitx,fity,lw=1,label = "fit offset:"+str(offset)+" width"+str(width))

			# plt.legend()
			# plt.loglog(fitx,mtsite.windowedCounts[2])
			# plt.show()
			# print('mtsite.windowedCounts')
			# print(mtsite.windowedCounts)
			mtsite.fitEratesPerYear()
			mtsite.plotEratesPerYear(fittype)

		plt.legend()
		plt.title('Rate geoelectric field is above threshold')
		plt.xlabel('Geoelectric Field (V/km)')
		plt.ylabel('Average rate per year distinct contiguous sample average is above E field (counts/year)')
		plt.show()
	
	def Estats(self,mtsites):
		for i in range(0,len(mtsites)):
			if( not Params.useMTsite[i]):
				continue
			mtsites[i].importSite()
			mtsites[i].createChunks()
			mtsites[i].loadEfields()
			mtsites[i].Estats()


	def plot1989Efields(self,mtsites):
		for i in range(0,len(mtsites)):
			if( not Params.useMTsite[i]):
				continue
			mtsite=mtsites[i]
			mtsite.importSite()
			mtsite.createChunks()
			mtsite.loadEfields()
			#for love, 2018 vaq58,FRD plot
			mtsite.plotEfields(i,732000-5*60-20+12*60,732000-5*60-20+108*60)

	def loadPreviousMTfits(self,mtsites):
		for i in range(0,len(mtsites)):
			if(not Params.useMTsite[i]):
				continue
			mtsites[i].importSite()
			mtsites[i].loadWindowedCounts()
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
			for i in range(0,len(self.allwindowperiods)):
				windowperiod=self.allwindowperiods[i]
				mean=self.combinedlogfits[i*5+1]
				std=self.combinedlogfits[i*5+2]
				loc=self.combinedlogfits[i*5+3]
				ratio=self.combinedlogfits[i*5+4]
				exponent=self.combinedpowerfits[i*5+1]
				refField=fits.powerlawxfromy(r,exponent)
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

	#adjusts ratesperyear for all the sites to the reference maglat and conductivity and saves the result as properties of this EarthModel instance
	def calcReferenceRateperyear(self,tfsites,mtsites,fittype,plot):
		self.refcumsum=[]
		self.refEfields=[]
		allcumsum=[]
		allE=[]
		allYears=[]
		allCounts=[]
		if(plot):
			plt.figure()
		for i in range(0,len(mtsites)):
			if(not Params.useMTsite[i]):
				continue
			tfsite=tfsites[i]
			mtsite=mtsites[i]
			plt.loglog()
			allEatwindow=[]
			allcumsumatwindow=[]
			allcountsatwindow=[]
			allyearsatwindow=[]

			for j in range(0,int(np.floor(len(mtsite.windowedCounts))/5)):
				durationindex = j*5
				efieldindex = j*5+1
				countsindex = j*5+2
				cumsumindex = j*5+3
				netyearsindex = j*5+4

				duration = mtsite.windowedCounts[durationindex]
				Efields = mtsite.windowedCounts[efieldindex]
				counts = mtsite.windowedCounts[countsindex]
				cumsum = mtsite.windowedCounts[cumsumindex]
				years = mtsite.windowedCounts[netyearsindex]
				rpy=cumsum/mtsite.cumulativeyears

				Efields = mtsite.windowedCounts[efieldindex]
				counts = mtsite.windowedCounts[countsindex]
				matchid=-1
				for k in range(0,len(self.allwindowperiods)):
					wp=self.allwindowperiods[k]
					if(wp==duration):
						if(matchid!=-1):
							print('error: there appear to be duplicate windows for a single MTsite!')
							quit()
						matchid=k
						adjustedEfields=self.averagedRatios[matchid]*Efields

				if(matchid==-1):
					print('error: no known adjustment for this windowperiod')
					print('be sure you\'ve calculated the E field vs duration already.')
					quit()

				apparentcond=tfsite.getApparentcCloseToWindowPeriod(duration)
				adjustedtosite=self.adjustEfieldsToRefCond(apparentcond,adjustedEfields,duration)
				adjustedtomaglat=np.real(np.array(adjustedtosite)/self.getMagLatDivisor(mtsite.maglat))
				allEatwindow.append(np.array(adjustedtomaglat))
				allcumsumatwindow.append(np.array(cumsum))
				allcountsatwindow.append(np.array(counts))
				allyearsatwindow.append(years)

			allcumsum.append(allcumsumatwindow)
			allE.append(allEatwindow)
			allCounts.append(allcountsatwindow)
			allYears.append(allyearsatwindow)

			if(plot):	
				duration = mtsite.windowedCounts[0]
				# Efields = mtsite.windowedCounts[1]

				counts = mtsite.windowedCounts[2]
				cumsum = mtsite.windowedCounts[3]
				rpy=cumsum/mtsite.cumulativeyears
				probtoRPYratio=np.max(cumsum)/mtsite.cumulativeyears

				[exponent]=fits.fitPower(allEatwindow[0],cumsum/np.max(cumsum))
				print('site '+str(mtsite.MTsitefn)+' window: '+str(duration))
				#use PDF to determine mean and standard deviation of underlying normal distribution
				[guessMean,guessStd]=fits.getGuesses(allEatwindow[0],counts,False)

				#fit to the datapoints (CDF has a probability of 1 at the first datapoint)
				[mean,std,loc]=fits.fitLognormalCDF(allEatwindow[0],cumsum/np.max(cumsum),guessMean,guessStd,False)

				fitplotted=False
				#if we've calculated power fits for all the windows
				plt.plot(allEatwindow[0],rpy,lw=1,label = "Efields averaged over "+str(duration)+' seconds, '+str(mtsite.sitename))
				if('power' in fittype or fittype=='all'):
					fitplotted=True
					# print('power coeffs')
					# print([exponent,probtoRPYratio])

					plt.plot(allEatwindow[0], fits.powerlaw(allEatwindow[0],exponent)*probtoRPYratio, lw=1,label = "Field averaged over "+str(duration)+" seconds, powerfit")

				if('lognormal' in fittype or fittype == 'all'):
					fitplotted=True

					# print('lognormal coeffs')
					# print([mean,np.abs(std),loc,probtoRPYratio])

					yfit=probtoRPYratio*fits.logcdf(np.array(allEatwindow[0]),mean,np.abs(std),loc)

					boundedfit=np.array(yfit[np.array(yfit)>10**-4])
					boundedEfields=np.array(allEatwindow[0])[np.array(yfit)>10**-4]
					plt.plot(boundedEfields,boundedfit,lw=1,label = "Field averaged over "+str(duration)+" seconds, lognormalfit")	

		if(plot):		
			plt.legend()
			plt.title('Rate geoelectric field is above threshold')
			plt.xlabel('Geoelectric Field (V/km)')
			plt.ylabel('Average rate per year distinct contiguous sample average is above E field (counts/year)')

			plt.show()

		self.refcumsum=allcumsum
		self.refEfields=allE
		self.refCounts=allCounts
		self.refYears=allYears

	#adjust mtsite fits to reference conductivity.
	#this assumes they all have
	def adjustEfieldsToRefCond(self,apparentcond,Efields,windowperiod):
		#find the closest TFsite frequency match to the given window period.
		#we want the apparent conductivity at the average frequency of the EMP, which can be approximated by 1/windowperiod

		return np.sqrt(apparentcond)/np.sqrt(self.refconductivity)*np.array(Efields)

	def calcCombinedRates(self,plot):
		self.combinedexponents=[]
		#combine the fields into one distribution
		for j in range(0,len(self.allwindowperiods)):
			windowperiod = self.allwindowperiods[j]
			combinedEfields=np.array([])
			combinedrates=np.array([])
			combinedcumsum=np.array([])
			combinedcounts=np.array([])
			combinedyears=0

			for i in range(0,len(self.refYears)):
				if(not Params.useMTsite[i]):
					continue
				combinedyears=combinedyears+self.refYears[i][j]
				combinedEfields=np.append(combinedEfields,self.refEfields[i][j])
				ratestmp=np.array(self.refcumsum[i][j])/self.refYears[i][j]
				combinedrates=np.append(combinedrates,ratestmp)
				combinedcumsum=np.append(combinedcumsum,self.refcumsum[i][j])
				combinedcounts=np.append(combinedcounts,self.refCounts[i][j])
			probtoRPYratio=np.max(combinedcumsum)/(combinedyears/len(self.refYears))
			
			#sort by descending E field
			counts = [x for _,x in sorted(zip(-np.array(combinedEfields),combinedcounts))]
			cumsum = [x for _,x in sorted(zip(-np.array(combinedEfields),combinedcumsum))]
			rates = [x for _,x in sorted(zip(-np.array(combinedEfields),combinedrates))]
			E = -np.sort(-np.array(combinedEfields))


			[exponent]=fits.fitPower(E,cumsum/np.max(cumsum))
			self.combinedpowerfits = self.combinedpowerfits + [[windowperiod,exponent,probtoRPYratio]]
			#use PDF to determine mean and standard deviation of underlying normal distribution
			[guessMean,guessStd]=fits.getGuesses(E,counts,False)


			#fit to the datapoints (CDF has a probability of 1 at the first datapoint)
			[mean,std,loc]=fits.fitLognormalCDF(E,cumsum/np.max(cumsum),guessMean,guessStd,False)

			self.combinedlogfits = self.combinedlogfits + [[windowperiod,mean,std,loc,probtoRPYratio]]
			if(plot):
				if(j!=0):
					continue
				plt.figure()
				plt.loglog()

				plt.plot(E,fits.powerlaw(E,exponent)*probtoRPYratio, lw=1,label = "Powerfit, field averaged over "+str(windowperiod)+" seconds")
				plt.plot(E,fits.logcdf(np.array(E),mean,np.abs(std),loc)*probtoRPYratio, lw=1,label = "Logfit, field averaged over "+str(windowperiod)+" seconds")
				plt.plot(E,rates,'.', lw=1,label = "Field averaged over "+str(windowperiod)+" seconds")

				plt.legend()
				plt.title('Rate geoelectric field is above threshold')
				plt.xlabel('Geoelectric Field (V/km)')
				plt.ylabel('Average rate per year distinct contiguous sample average is above E field (counts/year)')
			
				plt.show()

	def calcGCcoords(self):
		self.GClatitudes=[]
		self.GClongitudes=[]
		longInterval=int(np.floor(Params.longituderes/0.25))
		latInterval=int(np.floor(Params.latituderes/0.25))
		with h5py.File(Params.globalcondloc, 'r') as f:
			latitudes = f.get('LATGRID')[0]
			longitudes = f.get('LONGRID')[:,0]
			lenlat=len(latitudes)
			lenlong=len(longitudes)
			for longkey in range(0,lenlong-1,longInterval):
				for latkey in range(0,lenlat-1,latInterval):
					latitude = f['LATGRID'][longkey,latkey]
					longitude=f['LONGRID'][longkey,latkey]
					if(longkey==0):
						self.GClatitudes.append(latitude)
					if(latkey==0):
						self.GClongitudes.append(longitude)

	def getClosestCoordsIds(self,latitude,longitude):
		mindelta=np.inf
		latindex=np.nan

		for i in range(0,len(self.GClatitudes)):
			delta = abs(self.GClatitudes[i]-latitude)
			if(delta<mindelta):
				mindelta=delta
				latindex=i
		mindelta=np.inf
		longindex=np.nan;
		for i in range(0,len(self.GClongitudes)):
			delta = abs(self.GClongitudes[i]-longitude)
			if(delta<mindelta):
				mindelta=delta
				longindex=i

		return [latindex,longindex]

	def compareAllTFsites(self):
		loaddirectory=Params.allTFsitesdir
		allfiles=glob.glob(loaddirectory+'*.edi')
		self.loadApparentCond()
		self.calcGCcoords()
		print(len(allfiles))
		self.allTFsiteAppcs=[]
		self.allGCmodelAppcs=[]
		self.allmaglats=[]
		self.allmaglongs=[]
		for i in range(0,len(allfiles)):
			fn=allfiles[i]
			if(fn=='.' or fn=='..'):
				continue
			print(fn)

			tfsite=TFsite()
			tfsite.initFromfile(fn)
			if(not tfsite.initFromfile(fn)):
				print('error with site '+str(fn)+', skipping')
				continue
			magcoords=self.geotomag(tfsite.lat,tfsite.long)
			tfsite.maglat=magcoords.data[0][1]
			tfsite.maglong=magcoords.data[0][2]
			print('tfsite.maglat')
			print(tfsite.maglat)
			print('tfsite.maglong')
			print(tfsite.maglong)

				
			[tfapparentc,gcapparentc]=self.compareGCmodeltoTFsite(tfsite)#8mHz
			self.allTFsiteAppcs.append(tfapparentc)
			self.allGCmodelAppcs.append(gcapparentc)
			self.allmaglats.append(tfsite.maglat)
			self.allmaglongs.append(tfsite.maglong)
		np.save('AllTFSitesCorrelation',[np.array(self.allTFsiteAppcs),np.array(self.allGCmodelAppcs),np.array(self.allmaglats),np.array(self.allmaglongs)])

	def loadGCtoTFcomparison(self):
		[self.allTFsiteAppcs,self.allGCmodelAppcs,self.allmaglats,self.allmaglongs]=np.load('AllTFSitesCorrelation.npy')

	def plotGCtoTFcomparison(self):
		plt.figure()
		plt.loglog()
		plt.scatter(self.allTFsiteAppcs,self.allGCmodelAppcs)
		plt.xlabel('EMTF (accurate) conductivities, (Ohm m)^-1')
		plt.ylabel('GC model (inaccurate) conductivities (Ohm m)^-1')
		plt.show()
		plt.figure()
		df=pd.DataFrame({'longs':(self.allmaglongs).astype(np.float),'lats':self.allmaglats.astype(np.float),'TF':self.allTFsiteAppcs.astype(np.float),'GC':self.allGCmodelAppcs.astype(np.float)})
		# dfTFfiltered = df[np.multiply((df['TF'].T > 10**-6), (df['TF'].T < 10**2))]
		dfTFfiltered = df[np.multiply((1/df['TF'].T > 10**0), (1/df['TF'].T < 10**4))]
		dfGCfiltered = df[np.multiply((df['GC'].T > 10**-3), (df['GC'].T < 10**1))]
		dfResLogTF=dfTFfiltered
		dfResLogGC=dfGCfiltered
		dfResLogTF['TF']=1/dfTFfiltered['TF']#np.log(1/dfTFfiltered['TF'])
		dfResLogGC['GC']=np.log(dfGCfiltered['GC'])
		crs={'init':'epsg:3857'}#4326'}
		# plt.figure()
		geometryTF=[Point(xy) for xy in zip(dfResLogTF['lats'],dfResLogTF['longs'])]
		geometryGC=[Point(xy) for xy in zip(dfResLogGC['lats'],dfResLogGC['longs'])]
		geo_df_TF=gpd.GeoDataFrame(dfResLogTF,crs=crs,geometry=geometryTF)
		geo_df_GC=gpd.GeoDataFrame(dfResLogGC,crs=crs,geometry=geometryGC)
		# print(dfResLog['TF'])
		plt.figure()
		plt.plot(np.array(dfResLogTF['longs'].as_matrix()))
		plt.show()
		print('from TF')
		fig, ax = plt.subplots(1, 1)
		minTF=np.min(np.array(dfResLogTF['TF'].as_matrix()))
		maxTF=np.max(np.array(dfResLogTF['TF'].as_matrix()))

		# ctx.add_basemap(ax)

		geo_df_TF.plot(column='TF',ax=ax,legend=True,\
			#cmap='hot',\
			norm=colors.LogNorm(vmin=minTF, vmax=maxTF))
		plt.show()
		
		print('from gcmodel')
		fig, ax = plt.subplots(1, 1)
		geo_df_GC.plot(column='GC',ax=ax,legend=True)
		plt.show()
		# plt.show()

		# GCmodelMap = np.zeros((len(self.allmaglongs), len(self.allmaglats)))
		# comb=pd.DataFrame({'long':self.allmaglongs,'lats':self.allmaglats})

		# plt.imshow(np.flipud(GCmodelMap), cmap='hot', interpolation='nearest')
		# plt.title('Overheat fractions '+str(r)+' per year, max of all durations')
		# cbar = plt.colorbar()
		plt.show()


	def compareGCmodeltoTFsite(self,tfsite):
			[latid,longid]=self.getClosestCoordsIds(tfsite.maglat,tfsite.maglong)
			print(latid)
			print(longid)
			freqindex=tfsite.getClosestFreqIndex(self.f)
			print('freqindex')
			print(freqindex)
			tfsite.calcApparentc()
			tfapparentc=tfsite.apparentc[freqindex]
			print('tfapparentc')
			print(tfapparentc)
			gcapparentc=self.apparentCondMap[latid,longid]
			print('gcapparentc')
			print(gcapparentc)
			return [tfapparentc,gcapparentc]