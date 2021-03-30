import Params
import fits
from scipy.interpolate import griddata
import matplotlib
import matplotlib.colors as colors
import glob
import Model.TFsite
from Model.TFsite import TFsite
import Model.MTsite
from Model.MTsite import MTsite
import geopandas

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
import geoplot as gplt
from shapely.geometry import Point
import contextily as ctx
from scipy.stats import pearsonr
import matplotlib.ticker as mticker

class EarthModel:

	# initialize a site by providing a list of frequencies w to determine the transfer function.
	# if a path to a TF site EDI document is provided, the init function 
	def __init__(self):
		Params.importIfNotAlready()
		self.refcumsum=[]
		self.refEfields=[]
		self.refmatchedEfields=[]
		self.combinedmatchedpowerfits=[]
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
		self.apparentCondMap=np.load(Params.apparentcondloc,allow_pickle=True)
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

		np.save(Params.durationRatiosDir,self.averagedRatios)

		if(plot):
			# plt.show()
			# plt.figure()
			plt.plot(self.allwindowperiods,self.averagedRatios,lw=7,label = "Mean Ratio")
			plt.legend()
			plt.title('Peak Geoelectric Pulse Ratios')
			plt.xlabel('Pulse Duration (s)')
			plt.ylabel('Ratio of Peak Field to 60 Second Peak Field')

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

			# plt.f/=igure()
			# # plt.loglog(fitx,mtsite.windowedCounts[2])
			# # plt.show()
			# # mtsite.loadWindowedCounts()
			# # print(np.array(mtsite.windowedCounts[1]))

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

	def formatticklabels(self,minval,maxval,pp):
		print('formatting')
		colourbar = pp.get_figure().get_axes()[1]
		ticks_loc = colourbar.get_xticks().tolist()
		ticks_loc.insert(0,minval)
		ticks_loc.append(maxval)
		colourbar.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
		
		labels=[]
		for i in range(0,len(ticks_loc)):
			x=ticks_loc[i]
			exponent = int(np.log10(x))
			coeff = x/10**exponent
			if(i==0 or i==len(ticks_loc)):
				labels.append(r"${:2.1f} \times 10^{{ {:2d} }}$".format(coeff,exponent))					
			else:
				labels.append(r"${:2.1f} \times 10^{{ {:2d} }}$".format(coeff,exponent))
		colourbar.set_xticklabels(labels,rotation=33)


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

	# find E -fields at all locations on conductivity map for each duration
	# Correct for geomagnetic latitude and apparent conductivity
	def calcGlobalEfields(self, ratePerYears):
		self.windowedEmaps=[]
		for r in ratePerYears:
			print('rate per year: '+str(r))
			windowedEmapsatrate=[]
			#save a map of e fields for each duration at this rate per Year
			# for i in range(0,len(self.allwindowperiods)):
			i=0
			windowperiod=self.allwindowperiods[i]
			mean=self.combinedlogfits[i][1]
			std=self.combinedlogfits[i][2]
			loc=self.combinedlogfits[i][3]
			ratio=self.combinedlogfits[i][4]
			exponent=self.combinedpowerfits[i][1]
			slope=self.combinedpowerfits[i][2]

			# exponent=self.combinedpowerfits[i*5+1]
			refFieldLog=fits.logcdfxfromy(r/ratio,mean,std,loc)
			print('refFieldLog')
			print(refFieldLog)
			refFieldPower=fits.powerlawxfromy(r,slope,exponent)
			refField=refFieldLog
			print('refField level rest of map is proportional to this field (V/km), logfit:'+str(refFieldLog))
			print('refField level rest of map is proportional to this field (V/km), powerfit:'+str(refFieldPower))
			EArr = np.zeros((len(self.GClongitudes), len(self.GClatitudes)))
			minE=10000
			maxE=-1
			maxc=np.max(self.apparentCondMap)
			minc=np.min(self.apparentCondMap)

			allcs=[]
			landcs=[]
			latlist=[]
			longlist=[]
			Elist=[]
			Cadjust=[]
			mlds=[]
			for latkey in range(0,len(self.GClatitudes)):
				latitude = self.GClatitudes[latkey]
				for longkey in range(0,len(self.GClongitudes)):
					longitude = self.GClongitudes[longkey]
					c=self.apparentCondMap[latkey,longkey]
					magcoords = self.geotomag(latitude, longitude)
					maglat=magcoords.data[0][1]
					maglong=magcoords.data[0][2]
					mld = self.getMagLatDivisor(maglat)

					if(maglat>-70 and maglat <80):
		
					# if(c>=5*10**-1 or c<10**-4): #ocean, or messed up data
					# 	E=np.nan
					# else:
						allcs.append(c)
						E =refField*np.sqrt(self.refconductivity)/np.sqrt(abs(c))*mld

						if(np.sqrt(self.refconductivity)/np.sqrt(abs(c))==0):
							print('c')
							print(c)
							quit()
						if(mld==0):
							print('mld')
							print(mld)
							quit()
						if(refField==0):
							print('reffield')
							quit()
						if(E==0):
							print('E')
							quit()
						landcs.append(c)

						if(maxE < E):
							maxE = E

						if(minE > E):
							minE = E

						EArr[latkey][longkey]=E
						latlist.append(latitude)
						longlist.append(longitude)
						Elist.append(E)
						Cadjust.append(1/np.sqrt(abs(c)))
						mlds.append(mld)
			

			print('minE')
			print(minE)
			print('maxE')
			print(maxE)
			# plt.figure()
			# plt.xscale("log")
			# [xland,yland]=fits.binlognormaldist(landcs,[],1)
			# [xall,yall]=fits.binlognormaldist(allcs,[],3)
			# plt.plot(xland,yland)
			# plt.plot(xall,yall)
			# plt.show()

			world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
			df=pd.DataFrame({'longs':np.array(longlist),'lats':np.array(latlist),'E':np.array(Elist),'Cadjust':np.array(Cadjust),'mld':np.array(mlds),'allcs':np.array(allcs)})
			crs={'init':'epsg:3857'}
			geometry=[Point(xy) for xy in zip(df['longs'],df['lats'])]
			geo_df=gpd.GeoDataFrame(df,crs=crs,geometry=geometry)
			# ax=gplt.kdeplot(geo_df,column='E',cmap='Reds',shade=True)
			print('wierd')
			###########E Field Plot ####################
			print('abouttoplto')
			plt.close()
			ax=geo_df.plot(column='E',legend=True,
			cmap='viridis',\
			legend_kwds={'label': 'Field Level (V/km)','orientation': "horizontal"},\
			norm=colors.LogNorm(vmin=minE, vmax=maxE))
			pp=gplt.polyplot(world,ax=ax,zorder=1)
			self.formatticklabels(minE,maxE,pp)
			# plt.rcParams['text.usetex'] = True
			plt.title('Magnitude of Peak '+str(windowperiod) +' Second Geoelectric Field\n1 in '+str(1/r)+' Year Storm')			
			print('not s7')
			plt.savefig(Params.globalEfieldPlots+'Results'+str(r)+'perYearWindow'+str(windowperiod)+'s.png')
			plt.show()


			# Create colorbar as a legend
			# sm = plt.cm.ScalarMappable(cmap=’BuGn’, norm=plt.Normalize(vmin=minE, vmax=maxE))# empty array for the data range
			# sm._A = []# add the colorbar to the figure
			# cbar = fig.colorbar(sm)
			# plt.figure()
			# plt.imshow(np.flipud(np.log(EArr[12:,:])), cmap='hot', interpolation='nearest')
			# plt.title('E field levels, '+str(r)+' per year, '+str(
			# cbar = plt.colorbar()
			# current_cmap = matplotlib.cm.get_cmap()
			# current_cmap.set_bad(color='blue')
			# cbar.set_label('E field (V/km)')
			for j in range(0,len(self.allwindowperiods)):
				windowperiod = self.allwindowperiods[j]
			
				windowedEmapsatrate = windowedEmapsatrate + [windowperiod,self.averagedRatios[j]*EArr]
			
			self.windowedEmaps = self.windowedEmaps + [r,windowedEmapsatrate]
		np.save(Params.globalEfieldData,self.windowedEmaps)

		###########Conductivity Plot ####################

		ax=geo_df.plot(column='allcs',legend=True,legend_kwds={'label': 'Apparent Conductivity (Ohm m)^-1','orientation': "horizontal"}, cmap='viridis',norm=colors.LogNorm(vmin=minc, vmax=maxc))
		pp=gplt.polyplot(world,ax=ax,zorder=1)
		self.formatticklabels(minc,maxc,pp)
		plt.title('Global Conductivity at 1/120 Hz')


		plt.savefig(Params.figuresDir+'Conductivity.png')

		# gplt.kdeplot(geo_df,column='E',)
		plt.show()

		
		###########Conductivity Adjust Plot ####################

		mincadj=1/np.sqrt(maxc)
		maxcadj=1/np.sqrt(minc)
		ax=geo_df.plot(column='Cadjust',legend=True,legend_kwds={'label': 'Coefficient on Reference Field Level','orientation': "horizontal"}, cmap='viridis',norm=colors.LogNorm(vmin=mincadj, vmax=maxcadj))
		
		pp=gplt.polyplot(world,ax=ax,zorder=1)
		self.formatticklabels(mincadj,maxcadj,pp)

		plt.title('Geoelectric Field Multiplier\nGlobal Conductivity')
		plt.savefig(Params.figuresDir+'ConductivityMultiplier.png')

		# gplt.kdeplot(geo_df,column='E',)
		plt.show()

		

		###########Maglat Adjust Plot ####################
		mldmin=np.min(df['mld'].values)
		mldmax=np.max(df['mld'].values)
		ax=geo_df.plot(column='mld',legend=True,legend_kwds={'label': 'Coefficient on Reference Field Level','orientation': "horizontal"}, cmap='viridis',norm=colors.LogNorm(vmin=mldmin, vmax=mldmax))

		pp=gplt.polyplot(world,ax=ax,zorder=1)

		colourbar = pp.get_figure().get_axes()[1]
		ticks_loc = colourbar.get_xticks().tolist()
		ticks_loc.insert(0,np.min(df['mld'].values))
		ticks_loc.append(np.max(df['mld'].values))
		colourbar.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
		self.formatticklabels(mincadj,maxcadj,pp)

		plt.title('Geoelectric Field Multiplier\n Magnetic Latitude')
		plt.savefig(Params.figuresDir+'MagneticLatitudeMultiplier.png')


		print('s10')
		plt.show()
	#adjusts ratesperyear for all the sites to the reference maglat and conductivity and sets the result to properties of this EarthModel instance
	def calcReferenceRateperyear(self,tfsites,mtsites,fittype,plot):
		self.refcumsum=[]
		self.refEfields=[]
		allcumsum=[]
		allE=[]
		allmatchedE=[]#populated later on, will remain empty for now
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
			allmatchedEatwindow=[]
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
				# finalAdjustment=self.adjustEfieldsToMatch(adjustedtomaglat)

				# quit()
				allEatwindow.append(np.array(adjustedtomaglat))
				allmatchedEatwindow.append(1)
				allcumsumatwindow.append(np.array(cumsum))
				allcountsatwindow.append(np.array(counts))
				allyearsatwindow.append(years)

			allcumsum.append(allcumsumatwindow)
			allE.append(allEatwindow)
			allmatchedE.append(allmatchedEatwindow)
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

			print('s11')
			plt.show()

		self.refcumsum=allcumsum
		self.refEfields=allE
		self.refmatchedEfields=allmatchedE
		self.refCounts=allCounts
		self.refYears=allYears

	#adjust mtsite fits to reference conductivity.
	#this assumes they all have
	def adjustEfieldsToRefCond(self,apparentcond,Efields,windowperiod):
		#find the closest TFsite frequency match to the given window period.
		#we want the apparent conductivity at the average frequency of the EMP, which can be approximated by 1/windowperiod

		return np.sqrt(apparentcond)/np.sqrt(self.refconductivity)*np.array(Efields)

	def adjustEfieldsToMatch(self,plot):

		# for j in range(0,len(self.allwindowperiods)):
		j=0
		if(plot and (j==0)):
			plt.figure()
			plt.loglog()
		exponent=self.combinedpowerfits[j][1]
		scale=self.combinedpowerfits[j][2]

		# combinedEfields=np.array([])
		# combinedrates=np.array([])
		# combinedcumsum=np.array([])
		# combinedcounts=np.array([])
		# combinedyears=0


		for i in range(0,len(self.refYears)):
			if(not Params.useMTsite[i]):
				continue
			years=self.refYears[i][j]
			Efields=self.refEfields[i][j]
			rates=np.array(self.refcumsum[i][j])/self.refYears[i][j]
			cumsum=self.refcumsum[i][j]
			counts=self.refCounts[i][j]

			# linfitEfields=np.power((rates[rates>1]/scale),(1/exponent))
			linfitEfields=np.power((rates/scale),(1/exponent))

			#mean ratio between the average powerfit value and the E fields for this data.
			# adjustment=np.mean(linfitEfields/Efields[rates>1])
			adjustment=np.mean(linfitEfields/Efields)
			print('adjustment')
			print(adjustment)

			self.refmatchedEfields[i][j]=adjustment
			if(plot and (j==0)):
				plt.plot(Efields*adjustment,rates)
				plt.plot(Efields,np.array(Efields)**exponent*scale)

		if(plot and (j==0)):
			print('s12')
			plt.show()
	def loadCombinedFits(self):
		self.combinedlogfits=np.load(Params.combinedFitParamsLoc,allow_pickle=True)


	def calcMatchedCombinedRates(self,plot):
		self.combinedexponents=[]
		#combine the fields into one distribution
		# for j in range(0,len(self.allwindowperiods)):
		j=0
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
			combinedEfields=np.append(combinedEfields,self.refmatchedEfields[i][j]*self.refEfields[i][j])
			ratestmp=np.array(self.refcumsum[i][j])/self.refYears[i][j]
			combinedrates=np.append(combinedrates,ratestmp)
			combinedcumsum=np.append(combinedcumsum,self.refcumsum[i][j])
			combinedcounts=np.append(combinedcounts,self.refCounts[i][j])

		[EfinalAll,cumsumfinalAll,countsfinalAll]=fits.combinecounts(combinedEfields,combinedcounts)
		mask=EfinalAll>EfinalAll[0]/10
		Efinal=EfinalAll[mask]
		cumsumfinal=cumsumfinalAll[mask]
		countsfinal=countsfinalAll[mask]
		# Efinal=EfinalAll
		# cumsumfinal=cumsumfinalAll
		# countsfinal=countsfinalAll


		print('combinedyears')
		print(combinedyears)
		rates=cumsumfinal/combinedyears			
		probtoRPYratio=np.max(cumsumfinal)/combinedyears

		[exponent]=fits.fitPower(Efinal,cumsumfinal/np.max(cumsumfinal))
		# self.combinedmatchedpowerfits = self.combinedmatchedpowerfits + [[windowperiod,linearfit[0],np.exp(linearfit[1])]]
		#use PDF to determine mean and standard deviation of underlying normal distribution
		[guessMean,guessStd]=fits.getGuesses(Efinal,countsfinal,False)


		#fit to the datapoints (CDF has a probability of 1 at the first datapoint)
		[mean,std,loc]=fits.fitLognormalCDF(Efinal,cumsumfinal/np.max(cumsumfinal),guessMean,guessStd,False)


		self.combinedlogfits = self.combinedlogfits + [[windowperiod,mean,std,loc,probtoRPYratio]]
		np.save(Params.combinedFitParamsLoc,self.combinedlogfits,allow_pickle=True)
		if(plot):
			plt.figure()
			# plt.loglog()
			# plt.plot(rates,fits.logcdfxfromy(rates/,mean,std,loc))
			# plt.plot(cumsumfinal/max(cumsumfinal),fits.logcdfxfromy(cumsumfinal/max(cumsumfinal),mean,std,loc))
			print('s13')
			# plt.show()

			plt.figure()
			plt.loglog()
			plt.plot(Efinal,fits.powerlaw(Efinal,exponent)*probtoRPYratio, lw=1,label = "Powerfit, field averaged over "+str(windowperiod)+" seconds")
			# plt.plot(Efinal,E**linearfit[0]*np.exp(linearfit[1]), lw=1,label = "Powerfit, field averaged over "+str(windowperiod)+" seconds")
			plt.plot(Efinal,fits.logcdf(np.array(Efinal),mean,np.abs(std),loc)*probtoRPYratio, lw=1,label = "Logfit, field averaged over "+str(windowperiod)+" seconds")

			plt.plot(Efinal,rates,'.', lw=1,label = "Field averaged over "+str(windowperiod)+" seconds")
			plt.plot(fits.logcdfxfromy(rates/probtoRPYratio,mean,std,loc),rates,'.', lw=1,label = "xfromy")

			plt.legend()
			plt.title('Rate geoelectric field is above threshold')
			plt.xlabel('Geoelectric Field (V/km)')
			plt.ylabel('Average rate per year distinct contiguous sample average is above E field (counts/year)')
		
			print('s14')
			plt.show()



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
			
			#sort by descending E field
			counts = [x for _,x in sorted(zip(-np.array(combinedEfields),combinedcounts))]
			cumsum = [x for _,x in sorted(zip(-np.array(combinedEfields),combinedcumsum))]
			rates = [x for _,x in sorted(zip(-np.array(combinedEfields),combinedrates))]
			E = -np.sort(-np.array(combinedEfields))

			linearfit=np.polyfit(np.log(E),np.log(rates),1)
			# linfun=np.poly1d(linearfit)
			self.combinedpowerfits = self.combinedpowerfits + [[windowperiod,linearfit[0],np.exp(linearfit[1])]]
			if(plot):
				if(j!=0):
					continue
				plt.loglog()

				plt.plot(E,E**linearfit[0]*np.exp(linearfit[1]), lw=1,label = "Powerfit, field averaged over "+str(windowperiod)+" seconds")
				# plt.plot(E,fits.logcdf(np.array(E),mean,np.abs(std),loc)*probtoRPYratio, lw=1,label = "Logfit, field averaged over "+str(windowperiod)+" seconds")
				plt.plot(E,rates,'.', lw=1,label = "Field averaged over "+str(windowperiod)+" seconds")

				plt.legend()
				plt.title('Rate geoelectric field is above threshold')
				plt.xlabel('Geoelectric Field (V/km)')
				plt.ylabel('Average rate per year distinct contiguous sample average is above E field (counts/year)')
				
				print('s15')
				plt.show()

	def calcGCcoords(self):
		print('calculating GC coords')
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
				longitude=f['LONGRID'][longkey,0]
				self.GClongitudes.append(longitude)

			for latkey in range(0,lenlat-1,latInterval):
				latitude = f['LATGRID'][0,latkey]
				self.GClatitudes.append(latitude)

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
		self.allTFsiteAppcs=[]
		self.allGCmodelAppcs=[]
		self.alllats=[]
		self.alllongs=[]
		numsucceed=0

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
			print('tfsite')
			print(tfsite.lat)
			print('tfsite')
			print(tfsite.long)

			# magcoords=self.geotomag(tfsite.lat,tfsite.long)
			# tfsite.maglat=magcoords.data[0][1]
			# tfsite.maglong=magcoords.data[0][2]
			# print('tfsite.maglat')
			# print(tfsite.maglat)
			# print('tfsite.maglong')
			# print(tfsite.maglong)

				
			[tfapparentc,gcapparentc]=self.compareGCmodeltoTFsite(tfsite)#8mHz
			self.allTFsiteAppcs.append(tfapparentc)
			self.allGCmodelAppcs.append(gcapparentc)
			self.alllats.append(tfsite.lat)
			self.alllongs.append(tfsite.long)
			numsucceed=numsucceed+1
			# print('self.alllats')
			# print('self.alllongs')
			# print(self.alllats)
			# print(self.alllongs)
		np.save('AllTFSitesCorrelation',[np.array(self.allTFsiteAppcs),np.array(self.allGCmodelAppcs),np.array(self.alllats),np.array(self.alllongs)])
		print('len(allfiles)')
		print(len(allfiles))
		print('numfail')
		print(len(allfiles)-numsucceed)

	def loadGCtoTFcomparison(self):
		[self.allTFsiteAppcs,self.allGCmodelAppcs,self.alllats,self.alllongs]=np.load('Data/SmallData/AllTFSitesCorrelation.npy',allow_pickle=True)

	def findAverages(self,longitudes,latitudes,values,avgwin):
		allaveragedpoints=[]
		for i in range(0,len(values)):
			longitude=longitudes[i]
			latitude=latitudes[i]
			# print('longitude')
			# print(longitude)
			# print('latitude')
			# print(latitude)
			pointstoavg=np.array([])
			for j in range(0,len(values)):
				if((((longitudes[j]<longitude+avgwin) and (longitudes[j]>longitude-avgwin))) and ((latitudes[j]<latitude+avgwin) and (latitudes[j]>latitude-avgwin))):
					# print('values[j]')
					# print(values[j])
					pointstoavg=np.append(pointstoavg,values[j])
		
			allaveragedpoints.append(np.mean(pointstoavg))
		return np.array(allaveragedpoints)

	def findPredictivePower(self,TFcond,GCcond):
		ratios=np.array([])
		for i in range(0,len(TFcond)):
			tf=TFcond[i]
			gc=GCcond[i]
			ratios=np.append(ratios,tf/gc)
		plt.figure()
		plt.loglog()
		plt.scatter(TFcond,ratios)
		print('s16')
		plt.show()
		print('np.mean(ratio)')
		print(np.mean(ratios))
		print('np.mean(ratio)')
		print(np.median(ratios))

		print('sqrt(np.mean(ratio-**2)')
		print('np.sqrt(np.sum((ratios-np.mean(ratios))**2)/len(ratios))')
		print(np.sqrt(np.sum((ratios-np.mean(ratios))**2)/len(ratios)))
		print('np.sqrt(np.sum((ratios-np.mean(ratios))**2)/len(ratios))')
		print(np.sqrt(np.sum(ratios**2)/len(ratios)))

		rmsratio=np.sqrt(np.sum(ratios**2)/len(ratios))

		print('np.std(ratio)')
		print(np.std(ratios))
		print('2 sigma')
		print(2*np.std(ratios))
		print('rms ratio')
		print(np.std(ratios))
		print('2 rms ratio')
		print(2*rmsratio)

	def plotGCtoTFcomparison(self):
		plt.figure()
		# plt.loglog()
		# plt.scatter(self.allTFsiteAppcs,self.allGCmodelAppcs)
		# plt.xlabel('EMTF (accurate) conductivities, (Ohm m)^-1')
		# plt.ylabel('GC model (inaccurate) conductivities (Ohm m)^-1')
		print('s17')
		# plt.show()
		alllongs=np.array(self.alllongs).astype(np.float)
		alllats=np.array(self.alllats).astype(np.float)
		allTFsiteAppcs=np.real(np.array(self.allTFsiteAppcs)).astype(np.float)
		allGCmodelAppcs=np.real(np.array(self.allGCmodelAppcs)).astype(np.float)

		# target grid to interpolate to
		xi = yi = np.arange(0,1.01,0.01)
		xi,yi = np.meshgrid(xi,yi)

		# set mask
		# mask = (xi > 0.5) & (xi < 0.6) & (yi > 0.5) & (yi < 0.6)

		# interpolate
		# zi = griddata((x,y),z,(xi,yi),method='linear')


		# mask out the field
		# zi[mask] = np.nan

		# plot
		plt.figure()
		ax = fig.add_subplot(111)
		# plt.contourf(xi,yi,zi,np.arange(0,1.01,0.01))
		# plt.plot(alllongs,alllats,'k.')
		# plt.xlabel('xi',fontsize=16)
		# plt.ylabel('yi',fontsize=16)
		print('s18')
		# plt.show()
		# plt.savefig('interpolated.png',dpi=100)
		# plt.close(fig)
		# quit()
		df=pd.DataFrame({'longs':np.array(self.alllongs).astype(np.float),'lats':np.array(self.alllats).astype(np.float),'TF':np.real(np.array(self.allTFsiteAppcs).astype(np.float)),'GC':np.real(np.array(self.allGCmodelAppcs).astype(np.float))})

		# print(np.sum(np.isnan(np.array(self.allTFsiteAppcs).astype(np.float))))
		# print(np.sum(np.isnan(np.array(self.allGCmodelAppcs).astype(np.float))))
		# corr, _ = pearsonr(np.real(self.allTFsiteAppcs.astype(np.float)),np.real( self
		# dfTFfiltered = df[np.multiply((df['TF'].T > 10**-6), (df['TF'].T < 10**2))]
		# dfTFfiltered = df[np.multiply(np.multiply(np.multiply((1/df['TF'].T > 10**0), (1/df['TF'].T < 10**4)),df['GC'].T >3*10**-3),df['GC'].T <2*10**-1)]
		# dfTFfiltered = df[np.multiply(np.multiply((1/df['TF'].T > 10**0), (1/df['TF'].T < 10**4)),df['GC'].T >3*10**-3)]
		dfTFfiltered = df[np.multiply((1/df['TF'].T > 10**0), (1/df['TF'].T < 10**4))]

		plt.figure()
		# plt.loglog()
		# plt.scatter(dfTFfiltered['TF'],dfTFfiltered['GC'])
		# plt.xlabel('EMTF (accurate) conductivities, (Ohm m)^-1')
		# plt.ylabel('GC model (inaccurate) conductivities (Ohm m)^-1')
		print('s19')
		# plt.show()

		# dfGCfiltered = df[np.multiply((df['GC'].T > 10**-3), (df['GC'].T < 10**1))]
		dfResLogTF=dfTFfiltered
		dfResLogGC=dfTFfiltered
		dfTFfiltered['logTF']=np.log(dfTFfiltered['TF'])#1/dfTFfiltered['TF']#np.log(1/dfTFfiltered['TF'])
		dfTFfiltered['logGC']=np.log(dfTFfiltered['GC'])#1/dfTFfiltered['TF']#np.log(1/dfTFfiltered['TF'])
		# maxTF=np.max(dfTFfiltered['logTF'])#1/dfTFfiltered['TF']#np.log(1/dfTFfiltered['TF'])
		# minTF=np.min(dfTFfiltered['logTF'])#1/dfTFfiltered['TF']#np.log(1/dfTFfiltered['TF'])
		# maxGC=np.max(dfTFfiltered['logGC'])#1/dfTFfiltered['GC']#np.log(1/dfTFfiltered['GC'])
		# minGC=np.min(dfTFfiltered['logGC'])#1/dfTFfiltered['TF']#np.log(1/dfTFfiltered['TF'])
		# dfTFfiltered['logGCadjusted']=(dfTFfiltered['logGC']-minGC)*(maxTF-minTF)/(maxGC-minGC)+minTF#1/dfTFfiltered['TF']#np.log(1/dfTFfiltered['TF'])
		corr=dfTFfiltered['logTF'].corr(dfTFfiltered['logGC'])
		print('correlation: %.3f' % corr)
		# dfTFfiltered['GCadjusted']=np.exp(dfTFfiltered['logGCadjusted'])

		# dfResLogGC['GC']=dfGCfiltered['GC']#np.log(dfGCfiltered['GC'])
		crs={'init':'epsg:3857'}#4326'}
		plt.figure()
		geometryTF=[Point(xy) for xy in zip(dfTFfiltered['longs'],dfTFfiltered['lats'])]
		geometryGC=[Point(xy) for xy in zip(dfTFfiltered['longs'],dfTFfiltered['lats'])]
		geo_df_TF=gpd.GeoDataFrame(dfTFfiltered,crs=crs,geometry=geometryTF)
		geo_df_GC=gpd.GeoDataFrame(dfTFfiltered,crs=crs,geometry=geometryGC)
		# print(dfResLog['TF'])
		plt.figure()
		tf=np.array(dfTFfiltered['TF'].values)
		longs=np.array(dfTFfiltered['longs'].values)
		lats=np.array(dfTFfiltered['lats'].values)
		averages=np.exp(self.findAverages(longs,lats,np.log(tf),5))
		# averages=np.load('averages.npy',allow_pickle=True)
		# print(averages)
		# quit()
		# np.save('averages',averages)
		dfTFfiltered['TFavgd']=averages
		dfTFfiltered['TFavgdlog']=np.log(averages)

		plt.figure()
		# plt.loglog()
		plt.plot(dfTFfiltered['logGC'])
		plt.xlabel('EMTF (accurate) conductivities, (Ohm m)^-1')
		plt.ylabel('GC model (inaccurate) conductivities (Ohm m)^-1')
		print('s20')
		plt.show()
		
		plt.figure()
		# plt.loglog()
		plt.plot(dfTFfiltered['TFavgdlog'])
		plt.xlabel('EMTF (accurate) conductivities, (Ohm m)^-1')
		plt.ylabel('GC model (inaccurate) conductivities (Ohm m)^-1')
		print('s21')
		plt.show()
		TFavgd=np.array(dfTFfiltered['TFavgd'].values)
		GCarr=np.array(dfTFfiltered['GC'].values)
		# GCadjusted=GCarr*np.exp(np.log(np.mean(TFavgd))/np.log(np.mean(GCarr)))
		[xtf,ytf]=fits.binlognormaldist(TFavgd,[],4)

		[xgc,ygc]=fits.binlognormaldist(GCarr,[],2)
		
		plt.figure()
		plt.xscale("log")

		plt.plot(xtf,ytf,lw=1,label = "tf")
		plt.plot(xgc,ygc,lw=1,label = "gc")
		plt.legend()
		print('s22')
		plt.show()


		corr=dfTFfiltered['TFavgdlog'].corr(dfTFfiltered['logGC'])
		print('')
		corr, _ = pearsonr(np.log(TFavgd),np.log(GCarr))
		print('spatial average correlation: %.3f' % corr)
		print('')
		print('predictive power GC vs TF')
		self.findPredictivePower(TFavgd,GCarr)#*np.exp(np.log(np.mean(TFavgd))/np.log(np.mean(GCarr))))
		# self.findPredictivePower(TFavgd,np.ones(len(TFavgd))*np.mean(TFavgd))
				# gc=np.array(dfTFfiltered['GC'].values)
		print('s23')
		# plt.show()
		# calculate Pearson's correlation`
		print('from TF')
		print('from GC')

		fig, ax = plt.subplots(1, 1)
		minTF=np.min(np.array(dfTFfiltered['TF'].values))
		maxTF=np.max(np.array(dfTFfiltered['TF'].values))

		minGC=np.min(np.array(dfTFfiltered['GC'].values))
		maxGC=np.max(np.array(dfTFfiltered['GC'].values))
		minTFavgd=np.min(averages)
		maxTFavgd=np.max(averages)

		plt.figure()
		plt.loglog()
		plt.scatter(dfTFfiltered['TFavgd'],dfTFfiltered['GC'])
		plt.xlabel('EMTF (accurate) conductivities, (Ohm m)^-1')
		plt.ylabel('GC model (inaccurate) conductivities (Ohm m)^-1')
		print('s24')
		plt.show()

		# ctx.add_basemap(ax)

		# geo_df_TF.plot(column='TF',ax=ax,legend=True,\
		# 	cmap='viridis',\
		# norm=colors.LogNorm(vmin=minTF, vmax=maxTF))
		print('s25')
		# plt.show()

		fig, ax = plt.subplots(1, 1)
		geo_df_GC.plot(column='TFavgd',ax=ax,legend=True,
			cmap='viridis',\
		norm=colors.LogNorm(vmin=minTF, vmax=maxTF))
		print('s26')
		plt.show()
		
		fig, ax = plt.subplots(1, 1)
		geo_df_GC.plot(column='GC',ax=ax,legend=True,
			cmap='viridis',\
		norm=colors.LogNorm(vmin=minTF, vmax=maxTF))
		print('s27')
		plt.show()

		# print('from gcmodel')
		# fig, ax = plt.subplots(1, 1)
		# geo_df_GC.plot(column='GC',ax=ax,legend=True)
		print('s28')
		# plt.show()
		print('s29')
		# plt.show()

		# GCmodelMap = np.zeros((len(self.alllongs), len(self.alllats)))
		# comb=pd.DataFrame({'long':self.alllongs,'lats':self.alllats})

		# plt.imshow(np.flipud(GCmodelMap), cmap='hot', interpolation='nearest')
		# plt.title('Overheat fractions '+str(r)+' per year, max of all durations')
		# cbar = plt.colorbar()
		print('s30')
		# plt.show()


	def compareGCmodeltoTFsite(self,tfsite):
			[latid,longid]=self.getClosestCoordsIds(tfsite.lat,tfsite.long)
			print(latid)
			print(longid)
			print('frequency')
			print(self.f)
			freqindex=tfsite.getClosestFreqIndex(self.f)
			print('freqindex')
			print(freqindex)
			print('frequency from site')
			print(tfsite.f[freqindex])
			print('ratio')
			print(tfsite.f[freqindex]/self.f)

			print('len frequencies')
			print(len(tfsite.f))
			if(freqindex==0):
				print('lost the frequency of interest')
				quit()
			tfsite.calcApparentc()
			tfapparentc=tfsite.apparentc[freqindex]
			print('tfapparentc')
			print(tfapparentc)
			gcapparentc=self.apparentCondMap[latid,longid]
			print('gcapparentc')
			print(gcapparentc)
			return [tfapparentc,gcapparentc]

	def loadDurationRatios(self):
		self.averagedRatios=np.load(Params.durationRatiosDir,allow_pickle=True)
