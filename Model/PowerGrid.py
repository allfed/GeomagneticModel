import numpy as np
#calculate thermal rise from duration and GIC level
import Params
import h5py
import matplotlib.pyplot as plt
import geopandas
import pandas as pd
import matplotlib.colors as colors
import geopandas as gpd
import geoplot as gplt
from shapely.geometry import Point
import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import LogFormatter 
import matplotlib.ticker as mticker
import rasterio
import rasterio.features
import rasterio.warp
class PowerGrid:

	def __init__(self):
		Params.importIfNotAlready()
		self.windowedEmaps=[]
		self.HVTnames=Params.HVTnames
		self.GICperE=Params.GICperE

		self.temprise0=Params.temprise0
		self.temprise10=Params.temprise10
		self.temprise20=Params.temprise20
		self.temprise40=Params.temprise40
		self.temprise50=Params.temprise50
		self.temprise100=Params.temprise100
		self.temprise200=Params.temprise200
		self.pop25to40=Params.pop25to40
		self.pop0to25=Params.pop0to25
		self.pop40plus=Params.pop40plus
		self.tau=Params.tau
		self.latitudes=[]
		self.longitudes=[]

	def setWindowedEmaps(self,earthmodel):
		self.windowedEmaps=earthmodel.windowedEmaps
		self.latitudes=earthmodel.GClatitudes
		self.longitudes=earthmodel.GClongitudes
		self.geotomag=earthmodel.geotomag
		self.averagedRatios=earthmodel.averagedRatios
		self.combinedlogfits=earthmodel.combinedlogfits

	def loadEfieldMaps(self):
		loaddirectory=Params.globalEfieldData
		self.windowedEmaps=np.load(loaddirectory,allow_pickle=True)

	def calcOverheatMap(self):
		#whether to plot each duration
		plotintermediates=True
		print('calculating overheating')
		#for each rate per year calculated in earth model
		alloverheatfractions=np.array([])
		for i in range(0,int(np.floor(len(self.windowedEmaps))/2)):
			alldurations=[]
			allmeanoverheats=[]
			meanoverheat=np.array([])
			durations=np.array([])
			maxoverheatfractions=[]
			rateindex=i*2
			dataindex=i*2+1
			r=self.windowedEmaps[rateindex]
			print('rateperyear')
			print(r)
			EmapAtDuration=self.windowedEmaps[dataindex]
			latitudes=self.latitudes
			longitudes=self.longitudes

			lenlat=len(latitudes)
			lenlong=len(longitudes)

			# print('self.calcFractionOverTemp(.003,300)')
			# print(self.calcFractionOverTemp(0,10000))
			# # quit()
			# print('self.calcFractionOverTemp(3,100)')
			# print(self.calcFractionOverTemp(3,100))

			for j in range(0,int(np.floor(len(EmapAtDuration))/2)):
				durationindex = j*2
				Emapindex = j*2+1

				duration = EmapAtDuration[durationindex]
				durations=np.append(durations,duration)
				print('duration')
				print(duration)
				Emap = EmapAtDuration[1]#Emapindex]
				thisdurationoverheatmap = np.zeros((lenlong, lenlat))
				overheatFractionMap = np.zeros((lenlong,lenlat))
				latitudes=[]
				longitudes=[]
				allfractionsover=[]
				i=0
				for latkey in range(0,len(self.latitudes)):
					latitude = self.latitudes[latkey]
					for longkey in range(0,len(self.longitudes)):
						longitude = self.longitudes[longkey]
						E=Emap[latkey,longkey]
						fractionover=self.calcFractionOverTemp(self.averagedRatios[j]*E,duration)
						
						magcoords = self.geotomag(latitude, longitude)
						maglat=magcoords.data[0][1]
				
						if(maglat>-70 and maglat <80):
							allfractionsover.append(fractionover)
							latitudes.append(latitude)
							longitudes.append(longitude)
							thisdurationoverheatmap[latkey][longkey]=fractionover
							overheatFractionMap[latkey][longkey]=max(overheatFractionMap[latkey][longkey],fractionover)
							if(j==0):
								maxoverheatfractions.append(fractionover)
							else:
								maxoverheatfractions[i]=max(maxoverheatfractions[i],fractionover)
							i=i+1
				meanoverheat=np.append(meanoverheat,np.mean(allfractionsover))
				print('mean overheat fraction')
				print(np.mean(allfractionsover))
				print('max overheat fraction')
				print(np.max(allfractionsover))
				
				# [xall,yall]=fits.binnormaldist(allfractionsover,[],-1)

				# plt.figure()
				# plt.plot(xall,yall)
				# plt.title('distributions of overheat fractions,linear scale')
				# plt.show()

				if(plotintermediates):
					df=pd.DataFrame({'longs':np.array(longitudes),'lats':np.array(latitudes),'over':np.array(allfractionsover)})
					crs={'init':'epsg:3857'}
					geometry=[Point(xy) for xy in zip(df['longs'],df['lats'])]
					geo_df=gpd.GeoDataFrame(df,crs=crs,geometry=geometry)
					# ax=gplt.kdeplot(geo_df,column='over',cmap='Reds',shade=True)
					ax=geo_df.plot(column='over',legend=True,legend_kwds={'label': 'Fraction HVT Over','orientation': "horizontal"}, cmap='viridis',vmin=0, vmax=.5)
					# ax.
					# cbr = fig.colorbar(sm, cax=cax,)
					# cbr.ax.tick_params() 
					# cbr.ax.set_yticklabels(['{:.0f}'.format(x) for x in np.arange(.01, .4+.001, .001)], fontsize=16, weight='bold')
					world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
					pp=gplt.polyplot(world,ax=ax,zorder=1)
					
					plt.title('HVT Tie Bar Exceeding Short Term Temp Limit \n'+str(duration) +' Second Geoelectric Field\n 1 in '+str(1/r)+' Year Storm')
					# gplt.kdeplot(geo_df,column='E',)
					print(Params.overheatMapsDir+'ByDuration/'+str(1/r)+'/Overheat'+str(r)+'perYearWindow'+str(duration)+'s.png')
					plt.savefig(Params.overheatMapsDir+'ByDuration/'+str(1/r)+'/Overheat'+str(r)+'perYearWindow'+str(duration)+'s.png')

					# plt.show()
					# plt.close()
			alloverheatfractions=np.append(alloverheatfractions,[r,np.array(latitudes),np.array(longitudes),np.array(maxoverheatfractions)])
			print('np.len(alloverheatfractions)')
			print(len(alloverheatfractions))
			alldurations.append(durations)
			allmeanoverheats.append(meanoverheat)
			# plt.figure()
			# plt.plot(durations,meanoverheat)
			# plt.title('HVT Exceeding Short Term Temp Limit')
			# plt.ylabel('Fraction of HVT Population')
			# plt.xlabel('Duration (s)')
			# plt.show()
			print('maximum overheat fraction')
			print(np.max(maxoverheatfractions))
			df=pd.DataFrame({'longs':np.array(longitudes),'lats':np.array(latitudes),'over':np.array(maxoverheatfractions)})
			crs={'init':'epsg:3857'}
			geometry=[Point(xy) for xy in zip(df['longs'],df['lats'])]
			geo_df=gpd.GeoDataFrame(df,crs=crs,geometry=geometry)
			# ax=gplt.kdeplot(geo_df,column='over',cmap='Reds',shade=True)
			ax=geo_df.plot(column='over',legend=True,legend_kwds={'label': 'Fraction Exceeding','orientation': "horizontal"}, cmap='viridis',vmin=0, vmax=.5)
			# norm=colors.LogNorm(vmin=0.01, vmax=1))

			world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
			gplt.polyplot(world,ax=ax,zorder=1)
			plt.title('HVT Exceeding Short Term Temp Limit\n1 in '+str(1/r)+' Year Storm')
			plt.savefig(Params.overheatMapsDir+str(r)+'perYearMaxfromallWindows.png')
			# plt.show()
			# plt.close()
			# gplt.kdeplot(geo_df,column='E',)

		np.save(Params.overheatDataDir,alloverheatfractions)
		# plt.figure(1)
		# plt.imshow(np.flipud(overheatFractionMap[8:,:]), cmap='hot', interpolation='nearest',vmin=0,vmax=1)
		# plt.title('Overheat fractions '+str(r)+' per year, max of all durations')
		# cbar = plt.colorbar()
		# cbar.set_label('Fraction overheat (1=100%)')
		# plt.show()

	# #plots the relative propensity for HVT to overheat at the reference site by duration for a given rate per year storm
	# def plotOverheatByDuration(self):
	# 	plt.figure()
	# 	logfits=self.combinedlogfits
	# 	print(logfits)
	# 	durationratios=self.averagedRatios
	# 	for i in range(0,int(np.floor(len(self.windowedEmaps))/2)):
	# 		alldurations=[]
	# 		allfractionsover=[]
	# 		durations=np.array([])
	# 		rateindex=i*2
	# 		dataindex=i*2+1
	# 		r=self.windowedEmaps[rateindex]
	# 		print('rateperyear')
	# 		print(r)
	# 		EmapAtDuration=self.windowedEmaps[dataindex]

	# 		mean=logfits[0][1]
	# 		std=logfits[0][2]
	# 		loc=logfits[0][3]
	# 		ratio=logfits[0][4]
	# 		refFieldLog=fits.logcdfxfromy(r/ratio,mean,std,loc)
	# 		print('refFieldLog')
	# 		print(refFieldLog)

	# 		for j in range(0,int(np.floor(len(EmapAtDuration))/2)):
	# 			durationindex = j*2
	# 			print('durationratios[ij')
	# 			print(durationratios[j])
	# 			duration = EmapAtDuration[durationindex]
	# 			alldurations.append(duration)
	# 			fractionover=self.calcFractionOverTemp(refFieldLog*durationratios[i],duration)
	# 			allfractionsover.append(fractionover)
	# 		plt.plot(alldurations,allfractionsover,lw=1,label = "Once per "+str(1/r)+" year storm")
	# 	plt.title('HVT Exceeding Short Term Temp Limit')
	# 	plt.ylabel('Fraction of HVT Population')
	# 	plt.xlabel('Duration (s)')
	# 	plt.legend()
	# 	plt.show()


	#plots the relative propensity for HVT to as a function of E field and duration
	def plotOverheatByDuration(self):

		durationratios=self.averagedRatios
		i=0
		bucketsize=2
		nperbucket=10
		nbuckets=5
		meanatduration=[]
		allbuckets=np.array([])
		alldurations=[]
		dataindex=i*2+1
		EmapAtDuration=self.windowedEmaps[dataindex]
		ndurations=int(np.floor(len(EmapAtDuration))/2)
		allbuckets=np.zeros([nbuckets,ndurations])
		allfractionover=np.zeros([nbuckets,ndurations])
		for j in range(0,ndurations):
			durationindex = j*2
			duration = EmapAtDuration[durationindex]
			alldurations.append(duration)
			
			fractionovers=np.array([])
			buckets=np.array([])
			for l in range(0,nbuckets):
				bucket=np.array([])
				bucketfractions=np.array([])
				for k in range(0,nperbucket+1):
					E=bucketsize*(l+k/nperbucket)
					bucket=np.append(bucket,E)
					fractionover=self.calcFractionOverTemp(E*durationratios[j],duration)
					bucketfractions=np.append(bucketfractions,fractionover)
				allfractionover[l][j]=np.mean(bucketfractions)
				allbuckets[l][j]=np.mean(bucket)
		
		plt.figure()
		for j in range(0,nbuckets):
			plt.plot(alldurations,allfractionover[j,:],lw=1,label = str(allbuckets[j,1])+" V/km, 60 second E field")

		plt.title('HVT Exceeding Short Term Temp Limit')
		plt.ylabel('Fraction of HVT Population')
		plt.xlabel('Duration (s)')
		plt.legend()
		plt.savefig(Params.figuresDir+'OverheatByDuration.png')

		plt.show()

	

	# calculate fraction of population that exceeds its temperature limit for every HVT type
	def calcFractionOverTemp(self,E,duration):
		fraction=0
		for i in range(0,len(self.HVTnames)):
			temprise0=self.temprise0[i]
			temprise10=self.temprise10[i]
			temprise20=self.temprise20[i]
			temprise40=self.temprise40[i]
			temprise50=self.temprise50[i]
			temprise100=self.temprise100[i]
			temprise200=self.temprise200[i]
			GICperE=self.GICperE[i]
			maxextrapolation=100000
			Etransformer=np.array([0,10,20,40,50,100,200,maxextrapolation])/GICperE
			#linearly extrapolate last two datapoints for higher values
			longtermlimit=np.array([temprise0,temprise10,temprise20,temprise40,temprise50,temprise100,temprise200,temprise200+(temprise200-temprise100)*((maxextrapolation-200)/100)])
			# plt.figure()
			# plt.plot(Etransformer,longtermlimit)
			# plt.show()

			topoiltemp=90 #C
			temp=np.interp(E,Etransformer,longtermlimit)*(1-np.exp(-duration/self.tau[i]))+topoiltemp
			# temp=(E*self.alpha[i]+self.beta[i])*(1-np.exp(-duration/self.tau[i]))+90
			# print('temp')
			# print(temp)
			# temps = temps + [temp]
			if(temp>140):
				fraction = fraction + self.pop40plus[i]
			if(temp>160):
				fraction = fraction + self.pop25to40[i]
			if(temp>180):
				fraction = fraction + self.pop0to25[i]
		#we divide the results by 2, because we assume half of HVT are design 1, which will never overheat
		return fraction/2

	def calcPopulationAffected(self):
		#percent spares
		spares=.15

		pdata=rasterio.open(Params.popDensity15min)
		ldata=rasterio.open(Params.landArea15min)
		#latitude resolution in degrees divided by 30 minute resolution of population density data gives us the number of points below and above to sum the population for this transformer location.
		populationLatN=np.int(np.floor(Params.latituderes/.25))
		populationLongN=np.int(np.floor(Params.longituderes/.25))
		pArr=pdata.read(1)
		lArr=ldata.read(1)
		pArrZeroed = np.where(pArr<0, 0, pArr)
		lArrZeroed = np.where(lArr<0, 0, lArr)
		totPop=np.multiply(pArrZeroed,lArrZeroed)
		overheatdata=np.load(Params.overheatDataDir,allow_pickle=True)
		for i in range(0,int(np.floor(len(overheatdata))/4)):
			totalPopCELE=0
			rateindex=i*4
			latitudeindex=i*4+1
			longitudeindex=i*4+2
			dataindex=i*4+3
			r=overheatdata[rateindex]
			data=overheatdata[dataindex]
			longitudes=overheatdata[longitudeindex]
			latitudes=overheatdata[latitudeindex]
			populationaffected=[]
			for j in range(0,len(data)):
				longitude=longitudes[j]
				latitude=latitudes[j]
				overheat=data[j]
				if(overheat>spares):
					indexlat,indexlong=pdata.index(longitude,latitude)

					population=0
					#sum over the populations near the grid point of interest 
					for k in range(0,populationLatN):
						newLatIndex=indexlat-populationLatN//2+k
						for l in range(0,populationLongN):
							newLongIndex=indexlong-populationLongN//2+l
							population=population+totPop[newLatIndex,newLongIndex]
							totalPopCELE=totalPopCELE+totPop[newLatIndex,newLongIndex]
					populationaffected.append(population)
				else:
					populationaffected.append(0)


			#Catastrophic Electricity Loss assumed if spares<overheated transformers	
			regionsCELE= np.where(np.array(populationaffected)>0, 1, 0)
			df=pd.DataFrame({'longs':np.array(longitudes),'lats':np.array(latitudes),'pop':np.array(populationaffected),'regionsCELE':np.array(regionsCELE)})
			crs={'init':'epsg:3857'}
			geometry=[Point(xy) for xy in zip(df['longs'],df['lats'])]
			geo_df=gpd.GeoDataFrame(df,crs=crs,geometry=geometry)
			# ax=gplt.kdeplot(geo_df,column='over',cmap='Reds',shade=True)
			ax=geo_df.plot(column='regionsCELE')
			# ax.
			# cbr = fig.colorbar(sm, cax=cax,)
			# cbr.ax.tick_params() 
			# cbr.ax.sexzt_yticklabels(['{:.0f}'.format(x) for x in np.arange(.01, .4+.001, .001)], fontsize=16, weight='bold')
			world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
			pp=gplt.polyplot(world,ax=ax,zorder=1)
			
			plt.title('Predicted Regions HVT Overheated Exceed Spares \n'+str(1/r)+' Year Storm, Assuming '+str(100*spares)+'% Spares')
			# gplt.kdeplot(geo_df,column='E',)
			print(Params.figuresDir+'/RegionsCELE/regionsAffected'+str(r)+ 'peryear.png')
			plt.savefig(Params.figuresDir+'/RegionsCELE/regionsAffected'+str(r)+ 'peryear.png')
			plt.show()
			print('world population')
			print(np.sum(totPop))
			print('world population CELE, 1 in '+str(1/r)+' year storm')
			print(totalPopCELE)
			print('fraction world population CELE, 1 in '+str(1/r)+' year storm')
			print(totalPopCELE/np.sum(totPop))
