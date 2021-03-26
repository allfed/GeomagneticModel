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

	def loadEfieldMaps(self):
		loaddirectory=Params.globalEfieldData
		self.windowedEmaps=np.load(loaddirectory,allow_pickle=True)

	def calcOverheatMap(self):
		print('calculating overheating')
		#for each rate per year calculated in earth model
		for i in range(0,int(np.floor(len(self.windowedEmaps))/2)):
			meanoverheat=[]
			durations=[]
			maxoverheatfractions=[]
			maxallfractions=0
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
				durations.append(duration)
				print('duration')
				print(duration)
				Emap = EmapAtDuration[Emapindex]
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
						fractionover=self.calcFractionOverTemp(E,duration)
						
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
				meanoverheat.append(np.mean(allfractionsover))
				print('mean overheat fraction')
				print(np.mean(allfractionsover))
				print('max overheat fraction')
				print(np.max(allfractionsover))
				
				[xlog,ylog]=fits.binlognormaldist(allfractionsover,[],-1)
				[xall,yall]=fits.binnormaldist(allfractionsover,[],-1)

				plt.figure()
				plt.plot(xall,yall)
				plt.title('distributions of overheat fractions,linear scale')
				plt.show()

				df=pd.DataFrame({'longs':np.array(longitudes),'lats':np.array(latitudes),'over':np.array(allfractionsover)})
				crs={'init':'epsg:3857'}
				geometry=[Point(xy) for xy in zip(df['longs'],df['lats'])]
				geo_df=gpd.GeoDataFrame(df,crs=crs,geometry=geometry)
				# ax=gplt.kdeplot(geo_df,column='over',cmap='Reds',shade=True)
				ax=geo_df.plot(column='over',legend=True, cmap='rainbow',vmin=0.01, vmax=.4)#cnp.max(allfractionsover))
				# norm=colors.LogNorm(vmin=0.01, vmax=1))

				world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
				gplt.polyplot(world,ax=ax,zorder=1)
				plt.title('percent transformers overheat, '+str(r)+'per year, '+str(duration))
				# gplt.kdeplot(geo_df,column='E',)
				plt.savefig(Params.overheatMapsDir+'Overheat'+str(r)+'perYearWindow'+str(duration)+'s.png')
				plt.show()

				# plt.figure(1)
				# plt.imshow(np.flipud(thisdurationoverheatmap[:,:]), cmap='hot', interpolation='nearest',vmin=0,vmax=1)
				# plt.title('Overheat fractions '+str(r)+' per year, '+str(duration)+'s')
				# cbar = plt.colorbar()
				# cbar.set_label('Fraction overheat (1=100%)')
				# plt.show()

		plt.figure()
		plt.plot(durations,meanoverheat)
		plt.show()

		df=pd.DataFrame({'longs':np.array(longitudes),'lats':np.array(latitudes),'over':np.array(maxoverheatfractions)})
		crs={'init':'epsg:3857'}
		geometry=[Point(xy) for xy in zip(df['longs'],df['lats'])]
		geo_df=gpd.GeoDataFrame(df,crs=crs,geometry=geometry)
		# ax=gplt.kdeplot(geo_df,column='over',cmap='Reds',shade=True)
		ax=geo_df.plot(column='over',legend=True, cmap='rainbow',vmin=0.01, vmax=np.max(maxoverheatfractions))
		# norm=colors.LogNorm(vmin=0.01, vmax=1))

		world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
		gplt.polyplot(world,ax=ax,zorder=1)
		plt.title('percent transformers overheat, '+str(r)+'per year, '+str(duration))
		plt.savefig(Params.overheatMapsDir+'Overheat'+str(r)+'perYearMaxfromallWindows.png')
		# gplt.kdeplot(geo_df,column='E',)
		plt.show()

	
		# plt.figure(1)
		# plt.imshow(np.flipud(overheatFractionMap[8:,:]), cmap='hot', interpolation='nearest',vmin=0,vmax=1)
		# plt.title('Overheat fractions '+str(r)+' per year, max of all durations')
		# cbar = plt.colorbar()
		# cbar.set_label('Fraction overheat (1=100%)')
		# plt.show()


		# self.windowedEmaps = self.overheatmaps + [r,windowedEmapsatrate]

		# np.save(Params.globalEfieldData,self.windowedEmaps)


	# calculate fraction of population that exceeds its temperature limit for every transformer type
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
		return fraction