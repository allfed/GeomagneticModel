import numpy as np
#calculate thermal rise from duration and GIC level
import Params
import h5py
import matplotlib.pyplot as plt
import geopandas
import pandas as pd
import matplotlib.colors as colors
import geopandas as gpd
from geopandas.tools import sjoin
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
		self.voltageClass=Params.voltageClass
		self.voltage=Params.voltage
		self.phase=Params.phase

		self.pop25to40=Params.pop25to40
		self.pop0to25=Params.pop0to25
		self.pop40plus=Params.pop40plus
		self.tau=Params.tau
		self.temperatures=np.zeros(len(Params.tau))
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

	#calculate the temperatures for each transformer category tie bar
	#also show the mean temperature of the transformers for each duration
	def calcTemperatureMap(self):
		#whether to plot each duration
		plotintermediates=False
		plotAllrates=False
		print('calculating temperatures')
		#for each rate per year calculated in earth model
		allmaxtemps=np.array([])
		for i in range(0,int(np.floor(len(self.windowedEmaps))/2)):
			alldurations=[]
			allmeanoverheats=[]
			meanoverheat=np.array([])
			durations=np.array([])
			maxtemperatures=[]
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

			for j in range(0,int(np.floor(len(EmapAtDuration))/2)):
				durationindex = j*2
				Emapindex = j*2+1

				duration = EmapAtDuration[durationindex]
				durations=np.append(durations,duration)
				print('duration')
				print(duration)
				Emap = EmapAtDuration[1]#Emapindex]
				byDurationTempsTie1Map = np.zeros((lenlong, lenlat))
				tempsTie1Map = np.zeros((lenlong,lenlat))
				latitudes=[]
				longitudes=[]
				alltemps=[]
				alltempsTie1=[]
				i=0
				for latkey in range(0,len(self.latitudes)):
					latitude = self.latitudes[latkey]
					for longkey in range(0,len(self.longitudes)):
						longitude = self.longitudes[longkey]
						magcoords = self.geotomag(latitude, longitude)
						maglat=magcoords.data[0][1]
				
						if(maglat>-70 and maglat <80):
							E=Emap[latkey,longkey]
							temps=self.calcTemps(self.averagedRatios[j]*E,duration)
							alltemps.append(np.array(temps))
							alltempsTie1.append(temps[0])
							latitudes.append(latitude)
							longitudes.append(longitude)
							byDurationTempsTie1Map[latkey][longkey]=temps[0]
							tempsTie1Map[latkey][longkey]=max(tempsTie1Map[latkey][longkey],temps[0])
							if(j==0):
								maxtemperatures.append(np.array(temps))
							else:
								for k in range(0,len(temps)):
									maxtemperatures[i][k]=max(maxtemperatures[i][k],temps[k])
							i=i+1
				
				#plot the tie bar 1 temp for each durations at this rate per year
				if(plotintermediates):
					df=pd.DataFrame({'longs':np.array(longitudes),'lats':np.array(latitudes),'over':np.array(alltempsTie1)})
					crs={'init':'epsg:3857'}
					geometry=[Point(xy) for xy in zip(df['longs'],df['lats'])]
					geo_df=gpd.GeoDataFrame(df,crs=crs,geometry=geometry)
					ax=geo_df.plot(column='over',legend=True,legend_kwds={'label': 'Peak Temperature (C)','orientation': "horizontal"}, cmap='viridis',vmin=np.min(alltempsTie1),vmax=np.max(alltempsTie1))
					world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
					pp=gplt.polyplot(world,ax=ax,zorder=1)
					
					plt.title('526kV Single Phase GSU Tie Bar Peak Temp \n'+str(duration) +' second Geoelectric Field\n 1 in '+str(1/r)+' Year Storm')
					print(Params.overheatMapsDir+'ByDuration/'+str(1/r)+'/Overheat'+str(r)+'perYearWindow'+str(duration)+'s.png')
					plt.savefig(Params.overheatMapsDir+'ByDuration/'+str(1/r)+'/TempTie1_'+str(r)+'perYearWindow'+str(duration)+'s.png')

			#plot the tie bar 1 temp, max of all durations at this rate per year
			if(plotAllrates):
				alldurations.append(durations)
				allmeanoverheats.append(meanoverheat)
				print('maximum temperatures')
				print(np.max(maxtemperatures))
				tie1MaxTemps=np.array(maxtemperatures)[:,0]
				df=pd.DataFrame({'longs':np.array(longitudes),'lats':np.array(latitudes),'over':tie1MaxTemps})
				crs={'init':'epsg:3857'}
				geometry=[Point(xy) for xy in zip(df['longs'],df['lats'])]
				geo_df=gpd.GeoDataFrame(df,crs=crs,geometry=geometry)
				# ax=gplt.kdeplot(geo_df,column='over',cmap='Reds',shade=True)
				ax=geo_df.plot(column='over',legend=True,legend_kwds={'label': 'Peak Temperature (C)','orientation': "horizontal"}, cmap='viridis',vmin=np.min(tie1MaxTemps), vmax=np.max(tie1MaxTemps))

				world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
				gplt.polyplot(world,ax=ax,zorder=1)
				plt.title('526kV Single Phase GSU Tie Bar Peak Temp \n1 in '+str(1/r)+' Year Storm')
				plt.savefig(Params.overheatMapsDir+str(r)+'perYearMaxfromallWindows.png')
				plt.show()
			
			allmaxtemps=np.append(allmaxtemps,[r,np.array(latitudes),np.array(longitudes),np.array(maxtemperatures)])
		np.save(Params.tempDataDir,allmaxtemps)

	#calculates a pandas dataframe that saves the overheat fractions of the various categories of transformers along with true/false for whether the category exceeds the expected spare transformers. 
	#note all of the fractions are out of the total population!
	def calcOverheatMap(self):
		print('calculating overheat map')
		tempdata=np.load(Params.tempDataDir,allow_pickle=True)
		
		dataArr=np.array([])
		for i in range(0,int(np.floor(len(tempdata))/4)):
			rateindex=i*4
			latitudeindex=i*4+1
			longitudeindex=i*4+2
			dataindex=i*4+3
			r=tempdata[rateindex]
			data=tempdata[dataindex]
			longitudes=tempdata[longitudeindex]
			latitudes=tempdata[latitudeindex]
			stats=[]
			T230kvfractions=[]
			T345kvfractions=[]
			T500kvfractions=[]
			T765kvfractions=[]
			singlefractions=[]
			threefractions=[]
			CELEsingles=[]
			CELEthrees=[]
			isCELEs=[]
			for j in range(0,len(data)):
				longitude=longitudes[j]
				latitude=latitudes[j]
				# print('data')
				# print(data)
				# print('data1')
				# print(data[1])
				# print('lendata')
				# print(len(data))
				# print('j')
				# print(j)
				temperatures=data[j]
				T230kvfraction=0
				T345kvfraction=0
				T500kvfraction=0
				T765kvfraction=0
				singlefraction=0
				threefraction=0
				CELEsingle=False
				CELEthree=False
				isCELE=False
				for k in range(0,len(self.HVTnames)):
					# print('k')
					# print(k)
					# print('temperatures')
					# print('temperatures[0]')
					# print(temperatures)
					# print(temperatures[0])
					# print
					temp=temperatures[k]

					if(self.phase[k]=='single'):
						if(temp>140):
							singlefraction = singlefraction + self.pop40plus[i]
						if(temp>160):
							singlefraction = singlefraction + self.pop25to40[i]
						if(temp>180):
							singlefraction = singlefraction + self.pop0to25[i]

					if(self.voltageClass[k]=='345'):
						if(temp>140):
							T345kvfraction = T345kvfraction + self.pop40plus[i]
						if(temp>160):
							T345kvfraction = T345kvfraction + self.pop25to40[i]
						if(temp>180):
							T345kvfraction = T345kvfraction + self.pop0to25[i]
					if(self.voltageClass[k]=='500'):
						if(temp>140):
							T500kvfraction = T500kvfraction + self.pop40plus[i]
						if(temp>160):
							T500kvfraction = T500kvfraction + self.pop25to40[i]
						if(temp>180):
							T500kvfraction = T500kvfraction + self.pop0to25[i]
					if(self.voltageClass[k]=='765'):
						if(temp>140):
							T765kvfraction = T765kvfraction + self.pop40plus[i]
						if(temp>160):
							T765kvfraction = T765kvfraction + self.pop25to40[i]
						if(temp>180):
							T765kvfraction = T765kvfraction + self.pop0to25[i]
				if(singlefraction>.095):#1/3 of the 28% by population single transformers (see voltage and fractional population tab in transformers spreadsheet)
					CELEsingle=True
				if(threefraction>.235):#appropriate fraction such that total replacement transformers is 15% total spares available (according to at least one source in the literature, see paper)
					CELEthree=True
				isCELE= CELEsingle or CELEthree

				T230kvfractions.append(T230kvfraction)
				T345kvfractions.append(T345kvfraction)
				T500kvfractions.append(T500kvfraction)
				T765kvfractions.append(T765kvfraction)
				singlefractions.append(singlefraction)
				threefractions.append(threefraction)
				CELEsingles.append(CELEsingle)
				CELEthrees.append(CELEthree)
				isCELEs.append(isCELE)
			df=pd.DataFrame({'lats':np.array(latitudes),\
				'longs':np.array(longitudes),\
				'T230kvfractions':np.array(T230kvfractions),\
				'T345kvfractions':np.array(T345kvfractions),\
				'T500kvfractions':np.array(T500kvfractions),\
				'T765kvfractions':np.array(T765kvfractions),\
				'singlefractions':np.array(singlefractions),\
				'threefractions':np.array(threefractions),\
				'CELEsingles':np.array(CELEsingles),\
				'CELEthrees':np.array(CELEthrees),\
				'isCELEs':np.array(isCELEs)})

			dataArr=np.append(dataArr,[r,df])
		np.save(Params.overheatDataDir,dataArr)

	#calculate the temperatures for each transformer category tie bar
	# def calcOverheatMaps(self):
	# 	#whether to plot each duration
	# 	plotintermediates=True
	# 	print('calculating overheating')
	# 	#for each rate per year calculated in earth model
	# 	alltemps=np.array([])
	# 	for i in range(0,int(np.floor(len(self.windowedEmaps))/2)):
	# 		alldurations=[]
	# 		durations=np.array([])
	# 		maxtemperatures=[]
	# 		rateindex=i*2
	# 		dataindex=i*2+1
	# 		r=self.windowedEmaps[rateindex]
	# 		print('rateperyear')
	# 		print(r)
	# 		EmapAtDuration=self.windowedEmaps[dataindex]
	# 		latitudes=self.latitudes
	# 		longitudes=self.longitudes

	# 		lenlat=len(latitudes)
	# 		lenlong=len(longitudes)

	# 		for j in range(0,int(np.floor(len(EmapAtDuration))/2)):
	# 			durationindex = j*2
	# 			Emapindex = j*2+1

	# 			duration = EmapAtDuration[durationindex]
	# 			durations=np.append(durations,duration)
	# 			print('duration')
	# 			print(duration)
	# 			Emap = EmapAtDuration[1]#Emapindex]
	# 			thisDurationTempMap = np.zeros((lenlong, lenlat))
	# 			tempMap = np.zeros((lenlong,lenlat))
	# 			latitudes=[]
	# 			longitudes=[]
	# 			allfractionsover=[]
	# 			allfractions1over=[]
	# 			allfractions3over=[]
	# 			i=0
	# 			for latkey in range(0,len(self.latitudes)):
	# 				latitude = self.latitudes[latkey]
	# 				for longkey in range(0,len(self.longitudes)):
	# 					longitude = self.longitudes[longkey]
	# 					E=Emap[latkey,longkey]
	# 					temps=self.calcTemps(self.averagedRatios[j]*E,duration)
	# 					magcoords = self.geotomag(latitude, longitude)
	# 					maglat=magcoords.data[0][1]
				
	# 					if(maglat>-70 and maglat <80):
	# 						allfractionsover.append(fractionover)
	# 						allfractions1over.append(fractionSingleOver)

	# 						latitudes.append(latitude)
	# 						longitudes.append(longitude)
	# 						thisdurationoverheatmap[latkey][longkey]=fractionover
	# 						overheatFractionMap[latkey][longkey]=max(overheatFractionMap[latkey][longkey],fractionover)
	# 						if(j==0):
	# 							maxtemperatures.append(np.array(temperatures))
	# 						else:
	# 							for k in temperatures:
	# 								print(maxtemperatures[i])
	# 								quit()
	# 								maxtemperatures[k]=max(maxtemperatures[k],temperatures[k])
	# 						i=i+1
				
	# 		alltemps=np.append(alltemps,[r,np.array(latitudes),np.array(longitudes),np.array(maxtemperatures)])
	# 		print('np.len(alltemps)')
	# 		print(len(alltemps))
	# 		alldurations.append(durations)
	# 		allmeanoverheats.append(meanoverheat)
	# 		# plt.figure()
	# 		# plt.plot(durations,meanoverheat)
	# 		# plt.title('HVT Exceeding Short Term Temp Limit')
	# 		# plt.ylabel('Fraction of HVT Population')
	# 		# plt.xlabel('Duration (s)')
	# 		# plt.show()
	# 		print('maximum temperatures')
	# 		print(np.max(maxtemperatures))
	# 	np.save(Params.overheatDataDir,alloverheatfractions)
	# 	# plt.figure(1)
	# 	# plt.imshow(np.flipud(overheatFractionMap[8:,:]), cmap='hot', interpolation='nearest',vmin=0,vmax=1)
	# 	# plt.title('Overheat fractions '+str(r)+' per year, max of all durations')
	# 	# cbar = plt.colorbar()
	# 	# cbar.set_label('Fraction overheat (1=100%)')
	# 	# plt.show()


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


	# calculate fraction of population that exceeds its temperature limit for every HVT type
	#category can be single ('single'), three phase ('three'), or leave empty to calculate sum of both, giving total overheat fraction ('')
	def calcTemps(self,E,duration):
		ntransformers=len(self.HVTnames)
		temperatures=[0]*ntransformers
		for i in range(0,ntransformers):
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
			temperatures[i]=temp
		return temperatures

	def calcFractionOver(temperatures):
		ntransformers=len(self.HVTnames)
		fraction=0
		for i in range(0,ntransformers):		
			if(temp>140):
				fraction = fraction + self.pop40plus[i]
			if(temp>160):
				fraction = fraction + self.pop25to40[i]
			if(temp>180):
				fraction = fraction + self.pop0to25[i]
		#we divide the results by 2, because we assume half of HVT are design 1, which will never overheat
		return fraction/2


	# For estimating loss of electricity we use the net consumption by country, and estimate the fraction of people in that country that would lose electricity.
	# The algorithm is:
	#     sum over all countries:
	#         (country consumption)*(fraction population loses electricity in country)
	def calcElectricityAffected(self):

		#import electricity by nation data

		rawimport=pd.read_csv('Data/SmallData/ElectricityByNation/EleByCountry.csv',skiprows=4,header=0)
		eleByCountry=pd.DataFrame({'countrycode':rawimport.iloc[1::3,0].values,'country':rawimport.iloc[0::3,1].values,'consumption':rawimport.iloc[2::3,2].values,'generation':rawimport.iloc[1::3,2].values})
		eleByCountry['countrycode']=eleByCountry['countrycode'].str.slice(start=14,stop=17)
		codearr=rawimport.iloc[1::3,0].str.split(pat='-').values
		
		codes=['']*len(codearr)
		for i in range(0,len(codearr)):
			ar=codearr[i]
			if(len(ar)<2):
				codes[i]==''
			else:
				codes[i]=ar[2]

		eleByCountry['countrycode']=np.array(codes)

		eleByCountry['consumption']=eleByCountry['consumption'].replace('--','0').astype(float)
		eleByCountry['generation']=eleByCountry['generation'].replace('--','0').astype(float)

		#iterate through the population data, and use the fraction of population in the country to determine loss by net consumption 
		allPopData=np.load(Params.popCELEdir,allow_pickle=True)
		for k in range(0,len(allPopData)):
			r=allPopData[k][0]
			popCELE=allPopData[k][1]
			#import the country boundaries, so we can see which country coordinates fall into
			world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
			countryshapes=world['geometry'].values
			dfpoint=pd.DataFrame({'longs':popCELE['longs'].values,'lats':popCELE['lats'].values,'popCELE':popCELE['pop'].values})
			crs={'init':'epsg:3857'}
			geometry=[Point(xy) for xy in zip(dfpoint['longs'],dfpoint['lats'])]
			geo_df=gpd.GeoDataFrame(dfpoint,crs=crs,geometry=geometry)
			world['countrygeometry']=world['geometry']
			#get groups of countries which match the coordinates for population losing electricity
			pointInPolys = sjoin(geo_df, world, how='left')

			grouped = pointInPolys.groupby('index_right')

			eleLostArr=[]
			fractionEleLostArr=[]
			eleCountry=[]
			countrycode=[]
			geometries=[]

			#for each country
			netpop=0
			for key, values in grouped:
				code=values['iso_a3'].values[0]
				eleByCountry['countrycode'].str.match(code)
				countryele=[]
				for i in range(0,len(eleByCountry['countrycode'])):
					if(eleByCountry['countrycode'].values[i]==code):
						countryele=eleByCountry.iloc[i]
						break
				if(len(countryele)==0):
					continue

				countrypop=values['pop_est'].values[0]
				popCELE=values['popCELE'].sum()
				eleconsumption=countryele['consumption']
				# print('eleconsumption')
				# print(eleconsumption)
				# print('popCELE')
				# print(popCELE)
				# print('countrypop')
				# print(countrypop)
				eleLost=eleconsumption*(popCELE/countrypop)
				countrycode.append(code)
				eleCountry.append(eleconsumption)
				eleLostArr.append(eleLost)
				fractionEleLostArr.append(eleLost/eleconsumption)
				geometries.append(values['countrygeometry'].values[0])
			worldElectricity=gpd.GeoDataFrame({'iso_a3':np.array(countrycode),'geometry':np.array(geometries),'fraction':np.array(fractionEleLostArr)})

			fig, ax = plt.subplots(1, 1)
			worldElectricity.plot(column='fraction', ax=ax, legend=True,vmin=0,vmax=1)
			
			plt.title('Predicted Fraction Electricity By Country \n One in '+str(1/r)+' Year Storm')
			plt.show()

			totalEleLoss=np.sum(eleLostArr)
			totalEle=np.sum(eleCountry)
			print('totalElectricity (GW)')
			print(totalEle)
			print('totalEleLoss (GW)')
			print(totalEleLoss)
			print('fraction overall')
			print(totalEleLoss/totalEle)


	#high resolution population estimates (15 minute)
	def calcPopulationAffected(self):
		print('calculating population and regions with CELE')
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
		allPopData=[]
		for i in range(0,int(np.floor(len(overheatdata))/2)):
			totalPopCELE=0
			rateindex=i*2
			dataindex=i*2+1
			r=overheatdata[rateindex]
			print('rateperyear'+str(r))
			df=overheatdata[dataindex]
			populationaffected=[]
			for j in range(0,len(df)):
				isCELE=df['isCELEs'].values[j]
				latitude=df['lats'].values[j]
				longitude=df['longs'].values[j]
				if(isCELE):
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
			regionsCELE= np.where(np.array(populationaffected)>0, 1, 0) #population affected is always the full population in the region where spares<overheated transformers. Nonzero=>population in grid cell.
			df=pd.DataFrame({'longs':df['longs'].values,'lats':df['lats'].values,'pop':np.array(populationaffected),'regionsCELE':np.array(regionsCELE)})
			
			crs={'init':'epsg:3857'}
			geometry=[Point(xy) for xy in zip(df['longs'],df['lats'])]
			geo_df=gpd.GeoDataFrame(df,crs=crs,geometry=geometry)
			ax=geo_df.plot(column='regionsCELE')
			world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
			pp=gplt.polyplot(world,ax=ax,zorder=1)
			
			plt.title('Predicted Regions HVT Overheated Exceed Spares \n'+str(1/r)+' Year Storm')
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
			allPopData.append([r,df])
		np.save(Params.popCELEdir,allPopData)