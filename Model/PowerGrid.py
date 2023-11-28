import os
import glob
import numpy as np

# calculate thermal rise from duration and GIC level
import Params
import h5py
import matplotlib.pyplot as plt
import geopandas
import pandas as pd
import matplotlib.colors as colors
import geopandas as gpd
from geopandas.tools import sjoin
import geoplot as gplt
import shapely.geometry
from shapely.geometry import Point
import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import LogFormatter
import matplotlib.ticker as mticker
import rasterio
import rasterio as rio
import rasterio.features
import rasterio.warp
from Model.Network import Network
from matplotlib.colors import ListedColormap
from Plotter import Plotter


class PowerGrid:
    def __init__(self):
        Params.importIfNotAlready()
        self.windowedEmaps = []
        self.HVTnames = Params.HVTnames
        self.GICperE = Params.GICperE

        self.temprise0 = Params.temprise0
        self.temprise10 = Params.temprise10
        self.temprise20 = Params.temprise20
        self.temprise40 = Params.temprise40
        self.temprise50 = Params.temprise50
        self.temprise100 = Params.temprise100
        self.temprise200 = Params.temprise200
        self.voltageClass = Params.voltageClass
        self.voltage = Params.voltage
        self.phase = Params.phase

        self.pop25to40 = Params.pop25to40
        self.pop0to25 = Params.pop0to25
        self.pop40plus = Params.pop40plus
        self.tau = Params.tau
        self.temperatures = np.zeros(len(Params.tau))
        self.latitudes = []
        self.longitudes = []
        self.networks = []

    def setWindowedEmaps(self, earthmodel):
        self.windowedEmaps = earthmodel.windowedEmaps
        self.latitudes = earthmodel.GClatitudes
        self.longitudes = earthmodel.GClongitudes
        self.geotomag = earthmodel.geotomag
        self.averagedRatios = earthmodel.averagedRatios
        self.combinedlogfits = earthmodel.combinedlogfits

    def loadEfieldMaps(self):
        loaddirectory = Params.globalEfieldData
        self.windowedEmaps = np.load(loaddirectory, allow_pickle=True)

    # calculate the temperatures for each transformer category tie bar
    # also show the mean temperature of the transformers for each duration
    def calcTemperatureMap(self):
        # whether to plot each duration
        plotintermediates = True
        plotAllrates = True
        print("calculating temperatures")
        # for each rate per year calculated in earth model
        allmaxtemps = np.array([])
        for i in range(0, int(np.floor(len(self.windowedEmaps)) / 2)):
            alldurations = []
            allmeanoverheats = []
            meanoverheat = np.array([])
            durations = np.array([])
            maxtemperatures = []
            rateindex = i * 2
            dataindex = i * 2 + 1
            r = self.windowedEmaps[rateindex]
            print("rateperyear")
            print(r)
            EmapAtDuration = self.windowedEmaps[dataindex]
            latitudes = self.latitudes
            longitudes = self.longitudes

            lenlat = len(latitudes)
            lenlong = len(longitudes)

            for j in range(0, int(np.floor(len(EmapAtDuration)) / 2)):
                durationindex = j * 2
                Emapindex = j * 2 + 1

                duration = EmapAtDuration[durationindex]
                durations = np.append(durations, duration)
                print("duration")
                print(duration)
                Emap = EmapAtDuration[1]  # Emapindex]
                byDurationTempsTie1Map = np.zeros((lenlong, lenlat))
                tempsTie1Map = np.zeros((lenlong, lenlat))
                latitudes = []
                longitudes = []
                alltemps = []
                alltempsTie1 = []
                i = 0
                for latkey in range(0, len(self.latitudes)):
                    latitude = self.latitudes[latkey]
                    for longkey in range(0, len(self.longitudes)):
                        longitude = self.longitudes[longkey]
                        magcoords = self.geotomag(latitude, longitude)
                        maglat = magcoords.data[0][1]

                        if maglat > -70 and maglat < 80:
                            E = Emap[latkey, longkey]
                            temps = self.calcTempsFromE(
                                self.averagedRatios[j] * E, duration
                            )
                            alltemps.append(np.array(temps))
                            alltempsTie1.append(temps[0])
                            latitudes.append(latitude)
                            longitudes.append(longitude)
                            byDurationTempsTie1Map[latkey][longkey] = temps[0]
                            tempsTie1Map[latkey][longkey] = max(
                                tempsTie1Map[latkey][longkey], temps[0]
                            )
                            if j == 0:
                                maxtemperatures.append(np.array(temps))
                            else:
                                for k in range(0, len(temps)):
                                    maxtemperatures[i][k] = max(
                                        maxtemperatures[i][k], temps[k]
                                    )
                            i = i + 1

                # plot the tie bar 1 temp for each durations at this rate per year
                if plotintermediates:
                    df = pd.DataFrame(
                        {
                            "longs": np.array(longitudes),
                            "lats": np.array(latitudes),
                            "over": np.array(alltempsTie1),
                        }
                    )
                    crs = {"init": "epsg:3857"}
                    geometry = [Point(xy) for xy in zip(df["longs"], df["lats"])]
                    geo_df = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)
                    ax = geo_df.plot(
                        column="over",
                        legend=True,
                        legend_kwds={
                            "label": "Peak Temperature (C)",
                            "orientation": "horizontal",
                        },
                        cmap="viridis",
                        vmin=np.min(alltempsTie1),
                        vmax=np.max(alltempsTie1),
                    )
                    world = geopandas.read_file(
                        geopandas.datasets.get_path("naturalearth_lowres")
                    )
                    pp = gplt.polyplot(world, ax=ax, zorder=1)

                    plt.title(
                        "526kV Single Phase GSU Tie Bar Peak Temp \n"
                        + str(duration)
                        + " second Geoelectric Field\n 1 in "
                        + str(1 / r)
                        + " Year Storm"
                    )
                    print(
                        Params.overheatMapsDir
                        + "ByDuration/"
                        + str(1 / r)
                        + "/Overheat"
                        + str(r)
                        + "perYearWindow"
                        + str(duration)
                        + "s.png"
                    )
                    plt.savefig(
                        Params.overheatMapsDir
                        + "ByDuration/"
                        + str(1 / r)
                        + "/TempTie1_"
                        + str(r)
                        + "perYearWindow"
                        + str(duration)
                        + "s.png"
                    )

            # plot the tie bar 1 temp, max of all durations at this rate per year
            if plotAllrates:
                alldurations.append(durations)
                allmeanoverheats.append(meanoverheat)
                print("maximum temperatures")
                print(np.max(maxtemperatures))
                tie1MaxTemps = np.array(maxtemperatures)[:, 0]
                df = pd.DataFrame(
                    {
                        "longs": np.array(longitudes),
                        "lats": np.array(latitudes),
                        "over": tie1MaxTemps,
                    }
                )
                crs = {"init": "epsg:3857"}
                geometry = [Point(xy) for xy in zip(df["longs"], df["lats"])]
                geo_df = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)
                # ax=gplt.kdeplot(geo_df,column='over',cmap='Reds',shade=True)
                ax = geo_df.plot(
                    column="over",
                    legend=True,
                    legend_kwds={
                        "label": "Peak Temperature (C)",
                        "orientation": "horizontal",
                    },
                    cmap="viridis",
                    vmin=np.min(tie1MaxTemps),
                    vmax=np.max(tie1MaxTemps),
                )

                world = geopandas.read_file(
                    geopandas.datasets.get_path("naturalearth_lowres")
                )
                gplt.polyplot(world, ax=ax, zorder=1)
                plt.title(
                    "526kV Single Phase GSU Tie Bar Peak Temp \n1 in "
                    + str(1 / r)
                    + " Year Storm"
                )
                plt.savefig(
                    Params.overheatMapsDir + str(r) + "perYearMaxfromallWindows.png"
                )
                plt.show()

            allmaxtemps = np.append(
                allmaxtemps,
                [
                    r,
                    np.array(latitudes),
                    np.array(longitudes),
                    np.array(maxtemperatures),
                ],
            )
        np.save(Params.tempDataDir, allmaxtemps)

    # calculates a pandas dataframe that saves the overheat fractions of the various categories of transformers along with true/false for whether the category exceeds the expected spare transformers.
    # note all of the fractions are out of the total population!
    def calcOverheatMap(self):
        print("calculating overheat map")
        tempdata = np.load(Params.tempDataDir, allow_pickle=True)

        dataArr = np.array([])
        for i in range(0, int(np.floor(len(tempdata)) / 4)):
            rateindex = i * 4
            latitudeindex = i * 4 + 1
            longitudeindex = i * 4 + 2
            dataindex = i * 4 + 3
            r = tempdata[rateindex]
            data = tempdata[dataindex]
            longitudes = tempdata[longitudeindex]
            latitudes = tempdata[latitudeindex]
            stats = []
            T230kvfractions = []
            T345kvfractions = []
            T500kvfractions = []
            T765kvfractions = []
            singlefractions = []
            threefractions = []
            CELEsingles = []
            CELEthrees = []
            isCELEs = []
            for j in range(0, len(data)):
                longitude = longitudes[j]
                latitude = latitudes[j]
                # print('data')
                # print(data)
                # print('data1')
                # print(data[1])
                # print('lendata')
                # print(len(data))
                # print('j')
                # print(j)
                temperatures = data[j]
                T230kvfraction = 0
                T345kvfraction = 0
                T500kvfraction = 0
                T765kvfraction = 0
                singlefraction = 0
                threefraction = 0
                CELEsingle = False
                CELEthree = False
                isCELE = False
                for k in range(0, len(self.HVTnames)):
                    # print('k')
                    # print(k)
                    # print('temperatures')
                    # print('temperatures[0]')
                    # print(temperatures)
                    # print(temperatures[0])
                    # print
                    temp = temperatures[k]

                    if self.phase[k] == "single":
                        if temp > 140:
                            singlefraction = singlefraction + self.pop40plus[i]
                        if temp > 160:
                            singlefraction = singlefraction + self.pop25to40[i]
                        if temp > 180:
                            singlefraction = singlefraction + self.pop0to25[i]

                    if self.voltageClass[k] == "345":
                        if temp > 140:
                            T345kvfraction = T345kvfraction + self.pop40plus[i]
                        if temp > 160:
                            T345kvfraction = T345kvfraction + self.pop25to40[i]
                        if temp > 180:
                            T345kvfraction = T345kvfraction + self.pop0to25[i]
                    if self.voltageClass[k] == "500":
                        if temp > 140:
                            T500kvfraction = T500kvfraction + self.pop40plus[i]
                        if temp > 160:
                            T500kvfraction = T500kvfraction + self.pop25to40[i]
                        if temp > 180:
                            T500kvfraction = T500kvfraction + self.pop0to25[i]
                    if self.voltageClass[k] == "765":
                        if temp > 140:
                            T765kvfraction = T765kvfraction + self.pop40plus[i]
                        if temp > 160:
                            T765kvfraction = T765kvfraction + self.pop25to40[i]
                        if temp > 180:
                            T765kvfraction = T765kvfraction + self.pop0to25[i]
                if (
                    singlefraction > 0.095
                ):  # 1/3 of the 28% by population single transformers (see voltage and fractional population tab in transformers spreadsheet)
                    CELEsingle = True
                if (
                    threefraction > 0.235
                ):  # appropriate fraction such that total replacement transformers is 15% total spares available (according to at least one source in the literature, see paper)
                    CELEthree = True
                isCELE = CELEsingle or CELEthree

                T230kvfractions.append(T230kvfraction)
                T345kvfractions.append(T345kvfraction)
                T500kvfractions.append(T500kvfraction)
                T765kvfractions.append(T765kvfraction)
                singlefractions.append(singlefraction)
                threefractions.append(threefraction)
                CELEsingles.append(CELEsingle)
                CELEthrees.append(CELEthree)
                isCELEs.append(isCELE)
            df = pd.DataFrame(
                {
                    "lats": np.array(latitudes),
                    "longs": np.array(longitudes),
                    "T230kvfractions": np.array(T230kvfractions),
                    "T345kvfractions": np.array(T345kvfractions),
                    "T500kvfractions": np.array(T500kvfractions),
                    "T765kvfractions": np.array(T765kvfractions),
                    "singlefractions": np.array(singlefractions),
                    "threefractions": np.array(threefractions),
                    "CELEsingles": np.array(CELEsingles),
                    "CELEthrees": np.array(CELEthrees),
                    "isCELEs": np.array(isCELEs),
                }
            )

            dataArr = np.append(dataArr, [r, df])
        np.save(Params.overheatDataDir, dataArr)

    # calculate the temperatures for each transformer category tie bar
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

    # plots the relative propensity for HVT to as a function of E field and duration
    def plotOverheatByDuration(self):
        durationratios = self.averagedRatios
        i = 0
        bucketsize = 2
        nperbucket = 10
        nbuckets = 5
        meanatduration = []
        allbuckets = np.array([])
        alldurations = []
        dataindex = i * 2 + 1
        EmapAtDuration = self.windowedEmaps[dataindex]
        ndurations = int(np.floor(len(EmapAtDuration)) / 2)
        allbuckets = np.zeros([nbuckets, ndurations])
        allfractionover = np.zeros([nbuckets, ndurations])
        for j in range(0, ndurations):
            durationindex = j * 2
            duration = EmapAtDuration[durationindex]
            alldurations.append(duration)

            fractionovers = np.array([])
            buckets = np.array([])
            for l in range(0, nbuckets):
                bucket = np.array([])
                bucketfractions = np.array([])
                for k in range(0, nperbucket + 1):
                    E = bucketsize * (l + k / nperbucket)
                    bucket = np.append(bucket, E)
                    fractionover = self.calcFractionOverTemp(
                        E * durationratios[j], duration
                    )
                    bucketfractions = np.append(bucketfractions, fractionover)
                allfractionover[l][j] = np.mean(bucketfractions)
                allbuckets[l][j] = np.mean(bucket)

        plt.figure()
        for j in range(0, nbuckets):
            plt.plot(
                alldurations,
                allfractionover[j, :],
                lw=1,
                label=str(allbuckets[j, 1]) + " V/km, 60 second E field",
            )

        plt.title("HVT Exceeding Short Term Temp Limit")
        plt.ylabel("Fraction of HVT Population")
        plt.xlabel("Duration (s)")
        plt.legend()
        plt.savefig(Params.figuresDir + "OverheatByDuration.png")

        plt.show()

    def calcFractionOverTemp(self, E, duration):
        fraction = 0
        for i in range(0, len(self.HVTnames)):
            temprise0 = self.temprise0[i]
            temprise10 = self.temprise10[i]
            temprise20 = self.temprise20[i]
            temprise40 = self.temprise40[i]
            temprise50 = self.temprise50[i]
            temprise100 = self.temprise100[i]
            temprise200 = self.temprise200[i]
            GICperE = self.GICperE[i]
            maxextrapolation = 100000
            Etransformer = (
                np.array([0, 10, 20, 40, 50, 100, 200, maxextrapolation]) / GICperE
            )
            # linearly extrapolate last two datapoints for higher values
            longtermlimit = np.array(
                [
                    temprise0,
                    temprise10,
                    temprise20,
                    temprise40,
                    temprise50,
                    temprise100,
                    temprise200,
                    temprise200
                    + (temprise200 - temprise100) * ((maxextrapolation - 200) / 100),
                ]
            )
            # plt.figure()
            # plt.plot(Etransformer,longtermlimit)
            # plt.show()

            topoiltemp = 90  # C
            temp = (
                np.interp(E, Etransformer, longtermlimit)
                * (1 - np.exp(-duration / self.tau[i]))
                + topoiltemp
            )
            # temp=(E*self.alpha[i]+self.beta[i])*(1-np.exp(-duration/self.tau[i]))+90
            # print('temp')
            # print(temp)
            # temps = temps + [temp]
            if temp > 140:
                fraction = fraction + self.pop40plus[i]
            if temp > 160:
                fraction = fraction + self.pop25to40[i]
            if temp > 180:
                fraction = fraction + self.pop0to25[i]
        # we divide the results by 2, because we assume half of HVT are design 1, which will never overheat
        return fraction / 2

    # #calculate temperatures of all the transformer types from GIC given the duration and the GIC field level.
    # def calcProbabilityDamageFromGIC(self,gic,duration,voltage):
    # 	ntransformers=len(self.HVTnames)
    # 	# if(voltage==float(0)):
    # 	# 	# voltageClass=[230,345,500,735]
    # 	# 	# return 0
    # 	# if(voltage==float(275)):
    # 	# 	voltageClass='230'
    # 	# if(voltage==float(400)):
    # 	# 	voltageClass='345'
    # 	# print('voltageClass')
    # 	# print(voltageClass)
    # 	probabilityDamaged=0
    # 	numberTransformerTypesInVoltageClass=0
    # 	fraction=0
    # 	for i in range(0,ntransformers):
    # 		# numberTransformerTypesInVoltageClass=numberTransformerTypesInVoltageClass+1
    # 		print('i')
    # 		print(i)
    # 		print('gic')
    # 		print(gic)
    # 		print('duration')
    # 		print(duration)

    # 		temp=self.calcTemp(i,gic,duration)
    # 		print('temp')
    # 		print(temp)

    # 		if(temp>140):
    # 			fraction=fraction+Params.pop40plus[i]*.39
    # 		if(temp>160):
    # 			fraction=fraction+Params.pop25to40[i]*.25
    # 		if(temp>180):
    # 			fraction=fraction+Params.pop0to25[i]*.36

    # 	return fraction
    def calcProbabilityDamageFromGIC(self, gic, duration, voltage):
        ntransformers = len(self.HVTnames)
        if voltage < 290:
            voltageClass = "230"
        elif voltage < float(420):
            voltageClass = "345"
        elif voltage < float(620):
            voltageClass = "500"
        else:
            voltageClass = "735"

        # print('voltageClass')
        # print(voltageClass)
        probabilityDamaged = 0
        numberTransformerTypesInVoltageClass = 0
        fraction = 0
        for i in range(0, ntransformers):
            if float(self.voltageClass[i]) == float(voltageClass):
                numberTransformerTypesInVoltageClass = (
                    numberTransformerTypesInVoltageClass + 1
                )
                # print('i')
                # print(i)
                # print('gic')
                # print(gic)
                # print('duration')
                # print(duration)

                temp = self.calcTemp(i, gic, duration)
                # print('temp')
                # print(temp)

                if temp > 140:
                    fraction = fraction + Params.fractionTypeAtVoltage[i] * 0.39
                if temp > 160:
                    fraction = fraction + Params.fractionTypeAtVoltage[i] * 0.25
                if temp > 180:
                    fraction = fraction + Params.fractionTypeAtVoltage[i] * 0.36

        return fraction

    # calculate GIC given an E field value for a given duration
    # uses very simplistic assumptions about GIC levels induced from E fields
    def calcTempsFromE(self, E, duration):
        ntransformers = len(self.HVTnames)
        temperatures = [0] * ntransformers
        for i in range(0, ntransformers):
            gic = self.GICperE[i] * E
            temperatures[i] = self.calcTemp(i, gic, duration)

        return temperatures

    # calculate the temperature of a given transformer index using duration and GIC level
    def calcTemp(self, index, gic, duration):
        temprise0 = self.temprise0[index]
        temprise10 = self.temprise10[index]
        temprise20 = self.temprise20[index]
        temprise40 = self.temprise40[index]
        temprise50 = self.temprise50[index]
        temprise100 = self.temprise100[index]
        temprise200 = self.temprise200[index]
        GICperE = self.GICperE[index]

        maxextrapolation = 100000
        gics = np.array([0, 10, 20, 40, 50, 100, 200, maxextrapolation])
        # linearly extrapolate last two datapoints for higher values
        longtermlimit = np.array(
            [
                temprise0,
                temprise10,
                temprise20,
                temprise40,
                temprise50,
                temprise100,
                temprise200,
                temprise200
                + (temprise200 - temprise100) * ((maxextrapolation - 200) / 100),
            ]
        )
        # plt.figure()
        # plt.plot(Etransformer,longtermlimit)
        # plt.show()

        topoiltemp = 90  # C
        # the calculation of GIC co
        # we multiply the long term limit and the time constant by 2/3, because we assume the HVT are thermally "in the middle" between design 1 and design 2.
        # We multiply the tau value by 3/4 because design 1 is typically 8 minutes, and design 2 is typically 4 minutes, so we choose some average value.
        # This assumption could be improved by actually importing design 2 and performing a more data driven averaging scheme, and even further improved by including a distribution or running a monte carlo with a spread in the range between design 1 and design 2.
        temp = (
            np.interp(gic, gics, longtermlimit * (2 / 3))
            * (1 - np.exp(-duration / ((3 / 4) * self.tau[index])))
            + topoiltemp
        )
        return temp

    def calcFractionOver(temperatures):
        ntransformers = len(self.HVTnames)
        fraction = 0
        for i in range(0, ntransformers):
            if temp > 140:
                fraction = fraction + self.pop40plus[i]
            if temp > 160:
                fraction = fraction + self.pop25to40[i]
            if temp > 180:
                fraction = fraction + self.pop0to25[i]
        return fraction

    # # For estimating loss of electricity we use the net consumption by country, and estimate the fraction of people in that country that would lose electricity.
    # # The algorithm is:
    # #     sum over all countries:
    # #         (country consumption)*(fraction population loses electricity in country)
    # def calcElectricityAffected(self):

    # 	#import electricity by nation data

    # 	rawimport=pd.read_csv('Data/SmallData/ElectricityByNation/EleByCountry.csv',skiprows=4,header=0)
    # 	eleByCountry=pd.DataFrame({'countrycode':rawimport.iloc[1::3,0].values,'country':rawimport.iloc[0::3,1].values,'consumption':rawimport.iloc[2::3,2].values,'generation':rawimport.iloc[1::3,2].values})
    # 	eleByCountry['countrycode']=eleByCountry['countrycode'].str.slice(start=14,stop=17)
    # 	codearr=rawimport.iloc[1::3,0].str.split(pat='-').values

    # 	codes=['']*len(codearr)
    # 	for i in range(0,len(codearr)):
    # 		ar=codearr[i]
    # 		if(len(ar)<2):
    # 			codes[i]==''
    # 		else:
    # 			codes[i]=ar[2]

    # 	eleByCountry['countrycode']=np.array(codes)

    # 	eleByCountry['consumption']=eleByCountry['consumption'].replace('--','0').astype(float)
    # 	eleByCountry['generation']=eleByCountry['generation'].replace('--','0').astype(float)

    # 	#iterate through the population data, and use the fraction of population in the country to determine loss by net consumption
    # 	allPopData=np.load(Params.popCELEdir,allow_pickle=True)
    # 	for k in range(0,len(allPopData)):
    # 		r=allPopData[k][0]
    # 		popCELE=allPopData[k][1]
    # 		#import the country boundaries, so we can see which country coordinates fall into
    # 		world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
    # 		countryshapes=world['geometry'].values
    # 		dfpoint=pd.DataFrame({'longs':popCELE['longs'].values,'lats':popCELE['lats'].values,'popCELE':popCELE['pop'].values})
    # 		crs={'init':'epsg:3857'}
    # 		geometry=[Point(xy) for xy in zip(dfpoint['longs'],dfpoint['lats'])]
    # 		geo_df=gpd.GeoDataFrame(dfpoint,crs=crs,geometry=geometry)
    # 		world['countrygeometry']=world['geometry']
    # 		#get groups of countries which match the coordinates for population losing electricity
    # 		pointInPolys = sjoin(geo_df, world, how='left')

    # 		grouped = pointInPolys.groupby('index_right')

    # 		eleLostArr=[]
    # 		fractionEleLostArr=[]
    # 		elePerCapita=[]
    # 		eleCountry=[]
    # 		countrycode=[]
    # 		geometries=[]

    # 		#for each country
    # 		netpop=0
    # 		for key, values in grouped:
    # 			code=values['iso_a3'].values[0]
    # 			# eleByCountry['countrycode'].str.match(code)
    # 			countryele=[]
    # 			for i in range(0,len(eleByCountry['countrycode'])):
    # 				if(eleByCountry['countrycode'].values[i]==code):
    # 					countryele=eleByCountry.iloc[i]
    # 					break
    # 			if(len(countryele)==0):
    # 				continue

    # 			countrypop=values['pop_est'].values[0]
    # 			popCELE=values['popCELE'].sum()
    # 			eleconsumption=countryele['consumption']
    # 			# print('eleconsumption')
    # 			# print(eleconsumption)
    # 			# print('popCELE')
    # 			# print(popCELE)
    # 			# print('countrypop')
    # 			# print(countrypop)
    # 			eleLost=eleconsumption*(popCELE/countrypop)
    # 			countrycode.append(code)
    # 			eleCountry.append(eleconsumption)
    # 			elePerCapita.append(eleconsumption/countrypop)
    # 			eleLostArr.append(eleLost)
    # 			fractionEleLostArr.append(eleLost/eleconsumption)
    # 			geometries.append(values['countrygeometry'].values[0])

    # 		worldElectricity=gpd.GeoDataFrame({'iso_a3':np.array(countrycode),'geometry':np.array(geometries),'fraction':np.array(fractionEleLostArr),'total':np.array(eleCountry),'totalLost':np.array(eleLostArr),'elePerCapita':np.array(elePerCapita)},crs=crs)

    # 		#set area of each country
    # 		worldElectricity['area']=worldElectricity['geometry'].to_crs({'init': 'epsg:3395'}).map(lambda p: p.area / 10**6)

    # 		#electricity per km^2
    # 		worldElectricity['eleDensity']=np.log(worldElectricity['total']/worldElectricity['area'])

    # 		fig, ax = plt.subplots(1, 1)
    # 		worldElectricity.plot(column='fraction', ax=ax, legend=True,vmin=0,vmax=1)

    # 		plt.title('Predicted Fraction Electricity By Country \n One in '+str(1/r)+' Year Storm')
    # 		plt.show()

    # 		totalEleLoss=np.sum(eleLostArr)
    # 		totalEle=np.sum(eleCountry)
    # 		print('totalElectricity (GW)')
    # 		print(totalEle)
    # 		print('totalEleLoss (GW)')
    # 		print(totalEleLoss)
    # 		print('fraction overall')
    # 		print(totalEleLoss/totalEle)

    # 		# fig, ax = plt.subplots(1, 1)
    # 		# worldElectricity.plot(column='eleDensity', ax=ax, legend=True,vmin=np.min(worldElectricity['eleDensity']),vmax=np.max(worldElectricity['eleDensity']))

    # 		# plt.title('average Electricity Consumption By Country per unit kilometer')
    # 		# plt.show()

    # 		# fig, ax = plt.subplots(1, 1)
    # 		# worldElectricity.plot(column='totalLost', ax=ax, legend=True,vmin=0,vmax=1)

    # 		# plt.title('Predicted total Electricity lost By Country \n One in '+str(1/r)+' Year Storm')
    # 		# plt.title('Total Electricity Consumption  By Country')
    # 		# plt.show()

    # 		pdata=rasterio.open(Params.popDensity15min)
    # 		ldata=rasterio.open(Params.landArea15min)
    # 		pArr=pdata.read(1)
    # 		lArr=ldata.read(1)
    # 		pArrZeroed = np.where(pArr<0, 0, pArr)
    # 		lArrZeroed = np.where(lArr<0, 0, lArr)
    # 		totPop=np.multiply(pArrZeroed,lArrZeroed)
    # 		longs=[]
    # 		lats=[]
    # 		highResEle=[]
    # 		for latindex in range(0,len(totPop)):
    # 			for longindex in range(0,len(totPop[0])):
    # 				lats.append(90-latindex*.25)
    # 				longs.append(longindex*.25-180)
    # 				highResEle.append(totPop[latindex,longindex])
    # 				#figure out if point is inside country
    # 				# grouped = pointInPolys.groupby('index_right')

    # 				# for key, values in grouped:
    # 				# 	code=values['iso_a3'].values[0]
    # 				# 	countryele=[]
    # 				# 	for i in range(0,len(worldElectricity['iso_a3'])):
    # 				# 		if(worldElectricity['iso_a3'].values[i]==code):
    # 				# 			highResEle.append(worldElectricity['elePerCapita'].values[i]*totPop[latindex,longindex])

    # 				# 			break
    # 				# 	if(len(countryele)==0):
    # 				# 		continue
    # 		elemin=np.min(highResEle)+.1
    # 		elemax=np.max(highResEle)
    # 		df=pd.DataFrame({'longs':np.array(longs),'lats':np.array(lats),'highResEle':np.array(highResEle)})

    # 		crs={'init':'epsg:3857'}
    # 		geometry=[Point(xy) for xy in zip(df['longs'],df['lats'])]
    # 		latdiff=.25
    # 		longdiff=.25
    # 		geo_df=Plotter.calcGrid(df,latdiff,longdiff)
    # 		# geo_df=gpd.GeoDataFrame(df,crs=crs,geometry=geometry)

    # 		#loop through the 30 minute gridded population density data, and multiply the population in the country by the average electricity per person. Assuming everyone in a country uses about the same electricity, each data point gives the total electricity in that 30 minute square section.
    # 		ax=geo_df.plot(column='highResEle',legend=True,legend_kwds={'label': 'Coefficient on Reference Field Level','orientation': "horizontal"}, cmap='viridis',norm=colors.LogNorm(vmin=elemin, vmax=elemax))

    # 		pp=gplt.polyplot(world,ax=ax,zorder=1)
    # 		# print(df)
    # 		self.formatticklabels(elemin,elemax,pp)

    # 		plt.title('Geoelectric Field Multiplier\n Magnetic Latitude\n1 in '+str(1/r)+' Year Storm')

    # 		print('s10')
    # 		plt.show()

    def calcElectricityAffected(self, rateperyears, CELEpopAtRate):
        # import electricity by nation data

        rawimport = pd.read_csv(
            "Data/SmallData/ElectricityByNation/EleByCountry.csv", skiprows=4, header=0
        )
        eleByCountry = pd.DataFrame(
            {
                "countrycode": rawimport.iloc[1::3, 0].values,
                "country": rawimport.iloc[0::3, 1].values,
                "consumption": rawimport.iloc[2::3, 2].values,
                "generation": rawimport.iloc[1::3, 2].values,
            }
        )
        eleByCountry["countrycode"] = eleByCountry["countrycode"].str.slice(
            start=14, stop=17
        )
        codearr = rawimport.iloc[1::3, 0].str.split(pat="-").values

        codes = [""] * len(codearr)
        for i in range(0, len(codearr)):
            ar = codearr[i]
            if len(ar) < 2:
                codes[i] == ""
            else:
                codes[i] = ar[2]

        eleByCountry["countrycode"] = np.array(codes)

        eleByCountry["consumption"] = (
            eleByCountry["consumption"].replace("--", "0").astype(float)
        )
        eleByCountry["generation"] = (
            eleByCountry["generation"].replace("--", "0").astype(float)
        )

        # import the country boundaries, so we can see which country coordinates fall into
        world = geopandas.read_file(geopandas.datasets.get_path("naturalearth_lowres"))
        countryshapes = world["geometry"].values

        world["countrygeometry"] = world["geometry"]

        for r in rateperyears:
            CELEpop = CELEpopAtRate[r]
            CELEpop.to_pickle("CELEpop_%s.pkl" % str(r))
            del CELEpop["index_right"]
            # get groups of countries which match the coordinates for population losing electricity
            pointInPolys = sjoin(CELEpop, world)

            eleLostArr = []
            fractionEleLostArr = []
            elePerCapita = []
            eleCountry = []
            countrycode = []
            geometries = []

            # for each country
            netpop = 0
            for countryIndex in set(pointInPolys.index_right.values):
                row = world[world.index == countryIndex]

                code = row["iso_a3"].values[0]
                # eleByCountry['countrycode'].str.match(code)
                countryele = []
                for i in range(0, len(eleByCountry["countrycode"])):
                    if eleByCountry["countrycode"].values[i] == code:
                        countryele = eleByCountry.iloc[i]
                        break
                if len(countryele) == 0:
                    continue

                countrypop = row["pop_est"].values[0]
                countryCELEpop = pointInPolys[
                    pointInPolys["index_right"].values == countryIndex
                ]["Z"].sum()
                eleconsumption = countryele["consumption"]
                # print('eleconsumption')
                # print(eleconsumption)
                # print('popCELE')
                # print(countryCELEpop)
                # print('countrypop')
                # print(countrypop)
                eleLost = eleconsumption * (countryCELEpop / countrypop)
                countrycode.append(code)
                eleCountry.append(eleconsumption)
                elePerCapita.append(eleconsumption / countrypop)
                eleLostArr.append(eleLost)
                fractionEleLostArr.append(eleLost / eleconsumption)
                geometries.append(row["countrygeometry"].values[0])

            crs = {"init": "epsg:3857"}
            worldElectricity = gpd.GeoDataFrame(
                {
                    "iso_a3": countrycode,
                    "geometry": geometries,
                    "fraction": fractionEleLostArr,
                    "total": eleCountry,
                    "totalLost": eleLostArr,
                    "elePerCapita": elePerCapita,
                },
                crs=crs,
            )

            # set area of each country
            worldElectricity["area"] = (
                worldElectricity["geometry"]
                .to_crs({"init": "epsg:3395"})
                .map(lambda p: p.area / 10**6)
            )

            # electricity per km^2
            worldElectricity["eleDensity"] = np.log(
                worldElectricity["total"] / worldElectricity["area"]
            )

            totalEleLoss = np.sum(eleLostArr)
            totalEle = np.sum(eleCountry)
            print("rateperyear")
            print(r)
            print("totalElectricity (GW)")
            print(totalEle)
            print("totalEleLoss (GW)")
            print(totalEleLoss)
            print("fraction overall")
            print(totalEleLoss / totalEle)

            no_power_out = self.combinedDF[self.combinedDF[("powerOut" + str(r))] == 0]
            no_power_out = gpd.sjoin(self.rasterDF, no_power_out)
            no_power_out = no_power_out.drop(columns=["index_right"])
            no_power_out = no_power_out.sjoin(world)
            no_power_out = no_power_out[["iso_a3", "powerOut" + str(r)]]
            no_power_out = no_power_out.drop_duplicates()
            no_power_out = no_power_out[
                ~no_power_out["iso_a3"].isin(worldElectricity["iso_a3"])
            ]
            no_power_out["fraction"] = no_power_out["powerOut" + str(r)]

            worldElectricity = worldElectricity.append(no_power_out)[
                ["iso_a3", "fraction"]
            ]
            worldElectricity = world.merge(worldElectricity, on="iso_a3")
            worldElectricity.to_pickle("Data/SmallData/worldElectricity.pkl")

            fig, ax = plt.subplots()

            pp = gplt.polyplot(world, ax=ax, zorder=1)
            gplt.choropleth(
                worldElectricity,
                hue=worldElectricity["fraction"],
                cmap="viridis",
                ax=ax,
                legend=True,
                edgecolor="None",
            )

            plt.title(
                "Predicted Fraction Electricity Lost By Country \n One in "
                + str(1 / r)
                + " Year Storm"
            )
            plt.savefig(Params.figuresDir + str(r) + "fraction_electr_lost.png")
            plt.show()

    def formatticklabels(self, minval, maxval, pp):
        print("formatting")
        colourbar = pp.get_figure().get_axes()[1]
        ticks_loc = colourbar.get_xticks().tolist()
        ticks_loc.insert(0, minval)
        ticks_loc.append(maxval)
        colourbar.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))

        labels = []
        for i in range(0, len(ticks_loc)):
            x = ticks_loc[i]
            exponent = int(np.log10(x))
            coeff = x / 10**exponent
            if i == 0 or i == len(ticks_loc):
                labels.append(
                    r"${:2.1f} \times 10^{{ {:2d} }}$".format(coeff, exponent)
                )
            else:
                labels.append(
                    r"${:2.1f} \times 10^{{ {:2d} }}$".format(coeff, exponent)
                )
        colourbar.set_xticklabels(labels, rotation=33)

    # #high resolution population estimates (15 minute)
    # def calcPopulationAffected(self):
    # 	print('calculating population and regions with CELE')
    # 	pdata=rasterio.open(Params.popDensity15min)
    # 	ldata=rasterio.open(Params.landArea15min)
    # 	#latitude resolution in degrees divided by 30 minute resolution of population density data gives us the number of points below and above to sum the population for this transformer location.
    # 	populationLatN=np.int(np.floor(Params.latituderes/.25))
    # 	populationLongN=np.int(np.floor(Params.longituderes/.25))
    # 	pArr=pdata.read(1)
    # 	lArr=ldata.read(1)
    # 	pArrZeroed = np.where(pArr<0, 0, pArr)
    # 	lArrZeroed = np.where(lArr<0, 0, lArr)
    # 	totPop=np.multiply(pArrZeroed,lArrZeroed)
    # 	overheatdata=np.load(Params.overheatDataDir,allow_pickle=True)
    # 	allPopData=[]
    # 	for i in range(0,int(np.floor(len(overheatdata))/2)):
    # 		totalPopCELE=0
    # 		rateindex=i*2
    # 		dataindex=i*2+1
    # 		r=overheatdata[rateindex]
    # 		print('rateperyear'+str(r))
    # 		df=overheatdata[dataindex]
    # 		populationaffected=[]
    # 		for j in range(0,len(df)):
    # 			isCELE=df['isCELEs'].values[j]
    # 			latitude=df['lats'].values[j]
    # 			longitude=df['longs'].values[j]
    # 			if(not isCELE):
    # 				indexlat,indexlong=pdata.index(longitude,latitude)

    # 				population=0
    # 				#sum over the populations near the grid point of interest
    # 				for k in range(0,populationLatN):
    # 					newLatIndex=indexlat-populationLatN//2+k
    # 					for l in range(0,populationLongN):
    # 						newLongIndex=indexlong-populationLongN//2+l
    # 						population=population+totPop[newLatIndex,newLongIndex]
    # 						totalPopCELE=totalPopCELE+totPop[newLatIndex,newLongIndex]
    # 				populationaffected.append(population)
    # 			else:
    # 				populationaffected.append(0)

    # 		#Catastrophic Electricity Loss assumed if spares<overheated transformers
    # 		regionsCELE= np.where(np.array(populationaffected)>0, 1, 0) #population affected is always the full population in the region where spares<overheated transformers. Nonzero=>population in grid cell.
    # 		df=pd.DataFrame({'longs':df['longs'].values,'lats':df['lats'].values,'pop':np.array(populationaffected),'regionsCELE':np.array(regionsCELE)})
    # 		crs={'init':'epsg:3857'}
    # 		# nlatcells=(xmax-xmin)/latdiff
    # 		# nlongcells=(ymax-ymin)/longdiff

    # 		#turn the point values into a grid cell with the proper width and height
    # 		latdiff=Params.latituderes#
    # 		longdiff=Params.longituderes#
    # 		geo_df=Plotter.calcGrid(df,latdiff,longdiff)

    # 		ax=geo_df.plot(column='regionsCELE')
    # 		world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
    # 		pp=gplt.polyplot(world,ax=ax,zorder=1)

    # 		plt.title('Predicted Regions HVT Overheated Exceed Spares \n'+str(1/r)+' Year Storm')
    # 		# gplt.kdeplot(geo_df,column='E',)
    # 		print(Params.figuresDir+'/RegionsCELE/regionsAffected'+str(r)+ 'peryear.png')
    # 		plt.savefig(Params.figuresDir+'/RegionsCELE/regionsAffected'+str(r)+ 'peryear.png')
    # 		plt.show()
    # 		print('world population')
    # 		print(np.sum(totPop))
    # 		print('world population CELE, 1 in '+str(1/r)+' year storm')
    # 		print(totalPopCELE)
    # 		print('fraction world population CELE, 1 in '+str(1/r)+' year storm')
    # 		print(totalPopCELE/np.sum(totPop))
    # 		allPopData.append([r,df])
    # 	np.save(Params.popCELEdir,allPopData)

    def createNetwork(self):
        # add together the high voltage power grid of ... THE ENTIRE EARTH
        loaddirectory = Params.planetNetworkDir
        print(loaddirectory)
        allfiles = glob.glob(loaddirectory + "*")
        print(allfiles)
        networks = []
        alllines = []
        allvoltages = []
        for i in range(0, len(allfiles)):
            fn = allfiles[i]
            if fn == "." or fn == "..":
                continue
            network = Network()
            print("Processing " + fn)
            # if(i==1):
            # continue
            # leave country blank, to indicate the whole continent
            [voltages, lines] = network.importNetwork(fn.split("/")[-1], "")
            allvoltages = np.append(allvoltages, np.array(voltages))
            alllines = np.append(alllines, np.array(lines))
            networks.append(network)
        self.networks = networks
        np.save(
            Params.networkAnalysisDir + "voltageData", allvoltages, allow_pickle=True
        )
        np.save(Params.networkAnalysisDir + "linesData", alllines, allow_pickle=True)

    def createRegionNetwork(self, continent, country):
        # add together the high voltage power grid of ... THE ENTIRE EARTH
        loaddirectory = Params.countryNetworkDir + continent + "/" + str(country)
        print(loaddirectory)
        alllines = []
        allvoltages = []

        network = Network()
        print("Processing " + loaddirectory)

        [voltages, lines] = network.importNetwork(continent, country)
        allvoltages = np.append(allvoltages, np.array(voltages))
        alllines = np.append(alllines, np.array(lines))
        if not os.path.isdir(Params.networkAnalysisDir + network.region):
            os.mkdir(Params.networkAnalysisDir + network.region)
        np.save(
            Params.networkAnalysisDir + network.region + "/voltageData",
            allvoltages,
            allow_pickle=True,
        )
        print("allvoltages saved")
        np.save(
            Params.networkAnalysisDir + "/" + network.region + "/linesData",
            alllines,
            allow_pickle=True,
        )
        print("linedata saved")
        self.networks = [network]

    def splitIntoStationRegions(self, network, rateperyears):
        loadNodes = True

        if loadNodes:
            network.nodesDict = np.load(
                Params.networkAnalysisDir + network.region + "/analyzedNodes.npy",
                allow_pickle=True,
            ).item()
        # use the country boundary polygons

        # find the population for each substation node, by assuming all the nodes serve the areas closest to them
        indices = []
        lats = []
        longs = []
        failureProbs = {}
        for r in rateperyears:
            failureProbs[r] = []

        for n in network.nodesDict.items():
            node = n[1]
            if node["nodeType"] == "substation":
                if node["voltageClass"] > 320:
                    indices.append(node["index"])
                    lats.append(float(node["lat"]))
                    longs.append(float(node["long"]))
                    for r in rateperyears:
                        if "probFail" in node.keys():
                            failureProbs[r].append(node["probFail"][r])
                        else:
                            failureProbs[r].append(0)
        df = pd.DataFrame({"longs": np.array(longs), "lats": np.array(lats)})
        for r in rateperyears:
            df[str(r)] = np.array(failureProbs[r])
        crs = {"init": "epsg:3857"}

        geometry = [Point(xy) for xy in zip(df["longs"], df["lats"])]
        geo_df = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)
        # ax=geo_df.plot(column='over',legend=True,legend_kwds={'label': 'Peak Temperature (C)','orientation': "horizontal"}, cmap='viridis',vmin=np.min(alltempsTie1),vmax=np.max(alltempsTie1))

        allStationRegions = network.calcVoronoi(geo_df, network.region, rateperyears)
        print(allStationRegions)

        if not os.path.isdir(Params.networkAnalysisDir + network.region):
            os.mkdir(Params.networkAnalysisDir + network.region)

        allStationRegions.to_pickle(
            Params.networkAnalysisDir + network.region + "/allStationRegions.pkl"
        )
        # allStationRegions=pd.read_pickle(Params.networkAnalysisDir+network.region+"/allStationRegions.pkl")

        for r in rateperyears:
            Plotter.plotCombinedVoronoi(allStationRegions, r, False)
        for r in rateperyears:
            Plotter.plotCombinedVoronoi(allStationRegions, r, True)

    def plotNetwork(self):
        if len(self.networks) > 1 or len(self.networks) == 0:
            # load up the high voltage lines and put together for plotting
            allvoltages = np.load(
                Params.networkAnalysisDir + "voltageData.npy", allow_pickle=True
            )
            allines = np.load(
                Params.networkAnalysisDir + "linesData.npy", allow_pickle=True
            )
        else:
            region = self.networks[0].region
            if not os.path.isdir(Params.networkAnalysisDir + region):
                os.mkdir(Params.networkAnalysisDir + region)
            allvoltages = np.load(
                Params.networkAnalysisDir + region + "/voltageData.npy",
                allow_pickle=True,
            )
            allines = np.load(
                Params.networkAnalysisDir + region + "/linesData.npy", allow_pickle=True
            )

        indices = range(0, len(allines))
        f, ax = plt.subplots()
        maxi = 100000
        rawLinesbyVoltageindices = [x for _, x in sorted(zip(allvoltages, indices))]
        rawLinesbyVoltage = []
        for i in rawLinesbyVoltageindices:
            rawLinesbyVoltage.append(allines[i])
        sortedVoltages = sorted(allvoltages)
        alluniquevoltages = list(set(sortedVoltages))
        prevVoltage = -1
        linesByVoltage = []
        uniqueIndex = 0
        for i in range(0, len(rawLinesbyVoltage)):
            if sortedVoltages[i] == prevVoltage:
                if len(rawLinesbyVoltage[i]) < 2:
                    linesByVoltage[uniqueIndex - 1] = np.array(
                        list(linesByVoltage[uniqueIndex - 1]) + rawLinesbyVoltage[i]
                    )
                else:
                    linesByVoltage[uniqueIndex - 1] = np.append(
                        linesByVoltage[uniqueIndex - 1], rawLinesbyVoltage[i]
                    )
                continue
            if len(rawLinesbyVoltage[i]) < 2:
                linesByVoltage.append(np.array([]))
                linesByVoltage[uniqueIndex - 1] = rawLinesbyVoltage[i]
            else:
                linesByVoltage.append(np.array(rawLinesbyVoltage[i]))

            uniqueIndex = uniqueIndex + 1
            prevVoltage = sortedVoltages[i]
        assert len(linesByVoltage) == len(alluniquevoltages)

        n = len(alluniquevoltages)
        colornames = plt.cm.coolwarm(np.linspace(0.1, 0.9, n))
        mycolors = [
            ListedColormap(colornames[i]) for i in range(0, len(alluniquevoltages))
        ]
        for i in range(0, len(alluniquevoltages)):
            v = alluniquevoltages[i]
            if len(linesByVoltage[i]) < 2:
                maxi = i
                continue

            network = gpd.GeoDataFrame(
                {"voltage": np.array(v), "geometry": np.array(linesByVoltage[i])}
            )
            gdf = gpd.GeoDataFrame(geometry=linesByVoltage[i])

            gdf.plot(
                cmap=mycolors[i],
                ax=ax,
                legend=True,
                label=str(alluniquevoltages[i]) + " kV network",
            )
        ax.legend()
        leg = ax.get_legend()
        for i in range(0, len(alluniquevoltages)):
            if i == maxi:
                break
            leg.legendHandles[i].set_color(colornames[len(colornames) - i - 1])

        plt.show()

    def calcTransformerFailures(self, network, durationratios, durations, rateperyears):
        print("5")
        nodesarr = network.sortednodes
        for r in rateperyears:
            savedata = np.load(
                Params.networkAnalysisDir
                + network.region
                + "/gics"
                + str(r)
                + "perYear.npy",
                allow_pickle=True,
            )

            # savedata=np.load('gics'+str(r)+'perYear.npy',allow_pickle=True)
            [gics, lineLengths] = savedata
            # print(gics)
            transformerDamageProbs = []
            nnodes = 0
            for i in range(0, len(gics)):
                gic = gics[i]
                # if(gic>10):
                # 	print('didit')
                # 	quit()
                # if(i==800):
                # 	quit()
                # continue
                # print('i')
                # print(i)
                # print('len(nodesarr)')
                # print(len(nodesarr))
                # print('len(gic)')
                # print(len(gics))
                # print(nodesarr[i])
                nodetype = nodesarr[i]["nodeType"]  # node category
                voltage = nodesarr[i][
                    "voltageClass"
                ]  # max voltage at station, or voltage of winding
                maxprobabilityDamaged = 0
                if nodetype == "substation":
                    nnodes = nnodes + 1
                    for k in range(0, len(durationratios)):
                        durationratio = durationratios[k]
                        duration = durations[k]

                        # adjust to make average gic of given duration longer, and also increase GIC to reflect amount going to transformer rather than node.
                        adjustedGIC = gic * durationratio / 1.3

                        # print('duration')
                        # print(duration)
                        # print('durationratio')
                        # print(durationratio)
                        # print('gic')
                        # print(gic)
                        # if(adjustedGIC>10):
                        # 	print('wow')
                        # 	quit()
                        # quit()
                        prob = self.calcProbabilityDamageFromGIC(
                            adjustedGIC, duration, float(voltage)
                        )
                        if prob > maxprobabilityDamaged:
                            maxprobabilityDamaged = prob
                    transformerDamageProbs.append(maxprobabilityDamaged)
                    # save the probability of damage in the nodesDict for later analysis
                    if not (
                        "probFail"
                        in network.nodesDict[nodesarr[i]["substationName"]].keys()
                    ):
                        network.nodesDict[nodesarr[i]["substationName"]][
                            "probFail"
                        ] = {}
                    network.nodesDict[nodesarr[i]["substationName"]]["probFail"][
                        r
                    ] = maxprobabilityDamaged
            print("total nodes")
            print(nnodes)
            print("np.sum(maxprobabilityDamaged)")
            print(np.sum(transformerDamageProbs))
            print("fractiontotaldamaged")
            print(np.sum(transformerDamageProbs) / nnodes)
            print("r")
            print(r)
            # plt.plot(sorted(transformerDamageProbs))
            # plt.show()
        if not os.path.isdir(Params.networkAnalysisDir + network.region):
            os.mkdir(Params.networkAnalysisDir + network.region)
        np.save(
            Params.networkAnalysisDir + network.region + "/analyzedNodes.npy",
            network.nodesDict,
            allow_pickle=True,
        )

    def calcGICs(networks):
        gicsEachNetwork = []

        for n in networks:
            gicsEachNetwork.append(n.calcGICs())
            # by duration, calc temperature, and take maximum probability transformer will erheat at any given duration

    def combineRegions(self, regions, rateperyears):
        prevDF = []
        for region in regions:
            if not os.path.isfile(
                Params.networkAnalysisDir + region + "/allStationRegions.pkl"
            ):
                continue
            regionDF = pd.read_pickle(
                Params.networkAnalysisDir + region + "/allStationRegions.pkl"
            )
            if len(prevDF) > 0:
                combinedDF = pd.concat([prevDF, regionDF], ignore_index=True)
            else:
                combinedDF = regionDF
            prevDF = combinedDF
        for r in rateperyears:
            Plotter.plotCombinedVoronoi(combinedDF, r, True)
        self.combinedDF = combinedDF

    def specifyCombinedRegion(self, continent, country, rateperyears):
        if country:
            regionDF = pd.read_pickle(
                Params.networkAnalysisDir + country + "/allStationRegions.pkl"
            )
        else:
            regionDF = pd.read_pickle(
                Params.networkAnalysisDir + continent + "/allStationRegions.pkl"
            )
        for r in rateperyears:
            Plotter.plotCombinedVoronoi(regionDF, r, True)
        self.combinedDF = regionDF

    # https://gis.stackexchange.com/questions/384581/raster-to-geopandas
    def rasterToDF(self):
        ldata = rasterio.open(Params.landArea15min)
        lArr = ldata.read(1)
        lArrZeroed = np.where(lArr < 0, 0, lArr)
        with rio.Env():
            with rio.open(Params.popDensity15min) as src:
                crs = src.crs

                # create 1D coordinate arrays (coordinates of the pixel center)
                xmin, ymax = np.around(src.xy(0.00, 0.00), 9)  # src.xy(0, 0)
                xmax, ymin = np.around(
                    src.xy(src.height - 1, src.width - 1), 9
                )  # src.xy(src.width-1, src.height-1)
                x = np.linspace(xmin, xmax, src.width)
                y = np.linspace(
                    ymax, ymin, src.height
                )  # max -> min so coords are top -> bottom

                # create 2D arrays
                xs, ys = np.meshgrid(x, y)
                pArr = src.read(1)
                pArrZeroed = np.where(pArr < 0, 0, pArr)
                zs = np.multiply(pArrZeroed, lArrZeroed)

        data = {
            "X": pd.Series(xs.ravel()),
            "Y": pd.Series(ys.ravel()),
            "Z": pd.Series(zs.ravel()),
        }

        df = pd.DataFrame(data=data)
        geometry = gpd.points_from_xy(df.X, df.Y)
        gdf = gpd.GeoDataFrame(df, crs={"init": "epsg:3857"}, geometry=geometry)
        self.rasterDF = gdf  # assign to the object for repeat use
        return gdf

    # high resolution population estimates (15 minute)
    def calcPopulationPowerOut(self, rateperyears):
        populationpoints = self.rasterToDF()  # .to_crs('EPSG:3857')
        print("")
        print("calculating population and regions with CELE")

        print("total population")
        print(populationpoints["Z"].sum())
        outageRegionsAllRates = {}
        for rate in rateperyears:
            # determine power outage locations
            self.combinedDF["powerOut" + str(rate)] = 0
            for index, row in self.combinedDF.iterrows():
                if row[str(rate)] > 0.33:
                    self.combinedDF["powerOut" + str(rate)].iloc[index] = 1
                else:
                    self.combinedDF["powerOut" + str(rate)].iloc[index] = 0
            popsum = 0
            outageRegions = self.combinedDF[
                self.combinedDF[("powerOut" + str(rate))] == 1
            ]
            # outageRegions.to_pickle("outageRegions.pkl")
            # populationpoints.to_pickle("populationpoints.pkl")

            # this should be done with higher resolution population data at some point, to make sure it is accurate.It will tend to underestimate because it only uses grid points inside the polygon. The problem is the higher population at 2.5 arcsecond will take a very long time to calculate.

            pointsInPoly = gpd.sjoin(populationpoints, outageRegions)
            # print('Population affected at rate '+str(rate) +': ')
            populationAffectedThisRate = pointsInPoly["Z"].sum()
            # print(populationAffectedThisRate)

            print("world population CELE, 1 in " + str(1 / rate) + " year storm")
            print(populationAffectedThisRate)
            print(
                "fraction world population CELE, 1 in " + str(1 / rate) + " year storm"
            )
            print(populationAffectedThisRate / populationpoints["Z"].sum())
            outageRegionsAllRates[rate] = pointsInPoly
        return outageRegionsAllRates

    # compare the GIC calculated at the country and the continent
    # this only works for the simplistic model
    def compareGICresults(self, continent, country, rateperyears):
        an = np.load(
            "Data/SmallData/Networks/" + continent + "/analyzedNodes.npy",
            allow_pickle=True,
        )

        # get the nodes which were used for gics
        sortednodes = [
            v for k, v in sorted(an.item().items(), key=lambda item: item[1]["index"])
        ]
        consubnames = [
            x["substationName"] for x in sortednodes if ("savenode" in x.keys())
        ]
        conlats = [x["lat"] for x in sortednodes if ("savenode" in x.keys())]
        conlongs = [x["long"] for x in sortednodes if ("savenode" in x.keys())]
        contypes = [x["nodeType"] for x in sortednodes if ("savenode" in x.keys())]
        convolts = [x["voltageClass"] for x in sortednodes if ("savenode" in x.keys())]

        an = np.load(
            "Data/SmallData/Networks/" + country + "/analyzedNodes.npy",
            allow_pickle=True,
        )

        # get the nodes which were used for gics
        sortednodes = [
            v for k, v in sorted(an.item().items(), key=lambda item: item[1]["index"])
        ]
        cousubnames = [
            x["substationName"] for x in sortednodes if ("savenode" in x.keys())
        ]
        coulats = [x["lat"] for x in sortednodes if ("savenode" in x.keys())]
        coulongs = [x["long"] for x in sortednodes if ("savenode" in x.keys())]
        coutypes = [x["nodeType"] for x in sortednodes if ("savenode" in x.keys())]
        couvolts = [x["voltageClass"] for x in sortednodes if ("savenode" in x.keys())]
        print(cousubnames)
        print(coulats)
        # commented block below should also work, but may not have been updated at same time as analyzedNodes and gics[rate]perYear.npy

        # npf = np.float32
        # import the nodes from the continent
        # connetwork=open('Data/BigData/Networks/'+continent+'Nodes.txt','r')
        # netdata = connetwork.readlines()
        # nnodes = len(netdata)
        # consitenum = np.zeros(nnodes, dtype=np.int32)
        # congeolat, congeolon = np.zeros(nnodes, dtype=npf), np.zeros(nnodes, dtype=npf)
        # consitename=[]
        # for i in range(nnodes):
        # 	data = netdata[i].split("\t")
        # 	congeolat[i] = float(data[4])
        # 	congeolon[i] = float(data[5])
        # 	consitename.append(data[1])
        # 	consitenum[i] = int(data[0])

        # #import the nodes from the country
        # counetwork=open('Data/BigData/Networks/'+country+'Nodes.txt','r')
        # netdata = counetwork.readlines()
        # nnodes = len(netdata)
        # cousitenum = np.zeros(nnodes, dtype=np.int32)
        # cougeolat, cougeolon = np.zeros(nnodes, dtype=npf), np.zeros(nnodes, dtype=npf)
        # cousitename=[]
        # for i in range(nnodes):
        # 	data = netdata[i].split("\t")
        # 	cougeolat[i] = float(data[4])
        # 	cougeolon[i] = float(data[5])
        # 	cousitename.append(data[1])
        # 	cousitenum[i] = int(data[0])

        for rate in rateperyears:
            # import GICs from continent. GIC index is the same as the node index
            conGICs = np.load(
                "Data/SmallData/Networks/"
                + continent
                + "/gics"
                + str(rate)
                + "perYear.npy",
                allow_pickle=True,
            )[0]

            # import GICs from country. GIC index is the same as the node index
            couGICs = np.load(
                "Data/SmallData/Networks/"
                + country
                + "/gics"
                + str(rate)
                + "perYear.npy",
                allow_pickle=True,
            )[0]

            matchedConGICs = []
            matchedCouGICs = []

            matchedconlats = []
            matchedconlongs = []
            matchedcoulats = []
            matchedcoulongs = []
            matchedconvolts = []
            matchedcouvolts = []
            for sn in cousubnames:
                conindex = np.where(sn == np.array(consubnames))
                if (
                    len(conindex[0]) == 0
                ):  # not all continent substations match country substation
                    continue
                if len(conindex[0]) > 1:
                    print("ERROR: MULTIPLE CONTINENT SUBSTATION MATCH NAME")
                    break

                couindex = np.where(sn == np.array(cousubnames))
                # print(couindex)
                # print(conindex)
                if len(couindex[0]) == 0:
                    print("ERROR: NO COUNTRY SUBSTATIONS MATCH NAME")
                    break
                if len(couindex[0]) > 1:
                    print("ERROR: MULTIPLE COUNTRY SUBSTATIONS MATCH NAME")
                    break

                matchedConGICs.append(conGICs[conindex[0][0]])
                matchedCouGICs.append(couGICs[couindex[0][0]])

                matchedcoulats.append(coulats[couindex[0][0]])
                matchedcoulongs.append(coulongs[couindex[0][0]])
                matchedconlats.append(conlats[conindex[0][0]])
                matchedconlongs.append(conlongs[conindex[0][0]])

                matchedconvolts.append(convolts[conindex[0][0]])
                matchedcouvolts.append(couvolts[couindex[0][0]])
                # print(couGICs[couindex[0][0]])
                # print(conGICs[conindex[0][0]])
                # print('')
            plt.figure()
            plt.scatter(matchedconlats, matchedcoulats)
            plt.show()
            plt.figure()
            plt.scatter(matchedconlongs, matchedcoulongs)
            plt.show()
            plt.figure()
            plt.scatter(matchedconvolts, matchedcouvolts)
            plt.show()
            matchedconlats = []
            matchedconlats = []
            matchedcoulats = []
            matchedcoulongs = []
            assert len(matchedCouGICs) == len(matchedConGICs)

            plt.figure()
            plt.scatter(matchedCouGICs, matchedConGICs)
            plt.title("rate: " + str(rate))
            plt.show()

            #

            # import the nodes from the country
