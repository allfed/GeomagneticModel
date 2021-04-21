import fits
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import Params
import numpy as np
class ValidateModel:
	def __init__(self):

		#data comes from 100 year geoelectric fields from Love, et al
		# https://doi.org/10.1002/2017GL076042
		Params.importIfNotAlready()
		self.windowedEmaps=[]
		self.validationFields=[[-75.588,39.559,8.4333,7.4342],
		[-79.103,39.620,0.3435,0.2874],
		[-75.927,38.987,0.9627,0.7466],
		[-76.167,38.454,1.0162,1.1844],
		[-75.422,38.354,0.2585,0.2426],
		[-78.839,36.430,1.2189,0.9552],
		[-78.058,36.221,9.4080,7.7796],
		[-77.368,36.226,2.1478,1.6062],
		[-76.403,35.961,0.5847,0.7370],
		[-75.956,35.864,0.3643,0.3435],
		[-79.337,35.808,4.8083,3.7731],
		[-78.454,35.705,4.3451,3.6907],
		[-77.688,35.634,3.8815,2.9337],
		[-76.911,35.473,1.2574,1.0321],
		[-80.092,35.358,13.3659,11.1008],
		[-79.478,35.339,1.2897,0.8309],
		[-74.426,40.069,4.8865,3.9397],
		[-74.770,39.506,0.7780,0.7582],
		[-74.961,38.955,0.4623,0.3530],
		[-76.592,41.135,0.1655,0.0953],
		[-79.438,40.911,0.3881,0.2743],
		[-77.677,40.705,2.4462,1.4137],
		[-75.409,40.183,3.3923,3.1532],
		[-79.460,40.215,0.4031,0.2193],
		[-78.727,40.228,0.4078,0.3908],
		[-77.925,39.992,1.0990,1.1358],
		[-77.317,40.007,14.1253,12.8255],
		[-76.330,39.902,4.2121,4.3313],
		[-75.763,41.019,0.0765,0.0616],
		[-80.005,41.015,0.6823,0.8587],
		[-78.599,40.726,0.1655,0.1680],
		[-75.111,40.273,5.0176,4.5260],
		[-77.477,39.493,15.6494,12.1642],
		[-79.962,35.838,8.7096,6.9507],
		[-76.742,40.581,0.0641,0.0505],
		[-77.178,38.559,1.7559,1.2243],
		[-78.483,38.866,2.6546,2.5558],
		[-77.694,38.664,1.5381,1.2317],

		[-79.574,38.363,0.6823,0.5561],
		[-78.772,38.272,0.6584,0.6212],
		[-78.069,38.115,4.5551,3.3792],
		[-77.288,38.010,1.2545,0.9643],
		[-76.525,37.822,2.2387,1.9589],
		[-75.729,37.766,0.3999,0.4147],
		[-79.808,37.883,0.0467,0.0323],
		[-79.190,37.758,2.8641,2.3130],
		[-78.311,37.513,1.8923,1.2990],
		[-77.585,37.376,24.1824,19.9813],
		[-76.888,37.317,2.3577,1.9495],
		[-75.963,37.135,0.3706,0.4652],
		[-80.169,37.321,0.5035,0.3759],
		[-79.583,37.095,5.7147,4.5740],
		[-78.645,37.005,0.7780,0.8455],
		[-77.877,36.771,3.7110,2.3666],
		[-77.101,36.824,1.8556,1.5168],
		[-76.338,36.580,0.8053,0.9863],
		[-79.766,36.587,4.4004,4.6028],
		[-79.825,39.614,0.2648,0.2320],
		[-78.319,39.492,3.6307,2.7529],
		[-80.071,39.039,0.5854,0.5974],
		[-79.346,39.129,0.8433,0,5388]]

	def compareFields(self,earthmodel):
		loaddirectory=Params.globalEfieldData
		self.windowedEmaps=np.load(loaddirectory,allow_pickle=True)
		for i in range(0,int(np.floor(len(self.windowedEmaps))/2)):
			alldurations=[]
			durations=np.array([])
			rateindex=i*2
			dataindex=i*2+1
			r=self.windowedEmaps[rateindex]
			EmapAtDuration=self.windowedEmaps[dataindex]

			latitudes=earthmodel.GClatitudes
			longitudes=earthmodel.GClongitudes
			if(r==.01):
				ratios=[]
				ourEpredicteds=[]
				theirEpredicteds=[]
				for i in range(0,len(self.validationFields)):
					row=self.validationFields[i]
					longitude=row[0]
					latitude=row[1]
					theirEpredicted=row[2]

					[latid,longid]=earthmodel.getClosestCoordsIds(latitude,longitude)
					Emap60=EmapAtDuration[1]
					print('latid')
					print('longid')
					print(latid)
					print(longid)
					ourEpredicted=Emap60[latid,longid]
					ratio=ourEpredicted/theirEpredicted
					ratios.append(ratio)
					ourEpredicteds.append(ourEpredicted)
					theirEpredicteds.append(theirEpredicted)
				plt.figure()
				[xo,yo]=fits.binlognormaldist(ourEpredicteds,[],3,)
				[xt,yt]=fits.binlognormaldist(theirEpredicteds,[],3)
				plt.xscale("log")
				# plt.loglog()
				plt.title("Northeast US 100 year E field levels")
				plt.xlabel("Geoelectric Field (V/km)")
				plt.ylabel("PDF probability Site at Field Level")
				plt.plot(xo,yo, lw=1,label = " Our Predictions")
				plt.plot(xt,yt, lw=1,label = "Love,2017 Predictions")
				plt.legend()
				plt.show()

				plt.figure()
				[x,y]=fits.binlognormaldist(ratios,[],3)
				plt.xscale("log")
				plt.plot(x,y)
				plt.show()
				plt.figure()
				plt.loglog()
				plt.scatter(ourEpredicteds,theirEpredicteds)
				plt.title("Validation Against Love 100 year E fields")
				plt.xlabel("Geoelectric Field Global Model (V/km)")
				plt.ylabel("Geoelectric Field Love et. al. (V/km)")
				corr, _ = pearsonr(np.log(ourEpredicteds),np.log(theirEpredicteds))
				print(corr)

				plt.show()

	def calcGIC(self,E):
		# ntransformers=len(self.HVTnames)
		plt.title('Assumed GIC at E field of '+str(E)+'V/km')
		gic=np.array([9,9,34,37])*E
		fractions=[0.38502080443828,0.436061026352289,0.162829403606103,0.016088765603329]

		counts=np.floor(np.array(fractions)*100).astype(int)
		print('counts')
		print(counts)

		allcounts=np.array([gic[0]]*counts[0])
		allcounts=np.append(allcounts,[gic[1]]*counts[1])
		allcounts=np.append(allcounts,[gic[2]]*counts[2])
		allcounts=np.append(allcounts,[gic[3]]*counts[3])
		print('allcounts')
		print(allcounts)

		# plt.figure()
		# plt.plot(gic,fractions)
		plt.xlabel('Transformer Number')
		plt.ylabel('GIC (A)')
		# plt.figure()
		plt.bar(range(0,len(allcounts)),allcounts)
		plt.show()

				
		es=[]
		meangics=[]
		for tene in range(0,40):
			e=tene/10
			gic=np.array([9,9,34,37])*e
			meangic=np.sum(np.multiply(gic,np.array(fractions)))
			meangics.append(meangic)
			es.append(e)
		plt.plot(es,meangics)
		plt.xlabel("electric field (V/km)")
		plt.ylabel("Mean GIC over all transformers")
		plt.show()