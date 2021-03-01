import numpy as np
#calculate thermal rise from duration and GIC level
import Params
import h5py
import matplotlib.pyplot as plt

class PowerGrid:

	def __init__(self):
		Params.importIfNotAlready()
		self.windowedEmaps=[]
		self.HVTnames=Params.HVTnames
		self.alpha=Params.alpha
		self.beta=Params.beta
		self.pop25to40=Params.pop25to40
		self.pop0to25=Params.pop0to25
		self.pop40plus=Params.pop40plus
		self.tau=Params.tau

	def setWindowedEmaps(self,earthmodel):
		self.windowedEmaps=earthmodel.windowedEmaps

	def loadEfieldMaps(self):
		loaddirectory=Params.globalEfieldData
		self.windowedEmaps=np.load(loaddirectory)
		print('windowedEmaps')
		print(self.windowedEmaps)

	def calcOverheatMap(self):
		#for each rate per year calculated in earth model
		for i in range(0,int(np.floor(len(self.windowedEmaps))/2)):
			rateindex=i*2
			dataindex=i*2+1
			r=self.windowedEmaps[rateindex]
			EmapAtDuration=self.windowedEmaps[dataindex]
			windowedEmapsatrate = []

			longInterval=int(np.floor(Params.longituderes/0.25))
			latInterval=int(np.floor(Params.latituderes/0.25))
			with h5py.File(Params.globalcondloc, 'r') as f:
				latitudes = f.get('LATGRID')[0]
				longitudes = f.get('LONGRID')[:,0]

				lenlat=len(latitudes)
				lenlong=len(longitudes)
				for j in range(0,int(np.floor(len(EmapAtDuration))/2)):
				
					durationindex = j*2
					Emapindex = j*2+1

					duration = EmapAtDuration[durationindex]
					print(duration)
					Emap = EmapAtDuration[Emapindex]
					thisdurationoverheatmap = np.zeros((int(lenlong/longInterval), int(lenlat/latInterval)))
					overheatFractionMap = np.zeros((int(lenlong/longInterval), int(lenlat/latInterval)))
					col=0
					for longkey in range(0,lenlong-1,longInterval):
						row=0
						for latkey in range(0,lenlat-1,latInterval):
							E=Emap[row,col]
							fractionover=self.calcFractionOverTemp(E,duration)
							thisdurationoverheatmap = np.zeros((int(lenlong/longInterval), int(lenlat/latInterval)))
							overheatFractionMap[row][col]=max(overheatFractionMap[row][col],fractionover)
							row = row+1	
						col = col+1
					plt.figure(1)
					plt.imshow(np.flipud(overheatFractionMap[8:,:]), cmap='hot', interpolation='nearest')
					plt.title('Overheat fractions '+str(r)+' per year, '+str(duration)+'s')
					cbar = plt.colorbar()
					cbar.set_label('Fraction overheat (1=100%)')
					plt.savefig(Params.overheatMapsDir+'Overheat'+str(r)+'perYearWindow'+str(duration)+'s.png')
					plt.show()
			plt.figure(1)
			plt.imshow(np.flipud(overheatFractionMap[8:,:]), cmap='hot', interpolation='nearest')
			plt.title('Overheat fractions '+str(r)+' per year, max of all durations')
			cbar = plt.colorbar()
			cbar.set_label('Fraction overheat (1=100%)')
			plt.savefig(Params.overheatMapsDir+'Overheat'+str(r)+'perYearMaxfromallWindows.png')
			plt.show()
			# self.windowedEmaps = self.overheatmaps + [r,windowedEmapsatrate]

		# np.save(Params.globalEfieldData,self.windowedEmaps)


	# calculate fraction of population that exceeds its temperature limit for every transformer type
	def calcFractionOverTemp(self,E,duration):
		fraction=0
		for i in range(0,len(self.HVTnames)):
			temp=(E*self.alpha[i]+self.beta[i])*(1-np.exp(-duration/self.tau[i]))+90
			# temps = temps + [temp]
			if(temp>140):
				fraction = fraction + self.pop40plus[i]
			if(temp>160):
				fraction = fraction + self.pop25to40[i]
			if(temp>180):
				fraction = fraction + self.pop0to25[i]
		return fraction