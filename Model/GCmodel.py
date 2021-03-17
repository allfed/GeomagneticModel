import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import Params
import os

u0 = 1.25664e-6 #m kg s^-2 A^-2 (permeability of free space)

class GCmodel:
	def __init__(self):
		Params.importIfNotAlready()
		globalcondfn = Path(Params.globalcondloc)
		self.zgrid = []
		self.deltas=[]
		self.depths=[]
		self.downloaded=globalcondfn.is_file()
		self.w = 2*np.pi*(1.0/120) #s^-1 
		self.f = 1/120
		if not self.downloaded:
			currentdir=os.getcwd()

			print("Warning: Missing global conductivity model. Expected location in filesystem: "+currentdir+'/'+Params.globalcondloc)
			print("It should be the matlab version, with filename renamed from MODEL-VER-[version].mat to 'GCmodel.mat'. Dowload it from https://globalconductivity.ocean.ru/")
			print("Using pre-processed conductivity file instead.")

	def importGCmodel(self):

		#import from global conductivity model
		with h5py.File(Params.globalcondloc, 'r') as f:
			self.zgrid = f.get('ZGRID')[:]
			self.deltas = np.diff(np.transpose(self.zgrid)) #the difference between heights of layers
			self.depths = np.insert(deltas[0], 0, 100., axis=0) #the first layer is 100m thick


	# Calculate impedance recursively (see Chapter 1, Theory section of the NERC Application Guide)
	# we first calculate the bottom layer impedance, and repeat for each layer up.
	# Returns impedance at ground level
	def getImpedance(c):
		Z = [0]*len(self.depths) # Ohms (empty array of impedances)

		bottom=len(self.depths)-1

		# For each layer n starting at the bottom, calculate impedance. 
		# The final Z value at n=0 is the impedance at ground level.
		for n in range(bottom,-1,-1):
			k = np.sqrt(1j*self.w*u0*c[n]) #propagation constant
			
			if(n==bottom):
				Z[n]=1j*self.w*u0/k
			else:
				foo=k*Z[n+1]/(1j*self.w*u0)
				r=(1-foo)/(1+foo) #reflection coefficient
				Znumerator = 1-r*np.exp(-2*k*self.depths[n])
				Zdenominator = k*(1+r*np.exp(-2*k*self.depths[n]))

				Z[n] = 1j*self.w*u0*Znumerator/Zdenominator
		return Z[0]

	# here we loop through the map and determine apparent conductivity at each location, returning a 2d matrix
	# global conductivity model to determine local ground properties 
	# returns array of conductivities (depths specified with the ZGRID property from global conductivity model)
	def getImpedanceMap(latInterval,longInterval):
		maxcond = -1e19
		# maximp = -1e19
		mincond = 1e19
		# minimp = 1e19


		with h5py.File('../../Data/GeoelectricModel/MODEL-VER-5.mat', 'r') as f:
			row=0
			col=0
			lenlat=len(f['MODEL'][0,0,:])
			lenlong=len(f['MODEL'][0,:,0])

			# impArr = np.zeros((lenlong/longInterval,lenlat/latInterval))
			condArr = np.zeros((lenlong/longInterval, lenlat/latInterval))

			print('lenlong'+str(len(f['MODEL'][0,:,0])))
			print('lenlat'+str(len(f['MODEL'][0,0,:])))
			for longkey in range(0,lenlong-1,longInterval):
				print('longkey'+str(longkey-1))
				for latkey in range(0,lenlat,latInterval):

					c=f['MODEL'][:,longkey,latkey]
					compleximp=getImpedance(c)
					impedance = float(abs(compleximp))

					cond = u0*self.w/(impedance*impedance)

					if maxcond < cond:
						maxcond = cond

					if mincond > cond:
						mincond = cond

					# if maximp < impedance:
					# 	maximp = impedance

					# if minimp > impedance:
					# 	minimp = impedance

					# impArr[row,col] = impedance
					condArr[row,col] = cond

					row = row + 1

				row = 0
				col = col+1
			np.save('condreallybig',condArr)
			# np.save('imp',impArr)

			condForShow = np.flipud((np.log(condArr)-np.log(mincond))/(np.log(maxcond)-np.log(mincond)))
			# impForShow = np.flipud((np.log(impArr)-np.log(minimp))/(np.log(maximp)-np.log(minimp)))

			# print(impForShow[:,30])
			print(condForShow[:,30])


			plt.figure(1)
			plt.imshow(condForShow, cmap='hot', interpolation='nearest')
			plt.show()

			# plt.figure(2)
			# plt.imshow(impForShow, cmap='hot', interpolation='nearest')
			# plt.show()
		return 