import numpy as np
#calculate thermal rise from duration and GIC level
import Params
class PowerGrid:

	def __init__(self):
		Params.importIfNotAlready()

	def loadEfieldMap(self):


	def getTemp(deltat,E,tau,alpha):
		E*alpha*(1-np.exp(-deltat/tau))
		return E*alpha*(1-np.exp(-deltat/tau))


	def importODS(self,fn):
		print("doing things")
		print(Params.beta)