# A TF (transfer function) site refers to a measurement of the MT transfer function, aka the impedance tensor Z
# to download the TF data: http://ds.iris.edu/spud/emtf
# see https://code.usgs.gov/ghsc/geomag/emtf/fcu for impedance tensor transfer function documentation
# this site is useful as well https://ds.iris.edu/ds/products/emtf/
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import xml.etree.ElementTree as ET 

import numpy as np
import Params

u0 = 1.25664e-6 #N A^-2 = (permeability of free space, roughly same value inside the earth)

class TFsite:

	# initialize a site by providing a list of frequencies w to determine the transfer function.
	# if a path to a TF site EDI document is provided, the init function 
	def __init__(self, TFsiteindex):
		Params.importIfNotAlready()
		self.siteIndex=TFsiteindex
		folder=Params.tfsitesdir

		self.name=Params.tfsitenames[TFsiteindex]

		self.tfformat=Params.tfformats[TFsiteindex]

		self.TFsitefn = folder+self.name+'.'+self.tfformat
		print(self.TFsitefn)
		self.wInEDI = [] #the EDI file has its own list of frequencies, which we use as approximations to ours
		self.rows = [] #the EDI file frequencies are listed on specific rows and columns in the data.
		self.rowlen=0
		self.cols = [] #the EDI file frequencies are listed on specific rows and columns in the data.
		self.collen=0
		self.f=[]
		self.ZXX = [] # V m^-1 T^-1 = V A N^-1 (impedance tensor component XX)
		self.ZXY = [] # V m^-1 T^-1 = V A N^-1 (impedance tensor component XY)
		self.ZYX = [] # V m^-1 T^-1 = V A N^-1 (impedance tensor component YX)
		self.ZYY = [] # V m^-1 T^-1 = V A N^-1 (impedance tensor component YY)
		self.rhoXX = [] # V m A^-1 = Ohm m(apparent resistivity tensor component XX)
		self.rhoXY = [] # V m A^-1 = Ohm m(apparent resistivity tensor component XY)
		self.rhoYX = [] # V m A^-1 = Ohm m(apparent resistivity tensor component YX)
		self.rhoYY = [] # V m A^-1 = Ohm m(apparent resistivity tensor component YY)
		self.apparentc=[]
		if self.TFsitefn:
			if(self.tfformat=='xml'):
				self.importXML()
				self.wProvided = [x*2.0*np.pi for x in self.f] #rad s^-1 (angular frequency), list of frequencies provided
			elif(self.tfformat=='edi'):

				self.importFrequencies(self.TFsitefn)
				self.calcrows()
				self.calccols()
				self.T=[1.0/x for x in self.f] # s (period)
				self.wProvided = [x*2.0*np.pi for x in self.f] #rad s^-1 (angular frequency), list of frequencies provided
				self.importAllZ(self.TFsitefn)
			else:
				print('Error: TFsite format not supported. Choose EDI or XML.')
				quit()
			self.calcApparentc()
			self.calcApparentr()

	#calculate an array of all the columns at the TF site
	def calccols(self):
		cols=[]
		for i in range(0,self.collen):
			for j in range(0,self.rowlen):
				cols.append(j)
		self.cols=cols[0:len(self.f)]

	#calculate an array of all the rows at the TF site
	def calcrows(self):
		rows=[]
		for i in range(0,self.collen):
			for j in range(0,self.rowlen):
				rows.append(i)
		self.rows=rows[0:len(self.f)]

	# See Love "Geoelectric Hazard Maps for the
	# Mid-Atlantic United States: 100 Year Extreme Valuesand the 1989 Magnetic Storm" for the math to calculate apparent conductivity
	def printApparentc(self): 		
		for i in range(0,len(self.wInEDI)):
			print('Apparent conductivity at '+self.TFsitefn+ \
			' (S/m): '+str(self.apparentc[i])+ \
			' at w (rad/s):'+str(self.wInEDI[i]))
	
	def importXML(self):
		file=ET.parse(self.TFsitefn)
		root = file.getroot()
		self.T=[]
		self.f=[]
		for child in root:
			
			if(child.tag=='Data'):
				for period in child:
					self.T.append(np.float(period.attrib['value']))
					self.f.append(1/np.float(period.attrib['value']))
					for Zvals in period:
						if(Zvals.attrib['type']=='complex'):
							for value in Zvals:
								if(value.attrib['name']=='Zxx'):
									ZXX=value.text
									linevals = np.array(ZXX.split()).astype(np.float)
									# real and complex components
									self.ZXX.append(linevals[0]+1j*linevals[1])
								if(value.attrib['name']=='Zxy'):
									ZXY=value.text
									linevals = np.array(ZXY.split()).astype(np.float)
									# real and complex components
									self.ZXY.append(linevals[0]+1j*linevals[1])

								if(value.attrib['name']=='Zyx'):
									ZYX=value.text
									linevals = np.array(ZYX.split()).astype(np.float)
									# real and complex components
									self.ZYX.append(linevals[0]+1j*linevals[1])

								if(value.attrib['name']=='Zyy'):
									ZYY=value.text
									linevals = np.array(ZYY.split()).astype(np.float)
									# real and complex components
									self.ZYY.append(linevals[0]+1j*linevals[1])
	#imports the frequencies defined for this MT site. This defines rows and columns for the rest of the document.
	def importFrequencies(self,tfsitefn):
		file = open(tfsitefn, 'r')
		count=0
		while True: 
			count += 1
		  
			# Get next line from file 
			line = file.readline()

			# if line is empty 
			# end of file is reached 
			if not line: 
				break

			# find the frequency line and skip two lines down to get to the data
			if(line.find('****FREQUENCIES****')>-1):
				line = file.readline() #skip '>FREQ //[number of frequencies]' line

				#now we're at the frequencies, so we record all of them until we hit an empty line
				frequencies=[]
				collen=0
				while True:
					line = file.readline()
					linevals = np.array(line.split())
					if(collen==0):
						self.rowlen=len(linevals)

					# #the last row is shorter than the others in this case
					# if(linevals.size>0 and linevals.size<rowlen):
						
						
					#if frequencies are correct, return from function
					if(linevals.size==0):
						self.collen=collen
						self.f=frequencies
						return
					collen=collen+1
					freqsrow=list(linevals.astype(np.float))

					frequencies = frequencies + freqsrow
					if not line: 
						print('Error: Failed to import TFsite frequencies')
						quit()
		print('Error: Failed to find TFsite frequencies')
		quit()
		return

	#imports Z, the impedance (transfer function) for each frequency w of interest
	def importAllZ(self, TFsitefn):
		path = TFsitefn
		self.ZXX = []
		self.ZXY = []
		self.ZYX = []
		self.ZYY = []		
		for i in range(0,len(self.rows)):
			row=self.rows[i]
			col=self.cols[i]
			[ZXX,ZXY,ZYX,ZYY] = self.importZatw(path,row,col)

			self.ZXX = self.ZXX + [ZXX]
			self.ZXY = self.ZXY + [ZXY]
			self.ZYX = self.ZYX + [ZYX]
			self.ZYY = self.ZYY + [ZYY]
	

	#imports Z at only the frequency of interest specified by row, col
	def importZatw(self,path, row, column):
		file = open(path, 'r') 
		count = 0
		readnextline=False
		datanames=['ZXXR','ZXXI','ZXYR','ZXYI','ZYXR','ZYXI','ZYYR','ZYYI']

		u0 = 1.25664e-6 #m kg s^-2 A^-2 (permeability of free space, roughly same value inside the earth)


		#wait until get to impedances
		while True: 
			count += 1
		  
			# Get next line from file 
			line = file.readline()

			# if line is empty 
			# end of file is reached 
			if not line: 
				break

			# double check w is close to the frequency in the data 
			if(line.find('****FREQUENCIES****')>-1):
				line = file.readline() #skip '****FREQUENCIES****' line
				line = file.readline() #skip '>FREQ //30' line
				for r in range(row):
					line = file.readline() #skip rows until arrive at correct one
				
				allvals = line.split()

				f_in_file = float(allvals[column])
				# print('Frequency chosen in file '+path+' (Hz): '+str(f_in_file))
				self.wInEDI = self.wInEDI + [2*np.pi*f_in_file]

			if(line.find('****IMPEDANCES****')>-1):
				break
		  
		while True: 
			count += 1
		  
			# Get next line from file 
			line = file.readline()

			# if line is empty 
			# end of file is reached 
			if not line: 
				break

			if(readnextline):
				for r in range(row):
					line = file.readline() #skip rows until arrive at correct one

				readnextline=False
				allvals = line.split()

				Zi = float(allvals[column]) 
				
				# According to the documentation (https://ds.iris.edu/ds/products/emtf/):
				#   In practical magnetotellurics, the MT response tensor is colloqually known as the "MT impedance".
				#   Magnetotellurics employs the practical units for the MT response tensor, [mV/km]/[nT]. 

				# unit convert to kgs (standard SI units):  1e-3 mV/km/nT * (1V / 1e3mV) * (1km / 1e3m) * (1e9nT / 1T) = 1 V/m/T
				ziKGS = Zi/1e3

				if(matches[0]=='ZXXR'):
					ZXXR = ziKGS
				elif(matches[0]=='ZXXI'):
					ZXXI = ziKGS
				elif(matches[0]=='ZXYR'):
					ZXYR = ziKGS
				elif(matches[0]=='ZXYI'):
					ZXYI = ziKGS
				elif(matches[0]=='ZYXR'):
					ZYXR = ziKGS
				elif(matches[0]=='ZYXI'):
					ZYXI = ziKGS
				elif(matches[0]=='ZYYR'):
					ZYYR = ziKGS
				elif(matches[0]=='ZYYI'):
					ZYYI = ziKGS

			matches = [x for x in datanames if x in str(line)]

			if(len(matches)>0):
				readnextline=True

		file.close()

		return [ZXXR+1.0j*ZXXI,ZXYR+1.0j*ZXYI,ZYXR+1.0j*ZYXI,ZYYR+1.0j*ZYYI]

	def calcApparentc(self):
		self.apparentc = []
		ZXX = self.ZXX
		ZXY = self.ZXY
		ZYX = self.ZYX
		ZYY = self.ZYY
		w = self.wProvided
		for i in range(0,len(w)):
			# https://mathworld.wolfram.com/FrobeniusNorm.html
			# The frobenius norm involves absolute value squared of matrix elements, adding, then taking the square root.
			sumsquared = ZXX[i]*np.conj(ZXX[i])+ZXY[i]*np.conj(ZXY[i])+ZYX[i]*np.conj(ZYX[i])+ZYY[i]*np.conj(ZYY[i])

			# Love et al Equation 5: apparentc = w*u0/|Zi|^2.
			if(sumsquared==0):
				self.apparentc = self.apparentc + [np.inf]
			else:
				self.apparentc = self.apparentc + [w[i]*u0/sumsquared]

	#love et al 2019 equation 6
	def calcApparentr(self):
		self.rhoXX = []
		self.rhoXY = []
		self.rhoYX = []
		self.rhoYY = []
		ZXX = self.ZXX
		ZXY = self.ZXY
		ZYX = self.ZYX
		ZYY = self.ZYY
		w = self.wProvided
		for i in range(0,len(w)):
			# https://mathworld.wolfram.com/FrobeniusNorm.html
			# The frobenius norm involves absolute value squared of matrix elements, adding, then taking the square root.
			self.rhoXX = self.rhoXX + [ZXX[i]*np.conj(ZXX[i])/(2*np.pi*u0*self.f[i])]
			self.rhoXY = self.rhoXY + [ZXY[i]*np.conj(ZXY[i])/(2*np.pi*u0*self.f[i])]
			self.rhoYX = self.rhoYX + [ZYX[i]*np.conj(ZYX[i])/(2*np.pi*u0*self.f[i])]
			self.rhoYY = self.rhoYY + [ZYY[i]*np.conj(ZYY[i])/(2*np.pi*u0*self.f[i])]

	def plotApparentr(self):
		plt.figure(29)
		plt.loglog()
		# interpolatedfun=interp1d(self.f,self.rhoXX)
		# ynew = interpolatedfun(xnew)

		print(self.f[0])
		print(self.f[-1])

		# plt.plot(self.f,self.rhoXX,'k-', lw=1)
		plt.plot(self.f,self.ZXX*np.conj(self.ZXX), lw=1)
		plt.show()

	def getApparentcCloseToWindowPeriod(self,T):
		mindisti=1e19
		mini=0
		for i in range(0,len(self.T)):
			dist=np.abs(self.T[i]-T)
			if(dist<mindisti):
				mindisti=dist
				mini=i
		print('mindisti')
		print(mindisti)
		return self.apparentc[mini]