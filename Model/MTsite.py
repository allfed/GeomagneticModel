# A MTsite refers to a MT historical magnetic field record made by a magnetotelluric measurement site.
# This MT record of the B field over time is converted to a field record using the transfer functions determined in TFsite class
# Note: to find MT station fields over time, need to use https://supermag.jhuapl.edu/

from __future__ import print_function
import numpy as np
import scipy

from scipy import stats as scistats
from scipy.fftpack import fft, ifft
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from itertools import groupby
import recordtype
from recordtype import recordtype
import glob 
from datetime import datetime, timedelta

import fits
import Params

u0 = 1.25664e-6 #m kg s^-2 A^-2 (permeability of free space, roughly same value inside the earth)
#also has units N A^-2. 
secondsperyear=31556925.216 # s (seconds in a solar year)

class MTsite:
	def __init__(self,MTsitefn,siteIndex):
		Params.importIfNotAlready()
		self.MTsitefn = MTsitefn 
		self.chunks=[]
		self.windowedRates=[]
		self.windowedRatesTmp=[] #temporary stored here before saving
		self.sitename=Params.mtsitenames[siteIndex]
		self.maxchunksize=Params.maxchunksize
		self.siteIndex=siteIndex
		self.MTsiteChunk=recordtype('MTsiteChunk','chunkindex chunksize rawNS rawEW BfieldNS BfieldEW EfieldNS EfieldEW absE rateperyearxy')
		self.ds = []
		self.powerfits=[]
		self.maglat=0
		self.N=0
		self.nchunks=0
		self.polyFitFunNS=[]
		self.polyFitFunEW=[]

	def importSite(self):
		# chunk1 = MTsiteChunk([1,2,3],[2,3],[1,2,3],[2,3],[],[])
		# chunk2 = MTsiteChunk([1,2,3],[2,3],[1,2,3],[2,3,33],[],[])
		# self.chunks=[chunk1,chunk2]
		# print('chunk.absE')
		# print(self.chunks[0].BfieldE)
		self.ds = nc.Dataset(self.MTsitefn)	#  import the magnetic field record using the netcdf format as a big numpy array
		self.maglat=self.ds['mlat'][-1]
		
		self.N=len(self.ds['dbn_geo'])#size of this particular MTsite magnetic field  (B field) record (in this case north, but north and south are always the same length)
		self.sampleperiod = Params.sampleperiod[self.siteIndex]
		self.cumulativeyears=self.N*self.sampleperiod/secondsperyear #to be updated when processing ratesperyear
		if(Params.nchunks[self.siteIndex]):
			self.nchunks=Params.nchunks[self.siteIndex]
		else:
			self.nchunks=int(np.floor(self.N/self.maxchunksize))+1
			print('self.nchunks')
			print(self.nchunks)
		for i in range(0,self.nchunks):
			self.chunks=self.chunks+[self.MTsiteChunk(0,0,[],[],[],[],[],[],[],[])]
		self.windows=[]
		self.windowedRates=[]
		for i in range(0,len(Params.windows[self.siteIndex])):
			windowstring=Params.windows[self.siteIndex][i]
			if not windowstring:
				break
			window = int(windowstring)
			self.windows=self.windows+[window]
			self.windowedRates=self.windowedRates+[window*self.sampleperiod,[],[]]

		self.calcPolyFits()

	#MTsites are usually so large, one must process them in smaller chunks for the fourier transform and convolution with the TF site frequency dependent transfer function to succeed. We also need a TF site for each MT site to determine the transfer function and thus estimate the geoelectric field.
	def createChunks(self):

		for chunkindex in range(0,self.nchunks):
			newchunk=self.createChunk(chunkindex)
		print('')
		print('')

	def createChunk(self,chunkindex):
		minindex = chunkindex*self.maxchunksize
		print('minindex')
		print(minindex)
		maxindex = min((chunkindex+1)*self.maxchunksize-1,self.N-1)
		print('maxindex')
		print(maxindex)

		chunksize= maxindex-minindex+1
		print('chunksize')

		print(chunksize)
		print('importing chunkindex '+str(chunkindex)+', (chunk '+str(chunkindex+1)+' of '+str(self.nchunks)+')', end='\r')
		# print(minindex)
		# print('maxindex')
		# print(maxindex)
		# print('chunksize')
		# print(chunksize)

		#getBfield along maglat NS and EW
		rawNS = self.ds['dbn_geo'][minindex:maxindex+1]
		rawEW = self.ds['dbe_geo'][minindex:maxindex+1]
		print('len(rawNS)')
		print(len(rawNS))

		self.chunks[chunkindex] = self.MTsiteChunk(chunkindex,chunksize,rawNS,rawEW,[],[],[],[],[],[])

	# #MTsites are usually so large, one must process them in smaller chunks for the fourier transform and convolution with the TF site frequency dependent transfer function to succeed. We also need a TF site for each MT site to determine the transfer function and thus estimate the geoelectric field.
	# def calcEfields(self,tfsites):
	# 	for i in range(0,len(Params.windows[self.siteIndex])):
	# 		tfsite=tfsites[self.siteIndex]
	# 		window = Params.windows[self.siteIndex][i]

	# 		if not window: #if window is empty, dont't process the MTsite at this window averaging
	# 			continue

	# 		for chunk in self.chunks:

	# #find the chunk and index for a given time (time in seconds) 
	# findEfieldindices():


	def cleanBfields(self):
		for chunkindex in range(0,self.nchunks):
			self.cleanChunkBfields(chunkindex)

	def cleanChunkBfields(self,chunkindex):
		chunk=self.chunks[chunkindex]
		print('lenrawns')		
		rawNS=chunk.rawNS
		rawEW=chunk.rawEW
		print(len(rawNS))

		if(np.isnan(rawNS[0])):
			rawNS[0] = 0
		if(np.isnan(rawNS[chunk.chunksize-2])):
			rawNS[chunk.chunksize-2] = 0

		if(np.isnan(rawEW[0])):
			rawEW[0] = 0
		if(np.isnan(rawEW[chunk.chunksize-2])):
			rawEW[chunk.chunksize-2] = 0

		#linearly interpolate any missing "nan" values

		indices = np.arange(len(rawNS))
		print('indices')
		print(indices)

		maskNS = np.isfinite(rawNS)
		maskEW = np.isfinite(rawEW)

		print('indices[maskNS]')
		print(indices[maskNS])

		BfieldNS = np.interp(indices, indices[maskNS], rawNS[maskNS])
		BfieldEW = np.interp(indices, indices[maskEW], rawEW[maskEW])
		self.chunks[chunk.chunkindex].BfieldNS=BfieldNS
		self.chunks[chunk.chunkindex].BfieldEW=BfieldEW
		print('BfieldNS')
		print(BfieldNS)


	def plotEfields(self,chunkindex,startindex,endindex):
		# chunkE
		print('loading field from chunkindex '+str(chunkindex)+', chunk '+str(chunkindex+1)+' of '+str(self.nchunks), end='\r')
		print('')
		# chunkE=np.load(Params.mtEfieldsloc+str(self.sitename)+'c'+str(chunkindex)+'.npy')
		# print('chunksize')
		# print(len(chunkE))
		
		Efields=self.chunks[chunkindex].absE[startindex:endindex]
		startindexds=chunkindex*self.maxchunksize+startindex
		endindexds=chunkindex*self.maxchunksize+endindex

		dts=[str(datetime(2018, 1,1))]*(endindexds-startindexds)
		elapsedHours=[0]*(endindexds-startindexds)
		print('endindexds')
		print(endindexds)
		print('startindexds')
		print(startindexds)

		years=np.array(self.ds['time_yr'][startindexds:endindexds]).astype(int)
		months=np.array(self.ds['time_mo'][startindexds:endindexds]).astype(int)
		days=np.array(self.ds['time_dy'][startindexds:endindexds]).astype(int)
		hours=np.array(self.ds['time_hr'][startindexds:endindexds]).astype(int)
		minutes=np.array(self.ds['time_mt'][startindexds:endindexds]).astype(int)
		seconds=np.array(self.ds['time_sc'][startindexds:endindexds]).astype(int)
			
		for i in range(0,len(years)):
			dts[i]=datetime(years[i],months[i],days[i],hours[i],minutes[i],seconds[i])
			elapsedHours[i]=(dts[i]-dts[0]).total_seconds()/(60*60)+12

		print('elapsedHours '+str(elapsedHours[0]))
		print('start time '+str(dts[0]))
		print('end time '+str(dts[-1]))
		plt.figure()
		plt.yscale("log")
		plt.plot(elapsedHours,Efields)
		# beautify the x-labels
		plt.gcf().autofmt_xdate()
		plt.show()

	def saveChunkEfields(self,chunkindex):
		np.save(Params.mtEfieldsloc+str(self.sitename)+'c'+str(chunkindex),self.chunks[chunkindex].absE)

	def loadEfields(self):
		for i in range(0,self.nchunks):
			# chunkE
			print('loading field from chunkindex '+str(i)+', chunk '+str(i+1)+' of '+str(self.nchunks), end='\r')
			print('')
			chunkE=np.load(Params.mtEfieldsloc+str(self.sitename)+'c'+str(i)+'.npy')
			print('chunkstart: '+str(i*self.maxchunksize))
			# print('chunksize')
			# print(len(chunkE))
			self.chunks[i].absE=chunkE
			self.chunks[i].chunksize=len(chunkE)

	#MTsites are usually so large, one must process them in smaller chunks for the fourier transform and convolution with the TF site frequency dependent transfer function to succeed. We also need a TF site for each MT site to determine the transfer function and thus estimate the geoelectric field.
	def calcEfields(self,TFsite):
		for i in range(0,len(self.nchunks)):
			self.calcChunkEfields(TFsite,i)
		print('')
		print('')

	# calculate 2nd order polynomial fits to NS and EW B data (downsample to one in 1000 points), assign array of fitted values to poly fit property
	def calcPolyFits(self):
		downsampleratio=1000
		samplesNS=np.array(self.ds['dbn_geo'][0:-1:downsampleratio])
		samplesEW=np.array(self.ds['dbe_geo'][0:-1:downsampleratio])

		indices = np.arange(len(samplesNS))
		maskNS = np.isfinite(samplesNS)
		maskEW = np.isfinite(samplesEW)

		BfieldNS = np.interp(indices, indices[maskNS], samplesEW[maskNS])
		BfieldEW = np.interp(indices, indices[maskEW], samplesNS[maskEW])

		coeffsNS=np.polyfit(indices*downsampleratio,BfieldNS,2)
		coeffsEW=np.polyfit(indices*downsampleratio,BfieldEW,2)

		funNS=np.poly1d(coeffsNS)
		funEW=np.poly1d(coeffsEW)

		# plt.figure()
		# plt.plot(BfieldEW)
		# plt.plot(funEW(indices*downsampleratio))
		# plt.show()
		self.polyFitFunNS=funNS
		self.polyFitFunEW=funEW

	def calcChunkEfields(self,TFsite,chunkindex):
		chunk=self.chunks[chunkindex]
		chunksize=chunk.chunksize
		startindex=chunkindex*self.maxchunksize
		endindex=startindex+chunksize

		#see love, 2018 for the list of four corrections applied here.

		# first, subtract 2nd order polynomial fit

		indices = np.array(range(startindex,endindex))
		detrendedBN2=chunk.BfieldNS
		print('chunk.BfieldNSlen')
		print(len(chunk.BfieldNS))
		detrendedBNS=chunk.BfieldNS-self.polyFitFunNS(indices)
		detrendedBEW=chunk.BfieldEW-self.polyFitFunEW(indices)



		# see page 381, Love 2019 for an explanation of the math below

		#second, apply FFT with 64 bit precision (fast fourier transform the field into freqency space)
		ftNS = fft(detrendedBNS,chunksize)
		ftEW = fft(detrendedBEW,chunksize)

		halflength = int(np.floor(chunksize/2))
			
		#next, matrix multiply the tensor by the field
		#we do this by assigning each point on the b field freq array 
		#as some Z value determined by linear interpolation
		fs = 1.0/self.sampleperiod#sampling frequency of 1/60 Hz

		freqrange=[x*fs/chunksize for x in range(1, halflength)]
		periodrange=[chunksize/(x*fs) for x in range(1, halflength)]

		print('interpolating Z(w) for B(w), chunkindex '+str(chunk.chunkindex)+', '+str(chunk.chunkindex+1)+' of '+str(self.nchunks), end='\r')

		ZXXforB=self.getInterpolation(freqrange,TFsite.ZXX,TFsite.f)
		ZXYforB=self.getInterpolation(freqrange,TFsite.ZXY,TFsite.f)
		ZYXforB=self.getInterpolation(freqrange,TFsite.ZYX,TFsite.f)
		ZYYforB=self.getInterpolation(freqrange,TFsite.ZYY,TFsite.f)

 
		#frequency domain
		EfieldNS_fd = [x * y for x, y in zip(ZXXforB, ftNS)]+\
		[x * y for x, y in zip(ZXYforB, ftEW)]
		EfieldEW_fd = [x * y for x, y in zip(ZYXforB, ftNS)]+\
		[x * y for x, y in zip(ZYYforB, ftEW)]

		#time domain
		EfieldNS=ifft(EfieldNS_fd,chunksize)
		EfieldEW=ifft(EfieldEW_fd,chunksize)

		#finally, determine the magnitude of E field over time
		absE = np.real(np.sqrt(\
			EfieldNS*np.conj(EfieldNS)+\
			EfieldEW*np.conj(EfieldEW)))

		self.chunks[chunk.chunkindex].EfieldNS=EfieldNS
		self.chunks[chunk.chunkindex].EfieldEW=EfieldEW
		self.chunks[chunk.chunkindex].absE=absE
		self.chunks[chunk.chunkindex].chunksize=chunksize


	def calcEratesPerYear(self):
		self.windowedRates=[]
		savearr = []
		for w in self.windows:
			self.calcWindowEratesPerYear(w)


	#combine the E rates calculated for all the chunks by windowed average into an array of e fields and a y axis value for rate per year of those e fields. 
	def calcWindowEratesPerYear(self,window):
		allcumsumE=[]
		allcountsatE=[]
		samplecount=0
		for j in range(0,len(self.chunks)):
			# get the occurrence of High E rates per year vs E field level
			rateperyearxy=self.calcChunkEratesPerYear(window,j)
			allcumsumE = allcumsumE+rateperyearxy[0]
			allcountsatE =allcountsatE+ rateperyearxy[1]
			samplecount = samplecount + self.chunks[j].chunksize

		
		# int(np.floor(secondsperyear))
		# chunk.absE=np.array([4E-2])*int(np.floor(secondsperyear))
		# np.repeat(np.array([9E-3]),1)
		# np.repeat(np.array([4E-2]),1)

		#sort all the highest E fields highest to lowest for all the chunks combined
		sortedcountsbyE = [x for _,x in sorted(zip(-np.array(allcumsumE),allcountsatE))]
		sortedE = -np.sort(-np.array(allcumsumE))

		#total years measured
		self.cumulativeyears = samplecount*self.sampleperiod/secondsperyear #total records (once per minute) by minutes in a year
		# sortedE=np.array([3E-2,6E-3])
		# sortedcountsbyE=np.array([3E-2,1])*self.cumulativeyears
		#finally combine the counts and take the rate
		print('self.cumulativeyears')
		print(self.cumulativeyears)
		cumsum=0
		cumsumArr=[]
		uniqueEarr=[]
		prevE=0
		# plt.show()
		print('sortedE')
		print(sortedE)
		for i in range(0,len(sortedE)):
			E=sortedE[i] 
			countsthisE=sortedcountsbyE[i]
			cumsum = cumsum+countsthisE
			if(E==prevE):
				#increment the count for
				cumsumArr[-1]=cumsumArr[-1]+countsthisE
				continue
			cumsumArr = cumsumArr + [cumsum]
			uniqueEarr = uniqueEarr + [E] 
			prevE=E

			#if we are now looking at fields that occur 1 times a year or more, stop saving data
			if(cumsum/self.cumulativeyears >1):
				break
		rates=np.array(cumsumArr)/self.cumulativeyears
		self.windowedRates=self.windowedRates+[window*self.sampleperiod,uniqueEarr,rates]
		self.windowedRatesTmp=self.windowedRates

	# if E rates for year already exists for this site, add to it. Otherwise, create a new .npy to store Erates per year for this site
	def saveEratesPerYear(self):

		#add or replace the data for each window that has been created so far
		if(self.loadWindowedRates()):
			savearr=self.windowedRates
			print('self.windowedRates after loaded')
			print(self.windowedRates)
			
			loadeddurations=np.array([])
			for i in range(0,int(np.floor(len(self.windowedRates))/3)):
				durationindex = i*3
				efieldindex = i*3+1
				ratesindex = i*3+2

				duration = self.windowedRates[durationindex]
				print('duration')
				print(duration)
				loadeddurations=np.append(loadeddurations,duration)
				print('loadeddurations')
				print(loadeddurations)
				# Efields = self.windowedRates[efieldindex]
				# rates = self.windowedRates[ratesindex]

			print('self.windowedRatesTmp')
			print(self.windowedRatesTmp)
			for i in range(0,int(np.floor(len(self.windowedRatesTmp))/3)):
				durationindex = i*3
				efieldindex = i*3+1
				ratesindex = i*3+2

				duration = self.windowedRatesTmp[durationindex]
				Efields = self.windowedRatesTmp[efieldindex]
				rates = self.windowedRatesTmp[ratesindex]

				if(len(rates)!=0):#if data has been processed for this window
					loadeddurationindex=np.where(loadeddurations==duration)[0]
					
					if(len(loadeddurationindex)==0):
						print('nomatch')
						# print(type('windowedRates'))
						# print(type(self.windowedRates))
						# print(type('duration'))
						# print(type(duration))
						# print(type('Efields'))
						# print(type(Efields))
						# print(type('rates'))
						# print(type(rates))

						self.windowedRates=np.append(self.windowedRates,[duration,Efields,rates])
					else:
						print('matches')
						self.windowedRates[loadeddurationindex[0]*3] = duration
						print('Efields')
						print(Efields)
						print(type(Efields))
						print('self.windowedRates[1]')
						print(self.windowedRates[1])
						print('typeself.windowedRates')
						print(type(self.windowedRates[1]))
						self.windowedRates[loadeddurationindex[0]*3+1]=Efields
						self.windowedRates[loadeddurationindex[0]*3+2]=rates
		else:
			self.windowedRates=self.windowedRatesTmp
		print('savearr')
		print(self.windowedRates)
		np.save(Params.mtRepeatRatesDir+'MTsite'+str(self.siteIndex)+'EfieldRatesPerYear',self.windowedRates)


	def loadWindowedRates(self):
		loaddirectory=Params.mtRepeatRatesDir+'MTsite'
		allfiles=glob.glob(loaddirectory+'*.np[yz]')

		for f in allfiles:
			if(f==Params.mtRepeatRatesDir+'MTsite'+str(self.siteIndex)+'EfieldRatesPerYear.npy'):
				wr=np.load(Params.mtRepeatRatesDir+'MTsite'+str(self.siteIndex)+'EfieldRatesPerYear.npy')
				self.windowedRates=wr
				return True
		return False		
			

	def plotandFitEratesPerYear(self):
		plt.figure()
		plt.loglog()
		self.powerfits = []
		for i in range(0,int(np.floor(len(self.windowedRates))/3)):
			durationindex = i*3
			efieldindex = i*3+1
			ratesindex = i*3+2
			duration = self.windowedRates[durationindex]
			Efields = self.windowedRates[efieldindex]
			rates = self.windowedRates[ratesindex]
			[slope,exponent]=fits.fitPower(Efields,rates)
			print('coeffs')
			print([slope,exponent])
			plt.plot(Efields,rates, Efields, fits.powerlaw(Efields,slope,exponent), lw=1,label = "Field averaged over "+str(duration)+" seconds")
			# plt.plot(Efields,rates,lw=1,label = "Field averaged over "+str(duration)+" seconds")

			self.powerfits = self.powerfits + [[duration,slope,exponent]]
			
		# print('Efields')
		# print(Efields)
		# print('rates')
		# print(rates)

		plt.legend()
		plt.title('Rate geoelectric field is above threshold')
		plt.xlabel('Geoelectric Field (V/km)')
		plt.ylabel('Average rate per year distinct contiguous sample average is above E field (counts/year)')
	
		plt.show()

	def fitEratesPerYear(self):
		self.powerfits = []
		for i in range(0,int(np.floor(len(self.windowedRates))/3)):
			durationindex = i*3
			efieldindex = i*3+1
			ratesindex = i*3+2

			windowperiod = self.windowedRates[durationindex]
			Efields = self.windowedRates[efieldindex]
			rates = self.windowedRates[ratesindex]

			[slope,exponent]=fits.fitPower(Efields,rates)

			self.powerfits = self.powerfits + [[windowperiod,slope,exponent]]

			# plt.plot(Efields,rates, lw=1,label = "Field averaged over "+str(windowperiod)+" seconds")

	# def fitToRateperyear(self):
	# 	plt.figure()
	# 	plt.loglog()
	# 	for i in range(0,int(np.floor(len(self.windowedRates))/3)):
	# 		durationindex = i*3
	# 		efieldindex = i*3+1
	# 		ratesindex = i*3+2

	# 		windowperiod = self.windowedRates[durationindex]
	# 		Efields = self.windowedRates[efieldindex]
	# 		rates = self.windowedRates[ratesindex]
	# 	plt.legend()
	# 	plt.title('Rate geoelectric field is above threshold')
	# 	plt.xlabel('Geoelectric Field (V/km)')
	# 	plt.ylabel('Average rate per year distinct contiguous sample average is above E field (counts/year)')
	
	# 	plt.show()


	# determine number of times per year the field level or higher occurs in a chunk. Limit to a field which occurs less than 100 times per year for the chunk timespan
	def calcChunkEratesPerYear(self,window,chunkindex):
		chunk=self.chunks[chunkindex]
		if(window==1):
			windowedE = chunk.absE
		else: 
			windowedE = np.convolve(chunk.absE, np.ones(window)/window, mode='valid')
		# print('chunksize')
		# print(chunk.chunksize)
		# print('sampleperiod')
		# print(self.sampleperiod)
		chunkyears = chunk.chunksize*self.sampleperiod/secondsperyear #total recorded years
		# int(np.floor(secondsperyear))
		# chunk.absE=np.array([4E-2])*int(np.floor(secondsperyear))
		# np.repeat(np.array([9E-3]),1)
		# np.repeat(np.array([4E-2]),1)

		##this was some code to test it works properly
		# if(chunkindex==0):
		# 	sortedcountsbyE=(np.array([3E-2,1])*self.cumulativeyears).astype(int)
		# 	print('sortedcountsbyE')
		# 	print(sortedcountsbyE)
		# 	sortedE=np.concatenate((np.repeat(3E-2,sortedcountsbyE[0]),np.repeat(6E-3,sortedcountsbyE[1])))
		# 	print('sortedE')
		# 	print(sortedE)

		# else:
		# 	return [[],[]]


		sortedE=-np.sort(-windowedE)#sort largest E to smallest E

		cumsum=0
		cumsumarr = []
		countsatE = []
		uniqueEarr=[]
		thisEcount = 1
		prevE=0 
		numtimesexceedsprev = 0

		# print('sortedE')
		# print(sortedE)
		# print('chunkyears')
		# print(chunkyears)

		for E in sortedE:
			cumsum = cumsum+1

			#if this E field occurs more than 100 samples per year, there is no need to include more instances of it for the final combination of counts, as we're interested only in fields that occur once per year or fewer.
			if(cumsum/chunkyears >100):
				print('100 per year E for chunk'+str(chunkindex)+' and window of '+str(window)+':')
				print(E)
				break
			if(E==prevE):
				#increment the count for
				countsatE[-1]=countsatE[-1]+1
				thisEcount = thisEcount+1
				continue

			if(window!=1):

				#get groups of 1s and 0s where E field exceeds or meets threshold
				exceeds = (windowedE>=E)
				
				#count consecutive groups where windowed average is above threshold
				numtimesexceeds=sum(k for k,v in groupby(exceeds))
				cumsumarr = cumsumarr +[cumsum]

				countsatE = countsatE + [(numtimesexceeds-numtimesexceedsprev)]
				numtimesexceedsprev = numtimesexceeds
			else:
				countsatE = countsatE + [1]
			prevE=E
			uniqueEarr = uniqueEarr + [E] 
			thisEcount = 1
		# plt.figure()
		# plt.loglog()
		# plt.plot(uniqueEarr,countsatE,lw=1,label = "Field averaged over seconds")
		# plt.legend()
		# plt.title('Rate geoelectric field is above threshold')
		# plt.xlabel('Geoelectric Field (V/km)')
		# plt.ylabel('Average rate per year distinct contiguous sample average is above E field (counts/year)')
	
		# plt.show()
		# quit()
		return [uniqueEarr,countsatE]


	# The purpose of this function is to return an array of 
	# values of impedance for each frequency bin in the frequency
	# domain B field. This allows us to multiply the tensor 
	# components of the frequency domain B field with the 
	def getInterpolation(self,BfieldFreqs,Z,Zfreqs):
		sortedBFreqs=sorted(BfieldFreqs)
		sortedZFreqs=sorted(Zfreqs)
		# plt.figure()
		# plt.loglog()
		# print('sortedBfreqs')
		# plt.plot(sortedBFreqs)
		# plt.show()
		# plt.figure()
		# plt.loglog()
		# print('sortedZfreqs')
		# plt.plot(sortedZFreqs)
		# plt.show()

		sortedZbyFreq = [x for _,x in sorted(zip(Zfreqs,Z))]
		# plt.show()
		# plt.figure()
		# # plt.loglog()
		# print('sortedZfreqs')
		# plt.plot(sortedZFreqs,sortedZbyFreq)
		# plt.show()
		toInterpolateFreqs = []
		toInterpolateZ = []

		#we only use well-defined frequencies from 10^-1 to 10^-4. All frequencies outside of that are set to zero impedance.
		# toohigh=np.where(np.array(sortedZFreqs)>10**-1)
		# toolow=np.where(np.array(sortedZFreqs)<10**-4)
		
		bandpass=np.logical_and(10**-1>np.array(sortedZFreqs),np.array(sortedZFreqs)>10**-4)
		
		bandpassedZ=np.array(sortedZbyFreq)*bandpass
		
		print('bandpass')
		print(bandpass)
		if(sortedZFreqs[0]>10**-4 or sortedZFreqs[-1]<10**-1):
			print('ERROR!!!: the TFsite has an unsuitable freqency range') 
			quit()

		#Furthermore, we need to ensure there is no amplification of frequencies higher than the half sampling rate (nyquist limit). Any Z values higher than this frequency is set to zero impedance.

		lowpass=0.5/self.sampleperiod>np.array(sortedZFreqs)
		lowpassedZ=np.array(sortedZbyFreq)*lowpass

		#if the lowest B field frequency is lower than the lowest Z frequency, extrapolate the Z value of the lowest frequency Z to the same value at the lowest B field frequency
		# if(sortedBFreqs[0] < sortedZFreqs[0]):
		# 	print('islower')
		# 	toInterpolateFreqs = [sortedBFreqs[0]]
		# 	toInterpolateZ =  [sortedZbyFreq[0]]

		if(sortedBFreqs[0] < sortedZFreqs[0]):
			toInterpolateFreqs = [sortedBFreqs[0]]
			toInterpolateZ =  np.array([0])#[sortedZbyFreq[0]]
		print(len(sortedZFreqs))
		print(len(bandpassedZ))

		toInterpolateFreqs =np.append(toInterpolateFreqs, sortedZFreqs)
		toInterpolateZ = np.append(toInterpolateZ, lowpassedZ)
		print(len(toInterpolateZ))
		print(len(toInterpolateFreqs))

		# #if the highest B field frequency is higher than the highest Z frequency, extrapolate the Z value of the highest frequency Z to the same value at the highest B field frequency
		# if(sortedBFreqs[-1] > sortedZFreqs[-1]):
		# 	print('ishigher')
		# 	toInterpolateFreqs =toInterpolateFreqs + [sortedBFreqs[-1]]
		# 	toInterpolateZ = toInterpolateZ + [sortedZbyFreq[-1]]
		
		# if freq is higher,  
		if(sortedBFreqs[-1] > sortedZFreqs[-1]):
			print('ishigher')
			toInterpolateFreqs =np.append(toInterpolateFreqs, [sortedBFreqs[-1]])
			toInterpolateZ = np.append(oInterpolateZ, np.array([0]))
		
		# plt.show()
		# plt.figure()
		# # plt.loglog()
		# print('sortedZfreqs')
		# plt.plot(toInterpolateFreqs,toInterpolateZ)
		# plt.show()
		#return the interpolating function with Z values for every Bfield frequency
		interpolatedfun=interp1d(toInterpolateFreqs,toInterpolateZ)
		ZforBfield = interpolatedfun(BfieldFreqs)
		# plt.figure()
		# plt.loglog()
		# plt.plot(BfieldFreqs,ZforBfield)
		# print('plt.plot(BfieldFreqs,ZforBfield)')
		# plt.show()

		# plt.figure()
		# plt.loglog()
		# print('plt.plot(toInterpolateFreqs,toInterpolateZ)')
		# plt.plot(toInterpolateFreqs,np.real(toInterpolateZ))
		# plt.show()
		return ZforBfield