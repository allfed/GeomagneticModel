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
		self.maxchunksize=Params.maxchunksize
		self.siteIndex=siteIndex
		self.sitename=Params.mtsitenames[siteIndex]
		self.betaThreshold=Params.betaThreshold[siteIndex]
		self.MTsiteChunk=recordtype('MTsiteChunk','chunkindex chunksize rawNS rawEW BfieldNS BfieldEW EfieldNS EfieldEW absE rateperyearxy storms stormpeaks')
		self.ds = []
		self.powerfits=[]
		self.maglat=0
		self.N=0
		self.nchunks=0
		self.maxE=0
		self.polyFitFunNS=[]
		self.polyFitFunEW=[]

	def importSite(self):
		self.ds = nc.Dataset(self.MTsitefn)	#  import the magnetic field record using the netcdf format as a big numpy array
		self.maglat=self.ds['mlat'][-1]
		
		self.N=len(self.ds['dbn_geo'])#size of this particular MTsite magnetic field  (B field) record (in this case north, but north and south are always the same length)
		self.sampleperiod = Params.sampleperiod[self.siteIndex]
		self.hourWindow=int(np.floor((60*60)/self.sampleperiod))
		self.cumulativeyears=self.N*self.sampleperiod/secondsperyear #to be updated when processing ratesperyear
		if(Params.nchunks[self.siteIndex]):
			self.nchunks=Params.nchunks[self.siteIndex]
		else:
			self.nchunks=int(np.floor(self.N/self.maxchunksize))+1
			print('self.nchunks')
			print(self.nchunks)
		for i in range(0,self.nchunks):
			self.chunks=self.chunks+[self.MTsiteChunk(0,0,[],[],[],[],[],[],[],[], [],[])]

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

	# #MTsites are usually so large, one must process them in smaller chunks for the fourier transform and convolution with the TF site frequency dependent transfer function to succeed. We also need a TF site for each MT site to determine the transfer function and thus estimate the geoelectric field.
	def createChunk(self,chunkindex):
		minindex = chunkindex*self.maxchunksize
		maxindex = min((chunkindex+1)*self.maxchunksize-1,self.N-1)

		chunksize= maxindex-minindex+1

		print('importing chunkindex '+str(chunkindex)+', (chunk '+str(chunkindex+1)+' of '+str(self.nchunks)+')', end='\r')

		#getBfield along maglat NS and EW
		rawNS = self.ds['dbn_geo'][minindex:maxindex+1]
		rawEW = self.ds['dbe_geo'][minindex:maxindex+1]

		self.chunks[chunkindex] = self.MTsiteChunk(chunkindex,chunksize,rawNS,rawEW,[],[],[],[],[],[],[],[])

	def cleanBfields(self):
		for chunkindex in range(0,self.nchunks):
			self.cleanChunkBfields(chunkindex)

	def cleanChunkBfields(self,chunkindex):
		chunk=self.chunks[chunkindex]
		rawNS=chunk.rawNS
		rawEW=chunk.rawEW

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

		maskNS = np.isfinite(rawNS)
		maskEW = np.isfinite(rawEW)

		BfieldNS = np.interp(indices, indices[maskNS], rawNS[maskNS])
		BfieldEW = np.interp(indices, indices[maskEW], rawEW[maskEW])
		self.chunks[chunk.chunkindex].BfieldNS=BfieldNS
		self.chunks[chunk.chunkindex].BfieldEW=BfieldEW


	def plotEfields(self,chunkindex,startindex,endindex):
		print('loading field from chunkindex '+str(chunkindex)+', chunk '+str(chunkindex+1)+' of '+str(self.nchunks), end='\r')
		print('')
		
		Efields=self.chunks[chunkindex].absE[startindex:endindex]
		startindexds=chunkindex*self.maxchunksize+startindex
		endindexds=chunkindex*self.maxchunksize+endindex

		dts=[str(datetime(2018, 1,1))]*(endindexds-startindexds)
		elapsedHours=[0]*(endindexds-startindexds)

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
			print('loading field from chunkindex '+str(i)+', chunk '+str(i+1)+' of '+str(self.nchunks), end='\r')
			print('')
			chunkE=np.load(Params.mtEfieldsloc+str(self.sitename)+'c'+str(i)+'.npy')
			print('chunkstart: '+str(i*self.maxchunksize))
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
		print('calcEratesPerYear')
		self.calcStorms()
		self.windowedRates=[]
		savearr = []
		for windex in range(0,len(self.windows)):
			self.calcWindowEratesPerYear(windex)


	#combine the E rates calculated for all the chunks by windowed average into an array of e fields and a y axis value for rate per year of those e fields. This is a modification of (Love, 2018 page 9) which describes the overall process.
	def calcWindowEratesPerYear(self,windex):
		window=self.windows[windex]
		allcumsumE=[]
		allcountsatE=[]
		samplecount=0
		for j in range(0,len(self.chunks)):
			# get the occurrence of High E rates per year vs E field level
			rateperyearxy=self.calcChunkEratesPerYear(windex,j)
			allcumsumE = allcumsumE+rateperyearxy[0]
			allcountsatE =allcountsatE+ rateperyearxy[1]
			samplecount = samplecount + self.chunks[j].chunksize

		#sort all the highest E fields highest to lowest for all the chunks combined
		sortedcountsbyE = [x for _,x in sorted(zip(-np.array(allcumsumE),allcountsatE))]
		sortedE = -np.sort(-np.array(allcumsumE))

		#total years measured
		self.cumulativeyears = samplecount*self.sampleperiod/secondsperyear #total records (once per minute) by minutes in a year

		#finally combine the counts and take the rate
		cumsum=0
		cumsumArr=[]
		uniqueEarr=[]
		prevE=0
		
		for i in range(0,len(sortedE)):
			E=sortedE[i] 
			countsthisE=sortedcountsbyE[i]
			cumsum = cumsum+countsthisE
			if(E==prevE):
				#increment the count for this E field value
				cumsumArr[-1]=cumsumArr[-1]+countsthisE
				continue
			cumsumArr = cumsumArr + [cumsum]
			uniqueEarr = uniqueEarr + [E] 
			prevE=E

		rates=np.array(cumsumArr)/self.cumulativeyears


		windowalreadyexists=False

		#add this duration data to the existing windowedRates, or update if exists already
		for i in range(0,int(np.floor(len(self.windowedRates))/3)):
			durationindex = i*3
			efieldindex = i*3+1
			ratesindex = i*3+2

			duration=self.windowedRates[durationindex]
			if(duration==window*self.sampleperiod):
				self.windowedRates[efieldindex]=uniqueEarr
				self.windowedRates[ratesindex]=rates
				windowalreadyexists=True
				break

		if(not windowalreadyexists):
			self.windowedRates=self.windowedRates+[window*self.sampleperiod,uniqueEarr,rates]

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
			

	def plotEratesPerYear(self):
		plt.figure()
		plt.loglog()
		for i in range(0,int(np.floor(len(self.windowedRates))/3)):
			durationindex = i*3
			efieldindex = i*3+1
			ratesindex = i*3+2
			duration = self.windowedRates[durationindex]
			Efields = self.windowedRates[efieldindex]
			rates = self.windowedRates[ratesindex]
			# [slope,exponent]=fits.fitPower(Efields,rates)
			#if we've calculated power fits for all the windows
			if(len(self.powerfits)>0):

				slope=self.powerfits[i][1]
				exponent=self.powerfits[i][2]
				print('coeffs')
				print([slope,exponent])

				plt.plot(Efields,rates, Efields, fits.powerlaw(Efields,slope,exponent), lw=1,label = "Field averaged over "+str(duration)+" seconds")
			else:
				plt.plot(Efields,rates,lw=1,label = "Field averaged over "+str(duration)+" seconds")
			
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
	def calcStorms(self):
		for i in range(0,self.nchunks):
			self.calcChunkStorms(i)

		#combine any storms that occurred between chunks 
		for i in range(0,len(self.chunks)):
			if(i>0):
				#if first element is part of a group
				if(self.chunks[i].storms[0][0][0]==0):
				#combine this storm with the one at the end of the previous group if the last element of the previous group is also part of a storm
					if(self.chunks[i-1].storms[-1][0][-1]==(self.chunks[i-1].chunksize-self.hourWindow)):
						print('matches.')
						print('chunk'+str(i))
						stormextraminutes=len(self.chunks[i].storms[0][0])
						stormstart=self.chunks[i-1].storms[-1]		
						newindices=np.arange(stormstart[0][0],stormstart[0][-1]+stormextraminutes)
						stormtoappend=self.chunks[i].storms[0][1]
						self.chunks[i-1].storms[-1][0]=newindices
						self.chunks[i-1].storms[-1][1]=np.append(stormstart[1],stormtoappend)
						for j in range(0,len(self.nwindows)):
							self.chunks[i-1].stormpeaks[j][-1]=np.max([self.chunks[i-1].stormpeaks[j][-1],self.chunks[i].stormpeaks[j][0]])
						self.chunks[i].storms[0].pop()
						# print(len(newindices))
						# print(len(self.chunks[i-1].storms[-1][1]))
						# else:
						# print('no matches')
						# print('self.chunks[i-1].storms[-1][0][-1]')
						# print(self.chunks[i-1].storms[-1][0][-1])
						# print('self.chunks[i-1].chunksize-self.hourWindow')
						# print(self.chunks[i-1].chunksize-self.hourWindow)
		self.maxE=[0]*len(self.windows)
		for i in range(0,len(self.chunks)):
			for j in range(0,len(self.windows)):
				maxE=np.max(self.chunks[i].stormpeaks[j])
				if(maxE>self.maxE[j]):
					self.maxE[j]=maxE
		print('self.maxE')
		print(self.maxE)


		#remove any storms that are below 1/10 the peak E field for all measurements (remove different storms for each window), but the original storm E fields are only kept for window0.
		nstorms=0
		nwindows=len(self.windows)
		for i in range(0,len(self.chunks)):
			stormstmp=[]
			storms=self.chunks[i].storms
			peakstmp=[np.array([])]*nwindows
			for k in range(0,nwindows):
				sp=self.chunks[i].stormpeaks[k]
				for j in range(0,len(sp)):
					if(sp[j]>self.maxE[k]/10):
						peakstmp[k]=np.append(peakstmp[k],np.array(sp[j]))
						if(k==0):
							stormstmp.append(storms[j])
			self.chunks[i].stormpeaks=peakstmp
			self.chunks[i].storms=stormstmp
			print('num storms in chunk'+str(i))
			print(len(self.chunks[i].storms))
			nstorms=nstorms+len(self.chunks[i].storms)
		print('number storms per year: '+str(nstorms/self.cumulativeyears))

	def calcChunkStorms(self,chunkindex):
		chunk=self.chunks[chunkindex]
		#we lose nchunks hours out of entire dataset (one hour at the end of each chunk) by using windowed average, but this should be fine. 
		hourAvg = np.convolve(chunk.absE,np.ones(self.hourWindow)/self.hourWindow,mode='valid')
		# print('chunk.absE')
		# plt.figure()
		# plt.yscale("log")
		# plt.plot(chunk.absE[100000:1000100])
		# plt.plot(hourAvg[100000:1000100])
		# plt.show()
		# self.plotEfields(chunkindex,732000-5*60-20+12*60,732000-5*60-20+108*60)
		nwindows=len(self.windows)
		Efields=[np.array([])]*nwindows
		for i in range(0,nwindows):
			 smoothed = np.convolve(chunk.absE,np.ones(self.windows[i])/self.windows[i],mode='valid')
			 Efields[i] = smoothed[0:-self.hourWindow+1]

		exceeds=np.array(hourAvg)>self.betaThreshold
		data = zip(exceeds,Efields[0])
		Eindex=0
		storms=[]
		stormpeaks=[np.array([])]*nwindows

		print('getting all the storms in chunk '+str(chunkindex))

		for key, group in groupby(data, lambda x: x[0]):
			isstorm,field = next(group)
			elems = len(list(group)) + 1
			stormE=[np.array([])]*nwindows
			stormpeak=[0]*nwindows
			if(isstorm):
				Eindices=np.arange(Eindex,Eindex+elems)
				for i in range(0,nwindows):
					stormE[i]=Efields[i][Eindex:Eindex+elems]
					stormpeak[i]=np.max(np.array(stormE[i]))
					stormpeaks[i]=np.append(stormpeaks[i],stormpeak[i])

				storms.append([Eindices,stormE[0]])
			Eindex+=elems
		

		self.chunks[chunkindex].storms=storms
		self.chunks[chunkindex].stormpeaks=stormpeaks

	# determine number of times per year the field level or higher occurs in a chunk. Limit to a field which occurs less than 100 times per year for the chunk timespan
	def calcChunkEratesPerYear(self,windex,chunkindex):
		chunk=self.chunks[chunkindex]
		windowedE = np.array(chunk.stormpeaks[windex])
		chunkyears = chunk.chunksize*self.sampleperiod/secondsperyear #total recorded years

		##this was some code to test it works properly
		# if(chunkindex==0):
		# 	sortedcountsbyE=(np.array([3E-2,1])*self.cumulativeyears).astype(int)
		
		# 	print('sortedcountsbyE')
		# 	print(sortedcountsbyE)
		# 	sortedE=np.concatenate((np.repeat(3E-2,sortedcountsbyE[0]),np.repeat(6E-3,sortedcountsbyE[1])))
		# 	print('sortedE')
		# 	print(sortedE)


		sortedE=-np.sort(-windowedE)#sort largest E to smallest E

		cumsum=0
		cumsumarr = []
		countsatE = []
		uniqueEarr=[]
		thisEcount = 1
		prevE=0 
		numtimesexceedsprev = 0


		for E in sortedE:
			cumsum = cumsum+1

			#if this E field occurs more than 100 samples per year, there is no need to include more instances of it for the final combination of counts, as we're interested only in fields that occur once per year or fewer.
			if(E==prevE):
				#increment the count for this E value
				countsatE[-1]=countsatE[-1]+1
				thisEcount = thisEcount+1
				continue

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
		# plt.loglog()
		# print('sortedZfreqs')
		# plt.plot(sortedZFreqs,sortedZbyFreq)
		# plt.show()
		toInterpolateFreqs = []
		toInterpolateZ = []

		#we only use well-defined frequencies from 10^-1 to 10^-4. All frequencies outside of that are set to zero impedance.
		bandpass=np.logical_and(10**-1>np.array(sortedZFreqs),np.array(sortedZFreqs)>10**-4)
		
		bandpassedZ=np.array(sortedZbyFreq)*bandpass
		
		if(sortedZFreqs[0]>10**-4 or sortedZFreqs[-1]<10**-1):
			print('ERROR!!!: the TFsite has an unsuitable freqency range') 
			quit()

		#Furthermore, we need to ensure there is no amplification of frequencies higher than the half sampling rate (nyquist limit). Any Z values higher than this frequency is set to zero impedance.

		lowpass=0.5/self.sampleperiod>np.array(sortedZFreqs)
		lowpassedZ=np.array(sortedZbyFreq)*lowpass

		# if freq is lower,  
		if(sortedBFreqs[0] < sortedZFreqs[0]):
			toInterpolateFreqs = [sortedBFreqs[0]]
			toInterpolateZ =  np.array([0])#[sortedZbyFreq[0]]

		toInterpolateFreqs =np.append(toInterpolateFreqs, sortedZFreqs)
		toInterpolateZ = np.append(toInterpolateZ, lowpassedZ)

		# if freq is higher,  
		if(sortedBFreqs[-1] > sortedZFreqs[-1]):
			toInterpolateFreqs =np.append(toInterpolateFreqs, [sortedBFreqs[-1]])
			toInterpolateZ = np.append(oInterpolateZ, np.array([0]))
		
		# plt.show()
		# plt.figure()
		# plt.loglog()
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