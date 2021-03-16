#powerfit adapted from https://scipy-cookbook.readthedocs.io/items/FittingData.html
#lognormal fit adapted from https://machinelearningmastery.com/curve-fitting-with-python/
import pandas as pd
import numpy as np
from numpy import pi, r_
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.optimize import curve_fit
import scipy.stats
from scipy.stats import lognorm,norm
from scipy.stats import kurtosis, skew
from scipy.stats import skewnorm

#  log((exponent-1)*x**(-exponent))
# =log(exponent-1)+log(x)/exponent
# we know this equals 1 at x[0]
#  1-log(exponent-1)=log(x[0])/exponent
# 
def powerlaw(x,exponent):

	return (x**exponent)/np.min(x)**exponent

# y=s*(x^e)
# y/s=x^e
# (y/s)^(1/e)=x
# def powerlawxfromy(y,slope,exponent):
# 	return (y/slope)**(1/exponent)

#see Love, 2018 page 9 (equation 7)
# upsilon is the expected value of the PDF.
# epsilonsqd is epsilon squared, the variance of the distribution
def lognormal(x,upsilon,epsilon):
	# return np.exp(-np.power((np.log(x)-upsilon),2))
	# return x/(2*x)
	numerator=np.exp(-np.power((np.log(x)-upsilon),2)/(2*epsilon**2))
	denominator=x*np.sqrt(2*np.pi*epsilon**2)
	return numerator/denominator

# incorporated scaling
# upsilon is the expected value of the PDF.
# epsilonsqd is epsilon squared, the variance of the distribution
def lognormal(x,upsilon,epsilon):
	# return np.exp(-np.power((np.log(x)-upsilon),2))
	# return x/(2*x)
	numerator=np.exp(-np.power((np.log(x)-upsilon),2)/(2*epsilon**2))
	denominator=x*np.sqrt(2*np.pi*epsilon**2)
	return numerator/denominator

def loclognormal(x,upsilon,epsilon,loc):
	# return np.exp(-np.power((np.log(x)-upsilon),2))
	# return x/(2*x)
	numerator=np.exp(-np.power((np.log(x-loc)-upsilon),2)/(2*epsilon**2))
	denominator=x*np.sqrt(2*np.pi*epsilon**2)
	return numerator/denominator

def importedlognormal(x,upsilon,epsilon):
	s=np.abs(epsilon)
	scale=np.exp(upsilon)
	y=np.abs(x)/np.abs(scale)
	return lognorm.pdf(y, s)/np.abs(scale)

def locimportedlognormal(x,upsilon,epsilon,loc):
	s=np.abs(epsilon)
	scale=np.exp(upsilon)
	y=(np.abs(x-loc))/scale
	y=(x-loc)/scale
	return lognorm.pdf(y, s)/scale

def logcdf(x,mu,sigma,loc):
	scale=np.exp(mu)
	return 1-lognorm.cdf(x,sigma,loc,scale)

# def logcdf(x,mu,sigma,loc=0):
# 	return 1-lognorm.cdf(x,mu,loc,sigma)
def cdf(x):
	return 1-norm.cdf(x)
def logpdf(x,mu,sigma,loc=0):
	return lognorm.pdf(x,mu,loc,sigma)
# def logcdfskewkurt(x,mean, variance, skew, kurtosis):
# 	return skewnorm.cdf(mean,)

# def fitLognormal(xdata,ydata):
# 	guessupsilon=-11
# 	guessepsilonsqd=4.6#50000
# 	# for guessupsilon in range(-20,-5):

# 	# we help curve_fit out by providing a wide range of reasonable guess inputs
# 	plt.figure()
# 	# for guessupsilon in range(-20,-5,3):
# 	# for guessepsilonsqd in range(1,9,2):
# 	print('guessepsilonsqd')
# 	print(guessepsilonsqd)
# 	print('guessupsilon')
# 	print(guessupsilon)
# 	plt.loglog(xdata,lognormal(np.array(xdata),guessupsilon,guessepsilonsqd)*1.1,lw=1,label = "guess for guess upsilon:"+str(guessupsilon)+" guessepsilonsqd"+str(guessepsilonsqd))
# 	try:
# 		params=[]
# 		params, covar = curve_fit(lognormal,xdata,ydata,p0=[guessupsilon,guessepsilonsqd],maxfev=10000)
# 		upsilon, epsilonsqd = params
# 		print(params) 
# 		quit()
# 		if(covar[0][0]==inf):
# 			print('fit infinite covar: guessupsilon: '+str(guessupsilon)+'  epsilonsqd: '+str(guessepsilonsqd))
# 			# continue

# 	except:
# 		print('fit failed: guessupsilon: '+str(guessupsilon)+'  epsilonsqd: '+str(guessepsilonsqd))
# 		print('params')
# 		print(params)
# 		upsilon=guessupsilon
# 		epsilonsqd=guessepsilonsqd
# 		# continue
# 		fitted=lognormal(np.array(xdata),upsilon,epsilonsqd)
# 		print(fitted)
# 		plt.loglog(xdata,fitted,lw=1,label = "fit from guess upsilon:"+str(guessupsilon)+" guessepsilonsqd"+str(guessepsilonsqd))

# 	plt.loglog(xdata,ydata,lw=1,label = "data we're fitting")
# 	plt.legend()
# 	print('covar')
# 	print(covar)
# 	plt.show()		

# 	# plt.loglog(xdata, lognormal(np.array(xdata),upsilon,epsilonsqd))
# 	# plt.loglog(xdata, ydata)
# 	# plt.xlabel('X (log scale)')
# 	# plt.ylabel('Y (log scale)')
# 	# plt.show()

# 	return [upsilon,epsilonsqd]


def fitLognormalCDF(xdata,ydata,guessmu,guesssigma,plot):
	try:
		params, covar = curve_fit(logcdf,xdata,ydata,p0=[guessmu,guesssigma,0],maxfev=10000)
	except:
		print('fit failed')
		quit()
	# mean,std=params
	# loc=0
	mean,std,loc=params
	print('mean')
	print(mean)
	print('std')
	print(std)
	print('loc')
	print(loc)

	if(plot):
		plt.figure()
		plt.loglog(xdata,ydata,lw=1,label = "data")
		fitdata=logcdf(np.array(xdata),mean,np.abs(std),loc)
		yguess=logcdf(np.array(xdata),guessmu,guesssigma,0)
		plt.loglog(xdata,yguess,lw=1,label = "guess dist")
		plt.loglog(xdata,fitdata,lw=1,label = "fit dist")
		plt.legend()
		plt.show()
	return [mean,std,loc]

def fitLognormalwithloc(xdata,ydata):
	# plt.figure(45)
	# plt.subplot(2, 1, 2)
	# guessupsilon=.0005
	# guessepsilonsqd=50000 
	guesses=[\
	np.array([-1.8908820193414216e-05, 0.003528574384792193, 0.9972651549201421]),
	np.array([5.9588582, 4.37993949, 0.00277132933]),\
	np.array([-1.1966543898728443, 0.04531822610045514, 0.27575209465565304]),\
	np.array([-1.1966543898728443, 0.4, 0.02]),\
	np.array([-3.07289437, 0.02616835, 0.04285109]),\
	np.array([-6.43616259e-04, 4.75736705e-03, 9.89547338e-01]),\
	]

	for g in guesses:
		print('guessg')
		print(g)
		guessdata=locimportedlognormal(np.array(xdata),g[0],np.abs(g[1]),g[2])
		# plt.loglog(xdata,guessdata,lw=1,label = "guess ["+str(g[0])+', '+str(g[0])+', '+str(g[0])+']')
		try:
			logx=np.log(xdata)
			logy=np.log(ydata)
			params, covar = curve_fit(locimportedlognormal,xdata,ydata,\
				p0=g,\
				maxfev=10000)

			upsilon,epsilon,loc=params
			fitdata=locimportedlognormal(np.array(xdata),upsilon,np.abs(epsilon),loc)

			r = np.divide(np.array(ydata),np.array(fitdata))#ratios
			arezeros=np.array(1-r).astype(bool)
			for e in arezeros:
				if(not e):
					print('error with one of the fits: zero value occured!')
					continue
				if(e==np.inf):
					print('error with one of the fits: infinitely large log residual!')
					continue
			avgr=np.mean(np.log(np.abs(r)))
			print('avgr')
			print(avgr)
			plt.loglog(xdata, ydata)
			plt.loglog(xdata, fitdata)
			plt.xlabel('X (log scale)')
			plt.ylabel('Y (log scale)')
		except:
			print('fitfailed, trying again')
			params=[-5.8,.37,0.0001]	
			r=np.inf
			continue
		break
	plt.legend()
	# params, covar = curve_fit(locimportedlognormal,xdata,ydata,p0=[-5.8,.37,0.00000001],maxfev=10000)
	# print('covar')
	# print(covar)
	print('final average residuals')
	print(avgr)
	# paramsa, covar = curve_fit(lognormal,xdata,ydata,p0=[-5.8,.37],maxfev=10000)
	# plt.figure()
	# plt.loglog(xdata, locimportedlognormal(np.array(xdata),upsilon,epsilon,loc))
	# plt.loglog(xdata, ydata)
	# plt.xlabel('X (log scale)')
	# plt.ylabel('Y (log scale)')
	plt.show()
	return [upsilon,epsilon,loc]

def fitPower(xdata,ydata):
	# print(xdata)
	# print(ydata)
	# Note: all positive, non-zero data
	xdata=np.array(xdata)
	ydata=np.array(ydata)
	# Define function for calculating a power law
	xdata=np.flip(xdata)
	ydata=np.flip(ydata)

	# Note: all positive, non-zero data
	exponentguess=2
	# plt.plot(xdata, ydata)  # Fit
	##########
	# Fitting the data -- Least Squares Method
	##########

	# Power-law fitting is best done by first converting
	# to a linear equation and then fitting to a straight line.
	# Note that the `logyerr` term here is ignoring a constant prefactor.
	#
	#  y = a * x^b
	#  log(y) = log(a) + b*log(x)
	#

	logx = np.log10(xdata)
	x0=logx[0]
	logy = np.log10(ydata)

	# define our (line) fitting function
	fitfunc = lambda p, x: p[0] * (x-x0)
	errfunc = lambda p, x, y: (y - fitfunc(p, x))

	pinit = [exponentguess]
	out = optimize.leastsq(errfunc, pinit,args=(logx, logy), full_output=1)

	pfinal = out[0]
	# print(pfinal)
	# print(covar)

	exponent = pfinal[0]

	plot=False
	if(plot):
		##########
		# Plotting data
		##########
		plt.figure(45)
		plt.clf()
		plt.subplot(2, 1, 1)
		plt.plot(xdata, ydata)     # Fit
		plt.plot(xdata, powerlaw(xdata, exponent))     # Fit
		plt.title('Best Fit Power Law')
		plt.xlabel('X')
		plt.ylabel('Y')

		plt.subplot(2, 1, 2)
		plt.loglog(xdata, ydata)
		plt.loglog(xdata, powerlaw(xdata, exponent))
		plt.xlabel('X (log scale)')
		plt.ylabel('Y (log scale)')
		plt.show()
	return [exponent]

#thanks to https://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy:
def weighted_std(values, weights):
	average = np.average(values, weights=weights)
	# Fast and numerically precise:
	variance = np.average((values-average)**2, weights=weights)
	return np.sqrt(variance)

def binnormaldist(dist,counts,sigmas):
	comb=pd.DataFrame({'x':dist,'counts':counts})
	withduplicates=comb.loc[comb.index.repeat(comb.counts)]
	countswithduplicates=np.array(withduplicates.as_matrix(columns=withduplicates.columns[0:]))[:,0]
	kurt=kurtosis(countswithduplicates)
	sk=skew(countswithduplicates)

	print('kurtosis')
	print(kurt)
	print('skew')
	print(sk)
	print('')
	nbins=50
	bins=[0]*nbins
	bins2=[0]*nbins
	mdist=np.mean(dist)
	sdist=np.std(dist)

	minimum=-sigmas*sdist+mdist
	maximum=sigmas*sdist+mdist

	binvalslow=[]
	binvalshigh=[]
	for i in range(0,nbins):
		low=minimum+(maximum-minimum)/nbins*i
		high=minimum+(maximum-minimum)/nbins*(i+1)
		mask=np.logical_and(dist>low,dist<high)
		bins[i]=np.sum(counts[mask])
		binvalslow.append(low)
		binvalshigh.append(high)
		# bins2[i]=np.sum(mask2)

	# plt.figure()
	# plt.yscale("log")
	# plt.plot(binvalslow,'.')
	# plt.plot(binvalshigh,'.')
	# plt.plot(binvalslow,bins)
	# plt.plot(binvalshigh,bins)
	# plt.plot(bins2)
	# plt.show()
	countsnormalized=np.array(bins)/np.sum(bins)
	return [np.divide(np.array(binvalslow)+np.array(binvalshigh),2),countsnormalized]

def getGuesses(E,counts,plot):
	normalE=np.log(E)
	mean=np.average(normalE,weights=counts)
	std=weighted_std(normalE,counts)
	if(plot):
		centereddist=(normalE-mean)/std
		[stormx,stormy]=binnormaldist(centereddist,counts,5)

		normaldist=np.random.normal(0,1,100000)#len(E))
		normalcounts=np.ones(len(normaldist))
		[normalx,normaly]=binnormaldist(normaldist,normalcounts,5)

		plt.figure()
		plt.plot(stormx,stormy,lw=1,label = "storm")
		plt.plot(normalx,normaly,lw=1,label = "normal dist")
		plt.legend()
		plt.show()

	return [mean,std]

def combinecounts(E,allcountsatE):
	sortedcountsbyE = [x for _,x in sorted(zip(-np.array(E),allcountsatE))]
	sortedE = np.sort(-np.array(E))
	negE=np.array(sortedE)
	
	combined=pd.DataFrame({'counts':sortedcountsbyE,'negE':negE}).groupby(by=["negE"]).sum()
	cumsumcombined=combined.cumsum()
	cumsumfinal=np.transpose(np.array(cumsumcombined.as_matrix(columns=cumsumcombined.columns[:])))[0]
	countsfinal=np.transpose(np.array(combined.as_matrix(columns=combined.columns[:])))[0]
	Efinal=-np.array(cumsumcombined.index)
	return [Efinal,cumsumfinal,countsfinal]


def getEfieldCounts(peaks):
	Earr=[]
	for p in peaks:
		Earr.append(p)
	if(len(Earr)==0):
		return [np.array([]),np.array([])]
	negE=np.sort(-np.array(Earr))
	cumsum=0
	cumsumarr = []
	countsatE = []
	uniqueEarr=[]
	thisEcount = 1
	prevE=0 
	numtimesexceedsprev = 0


	combined=pd.DataFrame({'negE':negE,'counts':np.ones(len(negE))}).groupby(by=["negE"]).sum()
	countsfinal=np.transpose(np.array(combined.as_matrix(columns=combined.columns[:])))[0]
	Efinal=-np.array(combined.index)
	return [Efinal,countsfinal]
