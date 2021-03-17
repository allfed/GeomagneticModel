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
def powerlawxfromy(y,slope,exponent):
	return (y/slope)**(1/exponent)

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

def fitLognormalCDF(xdata,ydata,guessmu,guesssigma,plot):
	try:
		params, covar = curve_fit(logcdf,xdata,ydata,p0=[guessmu,guesssigma,0],maxfev=10000)
	except:
		print('log fit failed, using guess.')
		return [guessmu,guesssigma,0]
	# mean,std=params
	# loc=0
	mean,std,loc=params
	# print('mean')
	# print(mean)
	# print('std')
	# print(std)
	# print('loc')
	# print(loc)

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
	x0=np.min(logx)
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
	if(len(allcountsatE)==0):
		allcountsatE=np.ones(np.len(E))

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
