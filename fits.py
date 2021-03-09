#powerfit adapted from https://scipy-cookbook.readthedocs.io/items/FittingData.html
#lognormal fit adapted from https://machinelearningmastery.com/curve-fitting-with-python/
import numpy as np
from numpy import pi, r_
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.optimize import curve_fit
from scipy.stats import lognorm

def powerlaw(x,slope,exponent):
	return slope * (x**exponent)

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

def importedlognormal(x,upsilon,epsilon):
	s=np.abs(epsilon)
	scale=np.exp(upsilon)
	y=np.abs(x)/np.abs(scale)
	return lognorm.pdf(y, s)/np.abs(scale)

def locimportedlognormal(x,upsilon,epsilon,loc):
	s=np.abs(epsilon)
	scale=np.exp(upsilon)
	y=(np.abs(x-loc))/scale
	return lognorm.pdf(y, s)/scale

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

def fitLognormal(xdata,ydata):
	# plt.figure(45)
	# plt.subplot(2, 1, 2)
	# guessupsilon=.0005
	# guessepsilonsqd=50000
	params, covar = curve_fit(lognormal,xdata,ydata,p0=[-5.8,.37],maxfev=10000)
	print('covar')
	print(covar)
	print('params')
	print(params)
	upsilon, epsilonsqd = params
	# plt.figure()
	# plt.loglog(xdata, lognormal(np.array(xdata),upsilon,epsilonsqd))
	# plt.loglog(xdata, ydata)
	# plt.xlabel('X (log scale)')
	# plt.ylabel('Y (log scale)')
	# plt.show()

	return [upsilon,epsilonsqd]

def fitLognormalwithloc(xdata,ydata):
	# plt.figure(45)
	# plt.subplot(2, 1, 2)
	# guessupsilon=.0005
	# guessepsilonsqd=50000 
	try:
		params, covar = curve_fit(locimportedlognormal,xdata,ydata,\
			#p0=[5.9588582, 4.37993949, 0.00277132933],\
			p0=[-3.07289437, 0.02616835, 0.04285109],\
			#p0=[-6.43616259e-04, 4.75736705e-03, 9.89547338e-01],\
			maxfev=10000)
			# )
		upsilon,epsilon,loc=params
		fitdata=locimportedlognormal(np.array(xdata),upsilon,np.abs(epsilon),loc)
		# params, covar = curve_fit(locimportedlognormal,xdata,ydata,p0=[-5.8,.37,0.00000001],maxfev=10000)
		# print('covar')
		# print(covar)
		r = np.linalg.norm(ydata-fitdata)
		print('norm of residuals')
		print(r)
	except:
		print('fitfailed')
		params=[-5.8,.37,0.0001]
	# params, covar = curve_fit(lognormal,xdata,ydata,p0=[-5.8,.37],maxfev=10000)
	print('params')
	print(params)
	upsilon, epsilon, loc = params
	# plt.figure()
	# plt.loglog(xdata, locimportedlognormal(np.array(xdata),upsilon,epsilon,loc))
	# plt.loglog(xdata, ydata)
	# plt.xlabel('X (log scale)')
	# plt.ylabel('Y (log scale)')
	# plt.show()

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
	exponentguess=-1/1.10175936
	slopeguess=1/5

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
	logy = np.log10(ydata)

	# define our (line) fitting function
	fitfunc = lambda p, x: p[0] + p[1] * x
	errfunc = lambda p, x, y: (y - fitfunc(p, x))

	pinit = [slopeguess, exponentguess]
	out = optimize.leastsq(errfunc, pinit,args=(logx, logy), full_output=1)

	pfinal = out[0]
	covar = out[1]
	# print(pfinal)
	# print(covar)

	exponent = pfinal[1]
	slope = 10.0**pfinal[0]


	##########
	# Plotting data
	##########
	# plt.figure(45)
	# plt.clf()
	# plt.subplot(2, 1, 1)
	# plt.plot(xdata, ydata)     # Fit

	# plt.plot(xdata, powerlaw(xdata, slope, exponent))     # Fit
	# plt.title('Best Fit Power Law')
	# plt.xlabel('X')
	# plt.ylabel('Y')

	# plt.subplot(2, 1, 2)
	# plt.loglog(xdata, ydata)
	# plt.loglog(xdata, powerlaw(xdata, slope, exponent))
	# plt.xlabel('X (log scale)')
	# plt.ylabel('Y (log scale)')
	# plt.show()
	return [slope,exponent]
