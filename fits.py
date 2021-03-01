#adapted from https://scipy-cookbook.readthedocs.io/items/FittingData.html
import numpy as np
from numpy import pi, r_
import matplotlib.pyplot as plt
from scipy import optimize


def powerlaw(x,slope,exponent):
	return slope * (x**exponent)

# y=s*(x^e)
# y/s=x^e
# (y/s)^(1/e)=x
def powerlawxfromy(y,slope,exponent):
	return (y/slope)**(1/exponent)

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
