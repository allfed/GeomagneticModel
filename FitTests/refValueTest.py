import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import lognorm
from scipy.optimize import curve_fit, fsolve


def lognorm_wrapper(x, s, loc, scale):
    return lognorm.pdf(x, s, loc, scale)


def read_data():
    x = np.load("FitTests/combinedMatchedEfields.npy")
    y = np.load("FitTests/combinedMatchedRates.npy")
    data = np.vstack((x, y)).T
    data = data[data[:, 0].argsort()[::-1]]
    return data


data = read_data()
popt, pcov = curve_fit(lognorm_wrapper, data[:, 0], data[:, 1])
logfit = lognorm_wrapper(data[:, 0], *popt)

y = 1 / 50
z = fsolve(lambda x: lognorm_wrapper(x, *popt) - y, x0=25)
print(z)


fig, ax = plt.subplots()
ax.set_yscale("log")
ax.set_xscale("log")

plt.plot(data[:, 0], data[:, 1], ".")
plt.plot(data[:, 0], logfit, "-")
plt.plot(z, y, "ro", markersize=10)
plt.show()
