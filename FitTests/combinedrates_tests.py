import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import lognorm, powerlaw, pareto
from scipy.optimize import curve_fit, fsolve
from sklearn.metrics import r2_score


def lognorm_wrapper(x, s, loc, scale):
    return lognorm.pdf(x, s, loc, scale)


def powerlaw_wrapper(x, a, loc, scale):
    return pareto.pdf(x, a, loc, scale)


def read_data():
    prop_cycle = plt.rcParams["axes.prop_cycle"]
    colors = prop_cycle.by_key()["color"]

    # I'm assuming here that these are storms per year of strength E or higher
    # from magnetoalluric sites adjusted to the reference location
    series_length = [1311, 1409, 595, 316, 846, 103]
    x = np.load("FitTests/combinedEfields.npy")
    y = np.load("FitTests/combinedRates.npy")
    idx = 0
    for ii, sl in enumerate(series_length):
        a = x[idx : idx + sl]
        b = y[idx : idx + sl]
        b /= np.sum(b)
        plt.loglog(a, b, ".", color=colors[ii])
        # popt, pcov = curve_fit(lognorm_wrapper, a, b)
        # pl = lognorm_wrapper(a, *popt)
        # print("lognorm: ", r2_score(b, pl))
        # plt.loglog(a, pl, "-", color=colors[ii])
        # popt, pcov = curve_fit(powerlaw_wrapper, a, b)
        # pl = powerlaw_wrapper(a, *popt)
        # print("pareto: ", r2_score(b, pl))
        # plt.loglog(a, pl, "--", color=colors[ii])
        idx += sl
    plt.show()
    exit(3)
    data = np.vstack((x, y)).T
    data = data[data[:, 0].argsort()[::-1]]
    return data


data = read_data()
# popt, pcov = curve_fit(lognorm_wrapper, data[:, 0], data[:, 1])
# logfit = lognorm_wrapper(data[:, 0], *popt)

# y = 1 / 50
# z = fsolve(lambda x: lognorm_wrapper(x, *popt) - y, x0=25)
# print(z)


fig, ax = plt.subplots()
ax.set_yscale("log")
ax.set_xscale("log")

plt.plot(data[:, 0], data[:, 1], ".")
# plt.plot(data[:, 0], logfit, "-")
# plt.plot(z, y, "rx")
plt.show()
