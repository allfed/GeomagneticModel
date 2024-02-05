import fits
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score as sk_r2_score, mean_squared_error


def r2_score(y_true, y_pred):
    residuals = y_true - y_pred
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared


windowperiod = 60
combinedcounts = np.load("FitTests/FRN_combinedcounts.npy")
combinedrates = np.load("FitTests/FRN_combinedrates.npy")
probtoRPYratio = np.max(combinedrates)
combinedEfields = np.load("FitTests/FRN_combinedEfields.npy")

data = np.vstack((combinedEfields, combinedrates, combinedrates)).T
data = data[data[:, 0].argsort()]

[exponent] = fits.fitPower(data[:, 0], data[:, 1] / probtoRPYratio)
# self.combinedmatchedpowerfits = self.combinedmatchedpowerfits + [[windowperiod,linearfit[0],np.exp(linearfit[1])]]
# use PDF to determine mean and standard deviation of underlying normal distribution
[guessMean, guessStd] = fits.getGuesses(combinedEfields, combinedcounts, False)

# fit to the datapoints (CDF has a probability of 1 at the first datapoint)
[mean, std, loc] = fits.fitLognormalCDF(
    combinedEfields,
    combinedrates / probtoRPYratio,
    guessMean,
    guessStd,
    False,
)
fig, ax = plt.subplots()
ax.set_xscale("log")
ax.set_yscale("log")

plt.plot(
    data[:, 0],
    data[:, 1],
    ".",
    lw=1,
    label="Field averaged over " + str(windowperiod) + " seconds",
)

powerfit = fits.powerlaw(data[:, 0], exponent) * probtoRPYratio
print("r2 power: ", r2_score(data[:, 1], powerfit))
print("sk r2 power: ", sk_r2_score(data[:, 1], powerfit))
print("RMSE power: ", mean_squared_error(data[:, 1], powerfit, squared=False))

plt.plot(
    data[:, 0],
    powerfit,
    lw=1.5,
    label="Powerfit, field averaged over " + str(windowperiod) + " seconds",
)

logfit = fits.logcdf(data[:, 0], mean, np.abs(std), loc) * probtoRPYratio
print("r2 log: ", r2_score(data[:, 1], logfit))
print("sk r2 log: ", sk_r2_score(data[:, 1], logfit))
print("RMSE log: ", mean_squared_error(data[:, 1], logfit, squared=False))

plt.plot(
    data[:, 0],
    logfit,
    lw=1.5,
    label="Logfit, field averaged over " + str(windowperiod) + " seconds",
)

plt.legend()
plt.title("Rate geoelectric field is above threshold")
plt.xlabel("Geoelectric Field (V/km)")
plt.ylabel("Rate distinct storms above E field (storms/year)")

plt.show()
