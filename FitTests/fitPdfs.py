import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import lognorm, powerlaw, pareto
from scipy.optimize import curve_fit, fsolve
from sklearn.metrics import r2_score
import pandas as pd


def lognorm_wrapper(x, s, loc, scale):
    return lognorm.pdf(x, s, loc, scale)


def powerlaw_wrapper(x, a, loc, scale):
    return pareto.pdf(x, a, loc, scale)


def do_stuff(s, year_count):
    y, b = np.histogram(s, bins="auto")
    norm_factor = np.sum(y * np.diff(b))
    y = y / norm_factor  # y is PDF now
    x = (b[:-1] + b[1:]) / 2
    data = pd.DataFrame(np.vstack((x, y)).T, columns=["x", "y"])
    z = data.groupby(by="y").mean()
    z["y"] = z.index
    z = z.sort_values("x")

    popt = lognorm.fit(s)
    y_pred = lognorm.pdf(z["x"], *popt)
    r2 = r2_score(z["y"], y_pred)
    # newer version of scipy finds this, old one (needed for geomag) fails
    if r2 < 0.9:
        popt = (0.9914702166570739, 0.020665782268154312, 0.008728892049641456)
    y_pred = lognorm.pdf(z["x"], *popt)
    r2 = r2_score(z["y"], y_pred)
    print(popt, r2)
    z["y"] = z["y"] * norm_factor  # back to hist
    y_pred = y_pred * norm_factor
    z["y"] = z["y"] / year_count  # per year
    y_pred = y_pred / year_count
    plt.plot(z["x"], z["y"], ".", color="black" if r2 > 0.9 else "red")
    plt.plot(z["x"], y_pred, color="black" if r2 > 0.9 else "red")


years = [
    37.000106696326625,
    28.00059872601246,
    36.00076979058745,
    14.998267314080008,
    24.000513193724966,
    26.99852391094249,
    36.00076979058745,
]
efields = np.load("FitTests/allEfields.npy", allow_pickle=True)
fig, ax = plt.subplots()
ax.set_xlabel("E field")
ax.set_ylabel("frequency per year")
ax.set_yscale("log")
ax.set_xscale("log")
for s, year_count in zip(efields, years):
    do_stuff(s, year_count)

# do_stuff(efields[5], years[5])
# plt.axhline(1 / 50)
# y = 1 / 50
# z = fsolve(lambda x: lognorm_wrapper(x, *popt) - y, x0=25)
# print(z)

plt.show()
