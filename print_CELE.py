import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

# plt.rcParams.update({"font.size": 32})
sb.set_context("poster", font_scale=1.5)
# europe pop taken from https://www.worldometers.info/world-population/europe-population/
# usa from Data/SmallData/ElectricityByNation/PopByCountry.csv
# year=2018
x = pd.read_csv("CELE_region_values.csv")
# x = x[x["region"] != "usa"]
x["popRatio"] = (x["populationCELE"] / x["population"]) * 100
x["oneInXyears"] = 1 / x["rate"]
print(x)

_, axs = plt.subplots(1, 2, tight_layout=True, figsize=(22, 9))
for ax in axs:
    ax.set(xscale="log", yscale="linear")
sb.lineplot(
    x,
    x="oneInXyears",
    y="popRatio",
    hue="region",
    style="region",
    markers=True,
    dashes=False,
    markersize=28,
    linestyle="--",
    legend="full",
    ax=axs[0],
)
axs[0].set_xlabel("Strom recurrence [years]")
axs[0].set_ylabel("Population in CELE [%]")
sb.lineplot(
    x,
    x="oneInXyears",
    y="populationCELE",
    hue="region",
    style="region",
    markers=True,
    dashes=False,
    markersize=28,
    linestyle="--",
    legend=False,
    ax=axs[1],
)
axs[1].set_xlabel("Strom recurrence [years]")
axs[1].set_ylabel("Population in CELE [people]")
plt.savefig("CELE_values.eps", bbox_inches="tight", dpi=72)
# plt.show()
