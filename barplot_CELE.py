from operator import itemgetter
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb

# plt.rcParams.update({"font.size": 32})
sb.set_context("poster")
# europe pop taken from https://www.worldometers.info/world-population/europe-population/
# usa from Data/SmallData/ElectricityByNation/PopByCountry.csv
# year=2018
x = pd.read_csv("CELE_region_values.csv")
# x = x[x["region"] != "usa"]
x["popRatio"] = (x["populationCELE"] / x["population"]) * 100
x["oneInXyears"] = (1 / x["rate"]).astype(int)
x = x[x["popRatio"] > 0]
print(x)
palette = itemgetter(0, 3)(sb.color_palette())
_, axs = plt.subplots(1, 3, tight_layout=True, figsize=(22, 9))
sb.barplot(
    x,
    x="oneInXyears",
    y="popRatio",
    hue="region",
    legend=False,
    ax=axs[0],
    palette=palette,
)
axs[0].set_xlabel("Storm recurrence [years]")
axs[0].set_ylabel("Population in CELE [%]")
sb.barplot(
    x,
    x="oneInXyears",
    y="populationCELE",
    hue="region",
    legend=False,
    ax=axs[1],
    palette=palette,
)
axs[1].set_xlabel("Storm recurrence [years]")
axs[1].set_ylabel("Population in CELE [people]")
sb.barplot(
    x,
    x="oneInXyears",
    y="totalGWloss",
    hue="region",
    legend="full",
    ax=axs[2],
    palette=palette,
)
axs[2].set_xlabel("Storm recurrence [years]")
axs[2].set_ylabel("Electricity loss [GW]")
plt.savefig("CELE_bar.eps", bbox_inches="tight", dpi=72)
# plt.show()
