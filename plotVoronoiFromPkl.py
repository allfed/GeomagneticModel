import matplotlib.pyplot as plt
import geopandas
import geoplot as gplt
import pandas as pd
import mapclassify as mc


def plotCombinedVoronoi(sRegions, region, rate):
    _, ax = plt.subplots(figsize=(12, 10))

    sRegions["powerOut" + str(rate)] = sRegions[str(rate)] > 0.33

    gplt.choropleth(
        sRegions,
        hue=sRegions["powerOut" + str(rate)],
        cmap="viridis",
        ax=ax,
        legend=False,
        # legend_labels=["", "outage"],
        edgecolor="None",
        scheme=mc.UserDefined(sRegions["powerOut" + str(rate)].values, [0, 1]),
    )
    plt.title(
        "Substation at Risk of Electricity Loss, " + str(rate) + " per Year Storm"
    )
    plt.savefig(
        "Data/SmallData/Figures/Europe_USA/" + region + str(rate) + "peryearOutage.png",
        bbox_inches="tight",
    )
    plt.show()


file = "Data/SmallData/combinedVoronoi_europe_0.0001.pkl"
sRegions = pd.read_pickle(file)
rate = float(file[file.rfind("_") + 1 : file.rfind(".")])
region = file[file.find("_") + 1 : file.rfind("_")]
plotCombinedVoronoi(sRegions, region, rate)
