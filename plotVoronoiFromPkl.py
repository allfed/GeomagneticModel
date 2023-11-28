import matplotlib.pyplot as plt
import geopandas
import geoplot as gplt
import pandas as pd
import mapclassify as mc


def plotCombinedVoronoi(sRegions, rate, cutoff, continent=None):
    fig, ax = plt.subplots(figsize=(12, 10))

    sRegions["powerOut" + str(rate)] = sRegions[str(rate)] > 0.33

    if continent == "Europe":
        world = geopandas.read_file(geopandas.datasets.get_path("naturalearth_lowres"))
        world = world[world["name"] == "Russia"]
        inter = geopandas.sjoin(world, sRegions)
        sRegions = sRegions[~sRegions.index.isin(inter["index_right"])]

    gplt.choropleth(
        sRegions,
        hue=sRegions["powerOut" + str(rate)],
        cmap="viridis",
        ax=ax,
        legend=True,
        legend_labels=["", "outage"],
        edgecolor="None",
        scheme=mc.UserDefined(sRegions["powerOut" + str(rate)].values, [0, 1]),
    )
    plt.title("Substation at Risk of Electricity Loss, " + str(rate) + "per Year Storm")
    plt.show()
    # plt.savefig(Params.figuresDir + str(rate) + "peryearOutage.png")


sRegions = pd.read_pickle("Data/SmallData/combinedVoronoi_europe_0.0001.pkl")
plotCombinedVoronoi(sRegions, 0.0001, True, "Europe")
