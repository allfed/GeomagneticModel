import matplotlib.pyplot as plt
import geoplot as gplt
import pandas as pd
import mapclassify as mc
from shapely.geometry import Polygon
import geopandas


def make_bbox(long0, lat0, long1, lat1):
    bbox = Polygon([[long0, lat0], [long1, lat0], [long1, lat1], [long0, lat1]])
    bbox_gdf = geopandas.GeoDataFrame(index=[0], crs="epsg:4326", geometry=[bbox])
    return bbox_gdf


def plotCombinedVoronoi(sRegions, region, rate):
    _, ax = plt.subplots(figsize=(12, 10))

    sRegions["powerOut" + str(rate)] = sRegions[str(rate)] > 0.33

    if region == "usa":
        # http://bboxfinder.com/
        bbox = make_bbox(-129.682617, 24.527135, -64.863281, 49.866317)
        sRegions = sRegions.overlay(bbox, how="intersection")

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


# file = "Data/SmallData/combinedVoronoi_europe_0.0001.pkl"
# file = "Data/SmallData/combinedVoronoi_usa_0.0001.pkl"
file = "Data/SmallData/combinedVoronoi_northamerica_0.0001.pkl"

sRegions = pd.read_pickle(file)
rate = float(file[file.rfind("_") + 1 : file.rfind(".")])
region = file[file.find("_") + 1 : file.rfind("_")]
plotCombinedVoronoi(sRegions, region, rate)
