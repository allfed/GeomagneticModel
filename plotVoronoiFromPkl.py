import matplotlib.pyplot as plt
import geoplot as gplt
import pandas as pd
import mapclassify as mc
from shapely.geometry import Polygon
import geopandas
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
import numpy as np


def make_bbox(long0, lat0, long1, lat1):
    bbox = Polygon([[long0, lat0], [long1, lat0], [long1, lat1], [long0, lat1]])
    bbox_gdf = geopandas.GeoDataFrame(index=[0], crs="epsg:4326", geometry=[bbox])
    return bbox_gdf


def plotCombinedVoronoi(sRegions, region, rate, cutoff=True):
    _, ax = plt.subplots(figsize=(18, 22))

    if cutoff:
        sRegions["powerOut" + str(rate)] = sRegions[str(rate)] > 0.33
    else:
        sRegions["powerOut" + str(rate)] = sRegions[str(rate)]

    if region == "usa":
        # http://bboxfinder.com/
        bbox = make_bbox(-129.682617, 24.527135, -64.863281, 49.866317)
        sRegions = sRegions.overlay(bbox, how="intersection")

    no_cutoff_cmap = LinearSegmentedColormap.from_list("gr", ["grey", "red"], N=128)

    gplt.choropleth(
        sRegions,
        hue=sRegions["powerOut" + str(rate)],
        cmap="Set1_r" if cutoff else no_cutoff_cmap,
        ax=ax,
        legend=not cutoff,
        norm=PowerNorm(gamma=0.5, vmin=0, vmax=1),
        legend_kwargs={"shrink": 0.35, "label": "HVT failure rate"}
        if not cutoff
        else None,
        edgecolor="None",
        scheme=mc.UserDefined(sRegions["powerOut" + str(rate)].values, [0, 1])
        if cutoff
        else None,
    )
    no_cutoff_str = "" if cutoff else "no_cutoff"
    plt.savefig(
        "Data/SmallData/Figures/Europe_USA/"
        + region
        + str(rate)
        + no_cutoff_str
        + "peryearOutage.eps",
        bbox_inches="tight",
        dpi=72,
    )


plt.rcParams.update({"font.size": 24})

# file = "Data/SmallData/combinedVoronoi_europe_0.0001.pkl"
# file = "Data/SmallData/combinedVoronoi_usa_0.0001.pkl"
file = "Data/SmallData/combinedVoronoi_northamerica_0.0001.pkl"

sRegions = pd.read_pickle(file)
rate = float(file[file.rfind("_") + 1 : file.rfind(".")])
region = file[file.find("_") + 1 : file.rfind("_")]
plotCombinedVoronoi(sRegions, region, rate, False)
