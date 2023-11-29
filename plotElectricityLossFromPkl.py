import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import geoplot as gplt
from shapely.geometry import Polygon


def make_bbox(long0, lat0, long1, lat1):
    bbox = Polygon([[long0, lat0], [long1, lat0], [long1, lat1], [long0, lat1]])
    bbox_gdf = gpd.GeoDataFrame(index=[0], crs="epsg:4326", geometry=[bbox])
    return bbox_gdf


def plotElectricityLoss(worldElectricity, region, rate):
    _, ax = plt.subplots()

    if region == "europe":
        bbox = make_bbox(-36.386719, 29.228890, 60.292969, 74.543330)
        worldElectricity = worldElectricity.overlay(bbox, how="intersection")

    gplt.choropleth(
        worldElectricity,
        hue=worldElectricity["fraction"],
        cmap="viridis",
        ax=ax,
        legend=True,
        edgecolor="None",
    )

    plt.title(
        "Predicted Fraction Electricity Lost By Country \n One in "
        + str(1 / rate)
        + " Year Storm"
    )
    plt.savefig(
        "Data/SmallData/Figures/Europe_USA/"
        + region
        + str(rate)
        + "fraction_electr_lost.png",
        bbox_inches="tight",
    )
    plt.show()


file = "Data/SmallData/worldElectricity_europe_0.0001.pkl"
worldElectricity = pd.read_pickle(file)
rate = float(file[file.rfind("_") + 1 : file.rfind(".")])
region = file[file.find("_") + 1 : file.rfind("_")]
plotElectricityLoss(worldElectricity, region, rate)
