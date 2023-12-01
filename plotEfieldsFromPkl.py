import pandas as pd
import geopandas as gpd
import geoplot as gplt
import matplotlib.pyplot as plt
from shapely.geometry import Point
import matplotlib.colors as colors
from shapely.geometry import Polygon


def make_bbox(long0, lat0, long1, lat1):
    bbox = Polygon([[long0, lat0], [long1, lat0], [long1, lat1], [long0, lat1]])
    bbox_gdf = gpd.GeoDataFrame(index=[0], crs="epsg:4326", geometry=[bbox])
    return bbox_gdf


def plotEfield(Efields_df, region):
    world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    # world = world.to_crs(3857)
    geometry = [Point(xy) for xy in zip(Efields_df["longs"], Efields_df["lats"])]
    Efields_df = gpd.GeoDataFrame(Efields_df, crs=4326, geometry=geometry)

    # http://bboxfinder.com/
    bbox = None
    if region == "europe":
        bbox = make_bbox(-36.386719, 29.228890, 60.292969, 74.543330)
    if region == "northamerica":
        bbox = make_bbox(-169.277344, 3.864255, -48.339844, 74.354828)
    if region == "usa":
        bbox = make_bbox(-129.902344, 23.725012, -62.753906, 50.513427)
    if bbox is not None:
        Efields_df = Efields_df.overlay(bbox, how="intersection")
        world = world.overlay(bbox, how="intersection")

    ax = Efields_df.plot(
        column="E",
        legend=True,
        cmap="viridis",
        legend_kwds={
            "label": "Field Level (V/km)",
            "orientation": "horizontal",
        },
        norm=colors.LogNorm(vmin=Efields_df.E.min(), vmax=Efields_df.E.max()),
    )
    gplt.polyplot(world, ax=ax, zorder=1)
    plt.title(
        "Magnitude of Peak "
        + str(60)
        + " Second Geoelectric Field\n1 in "
        + str(1 / 0.01)
        + " Year Storm"
    )
    # # this is broken for some reason, the borders are shifted and don't match the E fields
    # plt.savefig(
    #     "Data/SmallData/Figures/Europe_USA/" + region + "Efield.png",
    #     bbox_inches="tight",
    # )
    # need to show and save to file manually
    plt.show()


f = pd.read_pickle("Data/SmallData/globalEfields.pkl")
plotEfield(f, "europe")
plotEfield(f, "usa")
plotEfield(f, "northamerica")
