import pandas as pd
import geopandas as gpd
import geoplot as gplt
import matplotlib.pyplot as plt
from shapely.geometry import Point
import matplotlib.colors as colors

world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
world.to_crs(3857)
df = pd.read_pickle("Data/SmallData/globalEfields.pkl")
crs = {"init": "epsg:3857"}
geometry = [Point(xy) for xy in zip(df["longs"], df["lats"])]
df = gpd.GeoDataFrame(df, geometry=geometry)

ax = df.plot(
    column="E",
    legend=True,
    cmap="viridis",
    legend_kwds={
        "label": "Field Level (V/km)",
        "orientation": "horizontal",
    },
    norm=colors.LogNorm(vmin=df.E.min(), vmax=df.E.max()),
)
gplt.polyplot(world, ax=ax, zorder=1)
plt.title(
    "Magnitude of Peak "
    + str(60)
    + " Second Geoelectric Field\n1 in "
    + str(1 / 0.01)
    + " Year Storm"
)
plt.show()
