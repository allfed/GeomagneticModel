import matplotlib.pyplot as plt
import geopandas
import geoplot as gplt
import pandas as pd
import mapclassify as mc


def plotCombinedVoronoi(sRegions, rate, cutoff):
    if cutoff:
        fig, ax = plt.subplots(figsize=(12, 10))
        world = geopandas.read_file(geopandas.datasets.get_path("naturalearth_lowres"))
        pp = gplt.polyplot(world, ax=ax, zorder=1)
        sRegions.geometry = sRegions["geometry"]
        sRegions["powerOut" + str(rate)] = 0

        for index, row in sRegions.iterrows():
            if row[str(rate)] > 0.33:
                sRegions["powerOut" + str(rate)].iloc[index] = 1
            else:
                sRegions["powerOut" + str(rate)].iloc[index] = 0

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
        # pp2=gplt.polyplot(df,ax=ax,zorder=1)
        plt.title(
            "Substation at Risk of Electricity Loss, " + str(rate) + "per Year Storm"
        )
        # df.plot(ax=ax)
        plt.show()
        # plt.savefig(Params.figuresDir + str(rate) + "peryearOutage.png")

    else:
        fig, ax = plt.subplots(figsize=(12, 10))
        world = geopandas.read_file(geopandas.datasets.get_path("naturalearth_lowres"))
        pp = gplt.polyplot(world, ax=ax, zorder=1)
        sRegions.geometry = sRegions["geometry"]
        gplt.choropleth(
            sRegions,
            hue=sRegions[str(rate)],
            cmap="viridis",
            ax=ax,
            legend=True,
            edgecolor="None",
        )
        # pp2=gplt.polyplot(df,ax=ax,zorder=1)
        # df.plot(ax=ax)
        plt.title(
            "Fraction of Transformers at Substations that Overheat, "
            + str(rate)
            + "per Year Storm"
        )

        plt.show()

        # plt.savefig(Params.figuresDir + str(rate) + "peryearOverheat.png")


sRegions = pd.read_pickle("Data/SmallData/combinedVoronoi_europe_0.0001.pkl")
plotCombinedVoronoi(sRegions, 0.0001, True)
