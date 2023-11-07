import os
import geoplot.crs as gcrs
from matplotlib.colors import ListedColormap
import numpy as np
import shapely
import geopandas
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import geoplot as gplt
from shapely.geometry import Point
import Params
from mpl_toolkits.axes_grid1 import make_axes_locatable
import mapclassify as mc
from geovoronoi.plotting import subplot_for_map, plot_voronoi_polys_with_points_in_area


class Plotter:
    def __init__(self):
        pass

    def plotGICsBubble(df, network, rate):
        world = geopandas.read_file(geopandas.datasets.get_path("naturalearth_lowres"))
        crs = {"init": "epsg:3857"}
        if network.country == "great-britain":
            shape = world[(world.name == "United Kingdom")]
            ax = gplt.webmap(shape, projection=gcrs.WebMercator())
        else:
            shape = gpd.GeoDataFrame(
                {"geometry": network.boundaryPolygon}, crs="epsg:3857"
            )
            ax = gplt.polyplot(shape, projection=gcrs.AlbersEqualArea())
        print(ax)
        geometry = [Point(xy) for xy in zip(df["longs"], df["lats"])]

        geo_df = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)
        print(np.min(geo_df["GIC"].values))
        print(np.max(geo_df["GIC"].values))
        print(np.min(geo_df["cutoffGIC"].values))
        print(np.max(geo_df["cutoffGIC"].values))
        scheme = mc.UserDefined(geo_df["GIC"], bins=[25, 50, 75, 100, 125, 150, 175])

        loadmap = False

        # f, ax = plt.subplots()

        # plt.show()
        # plt.show()
        # scheme = mc.EqualInterval(geo_df['GIC'], k=10)
        gplt.pointplot(
            geo_df,
            projection=gcrs.AlbersEqualArea(),
            scale="cutoffGIC",
            limits=(1, 30 * (np.max(geo_df["cutoffGIC"].values) / 150)),
            alpha=0.3,
            edgecolor="black",
            extent=shape.total_bounds,
            cmap="jet",
            ax=ax,
            # cmap='rainbow',
            legend=True,
            legend_var="hue",
            hue="GIC",
            scheme=scheme,
        )

        # print('plotNetwork')
        # Plotter.plotNetwork(network.voltages,network.lines,ax)
        plt.title("GIC in " + network.region + " for " + str(rate) + " Per Year Storm")

        if not os.path.isdir(Params.figuresDir + "Regions/" + network.region):
            os.mkdir(Params.figuresDir + "Regions/" + network.region)

        plt.savefig(
            Params.figuresDir
            + "Regions/"
            + network.region
            + "/"
            + str(rate)
            + "peryearGIC.png"
        )
        plt.show()

    def plotNetwork(voltages, lines, ax, region):
        n = len(voltages)
        colornames = plt.cm.coolwarm(np.linspace(0.1, 0.9, n))
        if len(voltages) == 1:
            mycolors = [ListedColormap(colornames[0])]
        mycolors = [ListedColormap(colornames[i]) for i in range(0, len(voltages))]
        # f, ax = plt.subplots()
        maxi = 100000
        # if(len(voltages)==1):
        # 	print('AAAh')
        # 	i=0
        # 	v=voltages[i]
        # 	network=gpd.GeoDataFrame({'voltage': np.array(v),'geometry':np.array(lines[i])})
        # 	gdf=gpd.GeoDataFrame(geometry=lines[i])

        # 	gdf.plot(cmap=mycolors[i],ax=ax,legend=True, label=str(voltages[i])+' kV network')

        for i in range(0, len(voltages)):
            v = voltages[i]
            print("v")
            print(v)
            # print(lines[i])
            print(len(lines[i]))
            if len(lines[i]) < 2:
                maxi = i
                continue
                # lines[i]=[lines[i],lines[i]]

            network = gpd.GeoDataFrame(
                {"voltage": np.array(v), "geometry": np.array(lines[i])}
            )
            gdf = gpd.GeoDataFrame(geometry=lines[i])

            gdf.plot(
                cmap=mycolors[i],
                ax=ax,
                legend=True,
                label=str(voltages[i]) + " kV network",
            )
        ax.legend()
        leg = ax.get_legend()
        for i in range(0, len(voltages)):
            if i == maxi:
                break
            leg.legendHandles[i].set_color(colornames[i])

        plt.title("High Voltage network in " + region + " above 100kV")

        if not os.path.isdir(Params.figuresDir + "Regions/" + region):
            os.mkdir(Params.figuresDir + "Regions/" + region)

        plt.savefig(Params.figuresDir + "Regions/" + region + "/Network")
        return ax

    def plotRegionEfields(df, polygon, network):
        minE = np.min(df["E"])
        maxE = np.max(df["E"])

        crs = {"init": "epsg:3857"}
        geometry = [Point(xy) for xy in zip(df["longs"], df["lats"])]
        # geo_df=gpd.GeoDataFrame(df,crs=crs,geometry=geometry)
        geo_df = Plotter.calcGridStandard(df)
        polyGeo_df = gpd.GeoDataFrame({"geometry": polygon})
        pp = gplt.polyplot(polyGeo_df, ax=ax, zorder=1)
        world = geopandas.read_file(geopandas.datasets.get_path("naturalearth_lowres"))

        if network.country == "great-britain":
            shape = world[(world.name == "United Kingdom")]
        else:
            shape = gpd.GeoDataFrame(
                {"geometry": network.boundaryPolygon}, crs="epsg:3857"
            )
        ax = geo_df.plot(
            column="E",
            legend=True,
            cmap="viridis",
            legend_kwds={"label": "Field Level (V/km)", "orientation": "horizontal"},
        )
        pp = gplt.polyplot(
            shape,
            zorder=1,
            projection=gcrs.AlbersEqualArea(),
            extent=shape.total_bounds,
        )
        # Plotter.plotNetwork(network.voltages,network.lines,ax)
        plt.show()

    def plotRegionEfields(df, rate, network):
        minE = np.min(df["E"])
        maxE = np.max(df["E"])

        crs = {"init": "epsg:3857"}
        geometry = [Point(xy) for xy in zip(df["longs"], df["lats"])]
        # geo_df=gpd.GeoDataFrame(df,crs=crs,geometry=geometry)

        geo_df = Plotter.calcGridStandard(df)
        ax = geo_df.plot(
            column="E",
            legend=True,
            cmap="viridis",
            legend_kwds={"label": "Field Level (V/km)", "orientation": "horizontal"},
        )
        geo_df.geometry = geometry
        [newnodes, worldpolys] = network.getPolysWithNodes(geo_df)
        polyGeo_df = gpd.GeoDataFrame({"geometry": network.boundaryPolygon})

        pp = gplt.polyplot(polyGeo_df, ax=ax, zorder=1)
        world = geopandas.read_file(geopandas.datasets.get_path("naturalearth_lowres"))

        # gb_shape=world[(world.name == "United Kingdom")]
        # pp=gplt.polyplot(gb_shape,ax=ax,zorder=1)
        # Plotter.plotNetwork(network.voltages,network.lines,ax)
        if not os.path.isdir(Params.figuresDir + "Regions/" + network.region):
            os.mkdir(Params.figuresDir + "Regions/" + network.region)
        plt.title(
            "E Field in " + network.region + " for " + str(rate) + " Per Year Storm"
        )

        if not os.path.isdir(Params.figuresDir + "Regions/" + network.region):
            os.mkdir(Params.figuresDir + "Regions/" + network.region)

        plt.savefig(
            Params.figuresDir
            + "Regions/"
            + network.region
            + "/"
            + str(rate)
            + "peryearE.png"
        )
        plt.show()

    def calcGridStandard(df):
        latdiff = Params.latituderes  #
        longdiff = Params.longituderes  #
        xmin = np.min(df["lats"])
        ymin = np.min(df["longs"])
        xmax = np.max(df["lats"])
        ymax = np.max(df["longs"])
        cells = []
        for index, row in df.iterrows():
            if (
                row["lats"] == xmin
                or row["lats"] == xmax
                or row["longs"] == ymin
                or row["longs"] == ymax
            ):
                cells.append(shapely.geometry.box(0, 0, 0, 0))
                continue

            cell = shapely.geometry.box(
                row["longs"] - longdiff / 2,
                row["lats"] + latdiff / 2,
                row["longs"] + longdiff / 2,
                row["lats"] - latdiff / 2,
            )
            cells.append(cell)
        crs = {"init": "epsg:3857"}

        geo_df = gpd.GeoDataFrame(df, crs=crs, geometry=cells)
        return geo_df

    def calcGrid(df, lats, longs, latdiff, longdiff):
        xmin = np.min(df["lats"])
        ymin = np.min(df["longs"])
        xmax = np.max(df["lats"])
        ymax = np.max(df["longs"])
        cells = []
        for index, row in df.iterrows():
            if (
                row["lats"] == xmin
                or row["lats"] == xmax
                or row["longs"] == ymin
                or row["longs"] == ymax
            ):
                cells.append(shapely.geometry.box(0, 0, 0, 0))
                continue

            cell = shapely.geometry.box(
                row["longs"] - longdiff / 2,
                row["lats"] + latdiff / 2,
                row["longs"] + longdiff / 2,
                row["lats"] - latdiff / 2,
            )
            cells.append(cell)
        crs = {"init": "epsg:3857"}

        geo_df = gpd.GeoDataFrame(df, crs=crs, geometry=cells)
        return geo_df

    # def plotVoronoiRegions(failureProbs,geometry,boundary_shape,poly_shapes,regionPointsDict,region,rate):
    # 	fig, ax = subplot_for_map()

    # 	colors=plt.cm.viridis(failureProbs)#/np.max(nodes['failureProbs']))

    # 	# print('np.max(nodes[failureProbs])')
    # 	# print(np.max(nodes['failureProbs']))
    # 	# print('falprobs')
    # 	# print(nodes['failureProbs'])
    # 	# print(len(colors))
    # 	# print(len(poly_shapes))
    # 	# print('colors')
    # 	# print(colors)
    # 	# np.save( 'boundary_shape',boundary_shape,allow_pickle=True)
    # 	# np.save( 'poly_shapes',poly_shapes,allow_pickle=True)
    # 	# np.save( 'nodes',nodes,allow_pickle=True)
    # 	# np.save( 'colors',colors,allow_pickle=True)
    # 	# nodes.to_pickle("nodes.pkl")
    # 	# output = pd. read_pickle("a_file.plkl")

    # 	plot_voronoi_polys_with_points_in_area(ax, boundary_shape,poly_shapes,geometry,voronoi_and_points_cmap='viridis',voronoi_color=colors)
    # 	ax.set_title('Voronoi regions of Substations in '+region+' rate:'+str(rate))
    # 	plt.tight_layout()
    # 	plt.show()

    def plotCombinedVoronoi(sRegions, rate, cutoff):
        if cutoff:
            fig, ax = plt.subplots(figsize=(12, 10))
            world = geopandas.read_file(
                geopandas.datasets.get_path("naturalearth_lowres")
            )
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
                legend_labels=[0, 1],
                edgecolor="None",
                scheme=mc.UserDefined(sRegions["powerOut" + str(rate)].values, [0, 1]),
            )
            # pp2=gplt.polyplot(df,ax=ax,zorder=1)
            plt.title(
                "Substation at Risk of Electricity Loss, "
                + str(rate)
                + "per Year Storm"
            )
            # df.plot(ax=ax)
            # plt.show()
            plt.savefig(Params.figuresDir + str(rate) + "peryearOutage.png")

        else:
            fig, ax = plt.subplots(figsize=(12, 10))
            world = geopandas.read_file(
                geopandas.datasets.get_path("naturalearth_lowres")
            )
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

            # plt.show()

            plt.savefig(Params.figuresDir + str(rate) + "peryearOverheat.png")

    # <AxesSubplot:>
    # pp1=gplt.polyplot(polygdf,ax=ax,zorder=2)
