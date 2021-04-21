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

class Plotter:
	def __init__(self):
		pass

	#calculates the appropriate geometries for each cell in a geopandas grid for plotting purposes
	@staticmethod
	def calcGrid(df):
		latdiff=Params.latituderes#
		longdiff=Params.longituderes#
		xmin=np.min(df['lats'])
		ymin=np.min(df['longs'])
		xmax=np.max(df['lats'])
		ymax=np.max(df['longs'])
		cells=[]
		for index,row in df.iterrows():
			if(row['lats']==xmin or row['lats']==xmax or row['longs']==ymin or row['longs']==ymax):
				cells.append(shapely.geometry.box(0,0,0,0))
				continue

			cell=shapely.geometry.box(row['longs']-longdiff/2, row['lats']+latdiff/2,row['longs']+longdiff/2 ,row['lats']-latdiff/2)
			cells.append(cell)
		crs={'init':'epsg:3857'}

		geo_df=gpd.GeoDataFrame(df,crs=crs,geometry=cells)
		return geo_df

	@staticmethod
	def plotRegionEfields(df,polygon):
		minE=np.min(df['E'])
		maxE=np.max(df['E'])

		crs={'init':'epsg:3857'}
		geometry=[Point(xy) for xy in zip(df['longs'],df['lats'])]
		# geo_df=gpd.GeoDataFrame(df,crs=crs,geometry=geometry)
		geo_df=Plotter.calcGrid(df)
		ax=geo_df.plot(column='E',legend=True,cmap='viridis',legend_kwds={'label': 'Field Level (V/km)','orientation': "horizontal"})
		polyGeo_df=gpd.GeoDataFrame({'geometry':polygon})
		pp=gplt.polyplot(polyGeo_df,ax=ax,zorder=1)
		plt.show()

	@staticmethod
	def formatticklabels(minval,maxval,pp):
		print('formatting')
		colourbar = pp.get_figure().get_axes()[1]
		ticks_loc = colourbar.get_xticks().tolist()
		ticks_loc.insert(0,minval)
		ticks_loc.append(maxval)
		colourbar.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
		
		labels=[]
		for i in range(0,len(ticks_loc)):
			x=ticks_loc[i]
			exponent = int(np.log10(x))
			coeff = x/10**exponent
			if(i==0 or i==len(ticks_loc)):
				labels.append(r"${:2.1f} \times 10^{{ {:2d} }}$".format(coeff,exponent))					
			else:
				labels.append(r"${:2.1f} \times 10^{{ {:2d} }}$".format(coeff,exponent))
		colourbar.set_xticklabels(labels,rotation=33)

	def plotNetwork(voltages,lines):

		n = len(voltages)
		colornames = plt.cm.coolwarm(np.linspace(0.1,0.9,n))
		if(len(voltages)==1):
			mycolors= [ListedColormap(colornames[0])]	
		mycolors = [ ListedColormap(colornames[i]) for i in range(0,len(voltages))]   
		f, ax = plt.subplots()
		maxi=100000
		# if(len(voltages)==1):
		# 	print('AAAh')
		# 	i=0
		# 	v=voltages[i]
		# 	network=gpd.GeoDataFrame({'voltage': np.array(v),'geometry':np.array(lines[i])})
		# 	gdf=gpd.GeoDataFrame(geometry=lines[i])

		# 	gdf.plot(cmap=mycolors[i],ax=ax,legend=True, label=str(voltages[i])+' kV network')

		for i in range(0,len(voltages)):
			v=voltages[i]
			# print('v')
			# print(v)
			# print(lines[i])
			# print(len(lines[i]))
			if(len(lines[i])<2):
				maxi=i
				continue
				# lines[i]=[lines[i],lines[i]]

			network=gpd.GeoDataFrame({'voltage': np.array(v),'geometry':np.array(lines[i])})
			gdf=gpd.GeoDataFrame(geometry=lines[i])

			gdf.plot(cmap=mycolors[i],ax=ax,legend=True, label=str(voltages[i])+' kV network')
		ax.legend()
		leg=ax.get_legend()
		# if(len(voltages)==1):
		# 	leg.legendHandles[i].set_color(colornames[0])
		for i in range(0,len(voltages)):
			if(i==maxi):
				break
			leg.legendHandles[i].set_color(colornames[i])

		# plt.show()