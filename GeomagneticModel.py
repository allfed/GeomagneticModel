# Copyright © Morgan Rivers, 2021 Alliance to Feed the Earth in Disasters 
# Email: morgan[at]allfed[dot]info>
# All Rights Reserved.
# This file is part of the GeomagneticModel package.
#
# The GeomagneticModel is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The GeomagneticModel is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the GeomagneticModel (in the LICENSE file).
import argparse
import numpy as np
import Model.PowerGrid
from Model.PowerGrid import PowerGrid

import Model.EarthModel
from Model.EarthModel import EarthModel

def get_args():
	arg_parser = argparse.ArgumentParser(
		description='''GeomagneticModel - assess dangers of geomagnetic storms''',
		usage='''GeomagneticModel.py [model] [-h] 
		[-f fit]
		[-p plots [plots ...]]
		[-m mtsites [mtsites ...]]
		[-t tfsites [tfsites ...]]
		[-r rateperyears [rateperyears ...]]''',
		formatter_class=argparse.RawTextHelpFormatter)

	arg_parser.add_argument(
		'Model', metavar='model', type=str,
		nargs='?',
		help='Specify part of geomagnetic model to run. If you do not require any of the model to be rerun, do not include this argument. WARNING: without a model argument included, ModelParams.ods is not updated in the model! Default: no model used.',
		choices=['TFsite','MTsite','GCmodel','PowerGrid','EarthModel'])

	arg_parser.add_argument(
		'Function', metavar='function', type=str,
		nargs='?',
		help='Specify subpart of geomagnetic model to run (a specific function as part of the model specified). If you want the entirety of this aspect of the model to be run, don\'t include any second argument. Default: entirety of model processed.',
		choices=['calcEfields','calcRecurrence','calcRecurrenceFit'])


	arg_parser.add_argument(
		'-f', '--fit', metavar='fit', type=str,
		help='Specify fit method for historical MT magnetic data. Default: lognormal',
		default='all',
		choices=['lognormal','power','all'], required=False)

	arg_parser.add_argument(
		'-p', '--plots', metavar='plots', type=str,nargs='+',
		help='Which plots to show. Default: GlobalEfield TransformersInDanger WhereNotEnoughSpares',
		default=['GlobalEfield','TransformersInDanger','WhereNotEnoughSpares'],
		choices=['StormRecurrence','GlobalEfield','GlobalConductivity','TransformersInDanger','WhereNotEnoumtsitesghSpares','1989Efields','Estats','ApparentResistivity','EvsDuration'], required=False) 
	arg_parser.add_argument(
		'-r', '--rate-per-year', metavar='rateperyears', type=float, nargs='+',
		help='Specify rates per year for analysis. Example: 0.1 .01 .001 (once per decade, once per century and once per millenia). Default: .01',
		default=.01,
		required=False)

	arg_parser.add_argument(
		'-t', '--TF-sites', metavar='tfsites', type=str, nargs='+',
		help='EMTF site names. Must match number of MT sites specified. Default: none (imported from ModelParams.ods)',
		required=False)

	arg_parser.add_argument(
		'-m', '--MT-sites', metavar='mtsites', type=str, nargs='+',
		help='Magnetotelluric site names. Must match number of FT sites specified. Default: none (imported from ModelParams.ods)',
		required=False)

	return vars(arg_parser.parse_args())


if __name__ == '__main__':
	args = get_args()
	print('Detected arguments: '+str(args))
	
	#this processing is run every time, regardless of input args
	earthmodel=EarthModel()
	powergrid=PowerGrid()

	tfsites=earthmodel.initTFsites()
	mtsites=earthmodel.initMTsites()
	gcmodel=earthmodel.initGCmodel()

	recurrencecalculated = False
	gcmodelprocessed = False
	durationratiosprocessed = False

	if(args['Model']=='TFsite'): 
		for tfs in tfsites:
			#tfs.printApparentc()
			i=tfs.getClosestFreqIndex(8**-3)
			print(tfs.name+' ratio to VAQ55')
			print(np.sqrt(.43814)/np.sqrt(tfs.apparentc[i]))
			if('ApparentResistivity' in args['plots']):
				tfs.plotApparentr()
	if(args['Model']=='MTsite'): 
		if(args['TF_sites'] or args['MT_sites']):
			if(args['TF_sites'] and args['MT_sites']):
				if(len(args['TF_sites']) != len(args['MT_sites'])):
					print("ERROR: when evaluating historic field data, must supply one TF site for each MT site. Modify your args to make sure they are equal in number.")
					quit()					
			else:
				print("ERROR: Must supply both MT and TF if MTsites are processed with MT sites and TF sites as command line arguments (as opposed to being auto-imported from ModelParams.ods).")
				quit()

		if(args['Function']=='calcEfields'):
			earthmodel.calcAndSaveEfields(mtsites,tfsites) 
		if(args['Function']=='calcRecurrence'):
			earthmodel.calcAndSaveRecurrence(mtsites,tfsites) 
			recurrencecalculated = True
		elif(not args['Function']):
			earthmodel.processMTsites(mtsites,tfsites)
			recurrencecalculated = True
	if(args['Model']=='GCmodel'): 
		gcmodelprocessed=True


	if('StormRecurrence' in args['plots']):
		earthmodel.calcandplotEratesPerYear(mtsites,args['fit'])
		# earthmodel.loadPreviousMTfits(mtsites)
		# mtsites[0].plotandFitEratesPerYear()
	if('Estats' in args['plots']):
		earthmodel.Estats(mtsites)
	if('1989Efields' in args['plots']):
		earthmodel.plot1989Efields(mtsites)
	if('EvsDuration' in args['plots']):
		earthmodel.peakEvsDuration(mtsites,True)
		durationratiosprocessed=True

	earthmodel.loadApparentCond()
	# Run earthmodel to first adjust mtsites to a consistent reference ground conductivity and geomagnetic latitude
	if(args['Model']=='EarthModel'): 
		#if no MT site was processed, load and use data from MTsite0 modeloutput 	for the plotting
		if(not recurrencecalculated):
			earthmodel.loadPreviousMTfits(mtsites)
		if(not durationratiosprocessed):
			earthmodel.peakEvsDuration(mtsites,False)

		earthmodel.calcReferenceRateperyear(tfsites,mtsites)
		earthmodel.plotCombinedRates()

		#apply the reference site back out across the earth by adjusting  for of reference ground conductivity, geomagnetic latitude, and determing the rate per year of each windowperiod across the earth. The output of this function is a series of maps with electric field levels at each duration for a given rate per year of interest.
		rateperyear=args['rate_per_year']
		if(type(rateperyear)==type(0.0)):
			earthmodel.calcGlobalEfields([rateperyear])
		else:
			earthmodel.calcGlobalEfields(rateperyear)

	if(args['Model']=='PowerGrid'): 
		powergrid=PowerGrid()
		# powergrid.setWindowedEmaps(earthmodel)
		powergrid.loadEfieldMaps()
		powergrid.calcOverheatMap()
