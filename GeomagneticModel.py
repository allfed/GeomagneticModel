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
 
import Model.ValidateModel
from Model.ValidateModel import ValidateModel

import Model.Network
from Model.Network import Network

# import Model.GIC_Model
# from Model.GIC_Model import GIC_Model

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
		choices=['TFsite','MTsite','GCmodel','PowerGrid','EarthModel','ValidateModel'])

	arg_parser.add_argument(
		'Function', metavar='function', type=str,
		nargs='?',
		help='Specify subpart of geomagnetic model to run (a specific function as part of the model specified). If you want the entirety of this aspect of the model to be run, don\'t include any second argument. Default: entirety of model processed.',
		choices=['calcEfields','calcRecurrence','calcRecurrenceFit','calcCombinedRecurrence','calcGlobalModel','evalGCmodel','calcEvsDuration','calcTimeBetweenStorms','WorldNetwork','Region','CalculateMagDict'])


	arg_parser.add_argument(
		'-f', '--fit', metavar='fit', type=str,
		help='Specify fit method for historical MT magnetic data. Default: lognormal',
		default='all',
		choices=['lognormal','power','all'], required=False)

	arg_parser.add_argument(
		'-p', '--plots', metavar='plots', type=str,nargs='+',
		help='Which plots to show. Default: GlobalEfield TransformersInDanger WhereNotEnoughSpares',
		default=['GlobalEfield','TransformersInDanger','WhereNotEnoughSpares'],
		choices=['StormRecurrence','GlobalEfield','GlobalConductivity','TransformersInDanger','WhereNotEnoughSpares','1989Efields','Estats','ApparentResistivity','EvsDuration','AdjustedRates','CombinedRates','CompareGCandTF','MapOfSites'], required=False) 
	arg_parser.add_argument(
		'-r', '--rate-per-year', metavar='rateperyears', type=float, nargs='+',
		help='Specify rates per year for analysis. Example: 0.1 .01 .001 (once per decade, once per century and once per millenia). Default: .01',
		default=.01,
		required=False)

	arg_parser.add_argument(
		'-con', '--continent', metavar='continent', type=str, nargs='?',
		help='Continent for network analysis.',
		required=False)

	arg_parser.add_argument(
		'-cou', '--country', metavar='country', type=str, nargs='?',
		help='Continent for network analysis. Optional ',
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
	print('done importing')
	print('Detected arguments: '+str(args))
	
	#this processing is run every time, regardless of input args
	earthmodel=EarthModel()
	powergrid=PowerGrid()

	rateperyears=args['rate_per_year']
	if(type(rateperyears)==type(0.0)):
		rateperyears=[rateperyears]

	print('initializing models')
	# tfsites=earthmodel.initTFsites()
	mtsites=earthmodel.initMTsites()
	gcmodel=earthmodel.initGCmodel()
	earthmodel.calcGCcoords()

	if(args['Model']=='EarthModel' and args['Function']=='CalculateMagDict'):
		earthmodel.generateGeomagneticArray()
		# magdictgenerated=True
	earthmodel.loadgeotomagDict()

	if('MapOfSites' in args['plots']):
		earthmodel.plotSites(tfsites,mtsites)
		quit()
	if(args['Model']=='PowerGrid'): 
		if(args['Function']=='WorldNetwork'):
			# powergrid.createNetwork()
			# powergrid.createRegionNetwork('europe','')

			# powergrid.createRegionNetwork('europe','')
			# powergrid.createRegionNetwork('south-america','')
			powergrid.createRegionNetwork('africa','')
			# powergrid.createRegionNetwork('north-america','')
			# powergrid.createRegionNetwork('australia-oceania','')
			# powergrid.createRegionNetwork('central-america','')
			# powergrid.createRegionNetwork('russia','')
			# powergrid.createRegionNetwork('asia','')
			powergrid.plotNetwork()
			# powergrid.calcGICs()
			quit()
		if(args['Function']=='Region'):
			continent=args['continent']
			continents=['europe','south-america','africa','north-america','australia-oceania','central-america','russia','asia']
			if(len(continent)==0 or not (continent in continents)):
				print('Error: Continent required to process region. Options are:')
				print(continents)
				quit()
			country=args['country']
			print('creating region network:'+continent +' '+str(country) )			
			powergrid.createRegionNetwork(continent,country)
			print('2')
			earthmodel.loadCombinedFits()
			earthmodel.loadDurationRatios()
			earthmodel.loadApparentCond()
			network=powergrid.networks[0]
			#the network contains information on the boundary of the region and the name.
			network.EfieldFiles=earthmodel.calcRegionEfields(rateperyears,network,plot=True)#the first and only network, as we're only calculating one region
			# earthmodel.plotRegionEfields()
			print('3')
			network.calcGICs()
			# powergrid.setWindowedEmaps(earthmodel)
			print('4')
			powergrid.calcTransformerFailures(network,earthmodel.allwindowperiods,earthmodel.averagedRatios,rateperyears)
			#load E fields for the 

			# powergrid.createRegionNetwork('europe','estonia')
			# powergrid.createRegionNetwork('south-america','')
			# powergrid.createRegionNetwork('africa','')
			# powergrid.createRegionNetwork('north-america','')
			# powergrid.createRegionNetwork('australia-oceania','')
			# powergrid.createRegionNetwork('central-america','')
			# powergrid.createRegionNetwork('russia','')
			# powergrid.createRegionNetwork('asia','')
			# powergrid.plotNetwork()
			# powergrid.calcGICs()
			quit()
		earthmodel.loadCombinedFits()
		earthmodel.loadDurationRatios()
		powergrid.setWindowedEmaps(earthmodel)
		powergrid.loadEfieldMaps()
		# powergrid.plotOverheatByDuration()
		# powergrid.calcTemperatureMap()
		# powergrid.calcOverheatMap()
		print('powergrid.pop_est')
		powergrid.calcPopulationAffected()
		print('powergrid.eleestimate')
		powergrid.calcElectricityAffected()
		quit()

	recurrencecalculated = False
	gcmodelprocessed = False
	durationratiosprocessed = False


	if(args['Model']=='TFsite'): 
		
		for tfs in tfsites:
			#tfs.printApparentc()
			i=tfs.getClosestFreqIndex(8**-3)
			# print(tfs.name+' apparent conductivity')
			# print(tfs.apparentc[i])
			if('ApparentResistivity' in args['plots']):
				tfs.plotApparentr()
			# print(np.sqrt(.43814)/np.sqrt(tfs.apparentc[i]))
			# if('ApparentResistivity' in args['plots']):
			# 	tfs.plotApparentr()

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
		if(args['Function']=='calcTimeBetweenStorms'):
			earthmodel.loadStorms(mtsites)
			earthmodel.calcTimeBetweenStorms(mtsites)
		elif(not args['Function']):
			earthmodel.processMTsites(mtsites,tfsites)
			recurrencecalculated = True
	if(args['Model']=='GCmodel'): 
		gcmodelprocessed=True
		gcmodel.importGCmodel()
		gcmodel.getImpedanceMap()
		quit()


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

	calccombinedrates=False
	calcglobalmodel=False
	evalgcmodel=False
	calcrecurrence=False
	calcpeakevsduration=False
	calcEvsDuration=True
	# Run earthmodel to first adjust mtsites to a consistent reference ground conductivity and geomagnetic latitude
	if(args['Model']=='EarthModel'):
		if(not recurrencecalculated):
			earthmodel.loadPreviousMTfits(mtsites)
		print(' previousmtfits loaded')
		if(not args['Function']):#do all the tasks
			calccombinedrates=True
			calcEvsDuration=True
			calcglobalmodel=True
			calcrecurrence=True
			calcpeakevsduration=True
		if('calcEvsDuration' in str(args['Function'])):
			calcpeakevsduration=True
			calcEvsDuration=True
		if('calcCombinedRecurrence' in str(args['Function'])):
			calcpeakevsduration=True
			calcEvsDuration=True
			calccombinedrates=True
			calcrecurrence=True
		if('calcGlobalModel' in str(args['Function'])):
			calccombinedrates=True
			calcEvsDuration=True
			calcrecurrence=True
			calcpeakevsduration=True
			calcglobalmodel=True
		if('evalGCmodel' in str(args['Function'])):
			evalgcmodel=True

	validatemodel=False
	if(args['Model']=='ValidateModel'):
		if(not args['Function']):#calc global model, then validate
			validatemodel=True
	if(validatemodel):
		validation=ValidateModel()
		validation.calcGIC(1)
		validation.calcGIC(5)
		validation.compareFields(earthmodel)
		quit()
	allTFandGCcompared=False

	if(evalgcmodel):
		earthmodel.compareAllTFsites()
		allTFandGCcompared=True
	#if no MT site was processed, load and use data from MTsite0 modeloutput 	for the plotting

	peakEvsDurationprocessed=False
	if((not durationratiosprocessed) and calcpeakevsduration):
		earthmodel.peakEvsDuration(mtsites,False)
		peakEvsDurationprocessed=True
		print('peakEvsDuration calculated')

	plotadjusted=False
	if('EvsDuration' in args['plots']):
		calcEvsDuration=True
		plotadjusted=True

	plotcombined=False
	if('CombinedRates' in args['plots']):
		calcEvsDuration=True
		calccombinedrates=True
		plotcombined=True

	if(calcEvsDuration):
		if(not peakEvsDurationprocessed):
			earthmodel.peakEvsDuration(mtsites,plotadjusted)
			peakEvsDurationprocessed=True

	if(calccombinedrates):
		earthmodel.calcReferenceRateperyear(tfsites,mtsites,args['fit'],plotadjusted)
		earthmodel.calcCombinedRates(plotcombined)
		earthmodel.adjustEfieldsToMatch(plotcombined)
		earthmodel.calcMatchedCombinedRates(plotcombined)
		print('calcCombinedRates ran')


	if(calcglobalmodel):
		#apply the reference site back out across the earth by adjusting  for of reference ground conductivity, geomagnetic latitude, and determing the rate per year of each windowperiod across the earth. The output of this function is a series of maps with electric field levels at each duration for a given rate per year of interest.
		earthmodel.calcGlobalEfields(rateperyears)

	if(validatemodel):
		validation=ValidateModel()
		validation.calcGIC(1)
		validation.calcGIC(5)
		validation.compareFields(earthmodel)

	if('CompareGCandTF' in args['plots']):
		if(not allTFandGCcompared):
			earthmodel.loadGCtoTFcomparison()
		else:
			earthmodel.compareAllTFsites()
		earthmodel.plotGCtoTFcomparison()

	if(args['Model']=='PowerGrid'): 
		powergrid.setWindowedEmaps(earthmodel)
		powergrid.loadEfieldMaps()
		powergrid.calcOverheatMap()