
def importIfNotAlready():
	global paramsinitialized
	try:
		paramsinitialized
	except NameError: #if the params haven't been initialized yet
		paramsinitialized = True
		importModelParams()
		importTransformers()
	else:
		return()


def importModelParams():
	from pyexcel_ods import get_data
	import numpy as np

	global paramsfilename
	global useMTsite
	global transformersfn
	global mtsitesdir
	global tfsitesdir
	global allTFsitesdir
	global mtRepeatRatesDir
	global mtRepeatRatePlotsDir
	global maxchunksize
	global apparentcondloc
	global globalcondloc
	global latituderes
	global longituderes
	global refconductivity
	global refmaglat
	global globalEfieldPlots
	global globalEfieldData
	global mtEfieldsloc
	global overheatMapsDir

	global tfsitenames
	global tfformats
	global mtsitenames
	global nchunks
	global sampleperiod
	global betaThreshold
	global toKeepThreshold
	global windows
	global window1
	global window2
	global window3
	global window4
	global window5
	global window6
	global window7
	global window8
	global window9
	global window10
	global window11
	global window12
	global window13
	global window14
	global window15
	global window16
	global window17
	global window18
	global window19
	global window20
	global window21
	global window22
	global window23
	global window24
	global window25

	paramsfilename = 'ModelParams.ods'
	data = get_data(paramsfilename)
	paramdata = data['Single Valued Params']

	column = []
	for coltitleindex in range(0,len(paramdata[1])):
		coltitle=paramdata[1][coltitleindex]
		if(coltitle == 'transformersfn'):
			transformersfn=paramdata[2][coltitleindex]
		if(coltitle == 'mtsitesdir'):
			mtsitesdir=paramdata[2][coltitleindex]
		if(coltitle == 'tfsitesdir'):
			tfsitesdir=paramdata[2][coltitleindex]
		if(coltitle == 'allTFsitesdir'):
			allTFsitesdir=paramdata[2][coltitleindex]
		if(coltitle == 'mtRepeatRatesDir'):
			mtRepeatRatesDir = paramdata[2][coltitleindex]
		if(coltitle == 'mtRepeatRatePlotsDir'):
			mtRepeatRatePlotsDir = paramdata[2][coltitleindex]
		if(coltitle == 'maxchunksize'):
			maxchunksize = paramdata[2][coltitleindex]
		if(coltitle == 'globalcondloc'):
			globalcondloc = paramdata[2][coltitleindex]
		if(coltitle == 'refconductivity'):
			refconductivity = paramdata[2][coltitleindex]
		if(coltitle == 'refmaglat'):
			refmaglat = paramdata[2][coltitleindex]
		if(coltitle == 'apparentcondloc'):
			apparentcondloc = paramdata[2][coltitleindex]
		if(coltitle == 'latituderes'):
			latituderes = paramdata[2][coltitleindex]
		if(coltitle == 'longituderes'):
			longituderes = paramdata[2][coltitleindex]
		if(coltitle == 'globalEfieldData'):
			globalEfieldData = paramdata[2][coltitleindex]
		if(coltitle == 'globalEfieldPlots'):
			globalEfieldPlots = paramdata[2][coltitleindex]
		if(coltitle == 'overheatMapsDir'):
			overheatMapsDir = paramdata[2][coltitleindex]
		if(coltitle == 'mtEfieldsloc'):
			mtEfieldsloc = paramdata[2][coltitleindex]
	
	paramdata = data['Params for each TF and MT site']

	column = []
	for coltitleindex in range(0,len(paramdata[1])):
		coltitle=paramdata[1][coltitleindex]
		lastvariableindex=len(paramdata[1])-1
		tmp = []
		if(coltitle): #if exists a column title
			for t in paramdata[2:]:
				if(not t):
					break 
				if(coltitleindex>len(t)-1):
					tmp = tmp + ['']
				else:
					tmp = tmp + [t[coltitleindex]]
		if(coltitle == 'useMTsite'):
			useMTsite=tmp
		if(coltitle == 'tfsitenames'):
			tfsitenames = tmp
		if(coltitle == 'tfformats'):
			tfformats = tmp
		if(coltitle == 'mtsitenames'):
			mtsitenames = tmp
		if(coltitle == 'nchunks'):
			nchunks = tmp
		if(coltitle == 'sampleperiod'):
			sampleperiod = tmp
		if(coltitle == 'betaThreshold'):
			betaThreshold = tmp
		if(coltitle == 'toKeepThreshold'):
			toKeepThreshold = tmp			
		if(coltitle == 'window1'):
			window1 = tmp
		if(coltitle == 'window2'):
			window2 = tmp
		if(coltitle == 'window3'):
			window3 = tmp
		if(coltitle == 'window4'):
			window4 = tmp
		if(coltitle == 'window5'):
			window5 = tmp
		if(coltitle == 'window6'):
			window6 = tmp
		if(coltitle == 'window7'):
			window7 = tmp
		if(coltitle == 'window8'):
			window8 = tmp
		if(coltitle == 'window9'):
			window9 = tmp
		if(coltitle == 'window10'):
			window10 = tmp
		if(coltitle == 'window11'):
			window11 = tmp
		if(coltitle == 'window12'):
			window12 = tmp
		if(coltitle == 'window13'):
			window13 = tmp
		if(coltitle == 'window14'):
			window14 = tmp
		if(coltitle == 'window15'):
			window15 = tmp
		if(coltitle == 'window16'):
			window16 = tmp
		if(coltitle == 'window17'):
			window17 = tmp
		if(coltitle == 'window18'):
			window18 = tmp
		if(coltitle == 'window19'):
			window19 = tmp
		if(coltitle == 'window20'):
			window20 = tmp
		if(coltitle == 'window21'):
			window21 = tmp
		if(coltitle == 'window22'):
			window22 = tmp
		if(coltitle == 'window23'):
			window23 = tmp
		if(coltitle == 'window24'):
			window24 = tmp
		if(coltitle == 'window25'):
			window25 = tmp
	windows=np.transpose([window1,window2,window3,window4,window5,window6,window7,window8,window9,window10,window11,window12,window13,window14,window15,window16,window17,window18,window19,window20,window21,window22,window23,window24,window25])

def importTransformers():
	from pyexcel_ods import get_data

	global HVTnames
	global temprise0
	global temprise10
	global temprise20
	global temprise40
	global temprise50
	global temprise100
	global temprise200
	global GICperE
	global pop25to40
	global pop0to25
	global pop40plus
	global tau

	data = get_data(transformersfn)
	thermaldata = data['Thermal']
	column = []
	for coltitleindex in range(0,len(thermaldata[1])):
		coltitle=thermaldata[1][coltitleindex]
		tmp = []
		if(coltitle): #if exists a column title
			for t in thermaldata[2:]:
				if(not t):
					break
				tmp = tmp + [t[coltitleindex]]
		if(coltitle == 'name'):
			HVTnames = tmp
		if(coltitle == 'temprise0'):
			temprise0 = tmp
		if(coltitle == 'temprise10'):
			temprise10 = tmp
		if(coltitle == 'temprise20'):
			temprise20 = tmp
		if(coltitle == 'temprise40'):
			temprise40 = tmp
		if(coltitle == 'temprise50'):
			temprise50 = tmp
		if(coltitle == 'temprise100'):
			temprise100 = tmp
		if(coltitle == 'temprise200'):
			temprise200 = tmp
		if(coltitle == 'GICperE'):
			GICperE = tmp
		if(coltitle == 'pop25to40'):
			pop25to40 = tmp
		if(coltitle == 'pop0to25'):
			pop0to25 = tmp
		if(coltitle == 'pop40plus'):
			pop40plus = tmp
		if(coltitle == 'tau'):
			tau = tmp