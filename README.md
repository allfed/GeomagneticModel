# GeoMagneticModel
The purpose of GeoMagneticModel is to report quantitative estimates of the likelihood of widespread power outages from severe geomagnetic storms.

## Description
Existing model outputs, including figures generated by the model, may be found in the 'Data/ModelOutput' directory. You may also run GeomagneticModel.py to recalculate and replot all aspects of the model, with custom command line parameters set directly as arguments to GeomagneticModel.py and with arguments imported from the ModelParameters.json file.

As an example, say you wanted to see a plot of the fraction of high voltage transformers disabled from a 'once per millenium' level storm. You would run: 

	$ python GeoMagneticModel.py -p TransformersInDanger -r .001

The -p / --plots flag allows you to easily choose which plots will show. The -r / --rate-per-year flag allows you to choose which rates per year to analyze, in units of years^-1.

Another example is if you wish to recalculate the E field over time and display a plot at a specific MT station. Such a plot requires specifying both an EMTF transfer function site, and a nearby MT magnetotelluric historical magnetic field record. You must supply an equal number of MT and TF sites for the program to succeed. They will be paired up in order. If neither MT or TF sites are specified, the arguments will be imported from the ModelParameters.json file. In this case although the default analysis is for once per century storms, the Efields plots include all rates per century up to once per year, and the --rate-per-year flag will not apply to this plot.

	$ python GeoMagneticModel.py MTsite -p EfieldFits -m VAQ55 -t FRD

To apply a power fit instead of a log normal fit you would run:

	$ python GeoMagneticModel.py MTsite -p EfieldFits -m VAQ55 -t FRD -f power
	
	-c, --conductivity
		Use the global conductivity dataset that you placed in  Data/BigData/GlobalConductivity/ and save the result as Data/Results/GlobalConductivity/userGC.npy.

		The script will always use userGC.npy if it is located in GlobalConductivityProcessed, and if userGC.npy doesn't exist defaultGC.npy will be used instead.

	-m, --mt-site
		Import data from magnetotelluric sites placed in BigData/MTsites and use as basis for data fit.

		A csv matching TFsites to MTsites is required for the program to run.  

#### Model modules

---

##### WARNING: out-of-date 

```GCmodelImport```, ```GCmodelProcess```, are not present in the codebase, and ```HVtransformers``` is neither used or usable (error in code).

---


| Module | Description |
| ----------- | ----------- |
| TFsite  | Allow estimation of geoelectric field near MTsites using electromagnetic transfer function (EMTF) ground measurements. |
| MTsite | Provide data over time for EW and NS magnetic field, and use TF sites to calculate corresponding geoelectric fields. |
| GCmodelImport | Importing the global conductivity model and determing apparent conductivity allows estimation of variation of geoelectric field over the earth by modelling global ground resistance. |
| GCmodelProcess | Process the global conductivity model into useful results. |
| HVtransformers | High voltage transformer properties including statistics on their relative prevalence, GIC thermal time constants and short term thermal sensitivities. |
| PowerGrid  | Very simple model of global power grid, countries and populations affected by outages. Provides fit function to the MT site data to with duration and repeat rate based field level at a site. |



#### Important functions
- TFsite 
	- ```importZ(path, row, column)```:
		
		imports impedance for the site at frequency specified by row, column
	- ```getApparentCond(Z,w)```:
	
		use impedance Z to calculate the apparent conductivity at the site, at frequency w
- MTsite: 
	- ```importFieldRecord(MTsiteName,minindex,maxindex)```:
	
		imports the magnetic field record from the netcdf in a sample set of indices of interest
	- ```calcEfields(TFsite)```:
	
		calculate the E field at the MT site given the transfer function site impedances as a function of frequency
	- ```getERepeatRates(nsamples,maxrate)```:
	
		for the particular MT site, what is the repeat rate for fields up to the maxrate specified (default 1 per year max rate)
	- ```fitERepeatRates(ERepeatRates)```:
	
		return a lognormal fit for a given plot
- GCmodelGenerate:
	- ```getGroundImpedance(layerconductivity)```:
	
		calculate ground impedance for a single location in a 1d layered earth model.
	- ```getApparentConductivity(impedance)```:
		
		calculate apparent conductivity from ground impedance
	- ```generateApparentConductivityMap(latInterval,longInterval, w, name=)```:
	
		generate apparent conductivity at each location for a frequency and save result into convenient compact file format
- GCmodelProcess:
	- Maglat

## Running the model with the PowerGrid parameter
### Environment set up

```
conda env create -n geomagmodel --file enivronment.yml
```
```
conda activate geomagmodel
```
- Get global conductivity model from https://globalconductivity.ocean.ru/matlabformat.html
- Put that file in Data/BigData/ (you need to make this directory yourself).
- Rename it to ```GCmodel.mat```.
- Make a Data/BigData/Networks directory.
- Get a poly file from https://download.geofabrik.de/.
- Rename it to pfile.poly and put in an appropriate folder (see below).
	- For example to run the model for Estonia we need this file https://download.geofabrik.de/europe/estonia.poly.
- Get transnet-model file(s) from https://github.com/OpenGridMap/transnet-models and put them in an appropriate folder (see below).
	- For Estonia, we need those files: https://github.com/OpenGridMap/transnet-models/tree/master/europe/estonia
- Make sure you have ```transnet``` and ```transnet-models``` directories one level UP from the root of geomagnetic models, with appropriate files in them, see example.
	- For the Estonia example you want this kind of directory structure:
	```
	|-- GeomagneticModelWrapperDirectory
	| |-- GeomagneticModel (this is the folder this README is in)
	| | |-- Data
	| | | |-- BigData
	| | | | |-- GCmodel.mat
	| | |-- ...
	| |-- transnet
	| | |-- data
	| | | |-- europe
	| | | | |-- estonia
	| | | | | |-- pfile.poly
	| | |-- configs
	| | | |-- estonia
	| |-- transnet-models
	| | |-- europe
	| | | |-- estonia
	| | | | |-- cim.xml
	| | | | |-- ...
	```

### Launch

Example 1.
```
python GeomagneticModel.py PowerGrid Region --continent europe --country estonia -r 0.01
```

Example 2.
```
python GeomagneticModel.py PowerGrid WorldNetwork -r 0.01
```

Example 3.
```
python GeomagneticModel.py PowerGrid compareGICresults --continent europe --country estonia -r 0.01
```

Example 4.
```
python GeomagneticModel.py PowerGrid LoadRegionE --continent europe --country estonia -r 0.01
```

### Result
Example 1.
- prints some info to stdout
- saves some data to pkl/npy files
- displays a plot of Estonia's power grid
- displays a plot of the E field over Estonia
- displays a plot of GICs in Estonia
- saves a power grid plot to Data/SmallData/Figures/Regions/estonia/Network.png
- saves a plot of E field to Data/SmallData/Figures/Regions/estonia/0.01peryearE.png
- saves a plot of GICs to Data/SmallData/Figures/Regions/estonia/0.01peryearGIC.png
- saves a plot of power outages to Data/SmallData/Figures/0.01peryearOutage.png, WARNING: file path collision
- saves a plot of overheating transformers to Data/SmallData/Figures/0.01peryearOverheat.png, WARNING: file path collision

Example 2.
- prints some info to stdout
- displays a plot of substations at risk of electricity loss for the whole world. Note: Colour bar is an overkill, result is binary: outage or no outage per sector.
- displays a plot of predicted electricity loss (as a fraction: amount lost over total consumption) by country. WARNING: doesn't display countries correctly.
- saves the first plot to Data/SmallData/Figures/0.01peryearOutage.png. WARNING: collision without other plots of this kind, see Example 1.
- creates CELEpop.pkl file

Example 3.

Compares GICs computed for the given continent (here, Europe) and country (here, Estonia).
- displays a plot of stations latitudes; x=continent, y=country
- displays a plot of stations longitudes; x=continent, y=country
- displays a plot of stations voltages; x=continent, y=country
- displays a plot of stations detected GICs; x=country, y=continent

Example 4.
- prints some info to stdout
- saves some data to pkl/npy files
- displays a plot of Estonia's power grid
- displays a plot of the E field over Estonia

### Notes/TODO

- overheat/outage plots are not saved to a region subfolder which can cause unwanted overwrites
- document PowerGrid printed output in more detail
- document other PowerGrid->compareGICresults, and PowerGrid->LoadRegionE
- the Estonia data mentioned in the manual is out of date
- what's the purpose of compareGICresults?