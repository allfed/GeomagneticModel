# A MTsite refers to a MT historical magnetic field record made by a magnetotelluric measurement site.
# This MT record of the B field over time is converted to a field record using the transfer functions determined in TFsite class
# Note: to find MT station fields over time, need to use https://supermag.jhuapl.edu/

from __future__ import print_function
import numpy as np
import scipy

from scipy import stats as scistats
from scipy.fftpack import fft, ifft
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from itertools import groupby
import recordtype
from recordtype import recordtype
import glob
from datetime import datetime, timedelta
from statsmodels.tsa.stattools import pacf
from statsmodels.tsa.ar_model import AutoReg
from statsmodels.tsa.ar_model import AutoReg, ar_select_order
import statsmodels.api as sm
import pandas as pd
from scipy.ndimage.filters import uniform_filter1d
import fits
from matplotlib.font_manager import FontProperties
import Params

u0 = 1.25664e-6  # m kg s^-2 A^-2 (permeability of free space, roughly same value inside the earth)
# also has units N A^-2.
secondsperyear = 31556925.216  # s (seconds in a solar year)


class MTsite:
    def __init__(self, MTsitefn, siteIndex):
        Params.importIfNotAlready()
        self.MTsitefn = MTsitefn
        self.chunks = []
        self.windowedCounts = []
        self.windowedCountsTmp = []  # temporary stored here before saving
        self.maxchunksize = Params.maxchunksize
        self.siteIndex = siteIndex
        self.sitename = Params.mtsitenames[siteIndex]
        self.betaThreshold = Params.betaThreshold[siteIndex]
        self.MTsiteChunk = recordtype(
            "MTsiteChunk",
            "chunkindex chunksize rawNS rawEW BfieldNS BfieldEW EfieldNS EfieldEW absE hourAvg rateperyearxy storms stormpeaks",
        )
        self.ds = []
        self.powerfits = []
        self.logfits = []
        self.maglat = 0
        self.N = 0
        self.nchunks = 0
        self.maxE = 0
        self.polyFitFunNS = []
        self.polyFitFunEW = []
        self.occurrenceRatio = []
        self.allstorms = []

    def importSite(self):
        self.ds = nc.Dataset(
            self.MTsitefn
        )  #  import the magnetic field record using the netcdf format as a big numpy array
        self.maglat = self.ds["mlat"][-1]
        self.lat = self.ds["glat"][-1].data
        self.long = self.ds["glon"][-1].data

        self.N = len(
            self.ds["dbn_geo"]
        )  # size of this particular MTsite magnetic field  (B field) record (in this case north, but north and south are always the same length)
        self.toKeepThreshold = Params.toKeepThreshold[self.siteIndex]
        self.sampleperiod = Params.sampleperiod[self.siteIndex]
        self.hourWindow = int(np.floor((60 * 60) / self.sampleperiod))
        self.cumulativeyears = (
            self.N * self.sampleperiod / secondsperyear
        )  # to be updated when processing ratesperyear
        if Params.nchunks[self.siteIndex]:
            self.nchunks = Params.nchunks[self.siteIndex]
        else:
            self.nchunks = int(np.floor(self.N / self.maxchunksize)) + 1
            print("self.nchunks")
            print(self.nchunks)
        for i in range(0, self.nchunks):
            self.chunks = self.chunks + [
                self.MTsiteChunk(0, 0, [], [], [], [], [], [], [], [], [], [], [])
            ]

        self.windows = []
        self.windowedCounts = []
        for i in range(0, len(Params.windows[self.siteIndex])):
            windowstring = Params.windows[self.siteIndex][i]
            if not windowstring:
                break
            window = int(windowstring)
            self.windows = self.windows + [window]
            self.windowedCounts = self.windowedCounts + [
                window * self.sampleperiod,
                [],
                [],
                [],
            ]
        self.calcPolyFits()

    # MTsites are usually so large, one must process them in smaller chunks for the fourier transform and convolution with the TF site frequency dependent transfer function to succeed. We also need a TF site for each MT site to determine the transfer function and thus estimate the geoelectric field.
    def createChunks(self):
        for chunkindex in range(0, self.nchunks):
            newchunk = self.createChunk(chunkindex)
        print("")
        print("")

    # #MTsites are usually so large, one must process them in smaller chunks for the fourier transform and convolution with the TF site frequency dependent transfer function to succeed. We also need a TF site for each MT site to determine the transfer function and thus estimate the geoelectric field.
    def createChunk(self, chunkindex):
        minindex = chunkindex * self.maxchunksize
        maxindex = min((chunkindex + 1) * self.maxchunksize - 1, self.N - 1)

        chunksize = maxindex - minindex + 1
        print(self.sitename + " MT site, c" + str(chunkindex))
        # print('importing chunkindex '+str(chunkindex)+', (chunk '+str(chunkindex+1)+' of '+str(self.nchunks)+')', end='\r')
        print(
            "importing chunkindex "
            + str(chunkindex)
            + ", (chunk "
            + str(chunkindex + 1)
            + " of "
            + str(self.nchunks)
            + ")"
        )

        # getBfield along maglat NS and EW
        rawNS = self.ds["dbn_geo"][minindex : maxindex + 1]
        rawEW = self.ds["dbe_geo"][minindex : maxindex + 1]

        self.chunks[chunkindex] = self.MTsiteChunk(
            chunkindex, chunksize, rawNS, rawEW, [], [], [], [], [], [], [], [], []
        )

    def cleanBfields(self):
        for chunkindex in range(0, self.nchunks):
            self.cleanChunkBfields(chunkindex)

    def cleanChunkBfields(self, chunkindex):
        chunk = self.chunks[chunkindex]
        rawNS = chunk.rawNS
        rawEW = chunk.rawEW

        if np.isnan(rawNS[0]):
            rawNS[0] = 0
        if np.isnan(rawNS[chunk.chunksize - 2]):
            rawNS[chunk.chunksize - 2] = 0

        if np.isnan(rawEW[0]):
            rawEW[0] = 0
        if np.isnan(rawEW[chunk.chunksize - 2]):
            rawEW[chunk.chunksize - 2] = 0

        # linearly interpolate any missing "nan" values

        indices = np.arange(len(rawNS))

        maskNS = np.isfinite(rawNS)
        maskEW = np.isfinite(rawEW)

        BfieldNS = np.interp(indices, indices[maskNS], rawNS[maskNS])
        BfieldEW = np.interp(indices, indices[maskEW], rawEW[maskEW])
        self.chunks[chunk.chunkindex].BfieldNS = BfieldNS
        self.chunks[chunk.chunkindex].BfieldEW = BfieldEW

    def plotEfields(self, chunkindex, startindex, endindex):
        # print('loading field from chunkindex '+str(chunkindex)+', chunk '+str(chunkindex+1)+' of '+str(self.nchunks), end='\r')
        print(
            "loading field from chunkindex "
            + str(chunkindex)
            + ", chunk "
            + str(chunkindex + 1)
            + " of "
            + str(self.nchunks)
        )
        print("")

        Efields = self.chunks[chunkindex].absE[startindex:endindex]
        startindexds = chunkindex * self.maxchunksize + startindex
        endindexds = chunkindex * self.maxchunksize + endindex

        dts = [str(datetime(2018, 1, 1))] * (endindexds - startindexds)
        elapsedHours = [0] * (endindexds - startindexds)

        years = np.array(self.ds["time_yr"][startindexds:endindexds]).astype(int)
        months = np.array(self.ds["time_mo"][startindexds:endindexds]).astype(int)
        days = np.array(self.ds["time_dy"][startindexds:endindexds]).astype(int)
        hours = np.array(self.ds["time_hr"][startindexds:endindexds]).astype(int)
        minutes = np.array(self.ds["time_mt"][startindexds:endindexds]).astype(int)
        seconds = np.array(self.ds["time_sc"][startindexds:endindexds]).astype(int)

        for i in range(0, len(years)):
            dts[i] = datetime(
                years[i], months[i], days[i], hours[i], minutes[i], seconds[i]
            )
            elapsedHours[i] = (dts[i] - dts[0]).total_seconds() / (60 * 60) + 12

        print("elapsedHours " + str(elapsedHours[0]))
        print("start time " + str(dts[0]))
        print("end time " + str(dts[-1]))
        print("mean(Efields)")
        print(np.mean(Efields))
        plt.figure()
        plt.yscale("log")
        plt.plot(elapsedHours, Efields)
        # beautify the x-labels
        plt.gcf().autofmt_xdate()
        plt.show()

    def saveChunkEfields(self, chunkindex):
        np.save(
            Params.mtEfieldsloc + str(self.sitename) + "c" + str(chunkindex),
            self.chunks[chunkindex].absE,
        )

    def loadEfields(self):
        for i in range(0, self.nchunks):
            # print('loading field from chunkindex '+str(i)+', chunk '+str(i+1)+' of '+str(self.nchunks), end='\r')
            print(
                "loading field from chunkindex "
                + str(i)
                + ", chunk "
                + str(i + 1)
                + " of "
                + str(self.nchunks)
            )
            print("")
            chunkE = np.load(
                Params.mtEfieldsloc + str(self.sitename) + "c" + str(i) + ".npy",
                allow_pickle=True,
            )
            print("chunkstart: " + str(i * self.maxchunksize))
            self.chunks[i].absE = chunkE
            self.chunks[i].chunksize = len(chunkE)

    # MTsites are usually so large, one must process them in smaller chunks for the fourier transform and convolution with the TF site frequency dependent transfer function to succeed. We also need a TF site for each MT site to determine the transfer function and thus estimate the geoelectric field.
    def calcEfields(self, TFsite):
        for i in range(0, len(self.nchunks)):
            self.calcChunkEfields(TFsite, i)
        print("")
        print("")

    # calculate 2nd order polynomial fits to NS and EW B data (downsample to one in 1000 points), assign array of fitted values to poly fit property
    def calcPolyFits(self):
        downsampleratio = 1000
        samplesNS = np.array(self.ds["dbn_geo"][0:-1:downsampleratio])
        samplesEW = np.array(self.ds["dbe_geo"][0:-1:downsampleratio])

        indices = np.arange(len(samplesNS))
        maskNS = np.isfinite(samplesNS)
        maskEW = np.isfinite(samplesEW)

        BfieldNS = np.interp(indices, indices[maskNS], samplesEW[maskEW])
        BfieldEW = np.interp(indices, indices[maskEW], samplesNS[maskNS])

        coeffsNS = np.polyfit(indices * downsampleratio, BfieldNS, 2)
        coeffsEW = np.polyfit(indices * downsampleratio, BfieldEW, 2)

        funNS = np.poly1d(coeffsNS)
        funEW = np.poly1d(coeffsEW)

        # plt.figure()
        # plt.plot(BfieldEW)
        # plt.plot(funEW(indices*downsampleratio))
        # plt.show()
        self.polyFitFunNS = funNS
        self.polyFitFunEW = funEW

    def calcChunkEfields(self, TFsite, chunkindex):
        chunk = self.chunks[chunkindex]
        chunksize = chunk.chunksize
        startindex = chunkindex * self.maxchunksize
        endindex = startindex + chunksize

        # see love, 2018 for the list of four corrections applied here.

        # first, subtract 2nd order polynomial fit

        indices = np.array(range(startindex, endindex))
        detrendedBN2 = chunk.BfieldNS
        detrendedBNS = chunk.BfieldNS - self.polyFitFunNS(indices)
        detrendedBEW = chunk.BfieldEW - self.polyFitFunEW(indices)

        # see page 381, Love 2019 for an explanation of the math below

        # second, apply FFT with 64 bit precision (fast fourier transform the field into freqency space)
        ftNS = fft(detrendedBNS, chunksize)
        ftEW = fft(detrendedBEW, chunksize)

        halflength = int(np.floor(chunksize / 2))

        # next, matrix multiply the tensor by the field
        # we do this by assigning each point on the b field freq array
        # as some Z value determined by linear interpolation
        fs = 1.0 / self.sampleperiod  # sampling frequency of 1/60 Hz

        freqrange = [x * fs / chunksize for x in range(1, halflength)]
        periodrange = [chunksize / (x * fs) for x in range(1, halflength)]

        # print('interpolating Z(w) for B(w), chunkindex '+str(chunk.chunkindex)+', '+str(chunk.chunkindex+1)+' of '+str(self.nchunks), end='\r')
        print(
            "interpolating Z(w) for B(w), chunkindex "
            + str(chunk.chunkindex)
            + ", "
            + str(chunk.chunkindex + 1)
            + " of "
            + str(self.nchunks)
        )

        ZXXforB = self.getInterpolation(freqrange, TFsite.ZXX, TFsite.f)
        ZXYforB = self.getInterpolation(freqrange, TFsite.ZXY, TFsite.f)
        ZYXforB = self.getInterpolation(freqrange, TFsite.ZYX, TFsite.f)
        ZYYforB = self.getInterpolation(freqrange, TFsite.ZYY, TFsite.f)

        # frequency domain
        EfieldNS_fd = [x * y for x, y in zip(ZXXforB, ftNS)] + [
            x * y for x, y in zip(ZXYforB, ftEW)
        ]
        EfieldEW_fd = [x * y for x, y in zip(ZYXforB, ftNS)] + [
            x * y for x, y in zip(ZYYforB, ftEW)
        ]

        # time domain
        EfieldNS = ifft(EfieldNS_fd, chunksize)
        EfieldEW = ifft(EfieldEW_fd, chunksize)

        # finally, determine the magnitude of E field over time
        absE = np.real(
            np.sqrt(EfieldNS * np.conj(EfieldNS) + EfieldEW * np.conj(EfieldEW))
        )

        self.chunks[chunk.chunkindex].EfieldNS = EfieldNS
        self.chunks[chunk.chunkindex].EfieldEW = EfieldEW
        self.chunks[chunk.chunkindex].absE = absE
        self.chunks[chunk.chunkindex].chunksize = chunksize

    def calcEratesPerYear(self, plotstorms):
        print("calcEratesPerYear")
        self.calcStorms(plotstorms)

        self.windowedCounts = []
        for windex in range(0, len(self.windows)):
            self.calcWindowEratesPerYear(windex)
        self.windowedCountsTmp = self.windowedCounts

    # combine the E rates calculated for all the chunks by windowed average into an array of e fields and a y axis value for rate per year of those e fields. This is a modification of (Love, 2018 page 9) which describes the overall process.
    def calcWindowEratesPerYear(self, windex):
        window = self.windows[windex]
        allE = []
        allcountsatE = []
        samplecount = 0
        for j in range(0, len(self.chunks)):
            peaks = self.chunks[j].stormpeaks[windex]
            # get the occurrence of High E rates per year vs E field level
            rateperyearxy = fits.getEfieldCounts(peaks)
            allE = np.append(allE, rateperyearxy[0])
            allcountsatE = np.append(allcountsatE, rateperyearxy[1])
            samplecount = samplecount + self.chunks[j].chunksize

        # total years measured
        self.cumulativeyears = (
            samplecount * self.sampleperiod / secondsperyear
        )  # total records (once per minute) by minutes in a year

        [Efinal, cumsumfinal, countsfinal] = fits.combinecounts(allE, allcountsatE)

        windowalreadyexists = False

        # add this duration data to the existing windowedCounts, or update if exists already
        for i in range(0, int(np.floor(len(self.windowedCounts)) / 5)):
            durationindex = i * 5
            efieldindex = i * 5 + 1
            countsindex = i * 5 + 2
            cumsumindex = i * 5 + 3
            netyearsindex = i * 5 + 4

            duration = self.windowedCounts[durationindex]
            if duration == window * self.sampleperiod:
                self.windowedCounts[efieldindex] = Efinal
                self.windowedCounts[countsindex] = countsfinal
                self.windowedCounts[cumsumindex] = cumsumfinal
                self.windowedCounts[netyearsindex] = self.cumulativeyears
                windowalreadyexists = True
                break

        if not windowalreadyexists:
            self.windowedCounts = self.windowedCounts + [
                window * self.sampleperiod,
                Efinal,
                countsfinal,
                cumsumfinal,
                self.cumulativeyears,
            ]

    # if E rates for year already exists for this site, add to it. Otherwise, create a new .npy to store Erates per year for this site
    def saveEratesPerYear(self):
        # add or replace the data for each window that has been created so far
        # if(self.loadWindowedCounts()):
        # 	loadeddurations=np.array([])
        # 	for i in range(0,int(np.floor(len(self.windowedCounts))/3)):
        # 		durationindex = i*3
        # 		efieldindex = i*3+1
        # 		ratesindex = i*3+2

        # 		#as long as the saved windows match up with modelparams specified windows, keep old values even if this window value was not just calculated.
        # 		if(durationindex<=len(self.windows) and (self.windows[i]*self.sampleperiod == self.windowedCounts[durationindex])):

        # 			duration = self.windowedCounts[durationindex]
        # 			loadeddurations=np.append(loadeddurations,duration)

        # 	for i in range(0,int(np.floor(len(self.windowedCountsTmp))/3)):
        # 		durationindex = i*3
        # 		efieldindex = i*3+1
        # 		ratesindex = i*3+2

        # 		duration = self.windowedCountsTmp[durationindex]
        # 		Efields = self.windowedCountsTmp[efieldindex]
        # 		rates = self.windowedCountsTmp[ratesindex]

        # 		if(len(rates)!=0):#if data has been processed for this window
        # 			loadeddurationindex=np.where(loadeddurations==duration)[0]

        # 			if(len(loadeddurationindex)==0):
        # 				self.windowedCounts=np.append(self.windowedCounts,[duration,Efields,rates])
        # 			else:
        # 				self.windowedCounts[loadeddurationindex[0]*3] = duration
        # 				self.windowedCounts[loadeddurationindex[0]*3+1]=Efields
        # 				self.windowedCounts[loadeddurationindex[0]*3+2]=rates
        # else:
        self.windowedCounts = self.windowedCountsTmp
        np.save(
            Params.mtRepeatRatesDir
            + "MTsite"
            + str(self.siteIndex)
            + "EfieldRatesPerYear",
            self.windowedCounts,
        )

    def loadWindowedCounts(self):
        loaddirectory = Params.mtRepeatRatesDir + "MTsite"
        allfiles = glob.glob(loaddirectory + "*.np[yz]")

        for f in allfiles:
            if (
                f
                == Params.mtRepeatRatesDir
                + "MTsite"
                + str(self.siteIndex)
                + "EfieldRatesPerYear.npy"
            ):
                wr = np.load(
                    Params.mtRepeatRatesDir
                    + "MTsite"
                    + str(self.siteIndex)
                    + "EfieldRatesPerYear.npy",
                    allow_pickle=True,
                )
                self.windowedCounts = wr
                return True
        return False

    def calcPeakEvsDuration(self, plot):
        durations = []
        highestEfields = []
        onceperyearE = []
        onceperfiveyearsE = []
        oncepersevenyearsE = []
        onceperdecadeE = []
        print("calcPeakEvsDuration")
        for i in range(0, int(np.floor(len(self.windowedCounts)) / 5)):
            durationindex = i * 5
            efieldindex = i * 5 + 1
            countsindex = i * 5 + 2
            cumsumindex = i * 5 + 3
            netyearsindex = i * 5 + 4

            duration = self.windowedCounts[durationindex]
            Efields = self.windowedCounts[efieldindex]
            counts = self.windowedCounts[countsindex]
            cumsum = self.windowedCounts[cumsumindex]
            rpy = cumsum / self.cumulativeyears

            highestEfields = np.append(highestEfields, np.max(Efields))
            durations = np.append(durations, duration)

            onceperyearE = np.append(onceperyearE, np.interp(10**0, rpy, Efields))
            onceperfiveyearsE = np.append(
                onceperfiveyearsE, np.interp(3 * 10**0, rpy, Efields)
            )
            oncepersevenyearsE = np.append(
                oncepersevenyearsE, np.interp(3 * 10**0, rpy, Efields)
            )
            onceperdecadeE = np.append(
                onceperdecadeE, np.interp(10**-1, rpy, Efields)
            )
        self.occurrenceRatio = onceperdecadeE / np.max(onceperdecadeE)

        if plot:
            print("self.occurrenceRatio")
            print(self.occurrenceRatio)
            plt.plot(
                durations,
                self.occurrenceRatio,
                lw=1,
                label="site " + str(self.sitename),
            )  # +"  beta "+str(self.betaThreshold))
        return self.occurrenceRatio

    def plotEratesPerYear(self, fittype):
        # plt.figure()
        plt.loglog()
        plots = np.array([])
        for i in range(0, int(np.floor(len(self.windowedCounts)) / 5)):
            if i != 0:
                break
            durationindex = i * 5
            efieldindex = i * 5 + 1
            countsindex = i * 5 + 2
            cumsumindex = i * 5 + 3
            netyearsindex = i * 5 + 4

            duration = self.windowedCounts[durationindex]
            Efields = self.windowedCounts[efieldindex]
            counts = self.windowedCounts[countsindex]
            cumsum = self.windowedCounts[cumsumindex]
            rpy = cumsum / self.cumulativeyears
            fitplotted = False
            # if we've calculated power fits for all the windows
            # plt.plot(Efields,rpy,lw=1,label = "Efields averaged over "+str(duration)+' seconds, '+str(self.sitename))
            if len(self.powerfits) > 0 and ("power" in fittype or fittype == "all"):
                fitplotted = True
                exponent = self.powerfits[i][1]
                ratio = self.powerfits[i][2]
                print("power coeffs")
                print([exponent, ratio])

                # plt.plot(Efields, fits.powerlaw(Efields,exponent)*ratio, lw=1,label = "Field averaged over "+str(duration)+" seconds, powerfit")

            if len(self.logfits) > 0 and ("lognormal" in fittype or fittype == "all"):
                fitplotted = True
                mean = self.logfits[i][1]
                std = self.logfits[i][2]
                loc = self.logfits[i][3]
                ratio = self.logfits[i][4]

                print("lognormal coeffs")
                print([mean, np.abs(std), loc, ratio])

                yfit = ratio * fits.logcdf(np.array(Efields), mean, np.abs(std), loc)

                boundedfit = np.array(yfit[np.array(yfit) > 10**-4])
                boundedEfields = np.array(Efields)[np.array(yfit) > 10**-4]
                plt.plot(boundedEfields, boundedfit, color="grey")
                p = plt.scatter(
                    Efields,
                    rpy,
                    lw=1,
                    label="Efields averaged over "
                    + str(duration)
                    + " seconds, "
                    + str(self.sitename),
                )
                plots = np.append(plots, p)
                # plt.plot(boundedEfields,boundedfit,lw=1,label = "Field lognormalfit")
            # if(not fitplotted):
            # plt.scatter(Efields,rpy,lw=1,label = "Field with lognormal fit")
            # plt.plot(Efields,rpy,lw=1,label = "Field averaged over "+str(duration)+" seconds")
        # print('Efields'-)
        # print(Efields)
        # print('rates')
        # print(rates)

        return plots
        # plt.show()

    def fitEratesPerYear(self):
        self.logfits = []
        self.powerfits = []
        for i in range(0, int(np.floor(len(self.windowedCounts)) / 5)):
            durationindex = i * 5
            efieldindex = i * 5 + 1
            countsindex = i * 5 + 2
            cumsumindex = i * 5 + 3
            netyearsindex = i * 5 + 4

            duration = self.windowedCounts[durationindex]
            Efields = self.windowedCounts[efieldindex]
            counts = self.windowedCounts[countsindex]
            cumsum = self.windowedCounts[cumsumindex]
            rpy = cumsum / self.cumulativeyears

            # get counts
            probtoRPYratio = np.max(cumsum) / self.cumulativeyears

            [exponent] = fits.fitPower(Efields, cumsum / np.max(cumsum))
            print("site " + str(self.MTsitefn) + " window: " + str(duration))
            # use PDF to determine mean and standard deviation of underlying normal distribution
            [guessMean, guessStd] = fits.getGuesses(Efields, counts, False)

            # fit to the datapoints (CDF has a probability of 1 at the first datapoint)
            [mean, std, loc] = fits.fitLognormalCDF(
                Efields, cumsum / np.max(cumsum), guessMean, guessStd, False
            )

            self.powerfits = self.powerfits + [[duration, exponent, probtoRPYratio]]
            self.logfits = self.logfits + [[duration, mean, std, loc, probtoRPYratio]]

            # plt.plot(Efields,rates, lw=1,label = "Field averaged over "+str(windowperiod)+" seconds")

    # convert counts at any field value to the cumulative rate per year (assumes counts are sorted by descending E field)
    def countstoRPY(self, counts):
        # dummy variable, does nothing
        dummy = range(0, len(counts))

        [_, cumsum, _] = fits.combinecounts(dummy, counts)

        return cumsum / self.cumulativeyears

    # def fitToRateperyear(self):
    # 	plt.figure()
    # 	plt.loglog()
    # 	for i in range(0,int(np.floor(len(self.windowedCounts))/3)):
    # 		durationindex = i*3
    # 		efieldindex = i*3+1
    # 		ratesindex = i*3+2

    # 		windowperiod = self.windowedCounts[durationindex]
    # 		Efields = self.windowedCounts[efieldindex]
    # 		rates = self.windowedCounts[ratesindex]
    # 	plt.legend()
    # 	plt.title('Rate geoelectric field is above threshold')
    # 	plt.xlabel('Geoelectric Field (V/km)')
    # 	plt.ylabel('Average rate per year distinct contiguous sample average is above E field (counts/year)')

    # 	plt.show()
    def calcStorms(self, plot):
        for i in range(0, self.nchunks):
            self.calcChunkStorms(i)

        # combine any storms that occurred between chunks
        for i in range(0, len(self.chunks)):
            if i > 0:
                if len(self.chunks[i].storms) == 0:
                    continue
                # if first element is part of a group
                if self.chunks[i].storms[0][0][0] == 0:
                    # combine this storm with the one at the end of the previous group if the last element of the previous group is also part of a storm
                    if self.chunks[i - 1].storms[-1][0][-1] == (
                        self.chunks[i - 1].chunksize - self.hourWindow
                    ):
                        print("chunk" + str(i))
                        stormextraminutes = len(self.chunks[i].storms[0][0])
                        stormstart = self.chunks[i - 1].storms[-1]
                        newindices = np.arange(
                            stormstart[0][0], stormstart[0][-1] + stormextraminutes
                        )
                        stormtoappend = self.chunks[i].storms[0][1]
                        self.chunks[i - 1].storms[-1][0] = newindices
                        self.chunks[i - 1].storms[-1][1] = np.append(
                            stormstart[1], stormtoappend
                        )

                        for j in range(0, len(self.windows)):
                            self.chunks[i - 1].stormpeaks[j][-1] = np.max(
                                [
                                    self.chunks[i - 1].stormpeaks[j][-1],
                                    self.chunks[i].stormpeaks[j][0],
                                ]
                            )
                        self.chunks[i].storms[0].pop()
        self.maxE = [0] * len(self.windows)
        for i in range(0, len(self.chunks)):
            for j in range(0, len(self.windows)):
                if len(self.chunks[i].stormpeaks[j]) == 0:
                    continue
                maxE = np.max(np.array(self.chunks[i].stormpeaks[j]))
                if maxE > self.maxE[j]:
                    self.maxE[j] = maxE

        # remove any storms that are below 1/10 the peak E field for all measurements (remove different storms for each window), but the original storm E fields are only kept for window0.
        nstormstotal = 0
        nwindows = len(self.windows)
        allstormdurations = []
        allstorms = [[], [], []]
        allE = np.array([])
        allHourAvg = np.array([])
        for i in range(0, len(self.chunks)):
            if plot:
                allE = np.append(allE, self.chunks[i].absE)
                allHourAvg = np.append(allHourAvg, self.chunks[i].absE)

            nstorms = 0
            stormstmp = []
            storms = self.chunks[i].storms
            peakstmp = [np.array([])] * nwindows
            for k in range(0, nwindows):
                sp = self.chunks[i].stormpeaks[k]
                for j in range(0, len(sp)):
                    if sp[j] > self.maxE[k] / self.toKeepThreshold:
                        peakstmp[k] = np.append(peakstmp[k], np.array(sp[j]))
                        if k == 0:
                            stormstmp.append(storms[j])
                            npoints = len(storms[j][1])
                            chunkoffset = self.maxchunksize * i
                            allstorms[0] = np.append(
                                np.array(allstorms[0]), storms[j][0] + chunkoffset
                            )
                            allstorms[1] = np.append(
                                np.array(allstorms[1]), storms[j][1]
                            )
                            hourAvg = self.chunks[i].hourAvg[storms[j][0]]
                            allstorms[2] = np.append(np.array(allstorms[2]), hourAvg)
                            allstormdurations.append(
                                npoints * self.sampleperiod / (60 * 60)
                            )
            self.chunks[i].stormpeaks = peakstmp
            self.chunks[i].storms = stormstmp
            if len(self.chunks[i].storms) > 0:
                nstorms = len(self.chunks[i].storms)
                nstormstotal = nstormstotal + nstorms

            print(str(nstorms) + " storms in chunk index " + str(i))

            if plot:
                if (
                    i != 0
                ):  # only take the first chunk for this plot, more is unnecessary
                    continue
                plt.figure()
                plt.plot(np.log(np.ones(len(allE)) * self.maxE[0]))
                plt.plot(np.log(np.ones(len(allE)) * self.maxE[0] / 10))
                plt.plot(np.log(allE))
                plt.plot(np.log(allHourAvg))
                plt.plot(np.log(uniform_filter1d(allE, size=self.hourWindow)))
                plt.plot(allstorms[0], np.log(allstorms[1]))
                plt.plot(allstorms[0], np.log(allstorms[2]))

                plt.show()

        meanduration = np.mean(allstormdurations)

        if plot:
            print("median E field")
            print(np.median(allE))
            print("log threshold")
            print(np.log(self.betaThreshold))
            print("log maxE")
            print(np.log(self.maxE[0] / 10))
            counts = np.ones(len(allstormdurations))
            [durationsx, bins] = fits.binnormaldist(allstormdurations, counts, 5)
            plt.figure()
            plt.xlabel("hours")
            plt.ylabel("percentage")
            plt.plot(durationsx, bins)
            plt.show()

        print("number storms per year: " + str(nstormstotal / self.cumulativeyears))
        print("average storms duration (hours): " + str(meanduration))
        self.allstorms = allstorms
        np.save(
            Params.mtStormsDir + "MTsite" + str(self.siteIndex) + "Storms", allstorms
        )
        print("saved storms")

    # calculates the peak value of storms given the windows and E fields
    # returns an array structured as:
    # [Storms,Peaks]
    # where
    # 	Storms=[[storm indices],[storm E fields],...]
    # 	Peaks=[[peaks window 0],[peaks window 1],...]
    # and storm indices, storm E fields, and peaks window [x] are numpy arrays, but the outer layer which combines them is the default array type.
    # each of these datatypes is stored for each chunk
    def calcChunkStorms(self, chunkindex):
        chunk = self.chunks[chunkindex]
        # we lose nchunks hours out of entire dataset (one hour at the end of each chunk) by using windowed average, but this should be fine.
        hourAvg = uniform_filter1d(chunk.absE, size=self.hourWindow)
        self.chunks[chunkindex].hourAvg = hourAvg
        # print('chunk.absE')
        # plt.figure()
        # plt.yscale("log")
        # plt.plot(chunk.absE[100000:1000100])
        # plt.plot(hourAvg[100000:1000100])
        # plt.show()
        # self.plotEfields(chunkindex,732000-5*60-20+12*60,732000-5*60-20+108*60)
        nwindows = len(self.windows)
        Efields = [np.array([])] * nwindows
        for i in range(0, nwindows):
            smoothed = uniform_filter1d(chunk.absE, size=self.windows[i])
            Efields[i] = np.array(smoothed)
        exceeds = np.array(hourAvg) > self.betaThreshold
        data = zip(exceeds, Efields[0])
        Eindex = 0
        storms = []
        stormpeaks = [np.array([])] * nwindows
        # print('chunkindex')
        # print(chunkindex)
        print("getting all the storms in chunk " + str(chunkindex))

        for key, group in groupby(data, lambda x: x[0]):
            isstorm, field = next(group)
            # print('nWindows')
            # print(nwindows)
            elems = len(list(group)) + 1
            stormE = [np.array([])] * nwindows
            stormpeak = [0] * nwindows
            if isstorm:
                Eindices = np.arange(Eindex, Eindex + elems)
                for i in range(0, nwindows):
                    stormE[i] = Efields[i][Eindex : Eindex + elems]
                    stormpeak[i] = np.max(np.array(stormE[i]))
                    stormpeaks[i] = np.append(stormpeaks[i], stormpeak[i])

                storms.append([Eindices, stormE[0]])
            Eindex += elems

        self.chunks[chunkindex].storms = storms
        self.chunks[chunkindex].stormpeaks = stormpeaks

    # The purpose of this function is to return an array of
    # values of impedance for each frequency bin in the frequency
    # domain B field. This allows us to multiply the tensor
    # components of the frequency domain B field with the
    def getInterpolation(self, BfieldFreqs, Z, Zfreqs):
        sortedBFreqs = sorted(BfieldFreqs)
        sortedZFreqs = sorted(Zfreqs)
        # plt.figure()
        # plt.loglog()
        # print('sortedBfreqs')
        # plt.plot(sortedBFreqs)
        # plt.show()

        # plt.figure()
        # plt.loglog()
        # print('sortedZfreqs')
        # plt.plot(sortedZFreqs)
        # plt.show()

        sortedZbyFreq = [x for _, x in sorted(zip(Zfreqs, Z))]
        # plt.show()
        # plt.figure()
        # plt.loglog()
        # print('sortedZfreqs')
        # plt.plot(sortedZFreqs,sortedZbyFreq)
        # plt.show()
        toInterpolateFreqs = []
        toInterpolateZ = []

        # we only use well-defined frequencies from 10^-1 (once per 10 seconds) to 10^-4. All frequencies outside of that are set to zero impedance.
        bandpass = np.logical_and(
            10**-1 > np.array(sortedZFreqs), np.array(sortedZFreqs) > 10**-4
        )

        bandpassedZ = np.array(sortedZbyFreq) * bandpass
        if sortedZFreqs[0] > 10**-4 or sortedZFreqs[-1] < 10**-2:
            print("ERROR!!!: the TFsite has an unsuitable freqency range")
            quit()

        # Furthermore, we need to ensure there is no amplification of frequencies higher than the half sampling rate (nyquist limit). Any Z values higher than this frequency is set to zero impedance.

        lowpass = 0.5 / self.sampleperiod > np.array(sortedZFreqs)
        lowpassedZ = np.array(sortedZbyFreq) * lowpass

        # if freq is lower,
        if sortedBFreqs[0] < sortedZFreqs[0]:
            toInterpolateFreqs = [sortedBFreqs[0]]
            toInterpolateZ = np.array([0])  # [sortedZbyFreq[0]]

        toInterpolateFreqs = np.append(toInterpolateFreqs, sortedZFreqs)
        toInterpolateZ = np.append(toInterpolateZ, lowpassedZ)

        # if freq is higher,
        if sortedBFreqs[-1] > sortedZFreqs[-1]:
            toInterpolateFreqs = np.append(toInterpolateFreqs, [sortedBFreqs[-1]])
            toInterpolateZ = np.append(oInterpolateZ, np.array([0]))

        # plt.show()
        # plt.figure()
        # plt.loglog()
        # print('sortedZfreqs')
        # plt.plot(toInterpolateFreqs,toInterpolateZ)
        # plt.show()

        # return the interpolating function with Z values for every Bfield frequency
        interpolatedfun = interp1d(toInterpolateFreqs, toInterpolateZ)
        ZforBfield = interpolatedfun(BfieldFreqs)

        # plt.figure()
        # plt.loglog()
        # plt.plot(BfieldFreqs,ZforBfield)
        # print('plt.plot(BfieldFreqs,ZforBfield)')
        # plt.show()

        # plt.figure()
        # plt.loglog()
        # print('plt.plot(toInterpolateFreqs,toInterpolateZ)')
        # plt.plot(toInterpolateFreqs,np.real(toInterpolateZ))
        # plt.show()
        return ZforBfield

    def Estats(self):
        E = self.chunks[1].absE
        normE = np.log(E)
        mean = np.mean(normE)
        std = np.std(normE)

        normBeta = (np.log(self.betaThreshold) - mean) / std

        print("mean")
        print(mean)
        print("std")
        print(std)
        # plt.plot(normal,'.')
        # plt.plot(normE)
        normData = (normE - mean) / std
        print("betaThreshold")
        print(normBeta)
        np.save("examplechunk", normData)
        # plt.figure()
        # # pcorrs=pacf(normData,15)
        # plt.plot(normData)
        # # plt.yscale("log")
        # plt.show()

        # model = AutoReg(normData, lags=15)
        # model_fit = model.fit()
        # coef = model_fit.params
        # print('coef')
        # print(coef)
        # print('pcorrs')
        # print(pcorrs)

        # sel = ar_select_order(normData[0:10000], 13, glob=True, old_names=False)
        # sel.ar_lags
        # res = sel.model.fit()
        # print(res.summary())

        # model = sm.tsa.ARMA(normData, (15, 5)).fit(trend='nc', disp=0)
        return
        # print('model.params')
        # print(model.params)
        # plt.figure()
        # plt.plot(model.params)
        # plt.show()

        nbins = 50
        bins = [0] * nbins
        minimum = -4
        maximum = 4
        for i in range(0, nbins):
            low = minimum + (maximum - minimum) / nbins * i
            high = minimum + (maximum - minimum) / nbins * (i + 1)
            mask = np.logical_and(normData > low, normData < high)
            # mask=np.logical_and(white>low,white<high)
            bins[i] = np.sum(mask)

        # gets autocorrelation of underlying normally distributed E fields
        # corrlag1=fits.compute_corr_lag_1(normData)
        # print('corrlag1')
        # print(corrlag1)

        plt.figure()
        # plt.yscale("log")
        # plt.plot(normal,'.')
        plt.plot(bins)
        plt.show()

    def loadStorms(self):
        self.allstorms = np.load(
            Params.mtStormsDir + "MTsite" + str(self.siteIndex) + "Storms.npy",
            allow_pickle=True,
        )

    def calcTimeBetweenStorms(self):
        # for i in range(0,self.chunks):
        i = 0
        # stormindices=self.chunks[i].storms[0][0]

        differences = np.diff(self.allstorms)
        differences0 = np.diff(self.allstorms[0])

        # differences>0

        # for i in range(0,len(stormindices)):
        # if(differences[i]>0):

        ranges = []

        # for x in xrange(1,10):
        # pass

        # for key, group in groupby(np.array(self.allstorms)>0):
        # print(key)
        # group = map(itemgetter(1), group)
        # if len(group) > 1:
        # ranges.append(xrange(group[0], group[-1]))
        # else:
        # ranges.append(group[0])

        # for group in groupby(stormindices, lambda x: x[0]):
        # 	isstorm,field = next(group)
        # 	# print('nWindows')
        # 	# print(nwindows)
        # 	elems = len(list(group)) + 1
        # 	stormE=[np.array([])]*nwindows
        # df=pd.DataFrame({'indices':np.array(self.allstorms[0]),'vals':np.array(self.allstorms[1]),'hourwin':np.array(self.allstorms[2])})
        tmp = np.append(differences0, 0) == 1
        mask = tmp
        print(mask)
        print(np.sum(mask))
        # for i in range(1,len(tmp)):
        # 	if(tmp[i]!=tmp[i-1]):
        # 		mask[i]=True

        storms = []
        # stormpeaks=[np.array([])]*nwindows
        # print('chunkindex')
        # print(chunkindex)
        stormstarts = []
        stormends = []
        stormdurations = []
        peaks = []
        data = zip(mask, self.allstorms[0])
        Eindex = 0
        for key, group in groupby(data, lambda x: x[0]):
            isstorm, stormstart = next(group)
            elems = len(list(group)) + 1
            if isstorm:
                stormE = self.allstorms[1][Eindex : Eindex + elems]
                stormindices = self.allstorms[0][Eindex : Eindex + elems]
                # print('stormstart')
                # print(stormstart)
                # print('len(stormE)')
                # print(len(stormE))
                peaks.append(np.max(np.array(stormE)))

                stormstarts.append(stormindices[0])
                stormends.append(stormindices[-1])
                stormdurations.append(stormindices[-1] - stormindices[0])
            Eindex += elems

        timesbetween = np.array(stormstarts[1:]) - np.array(stormends[0:-1])
        print("hoursbetween")
        hoursbetween = timesbetween / (60)
        print(hoursbetween)
        [x, y] = fits.binlognormaldist(hoursbetween, [], 3)
        # [xd,yd]=fits.binlognormaldist(np.array(stormdurations)/(60*60),[],3)
        plt.figure()

        plt.xscale("log")

        plt.plot(x, y)
        plt.show()

        # plt.figure()
        # plt.xscale("log")
        # plt.plot(xd,yd)
        # plt.show()

        #
