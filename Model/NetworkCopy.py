import Model.GIC_Model
from Model.GIC_Model import GIC_Model

from Plotter import Plotter
from PolyParser import PolyParser
import Params
import re
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, LineString
import matplotlib.pyplot as plt
import numpy as np
from numpy import savetxt
import PyCIM
from PyCIM import cimread
import xml.etree.ElementTree as ET

# import CIM14.IEC61970.Wires as Wires
import CIM15.IEC61970.Core.ConnectivityNode as ConnectivityNode
import CIM15.IEC61970.Wires.ACLineSegment as ACLineSegment
import CIM15.IEC61970.Core.Substation as SubStation
import CIM15.IEC61970.Wires.SynchronousMachine as SynchronousMachine
from shapely.geometry import LineString
from matplotlib.colors import ListedColormap
import cycler
import os
import sys
import subprocess


class Network:
    def __init__(self):
        # assumptions:
        # each transformer gets its own winding resistance, and it's assumed that transformers are all connected to the same grounding grid at a substation, with a given resistance per phase of the windings (transformertogroundR). These are all connected together, and there's an earthing resistance (groundingR) for the grounding grid to earth.
        self.lowestAllowedVoltage = 100  # kV
        self.groundingR = "1"  # ohms
        # self.transformertogroundR='.15'#ohms
        self.allvoltageset = set({})
        self.region = ""
        self.transformerVoltages = {}
        self.maxlat = -100000
        self.minlat = 100000
        self.maxlong = -100000
        self.minlong = 100000
        self.EfieldFiles = {}
        # https://www.nerc.com/pa/Stand/Project201303GeomagneticDisturbanceMitigation/ApplicableNetwork.pdf

        windingkV = [230, 345, 500, 735]
        windingConductance = 1 / np.array([0.692, 0.356, 0.195, 0.159])  # Ohms
        self.transformertogroundR = {}
        self.transformertogroundR[275] = 1 / np.interp(
            275, windingkV, windingConductance
        )
        self.transformertogroundR[400] = 1 / np.interp(
            400, windingkV, windingConductance
        )

        linekV = [110, 230, 345, 500, 735]  # kV
        lineConductance = 1 / np.array(
            [0.037 * 3, 0.072, 0.037, 0.013, 0.011]
        )  # Ohms/km. Later divided by 3 in GIC_Model
        self.linekV = linekV
        self.lineConductance = lineConductance
        # plt.plot(linekV,lineConductance)
        # plt.show()
        self.lineResistance = {}
        self.lineResistance[750] = 1 / np.interp(750, linekV, lineConductance)
        self.lineResistance[500] = 1 / np.interp(500, linekV, lineConductance)
        self.lineResistance[450] = 1 / np.interp(450, linekV, lineConductance)
        self.lineResistance[420] = 1 / np.interp(420, linekV, lineConductance)
        self.lineResistance[400] = 1 / np.interp(400, linekV, lineConductance)
        self.lineResistance[380] = 1 / np.interp(380, linekV, lineConductance)
        self.lineResistance[330] = 1 / np.interp(330, linekV, lineConductance)
        self.lineResistance[300] = 1 / np.interp(300, linekV, lineConductance)
        self.lineResistance[275] = 1 / np.interp(275, linekV, lineConductance)
        self.lineResistance[220] = 1 / np.interp(220, linekV, lineConductance)
        self.lineResistance[225] = 1 / np.interp(225, linekV, lineConductance)
        self.lineResistance[115] = 1 / np.interp(115, linekV, lineConductance)
        self.lineResistance[110] = 1 / np.interp(110, linekV, lineConductance)
        self.lineResistance[132] = 1 / np.interp(110, linekV, lineConductance)

    def importNetwork(self, continent, country):
        if len(country) == 0:
            self.region = continent
            self.regionDataDir = Params.planetNetworkDir + continent + "/"
            self.downloadedDataDir = Params.downloadedTransnetDataDir + continent + "/"
        else:
            self.region = country
            self.regionDataDir = (
                Params.countryNetworkDir + continent + "/" + country + "/"
            )
            self.downloadedDataDir = (
                Params.downloadedTransnetDataDir + continent + "/" + country + "/"
            )
        self.loadPolygon()
        self.continent = continent
        # self.regionDataDir='europe/luxembourg/'

        self.processedNetworkDir = Params.networkSaveDir

        # self.country='luxembourg'
        self.country = country
        self.continent = continent
        self.dbname = self.region.replace("-", "_")
        self.connections = []
        self.nodes = []
        self.nodesDict = {}
        self.cleanUpXML()
        self.rawNetworkData = cimread(self.regionDataDir + "cimClean.xml")

        self.conns = []
        self.gens = []
        self.voltages = []
        self.lines = []
        # import values from transnet data
        for rnd in self.rawNetworkData.values():
            #  the power line data
            if isinstance(rnd, ACLineSegment):
                self.conns.append(rnd)
                # print(self.conns)
            # power station data
            if isinstance(rnd, SynchronousMachine):
                self.gens.append(rnd)

        nodesDict = self.importNodes()
        self.sortednodes = {
            v[0]: v
            for k, v in sorted(self.nodesDict.items(), key=lambda item: item[1][0])
        }
        # self.sortednodes=[sortednodesDict.values()]

        self.importConnections(nodesDict)
        self.saveTextFiles()
        self.generateConfigFile()
        print("")
        print("")
        print("")
        print("")
        print("")
        # self.analyzeNetwork()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # Plotter.plotNetwork(self.voltages,self.lines,ax)
        return [self.voltages, self.lines]

    def loadPolygon(self):
        pfile = self.downloadedDataDir + "pfile.poly"
        self.boundaryPolygon = PolyParser.poly_to_polygon(pfile)

    def cleanUpXML(self):
        print("cleaning up the xml doc")
        f1 = open(self.regionDataDir + "cim.xml", "r")
        f2 = open(self.regionDataDir + "cimClean.xml", "w")

        def addunderscore(matchobj):
            return "_CN_" + matchobj.group(0)[-1]

        for line in f1:
            # f2.write(line.replace('>None<', '>NaN<'))
            f2.write(re.sub("_CN\d", addunderscore, line.replace(">None<", ">NaN<")))
        f1.close()
        f2.close()
        print("cleanupxml done")

    def importNodes(self):
        tree = ET.parse(self.regionDataDir + "cim.xml")
        root = tree.getroot()

        conns = []
        output = []
        windings = []
        nodes = []
        nodesDict = {}

        # NOTE: unfortunately PyCIM does not appear to import the Transformerwinding object, so I had to hack together my own import function for this property.
        index = 0
        for i in root:
            if i.tag.split("}")[-1] == "TransformerWinding":
                windings.append(i)

        # obtain transformer winding properties, build a node for each substation ground (_E) and voltage (_[voltage value])
        for w in windings:
            for info in w:
                # name of the winding
                if info.tag.split("}")[-1] == "IdentifiedObject.name":
                    windname = info.text
                if info.tag.split("}")[-1] == "TransformerWinding.ratedU":
                    voltage = float(info.text) / 1000  # kV
                if info.tag.split("}")[-1] == "TransformerWinding.r":
                    r = float(info.text)
                if info.tag.split("}")[-1] == "TransformerWinding.PowerTransformer":
                    tid = list(info.attrib.values())[0]
                    uuid = tid[1:]
            # 	if(info.tag.split('}')[-1]=='TransformerWinding.windingType'):
            # 		windtype=list(info.attrib.values())[0]

            # #we only care about primary windings, because we've already specified the winding by its voltage in the windname property
            # if(windtype.split('.')[-1]!="primary"):
            # 	continue

            print("here")
            # we only care about primary windings, because we've already specified the winding by its voltage in the windname property
            if voltage < self.lowestAllowedVoltage:
                continue

            # build a dictionary for finding the substation name from the transformer winding name
            transformer = self.rawNetworkData.get(uuid)
            substation = transformer.getEquipmentContainer()
            subname = substation.name
            lat = substation.getLocation().getPositionPoints()[0].yPosition
            long = substation.getLocation().getPositionPoints()[0].xPosition
            if float(lat) > self.maxlat:
                self.maxlat = float(lat)
            if float(long) > self.maxlong:
                self.maxlong = float(long)
            if float(long) < self.minlong:
                self.minlong = float(long)
            if float(lat) < self.minlat:
                self.minlat = float(lat)
            # every winding needs a connection to ground node of substation

            # every substation has one ground node
            # node 5 is resistance directly to ground from this node
            # node 6 is resistance directly to transformer from this node (I think, not really sure about this one)

            ## HERE WE CHECK IF SUBSTATION HAS APPEARED, IF SO USE MAX VOLTAGE
            if not (subname in nodesDict.keys()):
                index = index + 1
                # data format for ground node:
                # __0__ ____1____ ___2___ ___3___ _4_ __5__ _______6______ _7_
                # index node name  code   country lat long  gnd resistance  0
                nodes.append(
                    np.array(
                        [
                            index,
                            subname + "_E",
                            index,
                            self.region,
                            lat,
                            long,
                            self.groundingR,
                            "0",
                        ]
                    )
                )

                nodesDict[subname] = [
                    index,
                    lat,
                    long,
                    "substation",
                    "",
                    voltage,
                ]  # last index is voltage of maximum transformer winding voltage at substation
                subindex = index
            elif nodesDict[subname][5] < voltage:
                nodesDict[subname][5] = voltage

            # data format for transformer (winding) node:
            # __0__ ____1____ ___2___ ___3___ _4_ __5__ _______6______ _7_
            # index node name  code   country lat long  gnd resistance Inf
            index = index + 1
            nodes.append(
                np.array([index, windname, index, self.region, lat, long, "Inf", "Inf"])
            )

            # as it stands, there can only be one winding per voltage network for the substation. This may need to be adjusted if there were more complicated representations of networks.
            #
            # assert(not (windname in nodesDict.keys()))
            nodesDict[windname] = [index, lat, long, "winding", subname, voltage]

        # add power station nodes
        # for g in self.gens:
        # 	voltage=int(g.getBaseVoltage().nominalVoltage/1000)
        # 	genname=g.getGeneratingUnit().name
        # 	genlat=g.getGeneratingUnit().getLocation().getPositionPoints()[0].yPosition
        # 	genlong=g.getGeneratingUnit().getLocation().getPositionPoints()[0].xPosition

        # 	index=index+1

        # 	# data format for power station transformer node:
        # 	# __0__ ____1____ ___2___ ___2___ _3_ __4__ _______5______ _6_
        # 	# index node name  code   country lat long  gnd resistance Inf
        # 	nodes.append(np.array([index,genname,index,self.region,genlat,genlong,'Inf','Inf']))

        # 	nodesDict[genname]=[index,genlat,genlong,'generatingUnit',genname,voltage]
        self.nodes = nodes
        self.nodesDict = nodesDict
        return nodesDict

    def importConnections(self, nodesDict):
        conns = []

        # NOTE: again because of the bug with PyCIM, am hacking together a reference by using the names of the AC power lines, rather than the TransformerWinding(A) - Terminal - ConnectivityNode - Terminal - ACLineSegment -Terminal - ConnectivityNode - Terminal - TransformerWinding(B) UUID terminal chain from the cim.xml.
        outputConnections = []

        gdfs = []  # for plotting
        voltages = []
        lines = []
        startpoints = []
        stoppoints = []
        isConnectedToGround = {}
        index = 0
        notin = 0
        for c in self.conns:
            # print('self.conns')
            # print(c)
            v = int(c.getBaseVoltage().nominalVoltage / 1000)
            if v < self.lowestAllowedVoltage:
                continue
            # print(c)
            # print(c.name)
            winding1 = "TW_" + c.name.split("_CN_")[1]
            winding2 = "TW_" + c.name.split("_CN_")[2]

            if winding1 in nodesDict.keys():
                assert nodesDict[winding1][3] == "winding"
                # get the index for each winding
                index1 = nodesDict[winding1][0]
                subname1 = nodesDict[winding1][4]
                key1 = winding1
            else:
                continue
                # #there isn't a transformer for this connection. It must connect directly to a generating unit (aka power station). The gen units don't specify their voltage in the name.
                # gen1='G_'+c.name.split('_CN_')[1].split('_')[0]
                # # print('gen1')
                # # print(gen1)
                # if(gen1 in nodesDict.keys()):
                # 	assert(nodesDict[gen1][3]=='generatingUnit')
                # 	index1=nodesDict[gen1][0]
                # 	subname1=nodesDict[gen1][4]
                # 	key1=gen1
                # else:
                # 	print('Error: ACLineSegment'+str(c.name)+' does not connect to node!!!')

            if winding2 in nodesDict.keys():
                assert nodesDict[winding2][3] == "winding"
                # get the index for each winding
                index2 = nodesDict[winding2][0]
                subname2 = nodesDict[winding2][4]
                key2 = winding2
            else:
                continue
                # #there isn't a transformer for this connection. It must connect directly to a generating unit (aka power station). The gen units don't specify their voltage in the name.
                # gen2='G_'+c.name.split('_CN_')[2].split('_')[0]
                # if(gen2 in nodesDict.keys()):
                # 	assert(nodesDict[gen2][3]=='generatingUnit')
                # 	index2=nodesDict[gen2][0]
                # 	subname2=nodesDict[gen2][4]
                # 	key2=gen2
                # else:
                # 	print('Error: ACLineSegment'+str(c.name)+' does not connect to node!!!')
                # 	quit()
            self.allvoltageset.add(v)

            if not (v in self.lineResistance.keys()):
                self.lineResistance[v] = 1 / np.interp(
                    v, self.linekV, self.lineConductance
                )
            resistancePerKilometer = self.lineResistance[v]

            if not (v in self.transformertogroundR.keys()):
                self.transformertogroundR[v] = 1 / np.interp(
                    v, self.windingkV, self.windingConductance
                )
            windR = self.transformertogroundR[v]
            # print('resistancePerKilometer')
            # print(resistancePerKilometer)
            # print('windR')
            # print(windR)
            # print('voltage')
            # print(v)

            # print('c.length')
            # print(c.length)

            r = resistancePerKilometer * (
                c.length / 1000
            )  # Ohms (net resistance of the power line, length converted to kilometers)
            print("r")
            print(r)

            index = index + 1

            # data format for connection to other transformer node:
            # __0__ ____1____ ___2___ _3_ _4_ _____5____ ____6____
            # index nodefrom  nodeto   0   0  resistance  voltage
            outputConnections.append([index, index1, index2, "0", "0", str(r), str(v)])

            if not (index1 in isConnectedToGround.keys()):
                assert nodesDict[subname1][3] == "substation"
                sub1index = nodesDict[subname1][0]

                index = index + 1

                # data format for connection transformer winding node to ground node:
                # __0__ ____1____ ___2___ _3_ _4_ ________5_______ _6_
                # index nodefrom  nodeto   0   0  wind resistance   0
                outputConnections.append(
                    [index, index1, sub1index, "0", "0", windR, "0"]
                )

            # if we haven't connected winding2 to the substation ground yet
            if not (index2 in isConnectedToGround.keys()):
                sub2index = nodesDict[subname2][0]

                index = index + 1

                # data format for connection transformer winding node to ground node:
                # __0__ ____1____ ___2___ _3_ _4_ ________5_______ _6_
                # index nodefrom  nodeto   0   0  wind resistance   0
                outputConnections.append(
                    [index, index2, sub2index, "0", "0", windR, "0"]
                )

            # now we've connected the nodes connected to this power line to substation ground
            isConnectedToGround[index1] = True
            isConnectedToGround[index2] = True

            # for plotting
            lat1 = nodesDict[key1][1]
            long1 = nodesDict[key1][2]
            lat2 = nodesDict[key2][1]
            long2 = nodesDict[key2][2]
            p1 = Point(float(long1), float(lat1))
            p2 = Point(float(long2), float(lat2))
            line = LineString([[p1.x, p1.y], [p2.x, p2.y]])
            print("[[p1.x, p1.y], [p2.x, p2.y]]")
            print([[p1.x, p1.y], [p2.x, p2.y]])
            vindex = np.where(np.array(voltages) == v)[0]
            if len(vindex) == 0:
                voltages.append(v)
                lines.append([])
                i = len(voltages) - 1
            else:
                i = vindex[0]
            lines[i].append(line)

        # quit()

        self.connections = outputConnections
        self.lines = lines
        self.voltages = voltages
        self.nodesDict = nodesDict

    def analyzeNetwork(self):
        print("OOOO")
        lines = []
        points = []

        newnodes = {}
        newconnections = []
        # sortednodes={v[0]:v for k, v in sorted(self.nodesDict.items(), key=lambda item: item[1][0])}
        # print(sortednodes)
        for nk in self.nodesDict.keys():
            node = self.nodesDict[nk]
            point = Point([float(node[1]), float(node[2])])
            points.append(point)

            # print(node)
            if (
                float(node[2]) < 51.1
                and float(node[2]) > 50.2
                and float(node[1]) < -3
                and float(node[1]) > -5
            ):
                newnodes[nk] = node
                # print(node)
                for c in self.connections:
                    if c[1] == node[0] or c[2] == node[0]:
                        add = True

                        # print(c)
                        for nc in newconnections:
                            if nc[0] == c[0]:
                                add = False
                                break
                        if add:
                            newconnections.append(c)
                            # print(c[1])
                            # print('c[2]')
                            # print(c[2])
                            # print('sortednodes[c[1]]')
                            # print(sortednodes[c[1]])
                            # print('sortednodes[c[2]]')
                            # print(sortednodes[c[2]])

                            p1 = Point(
                                float(sortednodes[c[1]][2]), float(sortednodes[c[1]][1])
                            )
                            p2 = Point(
                                float(sortednodes[c[2]][2]), float(sortednodes[c[2]][1])
                            )
                            line = LineString([[p1.x, p1.y], [p2.x, p2.y]])
                            lines.append(line)
        # quit()
        print("newnodes")
        # print(newnodes)
        print("newconnections")
        # print(newconnections)
        print("OOOO")

        network = gpd.GeoDataFrame({"geometry": np.array(lines)})
        gdf = gpd.GeoDataFrame(geometry=lines)

        gdf.plot()
        plt.show()
        # quit()

    def saveTextFiles(self):
        # save the text files which will be processed

        savetxt(
            self.processedNetworkDir + self.region + "Nodes.txt",
            np.array(self.nodes),
            delimiter="\t",
            fmt="%s",
        )
        savetxt(
            self.processedNetworkDir + self.region + "Connections.txt",
            np.array(self.connections),
            delimiter="\t",
            fmt="%s",
        )

    def generateConfigFile(self):
        # save the .conf file so we can process the updated network with the correct voltages
        print("allvoltageset")
        print(self.allvoltageset)
        # save config file
        # np.save(self.regionDataDir+'voltageSet',self.allvoltageset)
        voltagestring = ""
        for voltage in self.allvoltageset:
            if voltage >= self.lowestAllowedVoltage:
                voltagestring = voltagestring + "|" + str(int(np.floor(voltage * 1000)))

        alreadydownloaded = True

        if alreadydownloaded:
            firstcomment = ""
            secondcomment = ""
            thirdcomment = "#"
            fourthcomment = "#"
        else:
            firstcomment = "#"
            secondcomment = "#"
            thirdcomment = ""
            fourthcomment = ""

        if len(self.country) == 0:
            destdir = self.continent
            urlsegment = self.continent
            continentspecifier = "continent=" + self.continent
        else:
            destdir = self.continent + "/" + self.country
            urlsegment = self.continent + "/" + self.country
            continentspecifier = ""

        conftxt = """
### Database
dname='%s'
duser='postgres'
dpassword='OpenGridMap'

### Data Source and Destination
## Data Source
# specify either the location of dump and poly file or specify the download link
%sddump='%sddump.pbf'
%spfile='%spfile.poly'
%sddump_url='http://download.geofabrik.de/%s-latest.osm.pbf'
%spfile_url='http://download.geofabrik.de/%s.poly'
## Destination Directory
destdir='%s'

%s

### Transnet
vlevels='%s'
## Transnet arguments
# -t plot topology
# -v verbose logging
# -l load estimation (does only work for Germany, Austria and Switzerland)
# -e evaluation of point-to-point connections (only makes sense for Germany, since coverage of OSM data is sufficient high)
trans_args='-t'
		""" % (
            self.dbname,
            firstcomment,
            Params.downloadedTransnetDataDir + self.region + "/",
            secondcomment,
            Params.downloadedTransnetDataDir + self.region + "/",
            thirdcomment,
            urlsegment,
            fourthcomment,
            urlsegment,
            destdir,
            continentspecifier,
            voltagestring[1:],
        )

        file = open(Params.networkConfigDir + self.region + ".conf", "w")
        file.write(conftxt)
        file.close()

    # this allows us to run the generation code so we can generate any country or continent's high voltage power grid network from open data
    # runs th prepare_db.sh script in transnet's bash/ folder, which uses the .conf file we  generate with "generateConfigFile" function
    def initTransnet(self):
        # cdir=os.getcwd()

        # print('make sure to modify shell script with cache size of your computer, if you want to reimport')
        # print('According to params, prepare db shell script location is:')
        # print(Params.prepareDBScriptLoc+'prepare_db.sh')

        # #make sure to modify script with cache size of your computer!
        # os.chdir(prepareDBScriptLoc)
        # subprocess.call(['sh','prepare_db.sh',])
        # os.chdir(cdir)
        # quit()

        # you must prepare the database for the country or continent of interest.

        # once in transnet/bash, command should look like:
        # $ ./prepare_db.sh ../configs/great-britain.conf
        # if that succeeds, we run again to actually infer the network:
        # $ ./run_country.sh ../configs/great-britain.conf

        # transnet generates files in the transnet/models directory which include a map of the network and the cim.xml data, which we use to generate our version of the network in GEOMAGICA format.
        pass

    def calcGICs(self):
        # The GEOMAGICA formatted data is run by our own version of GEOMAGICA, called GIC_Model.py.
        gicModel = GIC_Model()
        connsAndNodesPath = self.processedNetworkDir + self.region
        for r in self.EfieldFiles.keys():
            efilePath = self.EfieldFiles[r]
            networkFileDir = self.processedNetworkDir + self.region
            # [gics,lineLengths]=gicModel.runModel(efilePath,networkFileDir)
            # savedata=[gics,lineLengths]
            # np.save('gics'+str(r)+'perYear',savedata,allow_pickle=True)
            savedata = np.load("gics" + str(r) + "perYear.npy", allow_pickle=True)
            [gics, lineLengths] = savedata
            print("number of lines")
            print(len(lineLengths))
            print("median line length")
            print(np.median(lineLengths))
            print("max line length")
            print(np.max(lineLengths))
            points = []
            lats = []
            longs = []
            # print(gics)
            print("r")
            print(r)
            # print(len(gics))
            # print(len(self.sortednodes))
            voltages = []
            loggics = []
            cutoffGIC = []
            countedgics = []
            # for i in range(0,self.sortednodes):
            # 	self.sortednodes[i]
            nodesarr = list(self.sortednodes.values())
            for i in range(0, len(nodesarr)):
                node = nodesarr[i]

                if node[3] != "substation":
                    continue
                # print('node')
                # print(node)
                index = node[0] - 1
                # point=Point([float(node[1]),float(node[2])])
                # points.append(point)

                voltage = node[-1]
                # print('voltage')
                # print(voltage)
                # print(gics[index])
                voltages.append(voltage)
                lats.append(float(node[1]))
                longs.append(float(node[2]))
                nodeGIC = gics[index]
                transformerGIC = nodeGIC / 1.3

                if transformerGIC == 0:
                    loggics.append(np.nan)
                else:
                    loggics.append(np.log(transformerGIC))
                if transformerGIC > 175:
                    cutoffGIC.append(175)
                else:
                    cutoffGIC.append(transformerGIC)
                countedgics.append(transformerGIC)
            # plt.plot(sorted(gics[index]))
            # plt.show()

            df = pd.DataFrame(
                {
                    "lats": lats,
                    "longs": longs,
                    "GIC": countedgics,
                    "logGIC": loggics,
                    "cutoffGIC": cutoffGIC,
                }
            )
            df.dropna(subset=["GIC"], inplace=True)
            print("number of substations")
            print(len(df))
            # Plotter.plotGICsBubble(df,self)
        # print(sortednodes)

        # now that's done, we have network files located in the directory of
        # GEOMAGICA.

        # sys.path.insert(1, os.path.join(sys.path[0], '../transnet/app/'))
        # import Transnet as Transnet

        # transnet=Transnet()
        # transnet.
        # run transnet on all the countries
        # print('allvoltageset')
        # print(self.allvoltageset)
