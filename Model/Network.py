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
from geopandas.tools import sjoin
from geovoronoi import voronoi_regions_from_coords, points_to_coords
from shapely.ops import cascaded_union


class Network:
    def __init__(self):
        # assumptions:
        # each transformer gets its own winding resistance, and it's assumed that transformers are all connected to the same grounding grid at a substation, with a given resistance per phase of the windings (transformertogroundR). These are all connected together, and there's an earthing resistance (groundingR) for the grounding grid to earth.
        self.isSimplistic = True
        self.lowestAllowedVoltage = 10  # 200#100#kV
        self.groundingR = ".5"  # ohms
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

        self.windingkV = [230, 345, 500, 735]
        self.windingConductance = 1 / np.array([0.692, 0.356, 0.195, 0.159])  # Ohms
        self.transformertogroundR = {}
        self.transformertogroundR[275] = 1 / np.interp(
            275, self.windingkV, self.windingConductance
        )
        self.transformertogroundR[400] = 1 / np.interp(
            400, self.windingkV, self.windingConductance
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
        if not country:
            self.region = continent
            self.regionDataDir = Params.planetNetworkDir + continent + "/"
            self.downloadedDataDir = Params.planetTransnetDataDir + continent + "/"
        else:
            self.region = country
            self.regionDataDir = (
                Params.countryNetworkDir + continent + "/" + country + "/"
            )
            self.downloadedDataDir = (
                Params.countryTransnetDataDir + continent + "/" + country + "/"
            )
        self.continent = continent
        # self.regionDataDir='europe/luxembourg/'

        self.processedNetworkDir = Params.networkSaveDir

        # self.country='luxembourg'
        self.country = country
        self.loadPolygon()
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

        if self.isSimplistic:
            nodesDict = self.importNodesSimplistic()
        else:
            nodesDict = self.importNodes()

        # sort nodesdict by index
        self.sortednodes = [
            v
            for k, v in sorted(
                self.nodesDict.items(), key=lambda item: item[1]["index"]
            )
        ]

        if self.isSimplistic:
            self.importConnections(nodesDict)
        else:
            self.importConnectionsSimplistic(nodesDict)
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
        Plotter.plotNetwork(self.voltages, self.lines, ax, self.region)
        plt.show()
        print("network plotted")
        return [self.voltages, self.lines]

    def loadPolygon(self):
        pfile = self.downloadedDataDir + "pfile.poly"
        if not self.country:
            self.boundaryPolygon = PolyParser.poly_to_polygon(pfile)

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

    def importNodesSimplistic(self):
        tree = ET.parse(self.regionDataDir + "cim.xml")
        root = tree.getroot()

        conns = []
        output = []
        windings = []
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
            # every winding has a connection to the ground node of the substation. All such windings can be combined. The highest winding voltage is saved as the voltage class of the transformer at that node.

            # node 5 is resistance directly to ground from this node

            if not (subname in nodesDict.keys()):
                # data format for ground node:
                # __0__ ____1____ ___2___ ___3___ _4_ __5__ _______6______ _7_
                # index node name  code   country lat long  gnd resistance  0
                savenode = np.array(
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

                nodesDict[subname] = {
                    "savenode": savenode,
                    "index": index,
                    "lat": lat,
                    "long": long,
                    "nodeType": "substation",
                    "substationName": subname,
                    "voltageClass": voltage,
                }
                index = index + 1

            if windname in nodesDict.keys():
                continue

            savenode = np.array(
                [index, windname, index, self.region, lat, long, "Inf", "Inf"]
            )

            # we assume the voltage class of the substation is the highest voltage of the transformer winding at the substation
            if voltage > nodesDict[subname]["voltageClass"]:
                nodesDict[subname]["voltageClass"] = voltage

            # subname node points to the ground node. Do not include a savenode.
            nodesDict[windname] = {
                "index": nodesDict[subname]["index"],
                "lat": lat,
                "long": long,
                "nodeType": "winding",
                "substationName": subname,
                "voltageClass": voltage,
            }  # ground has voltage 0

        self.nodesDict = nodesDict
        return nodesDict

    def importNodes(self):
        tree = ET.parse(self.regionDataDir + "cim.xml")
        root = tree.getroot()

        conns = []
        output = []
        windings = []
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

            # every substation has one ground node
            # node 5 is resistance directly to ground from this node
            # node 6 is resistance directly to transformer from this node (I think, not really sure about this one)
            if not (subname in nodesDict.keys()):
                # data format for simplistic node:
                # __0__ ____1____ ___2___ ___3___ _4_ __5__ _______6______ _7_
                # index node name  code   country lat long  gnd resistance  0
                savenode = np.array(
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

                nodesDict[subname] = {
                    "savenode": savenode,
                    "index": index,
                    "lat": lat,
                    "long": long,
                    "nodeType": "substation",
                    "voltageClass": 0,
                }  # ground has voltage 0
                index = index + 1

            if windname in nodesDict.keys():
                continue
            # data format for transformer node:
            # __0__ ____1____ ___2___ ___3___ _4_ __5__ _______6______ _7_
            # index node name  code   country lat long  gnd resistance Inf

            savenode = np.array(
                [index, windname, index, self.region, lat, long, "Inf", "Inf"]
            )

            # as it stands, there can only be one winding per voltage network for the substation. This may need to be adjusted if there were more complicated representations of networks.
            #
            # assert(not (windname in nodesDict.keys()))

            nodesDict[windname] = {
                "savenode": savenode,
                "index": index,
                "lat": lat,
                "long": long,
                "nodeType": "winding",
                "substationName": subname,
                "voltageClass": voltage,
            }  # ground has voltage 0

            index = index + 1

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
                assert nodesDict[winding1]["nodeType"] == "winding"
                # get the index for each winding
                index1 = nodesDict[winding1]["index"]
                subname1 = nodesDict[winding1]["substationName"]
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
                assert nodesDict[winding2]["nodeType"] == "winding"
                # get the index for each winding
                index2 = nodesDict[winding2]["index"]
                subname2 = nodesDict[winding2]["substationName"]
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

            r = resistancePerKilometer * (
                c.length / 1000
            )  # Ohms (net resistance of the power line, length converted to kilometers)

            # data format for connection to other transformer node:
            # __0__ ____1____ ___2___ _3_ _4_ _____5____ ____6____
            # index nodefrom  nodeto   0   0  resistance  voltage
            outputConnections.append([index, index1, index2, "0", "0", str(r), str(v)])
            index = index + 1

            if not (index1 in isConnectedToGround.keys()):
                assert nodesDict[subname1]["nodeType"] == "substation"
                sub1index = nodesDict[subname1]["index"]

                # data format for connection transformer winding node to ground node:
                # __0__ ____1____ ___2___ _3_ _4_ ________5_______ _6_
                # index nodefrom  nodeto   0   0  wind resistance   0
                outputConnections.append(
                    [index, index1, sub1index, "0", "0", windR, "0"]
                )
                index = index + 1

            # if we haven't connected winding2 to the substation ground yet
            if not (index2 in isConnectedToGround.keys()):
                sub2index = nodesDict[subname2]["index"]

                # data format for connection transformer winding node to ground node:
                # __0__ ____1____ ___2___ _3_ _4_ ________5_______ _6_
                # index nodefrom  nodeto   0   0  wind resistance   0
                outputConnections.append(
                    [index, index2, sub2index, "0", "0", windR, "0"]
                )
                index = index + 1

            # now we've connected the nodes connected to this power line to substation ground
            isConnectedToGround[index1] = True
            isConnectedToGround[index2] = True

            # for plotting
            lat1 = nodesDict[key1]["lat"]
            long1 = nodesDict[key1]["long"]
            lat2 = nodesDict[key2]["lat"]
            long2 = nodesDict[key2]["long"]
            p1 = Point(float(long1), float(lat1))
            p2 = Point(float(long2), float(lat2))
            line = LineString([[p1.x, p1.y], [p2.x, p2.y]])
            vindex = np.where(np.array(voltages) == v)[0]
            if len(vindex) == 0:
                voltages.append(v)
                lines.append([])
                i = len(voltages) - 1
            else:
                i = vindex[0]
            lines[i].append(line)

        self.connections = outputConnections
        self.lines = lines
        self.voltages = voltages
        self.nodesDict = nodesDict

    # each voltage line naively connects to the other through no resistance
    # the combination of substation grounding and transformer winding is assumed to be 0.5 ohms
    def importConnectionsSimplistic(self, nodesDict):
        conns = []

        # NOTE: again because of the bug with PyCIM, am hacking together a reference by using the names of the AC power lines, rather than the TransformerWinding(A) - Terminal - ConnectivityNode - Terminal - ACLineSegment -Terminal - ConnectivityNode - Terminal - TransformerWinding(B) UUID terminal chain from the cim.xml.
        outputConnections = []

        gdfs = []  # for plotting
        voltages = []
        lines = []
        startpoints = []
        stoppoints = []
        index = 0
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
                assert nodesDict[winding1]["nodeType"] == "winding"
                index1 = nodesDict[winding1]["index"]
                subname1 = nodesDict[winding1]["substationName"]
                key1 = subname1  # connects directly to the substation
            else:
                continue
                # there isn't a transformer for this connection. It must connect directly to a generating unit (aka power station).

            if winding2 in nodesDict.keys():
                assert nodesDict[winding2]["nodeType"] == "winding"
                index2 = nodesDict[winding2]["index"]
                subname2 = nodesDict[winding2]["substationName"]
                key2 = subname2  # connects directly to the sublstation
            else:
                continue
                # #there isn't a transformer for this connection. It must connect directly to a generating unit (aka power station).
            self.allvoltageset.add(v)

            if not (v in self.lineResistance.keys()):
                self.lineResistance[v] = 1 / np.interp(
                    v, self.linekV, self.lineConductance
                )
            resistancePerKilometer = self.lineResistance[v]

            r = resistancePerKilometer * (
                c.length / 1000
            )  # Ohms (net resistance of the power line, length converted to kilometers)

            # data format for connection to other transformer node:
            # __0__ ____1____ ___2___ _3_ _4_ _____5____ ____6____
            # index nodefrom  nodeto   0   0  resistance  voltage
            outputConnections.append([index, index1, index2, "0", "0", str(r), str(v)])
            index = index + 1

            # now we've connected the nodes connected to this power line to substation ground

            # for plotting
            lat1 = nodesDict[key1]["lat"]
            long1 = nodesDict[key1]["long"]
            lat2 = nodesDict[key2]["lat"]
            long2 = nodesDict[key2]["long"]
            p1 = Point(float(long1), float(lat1))
            p2 = Point(float(long2), float(lat2))
            line = LineString([[p1.x, p1.y], [p2.x, p2.y]])
            vindex = np.where(np.array(voltages) == v)[0]
            if len(vindex) == 0:
                voltages.append(v)
                lines.append([])
                i = len(voltages) - 1
            else:
                i = vindex[0]
            lines[i].append(line)

        self.connections = outputConnections
        self.lines = lines
        self.voltages = voltages

    # def analyzeNetwork(self):
    # 	print('OOOO')
    # 	lines=[]
    # 	points=[]

    # 	newnodes={}
    # 	newconnections=[]
    # 	# sortednodes={v[0]:v for k, v in sorted(self.nodesDict.items(), key=lambda item: item[1][0])}
    # 	# print(sortednodes)
    # 	for nk in self.nodesDict.keys():
    # 		node=self.nodesDict[nk]
    # 		point=Point([float(node[1]),float(node[2])])
    # 		points.append(point)

    # 		# print(node)
    # 		if(float(node[2])<51.1 and float(node[2])>50.2 and float(node[1])<-3 and float(node[1])>-5):
    # 			newnodes[nk]=node
    # 			# print(node)
    # 			for c in self.connections:
    # 				if(c[1]==node[0] or c[2]==node[0]):
    # 					add=True

    # 					# print(c)
    # 					for nc in newconnections:
    # 						if(nc[0]==c[0]):
    # 							add=False
    # 							break
    # 					if(add):
    # 						newconnections.append(c)
    # 						# print(c[1])
    # 						# print('c[2]')
    # 						# print(c[2])
    # 						# print('sortednodes[c[1]]')
    # 						# print(sortednodes[c[1]])
    # 						# print('sortednodes[c[2]]')
    # 						# print(sortednodes[c[2]])

    # 						p1=Point(float(sortednodes[c[1]][2]),float(sortednodes[c[1]][1]))
    # 						p2=Point(float(sortednodes[c[2]][2]),float(sortednodes[c[2]][1]))
    # 						line=LineString([[p1.x, p1.y], [p2.x, p2.y]])
    # 						lines.append(line)
    # 	# quit()
    # 	print('newnodes')
    # 	# print(newnodes)
    # 	print('newconnections')
    # 	# print(newconnections)
    # 	print('OOOO')

    # 	network=gpd.GeoDataFrame({'geometry':np.array(lines)})
    # 	gdf=gpd.GeoDataFrame(geometry=lines)

    # 	gdf.plot()
    # 	plt.show()
    # quit()

    def saveTextFiles(self):
        # save the text files which will be processed by GIC_Model
        savenodes = [
            x["savenode"] for x in self.sortednodes if ("savenode" in x.keys())
        ]
        # savenodes=[x['savenode'] for x in self.sortednodes if (('savenode' in x.keys()) and (x['substationName'] in ['SS_516651650', 'SS_197846408', 'SS_180637986', 'SS_88462768', 'SS_254158424', 'SS_264275258', 'SS_183749097', 'SS_456199341', 'SS_264852019', 'SS_194575030', 'SS_255668124', 'SS_88144450', 'SS_266025855', 'SS_35150701', 'SS_186889669', 'SS_104388595', 'SS_207339053', 'SS_207331199', 'SS_42998577', 'SS_179121616', 'SS_194614470', 'SS_185902319', 'SS_191715059', 'SS_79803767', 'SS_150124101', 'SS_147525695', 'SS_185693432', 'SS_89760546', 'SS_288283900', 'SS_207728854', 'SS_27144164', 'SS_96278380', 'SS_104388599', 'SS_361261945', 'SS_208124159', 'SS_137050528', 'SS_166764802', 'SS_148266260', 'SS_137050527', 'SS_230134599', 'SS_142185552', 'SS_402783212', 'SS_227767019', 'SS_227767022', 'SS_339008021', 'SS_260836695']))]

        # savenodeindices=[x['index'] for x in self.sortednodes if (x['substationName'] in ['SS_516651650', 'SS_197846408', 'SS_180637986', 'SS_88462768', 'SS_254158424', 'SS_264275258', 'SS_183749097', 'SS_456199341', 'SS_264852019', 'SS_194575030', 'SS_255668124', 'SS_88144450', 'SS_266025855', 'SS_35150701', 'SS_186889669', 'SS_104388595', 'SS_207339053', 'SS_207331199', 'SS_42998577', 'SS_179121616', 'SS_194614470', 'SS_185902319', 'SS_191715059', 'SS_79803767', 'SS_150124101', 'SS_147525695', 'SS_185693432', 'SS_89760546', 'SS_288283900', 'SS_207728854', 'SS_27144164', 'SS_96278380', 'SS_104388599', 'SS_361261945', 'SS_208124159', 'SS_137050528', 'SS_166764802', 'SS_148266260', 'SS_137050527', 'SS_230134599', 'SS_142185552', 'SS_402783212', 'SS_227767019', 'SS_227767022', 'SS_339008021', 'SS_260836695'])]
        # print(len(savenodeindices))
        print(self.processedNetworkDir + self.region + "Nodes.txt")
        savetxt(
            self.processedNetworkDir + self.region + "Nodes.txt",
            np.array(savenodes),
            delimiter="\t",
            fmt="%s",
        )
        # print('self.connections')
        # print(self.connections)
        # onlyirelandconnections=[x for x in self.connections if ((x[1] in savenodeindices) and (x[2] in savenodeindices))]
        # print(len(onlyirelandconnections))
        # print(len(self.connections))
        savetxt(
            self.processedNetworkDir + self.region + "Connections.txt",
            np.array(self.connections),
            delimiter="\t",
            fmt="%s",
        )
        # savetxt(self.processedNetworkDir+self.region+'Connections.txt',np.array(onlyirelandconnections),delimiter='\t', fmt='%s')

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

        if not self.country:
            destdir = self.continent
            urlsegment = self.continent
            continentspecifier = "continent=" + self.continent
            tags = "-t -g"
        else:
            destdir = self.continent + "/" + self.country
            urlsegment = self.continent + "/" + self.country
            continentspecifier = ""
            tags = "-t"

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
trans_args='%s'
		""" % (
            self.dbname,
            firstcomment,
            Params.configTransnetDataDir + destdir + "/",
            secondcomment,
            Params.configTransnetDataDir + destdir + "/",
            thirdcomment,
            urlsegment,
            fourthcomment,
            urlsegment,
            destdir,
            continentspecifier,
            voltagestring[1:],
            tags,
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
        print("calcGICs")
        gicModel = GIC_Model()

        loadDataFromFiles = False

        for r in self.EfieldFiles.keys():
            efilePath = self.EfieldFiles[r]
            networkFileDir = self.processedNetworkDir + self.region
            if loadDataFromFiles:
                savedata = np.save(
                    Params.networkAnalysisDir
                    + self.region
                    + "/gics"
                    + str(r)
                    + "perYear.npy",
                    allow_pickle=True,
                )
                [gics, lineLengths] = savedata
            else:
                [gics, lineLengths] = gicModel.runModel(
                    efilePath, networkFileDir, self.isSimplistic
                )
                savedata = [gics, lineLengths]
                print("saved")
                print("len(gics)")
                # print(len(gics))
                np.save(
                    Params.networkAnalysisDir
                    + self.region
                    + "/gics"
                    + str(r)
                    + "perYear",
                    savedata,
                    allow_pickle=True,
                )
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
            # nodesarr=list(self.sortednodes.values())
            for i in range(0, len(self.sortednodes)):
                node = self.sortednodes[i]

                if node["nodeType"] != "substation":
                    continue
                # print('node')
                # print(node)
                # point=Point([float(node[1]),float(node[2])])
                # points.append(point)

                voltage = node["voltageClass"]
                # print('voltage')
                # print(voltage)
                # print(gics[index])
                voltages.append(voltage)
                lats.append(float(node["lat"]))
                longs.append(float(node["long"]))
                index = node["index"]
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
            # Plotter.plotGICsBubble(df,self,r)
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

    # https://towardsdatascience.com/how-to-create-voronoi-regions-with-geospatial-data-in-python-adbb6c5f2134
    def calcVoronoi(self, allnodes, region, rateperyears):
        # fig, ax = plt.subplots(figsize=(12, 10))
        # world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
        # pp=gplt.polyplot(world,zorder=1)

        # plt.show()
        # stationShapes=[]
        # sRegionFailureProbs=np.array([])
        sRegionFailureProbs = pd.DataFrame()
        [nodeArrList, polyList] = self.getPolysWithNodes(allnodes)
        prevprobfails = []
        # print('len(polyList)')
        # print(len(polyList))
        for i in range(0, len(polyList)):
            poly = polyList[i]
            nodes = nodeArrList[i]

            # print('poly')
            # print(poly)
            # print('nodes')
            # print(nodes)
            shape = gpd.GeoDataFrame({"geometry": [poly]}, crs="epsg:3857")
            boundary = shape.to_crs(epsg=3857)
            gdf_proj = nodes.to_crs(boundary.crs)
            boundary_shape = cascaded_union(boundary.geometry)
            coords = points_to_coords(gdf_proj.geometry)
            if len(coords) < 2:  # voronoi can't handle only one node.
                continue
            # Calculate Voronoi Regions (appears to scramble ID of polygons and nodes, which is fixed below)
            poly_shapes, poly_to_pt_assignment = voronoi_regions_from_coords(
                coords, boundary_shape
            )
            # for r in rateperyears:
            # 	Plotter.plotVoronoiRegions(nodes[str(r)],nodes['geometry'], boundary_shape, poly_shapes, poly_to_pt_assignment,region, r)

            # probfails=pd.DataFrame()
            nnodes = len(poly_to_pt_assignment)
            cols = [str(x) for x in rateperyears]
            cols.append("geometry")
            probfails = gpd.GeoDataFrame(
                0,
                index=len(prevprobfails)
                + np.linspace(0, nnodes - 1, nnodes).astype(int),
                columns=cols,
                crs="epsg:3857",
            )

            # we need to assign the polygon shapes to a dataframe with the appropriate rateperyear probability of failure of that region
            for key, item in poly_to_pt_assignment.items():
                for r in rateperyears:
                    probfails[str(r)].iloc[key] = nodes[str(r)].iloc[item[0]]
                probfails["geometry"].iloc[key] = poly_shapes[key]
                # NOTE: assertion below fails somewhere in Europe, not sure why.
                # assert(probfails['geometry'].iloc[key].contains(nodes['geometry'].iloc[item[0]]))
            if len(prevprobfails) > 0:
                newdataframe = pd.concat([prevprobfails, probfails])
            else:
                newdataframe = probfails
            prevprobfails = newdataframe
        return newdataframe

    # return a list of polygons (independent nations or sections of nations) and associated nodes which lie within those polygons
    # https://gis.stackexchange.com/questions/208546/check-if-a-point-falls-within-a-multipolygon-with-python
    def getPolysWithNodes(self, nodes):
        # get a set of country land boundaries
        world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))

        points = gpd.GeoSeries(nodes["geometry"])

        # get a list of polygons
        allpolys = []
        for index, row in world.iterrows():
            if row["geometry"].geom_type == "MultiPolygon":
                for polygon in row["geometry"]:
                    if points.within(polygon).any():
                        allpolys.append(polygon)
            else:
                if points.within(row["geometry"]).any():
                    allpolys.append(row["geometry"])
        # print(len(allpolys))
        crs = {"init": "epsg:3857"}
        # create the list of points by splitting the dataframe
        splitpoints = []
        for poly in allpolys:
            geo_df = gpd.GeoDataFrame({"geometry": [poly]}, crs="epsg:3857")
            pointInPolys = sjoin(nodes, geo_df, how="left")
            grouped = pointInPolys.groupby("index_right")
            # print(len(grouped))
            # print(grouped.head())
            splitpoints.append(grouped.get_group("index_right" == 0))
        return [splitpoints, allpolys]
