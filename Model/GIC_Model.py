#!/usr/bin/env python
"""
********************************************************************
Basic GIC calculation method following Lehtinen and Pirjola 1985.
Adapted from ComputeGIC.m from C. Beggan's (BGS) adaptation of
K. Turnbull's Fortran code.

For usage to calculate expected GIC values in the grid from the 
Horton et al. (2012) paper with a 1 V/km geoelectric field execute:
    $ python GIC_Model_Horton.py

NOTES:
-   Due to differences in distance calculations and Python Numpy
    methods, final values are off the exact Horton values by
    a few percent in most cases.
-   Lehtinen-Pirjola method needs more null points added in to 
    equate to nodal admission matrix method, hence the representation
    of lines between HV, LV, ground and switch nodes.
-   Note that in this case the transformer resistance is explicitly 
    given as a connection between nodes.

Created by R Bailey (ZAMG, Austria) on 2015-08-03.
Minor modifications to allow for arbitrary geolocation and modified 
input/output format by Morgan Rivers on 2021-04-12.
#********************************************************************
"""

import os
import Params
import sys
import getopt
import numpy as np
from scipy import interpolate
from math import radians, tan, atan, atan2, cos, sin, acos, asin
from math import sqrt, pi, log

try:
    import IPython
except:
    pass
import matplotlib.pyplot as plt


class GIC_Model:
    def __init__(self):
        pass

    #####################################################################
    #                       FUNCTIONS                                   #
    #####################################################################

    def grc_azimuth(self, lonlat1, lonlat2):
        """
        Function to compute the geographic distance on a prolate ellipsiod such
        as a planet or moon. This computes the distance accurately - not that it
        makes much difference in most cases, as the location of the points of interest is
        generally relatively poorly known.

        This function is based on the formula from the Wikipedia page, verified
        against the Geoscience Australia website calculator

        Author: Ciaran Beggan
        Rewritten from matlab into python by R Bailey, ZAMG.

        Returns azimuth between two points.
        """

        a, b = 6378.137, 6356.752
        f = (a - b) / a

        u1 = atan((1.0 - f) * tan(pi / 180.0 * (lonlat1[1])))
        u2 = atan((1.0 - f) * tan(pi / 180.0 * (lonlat2[1])))

        L = pi / 180.0 * (lonlat2[0] - lonlat1[0])
        Lambda, converge, iters = L, False, 0

        while not converge and iters < 20:
            sinsig = sqrt(
                (cos(u2) * sin(Lambda)) ** 2.0
                + (cos(u1) * sin(u2) - sin(u1) * cos(u2) * cos(Lambda)) ** 2.0
            )
            cossig = sin(u1) * sin(u2) + cos(u1) * cos(u2) * cos(Lambda)
            sig = atan2(sinsig, cossig)

            sinalpha = (cos(u1) * cos(u2) * sin(Lambda)) / sinsig
            cossqalpha = 1.0 - sinalpha**2.0
            cos2sigm = cossig - (2.0 * sin(u1) * sin(u2)) / cossqalpha

            C = (f / 16.0) * cossqalpha * (4.0 + f * (4.0 - 3.0 * cossqalpha))

            calclambda = L + (1.0 - C) * f * sinalpha * (
                sig + C * sinalpha * (cos2sigm + C * cossig * (-1.0 + 2.0 * cos2sigm))
            )

            if abs(Lambda - calclambda) < 10.0 ** (-12.0):
                converge = True
                Lambda = calclambda
            else:
                iters = iters + 1
                Lambda = calclambda

        usq = cossqalpha * ((a**2.0 - b**2.0) / b**2.0)
        A = 1.0 + usq / 16384.0 * (
            4096.0 + usq * (-768.0 + usq * (320.0 - 175.0 * usq))
        )
        B = usq / 1024.0 * (256.0 + usq * (-128.0 + usq * (74.0 - 47.0 * usq)))
        delsig = (
            B
            * sinsig
            * (
                cos2sigm
                + 0.25
                * B
                * (
                    cossig * (-1.0 + 2.0 * cos2sigm)
                    - (1.0 / 6.0)
                    * B
                    * cos2sigm
                    * (-3.0 + 4.0 * sinalpha**2.0)
                    * (-3.0 + 4.0 * cos2sigm**2.0)
                )
            )
        )
        s = b * A * (sig - delsig)
        a1 = atan2(
            cos(u2) * sin(Lambda), cos(u1) * sin(u2) - sin(u1) * cos(u2) * cos(Lambda)
        )
        # a2 = atan2(cos(u1)*sin(Lambda), -sin(u1)*cos(u2) + cos(u1)*sin(u2)*cos(Lambda) )
        # print usq, A, B, delsig, s, a1

        if np.isnan(a1):
            a1 = 0.0

        # a1 = -a1
        if a1 < 0.0:
            a1 = 2.0 * pi + a1

        return a1

    def grc_distance(self, lat1, lon1, lat2, lon2, result):
        """
        Calculate the great circle distance between two points
        on the earth (specified in decimal degrees) using the Haversine method.
        Combination of:
        http://stackoverflow.com/questions/4913349/haversine-formula-in-python-bearing-and-distance-between-two-gps-points
        &
        http://gis.stackexchange.com/questions/29239/calculate-bearing-between-two-decimal-gps-coordinates
        """

        # Convert decimal degrees to radians
        lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

        # Haversine formula for distance:
        # Source: Wiki https://en.wikipedia.org/wiki/Haversine_formula
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        # c =  2.* asin( sqrt( sin((dlat)/2.)**2. + cos(lat1)*cos(lat2)* sin((dlon)/2.)**2. ) )
        r = 6371.0  # Radius of earth in kilometers. Use 3956 for miles

        if dlat == 0.0 and dlon == 0.0:
            return 0.0

        # Great circle distance:
        c = acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(dlon))

        if result == "km":
            return c * r
        elif result == "rad":
            return c

    #####################################################################
    #                       MAIN PROGRAM                                #
    #####################################################################
    def runModel(self, efilePath, connsAndNodesPath, isSimplistic):
        print("efilePath")
        print(efilePath)
        print("connsAndNodesPath")
        print(connsAndNodesPath)
        testfield = False
        # ===============================================================
        # 1) DEFINE NETWORK CONFIGURATION
        # ===============================================================

        # networkpath =       "network/HortonGrid.txt"
        # connectionspath =   "network/HortonGrid_Connections.txt"
        networkpath = connsAndNodesPath + "Nodes.txt"  # "network/nodes.txt"
        connectionspath = (
            connsAndNodesPath + "Connections.txt"
        )  # "network/connections.txt"
        # Read station and transformer data:
        network = open(networkpath, "r")
        netdata = network.readlines()

        # Define number of nodes:
        nnodes = len(netdata)

        print("nnodes")
        print(nnodes)
        # Read data into arrays:
        npf = np.float32
        geolat, geolon = np.zeros(nnodes, dtype=npf), np.zeros(nnodes, dtype=npf)
        country = []
        sitecode, sitename = [], []
        sitenum = np.zeros(nnodes, dtype=np.int32)
        res_earth = np.zeros(nnodes, dtype=npf)
        res_trans = np.zeros(nnodes, dtype=npf)
        # Station-to-index and index-to-station dicts:
        s2i, i2s = {}, {}

        for i in range(nnodes):
            # READ NODE DATA:
            # ---------------
            data = netdata[i].split("\t")
            # Dictionary for station to index:
            s2i[data[2]] = int(data[0])
            # Dictionary for index to station:
            i2s[int(data[0])] = data[2]
            # Site number:
            sitenum[i] = int(data[0])
            # Site readable names:
            sitename.append(data[1])
            # Site code names:
            sitecode.append(data[2])
            # Site country:
            country.append(data[3])
            # Geographic latitude of node:
            geolat[i] = float(data[4])
            # Geographic longitude of node:
            geolon[i] = float(data[5])
            # Earthing resistance of each node: #multiplied by 3 to give per phase values
            res_earth[i] = float(data[6]) * 3
            # Transformer resistances (see notes):
            res_trans[i] = float(data[7])

        # READ CONNECTION DATA:
        # ---------------------
        connectionsfile = open(connectionspath, "r")
        conndata = connectionsfile.readlines()

        # Number of connections:
        nconnects = len(conndata)
        print("nconnects")
        print(nconnects)
        nodefrom = np.zeros(nconnects, dtype=np.int32)
        nodeto = np.zeros(nconnects, dtype=np.int32)
        res_line = np.zeros(nconnects, dtype=npf)
        voltage_lines = np.zeros(nconnects, dtype=np.int32)

        for i in range(0, nconnects):
            conns = conndata[i].split("\t")
            # Connection starts at site:
            nodefrom[i] = int(s2i[conns[1]])
            # Connection ends at site:
            nodeto[i] = int(s2i[conns[2]])
            # Line resistance between connecting sites:
            # (Divide by 3 for full transformer representation as paper values are given per phase) #What?? Per phase would be leaving unchanged...
            # res_line[i] = float(conns[5])/3.
            res_line[i] = float(conns[5])
            # Voltage level of line:
            voltage_lines[i] = float(conns[6])

        # Set inf values to very high resistance for computation:
        valInf = 1.0e9
        res_trans[res_trans == np.inf] = valInf
        res_earth[res_earth == np.inf] = valInf
        res_line[res_line == np.inf] = valInf

        # ===============================================================
        # 2) DEFINE LOCATION VARIABLES
        # ===============================================================

        # Boundaries of North American geographic box for Horton grid:
        # nbound_NA = 35.0
        # ebound_NA = -80.0
        # sbound_NA = 32.0
        # wbound_NA = -88.0

        # luxembourg
        # nbound_NA = 6.6
        # ebound_NA = 50
        # sbound_NA = 5.8
        # wbound_NA = 49.4

        # #ireland
        # nbound_NA = 56
        # ebound_NA = -5
        # sbound_NA = 51
        # wbound_NA = -10

        # # Option 1: Basic 1 V/km geoelectric fields, use NA coordinates for Horton example:
        # if testfield:
        #     nbound, ebound, sbound, wbound = nbound_NA, ebound_NA, sbound_NA, wbound_NA
        #     # Spacing of cells in box (roughly equidistant):
        #     # x_inc, y_inc = 6./60., 8./60.
        #     # x_inc, y_inc =(50-49.4)/60, (6.6- 5.8)/60
        #     x_inc, y_inc =(5)/60, (5)/60

        #     # Latitude and longitude range within boundaries:
        #     lat = np.arange(sbound, nbound, x_inc)
        #     lon = np.arange(wbound, ebound, y_inc)

        # TODO: Add in method to remove station

        # ===============================================================
        # 3) ESTABLISH MATRICES FOR LP1985 METHOD
        # ===============================================================

        # Define resistance matrices:
        resis = np.zeros((nnodes, nnodes), dtype=npf)
        connections = np.zeros((nnodes, nnodes), dtype=npf)
        for i in range(0, nconnects):
            x, y = int(nodefrom[i]), int(nodeto[i])
            # print('connection index')
            # print(i)
            # print('node from')
            # print(x)
            # print('node to')
            # print(y)

            if resis[x, y] > 0.0 and res_line[i] > 0.0:
                # If res already added to this line, add another line in parallel:
                resis[x, y] = 1.0 / (1.0 / resis[x, y] + 1.0 / res_line[i])
                resis[y, x] = resis[x, y]
            else:
                resis[x, y] = res_line[i]
                resis[y, x] = resis[x, y]
            if resis[x, y] > 0.0:
                connections[x, y] = 1.0 / resis[x, y]
                connections[y, x] = 1.0 / resis[x, y]

        # Calculate matrix of distance between each point:
        dists = np.zeros((nnodes, nnodes), dtype=npf)
        azi = np.zeros((nnodes, nnodes), dtype=npf)
        for i in range(0, nnodes):
            for j in range(0, nnodes):
                lati, loni = geolat[i], geolon[i]
                latj, lonj = geolat[j], geolon[j]
                if resis[i, j] != 0.0:
                    dists[i, j] = self.grc_distance(lati, loni, latj, lonj, "km")
                    dists[j, i] = dists[i, j]
                    try:
                        azi[i, j] = self.grc_azimuth([loni, lati], [lonj, latj])
                        if azi[i, j] >= 0.0 and azi[i, j] <= np.pi:
                            azi[j, i] = azi[i, j] + np.pi
                        if azi[i, j] > np.pi and azi[i, j] <= 2 * np.pi:
                            azi[j, i] = azi[i, j] - np.pi
                    except ZeroDivisionError:
                        azi[i, j] = np.nan
                        azi[j, i] = np.nan
        # Create earth impedance matrix:
        # LP1984 eq. (3): **Z**
        # (assuming that nodes are spaced far enough apart to not affect one another)
        if isSimplistic:
            earthimp = np.diag(res_earth)
        else:
            earthimp = np.diag(res_earth + res_trans)

        # Calculate network admittance matrix:
        # LP1984 eq. (10): **Y**
        netadmit = -1.0 * connections + np.diag(sum(connections))

        # Create system matrix (1+YZ), which needs to be inverted:
        # LP1984 eq. (part of 12): **1 + YZ**
        systemmat = np.dot(netadmit, earthimp) + np.identity(nnodes)

        # ===============================================================
        # 4) DEFINE GEOELECTRIC FIELD
        # ===============================================================

        layermodel = None
        surfacemodel = None

        # Option 1: Use 1 V/km field:
        # N1Vkm = True
        # E1Vkm = True
        # if testfield:
        #     en, ee = np.zeros((len(lat),len(lon)), dtype=npf), np.zeros((len(lat),len(lon)), dtype=npf)
        #     if N1Vkm:
        #         en.fill(1.)
        #     if E1Vkm:
        #         ee.fill(0)
        #     en_int = interpolate.interp2d(lon, lat, en)
        #     ee_int = interpolate.interp2d(lon, lat, ee)

        # # Option 2 (-e): Use external field file:
        # else:
        darray = np.loadtxt(efilePath, skiprows=0)
        # Reshape into usable array:
        lat = sorted(list(set(darray[:, 0])))
        lon = sorted(list(set(darray[:, 1])))
        nlat, nlon = len(lat), len(lon)
        ien, iee = 2, 5  # Indices of N and E field component (real parts only here)

        endata = np.reshape(darray[:, ien], (nlat, nlon), order="F")
        eedata = np.reshape(darray[:, iee], (nlat, nlon), order="F")

        # Interpolate E-field onto these points:
        en_int = interpolate.interp2d(lon, lat, endata)
        ee_int = interpolate.interp2d(lon, lat, eedata)

        # Can do quick check to make sure this matches up (direct comparison to random file lines):
        # print("Comparison to lines from file for column #s {} and {}".format(ien+1, iee+1))
        # print("            Lat        Lon         EN(R)      EN(C)      EN_tot    EE(R)      EE(C)       EE_tot")
        # print("Line #863:  48.0190    11.7933     0.0509     0.0650     0.0826    -0.0531    -0.0675     0.0859")
        # print("--> {:.4f} (N), {:.4f} (E)".format(en_int(11.7933, 48.0190)[0], ee_int(11.7933, 48.0190)[0]))
        # print("Line #2547: 48.3770    17.1800     0.1056     0.1363     0.1724    -0.0729    -0.0944     0.1192")
        # print("--> {:.4f} (N), {:.4f} (E)".format(en_int(17.1800, 48.3770)[0], ee_int(17.1800, 48.3770)[0]))

        # ===============================================================
        # 5) INTEGRATE FIELD ALONG LINES
        # ===============================================================

        # Number of steps in path integration (more is slower but more exact):
        steps = 200
        pathlatsteps, pathlonsteps = np.zeros((nconnects, steps), dtype=npf), np.zeros(
            (nconnects, steps), dtype=npf
        )
        Vn_tot, Ve_tot = np.zeros(nconnects, dtype=npf), np.zeros(nconnects, dtype=npf)
        E_e, E_n = np.zeros((nconnects, steps), dtype=npf), np.zeros(
            (nconnects, steps), dtype=npf
        )
        line_dists = []
        for i in range(0, nconnects):
            slat = geolat[nodefrom[i]]
            slon = geolon[nodefrom[i]]

            flat = geolat[nodeto[i]]
            flon = geolon[nodeto[i]]

            # Set number of steps in path integration:
            steplat = (flat - slat) / float(steps)
            steplon = (flon - slon) / float(steps)

            if steplat != 0.0:
                pathlatsteps[i, :] = np.arange(slat, flat - steplat * 0.1, steplat)
            else:
                pathlatsteps[i, :] = slat

            if steplon != 0.0:
                pathlonsteps[i, :] = np.arange(slon, flon - steplon * 0.1, steplon)
            else:
                pathlonsteps[i, :] = slon

            for p in range(0, len(pathlonsteps[i, :])):
                E_n[i, p] = en_int(pathlonsteps[i, p], pathlatsteps[i, p])
                E_e[i, p] = ee_int(pathlonsteps[i, p], pathlatsteps[i, p])

            # Integrate to get V = int(E*dL), use cylindrical coordinates:
            intline = dists[nodefrom[i], nodeto[i]]

            # we would like to return all the lines that are not substation grounding connections (which are at 0 volts)
            if voltage_lines[i] > 0:
                line_dists.append(intline)

            intazi = azi[nodefrom[i], nodeto[i]]
            for j in range(steps - 1):
                # For a North field:
                vnseg = (
                    0.5 * (E_n[i, j] + E_n[i, j + 1]) * cos(intazi) * (intline / steps)
                )
                # For an East field:
                veseg = (
                    0.5 * (E_e[i, j] + E_e[i, j + 1]) * sin(intazi) * (intline / steps)
                )

                Vn_tot[i] += vnseg
                Ve_tot[i] += veseg
        # ===============================================================
        # 6) USE V AND SYSTEM MATRIX TO DETERMINE CURRENT (J)
        # ===============================================================

        # Create Nvoltage matrix for northward field:
        Nvoltage = np.zeros((nnodes, nnodes), dtype=npf)
        for l in range(nconnects):
            Nvoltage[nodefrom[l], nodeto[l]] = Vn_tot[l]
            Nvoltage[nodeto[l], nodefrom[l]] = -Vn_tot[l]

        # Create Evoltage matrix for eastward field:
        Evoltage = np.zeros((nnodes, nnodes), dtype=npf)
        for l in range(nconnects):
            Evoltage[nodefrom[l], nodeto[l]] = Ve_tot[l]
            Evoltage[nodeto[l], nodefrom[l]] = -Ve_tot[l]

        # Use Ohm's law to calculate the current along each pathlength:
        # LP1985 eq. (14): J = V/R
        Nlinecurr = np.zeros((nnodes, nnodes), dtype=npf)
        for m in range(nnodes):
            for n in range(nnodes):
                if resis[m, n] > 0.0:
                    newval = Nvoltage[m, n] / resis[m, n]
                    if not np.isnan(newval):
                        Nlinecurr[m, n] = newval
                    else:
                        Nlinecurr[m, n] = 0.0

        Elinecurr = np.zeros((nnodes, nnodes), dtype=npf)

        for m in range(nnodes):
            for n in range(nnodes):
                if resis[m, n] > 0.0:
                    newval = Evoltage[m, n] / resis[m, n]
                    if not np.isnan(newval):
                        Elinecurr[m, n] = newval
                    else:
                        Elinecurr[m, n] = 0.0

        # Total line current from both components:
        Tlinecurr = Nlinecurr + Elinecurr

        # Total current at each node due to the E field:
        Nsourcevec = np.sum(Nlinecurr, axis=0)
        Esourcevec = np.sum(Elinecurr, axis=0)

        # GIC formed at each node due to the E field:
        # LP1985 eq. (part of 12): **(1 + YZ)^(-1)*J**
        netconN, resid, rank, s = np.linalg.lstsq(systemmat, np.transpose(Nsourcevec))
        netconE, resid, rank, s = np.linalg.lstsq(systemmat, np.transpose(Esourcevec))
        netconT = netconN + netconE

        # ****************************************************************
        results = [netconN, netconE]  # GIC PER NODE
        # ****************************************************************

        # Calculate current flowing through transformers (NOT line current):
        Vn_trans = np.multiply(netconN, (res_earth + res_trans))
        Ve_trans = np.multiply(netconE, (res_earth + res_trans))
        Iline_N, Iline_E = np.zeros(nconnects), np.zeros(nconnects)
        for l in range(0, nconnects):
            Iline_N[l] = (Vn_trans[nodefrom[l]] - Vn_trans[nodeto[l]]) / res_line[l]
            Iline_E[l] = (Ve_trans[nodefrom[l]] - Ve_trans[nodeto[l]]) / res_line[l]

        # ****************************************************************
        # this gives us the currents from the lines, but the currents through the nodes will give a closer approximation
        transresults = [Iline_N, Iline_E]  # GIC PER TRANSFORMER
        # ****************************************************************
        gic = []
        for i in range(nnodes):
            # print("{}  \t{:.2f}    \t{:.2f}".format(sitename[i], results[0][i], results[1][i]))
            gic.append(np.abs(results[0][i] + results[1][i]))
            # gic.append(np.sqrt(transresults[0][i]**2+transresults[1][i]**2))
        return [np.array(gic), line_dists]

    def printResults(self):
        # ===============================================================
        # 7) PRINT RESULTS
        # ===============================================================

        print("")
        print("-----------------------------------------------------------------")
        if testfield:
            print(
                "RESULTS FOR 1 V/KM NORTHWARD E field"
            )  # AND EASTWARD ELECTRIC FIELDS")
        else:
            print("RESULTS FOR COMPLEX GEOELECTRIC FIELD")
        print("-----------------------------------------------------------------")

        # Print GIC per node (only grounding currents are non-zero):
        print("")
        # print("Network node\tGIC_N    \tGIC_E")

        gic = []
        for i in range(nnodes):
            # print("{}  \t{:.2f}    \t{:.2f}".format(sitename[i], results[0][i], results[1][i]))
            gic.append(np.sqrt(results[0][i] ** 2 + results[1][i] ** 2))

        plt.plot(gic)

        plt.show()
        plt.plot(np.sort(np.array(gic)))
        plt.title("1 V/km field East, Ireland, 1 ohm grounding")
        plt.xlabel("Substation")
        plt.ylable("GIC into/out of substation (A)")
        plt.show()

        # For convenience, convert connections to transformer codes:
        transID = {
            "Sub1_LV-Sub1_E": "T1",
            "Sub2_LV-Sub2_E": "T3 /T4",
            "Sub3_HV-Sub3_LV": "T5 / T15 (series)",
            "Sub3_LV-Sub3_E": "T5 / T15 (common)",
            "Sub4_LV-Sub4_E": "T2 / T13 (LV), T12 / T14 (common)",
            "Sub4_HV-Sub4_LV": "T12 / T14 (series)",
            "Sub4_HV-Sub4_E": "T2 / T13 (HV)",
            "Sub5_LV-Sub5_G": "T8",
            "Sub5_HV-Sub5_G": "T9",
            "Sub6_HV-Sub6_E": "T6 / T7",
            "Sub8_HV-Sub8_E": "T10 / T11",
        }

        # Print GIC per transformer (transformers here being 'lines'):
        print("")
        maxslen = max(len(p) for p in transID.values())
        print("Transformer code\t\t\tGIC_N\tGIC_E")
        for i in range(nconnects):
            codes, codes_r = (
                sitename[nodefrom[i]] + "-" + sitename[nodeto[i]],
                sitename[nodeto[i]] + "-" + sitename[nodefrom[i]],
            )
            if codes in transID or codes_r in transID:
                try:
                    print(
                        (
                            "{}".ljust(1 + maxslen - len(transID[codes]))
                            + "\t{:.2f}\t{:.2f}"
                        ).format(transID[codes], Iline_N[i] / 3.0, Iline_E[i] / 3.0)
                    )
                except:
                    print(
                        (
                            "{}".ljust(1 + maxslen - len(transID[codes_r]))
                            + "\t{:.2f}\t{:.2f}"
                        ).format(transID[codes_r], Iline_N[i] / 3.0, Iline_E[i] / 3.0)
                    )

        # IPython.embed()
