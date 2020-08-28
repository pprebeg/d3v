from iohandlers import IOHandler
from signals import Signals
from geometry import Geometry
import openmesh as om
import os
import numpy as np
import csv
import math as Math


# import matplotlib.pyplot as plt

class HullFormMeshQuality:
    def __init__(self):
        self._numWL = 30
        self._numPnWLhalf = 30
        self._distPolyOrder=3

    @property
    def numPointsWLhalf(self):
        return self._numPnWLhalf

    def _getDistribution(self, maxv, minv, n, pot):
        x = [0.0] * n
        for i in range(n):
            fi = float(i)
            fn1 = float(n - 1)
            x[i] = fi ** pot / fn1 ** pot * (maxv - minv) + minv
        x.reverse()
        return x

    def genWLPositions(self,hWL_top,  hWL_bottom ):
        wlPos = self._getDistribution(hWL_top, hWL_bottom, self._numWL, self._distPolyOrder)
        return wlPos
    def genWLPositionsUsingObligatory(self,obligatoryLines:list ):
        testLines = self._getDistribution(obligatoryLines[0], obligatoryLines[-1], self._numWL, self._distPolyOrder)

        nol=len(obligatoryLines)
        wlPos=[]
        i1TL=0
        for iol in range(1,nol):
            hmax=obligatoryLines[iol-1]
            hmin = obligatoryLines[iol]
            numWL=0
            for iTL in range(i1TL,self._numWL):
                if testLines[iTL] < hmax:
                    if testLines[iTL] > hmin:
                        numWL =numWL + 1
                    else:
                        i1TL=iTL
                        break
            wlPosi = self._getDistribution(hmax, hmin, numWL+2, 1)
            for wl in wlPosi:
                if len[wlPos==0]:
                    wlPos.append(wl)
                elif wl < wlPos[-1]:
                    wlPos.append(wl)
        return wlPos


class HullForm(Geometry):
    def __init__(self, fileName):
        super().__init__()
        self.filename = fileName
        self.shipdata = {}
        self.pdecks =[]
        self.pbulkheads = []
        self.hfmq = HullFormMeshQuality()
        results = self.readShipData()
        self.shipdata = results[0]
        self.pdecks = results[1]
        self.pbulkheads = results[2]
        self.wlinesPos = []  # positive y waterlines
        self.wlinesNeg = []  # negative y waerlines
        self.wlKeel = []  # keel waterline (one waterline)

    def generateMesh(self):

        hmax=self.pdecks[0]
        #hmax=9.4
        wlPos=self.hfmq.genWLPositions(hmax, 0)
        #wlPos = self.hfmq.genWLPositionsUsingObligatory(self.pdecks)
        #results = self.hullGen(self.shipdata, wlPos, self.hfmq.numPointsWLhalf)
        lines = self.hullGen(self.shipdata, wlPos, self.hfmq.numPointsWLhalf)
        self.wlinesPos = lines[0]  # positive y waterlines
        self.wlinesNeg = lines[1]  # negative y waerlines
        self.wlKeel = lines[2]  # keel waterline (one waterline)
        self.mesh = self.genHullFormMeshPP(lines)
        pass

    def getResults(self):
        self.getVolume(9.4)
        pass

    def getVolume(self, hvl):
        mesh=self.mesh
        vol = 0
        h=hvl
        i=0

        for fh in mesh.faces():  # facet handle
            p = []
            for vh in mesh.fv(fh):  # vertex handle
                p.append(mesh.point(vh))
            # A
            Ax = p[0][0]
            Ay = p[0][1]
            Az = p[0][2]
            # B
            Bx = p[1][0]
            By = p[1][1]
            Bz = p[1][2]
            # C
            Cx = p[2][0]
            Cy = p[2][1]
            Cz = p[2][2]

            hA = h - Az
            hB = h - Bz
            hC = h - Cz
            area = 1 / 2 * (Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By))
            vol = vol + abs(area) * (hA + hB + hC) / 3.0
        print(vol)
        return vol


    # kod za formu
    def genHullFormMesh(self, lines: list):
        mesh = om.TriMesh()
        wlinesPos = lines[0]  # positive y waterlines
        wlinesNeg = lines[1]  # negative y waerlines
        wlKeel = lines[2]  # keel waterline (one waterline)
        n1 = np.array([0,0,0])
        m1 = np.array([0,0,0])
        m2 = np.array([0,0,0])
        pt1= np.array(3)
        pt2= np.array(3)
        n2 = np.array([0,0,0])
        n3 = np.array([0,0,0])
        n4 = np.array([0,0,0])

        for i in range(len(lines) - 2):  # lijevo desno kobilica
            for j in range(len(lines[i])):  # broji vodne linije
                for k in range(len(lines[i][j]) - 2):  # broj tocaka na vodnoj liniji
                    if j == len(lines[i]) - 1:          #kobilica
                        mpt1 = lines[i][j-1][k]
                        mpt2 = lines[i][j-1][k + 1]
                        mpt3 = lines[i][len(lines[i])-1][k]
                        mpt4 = lines[i][len(lines[i])-1][k + 1]

                        # volumenA = volumenA + Voltot(mpt1,mpt2,mpt3,mpt4)
                        mesh.add_face(mesh.add_vertex(mpt1), mesh.add_vertex(mpt2), mesh.add_vertex(mpt3))
                        mesh.add_face(mesh.add_vertex(mpt2), mesh.add_vertex(mpt4), mesh.add_vertex(mpt3))

                    if j != len(lines[i])-1:         #Sve ostale vodne linije
                        mpt1 = lines[i][j][k]
                        mpt2 = lines[i][j][k + 1]
                        mpt3 = lines[i][j+1][k]
                        mpt4 = lines[i][j+1][k + 1]
                        mesh.add_face(mesh.add_vertex(mpt1), mesh.add_vertex(mpt2), mesh.add_vertex(mpt3))
                        mesh.add_face(mesh.add_vertex(mpt2), mesh.add_vertex(mpt4), mesh.add_vertex(mpt3))


        for i in range(len(lines) - 2):  # lijevo desno kobilica
                for k in range(len(lines[i][0]) - 2):  # broj tocaka na vodnoj liniji
                        pt1 = lines[i][0][k]
                        pt2 = lines[i][0][k + 1]
                        m1[0] = pt1[0]
                        m1[2] = pt1[2]
                        m2[0] = pt2[0]

                        m2[2] = pt2[2]
                        mesh.add_face(mesh.add_vertex(pt1), mesh.add_vertex(m1), mesh.add_vertex(pt2))
                        mesh.add_face(mesh.add_vertex(pt2), mesh.add_vertex(m1), mesh.add_vertex(m2))

        return mesh

    def _genFaces(self,mesh:om.TriMesh,whs:list, doReverse:bool):
        nl=len(whs)
        npt=len(whs[0])
        for iL in range(1, nl):
            npt_iL = len(whs[iL])
            npt_iL_1 = len(whs[iL-1])
            dip=0
            if npt_iL > npt_iL_1:
                if doReverse:
                    mesh.add_face(whs[iL][0], whs[iL][1], whs[iL - 1][0])
                else:
                    mesh.add_face(whs[iL][1], whs[iL][0], whs[iL - 1][0])
                dip = 1
            for ipL_1 in range(1,npt_iL_1):
                ip = ipL_1+dip
                if doReverse:
                    mesh.add_face(whs[iL - 1][ipL_1 - 1], whs[iL][ip], whs[iL - 1][ipL_1])
                    mesh.add_face(whs[iL - 1][ipL_1 - 1], whs[iL][ip - 1], whs[iL][ip])
                else:
                    mesh.add_face(whs[iL - 1][ipL_1-1],   whs[iL - 1][ipL_1],whs[iL ][ip])
                    mesh.add_face(whs[iL - 1][ipL_1 - 1], whs[iL][ip],    whs[iL][ip-1])

    def genHullFormMeshPP(self, lines: list):
        mesh = om.TriMesh()
        wlinesPos = lines[0]  # positive y waterlines
        wlinesNeg = lines[1]  # negative y waerlines
        wlKeel = lines[2]  # keel waterline (one waterline)
        wlinesPos.reverse()
        wlinesNeg.reverse()

        whsPos = []
        whsNeg = []
        whsi = []
        whsPos.append(whsi)
        whsNeg.append(whsi)
        for p in wlKeel:
            whsi.append(mesh.add_vertex(p))


        for wl in wlinesPos:
            whsi = []
            whsPos.append(whsi)
            for p in wl:
                whsi.append(mesh.add_vertex(p))
        for wl in wlinesNeg:
            whsi = []
            whsNeg.append(whsi)
            for p in wl:
                whsi.append(mesh.add_vertex(p))

        self._genFaces(mesh,whsPos,True)
        self._genFaces(mesh, whsNeg,False)
        return mesh



    def hullGen(self, shipdata: dict, pdecks: list, nump):
        # gs is the grid size of a cell, in pixels
        # Reminder to make gridsize scaled to the screen width
        # Sets hullform data to slider values
        shipdata["loa_val"] = shipdata["loa_val"]
        shipdata["boa_val"] = shipdata["boa_val"]

        #
        midshipsM = shipdata["ms_val"]  # Constant m in JC equation
        bowRakeM = shipdata["bow_val"]  # Constant m in JC equation
        transomM = shipdata["tr_val"]  # Constant m in JC equation
        fwdDeckM = shipdata["deck_val"]  # Constant m in JC equation

        transomBeamMax = (shipdata["boa_val"] * shipdata["tb_val"]) / 2  # Transom half beam
        transomTip = shipdata["draft_val"] * shipdata["td_val"]
        ACU = shipdata["loa_val"] * shipdata["acu_val"]
        keelFwd = shipdata["loa_val"] * shipdata["kf_val"]
        slope = shipdata["sa_val"]

        midBeam = []  # Array with midships half beam per deck
        bowRake = []  # Array with bow rake per deck
        bowRakeS = []  # Array with bow rake per deck in superstructure
        TB = 0  # Transom half beam of a deck
        transomBeam = []  # Array with transom half beam per deck
        fwdDeckMArray = []  # Array with constants m in JC equation for deck outlines
        AE = 0  # Aft end of a deck
        aftEnd = []  # Array with aft end of each deck
        aftEndS = []  # Array with aft end of each deck in superstructure
        noseConeBaseRadius = []  # See excel tool
        ogiveRadius = []  # See excel tool
        pdecks2 = []  # Array with deck positions of hull decks
        pdecks3 = []  # Array with deck positions of superstructure decks

        for i in range(len(pdecks)):  # Assign values to variables above
            if pdecks[i] <= shipdata["draft_val"]:  # For each deck that is in the hull
                midBeam.append((Math.acosh(
                    (pdecks[i] / shipdata["draft_val"]) * (Math.cosh(midshipsM * Math.pi) - 1) + 1) / (
                                            midshipsM * Math.pi)) * (shipdata["boa_val"] / 2))
                bowRake.append((Math.acosh(
                    (pdecks[i] / shipdata["draft_val"]) * (Math.cosh(bowRakeM * Math.pi) - 1) + 1) / (
                                            bowRakeM * Math.pi)) * (shipdata["loa_val"] - keelFwd))
                if pdecks[i] > transomTip:
                    TB = ((Math.acosh(((pdecks[i] - transomTip) / (shipdata["draft_val"] - transomTip)) * (
                                Math.cosh(transomM * Math.pi) - 1) + 1) / (transomM * Math.pi)) * (transomBeamMax))

                else:
                    TB = 0

                transomBeam.append(TB)
                fwdDeckMArray.append(fwdDeckM * (pdecks[i] / (shipdata[
                    "draft_val"])) + 0.001)  # Changes constant m in JC equation to make deck outlines becomes slimmer with decreasing z position (see below)
                if (pdecks[i] >= transomTip):
                    AE = (shipdata["draft_val"] - pdecks[i]) * Math.tan(slope)

                else:
                    AE = (shipdata["draft_val"] - transomTip) * Math.tan(slope) + (transomTip - pdecks[i]) * (
                                (ACU - (shipdata["draft_val"] - transomTip) * Math.tan(slope)) / transomTip)

                aftEnd.append(AE)
                pdecks2.append(pdecks[i])

            else:  # For each deck that is in the superstructure
                aftEndS.append((pdecks[i] - shipdata["draft_val"]) * Math.tan(slope))
                bowRakeS.append(shipdata["loa_val"] - ((pdecks[i] - shipdata["draft_val"]) * Math.tan(slope)) - keelFwd)
                pdecks3.append(pdecks[i])

        for i in range(len(midBeam)):  # Assign values to variables above cont.
            noseConeBaseRadius.append(midBeam[i] - transomBeam[i])
            if noseConeBaseRadius[i] > 0:
                ogiveRadius.append(
                    (Math.pow(noseConeBaseRadius[i], 2) + Math.pow((shipdata["loa_val"] / 2) - aftEnd[i], 2)) / (
                                2 * noseConeBaseRadius[i]))

            else:
                ogiveRadius.append(0)

        deckOutlinesHull = []  # Array with hull deck outline x, y coordinates
        # Get y points for every x
        for idk in range(len(midBeam)):  # For each deck in hull
            deckOutlinesHull.append([])  # For each deck create array
            if pdecks2[idk] != 0:  # If not keel
                if transomBeam[idk] > 0:  # Add vertical hull line at transom
                    deckOutlinesHull[idk].append([aftEnd[idk], 0])
                kmin = aftEnd[idk]
                kmax = shipdata["loa_val"] / 2
                klist = np.linspace(kmin, kmax, nump)
                for xpt in klist:
                    deckOutlinesHull[idk].append([xpt, (
                                Math.sqrt(Math.pow(ogiveRadius[idk], 2) - Math.pow(xpt - shipdata["loa_val"] / 2, 2)) +
                                noseConeBaseRadius[idk] - ogiveRadius[idk] + transomBeam[idk])])

                kmin = shipdata["loa_val"] / 2
                kmax = keelFwd + bowRake[idk]
                klist = np.linspace(kmin, kmax, nump)
                for xpt in klist:
                    eqX = (xpt - shipdata["loa_val"] / 2) / (
                                keelFwd + bowRake[idk] - (shipdata["loa_val"] / 2))  # Value of x in JC equation
                    deckOutlinesHull[idk].append([xpt, (1 - ((Math.cosh(eqX * fwdDeckMArray[idk] * Math.pi) - 1) / (
                                Math.cosh(fwdDeckMArray[idk] * Math.pi) - 1))) * midBeam[idk]])


            else:  # If keel draw top
                kmin = aftEnd[idk]
                kmax = (keelFwd + bowRake[idk])
                klist = np.linspace(kmin, kmax, nump * 2)
                for xpt in klist:
                    deckOutlinesHull[idk].append([xpt, 0])  # Straight line

        deckOutlinesS = []  # Array with superstructure deck outline x, y coordinates
        tumblehome = []  # Superstructure tumblehome
        for n in range(len(aftEndS)):  # For each deck in superstructure
            deckOutlinesS.append([])  # For each deck create array
            tumblehome = (pdecks3[n] - shipdata["draft_val"]) * Math.tan(
                slope)  # Calculate tumblehome y offset to subtract below
            deckOutlinesS[n].append([aftEndS[n], 0])  # Add vertical hull line at transom

            kmin = aftEndS[n]
            kmax = shipdata["loa_val"] / 2
            klist = np.linspace(kmin, kmax, nump)
            for xpt in klist:
                deckOutlinesS[n].append([xpt, (
                            Math.sqrt(Math.pow(ogiveRadius[0], 2) - Math.pow(xpt - shipdata["loa_val"] / 2, 2)) +
                            noseConeBaseRadius[0] - ogiveRadius[0] + transomBeam[0] - tumblehome)])

            kmin = shipdata["loa_val"] / 2
            kmax = (keelFwd + bowRakeS[n])
            klist = np.linspace(kmin, kmax, nump)
            for xpt in klist:
                eqX = (xpt - shipdata["loa_val"] / 2) / (
                            keelFwd + bowRakeS[n] - (shipdata["loa_val"] / 2))  # Value of x in JC equation
                deckOutlinesS[n].append([xpt, (1 - ((Math.cosh(eqX * fwdDeckMArray[0] * Math.pi) - 1) / (
                            Math.cosh(fwdDeckMArray[0] * Math.pi) - 1))) * (midBeam[0] - tumblehome)])

        wlinesPos = []
        wlinesNeg = []
        wlKeel = []

        for ii in range(len(deckOutlinesS)):
            wlineP = list()
            wlineN = list()
            for item in deckOutlinesS[ii]:
                p = np.array([item[0], item[1], pdecks3[ii]])
                wlineP.append(p)
                p = np.array([item[0], -item[1], pdecks3[ii]])
                wlineN.append(p)
            wlinesPos.append(wlineP)
            wlinesNeg.append(wlineN)

        for ii in range(len(deckOutlinesHull)):

            if pdecks2[ii] != 0:
                wlineP = list()
                wlineN = list()
                for item in deckOutlinesHull[ii]:
                    p = np.array([item[0], item[1], pdecks2[ii]])
                    wlineP.append(p)
                    p = np.array([item[0], -item[1], pdecks2[ii]])
                    wlineN.append(p)
                wlinesPos.append(wlineP)
                wlinesNeg.append(wlineN)
            else:
                for item in deckOutlinesHull[ii]:
                    p = np.array([item[0], item[1], pdecks2[ii]])
                    wlKeel.append(p)

        return [wlinesPos, wlinesNeg, wlKeel]

    def readShipData(self):
        shipdata = {}
        pdecks = []
        pbulkheads = []
        with open(self.filename, newline='') as csvfile:
            f = csv.DictReader(csvfile)
            shipset = 0
            for row in f:  # there is only one row after header row!!!!!
                shipset = row

            shipdata["loa_val"] = float(shipset['LOA'])  # treba li učitavanje vrijednosti biti u petlji?
            shipdata["boa_val"] = float(shipset['BOA'])  # treba li učitavanje vrijednosti biti u petlji?
            shipdata['draft_val'] = float(shipset['D'])
            shipdata['shipname'] = shipset['Name']

            splitdata = str(shipset['HullData']).split(" ")
            shipdata["ms_val"] = float(splitdata[0])
            shipdata["ms_val"] = float(splitdata[0])
            shipdata["bow_val"] = float(splitdata[1])
            shipdata["tr_val"] = float(splitdata[2])
            shipdata["deck_val"] = float(splitdata[3])
            shipdata["tb_val"] = float(splitdata[4])
            shipdata["td_val"] = float(splitdata[5])
            shipdata["acu_val"] = float(splitdata[6])
            shipdata["kf_val"] = float(splitdata[7])
            shipdata["sa_val"] = float(splitdata[8])

            shipdata["sp_val"] = float(splitdata[9])
            shipdata["cwl_val"] = float(splitdata[10])
            shipdata["lcb_val"] = float(splitdata[11])
            shipdata["cb_val"] = float(splitdata[12])
            shipdata["mc_val"] = float(splitdata[13])
            shipdata["bb_val"] = float(splitdata[14])
            shipdata["tran_val"] = float(splitdata[15])
            shipdata["ab_val"] = float(splitdata[16])

            shipdata["lwl_val"] = float(splitdata[17])
            shipdata["bwl_val"] = float(splitdata[18])
            shipdata["tf_val"] = float(splitdata[19])
            shipdata["ta_val"] = float(splitdata[20])

            shipdata["app1"] = float(splitdata[21])
            shipdata["area_app1"] = float(splitdata[22])
            shipdata["app2"] = float(splitdata[23])
            shipdata["area_app2"] = float(splitdata[24])
            shipdata["area_app3"] = float(splitdata[25])
            shipdata["area_app4"] = float(splitdata[26])
            shipdata["app5"] = float(splitdata[27])
            shipdata["area_app5"] = float(splitdata[28])
            shipdata["area_app6"] = float(splitdata[29])
            shipdata["area_app7"] = float(splitdata[30])
            shipdata["app8"] = float(splitdata[31])
            shipdata["area_app8"] = float(splitdata[32])
            shipdata["area_app9"] = float(splitdata[33])
            shipdata["area_app10"] = float(splitdata[34])
            shipdata["area_app11"] = float(splitdata[35])

            shipdata["cg_val"] = float(splitdata[36])
            shipdata["heading_val"] = float(splitdata[37])
            shipdata["amplitude_val"] = float(splitdata[38])
            shipdata["roll_val"] = float(splitdata[39])
            shipdata["damping_val"] = float(splitdata[40])
            shipdata["plr_val"] = float(splitdata[41])
            shipdata["gmt_val"] = float(splitdata[42])

            draft = shipdata["draft_val"]
            splitdata = str(shipset['DeckPos']).split(" ")
            for dp in splitdata:
                pdecks.append(float(dp))
            splitdata = str(shipset['BHPos']).split(" ")
            for dp in splitdata:
                pbulkheads.append(float(dp))
        results = [shipdata,pdecks,pbulkheads]
        return  results
