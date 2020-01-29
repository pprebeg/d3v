from iohandlers import IOHandler
from signals import Signals
from geometry import Geometry
import openmesh as om
import os
import numpy as np
import csv
import math as Math
#import matplotlib.pyplot as plt
class HullFormImporter(IOHandler):
    def __init__(self):
        super().__init__()
        Signals.get().importGeometry.connect(self.importGeometry)

    def importGeometry(self, fileName):
        if len(fileName) < 1:
            return
        filename, file_extension = os.path.splitext(fileName)
        if file_extension != ".huf":
            return
        dbb=HullForm(fileName)
        g = Geometry()
        g.mesh = dbb.getmesh()
        Signals.get().geometryImported.emit(g)

    def getImportFormats(self):
        return []


def createIOHandler():
    return HullFormImporter()

class HullForm ():
    def __init__(self,fileName):
        self.filename = fileName
        self.testcalc()
    def getmesh(self):
        m=self.test()
        m=self.onCreateBox()
        return m
    def test(self):
        mesh= om.TriMesh()
        vhandle = []
        data = np.array([0, 1, 0])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([1, 0, 0])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([2, 1, 0])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([0, -1, 0])
        vhandle.append( mesh.add_vertex(data))
        data = np.array([2, -1, 0])
        vhandle.append( mesh.add_vertex(data))

        fh0 = mesh.add_face(vhandle[0], vhandle[1], vhandle[2])
        fh1 = mesh.add_face(vhandle[1], vhandle[3], vhandle[4])
        #fh2 = mesh.add_face(vhandle[0], vhandle[3], vhandle[1])

        vh_list = [vhandle[2], vhandle[1], vhandle[4]]
        fh3 = mesh.add_face(vh_list)

        return mesh
        pass
    def hullformmesh(self):
        mesh= om.TriMesh()
        #read self.filename
        return mesh
        pass
    def onCreateBox(self):
        mesh = om.TriMesh()
        # m --> min, M --> max
        # yapf: disable
        p0 = mesh.add_vertex([-1, -1, -1])
        p1 = mesh.add_vertex([-1, -1,  1])
        p2 = mesh.add_vertex([-1,  1, -1])
        p3 = mesh.add_vertex([-1,  1,  1])
        p4 = mesh.add_vertex([ 1, -1, -1])
        p5 = mesh.add_vertex([ 1, -1,  1])
        p6 = mesh.add_vertex([ 1,  1, -1])
        p7 = mesh.add_vertex([ 1,  1,  1])
        # yapf: enable

        mesh.add_face([p0, p6, p4])
        mesh.add_face([p0, p2, p6])
        mesh.add_face([p0, p4, p5])
        mesh.add_face([p0, p5, p1])
        mesh.add_face([p0, p3, p2])
        mesh.add_face([p0, p1, p3])
        mesh.add_face([p6, p2, p3])
        mesh.add_face([p6, p3, p7])
        mesh.add_face([p4, p7, p5])
        mesh.add_face([p4, p6, p7])
        mesh.add_face([p1, p5, p7])
        mesh.add_face([p1, p7, p3])

        return  mesh
    def testcalc(self):
        with open(self.filename, newline='') as csvfile:
            hfr = csv.reader(csvfile, delimiter='\t', quotechar='|')
            for row in hfr:
                rown=[]
                for x in row:
                    rown.append(float(x))
                print(row)
                print(rown)
                #print(', '.join(row))
        pass

    def hullGen(self,shipdata: dict, pdecks: list):
        # gs is the grid size of a cell, in pixels
        # Reminder to make gridsize scaled to the screen width
        gs = 5

        # gsm is the cell size in meters
        gsm = 0.5
        # Sets hullform data to slider values
        shipdata["loa_val"] = round(shipdata["loa_val"] / gsm) * gsm
        shipdata["boa_val"] = round(shipdata["boa_val"] / gsm) * gsm;

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
        for j in range(len(midBeam)):  # For each deck in hull
            deckOutlinesHull.append([])  # For each deck create array
            if pdecks2[j] != 0:  # If not keel
                if transomBeam[j] > 0:  # Add vertical hull line at transom
                    deckOutlinesHull[j].append([Math.ceil((aftEnd[j]) / gsm) * gsm * (gs / gsm), 0])

                # k=Math.ceil((aftEnd[j])/gsm)*gsm; k<Math.round((shipset[0].LOA/2)/gsm)*gsm; k+=gsm
                kmin = Math.ceil((aftEnd[j]) / gsm) * gsm
                kmax = round((shipdata["loa_val"] / 2) / gsm) * gsm
                step = gsm
                klist = np.arange(kmin, kmax, step)
                for k in klist:  # k=Math.ceil((aftEnd[j])/gsm)*gsm; k<round((shipdata["loa_val"]/2)/gsm)*gsm; k+=gsm:     #For aft half of each deck in the hull
                    deckOutlinesHull[j].append([k * (gs / gsm), (
                                Math.sqrt(Math.pow(ogiveRadius[j], 2) - Math.pow(k - shipdata["loa_val"] / 2, 2)) +
                                noseConeBaseRadius[j] - ogiveRadius[j] + transomBeam[j]) * (gs / gsm)])

                kmin = round((shipdata["loa_val"] / 2) / gsm) * gsm
                kmax = (keelFwd + bowRake[j])
                step = gsm
                klist = np.arange(kmin, kmax, step)
                for k in klist:  # (k=round((shipdata["loa_val"]/2)/gsm)*gsm; k<=(keelFwd + bowRake[j]); k+=gsm):  #For forward half of each deck in the hull
                    eqX = (k - shipdata["loa_val"] / 2) / (
                                keelFwd + bowRake[j] - (shipdata["loa_val"] / 2))  # Value of x in JC equation
                    deckOutlinesHull[j].append([k * (gs / gsm), (1 - (
                                (Math.cosh(eqX * fwdDeckMArray[j] * Math.pi) - 1) / (
                                    Math.cosh(fwdDeckMArray[j] * Math.pi) - 1))) * midBeam[j] * (gs / gsm)])


            else:  # If keel draw top
                kmin = Math.ceil((aftEnd[j]) / gsm) * gsm
                kmax = (keelFwd + bowRake[j])
                step = gsm
                klist = np.arange(kmin, kmax, step)
                for k in klist:  # k=Math.ceil((aftEnd[j])/gsm)*gsm; k<=(keelFwd + bowRake[j]) k+=gsm:
                    deckOutlinesHull[j].append([k * (gs / gsm), 0])  # Straight line

        deckOutlinesS = []  # Array with superstructure deck outline x, y coordinates
        tumblehome = []  # Superstructure tumblehome
        for n in range(len(aftEndS)):  # For each deck in superstructure
            deckOutlinesS.append([])  # For each deck create array
            tumblehome = (pdecks3[n] - shipdata["draft_val"]) * Math.tan(
                slope)  # Calculate tumblehome y offset to subtract below
            deckOutlinesS[n].append(
                [Math.ceil((aftEndS[n]) / gsm) * gsm * (gs / gsm), 0])  # Add vertical hull line at transom

            kmin = Math.ceil((aftEndS[n]) / gsm) * gsm
            kmax = round((shipdata["loa_val"] / 2) / gsm) * gsm
            step = gsm
            klist = np.arange(kmin, kmax, step)
            for k in klist:  # (k=Math.ceil((aftEndS[n])/gsm)*gsm; k<round((shipdata["loa_val"]/2)/gsm)*gsm; k+=gsm): #For aft half of each deck in the superstructure (same equation as above with tumblehome y offset subtracted)
                deckOutlinesS[n].append([k * (gs / gsm), (
                            Math.sqrt(Math.pow(ogiveRadius[0], 2) - Math.pow(k - shipdata["loa_val"] / 2, 2)) +
                            noseConeBaseRadius[0] - ogiveRadius[0] + transomBeam[0] - tumblehome) * (gs / gsm)])

            kmin = round((shipdata["loa_val"] / 2) / gsm) * gsm
            kmax = (keelFwd + bowRakeS[n])
            step = gsm
            klist = np.arange(kmin, kmax, step)
            for k in klist:  # (k=round((shipdata["loa_val"]/2)/gsm)*gsm; k<=(keelFwd + bowRakeS[n]); k+=gsm): #For forward half of each deck in the superstructure (same equation as above with tumblehome y offset subtracted)
                eqX = (k - shipdata["loa_val"] / 2) / (
                            keelFwd + bowRakeS[n] - (shipdata["loa_val"] / 2))  # Value of x in JC equation
                deckOutlinesS[n].append([k * (gs / gsm), (1 - ((Math.cosh(eqX * fwdDeckMArray[0] * Math.pi) - 1) / (
                            Math.cosh(fwdDeckMArray[0] * Math.pi) - 1))) * (midBeam[0] - tumblehome) * (gs / gsm)])

        wlinesPos = []
        wlinesNeg = []
        wlKeel = []

        for ii in range(len(deckOutlinesHull)):
            wlineP = list()
            wlineN = list()
            if pdecks2[ii] != 0:
                for item in deckOutlinesHull[ii]:
                    p = np.array([item[0], item[1], pdecks2[ii]])
                    wlineP.append(p)
                    p = np.array([item[0], -item[1], pdecks2[ii]])
                    wlineN.append(p)
            else:
                for item in deckOutlinesHull[ii]:
                    p = np.array([item[0], item[1], pdecks2[ii]])
                    wlKeel.append(p)
            wlinesPos.append(wlineP)
            wlinesNeg.append(wlineN)
        for ii in range(len(deckOutlinesS)):
            wlineP = list()
            wlineN = list()
            for item in deckOutlinesS[ii]:
                p = np.array([item[0], item[1], pdecks2[ii]])
                wlineP.append(p)
                p = np.array([item[0], -item[1], pdecks2[ii]])
            wlinesPos.append(wlineP)
            wlinesNeg.append(wlineN)

        return [wlinesPos,wlinesNeg,wlKeel]

    def readShipData(self):
        shipdata = {}
        pdecks = []
        pbulkheads = []
        with open('shipd.csv', newline='') as csvfile:
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

            splitdata = str(shipset['DeckPos']).split(" ")
            for dp in splitdata:
                pdecks.append(float(dp))
            splitdata = str(shipset['BHPos']).split(" ")
            for dp in splitdata:
                pbulkheads.append(float(dp))
        vlines = getDistribution(pdecks[0], pdecks[len(pdecks) - 1], 10, 3)
        # vlines= np.linspace(pdecks[0],pdecks[len(pdecks)-1],10)

        # hullgen.hullGen(shipdata,pdecks)
        self.hullGen(shipdata,vlines)


