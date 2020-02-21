from PySide2.QtWidgets import QApplication, QMenu, QMessageBox
from commands import Command
from iohandlers import IOHandler
import openmesh as om
import numpy as np
from signals import Signals
from geometry import Geometry
from oofem import OOFEM
import os


class OOFEMCommand(Command):
    def __init__(self):
        super().__init__()
        app = QApplication.instance()
        app.registerIOHandler(OOFEMImporter())

        #tools = app.mainFrame.menuTools
        mb = app.mainFrame.menuBar()

        self.menuOOFEM = QMenu("OOFEM")

        self.menuView = QMenu("&View")
        self.menuAnalysis = QMenu("&Analysis")
        self.menuOOFEM.addMenu(self.menuView)
        self.menuOOFEM.addMenu(self.menuAnalysis)

        self.menuResults = QMenu("&Results")
        self.menuView.addMenu(self.menuResults)
        self.menuResults.setEnabled(False)

        menuCSID = self.menuView.addAction("CSID")
        menuCSID.triggered.connect(self.onCSID)

        menuMatID = self.menuView.addAction("MatID")
        menuMatID.triggered.connect(self.onMatID)

        menuThID = self.menuView.addAction("Thickness")
        menuThID.triggered.connect(self.onThID)

        menuAnalyse = self.menuAnalysis.addAction("Analyse")
        menuAnalyse.triggered.connect(self.onAnalyse)

        menuStressX = self.menuResults.addAction("Sigma X")
        menuStressX.triggered.connect(self.onSx)

        menuStressY = self.menuResults.addAction("Sigma Y")
        menuStressY.triggered.connect(self.onSy)

        menuStressXY = self.menuResults.addAction("Tau XY")
        menuStressXY.triggered.connect(self.onTxy)

        menuStressVM = self.menuResults.addAction("Sigma VM")
        menuStressVM.triggered.connect(self.onSVM)

        menuDisp = self.menuResults.addAction("Tot. Disp.")
        menuDisp.triggered.connect(self.onDisp)

        menuDx = self.menuResults.addAction("Dx")
        menuDx.triggered.connect(self.onDx)

        menuDy = self.menuResults.addAction("Dy")
        menuDy.triggered.connect(self.onDy)

        menuDz = self.menuResults.addAction("Dz")
        menuDz.triggered.connect(self.onDz)

        #tools.addMenu(menu)
        mb.addMenu(self.menuOOFEM)
        self.menuOOFEM.setEnabled(False)
        Signals.get().geometryAdded.connect(self.registerOOFEM)
        self.oofem = 0

    def registerOOFEM(self, oofem):
        if isinstance(oofem,OOFEM):
            self.oofem=oofem
            self.menuOOFEM.setEnabled(True)

    def onSx(self):
        self.oofem.showElementStress(0)
        Signals.get().geometryAdded.emit(self.oofem)

    def onSy(self):
        self.oofem.showElementStress(1)
        Signals.get().geometryAdded.emit(self.oofem)

    def onTxy(self):
        self.oofem.showElementStress(2)
        Signals.get().geometryAdded.emit(self.oofem)


    def onSVM(self):
        self.oofem.showElementStress(10)
        Signals.get().geometryAdded.emit(self.oofem)


    def onDx(self):
        self.oofem.showNodeVertexColor(0)
        Signals.get().geometryAdded.emit(self.oofem)

    def onDy(self):
        self.oofem.showNodeVertexColor(1)
        Signals.get().geometryAdded.emit(self.oofem)

    def onDz(self):
        self.oofem.showNodeVertexColor(2)
        Signals.get().geometryAdded.emit(self.oofem)

    def onDisp(self):
        self.oofem.showNodeVertexColor(10)
        Signals.get().geometryAdded.emit(self.oofem)

    def onMatID(self):

        self.oofem.showMaterialID()
        Signals.get().geometryAdded.emit(self.oofem)

    def onCSID(self):

        self.oofem.showCrossSectID()
        Signals.get().geometryAdded.emit(self.oofem)
        pass

    def onThID(self):

        self.oofem.showThicknessID()
        Signals.get().geometryAdded.emit(self.oofem)
        pass

    def onAnalyse(self):
        self.oofem.ReadOutput()
        self.menuResults.setEnabled(True)


class OOFEMImporter(IOHandler):
    def __init__(self):
        super().__init__()
        Signals.get().importGeometry.connect(self.importGeometry)

    def importGeometry(self, fileName):
        if len(fileName) < 1:
            return
        filename, file_extension = os.path.splitext(fileName)
        if file_extension != ".in":
            return
        oofem=OOFEM(fileName)
        oofem.genMesh()
        Signals.get().geometryImported.emit(oofem)

    def getImportFormats(self):
        return []


def createCommand():
    return OOFEMCommand()
