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

        menuCSID = self.menuView.addAction("CSID")
        menuCSID.triggered.connect(self.onCSID)

        menuMatID = self.menuView.addAction("MatID")
        menuMatID.triggered.connect(self.onMatID)

        menuAnalyse = self.menuAnalysis.addAction("Analyse")
        menuAnalyse.triggered.connect(self.onAnalyse)

        #tools.addMenu(menu)
        mb.addMenu(self.menuOOFEM)
        self.menuOOFEM.setEnabled(False)
        Signals.get().geometryAdded.connect(self.registerOOFEM)
        self.oofem = 0

    def registerOOFEM(self, oofem):
        if isinstance(oofem,OOFEM):
            self.oofem=oofem
            self.menuOOFEM.setEnabled(True)

    def onMatID(self):

        self.oofem.showMaterialID()
        Signals.get().geometryAdded.emit(self.oofem)

    def onCSID(self):

        self.oofem.showCrossSectID()
        Signals.get().geometryAdded.emit(self.oofem)
        pass

    def onAnalyse(self):
        QMessageBox.information(None, "Analysis", "Run OOFEM Analysis")

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
