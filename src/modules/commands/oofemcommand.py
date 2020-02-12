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

        menu = QMenu("OOFEM")

        menuOpen = menu.addAction("&Open")
        menuAnalyse = menu.addAction("&Analyse")

        menuOpen.triggered.connect(self.onOpen)
        menuAnalyse.triggered.connect(self.onAnalyse)

        #tools.addMenu(menu)
        mb.addMenu(menu)

    def onOpen(self):
        mesh = om.TriMesh()
        g = Geometry()
        g.mesh = mesh
        Signals.get().geometryImported.emit(g)

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
