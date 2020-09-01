from PySide2.QtWidgets import QApplication, QMenu, QMessageBox
from PySide2.QtWidgets import QDialog, QPushButton,QGridLayout
from commands import Command
from iohandlers import IOHandler
import openmesh as om
import numpy as np
from signals import Signals
from geometry import Geometry
from hullultstrength import HullUltStrength
import os
from PySide2.QtCore import Slot

class HUSCommand(Command):
    def __init__(self):
        super().__init__()
        app = QApplication.instance()
        app.registerIOHandler(HullUltStrengthImporter())
        self.mainwin = app.mainFrame
        self.hf_prop = DialogHullFormModify(self.mainwin)
        self.hf_hydroCurves = DialogHullFormHydrostaticCurves(self.mainwin)
        self.hf=0
        self.si=0

        #tools = app.mainFrame.menuTools
        mb = app.mainFrame.menuBar()

        self.menuMain = QMenu("Ultimate strength")

        self.menuHullFormModify = self.menuMain.addAction("&Modify form")
        self.menuHullFormModify.triggered.connect(self.onModifyForm)

        self.menuResults = QMenu("&Results")
        self.menu2 = QMenu("&Menu2")
        self.menuMain.addMenu(self.menuResults)
        self.menuMain.addMenu(self.menu2)



        menuResulMomentCurvatureDiagram = self.menuResults.addAction("Moment- Curvature Diagram")
        menuResulMomentCurvatureDiagram.triggered.connect(self.onShowMomentCurvatureDiagram)



        #tools.addMenu(menu)
        mb.addMenu(self.menuMain)
        #self.menuMain.setEnabled(False)


        Signals.get().geometryAdded.connect(self.registerHullForm)
        Signals.get().selectionChanged.connect(self.registerSelection)

    @Slot()
    def registerHullForm(self, hullForm):
        if isinstance(hullForm, HullUltStrength):
            self.hf=hullForm
            self.menuMain.setEnabled(True)

    @Slot()
    def registerSelection(self, si):
        self.si=si

    def onModifyForm(self):
        if isinstance(self.hf, HullUltStrength):
            self.hf_prop.setCurrentHullForm(self.hf)
            self.hf_prop.exec()

    def onShowMomentCurvatureDiagram(self):
        if isinstance(self.hf,HullUltStrength):
            self.hf_hydroCurves.setCurrentHullForm(self.hf)
            self.hf_hydroCurves.exec()



class HullUltStrengthImporter(IOHandler):
    def __init__(self):
        super().__init__()
        Signals.get().importGeometry.connect(self.importGeometry)

    def importGeometry(self, fileName):
        if len(fileName) < 1:
            return
        filename, file_extension = os.path.splitext(fileName)
        if file_extension != ".hus":
            return
        hus=HullUltStrength(fileName)
        g = hus.readModel()
        Signals.get().geometryImported.emit(g)

    def getImportFormats(self):
        return [".hus"]


class DialogHullFormModify(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.mainwin = parent
        self.btnModify = self.createButton("&Modify", self.regenerateHullFormMesh)

        mainLayout = QGridLayout()
        mainLayout.addWidget(self.btnModify, 0, 0)
        self.setLayout(mainLayout)
        self.currentHullForm=0


    def createButton(self, text, member):
        button = QPushButton(text)
        button.clicked.connect(member)
        return button

    def regenerateHullFormMesh(self):
        self.currentHullForm.move(1, 0, 0)
        Signals.get().geometryRebuild.emit(self.currentHullForm)

    def setCurrentHullForm(self, currentHullForm):
        self.currentHullForm = currentHullForm
        self.setWindowTitle("Modify Hull Form")

class DialogHullFormHydrostaticCurves(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.mainwin = parent
        self.btnRefresh = self.createButton("&Refresh", self.refreshResults)

        mainLayout = QGridLayout()
        mainLayout.addWidget(self.btnRefresh, 0, 0)
        self.setLayout(mainLayout)
        self.currentHullForm=0


    def createButton(self, text, member):
        button = QPushButton(text)
        button.clicked.connect(member)
        return button

    def refreshResults(self):
        self.currentHullForm.getResults()
        pass

    def setCurrentHullForm(self, currentHullForm):
        self.currentHullForm = currentHullForm
        self.setWindowTitle("Hull Form Hydrostatic Curves")

def createCommand():
    return HUSCommand()
