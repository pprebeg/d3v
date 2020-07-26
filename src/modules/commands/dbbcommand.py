from PySide2.QtWidgets import QApplication, QMenu, QMessageBox
from PySide2.QtWidgets import QDialog, QPushButton,QGridLayout
from commands import Command
from iohandlers import IOHandler
import openmesh as om
import numpy as np
from signals import Signals
from geometry import Geometry
from dbb import DBBProblem,DBB
import os
from PySide2.QtCore import Slot

class DBBCommand(Command):
    def __init__(self):
        super().__init__()
        app = QApplication.instance()
        app.registerIOHandler(DBBImporter())
        self.mainwin = app.mainFrame
        self.dbbprop = DialogDBBProps(self.mainwin)
        self.dbbproblem=0
        self.si=0

        #tools = app.mainFrame.menuTools
        mb = app.mainFrame.menuBar()

        self.menuMain = QMenu("DBB")

        self.menuInitTestProblem = self.menuMain.addAction("TestProblem")
        self.menuInitTestProblem.triggered.connect(self.onGenerateTestProblerm)

        self.menuModifyBlock = QMenu("&Modify block")
        self.menu2 = QMenu("&Menu2")
        self.menuMain.addMenu(self.menuModifyBlock)
        self.menuMain.addMenu(self.menu2)



        menuMove = self.menuModifyBlock.addAction("Move")
        menuMove.triggered.connect(self.onMoveDBB)



        #tools.addMenu(menu)
        mb.addMenu(self.menuMain)
        #self.menuMain.setEnabled(False)


        Signals.get().geometryAdded.connect(self.registerDBB)
        Signals.get().selectionChanged.connect(self.registerSelection)
        self.dbb = 0

    @Slot()
    def registerDBB(self, dbbproblem):
        if isinstance(dbbproblem,DBBProblem):
            self.dbbproblem=dbbproblem
            self.menuMain.setEnabled(True)

    @Slot()
    def registerSelection(self, si):
        self.si=si

    def onGenerateTestProblerm(self):
        dbbproblem = DBBProblem("")
        Signals.get().geometryImported.emit(dbbproblem.hull)
        for deck in dbbproblem.decks:
            Signals.get().geometryImported.emit(deck)
        for dbb in dbbproblem.dbbs:
            Signals.get().geometryImported.emit(dbb)
        self.menuInitTestProblem.setEnabled(False)

    def onMoveDBB(self):
        if self.si.haveSelection():
            currDBB=self.si.getGeometry()
            if isinstance(currDBB,DBB):
                self.dbbprop.setCurrentDBB(currDBB)
                self.dbbprop.exec()



class DBBImporter(IOHandler):
    def __init__(self):
        super().__init__()
        Signals.get().importGeometry.connect(self.importGeometry)

    def importGeometry(self, fileName):
        if len(fileName) < 1:
            return
        filename, file_extension = os.path.splitext(fileName)
        if file_extension != ".dbb":
            return
        dbbproblem = DBBProblem(fileName)
        Signals.get().geometryImported.emit(dbbproblem.hull)
        for deck in dbbproblem.decks:
            Signals.get().geometryImported.emit(deck)
        for dbb in dbbproblem.dbbs:
            Signals.get().geometryImported.emit(dbb)

    def getImportFormats(self):
        return []

class DialogDBBProps(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.mainwin = parent
        self.btnMove = self.createButton("&Move", self.moveCurrentDBB)

        mainLayout = QGridLayout()
        mainLayout.addWidget(self.btnMove, 0, 0)
        self.setLayout(mainLayout)
        self.currentDBB=0


    def createButton(self, text, member):
        button = QPushButton(text)
        button.clicked.connect(member)
        return button

    def moveCurrentDBB(self):
        self.currentDBB.move(1, 0, 0)
        Signals.get().geometryRebuild.emit(self.currentDBB)

    def setCurrentDBB(self, currentDBB):
        self.currentDBB = currentDBB
        self.setWindowTitle("Move DBB")

def createCommand():
    return DBBCommand()
