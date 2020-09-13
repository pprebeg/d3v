from PySide2.QtWidgets import QApplication, QMenu, QMessageBox,QFormLayout,QWidget,QHeaderView, QSplitter
from PySide2.QtWidgets import QDialog, QPushButton,QGridLayout,QVBoxLayout,QHBoxLayout,QTableView,QTextEdit,QLabel
from commands import Command
from iohandlers import IOHandler
import openmesh as om
import numpy as np
from signals import Signals
from geometry import Geometry
from hullultstrength import HullUltStrength
import os
from PySide2.QtCore import Slot,Qt,SIGNAL
from PySide2.QtCore import QAbstractTableModel, QModelIndex, QRect
from PySide2.QtGui import QColor, QPainter
from PySide2.QtCharts import QtCharts

class HUSCommand(Command):
    def __init__(self):
        super().__init__()
        app = QApplication.instance()
        app.registerIOHandler(HullUltStrengthImporter())
        self.mainwin = app.mainFrame
        self.hf_prop = DialogHullFormModify(self.mainwin)
        self.hus_TableChart = DialogHullFormTableChart(self.mainwin)
        self.hus=0
        self.si=0

        #tools = app.mainFrame.menuTools
        mb = app.mainFrame.menuBar()

        self.menuMain = QMenu("Ultimate strength")

        self.menuHullFormModify = self.menuMain.addAction("&View Type")
        self.menuHullFormModify.triggered.connect(self.onModifyForm)

        self.menuResults = QMenu("&Results")
        self.menuMain.addMenu(self.menuResults)



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
            self.hus=hullForm
            self.menuMain.setEnabled(True)

    @Slot()
    def registerSelection(self, si):
        self.si=si

    def onModifyForm(self):
        if isinstance(self.hus, HullUltStrength):
            self.hf_prop.setCurrentHullForm(self.hus)
            self.hf_prop.exec()

    def onShowMomentCurvatureDiagram(self):
        if isinstance(self.hus, HullUltStrength):
            self.hus_TableChart.setCurrentHUS(self.hus)
            self.hus_TableChart.exec()



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
        Signals.get().geometryImported.emit(hus)

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

class CustomTableModel(QAbstractTableModel):
    def __init__(self):
        QAbstractTableModel.__init__(self)
        self.input_data = []
        self.input_names = []
        self.mapping = {}
        self.column_count = 0
        self.row_count = 0


    def setInputData(self,input_names, input_data):
        self.beginResetModel()
        self.input_data=input_data
        self.input_names= input_names
        self.column_count = len(self.input_names)
        self.row_count = len(self.input_data)
        self.endResetModel()



    def rowCount(self, parent=QModelIndex()):
        return self.row_count

    def columnCount(self, parent=QModelIndex()):
        return self.column_count

    def headerData(self, section, orientation, role):
        if role != Qt.DisplayRole:
            return None

        if orientation == Qt.Horizontal:
            return self.input_names[section]
        else:
            return "{}".format(section + 1)

    def data(self, index, role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            return self.input_data[index.row()][index.column()]
        elif role == Qt.EditRole:
            return self.input_data[index.row()][index.column()]
        elif role == Qt.BackgroundRole:
            for color, rect in self.mapping.items():
                if rect.contains(index.column(), index.row()):
                    return QColor(color)
            # cell not mapped return white color
            return QColor(Qt.white);
        return None

    def setData(self, index, value, role=Qt.EditRole):
        if index.isValid() and role == Qt.EditRole:
            self.input_data[index.row()][index.column()] = float(value)
            self.dataChanged.emit(index, index)
            return True
        return False

    def flags(self, index):
        return Qt.ItemIsEnabled | Qt.ItemIsEditable | Qt.ItemIsSelectable

    def add_mapping(self, color, area):
        self.mapping[color] = area

    def clear_mapping(self):
        self.mapping = {}



class DialogHullFormTableChart(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        sizetxt=25
        self.mainwin = parent


        self.model = CustomTableModel()
        self.table_view = QTableView()
        self.table_view.setModel(self.model)
        self.table_view.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table_view.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)

        self.chart = QtCharts.QChart()
        self.chart.setAnimationOptions(QtCharts.QChart.AllAnimations)


        self.chart_view = QtCharts.QChartView(self.chart)
        self.chart_view.setRenderHint(QPainter.Antialiasing)
        self.chart_view.setMinimumSize(400, 300)


        self.setWindowTitle("Hull Gitder Ultimate Strength")
        self.btnGenerate = self.createButton("&Generate", self.refreshResults)
        self.btnGenerate.setVisible(False)
        self.txtMultSagg = QTextEdit()
        self.txtMultSagg.setFixedHeight(sizetxt)
        self.txtMultSagg.setText('Postaviti granični moment')
        self.txtMultSagg.setAlignment(Qt.AlignRight)
        self.txtMultHogg = QTextEdit()
        self.txtMultHogg.setFixedHeight(sizetxt)
        self.txtMultHogg.setText('Postaviti granični moment')
        self.txtMultHogg.setAlignment(Qt.AlignRight)


        mainLayout = QVBoxLayout()
        mainLayout.setStretch(0,1)
        mainLayout.setStretch(1, 0)
        tablechartLayout = QGridLayout()
        tablechartSplitter = QSplitter()


        controlLayout = QHBoxLayout()
        controlWidget = QWidget()
        controlWidget.setFixedHeight(sizetxt*3)
        controlWidget.setLayout(controlLayout)

        inputLayout = QFormLayout()
        controlLayout.addLayout(inputLayout)
        #controlWidget.setVisible(False)
        controlLayout.addWidget(self.btnGenerate)

        inputLayout.addRow("Sagg Ultimate moment:", self.txtMultSagg)
        inputLayout.addRow("Hogg Ultimate moment", self.txtMultHogg)

        #tablechartLayout.addWidget(self.table_view, 0, 0)
        #tablechartLayout.addWidget(self.chart_view, 0, 1)
        tablechartSplitter.addWidget(self.table_view)
        tablechartSplitter.addWidget(self.chart_view)   
        #mainLayout.addLayout(tablechartLayout)
        mainLayout.addWidget(tablechartSplitter)
        mainLayout.addLayout(controlLayout)
        mainLayout.addWidget(controlWidget)

        self.setLayout(mainLayout)
        self.currentHUS=0


    def createButton(self, text, member):
        button = QPushButton(text)
        button.clicked.connect(member)
        return button

    def refreshResults(self):
        #self.currentHUS.getResults(9, 1.025)
        #return
        input_data = self.currentHUS.model_results['Lusa iteration results'].data
        input_names = self.currentHUS.model_results['Lusa iteration results'].column_names
        self.model.setInputData(input_names, input_data)
        self.chart.removeAllSeries()
        seriesColorHex = "#000000"
        ix=2
        iy=1

        series = QtCharts.QLineSeries()
        series.setName(input_names[iy]+ " " + input_names[ix])
        mapper = QtCharts.QVXYModelMapper(self)
        mapper.setXColumn(ix)
        mapper.setYColumn(iy)
        mapper.setSeries(series)
        mapper.setModel(self.model)
        self.chart.addSeries(series)
            # get the color of the series and use it for showing the mapped area
        seriesColorHex = "{}".format(series.pen().color().name())
        self.model.add_mapping(seriesColorHex, QRect(iy, 0, 1, self.model.rowCount()))
        self.model.add_mapping(seriesColorHex, QRect(ix, 0, 1, self.model.rowCount()))

        self.chart.createDefaultAxes()


    def setCurrentHUS(self, currentHUS):
        self.currentHUS = currentHUS
        self.setWindowTitle("Hull Ultimate Strength")
        self.refreshResults()

def createCommand():
    return HUSCommand()
