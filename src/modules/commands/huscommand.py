from PySide2.QtWidgets import QApplication, QMenu, QMessageBox,QFormLayout,QWidget,QHeaderView, QSplitter,QLineEdit
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
from PySide2.QtCore import QObject
class HUSCommand(Command):
    def __init__(self):
        super().__init__()
        app = QApplication.instance()
        app.registerIOHandler(HullUltStrengthImporter())
        self.mainwin = app.mainFrame
        self.hus_meshcontrol = DialogMeshControl(self.mainwin)
        self.hus_TableChart = DialogHullFormTableChart(self.mainwin)
        self.hus=0
        self.si=0

        #tools = app.mainFrame.menuTools
        mb = app.mainFrame.menuBar()

        self.menuMain = QMenu("Ultimate strength")

        self.menuViewType = QMenu("&View Type")
        self.menuMain.addMenu(self.menuViewType)
        self.menuResults = QMenu("&Results")
        self.menuMain.addMenu(self.menuResults)



        menuResulMomentCurvatureDiagram = self.menuResults.addAction("Moment- Curvature Diagram")
        menuResulMomentCurvatureDiagram.triggered.connect(self.onShowMomentCurvatureDiagram)

        menuMeshControl = self.menuMain.addAction("View Control")
        menuMeshControl.triggered.connect(self.onMeshControl)



        #tools.addMenu(menu)
        mb.addMenu(self.menuMain)
        #self.menuMain.setEnabled(False)


        Signals.get().geometryAdded.connect(self.registerHullForm)
        Signals.get().selectionChanged.connect(self.registerSelection)

    @Slot()
    def registerHullForm(self, hus):
        if isinstance(hus, HullUltStrength):
            self.hus=hus
            self.menuMain.setEnabled(True)
            self.addMenus()

    def addMenus(self):
        for key,atr in self.hus.attrib_val_functions.items():
            menuResul = self.menuViewType.addAction(key)
            menuResul.triggered.connect(self.onColorControlMenu)
        for key, res in self.hus.element_results.items():
            menuResul = self.menuResults.addAction(key)
            menuResul.triggered.connect(self.onColorControlMenu)
            pass
    def onColorControlMenu(self):
        action = QObject.sender(self)
        key=action.iconText()
        self.hus.prepareModelForVisualization(key)
        self.hus_meshcontrol.setCurrentHUS(self.hus)
        Signals.get().geometryRebuild.emit(self.hus)
        pass

    @Slot()
    def registerSelection(self, si):
        self.si=si

    def onMeshControl(self):
        if isinstance(self.hus, HullUltStrength):
            self.hus_meshcontrol.setCurrentHUS(self.hus)
            self.hus_meshcontrol.exec()

    def onShowMomentCurvatureDiagram(self):
        if isinstance(self.hus, HullUltStrength):
            self.hus_TableChart.setCurrentHUS(self.hus)
            self.hus_TableChart.exec()



class HullUltStrengthImporter(IOHandler):
    def __init__(self):
        super().__init__()

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


class DialogMeshControl(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.mainwin = parent


        self.mainLayout = QVBoxLayout()

        self.setLayout(self.mainLayout)
        self.currentHUS=0
        self.mc=0
        self.txtMCupptresh=0
        self.txtMClowtresh = 0


    def createButton(self, text, member):
        button = QPushButton(text)
        button.clicked.connect(member)
        return button

    def applyMeshControl(self):

        self.mc.lowertreshold=float(self.txtMClowtresh.text())
        self.mc.uppertreshold = float(self.txtMCupptresh.text())
        self.currentHUS.regenerateusingcolor()
        Signals.get().geometryRebuild.emit(self.currentHUS)

    def setCurrentHUS(self, currentHUS):
        self.currentHUS = currentHUS
        self.mc = currentHUS.mc
        if len(self.windowTitle()) < 10:
            self.setWindowTitle("Mesh Control")
            self.btnModify = self.createButton("&Modify", self.applyMeshControl)
            self.btnModify.setFixedWidth(50)
            #self.btnModify.setFixedHeight(20)

            propLayout=QFormLayout()
            self.mainLayout.addLayout(propLayout)
            self.mainLayout.addWidget(self.btnModify)
            self.txtMCupptresh = QLineEdit()

            propLayout.addRow("&Upper treshold:", self.txtMCupptresh)

            self.txtMClowtresh = QLineEdit()


            propLayout.addRow("&Lower treshold:", self.txtMClowtresh)

        self.txtMClowtresh.setText(str(self.mc.lowertreshold))
        self.txtMCupptresh.setText(str(self.mc.uppertreshold))

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
        input_data = self.currentHUS.model_results['Lusa iteration results Sagg'].data
        input_names = self.currentHUS.model_results['Lusa iteration results Sagg'].column_names
        input_data2 = self.currentHUS.model_results['Lusa iteration results Hogg'].data
        input_names2 = self.currentHUS.model_results['Lusa iteration results Hogg'].column_names

        for name in input_names2:
            input_names.append(name)
        for i in range(len(input_data)):
            input_data[i]= input_data[i]+input_data2[i]

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
        self.model.add_mapping(seriesColorHex, QRect(iy, 0, 2, self.model.rowCount()))

        ix = 6
        iy = 5
        series = QtCharts.QLineSeries()
        series.setName(input_names[iy] + " " + input_names[ix])
        mapper = QtCharts.QVXYModelMapper(self)
        mapper.setXColumn(ix)
        mapper.setYColumn(iy)
        mapper.setSeries(series)
        mapper.setModel(self.model)
        self.chart.addSeries(series)
        # get the color of the series and use it for showing the mapped area
        seriesColorHex = "{}".format(series.pen().color().name())
        self.model.add_mapping(seriesColorHex, QRect(iy, 0, 2, self.model.rowCount()))

        self.chart.createDefaultAxes()

        axX=(self.chart.axes(Qt.Horizontal))[0]
        axX.setTitleText("Ime x")
        axY = (self.chart.axes(Qt.Vertical))[0]
        axY.setTitleText("Ime Y")




    def setCurrentHUS(self, currentHUS):
        self.currentHUS = currentHUS
        self.setWindowTitle("Hull Ultimate Strength")
        self.refreshResults()

def createCommand():
    return HUSCommand()
