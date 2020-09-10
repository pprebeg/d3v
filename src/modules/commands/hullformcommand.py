from PySide2.QtWidgets import QApplication, QMenu, QMessageBox,QFormLayout,QWidget,QHeaderView
from PySide2.QtWidgets import QDialog, QPushButton,QGridLayout,QVBoxLayout,QHBoxLayout,QTableView,QTextEdit,QLabel
from commands import Command
from iohandlers import IOHandler
import openmesh as om
import numpy as np
from signals import Signals
from geometry import Geometry
from hullform import HullForm
import os
from PySide2.QtCore import Slot,Qt,SIGNAL
from PySide2.QtCore import QAbstractTableModel, QModelIndex, QRect
from PySide2.QtGui import QColor, QPainter
from PySide2.QtCharts import QtCharts

class HullFormCommand(Command):
    def __init__(self):
        super().__init__()
        app = QApplication.instance()
        app.registerIOHandler(HullFormImporter())
        self.mainwin = app.mainFrame
        self.hf_prop = DialogHullFormModify(self.mainwin)
        self.hf_hydroCurves = DialogHullFormHydrostaticCurves(self.mainwin)
        self.hf=0
        self.si=0

        #tools = app.mainFrame.menuTools
        mb = app.mainFrame.menuBar()

        self.menuMain = QMenu("Hull Form")

        self.menuHullFormModify = self.menuMain.addAction("&Modify form")
        self.menuHullFormModify.triggered.connect(self.onModifyForm)

        self.menuHullFormResults = QMenu("&Results")
        self.menu2 = QMenu("&Menu2")
        self.menuMain.addMenu(self.menuHullFormResults)
        self.menuMain.addMenu(self.menu2)



        menuResultHydrostaticCurves = self.menuHullFormResults.addAction("Hydrostatic Curves")
        menuResultHydrostaticCurves.triggered.connect(self.onShowHydrostaticCurves)



        #tools.addMenu(menu)
        mb.addMenu(self.menuMain)
        #self.menuMain.setEnabled(False)


        Signals.get().geometryAdded.connect(self.registerHullForm)
        Signals.get().selectionChanged.connect(self.registerSelection)

    @Slot()
    def registerHullForm(self, hullForm):
        if isinstance(hullForm, HullForm):
            self.hf=hullForm
            self.menuMain.setEnabled(True)

    @Slot()
    def registerSelection(self, si):
        self.si=si

    def onModifyForm(self):
        if isinstance(self.hf, HullForm):
            self.hf_prop.setCurrentHullForm(self.hf)
            self.hf_prop.exec()

    def onShowHydrostaticCurves(self):
        if isinstance(self.hf,HullForm):
            self.hf_hydroCurves.setCurrentHullForm(self.hf)
            self.hf_hydroCurves.exec()



class HullFormImporter(IOHandler):
    def __init__(self):
        super().__init__()

    def importGeometry(self, fileName):
        if len(fileName) < 1:
            return
        filename, file_extension = os.path.splitext(fileName)
        if file_extension != ".huf":
            return
        hf = HullForm(fileName)
        Signals.get().geometryImported.emit(hf)

    def getImportFormats(self):
        return (".huf")


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
            self.input_names[section]
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


class DialogHullFormHydrostaticCurves(QDialog):
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
        self.chart_view.setMinimumSize(640, 480)


        self.setWindowTitle("Hull Form Hydrostatic Curves")
        self.btnGenerate = self.createButton("&Generate", self.refreshResults)

        self.txtMaxWL = QTextEdit()
        self.txtMaxWL.setFixedHeight(sizetxt)
        self.txtMaxWL.setText('9.0')
        self.txtMaxWL.setAlignment(Qt.AlignRight)
        self.txtWLStep = QTextEdit()
        self.txtWLStep.setFixedHeight(sizetxt)
        self.txtWLStep.setText('0.5')
        self.txtWLStep.setAlignment(Qt.AlignRight)


        mainLayout = QVBoxLayout()
        mainLayout.setStretch(0,1)
        mainLayout.setStretch(1, 0)
        tablechartLayout = QGridLayout()


        controlLayout = QHBoxLayout()
        controlWidget = QWidget()
        controlWidget.setFixedHeight(sizetxt*3)
        controlWidget.setLayout(controlLayout)

        inputLayout = QFormLayout()
        controlLayout.addLayout(inputLayout)
        controlLayout.addWidget(self.btnGenerate)

        inputLayout.addRow("&Max. Waterline height:", self.txtMaxWL)
        inputLayout.addRow("&Waterline step:", self.txtWLStep)

        tablechartLayout.addWidget(self.table_view, 0, 0)
        tablechartLayout.addWidget(self.chart_view, 0, 1)
        mainLayout.addLayout(tablechartLayout)
        mainLayout.addLayout(controlLayout)
        mainLayout.addWidget(controlWidget)

        self.setLayout(mainLayout)
        self.currentHullForm=0


    def createButton(self, text, member):
        button = QPushButton(text)
        button.clicked.connect(member)
        return button

    def refreshResults(self):
        #self.currentHullForm.getResults(9, 1.025)
        #return
        input_data = []
        maxWL= float(self.txtMaxWL.toPlainText())
        stepWL = float(self.txtWLStep.toPlainText())
        h=maxWL
        while h > 0:
            result = self.currentHullForm.getResults(h, 1.025)
            input_data.append(result)
            h=h-stepWL
            if h <= 0:
                result = self.currentHullForm.getResults(0.001, 1.025)
                input_data.append(result)


        input_names = ['h', 'Volume', 'Awl', 'Xwl', 'KBz', 'KBx', 'Ib', 'Il', 'Lwl', 'Bwl', 'MF area',
                       'MoB','KMo','MlB','KMl','JZ','Cwl','CB','CP','CX']
        #input_names = ['h', 'Volume', 'Awl']
#        self.model.layoutAboutToBeChanged()
        self.model.setInputData(input_names, input_data)
        #self.model.layoutChanged()
        self.chart.removeAllSeries()
        seriesColorHex = "#000000"
        for i in range(1, len(input_names)):
            series = QtCharts.QLineSeries()
            series.setName(input_names[i])
            mapper = QtCharts.QVXYModelMapper(self)
            mapper.setXColumn(0)
            mapper.setYColumn(i)
            mapper.setSeries(series)
            mapper.setModel(self.model)
            self.chart.addSeries(series)
            # get the color of the series and use it for showing the mapped area
            seriesColorHex = "{}".format(series.pen().color().name())
            self.model.add_mapping(seriesColorHex, QRect(i, 0, 1, self.model.rowCount()))

        self.chart.createDefaultAxes()


    def setCurrentHullForm(self, currentHullForm):
        self.currentHullForm = currentHullForm
        self.setWindowTitle("Hull Form Hydrostatic Curves")



def createCommand():
    return HullFormCommand()
