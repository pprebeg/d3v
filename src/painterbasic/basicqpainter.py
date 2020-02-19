from PySide2.QtCore import Slot,Qt
from PySide2.QtGui import QOpenGLShaderProgram, QOpenGLShader
from PySide2.QtGui import QOpenGLVersionProfile, QOpenGLContext
from PySide2.QtGui import QSurfaceFormat,QFont
from PySide2.QtWidgets import QMessageBox
from painters import Painter, PainterSignals
from signals import Signals, DragInfo
from PySide2.QtCore import QPoint, QRect
from geometry import Geometry
import openmesh as om
import numpy as np
from selinfo import SelectionInfo
from PySide2.QtGui import QBrush, QPainter,QPen ,QPolygon,QColor


class BasicQPainter(Painter):
    def __init__(self):
        #super().__init__()
        super(Painter, self).__init__()

        self.drawLegend = False
        self.legendValues = []
        self.legendColors = []
        self.legendTitle=""
        self._doLegend = False
        self.width=0
        self.height=0
        self.pen = QPen(QColor(0, 10, 127, 255), 2)
        self.brush = QBrush()
        self.paintDevice=0

    def initializeGL(self):
        pass

    def setprogramvalues(self, proj, mv, normalMatrix, lightpos):
        pass
    def paintGL(self):
        super().paintGL()
        if not self.drawLegend:
            return
        painter = QPainter(self.paintDevice)
        font= QFont("Times", 10, QFont.Bold)
        painter.setFont(font)

        #painter.setBrush(self.brush)
        color = QColor()
        iCol = 0
        lqYpos = 10;
        lqXpos = 10 + self.width - 50;
        lqW = 30;
        lqH = 20;
        lqGap = 10;
        lqText = lqXpos + lqW + 10;
        lqYpos += lqH + lqGap;
        pts = []
        for legVal in self.legendValues:
            pts.clear()
            n = [0, 0, 1]
            c = self.legendColors[iCol]
            color.setRgbF(c[0],c[1],c[2],c[3])
            self.pen = QPen(QColor(0, 0, 0, 255), 1)
            self.brush=QBrush(color,Qt.BrushStyle(Qt.SolidPattern))
            painter.setPen(self.pen)
            painter.setBrush(self.brush)
            iCol = iCol + 1
            pts.append(QPoint(lqXpos, lqYpos));
            pts.append(QPoint(lqXpos + lqW, lqYpos));
            pts.append(QPoint(lqXpos + lqW, lqYpos + lqH));
            pts.append(QPoint(lqXpos, lqYpos + lqH));
            rect = QRect(lqXpos-lqW, lqYpos, lqW, lqH)
            lqYpos += lqH + lqGap;
            poly = QPolygon(pts)

            painter.drawPolygon(poly)
            painter.drawText(rect,Qt.AlignCenter,legVal)
        rect = QRect(lqXpos -100, 0, 130, lqH)
        painter.drawText(rect, Qt.AlignLeft, self.legendTitle)

    def resizeGL(self, w:int, h:int):
        super().resizeGL(w,h)
        self.width=w
        self.height=h

    def updateGL(self):
        super().updateGL()


    def addGeometry(self, geometry:Geometry):
        if hasattr(geometry, 'drawLegend'):
            self.drawLegend = geometry.drawLegend
            self.legendColors=geometry.legendColors
            self.legendValues = geometry.legendValues
            self.legendTitle = geometry.legendTitle
        #self.requestUpdateGL()







