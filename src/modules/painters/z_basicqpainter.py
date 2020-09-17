from PySide2.QtCore import Slot,Qt
from PySide2.QtGui import QOpenGLShaderProgram, QOpenGLShader
from PySide2.QtGui import QOpenGLVersionProfile, QOpenGLContext
from PySide2.QtGui import QSurfaceFormat,QFont
from PySide2.QtWidgets import QMessageBox
from painters import Painter
from signals import Signals, DragInfo
from PySide2.QtCore import QPoint, QRect
from geometry import Geometry
import openmesh as om
import numpy as np
from selinfo import SelectionInfo
from PySide2.QtGui import QBrush, QPainter,QPen ,QPolygon,QColor
import PySide2.QtCore
from OpenGL import GL
from PySide2.QtWidgets import QApplication


class BasicQPainter(Painter):
    def __init__(self):
        super().__init__()
        self.drawInfo = False
        self.info= ""
        self.drawInfoRect = True
        self.width=0
        self.height=0
        self.colorRect = QColor(100, 100, 255, 255)
        self.colorText = QColor(0, 0, 0, 255)
        self.infoFontName = "Times"
        self.infoFontSize = 10
        self.infoFontType = QFont.Bold
        self.paintDevice =0


        self.drawLegend = False
        self.legendValues = []
        self.legendColors = []
        self.legendTitle=""
        self._doLegend = False

        self.pen = QPen(QColor(0, 10, 127, 255), 2)
        self.brush = QBrush()


    def initializeGL(self):
        paintDevice = QApplication.instance().mainFrame.glWin
        super().initializeGL()
        self.paintDevice = paintDevice
        self.width = paintDevice.vport.width()
        self.height = paintDevice.vport.height()
        self.glf.initializeOpenGLFunctions()





    def setprogramvalues(self, proj, mv, normalMatrix, lightpos):
        pass
    def paintGL(self):
        painter = QPainter(self.paintDevice)
        if self.drawInfo:
            self.glf.glDisable(GL.GL_CULL_FACE)
            font= QFont(self.infoFontName, self.infoFontSize, self.infoFontType)
            painter.setFont(font)
            pen=QPen(self.colorText, 1)
            painter.setPen(pen)
            brush=QBrush(self.colorRect)
            painter.setBrush(brush)

            rect = QRect(5, 5, self.width-10, 20)

            if self.drawInfoRect:
                painter.drawRect(rect)
                #painter.fillRect(rect,self.colorRect)
            painter.drawText(rect, Qt.AlignCenter, self.info)



        if self.drawLegend:
            self.glf.glDisable(GL.GL_CULL_FACE)
            font= QFont("Times", 10, QFont.Normal)
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
                color =QColor(int(c[0]*255),int(c[1]*255),int(c[2]*255),int(c[3]*255))
                self.pen = QPen(QColor(0, 0, 0, 255), 1)
                self.brush=QBrush(color,Qt.BrushStyle(Qt.SolidPattern))
                painter.setPen(self.pen)
                painter.setBrush(self.brush)
                iCol = iCol + 1
                pts.append(QPoint(lqXpos, lqYpos));
                pts.append(QPoint(lqXpos + lqW, lqYpos));
                pts.append(QPoint(lqXpos + lqW, lqYpos + lqH));
                pts.append(QPoint(lqXpos, lqYpos + lqH));
                rect = QRect(lqXpos-2*lqW, lqYpos, 2*lqW, lqH)
                lqYpos += lqH + lqGap;
                poly = QPolygon(pts)

                painter.drawPolygon(poly)
                painter.drawText(rect,Qt.AlignCenter,str(legVal))
            rect = QRect(lqXpos -100, 0, 130, lqH)
            painter.drawText(rect, Qt.AlignCenter, self.legendTitle)

        painter.end()


    def resizeGL(self, w:int, h:int):
        super().resizeGL(w,h)
        self.width=w
        self.height=h

    def updateGL(self):
        super().updateGL()


    def addGeometry(self, geometry:Geometry):
        #self.drawInfo = True
        if isinstance(geometry.mesh, om.TriMesh):
            self.info = self.info+ "Tria Mesh"
        else:
            self.info = self.info+ "Poly Mesh"
        self.info = self.info + "; nf = " + str(geometry.mesh.n_faces())+"; "
        #self.requestGLUpdate()

        if hasattr(geometry, 'drawLegend'):
            self.drawLegend = geometry.drawLegend
            self.legendColors=geometry.legendColors
            self.legendValues = geometry.legendValues
            self.legendTitle = geometry.legendTitle

    def removeGeometry(self, geometry: Geometry):
        pass

    def rebuildGeometry(self, geometry: Geometry):
        if hasattr(geometry, 'drawLegend'):
            self.drawLegend = geometry.drawLegend
            self.legendColors = geometry.legendColors
            self.legendValues = geometry.legendValues
            self.legendTitle = geometry.legendTitle


def createPainter():
    return BasicQPainter()


