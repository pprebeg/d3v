from PySide2.QtCore import Slot
from PySide2.QtGui import QOpenGLShaderProgram, QOpenGLShader
from PySide2.QtGui import QOpenGLVersionProfile, QOpenGLContext
from PySide2.QtGui import QSurfaceFormat
from PySide2.QtWidgets import QMessageBox
from painters import Painter, PainterSignals
from signals import Signals, DragInfo
from painterbasic.glvertdatasforhaders import VertDataCollectorCoord3fNormal3fColor4f
from painterbasic.glhelp import GLEntityType
from OpenGL import GL
from PySide2.QtCore import QCoreApplication
from geometry import Geometry
import openmesh as om
import numpy as np
from selinfo import SelectionInfo

class BasicInfoPainter(Painter):
    def __init__(self):
        #super().__init__()
        super(Painter, self).__init__()
        self._dentsvertsdata = {}  # dictionary that holds vertex data for all primitive and  submodel combinations
        self._geo2Add = []
        self._doSelection=False
        self._si=SelectionInfo()
        self.program = 0
        self.projMatrixLoc = 0
        self.mvMatrixLoc = 0
        self.normalMatrixLoc = 0
        self.lightPosLoc = 0
        # self.vertexShader = self.vertexShaderSourceCore()
        # self.fragmentShader = self.fragmentShaderSourceCore()
        self.vertexShader = self.vertexShaderSource()
        self.fragmentShader = self.fragmentShaderSource()
        # model / geometry
        self.addGeoCount=0
        self.drawLegend = False
        self.legendValues = []
        self.legendColors = []
        self._doLegend = False
        self.width=0
        self.height=0

    def initializeGL(self):
        super().initializeGL()
        self.program = QOpenGLShaderProgram()
        # profile = QOpenGLVersionProfile()
        # profile.setVersion(2, 0)
        #context = QOpenGLContext.currentContext()
        #print("paintr init "+str(context))
        # self.glf = context.versionFunctions(profile)
        # if not self.glf:
        #     QMessageBox.critical(None, "Failed to Initialize OpenGL",
        #                          "Could not initialize OpenGL. This program requires OpenGL x.x or higher. Please check your video card drivers.")
        self.glf.initializeOpenGLFunctions()
        self.glf.glClearColor(0.0, 0.0, 0.0, 1)
        self.program.addShaderFromSourceCode(QOpenGLShader.Vertex, self.vertexShader)
        self.program.addShaderFromSourceCode(QOpenGLShader.Fragment, self.fragmentShader)
        self.program.link()
        self.program.bind()
        self.projMatrixLoc = self.program.uniformLocation("projMatrix")
        self.mvMatrixLoc = self.program.uniformLocation("mvMatrix")
        self.normalMatrixLoc = self.program.uniformLocation("normalMatrix")
        self.lightPosLoc = self.program.uniformLocation("lightPos")
        self.program.release()

    def setprogramvalues(self, proj, mv, normalMatrix, lightpos):
        self.program.bind()
        self.program.setUniformValue(self.lightPosLoc, lightpos)
        self.program.setUniformValue(self.projMatrixLoc, proj)
        self.program.setUniformValue(self.mvMatrixLoc, mv)
        self.program.setUniformValue(self.normalMatrixLoc, normalMatrix)
        self.program.release()

    def paintGL(self):
        self.updateGeometry()
        super().paintGL()
        #self.glf.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        #self.glf.glEnable(GL.GL_DEPTH_TEST)
        self.glf.glDisable(GL.GL_DEPTH_TEST)
        #self.glf.glEnable(GL.GL_CULL_FACE)
        self.glf.glDisable(GL.GL_CULL_FACE)
        self.program.bind()
        for key, value in self._dentsvertsdata.items():
            value.drawvao(self.glf)
        self.program.release()

    def resizeGL(self, w:int, h:int):
        super().resizeGL(w,h)
        self.width=w
        self.height=h
        self._doLegend=True

    def updateGL(self):
        super().updateGL()

    def resetmodel(self):
        """!
        Reset the model

        Cleans the dictionary
        """
        for key, value in self._dentsvertsdata.items():
            value.free()
        self._dentsvertsdata.clear()

    def initnewdictitem(self, key, enttype):
        """!
        Initialize a new dictionary item that holds data for rendering
        @param key: (\b str) item key
        @param enttype: (GLEntityType) primitive drawing entity type
        @retval None
        """

        self._dentsvertsdata[key] = VertDataCollectorCoord3fNormal3fColor4f(enttype)


    def appendlistdata_f3xyzf3nf4rgba(self, key, x, y, z, nx, ny, nz, r, g, b, a):
        """!
        Append Vertex collector dictionary item with new vertex data
        @param key: (\b str) dictonary key
        @param x: (\b float) x coordinate
        @param y: (\b float) y coordinate
        @param z: (\b float) z coordinate
        @param nx: (\b float) x normal coordinate
        @param ny: (\b float) y normal coordinate
        @param nz: (\b float) z normal coordinate
        @retval: (\b int) index of the added vertex
        """
        return self._dentsvertsdata[key].appendlistdata_f3xyzf3nf4rgba(x, y, z, nx, ny, nz,r,g,b,a)

    def appenddictitemsize(self, key, numents):
        """!
        Append dictionary item size with the specified number of entities
        :@param key:(str) key
        :@param numents:(\b int) number of entities to be added
        """
        self._dentsvertsdata[key].appendsize(numents)

    def allocatememory(self):
        """!
        Allocate memory for all dictionary items that holds data for rendering

        Allocation size is based on the information collected by client calls to appenddictitemsize()
        """
        for key, value in self._dentsvertsdata.items():
            value.allocatememory()
    def allocatememory(self,key):
        """!
        Allocate memory for all dictionary items that holds data for rendering

        Allocation size is based on the information collected by client calls to appenddictitemsize()
        """
        self._dentsvertsdata[key].allocatememory()

    def bindData(self, key):
        self._dentsvertsdata[key].setupVertexAttribs(self.glf)
        atrList = self._dentsvertsdata[key].GetAtrList()
        for ent in atrList:
            self.program.bindAttributeLocation(ent[0], ent[1])


# Shader code ********************************************************
    def vertexShaderSourceCore(self):
        return """#version 150
                in vec4 vertex;
                in vec3 normal;
                out vec3 vert;
                out vec3 vertNormal;
                out vec4 colorV;
                uniform mat4 projMatrix;
                uniform mat4 mvMatrix;
                uniform mat3 normalMatrix;
                void main() {
                   vert = vertex.xyz;
                   vertNormal = normalMatrix * normal;
                   gl_Position = vec4(vertex, 1.0);
                   colorV = color;
                }"""
    def fragmentShaderSourceCore(self):
        return """#version 150
                in highp vec3 vert;
                in highp vec3 vertNormal;
                in highp vec4 colorV; 
                out highp vec4 fragColor;
                uniform highp vec3 lightPos;
                void main() {
                   highp vec3 L = normalize(lightPos - vert);
                   highp float NL = max(dot(normalize(vertNormal), L), 0.0);
                   highp vec3 col = clamp(colorV.rgb * 0.2 + colorV.rgb * 0.8 * NL, 0.0, 1.0);
                   fragColor = vec4(col, colorV.a);
                }"""
    def vertexShaderSource(self):
        return """attribute vec4 vertex;
                attribute vec3 normal;
                attribute vec4 color;
                varying vec3 vert;
                varying vec3 vertNormal;
                varying vec4 colorV;
                uniform mat4 projMatrix;
                uniform mat4 mvMatrix;
                uniform mat3 normalMatrix;
                void main() {
                   vert = vertex.xyz;
                   vertNormal = normalMatrix * normal;
                   gl_Position = vec4(vertex, 1.0);
                   colorV = color;
                }"""
    #gl_Position = projMatrix * mvMatrix * vertex;
    #gl_Position = vec4(vertex, 1.0);
    #gl_Position = projMatrix * vertex;
    def fragmentShaderSource(self):
        return """varying highp vec3 vert;
                varying highp vec3 vertNormal;
                varying highp vec4 colorV; 
                uniform highp vec3 lightPos;
                void main() {
                   highp vec3 L = normalize(lightPos - vert);
                   highp float NL = max(dot(normalize(vertNormal), L), 0.0);
                   highp vec3 col = clamp(colorV.rgb * 0.2 + colorV.rgb * 0.8 * NL, 0.0, 1.0);
                   gl_FragColor = vec4(0,0,1,1);
                }"""

# Painter methods implementation code ********************************************************

    def addGeometry(self, geometry:Geometry):
        if hasattr(geometry, 'drawLegend'):
            self.drawLegend = geometry.drawLegend
            self._doLegend=True
            self.legendColors=geometry.legendColors
            self.legendValues = geometry.legendValues
        #self.requestUpdateGL()

    def addLegend(self):
        key= 0
        if key in self._dentsvertsdata:
            self._dentsvertsdata[key].free()
            self._dentsvertsdata.pop(key)
            pass
        if len(self.legendValues) > 0:
            self.initnewdictitem(key, GLEntityType.TRIA)
            nf = len(self.legendValues)*2
            self.appenddictitemsize(key, nf)
            self.allocatememory(key)
            self.addLegendData4oglmdl(key)
            self.bindData(key)

    def updateGeometry(self):
        if self._doLegend:
            self.addLegend()
            self._doLegend=False


    def addLegendData4oglmdl(self, key):
        iCol=0
        lqYpos = 10/self.height;
        lqXpos = 10/self.width+0.5;
        lqW = 20/self.width;
        lqH = 20/self.height;
        lqGap = 10/self.height;
        lqText = lqXpos + lqW + 10/self.height;
        lqYpos += lqH + lqGap;
        pts=[]
        for fh in self.legendValues:
            pts.clear()
            n = [0,0,1]
            c = self.legendColors[iCol]
            iCol = iCol+1
            pts.append([lqXpos, lqYpos,0]);
            pts.append([lqXpos + lqW, lqYpos,0]);
            pts.append([lqXpos + lqW, lqYpos + lqH,0]);
            pts.append([lqXpos + lqW, lqYpos + lqH, 0]);
            pts.append([lqXpos, lqYpos + lqH,0]);
            pts.append([lqXpos, lqYpos, 0]);
            lqYpos += lqH + lqGap;
            for p in pts:  # vertex handle
                self.appendlistdata_f3xyzf3nf4rgba(key,
                                                   p[0], p[1], p[2],
                                                   n[0], n[1], n[2],
                                                   c[0], c[1], c[2], c[3])
        return




