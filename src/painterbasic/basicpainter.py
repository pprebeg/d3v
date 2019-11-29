from PySide2.QtCore import Slot
from PySide2.QtGui import QOpenGLShaderProgram, QOpenGLShader

from painters import Painter
from signals import Signals, DragInfo
from painterbasic.glvertdatasforhaders import VertDataCollectorCoord3fNormal3fColor4f
from painterbasic.subdivsphere import Sphere
from painterbasic.glhelp import GLEntityType


class BasicPainter(Painter):
    def __init__(self):
        self._dentsvertsdata = {}  # dictionary that holds vertex data for all primitive and  submodel combinations
        Painter.__init__(self)
        Signals.get().dragging.connect(self.bgchanger)
        self.program = QOpenGLShaderProgram()
        #self.vertexShader = self.vertexShaderSourceCore()
        #self.fragmentShader = self.fragmentShaderSourceCore()
        self.vertexShader = self.vertexShaderSource()
        self.fragmentShader = self.fragmentShaderSource()
        #self.core = "--coreprofile" in QCoreApplication.arguments()
        self.projMatrixLoc = 0
        self.mvMatrixLoc = 0
        self.normalMatrixLoc = 0
        self.lightPosLoc = 0

        # model / geometry
        self.spheres = []
        self.gengeometry()

    def gengeometry(self):
        curSphere = Sphere(0, 0, 0, 0.5)  # pass X, Y, Z, radius
        self.spheres.append(curSphere)
        curSphere = Sphere(1.5, -1.5, 0, 2)
        self.spheres.append(curSphere)
        curSphere = Sphere(1, 2, 8, 0.4)
        self.spheres.append(curSphere)
        curSphere = Sphere(1, 1, 2, 0.4)
        self.spheres.append(curSphere)

        self.resetmodel()
        self.initnewdictitem("sphere", GLEntityType.TRIA)

        for i in range(0,len(self.spheres)):
            self.appenddictitemsize("sphere", self.spheres[i].getnumtria())

        self.allocatememory()
        self.adddata4oglmdl()

    def adddata4oglmdl(self):

        for isp in range(0,len(self.spheres)):
            nt = self.spheres[isp].getnumtria()
            for i in range(0, nt):  # interate through all of the triangles
                iimin = i * 3
                iimax = iimin + 3
                for ii in range(iimin, iimax):
                    iv = self.spheres[isp].indices[ii] * 3  # each vertex has xyz
                    ic = self.spheres[isp].indices[ii] * 4  # each vertex has xyz
                    sph=self.spheres[isp]
                    # self.appendlistdata_f3xyzf3n("sphere",
                    #     sph.vertices[iv],sph.vertices[iv + 1],sph.vertices[iv + 2],
                    #     sph.normals[iv], sph.normals[iv + 1], sph.normals[iv + 2])
                    self.appendlistdata_f3xyzf3nf4rgba("sphere",
                        sph.vertices[iv],sph.vertices[iv + 1],sph.vertices[iv + 2],
                        sph.normals[iv], sph.normals[iv + 1], sph.normals[iv + 2],
                        sph.colors[ic], sph.colors[ic + 1], sph.colors[ic + 2],sph.colors[ic + 3])
        pass

    def initializeGL(self):
        self.glf.initializeOpenGLFunctions()
        self.glf.glClearColor(0.0, 0.0, 0.0, 1)

        for key, value in self._dentsvertsdata.items():
            value.setupVertexAttribs(self.glf)

        self.program.addShaderFromSourceCode(QOpenGLShader.Vertex, self.vertexShader)
        self.program.addShaderFromSourceCode(QOpenGLShader.Fragment, self.fragmentShader)

        atrList = self._dentsvertsdata[list(self._dentsvertsdata.keys())[0]].GetAtrList()
        for ent in atrList:
            self.program.bindAttributeLocation(ent[0], ent[1])

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
        self.program.bind()
        for key, value in self._dentsvertsdata.items():
            value.drawvao(self.glf)
        self.program.release()

    def resizeGL(self, w:int, h:int):
        pass

    @Slot(DragInfo)
    def bgchanger(self, di: DragInfo):

        Signals.get().updateGL.emit()

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
                   gl_Position = projMatrix * mvMatrix * vertex;
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
                   gl_Position = projMatrix * mvMatrix * vertex;
                   colorV = color;
                }"""

    def fragmentShaderSource(self):
        return """varying highp vec3 vert;
                varying highp vec3 vertNormal;
                varying highp vec4 colorV; 
                uniform highp vec3 lightPos;
                void main() {
                   highp vec3 L = normalize(lightPos - vert);
                   highp float NL = max(dot(normalize(vertNormal), L), 0.0);
                   highp vec3 col = clamp(colorV.rgb * 0.2 + colorV.rgb * 0.8 * NL, 0.0, 1.0);
                   gl_FragColor = vec4(col, colorV.a);
                }"""
# ********************************************************************
#gl_FragColor = colorV;