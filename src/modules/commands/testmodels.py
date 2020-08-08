from PySide2.QtWidgets import QApplication, QMenu, QMessageBox
from commands import Command

import openmesh as om
import numpy as np
from signals import Signals
from geometry import Geometry
from PySide2.QtCore import Slot
import time

class TestModelsCommand(Command):
    def __init__(self):
        super().__init__()
        app = QApplication.instance()
        tools = app.mainFrame.menuTools
        menu = QMenu("DataModelsTest ...")
        lof = menu.addAction(" ... loop on faces")
        lovef = menu.addAction(" ... loop on vertices in each face")
        lof.triggered.connect(self.onTestLoopOnFaces)
        lovef.triggered.connect(self.onTestLoopOnVerticesInEachFace)
        tools.addMenu(menu)
        Signals.get().geometryAdded.connect(self.onGeometryAdded)
        self.geometries=[]
        self.dataModels = []
        self.dataModels.append(TestDataModelOpenMesh())
        self.dataModels.append(TestDataModelPythonListAndDict())
        pass

    @Slot()
    def onGeometryAdded(self, geometry:Geometry):
        self.geometries.append(geometry)
        for tm in self.dataModels:
            tm.initialize(geometry)

    def onTestLoopOnFaces(self):
        print("*** Test Loop on faces started! ***")
        for tm in self.dataModels:
            tStart = time.perf_counter()
            msg = tm.doTestLoopOnFaces()
            dt = time.perf_counter() - tStart
            print(tm.name+", s: ", dt)
            print(tm.name+" msg.: ",msg)
        print("*** Test Loop on faces ended! ***")

    def onTestLoopOnVerticesInEachFace(self):
        print("*** Test LoopOnVerticesInEachFace started! ***")
        for tm in self.dataModels:
            tStart = time.perf_counter()
            msg=tm.doTestLoopOnVerticesInEachFace()
            dt = time.perf_counter() - tStart
            print(tm.name + ", s: ", dt)
            print(tm.name + " msg.: ", msg)
        print("*** Test LoopOnVerticesInEachFace ended! ***")

class TestDataModel:
    def __init__(self):
        self._name=""

    @property
    def name(self):
        return self._name
    def initialize(self,geometry:Geometry):
        pass

    def doTestLoopOnFaces(self):
        pass
    def doTestLoopOnVerticesInEachFace(self):
       pass

class TestDataModelOpenMesh(TestDataModel):
    def __init__(self):
        super().__init__()
        self._meshes=[]
        self._name = "Open Mesh"

    def initialize(self,geometry:Geometry):
        self._meshes.append(geometry.mesh)

    def doTestLoopOnFaces(self):
        msg = ""
        i=0
        for mesh in self._meshes:
            for fh in mesh.faces():
                i = i + 1
        msg="Total number of faces in loop: "
        msg=msg+str(i)
        return msg

    def doTestLoopOnVerticesInEachFace(self):
        msg = ""
        i=0
        for mesh in self._meshes:
            for fh in mesh.faces():
                for vh in mesh.fv(fh):
                    p = mesh.point(vh)
                    i = i + 1
        msg = "Total number of points in loop: "
        msg = msg + str(i)
        return msg
class TestDataModelPythonListAndDict(TestDataModel):
    def __init__(self):
        super().__init__()
        self._faces=[]
        self._vertices = []
        self._name = "PythonListAndDict"


    def initialize(self,geometry:Geometry):
        nf=geometry.mesh.n_faces()
        nv=geometry.mesh.n_vertices()
        faceDict={}
        vertexDict = {}
        for fh in geometry.mesh.faces():
            vertList=[]
            for vh in geometry.mesh.fv(fh):
                vertList.append(vh.idx())
            faceDict[fh.idx()]=vertList
        for vh in geometry.mesh.vertices():
            vertCoord = []
            p = geometry.mesh.point(vh)
            vertCoord.append(p[0])
            vertCoord.append(p[1])
            vertCoord.append(p[2])
            vertexDict[vh.idx()] = vertCoord
        self._faces.append(faceDict)
        self._vertices.append(vertexDict)

    def doTestLoopOnFaces(self):
        msg = ""
        i=0
        for faceDict in self._faces:
            for fh in faceDict:
                i = i + 1
        msg="Total number of faces in loop: "
        msg=msg+str(i)
        return msg

    def doTestLoopOnVerticesInEachFace(self):
        msg = ""
        i=0
        ig=0
        for ig in range(len(self._faces)):
            faceDict=self._faces[ig]
            for fhKey in faceDict:
                vidxList=faceDict[fhKey]
                for vidx in vidxList:
                    p = self._vertices[ig][vidx]
                    i = i + 1
        msg = "Total number of points in loop: "
        msg = msg + str(i)
        return msg

def createCommand():
    return TestModelsCommand()
