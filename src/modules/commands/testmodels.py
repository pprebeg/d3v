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
        self.dataModels.append(TestDataModelOpenMeshIndices())
        self.dataModels.append(TestDataModelOpenMeshIndicesGlobal())
        self.dataModels.append(TestDataModelOpenMeshIndicesIndexLoop())
        self.dataModels.append(TestDataModelOpenMeshIndicesGlobalToList())
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
class TestDataModelOpenMeshIndices(TestDataModel):
    def __init__(self):
        super().__init__()
        self._meshes=[]
        self._name = "Open Mesh Indices"

    def initialize(self,geometry:Geometry):
        self._meshes.append(geometry.mesh)

    def doTestLoopOnFaces(self):
        msg = ""
        i=0

        for mesh in self._meshes:
            fv_indices = mesh.fv_indices()
            nfv=fv_indices.shape[0]
            for fv in fv_indices:
                i = i + 1
        msg="Total number of faces in loop: "
        msg=msg+str(i)
        return msg

    def doTestLoopOnVerticesInEachFace(self):
        msg = ""
        i=0
        for mesh in self._meshes:
            fv_indices = mesh.fv_indices()
            points = mesh.points()
            for fv in fv_indices:
                for iv in fv:
                   p = points[iv]
                   i = i + 1
        msg = "Total number of points in loop: "
        msg = msg + str(i)
        return msg
class TestDataModelOpenMeshIndicesIndexLoop(TestDataModel):
    def __init__(self):
        super().__init__()
        self._meshes=[]
        self._name = "Open Mesh Indices Index Loop"

    def initialize(self,geometry:Geometry):
        self._meshes.append(geometry.mesh)

    def doTestLoopOnFaces(self):
        msg = ""
        i=0

        for mesh in self._meshes:
            fv_indices = mesh.fv_indices()
            nfv=fv_indices.shape[0]
            for ifv in range(nfv):
                fv=fv_indices[ifv]
                i = i + 1
        msg="Total number of faces in loop: "
        msg=msg+str(i)
        return msg

    def doTestLoopOnVerticesInEachFace(self):
        msg = ""
        i=0
        for mesh in self._meshes:
            fv_indices = mesh.fv_indices()
            points = mesh.points()
            nfv = fv_indices.shape[0]
            for ifv in range(nfv):
                fv = fv_indices[ifv]
                nvp=fv.shape[0]
                for iiv in range(nvp):
                   iv=fv[iiv]
                   p = points[iv]
                   i = i + 1
        msg = "Total number of points in loop: "
        msg = msg + str(i)
        return msg
class TestDataModelOpenMeshIndicesGlobal(TestDataModel):
    def __init__(self):
        super().__init__()
        self._meshes=[]
        self._name = "Open Mesh Indices Global"
        self._points = []
        self._fv_indices = []

    def initialize(self,geometry:Geometry):
        self._meshes.append(geometry.mesh)
        self._fv_indices.append(geometry.mesh.fv_indices())
        self._points.append(geometry.mesh.points())

    def doTestLoopOnFaces(self):
        msg = ""
        i=0

        for fv_indices in self._fv_indices:
            for fv in fv_indices:
                i = i + 1
        msg="Total number of faces in loop: "
        msg=msg+str(i)
        return msg

    def doTestLoopOnVerticesInEachFace(self):
        msg = ""
        i=0
        for im in range(len(self._meshes)):
            fv_indices = self._fv_indices[im]
            points = self._points[im]
            for fv in fv_indices:
                for iv in fv:
                    p = points[iv]
                    i = i + 1
        msg = "Total number of points in loop: "
        msg = msg + str(i)
        return msg
class TestDataModelOpenMeshIndicesGlobalToList(TestDataModel):
    def __init__(self):
        super().__init__()
        self._meshes=[]
        self._name = "Open Mesh Indices Global To List"
        self._points = []
        self._fv_indices = []

    def initialize(self,geometry:Geometry):
        self._meshes.append(geometry.mesh)
        self._fv_indices.append(geometry.mesh.fv_indices().tolist())
        self._points.append(geometry.mesh.points().tolist())

    def doTestLoopOnFaces(self):
        msg = ""
        i=0

        for fv_indices in self._fv_indices:
            for fv in fv_indices:
                i = i + 1
        msg="Total number of faces in loop: "
        msg=msg+str(i)
        return msg

    def doTestLoopOnVerticesInEachFace(self):
        msg = ""
        i=0
        for im in range(len(self._meshes)):
            fv_indices = self._fv_indices[im]
            points = self._points[im]
            for fv in fv_indices:
                for iv in fv:
                    p = points[iv]
                    i = i + 1
        msg = "Total number of points in loop: "
        msg = msg + str(i)
        return msg
class TestDataModelPythonListAndDict(TestDataModel):
    def __init__(self):
        super().__init__()
        self._fv_dicts=[]
        self._points_dicts = []
        self._name = "PythonListAndDict"


    def initialize(self,geometry:Geometry):
        nf=geometry.mesh.n_faces()
        nv=geometry.mesh.n_vertices()
        fv_dict={}
        points_dict = {}
        for fh in geometry.mesh.faces():
            vertList=[]
            for vh in geometry.mesh.fv(fh):
                vertList.append(vh.idx())
            fv_dict[fh.idx()]=vertList
        for vh in geometry.mesh.vertices():
            vertCoord = []
            p = geometry.mesh.point(vh)
            vertCoord.append(p[0])
            vertCoord.append(p[1])
            vertCoord.append(p[2])
            points_dict[vh.idx()] = vertCoord
        self._fv_dicts.append(fv_dict)
        self._points_dicts.append(points_dict)

    def doTestLoopOnFaces(self):
        msg = ""
        i=0
        for faceDict in self._fv_dicts:
            for fh in faceDict:
                i = i + 1
        msg="Total number of faces in loop: "
        msg=msg+str(i)
        return msg

    def doTestLoopOnVerticesInEachFace(self):
        msg = ""
        i=0
        ig=0
        for ig in range(len(self._fv_dicts)):
            faceDict=self._fv_dicts[ig]
            for fhKey in faceDict:
                vidxList=faceDict[fhKey]
                for vidx in vidxList:
                    p = self._points_dicts[ig][vidx]
                    i = i + 1
        msg = "Total number of points in loop: "
        msg = msg + str(i)
        return msg

def createCommand():
    return TestModelsCommand()
