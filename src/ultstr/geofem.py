import openmesh as om
from enum import Enum
import numpy as np
from geometry import Geometry

class ViewType(Enum):
    """!
    ViewType
    """

    solid = 0
    eltype = 1
    propplate = 2
    propbeam = 3

class MeshControl():
    def __init__(self):
        self.viewtype=ViewType.solid
        pass
class Property ():
    def __init__(self):
        self.id = 0
        self.name= ''
        pass
    def init(self,id,name):
        self.id = id
        self.name = name
class GeoEntity():
    def __init__(self):
        self.id=0
        pass
    def appendMesh(self, mesh:om.TriMesh,mc:MeshControl):
        pass
    def init(self,id):
        self.id=id

class Group():
    def __init__(self):
        self.entities = []
        pass


class Material (Property):
    def __init__(self):
        super().__init__()
        self.ReH=0
        self.E=0
        self.ni=0
        pass

    def init(self, id, name):
        super().init(id,name)

class RodProperty (Property):
    def __init__(self):
        super().__init__()
        pass
    def init(self, id, name):
        super().init(id,name)

class BeamProperty (Property):
    def __init__(self):
        super().__init__()
        self.hw=0
        self.tw=0
        self.bf=0
        self.tf=0
        self.secType=''
        self.material=0
        pass
    def init(self, id, name):
        super().init(id,name)
class StiffLayoutProperty (Property):
    def __init__(self):
        self.beam = 0
        pass
    def init(self, id, name):
        super().init(id,name)
class PlateProperty (Property):
    def __init__(self):
        self.tp=0
        self.material =0
        pass
    def init(self, id, name):
        super().init(id,name)




class Node(GeoEntity):
    def __init__(self):
        super().__init__()
        self.p = np.array([0, 0, 0])
        pass
    def x(self):
        return self.p[0]
    def y(self):
        return self.p[1]
    def z(self):
        return self.p[2]
    def init(self,id,x,y,z):
        super().init(id)
        self.p[0]= x
        self.p[1] = y
        self.p[2] = z


class Element(GeoEntity):
    def __init__(self):
        super().__init__()
        self.property = 0
        self.nodes = []
        pass
    def addNode(self,node):
        self.nodes.append(node)
    def init(self,id):
        super().init(id)
    def updateMesh(self,mesh:om.TriMesh,mc:MeshControl):
        pass


class RodElement(Element):
    def __init__(self):
        super().__init__()
        pass


class BeamElement(Element):
    def __init__(self):
        super().__init__()
        self.wo = np.array([0, 0, 0]) #web orientation
        pass
    def updateMesh(self,mesh:om.TriMesh,mc:MeshControl):
        pass

class TriaElement(Element):
    def __init__(self):
        super().__init__()
        pass
class StiffTriaElement(TriaElement):
    def __init__(self):
        super().__init__()
        self.layout=0
        pass
    def init(self,id):
        self.id=id


class QuadElement(Element):
    def __init__(self):
        super().__init__()
        pass
    def init(self,id):
        self.id=id


class StiffQuadElement(QuadElement):
    def __init__(self):
        super().__init__()
        self.layout=0
        pass
    def init(self,id):
        self.id=id
    def updateMesh(self,mesh:om.TriMesh,mc:MeshControl):
        #print(self.nodes[0].x(),self.nodes[0].y(),self.nodes[0].z())
        #print(self.nodes[0].p)


        pass


class Units():
    def __init__(self):
        self.user2si_length=1
        self.user2si_force=1
        self.name_length='m'
        self.name_force = 'N'
        pass

class GeoFEM(Geometry):
    def __init__(self):
        super().__init__()
        self.nodes = {}
        self.elements = {}
        self.materials = {}
        self.properties = {}
        self.stiflayouts = {}
        self.groups = {}
        self.meshcontrol = 0
        self.units= Units()
        self.mc = MeshControl()
        pass
    def addNode(self,item):
        self.nodes[item.id]=item
    def getNode(self,id):
        return self.nodes[id]
    def getMaterial(self,id):
        return self.materials[id]
    def getProperty(self,id):
        return self.properties[id]
    def getStiffLayout(self,id):
        return self.stiflayouts[id]
    def addElement(self,item):
        self.elements[item.id]=item
    def addProperty(self,item):
        self.properties[item.id]=item
    def addMaterial(self,item):
        self.materials[item.id]=item
    def addStiffLayout(self, item):
        self.stiflayouts[item.id] = item
    def regenerate(self):
        self.mesh= om.TriMesh()
        for el in self.elements.values():
            el.updateMesh(self.mesh,self.mc)
        pass