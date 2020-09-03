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

    def updateMesh(self,mesh:om.TriMesh,mc:MeshControl):
        print("čvor 1")
        print(self.nodes[0].p)
        print(self.nodes[0].p[0])
        print(self.nodes[0].p[1])
        print(self.nodes[0].p[2])
        print("čvor 1a")
        print(self.nodes[0].p[0])
        print(self.nodes[0].p[1])
        print(self.nodes[0].p[2]-1)
        print("čvor 2")
        print(self.nodes[1].p)

        vhandle = []

        hw=self.property.hw
        tw=self.property.tw
        bf=self.property.bf
        tf = self.property.tf

        vhandle.append(mesh.add_vertex(self.nodes[0].p))
        vhandle.append(mesh.add_vertex(self.nodes[1].p))


        data = self.nodes[0].p + self.wo*hw  #2 točka 3
        vhandle.append(mesh.add_vertex(data))
        data = self.nodes[1].p + self.wo * hw  # 3 točka 4
        vhandle.append(mesh.add_vertex(data))

        # data = np.array([self.nodes[0].p]) + np.array(self.wo)*self.hw - np.cross(np.array([self.nodes[1].p]),np.array(self.wo))*(self.bf*0.5) #4 točak 5
        # vhandle.append(mesh.add_vertex(data))

        fh0 = mesh.add_face(vhandle[0], vhandle[1], vhandle[2]) #1-2-3
        fh1 = mesh.add_face(vhandle[2], vhandle[1], vhandle[3]) #3-2-4



        # ###TEST S RUCNIM UNOSOM TOCAKA###
        # data = np.array([self.nodes[0].p[0], self.nodes[0].p[1], self.nodes[0].p[2]]) #0 točka 1
        # vhandle.append(mesh.add_vertex(data))
        # data = np.array([self.nodes[1].p[0], self.nodes[1].p[1], self.nodes[1].p[2]]) #1 točka 2
        # vhandle.append(mesh.add_vertex(data))
        # data = np.array([self.nodes[1].p[0], self.nodes[1].p[1], self.nodes[1].p[2]-1]) #2 točka 3
        # vhandle.append(mesh.add_vertex(data))
        # data = np.array([self.nodes[0].p[0], self.nodes[0].p[1], self.nodes[0].p[2]-1]) #3 točak 4
        # vhandle.append(mesh.add_vertex(data))
        #
        #
        # data = np.array([self.nodes[1].p[0], self.nodes[1].p[1]-0.25, self.nodes[1].p[2]-1]) #4 točka 5
        # vhandle.append(mesh.add_vertex(data))
        # data = np.array([self.nodes[0].p[0], self.nodes[0].p[1]-0.25, self.nodes[0].p[2]-1]) #5 točak 6
        # vhandle.append(mesh.add_vertex(data))
        #
        #
        # data = np.array([self.nodes[1].p[0], self.nodes[1].p[1]+0.25, self.nodes[1].p[2]-1]) #6 točka 7
        # vhandle.append(mesh.add_vertex(data))
        # data = np.array([self.nodes[0].p[0], self.nodes[0].p[1]+0.25, self.nodes[0].p[2]-1]) #8 točak 7
        # vhandle.append(mesh.add_vertex(data))
        #
        #
        # fh0 = mesh.add_face(vhandle[0], vhandle[1], vhandle[2]) #1-2-3
        # fh1 = mesh.add_face(vhandle[0], vhandle[2], vhandle[3]) #1-3-4
        #
        # ####NERADI S 4 trokutra jer dolazi da moramo po treci puta prolaziti po istom pravcu
        # # fh2 = mesh.add_face(vhandle[3], vhandle[2], vhandle[4]) #4-3-5
        # # fh3 = mesh.add_face(vhandle[3], vhandle[4], vhandle[5]) #4-5-6
        # # fh4 = mesh.add_face(vhandle[3], vhandle[2], vhandle[6])  # 4-3-7
        # # fh5 = mesh.add_face(vhandle[3], vhandle[6], vhandle[7])  # 4-7-8
        #
        # fh2 = mesh.add_face(vhandle[4], vhandle[5], vhandle[7]) #5-6-8
        # fh3 = mesh.add_face(vhandle[4], vhandle[7], vhandle[6]) #5-8-7


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
    def updateMesh(self, mesh: om.TriMesh, mc: MeshControl):
        # print(self.nodes[0].x(),self.nodes[0].y(),self.nodes[0].z())
        # print(self.nodes[1].x(], self.nodes[1].y(], self.nodes[1].z())
        # print(self.nodes[2].x(], self.nodes[2].y(], self.nodes[2].z())
        # print(self.nodes[3].x(], self.nodes[3].y(], self.nodes[3].z())
        # print(self.nodes[4].x(], self.nodes[4].y(], self.nodes[4].z())
        # print(self.nodes[0].p[0], self.nodes[0].p[1], self.nodes[0].p[2], "Ivan")
        # print("čvor 1")
        # print(self.nodes[0].p)
        # print("čvor 2")
        # print(self.nodes[1].p)
        # print("čvor 3")
        # print(self.nodes[2].p)
        # print("čvor 4")
        # print(self.nodes[3].p)
        # mesh = om.TriMesh()
        vhandle = []
        data = np.array([self.nodes[0].p[0], self.nodes[0].p[1], self.nodes[0].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[1].p[0], self.nodes[1].p[1], self.nodes[1].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[2].p[0], self.nodes[2].p[1], self.nodes[2].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[3].p[0], self.nodes[3].p[1], self.nodes[3].p[2]])
        vhandle.append(mesh.add_vertex(data))

        fh0 = mesh.add_face(vhandle[0], vhandle[1], vhandle[2])
        fh1 = mesh.add_face(vhandle[0], vhandle[2], vhandle[3])
        # return mesh
        # pass

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
        mesh= om.TriMesh()
        for el in self.elements.values():
            el.updateMesh(mesh,self.mc)
        pass
        self.mesh = mesh

    # def showFaceColorP(self, propDict):
    #     colors = [[0, 0, 255, 255], [128, 0, 128, 255], [222, 184, 135, 255], [255, 165, 0, 255], [0, 255, 0, 255],
    #               [0, 128, 0, 255], [128, 0, 0, 255], [255, 0, 0, 255], [255, 192, 203, 255], [222, 184, 135, 255],
    #               [255, 165, 0, 255], [255, 127, 80, 255], [128, 128, 0, 255], [255, 255, 0, 255], [245, 245, 220, 255],
    #               [0, 255, 0, 255], [0, 128, 0, 255], [245, 255, 250, 255], [0, 128, 128, 255], [0, 255, 255, 255],
    #               [0, 0, 128, 255], [230, 230, 250, 255], [255, 0, 255, 255], [205, 133, 63, 255]]
    #
    #     floatColors = []
    #     for color in colors:
    #         floatColors.append([x / 255 for x in color])
    #     mesh = self.mesh
    #     mesh.request_face_colors()
    #     propColorDict = {}
    #     self.legendValues.clear()
    #     self.legendColors.clear()
    #     self.drawLegend = False
    #     for el in self.element2Face:
    #         idProp=propDict[el]
    #         nuc=len(propColorDict)
    #         indexColor = 0
    #         if idProp in propColorDict:
    #             indexColor=propColorDict[idProp]
    #         else:
    #             propColorDict[idProp]=nuc
    #             indexColor=nuc
    #             self.legendValues.append(str(idProp))
    #             self.legendColors.append(floatColors[indexColor])
    #         for fh in self.element2Face[el]:
    #             mesh.set_color(fh, floatColors[indexColor])
    #         pass
    #     if len(self.legendValues)> 0:
    #         self.drawLegend=True
    #     pass