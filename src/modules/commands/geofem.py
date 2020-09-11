import openmesh as om
from enum import Enum
import numpy as np
from geometry import Geometry
import pathlib

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
        self.useviewtreshold=False
        self.uppertreshold=0
        self.lowertreshold = 0
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
        self.face_value=0
        self.vertex_based_values=[]
        self.value_name=""
        pass

    def setFaceValueUsingElementID(self,result:dict):
        val= result.get(self.id)
        return val

    def addNode(self,node):
        self.nodes.append(node)
        self.vertex_based_values.append(0)
    def init(self,id):
        super().init(id)
    def updateMesh(self,mesh:om.TriMesh,mc:MeshControl):
        pass

    def onTPLValue(self):
        self.face_value = 0
        return self.face_value


class RodElement(Element):
    def __init__(self):
        super().__init__()
        pass


class BeamElement(Element):
    def __init__(self):
        super().__init__()
        self.wo = np.array([0, 0, 0]) #web orientation

    def onTPLValue(self):
        self.face_value = self.property.tw
        return self.face_value

    def updateMesh(self,mesh:om.TriMesh,mc:MeshControl):
        print("čvor 1")
        # print(self.nodes[0].p)
        # print(self.nodes[0].p[0])
        # print(self.nodes[0].p[1])
        # print(self.nodes[0].p[2])
        # print("čvor 1a")
        # print(self.nodes[0].p[0])
        # print(self.nodes[0].p[1])
        # print(self.nodes[0].p[2]-1)
        # print("čvor 2")
        # print(self.nodes[1].p)

        vhandle = []

        hw=self.property.hw
        tw=self.property.tw
        bf=self.property.bf
        tf = self.property.tf
        x = self.nodes[0].p - self.nodes[1].p
        y = self.wo
        v = np.cross(x, y)
        # z = self.nodes[0].p + self.wo*hw - (v*bf*0.5)
        print(self.nodes[0].p)
        print("x")
        print(x)
        # print("y")
        # print(y)
        print("v")
        print(v)
        print("wo")
        print(self.wo)
        print("hw")
        print(hw)
        print("wo*hw")
        print(self.wo*hw)
        # print("z")
        # print(z)






        vhandle.append(mesh.add_vertex(self.nodes[0].p)) #0 točka 1
        vhandle.append(mesh.add_vertex(self.nodes[1].p)) #1 točka 2


        data = self.nodes[0].p + self.wo*hw  #2 točka 3
        vhandle.append(mesh.add_vertex(data))
        data = self.nodes[1].p + self.wo * hw  # 3 točka 4
        vhandle.append(mesh.add_vertex(data))
        data = self.nodes[0].p + self.wo*hw - (v*bf*0.5)  #4 točka 5
        vhandle.append(mesh.add_vertex(data))
        data = self.nodes[1].p + self.wo*hw - (v*bf*0.5)  #5 točka 6
        vhandle.append(mesh.add_vertex(data))
        data = self.nodes[1].p + self.wo*hw + (v*bf*0.5)  #6 točka 7
        vhandle.append(mesh.add_vertex(data))
        data = self.nodes[0].p + self.wo*hw + (v*bf*0.5)  #7 točka 8
        vhandle.append(mesh.add_vertex(data))

        # data = np.array([self.nodes[0].p]) + np.array(self.wo)*self.hw - np.cross(np.array([self.nodes[1].p]),np.array(self.wo))*(self.bf*0.5) #4 točak 5
        # vhandle.append(mesh.add_vertex(data))

        fh0 = mesh.add_face(vhandle[0], vhandle[1], vhandle[2]) #1-2-3
        fh1 = mesh.add_face(vhandle[2], vhandle[1], vhandle[3]) #3-2-4
        fh1 = mesh.add_face(vhandle[4], vhandle[5], vhandle[6]) #5-6-7
        fh1 = mesh.add_face(vhandle[4], vhandle[6], vhandle[7])  # 5-7-8



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

    def onTPLValue(self):
        self.face_value = self.property.tp
        return self.face_value

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

    def onTPLValue(self):
        self.face_value = self.property.tp
        return self.face_value


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
class MaestroElementAssociation():
    def __init__(self):
        self.strakeGirder2fe = {}
        self.strakePlate2fe = {}
        self.strakeFrame2fe = {}
        self.endPoint2node = {}
        self.endPointStrakes = {}

    def addStrakePlate(self,key:int,elList:[]):
        self.strakePlate2fe[key]=elList

    def addStrakeGirder(self,key:int,elList:[]):
        self.strakeGirder2fe[key]=elList

    def addStrakeFrame(self,key:int,elList:[]):
        self.strakeFrame2fe[key]=elList

    def addEndPointNode(self,key:int,nodeFeTag:int):
        feTagList = self.endPoint2node.setdefault(key,[])
        feTagList.append(nodeFeTag)

    def addEndPointStrake(self, key:int, strakeID:int):
        strakeList = self.endPointStrakes.setdefault(key,[])
        strakeList.append(strakeID)

    def getPlateElsForEnpoint(self, endpointID):
        connectedElements=[]
        strakeList = self.endPointStrakes.get(endpointID)
        if strakeList is not None:
            for idStrake in strakeList:
                for idEl in self.strakePlate2fe[idStrake]:
                    connectedElements.append(idEl)
        return connectedElements

    def getPlateElemForStrake(self,strakeID):
        return self.strakePlate2fe.get(strakeID)

    def getGirderBeamElemForStrake(self,strakeID):
        return self.strakeGirder2fe.get(strakeID)


class LusaElementAssociation():
    def __init__(self):
        self.spc2fe = {}
        self.plate2fe = {}
        self.hc2fe = {}

    def addPlate(self, key:int, strakeID):
        self.plate2fe[key]=strakeID
    def addSPC(self,key:int,strakeID):
        self.spc2fe[key]=strakeID
    def addHC(self,key:int,strakeID):
        self.hc2fe[key]=strakeID

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
        self.mas = MaestroElementAssociation()
        self.element_results= {}
        self.element_vertex_results = {}
        self.vertex_results = {}
        self.model_results = {}
        self.result_name=""
        self.is_node_result=False
        self.is_element_result = False
        self.attrib_val_functions ={}
        self.populateAtribValFunctionsDictionary()
        self.minValue=0
        self.maxValue=0
        pass

    def prepareModelForVisualization(self,key):
        self.minValue = float("inf")
        self.maxValue = float("-inf")
        fatrib= self.attrib_val_functions.get(key)

        if fatrib == None:
            self.doResultValue(key)
        else:
            fatrib()


    # region Attribute Value Functions

    def populateAtribValFunctionsDictionary(self):
        self.addAttValFunc('TPL', self.doTPLValue)
        self.addAttValFunc('Material ID', self.doMaterialIDValue)
        self.addAttValFunc('Property ID', self.doPropertyIDValue)

    def addAttValFunc(self, key, f):
        self.attrib_val_functions[key] = f

    def doTPLValue(self,key):
        for el in self.elements:
            val= el.onTPLValue(self)
            self.checkMinMax(val)


    def doMaterialIDValue(self,key):
        pass

    def doPropertyIDValue(self, key):
        pass


    def doResultValue(self,key):
        self.setValueToItemResults(key)

    # endregion


    def checkMinMax(self,val):
        if val > self.maxValue:
            self.maxValue = val
        if val < self.minValue:
            self.minValue = val

    def isElementFaceResult(self):
        return self.is_element_result and (not self.is_node_result)

    def isElementNodeResult(self):
        return self.is_element_result and self.is_node_result

    def isNodeResult(self):
        return (not self.is_element_result) and self.is_node_result

    def setValueToItemResults(self, resultName):
        result = self.element_results.get(resultName)
        if result != None:
            self.is_element_result= True
            self.is_node_result = False
            for el in self.elements:
                val = el.setFaceValueUsingElementID(result)
                self.checkMinMax(val)
            return

        result = self.element_vertex_results.get(resultName)
        if result != None:
            self.is_element_result= True
            self.is_node_result = True
            return

        result = self.vertex_results.get(resultName)
        if result != None:
            self.is_element_result = False
            self.is_node_result = True
            return


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
    def setResultValuesOnElements(self):
        pass



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


class Result():
    def __init__(self, name):
        self.name = name
        pass


class GeneralResultsDictionary(Result):
    def __init__(self, name):
        super().__init__(name)
        self.results = {}
        self.valueNames = []
        self.keyName = ""
        pass

    def initilizeResultDictionary(self, keyName, valueNames):
        self.keyName = keyName
        self.valueNames = [valueNames]

    def addValues(self, key, values: []):
        self.results[key] = values

    def addValues2(self, key, values: list):
        self.results[key] = values.copy()

    def getValues(self, key):
        return self.results[key]

    def getValue(self, key, index: int):
        return self.results[key][index]

    def appdendListwithResultData(self, x: list, y: list, iresx: int, iresy: int):
        for res in self.results.values():
            x.append(res[iresx])
            y.append(res[iresy])

    def appdendListwithKeyPairedResultData(self, x: list, y: list, iresy: int):
        for key, res in self.results:
            x.append(key)
            y.append(res[iresy])


class ElementResult(Result):
    def __init__(self, name):
        super().__init__(name)
        self.feres = {}
        pass

    def getValue(self, feID):
        return self.feres[feID]

    def setValue(self, feID, value):
        self.feres[feID]=value

    def keyExist(self,feID):
        return feID in self.feres


class FEMModelResults:
    def __init__(self, name):
        self.name = name
        pass

    def readOutput(self, path):
        pass

    def setResultsToModel(self, fem: GeoFEM):
        pass


class LusaResults(FEMModelResults):
    def __init__(self, name, mas:MaestroElementAssociation):
        super().__init__(name)
        self.las = LusaElementAssociation()
        self.mas=mas
        self.lers = {} #Lusa element results
        pass

    def readOutput(self, path):
        abspath1 = '\\'.join(path.split('\\')[0:-1])
        abspath2 = '/'.join(path.split('/')[0:-1])
        if len(abspath2) > len(abspath1):
            abspath=abspath2 + '/'
        else:
            abspath=abspath1 + '\\'

        abspath_hoggCSD = abspath + 'LUSAhoggCSD.OUT'
        abspath_saggCSD=  abspath + 'LUSAsaggCSD.OUT'

        self.readCSDFile(abspath_hoggCSD)
        self.readCSDFile(abspath_saggCSD)
        pass

    def readCSDFile(self,path):
        file=pathlib.Path(path)
        if not file.exists():
            return
        f = open(path, "r")
        nlines2skip=0
        isSPCdata=False
        isGPCdata = False
        isHCdata = False
        collapse_stress = ElementResult('Collapse Stress')
        self.lers[collapse_stress.name]=collapse_stress
        collapse_mod = ElementResult('Collapse Mod')
        self.lers[collapse_mod.name] = collapse_mod
        collapse_cycle = ElementResult('Collapse Cycle')
        self.lers[collapse_cycle.name] = collapse_cycle
        for line in f:
            if  nlines2skip > 0:
                nlines2skip=nlines2skip-1
                continue
            line = ' '.join(line.split())
            if line.startswith('*'):
                continue
            if line == "" or line == " ":
                continue

            if 'Stiffener - Plate Combinations (SPCs)' in line:
                nlines2skip=4
                isSPCdata=True
                continue
            if 'Girder - Plate Combinations (GPCs)' in line:
                nlines2skip=4
                isGPCdata=True
                isSPCdata=False
                continue
            if 'Hard Corners (HCs)' in line:
                nlines2skip = 4
                isGPCdata = False
                isHCdata = True
                continue
            sline = line.split(" ")
            if  len(sline)== 0:
                continue
            el_no_lusa=-1
            if isSPCdata:
                if len(sline) > 5:
                    strakeNo    = int(sline[0])
                    el_no_lusa        = int(sline[1])
                    self.las.addPlate(el_no_lusa, strakeNo)

                    elIDs=self.mas.getPlateElemForStrake(strakeNo)
                    cc_new = int(sline[4])
                    for id_el in elIDs:
                        bAddNew=True
                        cc_old= collapse_cycle.getValue(id_el)
                        if cc_old != None and cc_old < cc_new:
                            pass
                        else:
                            collapse_stress.setValue(id_el, float(sline[2]))
                            collapse_mod.setValue(id_el, float(sline[3]))
                            collapse_cycle.setValue(id_el, cc_new)
            elif isGPCdata:
                if len(sline) > 5:
                    strakeNo    = int(sline[0])
                    el_no_lusa        = int(sline[1])
                    self.las.addSPC(el_no_lusa, strakeNo)

                    elIDs = self.mas.getGirderBeamElemForStrake(strakeNo)
                    cc_new = int(sline[4])
                    for id_el in elIDs:
                        bAddNew = True
                        cc_old = collapse_cycle.getValue(id_el)
                        if cc_old != None and cc_old < cc_new:
                            pass
                        else:
                            collapse_stress.setValue(id_el, float(sline[2]))
                            collapse_mod.setValue(id_el, float(sline[3]))
                            collapse_cycle.setValue(id_el, cc_new)
            elif isHCdata:
                if len(sline) > 5:
                    endPtNo    = int(sline[0])
                    el_no_lusa        = int(sline[1])
                    self.las.addHC(el_no_lusa, endPtNo)
                    elIDs = self.mas.getPlateElsForEnpoint(endPtNo)
                    cc_new = int(sline[4])
                    for id_el in elIDs:
                        bAddNew = True
                        cc_old = collapse_cycle.getValue(id_el)
                        if cc_old != None and cc_old < cc_new:
                            pass
                        else:
                            collapse_stress.setValue(id_el, float(sline[2]))
                            collapse_mod.setValue(id_el, float(sline[3]))
                            collapse_cycle.setValue(id_el, cc_new)

        f.close()
        pass

    def setResultsToModel(self, fem: GeoFEM):
        pass