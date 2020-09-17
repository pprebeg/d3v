import openmesh as om
from enum import Enum
import numpy as np
from geometry import Geometry
import pathlib

class ViewType(Enum):
    """!
    ViewType
    """

    constant_color = 0
    face_colors = 1
    face_vertex_colors = 2

class MeshControl():
    def __init__(self):
        self.viewtype=ViewType.constant_color
        self.useviewtreshold=False
        self.uppertreshold=0
        self.lowertreshold = 0
        pass
    def getUpperTresholdColor(self):
        return [0.1,0.1,0.1,1.0]
    def getLowerTresholdColor(self):
        return [0.3,0.3,0.3,1.0]
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
        self.face_value=result.get(self.id)
        return self.face_value

    def getColorForFaceResult(self,fun_getcolor, minvalue, maxvalue):
        # define color

        if fun_getcolor != None:
            color = fun_getcolor(self.face_value, minvalue, maxvalue)
        else:
            color=[]
        return color

    def getColorForFaceVertexResult(self,fun_getcolor, minvalue, maxvalue):
        # define color
        colors = []
        if fun_getcolor != None:
            for i in range(len(self.vertex_based_values)):
                colors.append(fun_getcolor(self.vertex_based_values[i], minvalue, maxvalue))

    def addNode(self,node):
        self.nodes.append(node)
        self.vertex_based_values.append(0)
    def init(self,id):
        super().init(id)
    def updateMesh(self,mesh:om.TriMesh,mc:MeshControl, const_color = [0.4, 1.0, 1.0, 1.0],
                   fun_getcolor=None,minvalue = 0,maxvalue = 1 ):
        pass

    def onTPLValue(self):
        self.face_value = 0
        return self.face_value

    def onMatID(self):
        self.face_value = self.property.material.id
        return self.face_value

    def onPropID(self):
        self.face_value = self.property.id
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

    def updateMesh(self,mesh:om.TriMesh,mc:MeshControl, const_color = [0.4, 1.0, 1.0, 1.0],
                   fun_getcolor=None,minvalue = 0,maxvalue = 1 ):
        color = const_color
        vhandle = []
        handleColorIndex = []
        fhs = []

        hw=self.property.hw
        tw=self.property.tw
        bf=self.property.bf
        tf = self.property.tf
        x = self.nodes[0].p - self.nodes[1].p
        y = self.wo
        v = np.cross(x, y)
        # z = self.nodes[0].p + self.wo*hw - (v*bf*0.5)

        vhandle.append(mesh.add_vertex(self.nodes[0].p)) #0 točka 1
        handleColorIndex.append(0)
        vhandle.append(mesh.add_vertex(self.nodes[1].p)) #1 točka 2
        handleColorIndex.append(1)

        data = self.nodes[0].p + self.wo*hw  #2 točka 3
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(0)
        data = self.nodes[1].p + self.wo * hw  # 3 točka 4
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(1)
        data = self.nodes[0].p + self.wo*hw - (v*bf*0.5)  #4 točka 5
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(0)
        data = self.nodes[1].p + self.wo*hw - (v*bf*0.5)  #5 točka 6
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(1)
        data = self.nodes[1].p + self.wo*hw + (v*bf*0.5)  #6 točka 7
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(1)
        data = self.nodes[0].p + self.wo*hw + (v*bf*0.5)  #7 točka 8
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(0)

        # data = np.array([self.nodes[0].p]) + np.array(self.wo)*self.hw - np.cross(np.array([self.nodes[1].p]),np.array(self.wo))*(self.bf*0.5) #4 točak 5
        # vhandle.append(mesh.add_vertex(data))

        fhs.append(mesh.add_face(vhandle[0], vhandle[1], vhandle[2])) #1-2-3
        fhs.append(mesh.add_face(vhandle[2], vhandle[1], vhandle[3])) #3-2-4
        fhs.append(mesh.add_face(vhandle[4], vhandle[5], vhandle[6])) #5-6-7
        fhs.append(mesh.add_face(vhandle[4], vhandle[6], vhandle[7])) # 5-7-8

        # define color
        if fun_getcolor != None:
            if mc.viewtype == ViewType.face_colors:
                color = self.getColorForFaceResult(fun_getcolor, minvalue, maxvalue)
                for fh in fhs:
                    mesh.set_color(fh, color)
            if mc.viewtype == ViewType.face_vertex_colors:
                colors = self.getColorForFaceVertexResult(fun_getcolor, minvalue, maxvalue)
                for ivh in range(len(vhandle)):
                    mesh.set_color(vhandle[ivh], color[handleColorIndex[ivh]])

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
    def updateMesh(self,mesh:om.TriMesh,mc:MeshControl, const_color = [0.4, 1.0, 1.0, 1.0],
                   fun_getcolor=None,minvalue = 0,maxvalue = 1 ):

        color = const_color
        vhandle = []
        fhs = []
        data = np.array([self.nodes[0].p[0], self.nodes[0].p[1], self.nodes[0].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[1].p[0], self.nodes[1].p[1], self.nodes[1].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[2].p[0], self.nodes[2].p[1], self.nodes[2].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[3].p[0], self.nodes[3].p[1], self.nodes[3].p[2]])
        vhandle.append(mesh.add_vertex(data))

        fhs.append(mesh.add_face(vhandle[0], vhandle[1], vhandle[2]))
        fhs.append(mesh.add_face(vhandle[0], vhandle[2], vhandle[3]))

        # define color
        if fun_getcolor != None:
            if mc.viewtype == ViewType.face_colors:
                color =self.getColorForFaceResult(fun_getcolor, minvalue, maxvalue)
                for fh in fhs:
                    mesh.set_color(fh, color)
            elif mc.viewtype == ViewType.face_vertex_colors:
                colors =self.getColorForFaceVertexResult(fun_getcolor, minvalue, maxvalue)
                for ivh in range(len(vhandle)):
                    mesh.set_color(vhandle[ivh], color[ivh])




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
        self.numDiffValues =0
        self.valueIndexColor = {}

        self.drawLegend = False
        self.legendValues = []
        self.legendColors = []
        self.legendTitle = ""
        self.fixed_color_list = self.initColorList()
        self.max_fixed_colors= len(self.fixed_color_list)

        pass
    def initColorList(self):
        colors = [[0, 0, 255,255],[128, 0, 128,255],[222, 184, 135,255],[255, 165, 0,255],[0, 255, 0,255],
                  [ 0, 128, 0,255],[128, 0, 0,255],[255, 0, 0,255],[255, 192, 203,255],[222, 184, 135,255],
                  [255, 165, 0,255],[255, 127, 80,255],[128, 128, 0,255],[255, 255, 0,255],[245, 245, 220,255],
                  [0, 255, 0,255],[ 0, 128, 0,255],[245, 255, 250,255],[0, 128, 128,255],[0, 255, 255,255],
                  [0, 0, 128,255],[230, 230, 250,255],[255, 0, 255,255],[205, 133, 63,255]]
        floatColors = []
        for color in colors:
            floatColors.append([x / 255 for x in color])
        return floatColors

    def prepareModelForVisualization(self,key):
        self.minValue = float("inf")
        self.maxValue = float("-inf")
        self.numDiffValues = 0
        self.valueIndexColor.clear()

        fatrib= self.attrib_val_functions.get(key)

        if fatrib == None:
            self.doResultValue(key)
        else:
            fatrib(key)
        self.mc.lowertreshold=self.minValue
        self.mc.uppertreshold=self.maxValue
        self.mc.viewtype = ViewType.face_colors
        self.legendTitle=key
        self.regenerateusingcolor()

    def getContinuousColor(self, v, vmin, vmax):
        color = [1.0, 1.0, 1.0, 1.0]
        if v > self.mc.uppertreshold:
            return self.mc.getUpperTresholdColor()
        elif v < self.mc.lowertreshold:
            return self.mc.getLowerTresholdColor()
        vmin   = max(vmin,self.mc.lowertreshold)
        vmax = min(vmax, self.mc.uppertreshold)

        if v < vmin:
            v = vmin
        if v > vmax:
            v = vmax
        dv = vmax - vmin

        if (v < (vmin + 0.25 * dv)):
            color[0] = 0
            color[1] = 4 * (v - vmin) / dv
        elif (v < (vmin + 0.5 * dv)):
            color[0] = 0
            color[2] = 1 + 4 * (vmin + 0.25 * dv - v) / dv
        elif (v < (vmin + 0.75 * dv)):
            color[0] = 4 * (v - vmin - 0.5 * dv) / dv
            color[2] = 0
        else:
            color[1] = 1 + 4 * (vmin + 0.75 * dv - v) / dv
            color[2] = 0
        return color

    def getColorFromList(self, v, vmin, vmax):
        if v > self.mc.uppertreshold:
            return self.mc.getUpperTresholdColor()
        elif v < self.mc.lowertreshold:
            return self.mc.getLowerTresholdColor()

        index=self.getValueColorIndex(v)
        color = self.fixed_color_list[index]
        return color


    def prepContColorLegend(self,fun_getcolor,minVal, maxVal,nColor):
        self.legendValues.clear()
        self.legendColors.clear()
        self.drawLegend = True
        minVal = max(minVal, self.mc.lowertreshold)
        maxVal = min(maxVal, self.mc.uppertreshold)
        legendValues=np.linspace(minVal,maxVal,nColor)
        for x in legendValues:
            self.legendValues.append(f"{x:.4g}")
        for val in legendValues:
            color = fun_getcolor(val, minVal, maxVal)
            self.legendColors.append(color)

    def prepListColorLegend(self,fun_getcolor):
        self.legendValues.clear()
        self.legendColors.clear()
        self.drawLegend = True
        for key,index in self.valueIndexColor.items():
            if  key < self.mc.lowertreshold or key > self.mc.uppertreshold:
                continue
            self.legendValues.append(f"{key:.4g}")
            color = fun_getcolor(key, 0, self.numDiffValues)
            self.legendColors.append(color)

    # endregion
    # region Attribute Value Functions

    def populateAtribValFunctionsDictionary(self):
        self.addAttValFunc('TPL', self.doTPLValue)
        self.addAttValFunc('Material ID', self.doMaterialIDValue)
        self.addAttValFunc('Property ID', self.doPropertyIDValue)

    def addAttValFunc(self, key, f):
        self.attrib_val_functions[key] = f

    def doTPLValue(self,key):
        for el in self.elements.values():
            val= el.onTPLValue()
            self.checkMinMax(val)


    def doMaterialIDValue(self,key):
        for el in self.elements.values():
            val = el.onMatID()
            self.checkMinMax(val)

    def doPropertyIDValue(self, key):
        for el in self.elements.values():
            val = el.onPropID()
            self.checkMinMax(val)


    def doResultValue(self,key):
        self.setValueToItemResults(key)

    # endregion


    def checkMinMax(self,val):
        if val > self.maxValue:
            self.maxValue = val
        if val < self.minValue:
            self.minValue = val

        if self.numDiffValues <= self.max_fixed_colors:
            index = self.getValueColorIndex(val)
            if index == None:
                self.valueIndexColor[val]=self.numDiffValues
                self.numDiffValues=self.numDiffValues+1


    def getValueColorIndex(self,value):
        return self.valueIndexColor.get(value)


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
            for key, el in self.elements.items():
                val = el.setFaceValueUsingElementID(result.feres)
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

        self.mc.viewtype = ViewType.constant_color
       # mesh.request_face_colors()
        for el in self.elements.values():
            el.updateMesh(mesh,self.mc)
        pass
        self.mesh = mesh

    def regenerateusingcolor(self):
        fun_getcolor = self.getColorFromList
        if self.numDiffValues > self.max_fixed_colors:
            fun_getcolor= self.getContinuousColor

        mesh= om.TriMesh()
        if self.mc.viewtype == ViewType.constant_color or self.mc.viewtype == ViewType.face_colors:
            #mesh.release_vertex_colors()
            mesh.request_face_colors()
        elif self.mc.viewtype == ViewType.face_vertex_colors:
            #mesh.release_face_colors()
            mesh.request_vertex_colors()

        const_color = [0.4, 1.0, 1.0, 1.0]

        for el in self.elements.values():
            el.updateMesh(mesh,self.mc, const_color,fun_getcolor,self.minValue,self.maxValue)

        self.mesh = mesh

        if self.numDiffValues > self.max_fixed_colors:
            self.prepContColorLegend(fun_getcolor,self.minValue, self.maxValue, 12)
        else:
            self.prepListColorLegend(fun_getcolor)



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
        pass

    def appendValue(self, key, value):
        resultList = self.results.setdefault(key, [])
        resultList.append(value)

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

class GeneralResultsTableModel(Result):
    def __init__(self, name):
        super().__init__(name)
        self.data =[]
        self.column_names =[]
        pass

    def appendName(self, name: str):
        self.column_names.append(name)

    def addRow(self, values: list):
        self.data.append(values)

    def getRowValues(self, row_index):
        return self.data[row_index]

    def getValue(self, row_index, column_index: int):
        return self.data[row_index][column_index]

class ElementResult(Result):
    def __init__(self, name):
        super().__init__(name)
        self.feres = {}
        pass

    def getValue(self, feID):
        return self.feres.get(feID)

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
        self.modres={}
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

        self.readCSDFile(abspath_hoggCSD,False)
        self.readCSDFile(abspath_saggCSD,True)

        abspath_hogg = abspath + 'LUSAhogg.OUT'
        abspath_sagg = abspath + 'LUSAsagg.OUT'

        iterationResults = GeneralResultsTableModel('Lusa iteration results Sagg')
        self.modres[iterationResults.name]=iterationResults


        iterationResults.appendName('CycleNo Sagg')
        iterationResults.appendName('Moment Sagg, kNm')
        iterationResults.appendName('Curvature Sagg, 1/m')
        iterationResults.appendName('y_NA Sagg, m')

        self.readMainLusaFile(abspath_sagg,True,iterationResults)

        iterationResults = GeneralResultsTableModel('Lusa iteration results Hogg')
        self.modres[iterationResults.name] = iterationResults
        iterationResults.appendName('CycleNo Hogg')
        iterationResults.appendName('Moment Hogg, kNm')
        iterationResults.appendName('Curvature Hogg, 1/m')
        iterationResults.appendName('y_NA Hogg, m')

        self.readMainLusaFile(abspath_hogg, False, iterationResults)

        pass

    def readMainLusaFile(self, path, isSagg, tableResult:GeneralResultsTableModel):
        file = pathlib.Path(path)
        if not file.exists():
            return
        f = open(path, "r")
        nlines2skip = 0
        isCycleData = False

        for line in f:
            if nlines2skip > 0:
                nlines2skip = nlines2skip - 1
                continue
            line = ' '.join(line.split())
            if line.startswith('*'):
                continue
            if line == "" or line == " ":
                continue

            if 'HULL MODULE RESPONSE DATA' in line:
                nlines2skip = 4
                isCycleData = True
                continue
            if 'ULTIMATE CAPACITY IS' in line:
                f.close()
#                if isSagg:
#                    tableResult.data.reverse()
                return
            sline = line.split(" ")
            if len(sline) == 0:
                continue
            if isCycleData:
                if len(sline) == 4:
                    rowValues=[0]*4
                    rowValues[0]= int(sline[0])
                    rowValues[1] = float(sline[1])
                    rowValues[2] = float(sline[2])
                    rowValues[3] = float(sline[3])
                    tableResult.addRow(rowValues)

        f.close()

    def readCSDFile(self,path,isSagg):
        file=pathlib.Path(path)
        if not file.exists():
            return
        f = open(path, "r")
        nlines2skip=0
        isSPCdata=False
        isGPCdata = False
        isHCdata = False
        if isSagg:
            collapse_stress = ElementResult('Collapse Stress Sagg')
            collapse_mod = ElementResult('Collapse Mod Sagg')
            collapse_cycle = ElementResult('Collapse Cycle Sagg')
        else:
            collapse_stress = ElementResult('Collapse Stress Hogg')
            collapse_mod = ElementResult('Collapse Mod Hogg')
            collapse_cycle = ElementResult('Collapse Cycle Hogg')

        self.lers[collapse_stress.name]=collapse_stress
        self.lers[collapse_mod.name] = collapse_mod
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
                    elIDsPlate = self.mas.getPlateElemForStrake(strakeNo)
                    if elIDs == None:
                        elIDs = elIDsPlate
                    elif elIDsPlate != None:
                        for el in elIDsPlate:
                            elIDs.append(el)
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