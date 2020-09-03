from geofem import GeoFEM,Node,StiffQuadElement, Material, BeamProperty,StiffLayoutProperty,PlateProperty
from geofem import StiffTriaElement, BeamElement, RodElement,Element
import  numpy as np
import uuid
import xml.etree.ElementTree as ET
def getFloat(item, key):
    return  float(item.attrib[key])
def getInt(item,key):
    return  int(item.attrib[key])
def getGuid(item,key):
    return  uuid.UUID(item.attrib[key])
def getStr(item, key):
    return  (item.attrib[key])

def getNPVector(item, key):
    npVec=  (np.array(item.attrib[key].split(','))).astype(np.float)
    return  npVec


class MaestroXML:

    def __init__(self, xmlPath):
        self.xmlpath = xmlPath
    def addElement(self,fem:GeoFEM,elem:Element,id,idProp,nodeIds):
        sNodeIds = nodeIds.split(' ')
        elem.init(id)
        elem.property = fem.getProperty(idProp)
        for idNod in sNodeIds:
            node = fem.getNode(int(idNod))
            elem.addNode(node)
        fem.addElement(elem)
    def readModelToGeoFEM(self,fem:GeoFEM):
        fname = self.xmlpath
        tree = ET.parse(fname)
        root = tree.getroot()
        jobctrl = 0
        nodes = 0
        elements = 0
        properties = 0
        stiflayouts = 0
        materials=0
        bars=0
        plates=0
        groups = 0
        evaluation = 0
        units =0
        for child in root:
            for ch in child:
                #print(ch.tag, ch.attrib)
                sname=(ch.attrib['sName'])
                if sname == 'Maestro':
                    units = ch[0][0]
                elif sname=='FeModel':
                    jobctrl=ch[0][0][0]
                    nodes = ch[0][0][1]
                    elements = ch[0][0][2]
                    properties = ch[0][0][3]
                    #stiflayouts = ch[0][0][4]
                    groups = ch[0][0][5]
                    #evaluation = ch[0][0][5][0]
                if sname=='FeStifLayout':
                    stiflayouts = ch[0]
                if sname=='MaterialIso':
                    materials = ch[0]
                if sname == 'FePropBar':
                    bars = ch[0]
                if sname == 'FePropPlate':
                    plates = ch[0]


        dbars = {}
        dstifflayprops ={}
        dmaterials = {}
        dproperties={}
        for item in units:
            unit = Material()
            if item.attrib['sTag']== "length":
                fem.units.name_length=item.attrib['sLabel']
                fem.units.user2si_length = getFloat(item,"dUserToSI")
            if item.attrib['sTag']== "force":
                fem.units.name_force=item.attrib['sLabel']
                fem.units.user2si_force = getFloat(item,"dUserToSI")
        for item in materials:
            material = Material()
            material.init(getInt(item, 'iTag'),item.attrib['sName'])
            guid=getGuid(item,'idSelf')
            dmaterials[guid]=material
            fem.addMaterial(material)
            material.E=getFloat(item, 'dYoungs')
            material.ni = getFloat(item, 'dPoisson')
            material.ReH = getFloat(item, 'dYield')
        for item in bars:
            barprop = BeamProperty()
            barprop.init(getInt(item, 'iId'),item.attrib['sName'])
            guid=getGuid(item,'idSelf')
            dbars[guid]=barprop
            fem.addProperty(barprop)
            barprop.hw = getFloat(item,"dWebHeight")
            barprop.tw = getFloat(item,"dWebThick");
            barprop.bf = getFloat(item,"dFlangeBreadth");
            barprop.tf = getFloat(item,"dFlangeThick");
            barprop.secType = getStr(item,"sSecType");
            if item[0].tag == 'Moniker':
                guid=getGuid(item[0],'idrefObj')
                barprop.material=dmaterials[guid]
        for item in stiflayouts:
            stifflay = StiffLayoutProperty()
            stifflay.init(getInt(item, 'iTag'), item.attrib['sName'])
            fem.addStiffLayout(stifflay)
            if item[0].tag == 'Moniker':
                guid = getGuid(item[0], 'idrefObj')
                stifflay.beam = dbars[guid]

        for item in plates:
            plate = PlateProperty()
            plate.init(getInt(item, 'iId'),item.attrib['sName'])
            fem.addProperty(plate)
            plate.tp = getFloat(item[0][0],"dThick")
            if item[0][0][0].tag == 'Moniker':
                guid=getGuid(item[0][0][0],'idrefObj')
                plate.material=dmaterials[guid]

        for item in nodes:
            node=Node()
            node.init(getInt(item, 'iId'), getFloat(item, 'dX'), getFloat(item, 'dY'), getFloat(item, 'dZ'))
            fem.addNode(node)
        for item in elements:
            if  item.tag == 'Quad':
                elem = StiffQuadElement()
                elem.layout     = fem.getStiffLayout(getInt(item, "iIdLayout"))
                self.addElement(fem, elem, getInt(item, 'iId'), getInt(item, "iIdProp"), getStr(item, 'sNodeIds'))
            elif  item.tag == 'Tri':
                elem = StiffTriaElement()
                elem.layout     = fem.getStiffLayout(getInt(item, "iIdLayout"))
                self.addElement(fem, elem, getInt(item, 'iId'), getInt(item, "iIdProp"), getStr(item, 'sNodeIds'))
            elif  item.tag == 'Bar':
                elem = BeamElement()
                sNodeIds = item.attrib['sNodeIds'].split(' ')
                elem.wo= getNPVector(item,'sWebVec')
                self.addElement(fem, elem, getInt(item, 'iId'), getInt(item, "iIdProp"),getStr(item, 'sNodeIds'))
            elif  item.tag == 'Rod':
                elem = RodElement()
                self.addElement(fem, elem, getInt(item, 'iId'), getInt(item, "iIdProp"), getStr(item, 'sNodeIds'))

        # for item in groups:
        #     print(item.tag)
            # print(item.attrib)
        # for item in evaluation:
        #     print(item.tag)
            # print(item.attrib)

pass