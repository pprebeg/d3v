import random

from geometry import Geometry
import openmesh as om
import numpy as np



class OOFEM (Geometry):
    def __init__(self,fileName, guid = None):
        super().__init__(guid)
        self.drawLegend = False
        self.legendValues = []
        self.legendColors = []
        self.legendTitle=""
        self.filename=fileName
        self.all_face_handles = {} # key=faceHandle, value = element ID
        self.element2Face = {}     # key=element ID, value = [faceHandle1,...faceHandlex]
        self.elementId_color_dict = {
            1: 0,
            2: 3,
            3: 13,
            4: 7,
            5: 0,
            6: 11,
            7: 13,
            8: 15,
            9: 0,
            10: 0
        }
        self.crosssectdict={}
        self.materialdict = {}
        self.elementcrosssectdict = {}


    def genMesh(self):
        self.mesh = self.oofemmesh()

    def pokusaj2(self):
        #prvaTocka=[0,0,0]
        #drugaTocka=[1,0,0]
        #orient=[0,1,0] # orijentacija struka
        #hw=0.3
        #tw=0.03
        #bf=0.2
        #tf=0.02
        #bp=0.25
        #tp=0.025

        prvaTocka = []
        for i in range(0, 3):
            ele = float(input("Unesite koordinatu " + str(i + 1) + " prve točke: "))
            prvaTocka.append(ele)

        drugaTocka = []
        for i in range(0, 3):
            ele = float(input("Unesite koordinatu " + str(i + 1) + " druge točke: "))
            drugaTocka.append(ele)

        orijent = []
        for i in range(0, 3):
            ele = float(input("Unesite koordinatu " + str(i + 1) + " orijentacije struka: "))
            orijent.append(ele)

        hw = float(input("Unesite visinu struka u (m) : "))
        tw = float(input("Unesite debljinu struka u (m): "))
        bf = float(input("Unesite širinu flanže u (m): "))
        tf = float(input("Unesite debljinu flanže u (m): "))
        bp = float(input("Unesite širinu pokrova u (m): "))
        tp = float(input("Unesite debljinu pokrova u (m): "))


        return self.zadatak2(prvaTocka, drugaTocka,orijent,hw,tw,bf,tf,bp,tp)

        #return self.zadatak2(prvaTocka, drugaTocka, hw, tw, bf, tf, bp, tp)

    def zadatak2(self, prvaTocka, drugaTocka,orijent, hw,tw,bf,tf,bp,tp):
        mesh= om.TriMesh()

        vertexHandles = []
        faceHandles = []
        faceHandlesDict = {}

        Tocka1 = np.array(prvaTocka)
        vertexHandles.append(mesh.add_vertex(Tocka1))       #0
        Tocka2 = np.array(drugaTocka)
        vertexHandles.append(mesh.add_vertex(Tocka2))       #1
        Tocka1_ = Tocka1 + np.array([0, 0, 1])
        vertexHandles.append(mesh.add_vertex(Tocka1_))      #2
        Tocka2_ = Tocka2 + np.array([0, 0, 1])
        vertexHandles.append(mesh.add_vertex(Tocka2_))      #3
        Tocka1a = Tocka1 + np.array([0, 0, 2])
        vertexHandles.append(mesh.add_vertex(Tocka1a))      #4
        Tocka2a = Tocka2 + np.array([0, 0, 2])
        vertexHandles.append(mesh.add_vertex(Tocka2a))      #5

        x = np.array(Tocka2 - Tocka1)
        y = np.array(orijent)
        v = np.cross(x, y)
        v1 = v / np.linalg.norm(v)

        Tocka3 = Tocka2 + np.array(orijent)*hw
        vertexHandles.append(mesh.add_vertex(Tocka3))       #6
        Tocka4 = Tocka1 + np.array(orijent)*hw
        vertexHandles.append(mesh.add_vertex(Tocka4))       #7
        Tocka3_ = Tocka2_ + np.array(orijent)*hw
        vertexHandles.append(mesh.add_vertex(Tocka3_))      #8
        Tocka4_ = Tocka1_ + np.array(orijent)*hw
        vertexHandles.append(mesh.add_vertex(Tocka4_))      #9
        Tocka5_ = Tocka1_ + np.array(orijent)*hw - np.array(v1)*(bf * 0.5)
        vertexHandles.append(mesh.add_vertex(Tocka5_))      #10
        Tocka6_ = Tocka1_ + np.array(orijent)*hw + np.array(v1)*(bf * 0.5)
        vertexHandles.append(mesh.add_vertex(Tocka6_))      #11
        Tocka7_ = Tocka2_ + np.array(orijent)*hw + np.array(v1)*(bf * 0.5)
        vertexHandles.append(mesh.add_vertex(Tocka7_))      #12
        Tocka8_ = Tocka2_ + np.array(orijent)*hw - np.array(v1)*(bf * 0.5)
        vertexHandles.append(mesh.add_vertex(Tocka8_))      #13
        Tocka3a = Tocka2a + np.array(orijent)*hw
        vertexHandles.append(mesh.add_vertex(Tocka3a))      #14
        Tocka4a = Tocka1a + np.array(orijent)*hw
        vertexHandles.append(mesh.add_vertex(Tocka4a))      #15
        Tocka5a = Tocka1a + np.array(orijent)*hw - np.array(v1)*(bf * 0.5)
        vertexHandles.append(mesh.add_vertex(Tocka5a))      #16
        Tocka6a = Tocka1a + np.array(orijent)*hw + np.array(v1)*(bf * 0.5)
        vertexHandles.append(mesh.add_vertex(Tocka6a))      #17
        Tocka7a = Tocka2a + np.array(orijent)*hw + np.array(v1)*(bf * 0.5)
        vertexHandles.append(mesh.add_vertex(Tocka7a))      #18
        Tocka8a = Tocka2a + np.array(orijent)*hw - np.array(v1)*(bf * 0.5)
        vertexHandles.append(mesh.add_vertex(Tocka8a))      #19
        Tocka9a = Tocka1a - np.array(v1)*(bp * 0.5)
        vertexHandles.append(mesh.add_vertex(Tocka9a))      #20
        Tocka10a = Tocka1a + np.array(v1)*(bp * 0.5)
        vertexHandles.append(mesh.add_vertex(Tocka10a))     #21
        Tocka11a = Tocka2a + np.array(v1)*(bp * 0.5)
        vertexHandles.append(mesh.add_vertex(Tocka11a))     #22
        Tocka12a = Tocka2a - np.array(v1)*(bp * 0.5)
        vertexHandles.append(mesh.add_vertex(Tocka12a))     #23

        #faceHandles.append(mesh.add_face(vertexHandles[0],vertexHandles[1],vertexHandles[6]))
        faceHandlesDict[mesh.add_face(vertexHandles[0],vertexHandles[1],vertexHandles[6])] = 0
        #faceHandles.append(mesh.add_face(vertexHandles[6],vertexHandles[7],vertexHandles[0]))
        faceHandlesDict[mesh.add_face(vertexHandles[6], vertexHandles[7], vertexHandles[0])] = 6

        #faceHandles.append(mesh.add_face(vertexHandles[2], vertexHandles[3],vertexHandles[8]))
        faceHandlesDict[mesh.add_face(vertexHandles[2], vertexHandles[3], vertexHandles[8])] = 2
        #faceHandles.append(mesh.add_face(vertexHandles[8], vertexHandles[9],vertexHandles[2]))
        faceHandlesDict[mesh.add_face(vertexHandles[8], vertexHandles[9], vertexHandles[2])] = 8
        #faceHandles.append(mesh.add_face(vertexHandles[10],vertexHandles[11],vertexHandles[12]))
        faceHandlesDict[mesh.add_face(vertexHandles[10], vertexHandles[11], vertexHandles[12])] = 10
        #faceHandles.append(mesh.add_face(vertexHandles[12],vertexHandles[13],vertexHandles[10]))
        faceHandlesDict[mesh.add_face(vertexHandles[12], vertexHandles[13], vertexHandles[10])] = 12

        #faceHandles.append(mesh.add_face(vertexHandles[4],vertexHandles[5],vertexHandles[14]))
        faceHandlesDict[mesh.add_face(vertexHandles[4], vertexHandles[5], vertexHandles[14])] = 4
        #faceHandles.append(mesh.add_face(vertexHandles[14],vertexHandles[15], vertexHandles[4]))
        faceHandlesDict[mesh.add_face(vertexHandles[14], vertexHandles[15], vertexHandles[4])] = 14
        #faceHandles.append(mesh.add_face(vertexHandles[16],vertexHandles[17],vertexHandles[18]))
        faceHandlesDict[mesh.add_face(vertexHandles[16], vertexHandles[17], vertexHandles[18])] = 16
        #faceHandles.append(mesh.add_face(vertexHandles[18],vertexHandles[19],vertexHandles[16]))
        faceHandlesDict[mesh.add_face(vertexHandles[18], vertexHandles[19], vertexHandles[16])] = 18
        #faceHandles.append(mesh.add_face(vertexHandles[20],vertexHandles[21],vertexHandles[22]))
        faceHandlesDict[mesh.add_face(vertexHandles[20], vertexHandles[21], vertexHandles[22])] = 20
        #faceHandles.append(mesh.add_face(vertexHandles[22],vertexHandles[23],vertexHandles[20]))
        faceHandlesDict[mesh.add_face(vertexHandles[22], vertexHandles[23], vertexHandles[20])] = 22



        #print(vertexHandles)
        #print(faceHandles)
        print(faceHandlesDict)

        return mesh
        pass


    def oofemmesh(self):

        mesh = om.TriMesh()
        self.mesh=mesh
        f = open(self.filename, newline='')
        all_vertices = {}
        all_face_handles = {}
        plateElementTypes={"dktplate","mitc4shell","planestress2d","trplanestress2d","trplanestressrotallman","trplanestrrot"}

        for line in f:
            line = ' '.join(line.split())
            sline = line.split(" ")
            if len(sline) < 3:
                continue
            #print(sline[0].lower())
            if sline[0].lower() == "node":
                d = [float(sline[4]), float(sline[5]), float(sline[6])]
                all_vertices[sline[1]] = d
            elif  sline[0].lower() in plateElementTypes:
                id=int(sline[1])
                elementFaceHandles=[]
                self.element2Face[id]=elementFaceHandles
                self.elementcrosssectdict[id]=int(sline[-1])
                numNodes=int(sline[3])
                if numNodes >= 3:
                    vh_list = [mesh.add_vertex(all_vertices.get(sline[4])), mesh.add_vertex(all_vertices.get(sline[5])),
                               mesh.add_vertex(all_vertices.get(sline[6]))]
                    fh=mesh.add_face(vh_list)
                    self.all_face_handles[fh] = int(sline[1])
                    elementFaceHandles.append(fh)
                    pass
                    if numNodes == 4:
                        vh_list = [mesh.add_vertex(all_vertices.get(sline[6])), mesh.add_vertex(all_vertices.get(sline[7])),
                                   mesh.add_vertex(all_vertices.get(sline[4]))]
                        fh = mesh.add_face(vh_list)
                        self.all_face_handles[fh] = int(sline[1])
                        elementFaceHandles.append(fh)
                        pass
                else:
                    print ("unhandled type of the element")
                    pass
            elif sline[0].lower() == "simplecs":
                self.crosssectdict[int(sline[1])]=sline
            elif sline[0].lower() == "isole":
                self.materialdict[int(sline[1])] = sline
        f.close()
        #self.showFaceColorC()
        #self.showFaceColorP()
        #self.showVertexColor()
        return mesh
    def getPropIDforElID(self,elID):
        n=elID
        while n > 10:
            n=n-10
        propID = self.elementId_color_dict[n]
        return propID

    def getElResult(self,elID,numEl):
        retVal=0
        retVal=random.uniform(0, 1)
        return  retVal

    def showSolidColor(self):
        self.legendValues.clear()
        self.legendColors.clear()
        self.drawLegend = False

    def showCrossSectID(self):
        self.legendTitle="Cross section ID"
        self.showFaceColorP(self.elementcrosssectdict)

    def showMaterialID(self):
        self.legendTitle = "Material ID"
        propDict={}
        for el in self.elementcrosssectdict:
            csID= self.elementcrosssectdict[el]
            sline= self.crosssectdict[csID]
            propDict[el]=int(sline[-1])
        self.showFaceColorP(propDict)

    def showFaceColorP(self,propDict):
        colors = [[139,0,0,255],[220,20,60,255],[255,0,0,255],[255,20,147,255],[255,105,180,255],[255,192,203,255],[255,182,193,255],[0,100,0,255],[46,139,87,255],[143,188,143,255],[50,205,50,255],[0,255,0,255],[152,251,152,255],[0,0,139,255],[0,0,255,255],[65,105,225,255],[30,144,255,255],[0,191,255,255],[135,206,235,255],[173,216,230,255]]
        floatColors = []
        for color in colors:
            floatColors.append([x / 255 for x in color])
        mesh= self.mesh
        mesh.request_face_colors()
        propColorDict={}
        self.legendValues.clear()
        self.legendColors.clear()
        self.drawLegend = False
        for el in self.element2Face:
            idProp=propDict[el]
            nuc=len(propColorDict)
            indexColor = 0
            if idProp in propColorDict:
                indexColor=propColorDict[idProp]
            else:
                propColorDict[idProp]=nuc
                indexColor=nuc
                self.legendValues.append(str(idProp))
                self.legendColors.append(floatColors[indexColor])
            for fh in self.element2Face[el]:
                mesh.set_color(fh, floatColors[indexColor])
            pass
        if  len(self.legendValues)> 0:
            self.drawLegend=True
        pass

    def prepContColorLegend(self,minVal, maxVal,nColor):
        self.legendValues.clear()
        self.legendColors.clear()
        self.drawLegend = True
        self.legendValues=np.linspace(minVal,maxVal,nColor)
        for val in self.legendValues:
            color = self.getContinuousColor(val, minVal, maxVal)
            self.legendColors.append(color)


    def showFaceColorC(self,dictElVals):
        mesh= self.mesh
        mesh.request_face_colors()
        minVal=float('inf')
        maxVal = float('-inf')
        for el in dictElVals:
            val=dictElVals[el]
            if val > maxVal:
                maxVal=val
            if val < minVal:
                minVal = val
        index = 0
        self.prepContColorLegend(minVal,maxVal,10)
        for el in self.element2Face:
            index = index + 1
            val = dictElVals[el]
            for fh in self.element2Face[el]:
                color = self.getContinuousColor(val, minVal, maxVal)
                mesh.set_color(fh, color)
            pass
        pass

    def showVertexColor(self):
        mesh = self.mesh
        mesh.request_face_colors()
        for el in self.element2Face:
            for fh in self.element2Face[el]:
                iv=0
                for vh in mesh.fv(fh):
                    mesh.set_color(vh, [0, 0.5, 1, 1])
                    iv=iv+1
            pass
        pass

    def getContinuousColor(self, v, vmin, vmax):
        color = [1, 1, 1, 1]
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




