from iohandlers import IOHandler
from signals import Signals
from geometry import Geometry
import openmesh as om
import os
import numpy as np

class OOFEMImporter(IOHandler):
    def __init__(self):
        super().__init__()
        Signals.get().importGeometry.connect(self.importGeometry)

    def importGeometry(self, fileName):
        if len(fileName) < 1:
            return
        filename, file_extension = os.path.splitext(fileName)
        if file_extension != ".in":
            return
        oofem=OOFEM(fileName)
        #g = Geometry()
        #g.mesh =
        oofem.getmesh()
        Signals.get().geometryImported.emit(g)

    def getImportFormats(self):
        return []


def createIOHandler():
    return OOFEMImporter()

class OOFEM (Geometry):
    def __init__(self,fileName):
        self.filename=fileName
        self.all_face_handles = {} # key=faceHandle, value = element ID
        self.element2Face = {}     # key=element ID, value = [faceHandle1,...faceHandlex]


    def getmesh(self):
        #mesh=self.test()
        #mesh = self.test2()

        #mesh = self.pokusaj1()
        #mesh = self.pokusaj2()

        # mesh = self.oofemmesh()
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

    def test(self):
        mesh= om.TriMesh()
        vhandle = [0]*5
        data = np.array([0, 1, 0])
        vhandle[0] = mesh.add_vertex(data)
        #vhandle.append(mesh.add_vertex(data))
        data = np.array([1, 0, 0])
        vhandle[1] = mesh.add_vertex(data)
        data = np.array([2, 1, 0])
        vhandle[2] = mesh.add_vertex(data)
        data = np.array([0, -1, 0])
        vhandle[3] = mesh.add_vertex(data)
        data = np.array([2, -1, 0])
        vhandle[4] = mesh.add_vertex(data)

        fh0 = mesh.add_face(vhandle[0], vhandle[1], vhandle[2])
        fh1 = mesh.add_face(vhandle[1], vhandle[3], vhandle[4])
        fh2 = mesh.add_face(vhandle[0], vhandle[3], vhandle[1])

        vh_list = [vhandle[2], vhandle[1], vhandle[4]]
        fh3 = mesh.add_face(vh_list)


        return mesh
        pass
    def test2(self):
        mesh = om.TriMesh()
        # m --> min, M --> max
        # yapf: disable
        p0 = mesh.add_vertex([-1, -1, -1])
        p1 = mesh.add_vertex([-1, -1,  1])
        p2 = mesh.add_vertex([-1,  1, -1])
        p3 = mesh.add_vertex([-1,  1,  1])
        p4 = mesh.add_vertex([ 1, -1, -1])
        p5 = mesh.add_vertex([ 1, -1,  1])
        p6 = mesh.add_vertex([ 1,  1, -1])
        p7 = mesh.add_vertex([ 1,  1,  1])
        # yapf: enable

        mesh.add_face([p0, p6, p4])
        mesh.add_face([p0, p2, p6])
        mesh.add_face([p0, p4, p5])
        mesh.add_face([p0, p5, p1])
        mesh.add_face([p0, p3, p2])
        mesh.add_face([p0, p1, p3])
        mesh.add_face([p6, p2, p3])
        mesh.add_face([p6, p3, p7])
        mesh.add_face([p4, p7, p5])
        mesh.add_face([p4, p6, p7])
        mesh.add_face([p1, p5, p7])
        mesh.add_face([p1, p7, p3])

    def oofemmesh(self):
        mesh = om.TriMesh()
        self.mesh=mesh
        f = open(self.filename, newline='')
        all_vertices = {}
        all_face_handles = {}
        plateElementTypes={"dktplate","mitc4shell"}

        for line in f:
            line = ' '.join(line.split())
            sline = line.split(" ")
            if len(sline) < 3:
                continue
            #print(sline[0].lower())
            if sline[0] == "node":
                d = [float(sline[4]), float(sline[5]), float(sline[6])]
                all_vertices[sline[1]] = mesh.add_vertex(d)
            elif  sline[0].lower() in plateElementTypes:
                id=int(sline[1])
                elementFaceHandles=[]
                self.element2Face[id]=elementFaceHandles
                numNodes=int(sline[3])
                if numNodes >= 3:
                    #fh=mesh.add_face(vh_list)
                    #self.all_face_handles[fh] = int(sline[1])
                    #elementFaceHandles.append(fh)
                    pass
                    if numNodes == 4:
                        pass
                    else:
                        # unhandled type of the element
                        pass

            elif sline[0].lower() == "dktplate":
                vh_list = [all_vertices.get(sline[4]), all_vertices.get(sline[5]),
                           all_vertices.get(sline[6])]
                all_face_handles[mesh.add_face(vh_list)] = sline[1]
                #mesh.add_face(vh_list)
            elif sline[0].lower() == "mitc4shell":
                vh_list = [all_vertices.get(sline[4]), all_vertices.get(sline[5]),
                           all_vertices.get(sline[6])]
                all_face_handles[mesh.add_face(vh_list)] = sline[1]
                #mesh.add_face(vh_list)
                vh_list = [all_vertices.get( sline[6]), all_vertices.get(sline[7]),
                           all_vertices.get(sline[4])]
                all_face_handles[mesh.add_face(vh_list)] = sline[1]
                #mesh.add_face(vh_list)
            elif sline[0].lower() == "planestress2d":
                vh_list = [all_vertices.get(sline[4]), all_vertices.get(sline[5]),
                           all_vertices.get(sline[6])]
                all_face_handles[mesh.add_face(vh_list)] = sline[1]
                #mesh.add_face(vh_list)
                vh_list = [all_vertices.get(sline[6]), all_vertices.get(sline[7]),
                           all_vertices.get(sline[4])]
                all_face_handles[mesh.add_face(vh_list)] = sline[1]
            elif sline[0].lower() == "trplanestress2d":
                vh_list = [all_vertices.get(sline[4]), all_vertices.get(sline[5]),
                           all_vertices.get(sline[6])]
                all_face_handles[mesh.add_face(vh_list)] = sline[4]
                #mesh.add_face(vh_list)
            elif sline[0].lower() == "trplanestressrotallman":
                vh_list = [all_vertices.get(sline[4]), all_vertices.get(sline[5]),
                           all_vertices.get(sline[6])]
                all_face_handles[mesh.add_face(vh_list)] = sline[4]
                #mesh.add_face(vh_list)
            elif sline[0].lower() == "trplanestrrot":
                vh_list = [all_vertices.get(sline[4]), all_vertices.get(sline[5]),
                           all_vertices.get(sline[6])]
                all_face_handles[mesh.add_face(vh_list)] = sline[4]
                #mesh.add_face(vh_list)
        print(all_face_handles)
        return mesh
        pass


