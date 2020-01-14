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
        g = Geometry()
        g.mesh = oofem.getmesh()
        Signals.get().geometryImported.emit(g)

    def getImportFormats(self):
        return []


def createIOHandler():
    return OOFEMImporter()

class OOFEM ():
    def __init__(self,fileName):
        self.filename=fileName

    def getmesh(self):
        #mesh=self.test()
        #mesh = self.test2()
       #mesh = self.oofemmesh()

        #mesh = self.pokusaj1()
        mesh = self.pokusaj2()


        return mesh

    def pokusaj2(self):
        prvaTocka=[0,0,0]
        drugaTocka=[0,0,0]
        orient=[0,1,0] # orijentacija struka
        hw=0.3
        tw=0.03
        bf=0.2
        tf=0.02
        bp=0.25
        tp=0.025
        self.zadatak2(prvaTocka, drugaTocka, orient, hw,tw,bf,tf,bp,tp)

    def zadatak2(prvaTocka, drugaTocka, orient, hw,tw,bf,tf,bp,tp):
        mesh= om.TriMesh()
        #implementacija
        return mesh


    def pokusaj1(self):
        prvaTocka = []
        for i in range(0, 3):
            ele = float(input("Unesite koordinatu " + str(i + 1) + " prve tocke: "))
            prvaTocka.append(ele)

        drugaTocka = []
        for i in range(0, 3):
            ele = float(input("Unesite koordinatu " + str(i + 1) + " druge tocke: "))
            drugaTocka.append(ele)

        visinaS = []
        for i in range(0, 3):
            ele = float(input("Unesite koordinatu " + str(i + 1) + " visine struka: "))
            visinaS.append(ele)

        polaSirineF = []
        for i in range(0, 3):
            ele = float(input("Unesite koordinatu " + str(i + 1) + " sirine flanze: "))
            polaSirineF.append(ele)

        debljinaF = []
        for i in range(0, 3):
            ele = float(input("Unesite koordinatu " + str(i + 1) + " debljine flanze: "))
            debljinaF.append(ele)

        debljinaS = []
        for i in range(0, 3):
            ele = float(input("Unesite koordinatu " + str(i + 1) + " debljine struka: "))
            debljinaS.append(ele)

        sirinaP = []
        for i in range(0, 3):
            ele = float(input("Unesite koordinatu " + str(i + 1) + " sirine pokrova: "))
            sirinaP.append(ele)

        debljinaP = []
        for i in range(0, 3):
            ele = float(input("Unesite koordinatu " + str(i + 1) + " debljine pokrova: "))
            debljinaP.append(ele)

        print(prvaTocka)
        print(drugaTocka)
        print(visinaS)
        print(polaSirineF)
        print(debljinaF)
        print(debljinaS)
        print(sirinaP)
        print(debljinaP)
        mesh=self.zadatak(prvaTocka, drugaTocka, [visinaS,polaSirineF,debljinaF,debljinaS,sirinaP,debljinaP])
        #mesh=self.zadatak([0, 0, 0], [1, 0, 0], [[0,0.12,0],[0,0,0.03],[0.004,0,0],[0,0.004,0],[0,0,0.04],[0.004,0,0]])
        return mesh
    def zadatak(self, tocka1, tocka2, skalari, mesh=om.TriMesh()):
        print(tocka1)
        numpyTocka1 = np.array(tocka1)
        numpyTocka2 = np.array(tocka2)
        numpyTocka1_ = numpyTocka1 + np.array([0, 0, 1])
        numpyTocka2_ = numpyTocka2 + np.array([0, 0, 1])
        numpyTocka1a = numpyTocka1 + np.array([0, 0, 2])
        numpyTocka2a = numpyTocka2 + np.array([0, 0, 2])
        if(len(skalari) == 2):
            mesh.add_face(mesh.add_vertex(numpyTocka1), mesh.add_vertex(numpyTocka2), mesh.add_vertex(numpyTocka2+skalari[0]))
            mesh.add_face(mesh.add_vertex(numpyTocka2+skalari[0]), mesh.add_vertex(numpyTocka1 + skalari[0]), mesh.add_vertex(numpyTocka1))
        elif(len(skalari) == 4):
            mesh.add_face(mesh.add_vertex(numpyTocka1), mesh.add_vertex(numpyTocka2),mesh.add_vertex(numpyTocka2 + skalari[0]))
            mesh.add_face(mesh.add_vertex(numpyTocka2 + skalari[0]), mesh.add_vertex(numpyTocka1 + skalari[0]),mesh.add_vertex(numpyTocka1))

            mesh.add_face(mesh.add_vertex(numpyTocka1_), mesh.add_vertex(numpyTocka2_),mesh.add_vertex(numpyTocka2_ + skalari[0]))
            mesh.add_face(mesh.add_vertex(numpyTocka2_ + skalari[0]), mesh.add_vertex(numpyTocka1_ + skalari[0]),mesh.add_vertex(numpyTocka1_))
            mesh.add_face(mesh.add_vertex(numpyTocka1_ + skalari[0] - skalari[1]), mesh.add_vertex(numpyTocka1_ + skalari[0] + skalari[1]), mesh.add_vertex(numpyTocka2_ + skalari[0] + skalari[1]))
            mesh.add_face(mesh.add_vertex(numpyTocka2_ + skalari[0] + skalari[1]), mesh.add_vertex(numpyTocka2_ + skalari[0] - skalari[1]), mesh.add_vertex(numpyTocka1_ + skalari[0] - skalari[1]))
        elif(len(skalari) == 6):
            mesh.add_face(mesh.add_vertex(numpyTocka1), mesh.add_vertex(numpyTocka2),mesh.add_vertex(numpyTocka2 + skalari[0]))
            mesh.add_face(mesh.add_vertex(numpyTocka2 + skalari[0]), mesh.add_vertex(numpyTocka1 + skalari[0]),mesh.add_vertex(numpyTocka1))

            mesh.add_face(mesh.add_vertex(numpyTocka1_), mesh.add_vertex(numpyTocka2_),mesh.add_vertex(numpyTocka2_ + skalari[0]))
            mesh.add_face(mesh.add_vertex(numpyTocka2_ + skalari[0]), mesh.add_vertex(numpyTocka1_ + skalari[0]),mesh.add_vertex(numpyTocka1_))
            mesh.add_face(mesh.add_vertex(numpyTocka1_ + skalari[0] - skalari[1]),mesh.add_vertex(numpyTocka1_ + skalari[0] + skalari[1]),mesh.add_vertex(numpyTocka2_ + skalari[0] + skalari[1]))
            mesh.add_face(mesh.add_vertex(numpyTocka2_ + skalari[0] + skalari[1]),mesh.add_vertex(numpyTocka2_ + skalari[0] - skalari[1]),mesh.add_vertex(numpyTocka1_ + skalari[0] - skalari[1]))

            mesh.add_face(mesh.add_vertex(numpyTocka1a), mesh.add_vertex(numpyTocka2a),mesh.add_vertex(numpyTocka2a + skalari[0]))
            mesh.add_face(mesh.add_vertex(numpyTocka2a + skalari[0]), mesh.add_vertex(numpyTocka1a + skalari[0]),mesh.add_vertex(numpyTocka1a))
            mesh.add_face(mesh.add_vertex(numpyTocka1a + skalari[0] - skalari[1]),mesh.add_vertex(numpyTocka1a + skalari[0] + skalari[1]),mesh.add_vertex(numpyTocka2a + skalari[0] + skalari[1]))
            mesh.add_face(mesh.add_vertex(numpyTocka2a + skalari[0] + skalari[1]),mesh.add_vertex(numpyTocka2a + skalari[0] - skalari[1]), mesh.add_vertex(numpyTocka1a + skalari[0]- skalari[1]))
            mesh.add_face(mesh.add_vertex(numpyTocka1a - skalari[4]), mesh.add_vertex(numpyTocka1a + skalari[4]), mesh.add_vertex(numpyTocka2a + skalari[4]))
            mesh.add_face(mesh.add_vertex(numpyTocka2a + skalari[4]), mesh.add_vertex(numpyTocka2a - skalari[4]), mesh.add_vertex(numpyTocka1a - skalari[4]))

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
        f = open(self.filename, newline='')
        all_vertices = {}
        for line in f:
            line = ' '.join(line.split())
            sline = line.split(" ")
            if len(sline) < 3:
                continue
            #print(sline[0].lower())
            if sline[0] == "node":
                d = [float(sline[4]), float(sline[5]), float(sline[6])]
                all_vertices[sline[1]] = mesh.add_vertex(d)
            elif sline[0].lower() == "dktplane":
                vh_list = [all_vertices.get(sline[4]), all_vertices.get(sline[5]),
                           all_vertices.get(sline[6])]
                fh=mesh.add_face(vh_list)
                #if fh==None:
                    #print('nije dodan dface')
            elif sline[0].lower() == "mitc4shell":
                vh_list = [all_vertices.get(sline[4]), all_vertices.get(sline[5]),
                           all_vertices.get(sline[6])]
                mesh.add_face(vh_list)
                vh_list = [all_vertices.get( sline[6]), all_vertices.get(sline[7]),
                           all_vertices.get(sline[4])]
                mesh.add_face(vh_list)
            elif sline[0].lower() == "planestress2d":
                vh_list = [all_vertices.get(sline[4]), all_vertices.get(sline[5]),
                           all_vertices.get(sline[6]), all_vertices.get(sline[7])]
                mesh.add_face(vh_list)
            elif sline[0].lower() == "trplanestress2d":
                vh_list = [all_vertices.get(sline[4]), all_vertices.get(sline[5]),
                           all_vertices.get(sline[6])]
                mesh.add_face(vh_list)
            elif sline[0].lower() == "trplanestressrotallman":
                vh_list = [all_vertices.get(sline[4]), all_vertices.get(sline[5]),
                           all_vertices.get(sline[6])]
                mesh.add_face(vh_list)
            elif sline[0].lower() == "trplanestrrot":
                vh_list = [all_vertices.get(sline[4]), all_vertices.get(sline[5]),
                           all_vertices.get(sline[6])]
                mesh.add_face(vh_list)
        return mesh
        pass


