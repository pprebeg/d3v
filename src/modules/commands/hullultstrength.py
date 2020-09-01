from iohandlers import IOHandler
from signals import Signals
from geometry import Geometry
import openmesh as om
import os
import numpy as np
import readxml
import geofem
import csv





class HullUltStrength ():
    def __init__(self,fileName):
        self.filename = fileName
    def getmesh(self):
        #m=self.test
        m=self.readModel()
        return m

    def test(self):
        mesh= om.TriMesh()
        vhandle = [0]*100
        data = np.array([0, 1, 0])
        #vhandle[0] = mesh.add_vertex(data)
        #vhandle.append(mesh.add_vertex(data))
        data = np.array([1, 0, 0])
        #vhandle[1] = mesh.add_vertex(data)
        data = np.array([2, 1, 0])
        #vhandle[2] = mesh.add_vertex(data)
        data = np.array([0, -1, 0])
        #vhandle[3] = mesh.add_vertex(data)
        data = np.array([2, -1, 0])
        #vhandle[4] = mesh.add_vertex(data)
        #data = np.array([1, 0, 1

        #fh0 = mesh.add_face(vhandle[0], vhandle[1], vhandle[2])
        #fh1 = mesh.add_face(vhandle[1], vhandle[3], vhandle[4])
        #fh2 = mesh.add_face(vhandle[0], vhandle[3], vhandle[1])

        #vh_list = [vhandle[2], vhandle[1], vhandle[4]]
        #fh3 = mesh.add_face(vh_list)

        data = np.array([-12, 0, 0])
        vhandle[5] = mesh.add_vertex(data)
        data = np.array([12, 0, 0])
        vhandle[6] = mesh.add_vertex(data)
        data = np.array([12, 2.8, 0])
        vhandle[7] = mesh.add_vertex(data)
        data = np.array([-12, 2.8, 0])
        vhandle[8] = mesh.add_vertex(data)

        fh4 = mesh.add_face(vhandle[5], vhandle[6], vhandle[7])
        fh5 = mesh.add_face(vhandle[5], vhandle[7], vhandle[8])

        data = np.array([-12, 0, 4.8])
        vhandle[9] = mesh.add_vertex(data)
        data = np.array([12, 0, 4.8])
        vhandle[10] = mesh.add_vertex(data)
        data = np.array([12, 2.8, 4.8])
        vhandle[11] = mesh.add_vertex(data)
        data = np.array([-12, 2.8, 4.8])
        vhandle[12] = mesh.add_vertex(data)

        fh6 = mesh.add_face(vhandle[9], vhandle[10], vhandle[11])
        fh7 = mesh.add_face(vhandle[9], vhandle[11], vhandle[12])


        data = np.array([-12, 0, 9.6])
        vhandle[13] = mesh.add_vertex(data)
        data = np.array([12, 0, 9.6])
        vhandle[14] = mesh.add_vertex(data)
        data = np.array([12, 2.8, 9.6])
        vhandle[15] = mesh.add_vertex(data)
        data = np.array([-12, 2.8, 9.6])
        vhandle[16] = mesh.add_vertex(data)

        fh8 = mesh.add_face(vhandle[13], vhandle[14], vhandle[15])
        fh9 = mesh.add_face(vhandle[13], vhandle[15], vhandle[16])

        data = np.array([-12, 0, 13.6])
        vhandle[17] = mesh.add_vertex(data)
        data = np.array([12, 0, 13.6])
        vhandle[18] = mesh.add_vertex(data)
        data = np.array([12, 2.8, 13.6])
        vhandle[19] = mesh.add_vertex(data)
        data = np.array([-12, 2.8, 13.6])
        vhandle[20] = mesh.add_vertex(data)

        fh10 = mesh.add_face(vhandle[17], vhandle[19], vhandle[18])
        fh11 = mesh.add_face(vhandle[17], vhandle[20], vhandle[19])


        fh12 = mesh.add_face(vhandle[5], vhandle[8], vhandle[20])    #ZAŠTO OVO NERADI???
        fh13 = mesh.add_face(vhandle[5], vhandle[20], vhandle[17])

        fh14 = mesh.add_face(vhandle[6], vhandle[7], vhandle[19])
        fh15 = mesh.add_face(vhandle[6], vhandle[19], vhandle[18])

        # data = np.array([-4, 0, 0])
        # vhandle[21] = mesh.add_vertex(data)
        # data = np.array([-4, 2.8, 0])
        # vhandle[22] = mesh.add_vertex(data)
        # data = np.array([-4, 2.8, 4.8])
        # vhandle[23] = mesh.add_vertex(data)
        # data = np.array([-4, 0, 4.8])
        # vhandle[24] = mesh.add_vertex(data)
        #
        # fh16 = mesh.add_face(vhandle[21], vhandle[22], vhandle[23])
        # fh17 = mesh.add_face(vhandle[21], vhandle[23], vhandle[24])
        #
        # data = np.array([4, 0, 0])
        # vhandle[25] = mesh.add_vertex(data)
        # data = np.array([4, 2.8, 0])
        # vhandle[26] = mesh.add_vertex(data)
        # data = np.array([4, 2.8, 4.8])
        # vhandle[27] = mesh.add_vertex(data)
        # data = np.array([4, 0, 4.8])
        # vhandle[28] = mesh.add_vertex(data)
        #
        # fh18 = mesh.add_face(vhandle[25], vhandle[26], vhandle[27])
        # fh19 = mesh.add_face(vhandle[25], vhandle[27], vhandle[28])
        #
        #
        # data = np.array([-12, 0, 0])
        # vhandle[5] = mesh.add_vertex(data)
        # data = np.array([-12, 2.8, 0])
        # vhandle[8] = mesh.add_vertex(data)
        # data = np.array([-12, 2.8, 13.6])
        # vhandle[20] = mesh.add_vertex(data)
        # data = np.array([-12, 0, 13.6])
        # vhandle[17] = mesh.add_vertex(data)
        #
        # fh16 = mesh.add_face(vhandle[5], vhandle[8], vhandle[20])
        # fh17 = mesh.add_face(vhandle[5], vhandle[20], vhandle[17])
        #
        # data = np.array([12, 0, 0])
        # vhandle[6] = mesh.add_vertex(data)
        # data = np.array([12, 2.8, 0])
        # vhandle[7] = mesh.add_vertex(data)
        # data = np.array([12, 2.8, 13.6])
        # vhandle[19] = mesh.add_vertex(data)
        # data = np.array([12, 0, 13.6])
        # vhandle[18] = mesh.add_vertex(data)
        #
        # fh18 = mesh.add_face(vhandle[6], vhandle[7], vhandle[19])
        # fh19 = mesh.add_face(vhandle[6], vhandle[19], vhandle[18])

        print(mesh.n_faces())
        return mesh
        pass

    def readModel(self):
        with open(self.filename, newline='') as csvfile:
            hfr = csv.reader(csvfile, delimiter='\t', quotechar='|')
            data=[]
            for row in hfr:
                rown=[]
                for x in row:
                    rown.append(x)
                data.append(rown)
        xmlfile = data[0][1]
        abspath1 = '\\'.join(self.filename.split('\\')[0:-1])
        abspath2 = '/'.join(self.filename.split('/')[0:-1])
        if len(abspath2) > len (abspath1):
            abspath = abspath2 + '/' + xmlfile
        else:
            abspath = abspath1 + '\\' + xmlfile

        fem = geofem.GeoFEM()
        m = readxml.MaestroXML(abspath)
        m.readModelToGeoFEM(fem)
        fem.regenerate()
        return fem


