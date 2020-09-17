from iohandlers import IOHandler
from signals import Signals
import os
import numpy as np
import readxml
import geofem
import csv





class HullUltStrength (geofem.GeoFEM):
    def __init__(self,fileName):
        self.filename = fileName
        super().__init__()
        self.readModel()

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

        m = readxml.MaestroXML(abspath)
        m.readModelToGeoFEM(self)
        self.regenerate()
        return


