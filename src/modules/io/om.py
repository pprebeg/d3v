from iohandlers import IOHandler
from signals import Signals
from geometry import Geometry
import openmesh as om
import os

class OpenMeshImporter(IOHandler):
    def __init__(self):
        super().__init__()
        Signals.get().importGeometry.connect(self.importGeometry)

    def importGeometry(self, fileName):
        g = Geometry()
        fnam, fext = os.path.splitext(fileName)
        try:
            if fext == ".ply":
                m = om.read_trimesh(fileName,vertex_normal=True)
            else:
                m = om.read_trimesh(fileName)
            pass
        except:
            print("File not supported for read with openmesh")
            return
        g.mesh = m
        Signals.get().geometryImported.emit(g)
    def getImportFormats(self):
        return []
def createIOHandler():
    return OpenMeshImporter()