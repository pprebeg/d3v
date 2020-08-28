from iohandlers import IOHandler
from signals import Signals
from geometry import Geometry
import openmesh as om
import os
import numpy as np
import copy
import myfunctions as mf
from hullform import HullForm

class DBBGeometry (Geometry):
	def __init__(self):
		super().__init__()
		self.subdivboxes=[]
		self.nsubdivbox=0
		self.minsubdivbox = 0


class DBBProblem ():
	def __init__(self,fileName):
		self.filename=fileName
		self.hull =0
		self.decks =[]
		self.dbbs=[]


		if (fileName != ""):
			self.readProblem()


	def readProblem(self):
		fnhull = "dsfdsaf" # iz dbb datoteke učitati
		self.hull= DBBHullForm(fnhull)
		for deckIndex in range(len(self.hull.pdecks)):
			self.decks.append(DBBDeck(self.hull,deckIndex))
		pass

	def testProblem(self, scale, block_dims):
		self.hull= DBBHullForm("", scale)
		self.dbbs.append(DBB(self.hull,0))
		#self.dbbs[-1].setPosition(np.array([0,0,0]))
		self.dbbs[-1].testMesh(block_dims)

		#self.dbbs.append(DBB(self.hull, 0))
		#self.dbbs[-1].setPosition(0, -1, -4)
		#self.dbbs[-1].testMesh()

class DBBHullForm (HullForm):
	def __init__(self, fileName, scale):
		super().__init__(fileName)
		self.position = np.array([0.,0.,0.])
		if (fileName != ""):
			self.readHullForm()
		else:
			self.testMesh(scale)

	def readHullForm(self):
		pass

	def regenerateMesh(self):
		self.mesh = mf.make_form(scale = self.scale, move_vector = self.position)

	def move(self, move_vector):
		self.position += move_vector
		self.mesh = mf.move_mesh(self.mesh, move_vector)
		#self.regenerateMesh()

	def setPosition(self, new_position):
		old_position = self.position
		self.position = new_position
		self.mesh = mf.move_mesh(self.mesh, new_position - old_position)
		
	def testMesh(self, scale):
		self.scale = scale
		self.mesh = mf.make_form(scale = self.scale, move_vector = self.position)

		
class DBBDeck (Geometry):
	def __init__(self, hullform, deckIndex):
		super().__init__()
		self.hullform=hullform
		self.z=hullform.pdecks[deckIndex]
		self.deckIndex = deckIndex

	def regenerateMesh(self):
		self.mesh = om.TriMesh()

		pass


class DBB (Geometry):
	def __init__(self, hullform, deck):
		super().__init__()
		self.hullform= hullform
		self.deck=deck
		self.position = np.array([0.,0.,0.])

	def regenerateMesh(self):
		self.mesh= mf.make_block(block_dims = self.block_dims, move_vector = self.position)

	def move(self, move_vector):
		self.position += move_vector
		#self.regenerateMesh()
		self.mesh = mf.move_mesh(self.mesh, move_vector)
		
		#self.position[0] = self.position[0] + dx
		#self.position[1] = self.position[0] + dy
		#self.position[2] = self.position[0] + dz
		#self.regenerateMesh()

	def setPosition(self, new_position):
		old_position = self.position
		self.position = new_position
		self.mesh = mf.move_mesh(self.mesh, new_position - old_position)
		
	def testMesh(self, block_dims):
		self.block_dims = block_dims
		self.position = -block_dims / 2 #centrira block u 0,0,0
		self.position[1] = 0	#centrira block s obzirom na xz ravninu 
		self.mesh= mf.make_block(block_dims = self.block_dims, move_vector = self.position)

	def cutMesh(self):
		form_mesh = self.hullform.mesh
		block_mesh = self.mesh
		self.mesh = mf.cut_meshes(block_mesh, form_mesh)
	
	def calcVolume(self):
		print(mf.calc_mesh_volume(self.mesh))
		
	def IsClosed(self):
		print(mf.is_mesh_closed(self.mesh))
		
		
		
		
		
		
		