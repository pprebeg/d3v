from iohandlers import IOHandler
from signals import Signals
from geometry import Geometry
import openmesh as om
import os
import numpy as np
import copy
import myfunctions as mf

class DBBProblem ():
	def __init__(self,fileName):
		self.filename=fileName
		self.hull =0
		self.decks =[]
		self.dbbs=[]
		if (fileName != ""):
			self.readProblem()
		else:
			self.testProblem()

	def readProblem(self):
		pass
	def testProblem(self):
		self.hull= DBBHullForm("")
		self.dbbs.append(DBB(self.hull,0))
		self.dbbs[-1].setPosition(0,0,0)
		self.dbbs[-1].testMesh()

		#self.dbbs.append(DBB(self.hull, 0))
		#self.dbbs[-1].setPosition(0, -1, -4)
		#self.dbbs[-1].testMesh()

class DBBHullForm (Geometry):
	def __init__(self, fileName):
		self.filename = fileName
		super().__init__()
		if (fileName != ""):
			self.readHullForm()
		else:
			self.testMesh()

	def readHullForm(self):
		pass

	def regenerateMesh(self):
		pass

	def testMesh(self):
		self.mesh = mf.make_form()

class DBBDeck (Geometry):
	def __init__(self, hullform):
		super().__init__()
		self.hullform=hullform

	def regenerateMesh(self):
		self.testMesh()

	def testMesh(self):
		pass

class DBB (Geometry):
	def __init__(self, hullform, deck):
		super().__init__()
		self.hullform= hullform
		self.deck=deck
		self.position = np.array([0,0,0])

	def regenerateMesh(self):
		self.mesh= mf.make_block(move_vector = self.position)

	def move(self,dx,dy,dz):
		self.position[0] = self.position[0] + dx
		self.position[1] = self.position[0] + dy
		self.position[2] = self.position[0] + dz
		self.regenerateMesh()

	def setPosition(self,x,y,z):
		self.position[0] = x
		self.position[1] = y
		self.position[2] = z

	def testMesh(self):
		self.mesh= mf.make_block(move_vector = self.position)

		
		
		
		
		
		
		



