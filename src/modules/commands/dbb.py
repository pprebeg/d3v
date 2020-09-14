from iohandlers import IOHandler
from signals import Signals
from geometry import Geometry
import openmesh as om
import os
import numpy as np
import copy
import myfunctions as mf
from hullform import HullForm
import csv

class DBBGeometry (Geometry):
	def __init__(self):
		super().__init__()
		self.subdivboxes=[]
		self.nsubdivbox=0
		self.minsubdivbox = 0


class DBBProblem ():
	def __init__(self,fileName):
		self.filename=fileName
		self.hull = 0
		self.decks = []
		self.dbbs = []


		if (fileName != ""):
			self.readProblem()


	def readProblem(self):
		with open(self.filename, newline='') as csvfile:
			hfr = csv.reader(csvfile, delimiter='\t', quotechar='|')
			data = []
			for row in hfr:
				rown = []
				for x in row:
					rown.append(x)
				data.append(rown)


		abspath1 = '\\'.join(self.filename.split('\\')[0:-1])
		abspath2 = '/'.join(self.filename.split('/')[0:-1])
		if len(abspath2) > len(abspath1):
			abspath = abspath2 + '/'
		else:
			abspath = abspath1 + '\\'
		
		huf_path = abspath + data[0][1]
		block_data_path = abspath + data[1][1]
		#make hull from huf file:
		self.hull= DBBHullForm(huf_path)
		deckz_list= self.hull.pdecks[:-1] #bez keela
		#deckz_list = []
		#deckz_list.append(self.hull.pdecks[2])
		#deckz_list.append(self.hull.pdecks[3])
		#deckz_list.append(self.hull.pdecks[4])
		for deckIndex in range(len(deckz_list)):
			self.decks.append(DBBDeck(self.hull, deckz_list[deckIndex], deckIndex)) 
			
		#make blocks
		with open(block_data_path, "r") as csv_file:
			csv_reader = csv.DictReader(csv_file)
			for row in csv_reader:	#each row contains 1 block data
				deck  = int(row["deck"])
				Ax = float(row["Ax"])
				Ay = float(row["Ay"])	
				Az = self.decks[deck].z
				x = float(row["x"])
				y = float(row["y"])
				
				try:
					z = self.decks[deck - 1].z - self.decks[deck].z
				except IndexError:
					print("Invalid deck position")
				else:
					block_dims = np.array([x,y,z])
					position_A = np.array([Ax,Ay])
				
					block = DBB(self.hull, self.decks[deck], block_dims, position_A, abspath)
					self.dbbs.append(block)
				
		#print(self.hull.wlinesNeg)
		#print(self.hull.wlinesPos[1])
		#print(self.hull)
		#print(len(self.decks))
		#print(self.dbbs)
		#print(self.filename)
	def testProblem(self, scale, block_dims):
		self.hull= DBBHullForm("", scale)
		self.dbbs.append(DBB(self.hull,0))
		#self.dbbs[-1].setPosition(np.array([0,0,0]))
		self.dbbs[-1].genMesh(block_dims)

		#self.dbbs.append(DBB(self.hull, 0))
		#self.dbbs[-1].setPosition(0, -1, -4)
		#self.dbbs[-1].testMesh()

		
		
		
class DBBHullForm (HullForm):
	def __init__(self, fileName, scale = 1):
		super().__init__(fileName)		#vec u inicijalizaciji stvara formu
		self.position = np.array([0.,0.,0.])
		self.centroid = np.array([0.,0.,0.])
		self.LOA = self.shipdata["loa_val"]
		self.centroid = np.array([self.LOA / 2, 0, 0])  # sredina samo za x za sada
		#self.mesh = mf.hard_merge_meshes([self.mesh])

	#def readHullForm(self):
		#HullForm.__init__()


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

		
class DBBDeck(Geometry):
	def __init__(self, hullform, z, deckIndex):
		super().__init__()
		self.hullform = hullform
		self.z = z
		self.deckIndex = deckIndex
		self.genMesh()

	
	def regenerateMesh(self):
		self.mesh = om.TriMesh()

	def genMesh(self):	#za keel koji je na 0 nema wl
		for wline in self.hullform.wlinesPos:
			if np.isclose(self.z, wline[0][2]):
				deck_points = np.asarray(wline)
				break
		self.mesh = mf.make_deck(deck_points, subdivide = False)
		#for i in range(deck_points.shape[0]):
		#	bad_i = []
		#	point = deck_points[i]
		#	next_point = deck_points[(i+1) % deck_points.shape[0]]
		#	if np.allclose(point, next_point):
		#		bad_i.append(i)
		#deck_points = np.delete(deck_points, bad_i, 0)

		
		#new_points = []
		#for point in deck_points:
		#	if point[1] != 0: #ako point nije na osi x 
		#		new_point = copy.copy(point)
		#		new_point[1] = 0
		#		new_points.append(new_point)
		
		#deck_points = np.append(deck_points, np.asarray(new_points), axis = 0)
		#deck_points = np.asarray(deck_points + new_points)
		#deck_points = np.unique(deck_points, axis = 0)		#duplikat na tocki x=50? nakon uniqua, arrejevi su close; ne equal pa ih unique ne reze?
		

		
		#print(deck_points)
		
		#print(deck_points[18],deck_points[19],deck_points[20],)
		#print(np.allclose(deck_points[18],deck_points[19]))
		#self.mesh = mf.cut_meshes(self.mesh, self.hullform.mesh) 	#predugo traje

	def move(self, move_vector):
		self.z += move_vector[2]
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

		
		
		
		
		
		
		
		
		
class DBB(Geometry):
	def __init__(self, hullform, deck:DBBDeck,block_dims, position, abspath):
		super().__init__()
		self.folderpath = abspath
		self.hullform= hullform
		self.deck=deck
		self.position = np.array([position[0],position[1],deck.z])
		self.block_dims=block_dims
		self.genMesh()
		self.cutMesh()
		#print(self.hullform.filename)
		
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

		self.position = new_position
		self.genMesh()

		old_position = self.position
		self.position = new_position
		self.mesh = mf.move_mesh(self.mesh, new_position - old_position)
		
	def genMesh(self):
		self.mesh = mf.make_block_from_unit_csv(self.block_dims, self.position, self.folderpath)
		#self.mesh= mf.make_block(block_dims = self.block_dims, move_vector = self.position)

	def cutMesh(self):
		mf.fit_block_to_form(self.mesh, self.block_dims, self.position, self.hullform.mesh)
		#mf.cut_meshes(self.mesh, self.hullform.mesh)
		pass
		
	def calcVolume(self):
		print(mf.calc_mesh_volume(self.mesh))
		
	def IsClosed(self):
		print(mf.is_mesh_closed(self.mesh))
		
		
		
		
		
		
		