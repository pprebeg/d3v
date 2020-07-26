import openmesh as om
import numpy as np
import copy

def print_data(mesh, report):

	for faceh in mesh.faces():
		print("Face index: ", faceh.idx(), "\n")
		report.write("Face index: " + str(faceh.idx()) + "\n")
		
		vertexi = mesh.face_vertex_indices()[faceh.idx()]
		face_points = mesh.points()[vertexi]
		print("Face points: \n", face_points)
		report.write("Face points: \n" + str(face_points) + "\n\n")
		
		
		for edgeh in mesh.fe(faceh): 
			print("\nEdge index: \n", edgeh.idx())
			report.write("Edge index: " + str(edgeh.idx()) + "\n\n")
			
			
			edgei = mesh.edge_vertex_indices()[edgeh.idx()]
			edge_points = mesh.points()[edgei]
			edge_vector = mesh.calc_edge_vector(edgeh)
			
			print("Edge points: \n", edge_points)
			print("Edge vector: \n", edge_vector)
			
			report.write("Edge points: \n" + str(edge_points) + "\n\n")
			report.write("Edge vector: \n" + str(edge_vector) + "\n\n")
			

def test_split(mesh):
	for fh in mesh.faces():
		center = mesh.calc_face_centroid(fh)
		vh = mesh.add_vertex(center)
		#mesh.split(fh,vh)
		mesh.triangulate(fh)
			
			
			
#---------------------------------------			
def sort_by_dist(points):
	aritm_middle = np.sum(points, axis = 0)/points.shape[0]
	start_point = aritm_middle
	sorted_array_shape = list(points.shape)
	sorted_array_shape[0] = 0
	sorted_array = np.empty(sorted_array_shape)

	
	for p in range(points.shape[0]):
		print(start_point)
		print("-------------")
		print(points)
		d = start_point - points
		distance = np.sum(d**2, axis = 1)**0.5
		i = np.array(np.where(distance == distance.max())).flatten()[0]	#sto ako su 2 tocke jednako udaljene?	[0] da odabere samo 1
		print("index" + str(i))
		start_point = points[i]
		points = np.delete(points, i, axis = 0)
		print("star point"+str(start_point))
		sorted_array = np.append(sorted_array, np.expand_dims(start_point, axis = 0), axis = 0)
		
	return sorted_array
	


		
		
def merge_meshes(*args):			#problem ako imamo 2 ista pointa u oba mesha: radi 2 verteksa
	points = np.empty((0,3))
	face_vertex_indices = np.empty((0,3))
	
	for mesh in args:
		if len(mesh.face_vertex_indices().shape) == 1:		#za slucaj da nema fvi u mesh
			add_fvi = np.empty((0,3))
		else:
			add_fvi = mesh.face_vertex_indices()
		face_vertex_indices = np.append(face_vertex_indices, add_fvi + points.shape[0], axis = 0)				#+points.chape[0] je tu da poreda face_vertex_indices sa njihovim indexom u novom arrayu
		points = np.append(points, mesh.points(), axis = 0)
	
	combined_mesh = om.TriMesh(points, face_vertex_indices)	
	
	return combined_mesh

	
	
def rotate_vector(vector, theta):
	theta = theta * np.pi/180		#pretvara u radijane
	s = np.sin(theta)
	c = np.cos(theta)
	rot_matrix = np.array([[c, -s],[s, c]])
	return np.matmul(rot_matrix, vector.T).T
			
def triangle_surface(triangle_points):
	AB = triangle_points[1]-triangle_points[0]
	AC = triangle_points[2]-triangle_points[0]
	surface = np.linalg.norm(np.cross(AB, AC))/2
	return surface


	
def is_inside_triangle(point, triangle_points, include_edges = True):
	ABC = triangle_surface(triangle_points)
	Ai = np.empty((0))
	for i in range(triangle_points.shape[0]):
		points = copy.copy(triangle_points)
		points[i] = point
		Ai = np.append(Ai, triangle_surface(points))
	
	if np.isclose(np.sum(Ai), ABC) == True:
		if include_edges == True:
			return True
		elif include_edges == False:
			if np.isclose(np.any(Ai), 0) == True:
				return False
			else:
				return True
	
	else:
		return False

def is_inside_sphere(point, c, R):		#point je tocka koju provjeravamo, c je radij vector centra sfere, R je radius sfere , r je radius od centra do pointa
	r = np.sum((point - c) ** 2) ** 0.5
	if r <= R:
		return True
	else: 
		return False
		
def calc_triangle_circumradius(triangle_points):
	a = triangle_points[0] - triangle_points[1]	#vector
	a = np.sum(a**2)**0.5
	b = triangle_points[0] - triangle_points[2]
	b = np.sum(b**2)**0.5
	c = triangle_points[1] - triangle_points[2]
	c = np.sum(c**2)**0.5
	
	R = a*b*c/(((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c))**0.5)
	return R
	
def calc_triangle_circumcenter(triangle_points):
	a = triangle_points[0]
	b = triangle_points[1]
	c = triangle_points[2]
	ca = c-a
	ba = b-a
	baXca = np.cross(ba, ca)
	
	v1 = np.sum(ca**2) * np.cross(baXca, ba)
	v2 = np.sum(ba**2) * np.cross(-baXca, ca)
	v3 = 2 * np.sum(baXca**2)
	
	m = a + (v1+v2)/v3
	return m
	
def make_block(block_dims = np.array([20,6,3]), move_vector = np.array([0,0,0])):
	mesh = om.TriMesh()
	axes = []
	#stvara 2 tocke na svakoj osi
	for dim in block_dims:
		axes.append(np.linspace(0, dim, 2))
		
	block_corners = np.asarray(np.meshgrid(*axes)).T.reshape(8,3)
	block_corners += move_vector
	
	corner_vertices = []
	for corner in block_corners:
		corner_vertices.append(mesh.add_vertex(corner))

		
		
	#x+face
	mesh.add_face(corner_vertices[2],corner_vertices[3],corner_vertices[6])
	mesh.add_face(corner_vertices[3],corner_vertices[7],corner_vertices[6])
		
	#x-face
	mesh.add_face(corner_vertices[0],corner_vertices[4],corner_vertices[1])
	mesh.add_face(corner_vertices[1],corner_vertices[4],corner_vertices[5])
		
	#y+face
	mesh.add_face(corner_vertices[3],corner_vertices[1],corner_vertices[5])
	mesh.add_face(corner_vertices[3],corner_vertices[5],corner_vertices[7])
		
	#y-face
	mesh.add_face(corner_vertices[0],corner_vertices[2],corner_vertices[4])
	mesh.add_face(corner_vertices[2],corner_vertices[6],corner_vertices[4])
		
	#z+face
	mesh.add_face(corner_vertices[4],corner_vertices[6],corner_vertices[5])
	mesh.add_face(corner_vertices[6],corner_vertices[7],corner_vertices[5])
		
	#z-face
	mesh.add_face(corner_vertices[2],corner_vertices[0],corner_vertices[1]) 
	mesh.add_face(corner_vertices[2],corner_vertices[1],corner_vertices[3])
		
	#print("block")
	#print ("points",block_corners)
	#print ("indices",mesh.face_vertex_indices())
		
		
	return mesh

def make_form():
	mesh = om.TriMesh()
	points = np.array([[1,2,-1],[2,1,-1]])*5
	vhandles = []
	fhandles = []
		
	#x mirrored points
	mirrored_points = copy.deepcopy(points)
	mirrored_points[:,1:2] *= -1
		
	points = np.append(points, mirrored_points, axis = 0)
	
	#z mirrored points
	mirrored_points = copy.deepcopy(points)
	mirrored_points[:,2:] *= -1
		
	points = np.append(points, mirrored_points, axis = 0)
			
	for point in points:
		vhandles.append(mesh.add_vertex(point))
		

	fhandles.append(mesh.add_face(vhandles [0], vhandles [1], vhandles [4]))
	fhandles.append(mesh.add_face(vhandles [4], vhandles [1], vhandles [5]))
	fhandles.append(mesh.add_face(vhandles [5], vhandles [1], vhandles [7]))
	fhandles.append(mesh.add_face(vhandles [7], vhandles [1], vhandles [3]))
	fhandles.append(mesh.add_face(vhandles [7], vhandles [3], vhandles [2]))
	fhandles.append(mesh.add_face(vhandles [7], vhandles [2], vhandles [6]))
	
	return mesh
	
def make_supertetrahedron(points):
	supertetrahedron = om.TriMesh()
	supertetrahedron_points = np.array([[0,0,1],[(8/9)**0.5,0,-1/3],[-(2/9)**0.5,(2/3)**0.5,-1/3],[-(2/9)**0.5,-(2/3)**0.5,-1/3]])		#tocke za R = 1
	vhandles = []
	fhandles = []
	
	arithmetic_middle = np.sum(points, axis = 0) / points.shape[0]
	r = np.max(np.sum((points - arithmetic_middle)**2, axis = 1) ** 0.5)  #radius je najveci delta r or svih pointova  
	R = 3.05*r			#normalni radius opisane sfere je 3r ali ja stavljam 3.05 zbog sigurnosti () za slucaj da je neka tocka bas na sjecistu upisane sfere i tetrahedrona 
	supertetrahedron_points = arithmetic_middle + supertetrahedron_points * R		#scaleanje osnovnih tocaka za pravi R

	for point in supertetrahedron_points:
		vhandles.append(supertetrahedron.add_vertex(point))
	

	fhandles.append(supertetrahedron.add_face(vhandles [0], vhandles [1], vhandles [2]))
	fhandles.append(supertetrahedron.add_face(vhandles [0], vhandles [2], vhandles [3]))
	fhandles.append(supertetrahedron.add_face(vhandles [0], vhandles [3], vhandles [1]))
	fhandles.append(supertetrahedron.add_face(vhandles [3], vhandles [2], vhandles [1]))
	
	return supertetrahedron
	
def make_supertriangle(points):	
	supertriangle= om.TriMesh()
	supertriangle_points = np.array([[1,0,0],[-1/2,3**0.5/2,0],[-1/2,-3**0.5/2,0]])
	fhandles = []
	vhandles = []
	
	arithmetic_middle = np.sum(points, axis = 0) / points.shape[0]
	r = np.max(np.sum((points - arithmetic_middle)**2, axis = 1) ** 0.5)  #radius je najveci delta r or svih pointova  
	R = 2*r*1.05	#sigurnost od 1.05
	supertriangle_points = arithmetic_middle + supertriangle_points * R
	
	for point in supertriangle_points:
		vhandles.append(supertriangle.add_vertex(point))
	fhandles.append(supertriangle.add_face(vhandles[0],vhandles[1],vhandles[2]))
	
	return supertriangle
	
	
def get_intersection(radius, vector, face_points):
	face_radius = face_points[0]
	face_vector1 = face_points[1] - face_points[0]
	face_vector2 = face_points[2] - face_points[0]
	
	face_radius, face_vector1, face_vector2,
	
	radius_matrix = (radius - face_radius).T
	vector_matrix = np.empty((3,3))
	vector_matrix[:,0] = face_vector1
	vector_matrix[:,1] = face_vector2
	vector_matrix[:,2] = -vector
		
	try:
		edge_parameter = np.linalg.solve(vector_matrix, radius_matrix)[2]
	except:
		return None
	
	intersection_point = radius + (vector * edge_parameter)
	if is_inside_triangle(intersection_point, face_points) == True:
		return intersection_point
	else:
		return None

def get_meshes_intersections(mesh1, mesh2):
	intersection_points = np.empty((0,3))
	for eh in mesh1.edges():
		edge_vector = mesh1.calc_edge_vector(eh)
		evi = mesh1.edge_vertex_indices()[eh.idx()]		#evi je edge vertex indices (array sa idovima pointova za taj edge)
		edge_radius = mesh1.points()[evi][0]			#za radius nam treba samo jedna tocka 
		
		for fh in mesh2.faces():
			fvi = mesh2.face_vertex_indices()[fh.idx()]
			face_points = mesh2.points()[fvi]
			
			intersection_point = get_intersection(edge_radius, edge_vector, face_points)
			if intersection_point is not None:		#uvijet ako intersection point nije u trokutu ili edge i face su paralelni:
				intersection_points = np.append(intersection_points, intersection_point.reshape(1,3), axis = 0)
			
	for eh in mesh2.edges():
		edge_vector = mesh2.calc_edge_vector(eh)
		evi = mesh2.edge_vertex_indices()[eh.idx()]		
		edge_radius = mesh2.points()[evi][0]			
		
		for fh in mesh1.faces():
			fvi = mesh1.face_vertex_indices()[fh.idx()]
			face_points = mesh1.points()[fvi]
			
			intersection_point = get_intersection(edge_radius, edge_vector, face_points)
			if intersection_point is not None:		
				intersection_points = np.append(intersection_points, intersection_point.reshape(1,3), axis = 0)
			
	return np.unique(intersection_points, axis = 0) 					#mice sve duplikate tocaka


def BowyerWatson_triangulation_algorithm(points, no_boundary = True):
	points = sort_by_dist(points)
	mesh = make_supertriangle(points)
	start_points = copy.copy(mesh.points())
	
	
	for point in points:	#stavlj vertexe u mesh i identificira badTriangles
		vh = mesh.add_vertex(point)		
		bad_fvi = np.empty((0,3))
		bad_fh = []
		for fh in mesh.faces():
			fvi = mesh.face_vertex_indices()[fh.idx()]
			triangle_points = mesh.points()[fvi]
			c = calc_triangle_circumcenter(triangle_points)
			R = calc_triangle_circumradius(triangle_points)
			if is_inside_sphere(point, c, R) == True:
				bad_fvi = np.append(bad_fvi, fvi.reshape(1,3), axis = 0)
				bad_fh.append(fh)
		
		#badTriangles mesh da mogu lagano traziti boundary edgeve
		badTriangles = om.TriMesh(mesh.points(), bad_fvi)		#inhereta pointove iz mesha!
		polygon = []		#evi boundary edgeva koje kasnije retruanguliramo
		
		#find the boundary of the polygonal hole
		for fh in badTriangles.faces():					
			for eh in badTriangles.fe(fh):
				if badTriangles.is_boundary(eh) == True:
					evi = badTriangles.edge_vertex_indices()[eh.idx()]
					polygon.append(evi) 
		
		for fh in bad_fh:
			mesh.delete_face(fh, False)
		mesh.garbage_collection()
		
		#retriangualte
		for evi in polygon:								#tu dodaj uvijet za complex edge
			vh0 = mesh.vertex_handle(evi[0])			#jeli vertex indice isti kao vh.idx()
			vh1 = mesh.vertex_handle(evi[1])
			vh2 = vh
			n_faces_before = len(mesh.faces())
			mesh.add_face(vh0, vh1, vh2) 			#tako da zadrzi orijentaciju
			n_faces_after = len(mesh.faces())
							
			if n_faces_after == n_faces_before:			#ako je broj faceva isti face se nije dodao pa mijenjamo orijentaciuju		
				mesh.add_face(vh2, vh1, vh0)		#reverse orientation
				

	
	#delete start vertices
	for vh in mesh.vertices():
		vpoint = mesh.points()[vh.idx()]
		for start_point in start_points:
			if np.array_equal(vpoint, start_point) == True:
				mesh.delete_vertex(vh)
	mesh.garbage_collection()	


	return mesh


   
   
def BowyerWatson_triangulation_algorithm_backup(points, report, no_boundary = True):	#jos ga muci nes; ako su tocke kao u primjeru badT su svi facevi i brise ih ali ne ostavlja edgeve sa prosim pointom iteracije (svi facevi su u bad Tri)
	#da probam max r po max r od pointa do pointa?
	points = sort_by_dist(points)
	report.write("input points: \n" + str(points) + "\n\n")
	mesh = make_supertriangle(points)
	start_points = copy.copy(mesh.points())
	
	report.write("start_vert indices: \n" + str(start_points) + "\n\n")
	report.write("all points: \n" + str(np.append(mesh.points(), points, axis = 0)) + "\n\n")
	
	
	for point in points:
		print(point) 
		report.write("iteration for point:  \n" + str(point) + "\n\n")
		vh = mesh.add_vertex(point)		
		bad_fvi = np.empty((0,3))
		bad_fh = []
		for fh in mesh.faces():
			fvi = mesh.face_vertex_indices()[fh.idx()]
			triangle_points = mesh.points()[fvi]
			c = calc_triangle_circumcenter(triangle_points)
			R = calc_triangle_circumradius(triangle_points)
			if is_inside_sphere(point, c, R) == True:
				bad_fvi = np.append(bad_fvi, fvi.reshape(1,3), axis = 0)
				bad_fh.append(fh)
				
		report.write("bad fvi: \n" + str(bad_fvi) + "\n\n")
	
		
		badTriangles = om.TriMesh(mesh.points(), bad_fvi)		#inhereta pointove iz mesha!
		polygon = []		#evi
		
		for fh in badTriangles.faces():					#find the boundary of the polygonal hole
			for eh in badTriangles.fe(fh):
				print(badTriangles.edge_vertex_indices()[eh.idx()])
				print(badTriangles.is_boundary(eh))
				if badTriangles.is_boundary(eh) == True:
					evi = badTriangles.edge_vertex_indices()[eh.idx()]#neznam dali je u meshu i badtrianglesu eh isti za svaki edge najvj ne jer je generiran iz tocaka na pocetku, ali ako je generiran iz istih tocaka mozda su vertici isti(por. tocaka)
					print(evi)
					polygon.append(evi) 
					
		
		
		report.write("polygon: \n" + str(polygon) + "\n\n")
		report.write("mesh points: \n" + str(mesh.points()) + "\n\n")
		report.write("mesh face vertex indices: \n" + str(mesh.face_vertex_indices()) + "\n\n")
		
		
		for fh in bad_fh:
			mesh.delete_face(fh, False)
		mesh.garbage_collection()
		
		
		
		#for fh in badTriangles.faces():				#nisu isti fh!!
		#	mesh.delete_face(fh, False)			#false da ne deleta izolirane vertexe
		#mesh.garbage_collection()
		
		report.write("mesh face vertex indices after deleting bad triangles: \n" + str(mesh.face_vertex_indices()) + "\n\n")
		#retriangualte
		for vh in mesh.vertices():
			print("aaaaa"+str(mesh.points()[vh.idx()])+"aaaaaa")
		for evi in polygon:								#tu dodaj uvijet za complex edge
			vh0 = mesh.vertex_handle(evi[0])			#jeli vertex indice isti kao vh.idx()
			vh1 = mesh.vertex_handle(evi[1])
			vh2 = vh
			n_faces_before = len(mesh.faces())
			mesh.add_face(vh0, vh1, vh2) #tako da zadrzi orijentaciju
			n_faces_after = len(mesh.faces())
			
			report.write("edge points to be merged with point \n" + str(mesh.points()[evi]) + "\n\n")
			
			
			
			if n_faces_after == n_faces_before:			#ako je broj faceva isti face se nije dodao pa mijenjamo orijentaciuju		
				report.write("same N of faces; executing first reverse orientation: \n" + str(n_faces_after) + "\n\n")
				mesh.add_face(vh2, vh1, vh0)		#reverse orientation
				n_faces_after = len(mesh.faces())
				report.write("N faces after first reverse orientation: \n" + str(n_faces_after) + "\n\n")
			else:
				report.write("different N of faces: \n" + str(n_faces_after) + "\n\n")
			#
			#n_faces_after = len(mesh.faces())
			#if n_faces_after == n_faces_before:			
			#	report.write("STILL the same n of faces after reverse orientation: \n" + str(n_faces_after) + "\n\n")

			
		report.write("mesh face vertex indices after adding faces from polygon: \n" + str(mesh.face_vertex_indices()) + "\n\n")
	
	report.write("all points before deleting all vertices: \n" + str(mesh.points()) + "\n\n")		
	
	#delete start vertices	: brise i stvari koje nisu start pointovi
	#vh iz starta nije isti na kraju
	#algoritam ima problema: izbacuje izolirane vertekse(ne popunjava ih sve) najvj jer u nekeom koraku se brisu sva tri facea na koji je vertex spojen; provjeri jos u reportu!!!!!!!!!!!!!!!
	for vh in mesh.vertices():
		vpoint = mesh.points()[vh.idx()]
		for start_point in start_points:
			if np.array_equal(vpoint, start_point) == True:
				report.write("start point koji je algoritam uzeo za brisanje: \n" + str(vpoint) + "\n\n")
				mesh.delete_vertex(vh, False)
	mesh.garbage_collection()	
		
	#for id in start_vi:			#brise sve vertekse zbog nekog razloga takodjer baca complex edgeve!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! na drugoj tocki daje 2 complex edga testiraj reference za svaki slucaj 	
	#	vh = mesh.vertex_handle(id)
	#	mesh.delete_vertex(vh)			#vh iz start vh se promijenio; vise nisu tocke supertrokuta
	#mesh.garbage_collection()	
	report.write("mesh face vertex indices after deleting start vertices: \n" + str(mesh.face_vertex_indices()) + "\n\n")
	report.write("all points after deleting start vertices: \n" + str(mesh.points()) + "\n\n")
	#isloated_points = np.empty((0,3))
	#for vh in mesh.vertices():
	#	i = 0
	#	for fh in mesh.vf(vh):
	#		print(fh.idx())
	#		vi = mesh.vertex_vertex_indices()[vh.idx()]
	#		isloated_points = np.append(isloated_points, mesh.points()[vi], axis = 0)
		
	#report.write("isolated points: \n" + str(points) + "\n\n")
	
	return mesh

   
   
			
def cut_meshes2(block, form , report):		
	intersection_points = get_meshes_intersections(block, form)
	mesh = BowyerWatson_triangulation_algorithm_backup(intersection_points, report)
	print(intersection_points)
	return mesh
	

	
	
	
#tests:
#edge_radius = np.array([0,0,0])
#edge_vector = np.array([5,0,0])
#face_points = np.array([[2.5,-2.5,0],[2.5,2.5,2.5],[2.5,2.5,-2.5]])
#print (get_intersection(edge_radius, edge_vector, face_points))



#points = np.array([[1,0,0],[0,2,0],[0,0,3]])
#print(make_supertetrahedron(points).points())



#point = np.array([3,3,4.000111])
#c = np.array([3,3,3])
#R = 1
#print(is_inside_sphere(point, c , R))

#triangle_points = np.array([[5,99,3],[2,-3,-2],[9,11,3]])
#print(calc_triangle_circumcenter(triangle_points), calc_triangle_circumradius(triangle_points))









#sort_by_dist(np.arange(30).reshape(10,3))

































