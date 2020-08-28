import openmesh as om
import numpy as np
import copy
import rotm as rm

from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

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

def get_best_fitting_plane_normal(points):
	#points = np.array([[-4,-1,0],[3,-3,0],[-2,-2,0],[1,2,0],[0,0,0]]).T
	points = points.T		
	points = (points.T - np.sum(points,1) / points.shape[1]).T 
	svd = np.linalg.svd(points)

	#print(svd[0])
	plane_base_vector1 = svd[0][:,0]
	plane_base_vector2 = svd[0][:,1]
	plane_normal = svd[0][:,2]
	#calc_normal = np.cross(plane_base_vector2, plane_base_vector1)
	#print(calc_normal)
	#print(np.sum(calc_normal**2, axis = 0)**0.5) #length
	#print(calc_normal / (np.sum(calc_normal**2, axis = 0)**0.5))
	#print(calc/(np.sum(calc**2, axis = )))
	return plane_normal

#print(get_best_fitting_plane_normal(np.array([[0,0,0],[0,1,0],[1,1,0],[1,0,0]])))

def array_where_equal(a,b, bool = True):		#ako hocemo stupce samo stavi a.T , ako hocemo i gdje je razlicito stavimo bool = False
	i_array = np.empty(0, dtype = "int64")
	for i in range(a.shape[0]):
		if np.array_equal(a[i], b) == bool:
			i_array = np.append(i_array, i)
			
	return i_array


def flip_mesh_face_orientation(mesh):
	points = mesh.points()
	flipped_fvi = np.flip(mesh.face_vertex_indices(), axis = 1)
	return om.TriMesh(points, flipped_fvi)
	



def sort_by_dist(points, sort_by_min_dist = True):
	aritm_middle = np.sum(points, axis = 0)/points.shape[0]
	start_point = aritm_middle
	sorted_array_shape = list(points.shape)
	sorted_array_shape[0] = 0
	sorted_array = np.empty(sorted_array_shape)

	if sort_by_min_dist == True:
		for p in range(points.shape[0]):
			#print(start_point)
			#print("-------------")
			#print(points)
			d = start_point - points
			distance = np.sum(d**2, axis = 1)**0.5
			i = np.array(np.where(distance == distance.min())).flatten()[0]	#sto ako su 2 tocke jednako udaljene?	[0] da odabere samo 1
			#print("index" + str(i))
			start_point = points[i]
			points = np.delete(points, i, axis = 0)
			#print("star point"+str(start_point))
			sorted_array = np.append(sorted_array, np.expand_dims(start_point, axis = 0), axis = 0)
	
	elif sort_by_min_dist == False:		#trazimo po max duljini
		for p in range(points.shape[0]):
			#print(start_point)
			#print("-------------")
			#print(points)
			d = start_point - points
			distance = np.sum(d**2, axis = 1)**0.5
			i = np.array(np.where(distance == distance.max())).flatten()[0]	#sto ako su 2 tocke jednako udaljene?	[0] da odabere samo 1
			#print("index" + str(i))
			start_point = points[i]
			points = np.delete(points, i, axis = 0)
			#print("star point"+str(start_point))
			sorted_array = np.append(sorted_array, np.expand_dims(start_point, axis = 0), axis = 0)
		
	return sorted_array
	


		
		
def soft_merge_meshes(meshes):	#meshes je lista sa meshevima		#problem ako imamo 2 ista pointa u oba mesha: radi 2 verteksa
	points = np.empty((0,3))
	face_vertex_indices = np.empty((0,3))

	for mesh in meshes:
		if len(mesh.face_vertex_indices().shape) == 1:		#za slucaj da nema fvi u mesh
			add_fvi = np.empty((0,3))
		else:
			add_fvi = mesh.face_vertex_indices()
		face_vertex_indices = np.append(face_vertex_indices, add_fvi + points.shape[0], axis = 0)				#+points.chape[0] je tu da poreda face_vertex_indices sa njihovim indexom u novom arrayu
		points = np.append(points, mesh.points(), axis = 0)
	
	combined_mesh = om.TriMesh(points, face_vertex_indices)	
	
	return combined_mesh



#test divide points
#pretpostavka je da su im orijentacije prije merga dobre
def hard_merge_meshes(meshes):   #meshes je lista sa meshevima
	merged_mesh = soft_merge_meshes(meshes)
	data = np.unique(merged_mesh.points(), return_counts = True, axis = 0)	
	unique_points = data[0]
	duplicate_counter = data[1]
	points_with_duplicates = unique_points[np.where(duplicate_counter > 1)]
	
	#duplicate_point_idx = []
	#for idx in range(merged_mesh.points().shape[0]):
	#	for dpoint in points_with_duplicates:
	#		if np.array_equal(dpoint, merged_mesh.points()[idx]) == True:
	#			duplicate_point_idx.append(idx)
	
	new_vh = []
	for dpoint in points_with_duplicates:
		new_vh.append(merged_mesh.add_vertex(dpoint))
			
	bad_fh_list = []
	new_fvi_list = []
	for fh in merged_mesh.faces():		#trazimo bad fh i mijenjamo njihov fvi u novi
		fvi = merged_mesh.face_vertex_indices()[fh.idx()]		#fvi je array sa shape(3)
		#print("for fvi: "+ str(fvi) + "\n")
		face_points = merged_mesh.points()[fvi]
		new_fvi = copy.copy(fvi)
		for nvh in new_vh:		#trazi jeli novi vh i face imaju istu tocku
			new_point = merged_mesh.points()[nvh.idx()] 
			
			new_fvi[array_where_equal(face_points, new_point)] = nvh.idx()
			#print(new_fvi)
			#for i in range(face_points.shape[0]):
			#	if np.array_equal(face_points[i], new_point) == True:
			#		new_fvi[i] = nvh.idx()
			#		print(new_fvi)
		if np.array_equal(fvi, new_fvi) == False:		#ako originalni i novi fvi nisu isti dodaje novi fvi u listu
			#print("ping")
			bad_fh_list.append(fh)
			new_fvi_list.append(new_fvi)
			
		
	
	#for fh in bad_fh_list:
		#bad_fvi = merged_mesh.face_vertex_indices()[fh.idx()]
		#print(fh.idx())
	
	#print(new_fvi_list)					
	#print(merged_mesh.points())
	
	#print("points_with_duplicates:  " + str(points_with_duplicates) + "\n\n")
			
			
	for bad_fh in bad_fh_list:					#delete bad faces:
		merged_mesh.delete_face(bad_fh, False)  	#false da ne deletea izolirane vertexe
	
	merged_mesh.garbage_collection()
		
	
	for new_fvi in new_fvi_list:				#retriangularizacija sa novim fvi	
		new_face_vhandles = []
		for vi in new_fvi:
			new_vh = merged_mesh.vertex_handle(vi)
			new_face_vhandles.append(new_vh)
		
		merged_mesh.add_face(new_face_vhandles)
			
	delete_isolated_vertices(merged_mesh)
	merged_mesh.garbage_collection()			
			
	return merged_mesh	
							
def delete_isolated_vertices(mesh):
	for vh in mesh.vertices():	#delete isolated verices
		i = 0
		for fh in mesh.vf(vh):
			i+=1
		else:
			if i == 0: #nema neighbouring faceva
				mesh.delete_vertex(vh)

							
	#print("points_with_duplicates:  " + str(points_with_duplicates) + "\n\n")
		
	#trebam povezati point sa verteksima koje treba zamijeniti
	#svaki face koji sadrzi duplicate point dodajemo u bad triangles, imamo vec duplicate points, deletamo ih; retrianguliramo ih sa novim verteksima
	
	#for dpoint in points_with_duplicates: #za svaki duplicate point addamo novi vertex i njegov id, nalazimo id od verteksa koji dijele point sa duplicate pointom
	#	new_vh = merged_mesh.add_vertex(dpoint)
	#	new_idx = new_vh.idx()
	#	
	#	dpoint_idx = []
	#	for i in range(merged_mesh.points().shape[0])		:
	#		if np.array_equal(merged_mesh.points()[i], dpoint) ==  True:
	#			dpoint_idx.append(i)
		#dpoint_idx = np.where(merged_mesh.points() == dpoint) #tu greska, ne trazi ih dobro
		#print(str(dpoint_idx) + "aaaaaaaaaaaaaaaaaaaaa")
		#for old_idx in dpoint_idx:		#za svaki stari idx od duplicate pointova uzimamo njegov stari vertex handle
		#	old_vh = merged_mesh.vertex_handle(old_idx)
			
		#	for old_fh in merged_mesh.vf(old_vh): 	#za svaki face koji dijeli stari vertex handle; radimo novi fvi tako da starom zamijenimo stari idx duplicate pointa sa novim
		#		old_fvi = merged_mesh.face_vertex_indices()[old_fh.idx()]
		#		new_fvi = old_fvi
		#		new_fvi = np.where(old_fvi == old_idx, new_idx, new_fvi)	#tu greska	#mogu li zadat face preko fvi? 
		#		
		#		new_face_vhandles = []
		#		for vi in new_fvi:  #pretvaramo new fvi u new vertex list za unos u add face
		#			new_face_vhandles.append(merged_mesh.vertex_handle(vi))
		#		merged_mesh.delete_face(old_fh, False) #ne deleta izolirane vertekse; njih na samom kraju
		#		merged_mesh.add_face(new_face_vhandles)		#mogu li kao listu? ili [0][1][2]
		#	merged_mesh.garbage_collection()
		#merged_mesh.garbage_collection()
	
		
	

				
				
	
	
	
	
	
	
	
	
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


	
def is_inside_triangle(point, triangle_points):
	ABC = triangle_surface(triangle_points)
	Ai = np.empty((0))
	for i in range(triangle_points.shape[0]):
		points = copy.copy(triangle_points)
		points[i] = point
		Ai = np.append(Ai, triangle_surface(points))
	
	if np.isclose(np.sum(Ai), ABC) == True:
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

def make_help_axes(scale = 5.):	# longest is x, 2nd longest is y, shortest is z
	axes_mesh = om.TriMesh()
	x_points = (np.array([[0,0,0],[2,0,0],[0,0.1,0]]) * scale) - np.array([1,1,1])
	y_points = (np.array([[0,0,0],[0,1,0],[0.1,0,0]]) * scale) - np.array([1,1,1])
	z_points = (np.array([[0,0,0],[0,0,0.5],[0.1,0,0]]) * scale) - np.array([1,1,1])
	x_vhandles = []
	y_vhandles = []
	z_vhandles = []
	
	for x in x_points:
		x_vhandles.append(axes_mesh.add_vertex(x))
	for y in y_points:
		y_vhandles.append(axes_mesh.add_vertex(y))	
	for z in z_points:
		z_vhandles.append(axes_mesh.add_vertex(z))
	
	x_fh = axes_mesh.add_face(x_vhandles[0], x_vhandles[1], x_vhandles[2])
	y_fh = axes_mesh.add_face(y_vhandles[0], y_vhandles[1], y_vhandles[2])	
	z_fh = axes_mesh.add_face(z_vhandles[0], z_vhandles[1], z_vhandles[2])
	
	return axes_mesh
	
	
	
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

def make_form(scale = 5, move_vector = np.array([0.,0.,0.])):
	mesh = om.TriMesh()
	points = np.array([[2.,1.,2.],[2.,2.,1.]]) * scale
	vhandles = []
	fhandles = []
		
	#xy plane mirrored points
	mirrored_points = copy.deepcopy(points)
	mirrored_points[:,2] *= -1
		
	points = np.append(points, mirrored_points, axis = 0)
	
	#yz plane mirrored points
	mirrored_points = copy.deepcopy(points)
	mirrored_points[:,0] *= -1
		
	points = np.append(points, mirrored_points, axis = 0)
	points += move_vector
	
	for point in points:
		vhandles.append(mesh.add_vertex(point))
		

	fhandles.append(mesh.add_face(vhandles [1], vhandles [4], vhandles [0]))
	fhandles.append(mesh.add_face(vhandles [1], vhandles [5], vhandles [4]))
	fhandles.append(mesh.add_face(vhandles [1], vhandles [3], vhandles [5]))
	fhandles.append(mesh.add_face(vhandles [3], vhandles [7], vhandles [5]))
	fhandles.append(mesh.add_face(vhandles [3], vhandles [2], vhandles [7]))
	fhandles.append(mesh.add_face(vhandles [2], vhandles [6], vhandles [7]))
	
	return mesh
	
	

#def is_outside_mesh(point, mesh):
#	outside_check_array = np.empty((0))
#	for fh in mesh.faces():
#		fvi = mesh.face_vertex_indices()[fh.idx()]
#		face_points = mesh.points()[fvi]
#		face_radius = face_points[0]
#		face_vector1 = face_points[1] - face_points[0]
#		face_vector2 = face_points[2] - face_points[0]
#		face_normal = mesh.calc_face_normal(fh)
#		
#		radius_matrix = (point - face_radius).T
#		vector_matrix = np.empty((3,3))
#		vector_matrix[:,0] = face_vector1
#		vector_matrix[:,1] = face_vector2
#		vector_matrix[:,2] = face_normal
#		
#		try:			#sto ako je null face?
#			normal_parameter = np.linalg.solve(vector_matrix, radius_matrix)[2]
#			if normal_parameter > 0:	#if point is outside mesh face
#				outside_check_array = np.append(outside_check_array, True)
#			else:
#				outside_check_array = np.append(outside_check_array, False)
#		
#		
#		except:		#ako su v1 i v2 na istom pravcu
#			outside_check_array = np.append(outside_check_array, False)
#		
#	if np.all(outside_check_array) == True:
#		return True
#	else:
#		return False
		
#def get_outside_points(mesh, boundary_mesh): #izbacuje pointove mesha koji su sa pozitivne strane svih faceva u boundary meshu#
#	outside_points = np.empty((3,3))
#	for point in mesh.points():
#		if is_outside_mesh(point, boundary_mesh) == True:#
#			outside_points = np.append(outside_points, np.expand_dims(point, axis = 0), axis = 0)
#	
#	return outside_points
	
def move_mesh(mesh, move_vector):
	for vh in mesh.vertices():
		new_point = mesh.points()[vh.idx()] + move_vector
		mesh.set_point(vh, new_point)
	return mesh
		
	
def get_intersection(radius, vector, face_points):
	face_radius = face_points[0]
	face_vector1 = face_points[1] - face_points[0]
	face_vector2 = face_points[2] - face_points[0]
	
	#face_radius, face_vector1, face_vector2,
	
	radius_matrix = (radius - face_radius).T
	vector_matrix = np.empty((3,3))
	vector_matrix[:,0] = face_vector1
	vector_matrix[:,1] = face_vector2
	vector_matrix[:,2] = -vector
		
	try:
		edge_parameter = np.linalg.solve(vector_matrix, radius_matrix)[2]
	except:
		return (None, None)
	
	intersection_point = radius + (vector * edge_parameter)
	if is_inside_triangle(intersection_point, face_points) == True:
		return (intersection_point, edge_parameter)
	else:
		return (None, edge_parameter)

		
		
		
		
def get_meshes_intersections(mesh1, mesh2):
	intersection_points = np.empty((0,3))
	for eh in mesh1.edges():
		edge_vector = mesh1.calc_edge_vector(eh)
		evi = mesh1.edge_vertex_indices()[eh.idx()]		#evi je edge vertex indices (array sa idovima pointova za taj edge)
		edge_radius = mesh1.points()[evi][0]			#za radius nam treba samo jedna tocka 
		
		for fh in mesh2.faces():
			fvi = mesh2.face_vertex_indices()[fh.idx()]
			face_points = mesh2.points()[fvi]
			
			int_data= get_intersection(edge_radius, edge_vector, face_points)
			intersection_point = int_data[0]
			int_parameter = int_data[1]
			if (intersection_point is not None) and (0 <= int_parameter <= 1):		#uvijet ako intersection point nije u trokutu ili edge i face su paralelni, int point je izvan edga:
				intersection_points = np.append(intersection_points, intersection_point.reshape(1,3), axis = 0)
			
	for eh in mesh2.edges():
		edge_vector = mesh2.calc_edge_vector(eh)
		evi = mesh2.edge_vertex_indices()[eh.idx()]		
		edge_radius = mesh2.points()[evi][0]			
		
		for fh in mesh1.faces():
			fvi = mesh1.face_vertex_indices()[fh.idx()]
			face_points = mesh1.points()[fvi]
			
			int_data= get_intersection(edge_radius, edge_vector, face_points)
			intersection_point = int_data[0]
			int_parameter = int_data[1]
			if (intersection_point is not None) and (0 <= int_parameter <= 1):		#uvijet ako intersection point nije u trokutu ili edge i face su paralelni, int point je izvan edga:	
				intersection_points = np.append(intersection_points, intersection_point.reshape(1,3), axis = 0)
			
	return np.unique(intersection_points, axis = 0) 					#mice sve duplikate tocaka



	#rot_axis vector nam moze bit 0,0,0, face subdividane forme ima 2 iste tocke!!!!!!! provijeri
def make_supertriangle(points, report):			#problem u best fitting plane normali; rpevise dimenzija
	report.write("\n\nmake_supertriangle_input points:\n\n" + str(points) + "\n\n")
	supertriangle= om.TriMesh()
	supertriangle_points = np.array([[1,0,0],[-1/2,3**0.5/2,0],[-1/2,-3**0.5/2,0]])
	supertriangle_centroid =  np.sum(supertriangle_points, axis = 0) / supertriangle_points.shape[0]
	supertriangle_normal = np.cross(supertriangle_points[1]-supertriangle_points[0], supertriangle_points[2]-supertriangle_points[0])
	supertriangle_normal = supertriangle_normal / (np.sum(supertriangle_normal**2))**0.5	#jed vektor
	points_normal = get_best_fitting_plane_normal(points)
	points_normal = points_normal / (np.sum(points_normal**2))**0.5		#jed vektor
	rot_axis_vector = np.cross(supertriangle_normal, points_normal)		#ako je rot axis vector 0,0,0 onda su normale kolinearne	i ne treba rotirat 
	
	report.write("\n\nsupertriangle points:\n\n" + str(supertriangle_points) + "\n\n")
	report.write("\n\nsupertr centroid:\n\n" + str(supertriangle_centroid) + "\n\n")
	report.write("\n\nsupertr normal:\n\n" + str(supertriangle_normal) + "\n\n")
	report.write("\n\npoints normal:\n\n" + str(points_normal) + "\n\n")
	report.write("\n\nrot_axis_vector:\n\n" + str(rot_axis_vector ) + "\n\n")
	
	if np.allclose(rot_axis_vector, 0) == True:
		t = 0
	else:
		rot_axis_vector = rot_axis_vector / (np.sum(rot_axis_vector**2))**0.5	#jedinicni vektor	#ako rot axis vektor 0,0,0; jednicni vektor je nan; zato ga racunamo nakon else
		t = np.arccos(np.dot(supertriangle_normal, points_normal))	#abs kut izmedju vektora
		#provjera predznaka kuta t
		test_normal = copy.copy(supertriangle_normal)
		test_normal = rm.rotate_point(test_normal, supertriangle_centroid, rot_axis_vector, t)
		if np.allclose(test_normal, points_normal) == False:	#ako se normale nisu poklopile
			t = -t
		#rotacija supertrokuta:
		supertriangle_points[0] = rm.rotate_point(supertriangle_points[0], supertriangle_centroid, rot_axis_vector, t)
		supertriangle_points[1] = rm.rotate_point(supertriangle_points[1], supertriangle_centroid, rot_axis_vector, t)
		supertriangle_points[2] = rm.rotate_point(supertriangle_points[2], supertriangle_centroid, rot_axis_vector, t)
	
	report.write("\n\ntheta:\n\n" + str(t) + "\n\n")
	report.write("\n\nrotated supertriangle points:\n\n" + str(supertriangle_points) + "\n\n")
	#supertriangle_points = rm.rotate_points_to_fit_plane_normal(supertriangle_points, points_normal)	ovo kasnije dodaj za sad je dobro ovak 
	
	
	arithmetic_middle = np.sum(points, axis = 0) / points.shape[0]
	r = np.max(np.sum((points - arithmetic_middle)**2, axis = 1) ** 0.5)  #radius je najveci delta r or svih pointova  
	R = 2*r*1000	#sigurnost od 1.05		## ako se pogodi tako da je tocka na pravcu izedju supertr pointova algoritam nece dobro funkcionirat
	#R = 2*r*1.0
	supertriangle_points = arithmetic_middle + supertriangle_points * R
	report.write("\n\ntranslated and scaled supertr points:\n\n" + str(supertriangle_points) + "\n\n")
	
	fhandles = []
	vhandles = []
	for point in supertriangle_points:
		vhandles.append(supertriangle.add_vertex(point))
	fhandles.append(supertriangle.add_face(vhandles[0],vhandles[1],vhandles[2]))
	
	return supertriangle
	
   
   
def BowyerWatson_triangulation_algorithm_backup(points, report):	#jos ga muci nes; ako su tocke kao u primjeru badT su svi facevi i brise ih ali ne ostavlja edgeve sa prosim pointom iteracije (svi facevi su u bad Tri)
	#da probam max r po max r od pointa do pointa?
	points = sort_by_dist(points, True)
	report.write("input points: \n" + str(points) + "\n\n")
	mesh = make_supertriangle(points, report)
	start_points = copy.copy(mesh.points())
	
	
	
	for point in points:
		#print(point) 
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
				#print(badTriangles.edge_vertex_indices()[eh.idx()])
				#print(badTriangles.is_boundary(eh))
				if badTriangles.is_boundary(eh) == True:
					evi = badTriangles.edge_vertex_indices()[eh.idx()]#neznam dali je u meshu i badtrianglesu eh isti za svaki edge najvj ne jer je generiran iz tocaka na pocetku, ali ako je generiran iz istih tocaka mozda su vertici isti(por. tocaka)
					#print(evi)
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
		#for vh in mesh.vertices():
			#print("aaaaa"+str(mesh.points()[vh.idx()])+"aaaaaa")
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


def is_mesh_closed(mesh):	
	for eh in mesh.edges():		#check if mesh has any boundary edges if not mesh is closed, if yes mesh is open
		if mesh.is_boundary(eh) == True:
			return False
	else:
		return True

def is_point_on_edge_or_vertex(point, mesh):	# fh provjerava hogu li dobiti taj point samo pomocu jednog vektora tog facesa
	eh_list = []
	fh_list = []
	vh_list = []
	for fh in mesh.faces():
		for eh in mesh.fe(fh):
			evi = mesh.edge_vertex_indices()[eh.idx()]
			edge_points = mesh.points()[evi]
			vhA = mesh.vertex_handle(evi[0])	#verteksi edga
			vhB = mesh.vertex_handle(evi[1])
			AB = edge_points[1] - edge_points[0]
			AC = point - edge_points[0]
			if np.isclose(np.sum(np.cross(AB,AC)**2)**0.5, 0) == True: #vektori su kolinearni ako im je len(cross) produkt 0
				DAC = np.dot(AB, AC)
				DAB = np.dot(AB, AB)
				if 0 <= DAC <= DAB:	#ako je dot produkt vektora AC izmedju 0  i DAB point je na edgu 
					if fh not in fh_list:
						fh_list.append(fh)		# ako je point neki verteks
					if eh not in eh_list:
						eh_list.append(eh)		 #ako eh nije u listi(sprijecavanje duplikata)
					if np.isclose(DAC, 0):		#ako je 0 A = C, ako je AB; C = B 
						if vhA not in vh_list:		#dodavanje verteksa u vh_list
							vh_list.append(vhA)
					if np.isclose(DAC, DAB):
						if vhB not in vh_list:
							vh_list.append(vhB)
					
					
	if len(eh_list) > 0:
		return (True, fh_list, eh_list, vh_list)
	else:
		return (False, fh_list, eh_list, vh_list)
	
def are_all_points_colinear(points):	#min 2 pointa u arrayu					#cross od dva 0,0,0 vectora je 0,0,0
	start_vector = points[1] - points[0]
	for point in points:	#svi cross produkti moraju biti 0,0,0
		vector = point - points[0]	#svi vektori su gledani od 1.pointa
		if np.allclose(np.cross(start_vector, vector), 0) == False:
			return False
	else:
		return True
				
				
			
#	for i in range(face_points.shape[0]): #za svai vektor u faceu
#		face_vector = face_points[i] - face_points[(i+1) % face_points.shape[0]] #% da loopa natrag na 0 
#		face_vector_radius = face_points[i]
#		vector_matrix = face_vector
#		radius_matrix = point - face_vector_radius
#		parameters = radius_matrix / vector_matrix		#ako su svi parametri isti i između 0 i 1 onda je unutar edga; nan ne ubrajamo
#		parameters = np.where(np.isnan(parameters) == True, parameters[np.where(np.isnan(parameters) == False)[0]], parameters) #sve nan-ove u parametrima pretvaramo u prvi ne nan broj i onda gledamo jesu li svi isti; nan predstavlja slucaj kada bilo  koji broj zadovoljava jdb
#		if (np.all(parameters == parameters[0]) and np.all(parameters >= 0) and np.all(parameters <= 1)) == True:					#jesu li svi parametri jednaki i između 0 i 1(ako us 0 ili 1 onda su na verteksu)
#			return True
#	else:	#ako se niti jedan vektor ne poklapa sa tockom vraca False
#		return False
			
#def is_point_a_vertex_in_mesh(point, mesh):	
#	vindices = array_where_equal(mesh.points(), point)
#	vh = []
#	if vindices.shape[0] == 0:	#ako nista nije dodano: return False
#		return (False,vh)
#	else:
#		for vi in vindices:
#			vh.append(mesh.vertex_handle(vi))	
#		return (True,vh)
			
	
	
	
	#vhandles = []	# za slucaj vise verteksa
	#for vh in mesh.vertices():
	#	mvp = mesh.points()[vh.idx()]
	#	if np.allclose(point, mvp):	#ako su oba arraya jako blizu
	#		vhandles.append(vh)
	#
	#if len(vhandles) == 0:		#ako nema vh koji se podudaraju
	#	return (False,)
	
	#else:
	#	if return_vh == False:
	#		return (True,)
	#	elif return_vh == True:
	#		return (True, vhandles)

		#ne funkcionira ako su intersectioni na edgu!!! ipak treba preko is on edge, vertex
	#pogreske ako intersecta u edgu ili verteksu
	#moguce rjesenje: ako se desi da je ip u nekom edgu ili verteksu; zarotiramo vektor i sve ispocetka?
	#jos bolje: postavimo vektor presijecanja da ide kroz centroid nekog facea; sanse da je mesh takav da cak i onda ga presijece na vektoru ili edgu je mala
def is_point_inside_closed_mesh(point, form_mesh):		#prvo gleda jeli point na facevima(ako je ignoriraj ostatak algoritma i dodaj point u inside) nakon toga broji intersectione sa facevima mesha. ( problemi ako se pogodi vertex ili edge jer se int point broji face po face i moze doci do previse int_countera (moguce rjesenje sa is in edge/vertex i dali ga point intersecta sa pozitivne stane normale ili negative(tesko za implementirat)))) 
	#arithmetic_middle =	np.sum(mesh.points(), axis = 0) / (mesh.points().shape[0])	
	#point_vector = point - arithmetic_middle
	#point_vector = np.random.rand(1,3)
	checked_fh = []
	intersection_list = [] #intersection list je lista sa ip koji su se vec desili
	int_counter = 0
	face_centroid = form_mesh.calc_face_centroid(form_mesh.face_handle(0))		#face centroid je centar facea sa indeksom 0
	point_vector = face_centroid - point
	#report.write("face0:\n " + str(form_mesh.points()[form_mesh.face_vertex_indices()[0]]))
	#report.write("\n\nface0 centroid:\n " + str(face_centroid))
	#prva provijera je jeli point vec lezi na jednom od faceva; ovo je odvojena petlja jel se nekad zna desit da face na kojem point lezi je dodan u checked_fh ako j pogodjen u edge koji ga granici
	for fh in form_mesh.faces():
		fvi = form_mesh.face_vertex_indices()[fh.idx()]
		face_points = form_mesh.points()[fvi]
		if is_inside_triangle(point, face_points):	#jeli point na na faceu; ako je automatski je point unutar mesha		#t mora bit pozitivan!!!!!!!!!
			return True	
	
	for fh in form_mesh.faces():
		#report.write("\n\nfor fh idx:\n " + str(fh.idx()))
		if fh not in checked_fh:
			checked_fh.append(fh)
			fvi = form_mesh.face_vertex_indices()[fh.idx()]
			face_points = form_mesh.points()[fvi]
			data1 = get_intersection(point, point_vector, face_points) 
			intersection_point = data1[0]
			intersection_parameter = data1[1]
			if (intersection_point is not None) and (intersection_parameter > 0) and (intersection_point.tolist() not in intersection_list):	#ako se taj intersection vec jedamput desio
				intersection_list.append(intersection_point.tolist())
				#report.write("\nintersection at:\n " + str(intersection_point))
				int_counter += 1
				data2 = is_point_on_edge_or_vertex(point, form_mesh)		#provijera jeli point na edgu ili verteksu
				if data2[0] == True:
					checked_fh += (data2[1])		#dodajemo checked_fh faceve koji granice sa tim pointom

					
	#print(int_counter)	
	if (int_counter % 2) == 1: # ako je counter neparan
		return True
	else:
		return False		
		
		
		
		
		
		
		
		
#		#ako je int point sve isponova sa novim random vectorom u nadi da nece pogodit vertex ili edge(is inside edge vec proverava jel u verteksu (parametar 0)
#		if intersection_point is not None:
#			if is_point_on_edge(intersection_point, face_points) == True:
#				return (is_point_inside_closed_mesh(point, form_mesh))	#rekurzija sa novim random vektorom
#			else:
#				int_counter += 1
#			
#	if (int_counter % 2) == 1: # ako je counter neparan
#		return True
#	else:
#		return False

#!!!!!!!!!!!!! pretpostavka je forma broda koja je simetricna s obziron na xz os!!!!!!!!! ne funkcionira pravilno za ostale slucajeve
def is_point_inside_open_mesh(point, form_mesh): #ovdje je moguc samo jedan intersection sa formom( ako idemo od xz ravnine)
	point_vector = np.array([0,1,0])	#1 j
	
	for fh in form_mesh.faces():
		fvi = form_mesh.face_vertex_indices()[fh.idx()]
		face_points = form_mesh.points()[fvi]
		if is_inside_triangle(point, face_points):	#jeli point na nekom faceu
			return True
		intersection_point = get_intersection(point, point_vector, face_points)[0]
		
		if intersection_point is not None:
			if abs(point[1]) <= abs(intersection_point[1]):  		#ako je abs(y) manji ili jednak od abs(IPy) : point lezi u meshu
				return True
			else:
				return False
		
		#ako face ima vise od 2 ip u njemu treba posebno ili jednostavno lupi BW algoritam(probaj sa ogromnim supertrokutom (rpovjeri ovo)) dodaj uvijet u supertrokutu da mu se orijentacija flipa tako da moze za svaku stranu napraviti mesh prave orijentacije
		#daje complex edgeve; pogledaj unos fvi, pogledaj uvijet za 2 inside pointa ne dijeli ih dobro!!!!!!!!!!!!!!! 2 pointa unutra

		
		
#funkcionira ;ali ako je inside points empty nema garancije da ce BW pogodit tocnu strukturu mesha; u tom slucaju treba subdividat mesh		
#ako je inside points empty; new face od intersection pointova? kod form inside pointova mislim da nije dobar is inside closed mesh provijeri
def slice_mesh(mesh, intersection_points, inside_points, outside_points, report):
	badTri = []
	badTri_normals = []
	new_face_points_list = []
	new_faces_meshes = []
	if inside_points.shape[0] == 0:   #specijalan slucaj kada nema inside pointova(npr. kada su svu intersection pointovi na jednom faceu) radimo novi mesh samo od intersection pointova / subdividamo mesh tako da ima inside pointova
		mesh = subdivide_mesh([mesh])
		for fh in mesh.faces(): #svi facevi su bad facevi
			badTri.append(fh)
		new_face_mesh = BowyerWatson_triangulation_algorithm_backup(intersection_points, report)
		new_faces_meshes.append(new_face_mesh)
	else:	#ako ima inside pointova:
		for fh in mesh.faces():		#kategorizacija fh i pointova
			fvi = mesh.face_vertex_indices()[fh.idx()]
			face_points = mesh.points()[fvi]
			report.write("\ncheck for face:\n " + str(fh.idx()) + "\nwith points:\n" + str(face_points))
			face_intersection_points = np.empty((0,3))
			face_inside_points = np.empty((0,3))
			face_outside_points = np.empty((0,3))
		
			for intersection_point in intersection_points:		#trazi jeli face los(ako ima intersection pointove u sebi)
				if is_inside_triangle(intersection_point, face_points) == True:
					face_intersection_points = np.append(face_intersection_points, np.expand_dims(intersection_point, axis = 0), axis = 0)
					report.write("\nintersection points in face: \n" + str(face_intersection_points) + "\n")
					
			if face_intersection_points.shape[0] != 0:	#ako face u sebi ima intersection point: posebni slucaj za samo 1 intersection point: mesh presijecen u verteksu pa ga ne ubrajamo u BadTri	and (face_intersection_points.shape[0] != 1)
				for face_point in face_points:		#trazimo koji pointovi facea su outside a koji inside 
					face_inside_points = np.append(face_inside_points, inside_points[array_where_equal(inside_points, face_point)], axis = 0)
					face_outside_points = np.append(face_outside_points, outside_points[array_where_equal(outside_points, face_point)], axis = 0)
				
				new_face_points = np.unique(np.append(face_intersection_points, face_inside_points, axis = 0), axis = 0)	#novi unique pointovi za face
				report.write("\nface int points are more than 1:------------------------------------------------------ \n" )
				report.write("\nface_inside_points\n: " + str(face_inside_points) + "\n")
				report.write("\nface outside_points\n: " + str(face_outside_points) + "\n")
				report.write("\nnew_face_points: \n" + str(new_face_points) + "\n")
							
				
				
				
				#provijera uvijeta: 1) min 3 pointa; nesmiju biti kolinearni 
				report.write("\nuvijeti:\n " )
				report.write("\nshape uvijet:\n " + str(new_face_points.shape[0] < 3 == False) + "\n")
				report.write("\nnew face points shape:\n " + str(new_face_points.shape[0]) + "\n")
				#report.write("\nkolinearity:\n " + str(are_all_points_colinear(new_face_points) == False) + "\n")
				
				if ((new_face_points.shape[0] < 3) == False) and (are_all_points_colinear(new_face_points) == False):
					bad_face_normal = mesh.calc_face_normal(fh)
					bad_face_normal = bad_face_normal / (np.sum(bad_face_normal**2))**0.5		#jedinicni vektor
					badTri.append(fh)
					new_face_points_list.append(new_face_points)
					badTri_normals.append(bad_face_normal)
			
			

		for i in range(len(badTri)): #izrada novih faceva
			face_mesh = BowyerWatson_triangulation_algorithm_backup(new_face_points_list[i], report)
			test_fh = face_mesh.face_handle(0) #zelimo neki f da mozemo izracunat normalu
			new_normal = face_mesh.calc_face_normal(test_fh)
			new_normal = new_normal / (np.sum(new_normal**2))**0.5	#jedninicni vektor
			if np.allclose(new_normal, badTri_normals[i]) == True:
				new_faces_meshes.append(face_mesh)
			else:		#ako orijentacija nije dobra flip orientation
				new_faces_meshes.append(flip_mesh_face_orientation(face_mesh))
		
		report.write("\nbadTri: " + str(badTri) + "\n")
		report.write("\nnew_face_points_list:"  + str(new_face_points_list) + "\n")
		report.write("\nnew_faces_meshes:"  + str(new_faces_meshes) + "\n")
					

	for fh in badTri:
		fvi = mesh.face_vertex_indices()[fh.idx()]
		face_points = mesh.points()[fvi]
		report.write("badface points\n" + str(fh.idx()) + "\n" + str(face_points) + "\n\n")
	
	
	#delete bad faces
	for fh in badTri:
		mesh.delete_face(fh, False)
	mesh.garbage_collection()
	
	#merge new faces with old mesh
	#new_faces_meshes.append(mesh)	#probaj promijeniti redoslied liste
	new_faces_meshes.insert(0, mesh)#	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	print("new faces meshesssssssssss" + str(new_faces_meshes)) 
	sliced_mesh = hard_merge_meshes(new_faces_meshes)
	

#ne cisti dobro outside faceve
	
	#delete outside vertices
	if inside_points.shape[0] != 0:	
		for outside_point in outside_points:
			vi = array_where_equal(sliced_mesh.points(), outside_point)
			if len(vi) != 0:
				vh = sliced_mesh.vertex_handle(vi)
				sliced_mesh.delete_vertex(vh)
		sliced_mesh.garbage_collection()
		
	delete_isolated_vertices(sliced_mesh)
	sliced_mesh.garbage_collection()
	print(len(new_faces_meshes))
	
	return sliced_mesh
	#return new_faces_meshes[0]
	
	
	
	
	
#	intersection_points = sort_by_dist(intersection_points, True)   #ovo mozda ne fukcionira u najopcenition slucaju; treba pri cuttanju sortat da su po redu ali nezz kako
#	sliced_mesh = soft_merge_meshes([mesh]) # sa samo jednim meshem je ovo efektivno copy mesh/ mislim da copy.copy ne funkcionira
#	badTriangles = []
#	new_fvi_list = []
#	
#	first_vh = sliced_mesh.add_vertex(intersection_points[0]) #inicijalizacija tako da vertekse ne dodaje dvaput nije dobra jel ce mi prvi vertex dodat dvaput
#	
#	report.write("sliced_mesh points: " + str(sliced_mesh.points()) + "\n\n")
#	for ip_index in range(intersection_points.shape[0]):
#		ip0 = intersection_points[ip_index]
#		ip1 = intersection_points[(ip_index+1) % intersection_points.shape[0]]
#		ipvh0 = sliced_mesh.vertex_handle(array_where_equal(sliced_mesh.points(), ip0)) #vraca vh za vh koji ima istu tocku kao ip[index]
#		
#		if ip_index == (intersection_points.shape[0] - 1):		#ako je loopao do kraja i ovo mu je sad prva tocka
#			ipvh1 = first_vh
#		else:
#			ipvh1 = sliced_mesh.add_vertex(ip1) #dodaje slijedecu tocku tao za se izbjegne dodavanje dvaput
#		
#		
#		
#		
#		ipvi0 = ipvh0.idx()
#		ipvi1 = ipvh1.idx()
#		
#		report.write("\n--------------------------------------------")
#		report.write("\nfor ip_index: " + str(ip_index) + "\n\n")
#		report.write("ip0 " + str(ip0) + " at index ipvi0 " + str(ipvi0) + "\n")
#		report.write("ip1 " + str(ip1) + " at index ipvi1 " + str(ipvi1) + "\n\n")
#		
#		
#		for fh in sliced_mesh.faces():
#			fvi = sliced_mesh.face_vertex_indices()[fh.idx()]
#			face_points = sliced_mesh.points()[fvi]
#
#		
#			
#			
#			if (is_inside_triangle(ip0, face_points)) == True and (is_inside_triangle(ip1, face_points)) == True:	#ako su u faceu obije tocke; dodajemo fh u badTri, trazimo koje su tocke u inside, a koje outisde i radimo novi fvi
#				if fh not in badTriangles: #za slucaj da je isti fh pogodio dva puta   neznam mogu li sa handlovima ovo; provijeri doma----mogu
#					badTriangles.append(fh)   
#				bad_face_normal = sliced_mesh.calc_face_normal(fh)
#				bad_face_normal = bad_face_normal / (np.sum(bad_face_normal**2))**0.5		#jedinicni vektor
#				report.write("\ntwo intersection points are inside in fh_index: " + str(fh.idx()) + "\n\n")
#				report.write("fvi " + str(fvi) + " with points:\n" + str(face_points) + "\n")
#				report.write("bad face normal " + str(bad_face_normal) + "\n")
#				#dijeli face pointove na outside i inside
#				face_inside_points = np.empty((0,3))
#				face_outside_points = np.empty((0,3))
#				for face_point in face_points: 
#					face_inside_points = np.append(face_inside_points, inside_points[array_where_equal(inside_points, face_point)], axis = 0)
#					face_outside_points = np.append(face_outside_points, outside_points[array_where_equal(outside_points, face_point)], axis = 0)
#			
#				report.write("inside points \n" + str(face_inside_points) + "\noutside points:\n" + str(face_outside_points) + "\n")
#			
#				if face_inside_points.shape[0] == 1:	#ako face ima samo jednu tocku inside
#					report.write("\nFace has only 1 point inside mesh\n") 
#					fp0 = face_inside_points[0]
#					fp0vi = array_where_equal(sliced_mesh.points(), fp0)[0]
#					new_points = [ip0, ip1, fp0]
#					new_fvi = [ipvi0, ipvi1, fp0vi]
#					vector1 = new_points[1] - new_points[0]
#					vector2 = new_points[1] - new_points[2]
#					new_normal = np.cross(vector1, vector2)
#					new_normal = new_normal / (np.sum(new_normal**2))**0.5	#jedninicni vektor
#					report.write("new_points " + str(new_points) + "\n with new fvi:\n" + str(new_fvi) + "\n")
#					report.write("\nNew normal: " + str(new_normal) + "\n")
#					if np.allclose(bad_face_normal, new_normal) == True:		#ako su im jedinicne normale iste znaci da su im i orijentacije tockaka iste
#						new_fvi_list.append(new_fvi)
#						report.write("\nNormals are the same; appending: \n" + str(new_fvi) + "\n")
#					else:
#						new_fvi_list.append(new_fvi[::-1])		#ako normale nisu jednake; mijenjamo redoslijed unosa tocki
#						report.write("\nNormals are	NOT the same; appending: \n" + str(new_fvi[::-1]) + "\n")
#					
#					report.write("\nnew fvi list after appending:\n" + str(new_fvi_list) + "\n")
#					
#					
#				elif face_inside_points.shape[0] == 2:	#ako face ima dvije tocke inside
#					report.write("\nFace has 2 points inside mesh")
#					ip_array = np.array([ip0, ip1])
#					ipvi_array = np.array([ipvi0, ipvi1])
#					fp0 = face_inside_points[0]
#					fp1 = face_inside_points[1]
#					fp0vi = array_where_equal(sliced_mesh.points(), fp0)[0]
#					fp1vi = array_where_equal(sliced_mesh.points(), fp1)[0]
#					
#					d = np.sum((ip_array - fp0)**2, axis = 1)**0.5
#					nearest_ip_id = np.where(d == d.min())[0]	#[0] ako su oboje jednako udaljeni
#					farthest_ip_id = (nearest_ip_id + 1) % ip_array.shape[0]	#vraca id druge tocke 
#					
#					new_points1 = np.array([fp0, ip_array[farthest_ip_id], ip_array[nearest_ip_id]])
#					new_points2 = np.array([fp0, fp1, ip_array[farthest_ip_id]])
#					new_fvi1 = [fp0vi, ipvi_array[farthest_ip_id][0], ipvi_array[nearest_ip_id][0]]
#					new_fvi2 = [fp0vi, fp1vi, ipvi_array[farthest_ip_id][0]]
#					report.write("CHECK THIS")
#					report.write("new_points1 " + str(new_points1) + "\n with new fvi1:\n" + str(new_fvi1) + "\n")
#					report.write("new_points2 " + str(new_points2) + "\n with new fvi2:\n" + str(new_fvi2) + "\n")
					
					

					
					#provjerava jeli odabir tocaka tocan	### tu pocinju problemi!!!!!!!		#ako su 2 ip na sredini a ne na edgevima moguc problem
					#A = triangle_surface(np.array([face_outside_points[0], fp0, fp1]))
					#A1 = triangle_surface(np.array([face_outside_points[0], ip0, ip1]))
					#A2 = triangle_surface(new_points1)
					#A3 = triangle_surface(new_points2)

					#if np.isclose(A,(A1+A2+A3)) == False:
					#	report.write("\nfvi were in the wrong order with: " + str(A) + " not equal to " + str(A1+A2+A3))
					#	new_points2 = np.array([fp0, fp1, ip1])
					#	new_fvi2 = [fp0vi, fp1vi, ipvi1]
					#	report.write("\nnew points 2 " + str(new_points2) + " new_fvi " + str(new_fvi2))
						
						
					#provijera orijentacije preko normala face1
#					vector11 = new_points1[1] - new_points1[0]
#					vector12 = new_points1[1] - new_points1[2]
#					new_normal1 = np.cross(vector11, vector12)
#					new_normal1 = new_normal1 / (np.sum(new_normal1**2))**0.5	#jedninicni vektor
#					report.write("\nNew normal1: " + str(new_normal1) + "\n")
#					if np.allclose(bad_face_normal, new_normal1) == True:		#ako su im jedinicne normale iste znaci da su im i orijentacije tockaka iste
#						new_fvi_list.append(new_fvi1)
#						report.write("\nNormals1 are the same; appending: \n" + str(new_fvi1) + "\n")
#					
#					else:
#						new_fvi_list.append(new_fvi1[::-1])		#ako normale nisu jednake; mijenjamo redoslijed unosa tocki 
#						report.write("\nNormals1 are NOT the same; appending: \n" + str(new_fvi1[::-1]) + "\n")
#					
#					#provijera orijentacije preko normala face2
#					vector21 = new_points2[1] - new_points2[0]
#					vector22 = new_points2[1] - new_points2[2]
#					new_normal2 = np.cross(vector21, vector22)
#					new_normal2 = new_normal2 / (np.sum(new_normal2**2))**0.5	#jedninicni vektor
#					report.write("\nNew normal2: " + str(new_normal2) + "\n")
#					if np.allclose(bad_face_normal, new_normal2) == True:		#ako su im jedinicne normale iste znaci da su im i orijentacije tockaka iste
#						new_fvi_list.append(new_fvi2)
#						report.write("\nNormals2 are the same; appending: \n" + str(new_fvi2) + "\n")
#					else:
#						new_fvi_list.append(new_fvi2[::-1])		#ako normale nisu jednake; mijenjamo redoslijed unosa tocki 
#						report.write("\nNormals2 are NOT the same; appending: \n" + str(new_fvi2) + "\n")

		#ipvh0 = mesh.add_vertex(intersection_points[ip_index])
		#ipvh1 = mesh.add_vertex(intersection_points[(ip_index+1) % intersection_points.shape[0]]) #dodaje prvu i slijedecu tocku kao vertex
		
#retriangulate
#	report.write("\nretriangularizaton:\n")
#	for fvi in new_fvi_list:#
#		vhandles = []
#		for vi in fvi:
#			vhandles.append(sliced_mesh.vertex_handle(vi))		#moze li vertex_handle podrzat array input ili mora jeno po jedno?
#		n_faces_before = len(sliced_mesh.faces())
#		sliced_mesh.add_face(vhandles)
#		n_faces_after = len(sliced_mesh.faces())
#		if n_faces_after == n_faces_before:
#			report.write("\n bad fvi with complex edge: " + str(fvi) + "\nwith points\n" + str(sliced_mesh.points()[fvi]))
#		else:
#			report.write("\ngood triangularization with: " + str(fvi) + "\n")
			

def subdivide_mesh(mesh_list, c = 0 ,n = 1): #face po face subdividamo n puta,c je counter
	if c < n:
		new_meshes = []
		for mesh in mesh_list:
			for fh in mesh.faces():
				face_points = np.empty((0,3))
				midpoints = np.empty((0,3))
				#i=0
				
				for heh in mesh.fh(fh):
					#i+=1
					#print(i)
					hevi = mesh.halfedge_vertex_indices()[heh.idx()]
					halfedge_points = mesh.points()[hevi] 
					face_points = np.append(face_points, np.expand_dims(halfedge_points[0], axis = 0), axis = 0)	# da se zadrzi orijentacija
					midpoint = halfedge_points[0] + (halfedge_points[1] - halfedge_points[0]) * 0.5
					#print(midpoint)
					midpoints =  np.append(midpoints, np.expand_dims(midpoint, axis = 0), axis = 0)
				#print(face_points)
				#print(midpoints)
				new_mesh = om.TriMesh()
				vhandles = []
				fhandles = []
				for point in np.append(face_points, midpoints, axis = 0):
					vhandles.append(new_mesh.add_vertex(point))
				
				
				#vhandles(0-2) = face points, vhandles(3-5) midpoints
				#print(np.append(face_points, midpoints, axis = 0))
				#print(len(vhandles))
				fhandles.append(new_mesh.add_face(vhandles[0], vhandles[3], vhandles[5]))
				fhandles.append(new_mesh.add_face(vhandles[3], vhandles[1], vhandles[4]))
				fhandles.append(new_mesh.add_face(vhandles[5], vhandles[3], vhandles[4]))
				fhandles.append(new_mesh.add_face(vhandles[5], vhandles[4], vhandles[2]))
				new_meshes.append(new_mesh)
	
		return subdivide_mesh(new_meshes, c = c + 1, n = n)
	
	else:
		return hard_merge_meshes(mesh_list)
		
		
		
		
#block mesh je mesh za koji gledamo jesu li pointovi izvan ili unutar form mesha
def divide_points_by_outside_mesh(block_mesh, form_mesh):	#2 slucaja: kada je mesh zatvoren i kada je mesh otvoren, pretpostavljamo da je forma i gledamo jeli unutra projekcije na xy os(centar broda)
	return_dict = {}
	outside_points = np.empty((0,3))
	inside_points = np.empty((0,3)) # ako je point na granici mesha onda se smatra da je unutar mesha
	closed_mesh = is_mesh_closed(form_mesh)
	if closed_mesh == True:
		for vh in block_mesh.vertices():
			point = block_mesh.points()[vh.idx()]
			if is_point_inside_closed_mesh(point, form_mesh) == True:
				inside_points = np.append(inside_points, np.expand_dims(point, axis = 0), axis = 0)
			else:
				outside_points = np.append(outside_points, np.expand_dims(point, axis = 0), axis = 0)
 
	
	
	elif closed_mesh == False:	#forma broda
		for vh in block_mesh.vertices():
			point = block_mesh.points()[vh.idx()]
			if is_point_inside_open_mesh(point, form_mesh) == True:
				inside_points = np.append(inside_points, np.expand_dims(point, axis = 0), axis = 0)
			else:
				outside_points = np.append(outside_points, np.expand_dims(point, axis = 0), axis = 0)

	return_dict["inside points"] = inside_points
	return_dict["outside points"] = outside_points
	return return_dict
	
				
def paint_mesh_faces(meshes_list):
	for mesh in meshes_list:
		print(mesh.has_face_colors())
		mesh.request_face_colors()
		for fh in mesh.faces():
			#n_faces = len(mesh.faces())
			#mesh.request_edge_colors()
			#mesh.request_face_colors()
			color = np.array([100,0,0,1])
			mesh.set_color(fh, color)
		print(mesh.has_face_colors())
	return meshes_list

def get_scatter(axes, x, y, z, color):	
	return axes.scatter3D(x, y, z, color = color);

def plot3D(points_list):
	scatter_points = []
	colors = ["red","green","yellow","blue","orange","purple","black","cyan","magneta"]
	fig=plt.figure()
	axes = plt.axes(projection='3d')
	axes.set_xlabel("x")
	axes.set_ylabel("y")
	axes.set_zlabel("z")
	
	for i in range(len(points_list)):
		points = points_list[i]
		color = colors[i % (len(colors)-1)]
		x = points[:,0]
		y = points[:,1]
		z = points[:,2]
		scatter_points.append(get_scatter(axes, x, y, z, color))

	#points1 = axes.scatter3D(points[:,0],points[:,1],points[:,2], color = "green"); 
	#points2 = axes.scatter3D(supertr_points[:,0],supertr_points[:,1],supertr_points[:,2], color = "red"); 
	#point_center = axes.scatter3D(points_centroid[0],points_centroid[1],points_centroid[2], color = "blue");
	#supertr_center = axes.scatter3D(supertr_centroid[0],supertr_centroid[1],supertr_centroid[2], color = "black");		
	plt.show()
	
   
			
#def cut_meshes2(block, form , report):					#outside points ne radi dobro
#	intersection_points = get_meshes_intersections(block, form)
#	outside_points = get_outside_points(block, form)
#	mesh = BowyerWatson_triangulation_algorithm_backup(intersection_points, report)
#	print(intersection_points)
#	print(outside_points)
#	print(block.points())
#	return mesh
	
#specijalizirano za faceve; unosimo 3 pointa; 4 point je 0,0,0 -> jednostavnija formula
def calc_face_volume(face_points):
	a = face_points[0]
	b = face_points[1]
	c = face_points[2]
	volume = np.abs(np.dot(a, np.cross(b,c))) / 6
	return volume
	
def calc_mesh_volume(mesh):	
	mesh_volume = 0
	for fh in mesh.faces():
		face_normal = mesh.calc_face_normal(fh)
		face_centroid = mesh.calc_face_centroid(fh)
		fvi = mesh.face_vertex_indices()[fh.idx()]
		face_points = mesh.points()[fvi]
		face_volume = calc_face_volume(face_points)
		face_sign = np.sign(np.dot(face_centroid, face_normal))
		mesh_volume += face_volume * face_sign
	return mesh_volume

def cut_meshes(block_mesh, form_mesh):
		report = open("report.txt", "w")
		objects = {}
		max_subdiv = 3
		for i in range(max_subdiv):
			intersection_points = get_meshes_intersections(block_mesh, form_mesh)
			#print(intersection_points)
			block_data = divide_points_by_outside_mesh(block_mesh, form_mesh)
			block_inside_points = block_data["inside points"]
			block_outside_points = block_data["outside points"]
			form_data = divide_points_by_outside_mesh(form_mesh, block_mesh)			
			form_inside_points = form_data["inside points"]
			form_outside_points = form_data["outside points"]
			if form_inside_points.shape[0] == 0 or block_inside_points.shape[0] == 0:
				block_mesh = subdivide_mesh([block_mesh])
				form_mesh = subdivide_mesh([form_mesh])
			else:
				break
			
			
		#print(form_inside_points)
		mesh1 = slice_mesh(block_mesh, intersection_points, block_inside_points, block_outside_points, report)		#complex edgevi
		mesh2 = slice_mesh(form_mesh, intersection_points, form_inside_points, form_outside_points, report)
		help_axes = make_help_axes()
		#return soft_merge_meshes([block_mesh, help_axes, form_mesh])
		#return block_mesh
		#return mesh2
		#plot3D([block_mesh.points(), form_mesh.points(), intersection_points])
		mesh = hard_merge_meshes([mesh1, mesh2])
		#for eh in mesh.edges():	
		#	if mesh.is_boundary(eh):
		#		mesh.delete_edge(eh, False)
		#		print("ping")
		#mesh.garbage_collection()
		return mesh
		#return mesh
	



	
	
#tests:
#mesh = om.TriMesh()
#mesh = make_block()
#points = np.array([[0,0,0],[1,0,0],[0,1,0]])
#test_point = np.array([1,1,1])
#vhandles = []
#for point in points:
#	vhandles.append(mesh.add_vertex(point))
#fh1 = mesh.add_face(vhandles[0],vhandles[1],vhandles[2])
#for vh in mesh.vertices():
#	print("idx \n" + str(vh.idx()) + "\n\n")
#	print("idx point\n" + str(mesh.points()[vh.idx()]) + "\n\n")
#	print("vvi\n" + str(mesh.vertex_vertex_indices()[vh.idx()]) + "\n\n")
#	vi = mesh.vertex_vertex_indices()[vh.idx()]		#vertex vertex indices izbacuje array sa dva inta
#	print("vvi points\n" + str(mesh.points()[vi]) + "\n\n")
	
#normal = mesh.calc_face_normal(fh1)
#print("normal:\n\n" + str(normal) + "\n\n")
#print(is_outside_mesh(test_point, mesh))
	

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



#print(sort_by_dist(np.arange(30).reshape(30,1), False))


#a = np.arange(16).reshape(4,4)
#b = a[2]
#print(array_where_equal(a,b,))

#points = np.array([[0,0,0],[1,0,0],[2,0,0],[-100,0,0]])
#print(are_all_points_colinear(points))




#mesh = make_block()
#print(mesh)
#print(calc_mesh_volume(mesh))















