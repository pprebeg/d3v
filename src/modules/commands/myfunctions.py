import openmesh as om
import numpy as np
import copy
import sys
import csv

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
	




		
		
def soft_merge_meshes(meshes):	#meshes je lista sa meshevima		
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




def hard_merge_meshes(meshes):   #meshes je lista sa meshevima
	merged_mesh = soft_merge_meshes(meshes)
	data = np.unique(merged_mesh.points(), return_counts = True, axis = 0)	
	unique_points = data[0]
	duplicate_counter = data[1]
	points_with_duplicates = unique_points[np.where(duplicate_counter > 1)]
	new_vh = []
	for dpoint in points_with_duplicates:
		new_vh.append(merged_mesh.add_vertex(dpoint))
			
	bad_fh_list = []
	new_fvi_list = []
	
	
	merged_mesh_fvi = merged_mesh.face_vertex_indices().tolist()
	merged_mesh_points = merged_mesh.points()
	for i in range(len(merged_mesh_fvi)):	#trazimo bad fh i mijenjamo njihov fvi u novi
		fvi = np.asarray(merged_mesh_fvi[i])		
		face_points = merged_mesh_points[fvi]
		new_fvi = copy.copy(fvi)
		for nvh in new_vh:		#trazi jeli novi vh i face imaju istu tocku
			new_point = merged_mesh_points[nvh.idx()] 
			
			new_fvi[array_where_equal(face_points, new_point)] = nvh.idx()
		if np.array_equal(fvi, new_fvi) == False:		#ako originalni i novi fvi nisu isti dodaje novi fvi u listu
			fh = merged_mesh.face_handle(i)
			bad_fh_list.append(fh)
			new_fvi_list.append(new_fvi)
			
			
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
			
	return merged_mesh	
							 
def delete_isolated_vertices(mesh): 
	mesh_points = mesh.points().tolist()
	mesh_vertex_face_indices = mesh.vertex_face_indices().tolist() #kod vertex_face_indices izolirani su oni kojima je svima -1 (arrajevi moraju biti svi istevelicine pa je null -1)
	for vh_idx in range(len(mesh_points)):
		
		neighbouring_faces_fh_idx = np.asarray(mesh_vertex_face_indices[vh_idx])
		if np.all(neighbouring_faces_fh_idx == -1):
			vh = mesh.vertex_handle(vh_idx)
			mesh.delete_vertex(vh)
	mesh.garbage_collection()
	
			
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
		
		
		
	return mesh

def make_deck(wline_points, subdivide = False):
	#clean duplicates
	wline_points = np.unique(wline_points, axis = 0) 
	central_points = np.empty((0,3))
	for point in wline_points:
		if np.isclose(point[1], 0) ==  False:
			central_point = copy.copy(point)
			central_point[1] = 0
			central_points = np.append(central_points, np.expand_dims(central_point, 0), axis = 0)
	
	deck_points = np.append(wline_points, central_points, axis = 0)
	
	
	w_max_index = wline_points.shape[0]
	c_max_index = central_points.shape[0]
	
	#pocetni i zadnji fvi koji nemogu u petlju
	deck_fvi = np.array([[0,0+w_max_index,1],[deck_points.shape[0]-1,w_max_index-1,w_max_index-2]])
	 
	
	for interval_index in range(len(central_points-1)):
		fvi = np.array([[w_max_index,w_max_index+1,1],[w_max_index+1,2,1]]) + interval_index
		deck_fvi = np.append(deck_fvi, fvi, axis = 0)
	
	if subdivide == False:
		return om.TriMesh(deck_points, deck_fvi)
	elif subdivide == True:
		return subdivide_mesh([om.TriMesh(deck_points, deck_fvi)])
	
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



def is_mesh_closed(mesh):	
	for eh in mesh.edges():		#check if mesh has any boundary edges if not mesh is closed, if yes mesh is open
		if mesh.is_boundary(eh) == True:
			return False
	else:
		return True


def subdivide_mesh(mesh_list, c = 0 ,n = 1): #face po face subdividamo n puta,c je counter
	if c < n:
		new_meshes = []
		for mesh in mesh_list:
			mesh_points =  mesh.points()
			mesh_fvi = mesh.face_vertex_indices().tolist()
			mesh_hei = mesh.face_halfedge_indices().tolist() #lista sa 3 vrijednosti unutra
			face_hevi = mesh.halfedge_vertex_indices().tolist() #heh idx -> vertex indices	#lista [.........] velika sa slistama od 2 point = points[vindices]
			for i in range(len(mesh_fvi)): #i je idx od fh
				face_points = np.empty((0,3))
				midpoints = np.empty((0,3))
				for j in mesh_hei[i]: #j je idx od heh za halfedgeve na tom faceu /za svaki halfedge handleidx u faceu
					
					hevi = (face_hevi[j])	#tu se vrti
					halfedge_points = mesh_points[hevi] #array 
					face_points = np.append(face_points, np.expand_dims(halfedge_points[0], axis = 0), axis = 0)	# da se zadrzi orijentacija
					midpoint = halfedge_points[0] + (halfedge_points[1] - halfedge_points[0]) * 0.5
					midpoints =  np.append(midpoints, np.expand_dims(midpoint, axis = 0), axis = 0)
				new_mesh = om.TriMesh()
				vhandles = []
				fhandles = []
				for point in np.append(face_points, midpoints, axis = 0):
					vhandles.append(new_mesh.add_vertex(point))
				
				
				fhandles.append(new_mesh.add_face(vhandles[0], vhandles[3], vhandles[5]))
				fhandles.append(new_mesh.add_face(vhandles[3], vhandles[1], vhandles[4]))
				fhandles.append(new_mesh.add_face(vhandles[5], vhandles[3], vhandles[4]))
				fhandles.append(new_mesh.add_face(vhandles[5], vhandles[4], vhandles[2]))
				new_meshes.append(new_mesh)
	
		return subdivide_mesh(new_meshes, c = c + 1, n = n)
	
	else:
		return hard_merge_meshes(mesh_list)
		
		
		

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
	
def make_block_csv_file():
	block_dims = np.array([1,1,1])
	mesh = make_block(block_dims)
	mesh = subdivide_mesh([mesh], n = 1)
	points = mesh.points().tolist()
	fvi = mesh.face_vertex_indices().tolist()

	with open("unit_block_points.csv", "w", newline = "") as csv_file:
		csv_writer = csv.writer(csv_file)
		for point in points:
			csv_writer.writerow([point[0],point[1],point[2]])
	
	with open("unit_block_fvi.csv", "w", newline = "") as csv_file:
		csv_writer = csv.writer(csv_file)
		for f in fvi:
			csv_writer.writerow([f[0],f[1],f[2]])
	
def	make_block_from_unit_csv(block_dims = np.array([1,1,1]), move_vector = np.array([0,0,0]), path = ""):
	with open(path + "unit_block_points.csv", "r", newline ="") as csv_file:
		csv_reader = csv.reader(csv_file)
		points = np.asarray([line for line in csv_reader]).astype(float)
		
	with open(path + "unit_block_fvi.csv", "r", newline ="") as csv_file:
		csv_reader = csv.reader(csv_file)
		fvi = np.asarray([line for line in csv_reader]).astype(int)
	

	return om.TriMesh(points * block_dims + move_vector, fvi)
	
def is_point_inside_form_mesh_y(point, form_fh_idx_to_check, form_mesh): #ovdje je moguc samo jedan intersection sa formom( ako idemo od xz ravnine)
	point_vector_j = np.array([0,1,0])
	form_mesh_points = form_mesh.points()
	form_mesh_fvi = form_mesh.face_vertex_indices().tolist()
	
	for fh_idx in form_fh_idx_to_check:
		fvi = form_mesh_fvi[fh_idx]
		face_points = form_mesh_points[fvi]
		if is_inside_triangle(point, face_points):	#jeli point na nekom faceu
			return (True, point)
		
		#print(point, face_points)
		intersection_point = get_intersection(point, point_vector_j, face_points)[0]
		#print(intersection_point)
		if intersection_point is not None:
			if abs(point[1]) <= abs(intersection_point[1]):  		#ako je abs(y) manji ili jednak od abs(IPy) : point lezi u meshu
				return (True, intersection_point) 
			else:
				return (False, intersection_point)
		
	else:
		return (None, point) 

	
def is_point_inside_form_mesh_x(point, form_fh_idx_to_check, form_mesh): #ovdje je moguc samo jedan intersection sa formom( ako idemo od xz ravnine)
	point_vector_x = np.array([1,0,0])
	form_mesh_points = form_mesh.points()
	form_mesh_fvi = form_mesh.face_vertex_indices().tolist()
	
	for fh_idx in form_fh_idx_to_check:
		fvi = form_mesh_fvi[fh_idx]
		face_points = form_mesh_points[fvi]
		if is_inside_triangle(point, face_points):	#jeli point na nekom faceu
			return (True, point)
		
		intersection_point = get_intersection(point, point_vector_x, face_points)[0]
		if intersection_point is not None:
			xi = intersection_point[0] 
			xp = point[0]
			delta = xi - xp
			fh = form_mesh.face_handle(fh_idx)
			f_normal = form_mesh.calc_face_normal(fh)
			if f_normal[0] < 0:	#face je na krmi 
				if delta < 0:
					return (True, point)
				else:
					return (False, intersection_point)
				
			elif f_normal[0] > 0: #face je na palubi
				if delta > 0:
					return (True, point)
				else:
					return (False, intersection_point)
			
	return (True, point)

	
#ekskluzivno za formu broda
def fit_block_to_form(block_mesh, block_dims, block_position, form_mesh):
	#1) koji facevi su blizu blocka
	form_mesh_fvi = form_mesh.face_vertex_indices().tolist()
	form_mesh_vfi = form_mesh.vertex_face_indices().tolist()
	form_mesh_ffi = form_mesh.face_face_indices().tolist()
	form_mesh_points = form_mesh.points()
	xmin = block_position[0] 
	xmax = block_position[0] + block_dims[0]
	ymin = block_position[1] 
	ymax = block_position[1] + block_dims[1]
	zmin = block_position[2]
	zmax = block_position[2] + block_dims[2]
	
	near_fh_idx_list = []
	
	#trazi fh_idx od svih faceva koji imaju tocku u projekciji sa blokm na xz ravninu
	for vh_idx in range(form_mesh_points.shape[0]):
		point = form_mesh_points[vh_idx]
		if (xmin <= point[0] <= xmax) and (ymin <= point[1] <= ymax) and (zmin <= point[2] <= zmax): #ako je point izmedju projekcije
			neighbouring_fh_idx = form_mesh_vfi[vh_idx]
			near_fh_idx_list += neighbouring_fh_idx 

	#micanje duplikata i -1 u listi
	near_fh_idx_list = set(near_fh_idx_list)
	near_fh_idx_list = list(near_fh_idx_list)
	try:
		near_fh_idx_list.remove(-1)
	except:
		pass
	#trazenje susjednih facea susjednim facevima
	extra_fh_idx = []
	for fh_idx in near_fh_idx_list:
		neighbouring_fh_idx = form_mesh_ffi[fh_idx]
		extra_fh_idx += neighbouring_fh_idx

	near_fh_idx_list += extra_fh_idx
	#micanje duplikata i -1 u listi
	near_fh_idx_list = set(near_fh_idx_list)
	near_fh_idx_list = list(near_fh_idx_list)
	try:
		near_fh_idx_list.remove(-1)
	except:
		pass
	
	#trazenje jos susjednih faceva
	extra_fh_idx = []
	for fh_idx in near_fh_idx_list:
		neighbouring_fh_idx = form_mesh_ffi[fh_idx]
		extra_fh_idx += neighbouring_fh_idx

	near_fh_idx_list += extra_fh_idx
	#micanje duplikata i -1 u listi
	near_fh_idx_list = set(near_fh_idx_list)
	near_fh_idx_list = list(near_fh_idx_list)
	try:
		near_fh_idx_list.remove(-1)
	except:
		pass
		
	#jesu li facevi na krmi ili palubi?
	fh = form_mesh.face_handle(near_fh_idx_list[0])
	normal = form_mesh.calc_face_normal(fh)
	normal = np.sign(normal[0])

	
	#micanje tocaka na formu po y:
	block_mesh_points = block_mesh.points()
	x_outside_points_vh_idx = []
	x_inside_points_vh_idx = []	
	for vh_idx in range(block_mesh_points.shape[0]):
		block_point = block_mesh_points[vh_idx]
		data = is_point_inside_form_mesh_y(block_point, near_fh_idx_list, form_mesh)	#ispitivanje po y osi

		if data[0] == False:
			intersection_point = data[1]
			vh = block_mesh.vertex_handle(vh_idx)
			block_mesh.set_point(vh, intersection_point)
			x_inside_points_vh_idx.append(vh_idx)
		elif data[0] == None:
			x_outside_points_vh_idx.append(vh_idx)
	
	if len(x_outside_points_vh_idx) != 0:	#ako ima outside tocaka
		block_mesh_points = block_mesh.points()
		#block_inside_points = block_mesh_points[x_inside_points_vh_idx]
		block_outside_points = block_mesh_points[x_outside_points_vh_idx]
		
		if normal > 0:	#pramac	min y na max x
			for vh_idx in x_outside_points_vh_idx:
				outside_point =  block_mesh_points[vh_idx]
				vh = block_mesh.vertex_handle(vh_idx)
				point_to_set_to = np.array([outside_point[0],block_position[1],outside_point[2]])
				block_mesh.set_point(vh, point_to_set_to)
		
		elif normal < 0: #krma	max y na min x
			block_inside_points = block_mesh_points[x_inside_points_vh_idx]
			block_zs = np.unique(block_mesh_points[:,2]) #trebam li ovo refreshat?
			for block_z in block_zs:
				block_inside_points_with_same_z = np.unique(block_inside_points[np.where(block_inside_points[:,2] == block_z)], axis = 0)
				block_outside_points_with_same_z = np.unique(block_outside_points[np.where(block_outside_points[:,2] == block_z)], axis = 0)
				
				min_x_points = block_inside_points_with_same_z[np.where(block_inside_points_with_same_z[:,0] == np.min(block_inside_points_with_same_z[:,0]))]
				max_y_point = min_x_points[np.where(min_x_points[:,1] == np.max(min_x_points[:,1]))][0]	#0 na kraju za slucaj da ih je vise na istoj tocki
				
				for vh_idx in x_outside_points_vh_idx:
					outside_point = block_mesh_points[vh_idx]
					if abs(outside_point[1]) > abs(max_y_point[1]): 
						vh = block_mesh.vertex_handle(vh_idx)
						point_to_set_to = np.array([outside_point[0],max_y_point[1],outside_point[2]])
						block_mesh.set_point(vh, point_to_set_to)
				
	
		block_mesh_points = block_mesh.points()
		for vh_idx in x_outside_points_vh_idx:
			block_point = block_mesh_points[vh_idx]
			data = is_point_inside_form_mesh_x(block_point, near_fh_idx_list, form_mesh)	#ispitivanje po x osi
			if data[0] == False:
				intersection_point = data[1]
				vh = block_mesh.vertex_handle(vh_idx)
				block_mesh.set_point(vh, intersection_point)		

			
			
			



