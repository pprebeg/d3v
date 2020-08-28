import numpy as np
import copy 
#rm = np.array([[],[],[]])
	
def Rx(a):		#matrica rotacije oko osi x
	return np.array([[1,0,0],[0, np.cos(a), -np.sin(a)],[0, np.sin(a), np.cos(a)]])
	
def Ry(b):		#matrica rotacije oko osi y
	return np.array([[np.cos(b), 0, np.sin(b)],[0,1,0],[-np.sin(b), 0, np.cos(b)]])
	
def Rz(g):		#matrica rotacije oko osi z
	return np.array([[np.cos(g), -np.sin(g), 0],[np.sin(g), np.cos(g), 0],[0,0,1]])
		
def get_abs_angle(vector1, vector2):	
	angle = np.arccos(np.dot(vector1, vector2)/ (np.sum(vector1**2)**0.5 * np.sum(vector2**2)**0.5))
	#print(angle*180/np.pi)
	return angle
	
def get_vector_angles(vector):
	angles = np.full(3, 0.0)	#(a,b,g)
	i = np.array([1,0,0])
	j = np.array([0,1,0])
	k = np.array([0,0,1])
	#projekcije na osi
	vyz = copy.copy(vector)		#a
	vyz[0] = 0
	
	vxz = copy.copy(vector)		#b
	vxz[1] = 0
	
	vxy = copy.copy(vector)		#g
	vxy[2] = 0
	
	a = get_abs_angle(vyz, j)	#kut mjeren od j
	b = get_abs_angle(vxz, k)	#kut mjeren od k
	g = get_abs_angle(vxy, i)	#kut mjeren od i
	#print(vyz,vxz,vxy)
	#print(a*180/np.pi,b*180/np.pi,g*180/np.pi)
	#print("--------------")
	#print(j, vyz)
	#print(np.cross(j, vyz))		## nes tu ne valja
	normal_vyz = np.sign(np.cross(j, vyz))[0]
	#print(normal_vyz)
	
	normal_vxz = np.sign(np.cross(k, vxz))[1]
	#print(normal_vxz)
	
	normal_vxy = np.sign(np.cross(i, vxy))[2]
	#print(normal_vxy)
	#--
	if np.isclose(a, 0) or np.isclose(a, np.pi):
		angles[0] = a
	else:
		angles[0] = a * normal_vyz

	if np.isclose(b, 0) or np.isclose(b, np.pi):
		angles[1] = b
	else:
	#	print(b, normal_vxz,b * normal_vxz)
		angles[1] = b * normal_vxz		

	if np.isclose(g, 0) or np.isclose(g, np.pi):
		angles[2] = g
	else:

		angles[2] = g * normal_vxy
	#eliminacija negativnih kuteva
	#print(angles*180/np.pi)
	angles = np.where(angles < 0, angles + 2*np.pi, angles )		#where(condition, x, y) x stavlja ako je condition true, a y ako nije
	#print(angles*180/np.pi)
	return angles


#def rotate_points_to_fit_plane_normal(points, plane_normal):
#	centroid = np.sum(points, axis = 0) / points.shape[0]
#	points_normal = np.cross(points[1] - points[0], points[2] - points[0])
#	delta_angles = get_vector_angles(plane_normal) - get_vector_angles(points_normal)
#
#	#trebaju nam samo 2 rotacije oko drugacije osi koje nisu nan ili 0
#	#transformacija delta angla da im se nanovi zamijene sa 0 (nema rotacije) i s
#	delta_angles = np.where(np.isnan(delta_angles) == True, 0, delta_angles)
#	
#	
#	
#	
#	
#	while np.all(delta_angles != 0):		#dodaj uvijet ta se ne vrti u krug zauvjek break uvijet
#		Rotm = np.array([Rx(delta_angles[0]), Ry(delta_angles[1]), Rz(delta_angles[2])]) #matrice rotacije rx(a),ry(b),rz(g)
#		i = np.where(np.isnan(delta_angles) == False)
#
#		for j in i:
#			points = (np.squeeze(Rotm[j]) @ points.T).T
#			
#		points_normal = np.cross(points[1] - points[0], points[2] - points[0])
#		delta_angles = get_vector_angles(plane_normal) - get_vector_angles(points_normal)
#		print(delta_angles)
#	points = points + centroid
#	return points






	
	
	#angles = vector/
#ovo sve jos naruk provjeri
#print("---------------------")
#points = np.array([[1,1,0],[1,0,0],[0,0,0],[0,1,0]])	#.T zbog mnozenja matrica
#centroid = np.sum(points, axis = 0) / points.shape[0]
#triangle_normal = np.cross(points[1] - points[0], points[2] - points[0])
#plane_normal = np.array([1,0,0])	
#delta_angles = get_vector_angles(plane_normal) - get_vector_angles(triangle_normal)

#while np.all(delta_angles != 0):
#	Rotm = np.array([Rx(delta_angles[0]), Ry(delta_angles[1]), Rz(delta_angles[2])]) #matrice rotacije rx(a),ry(b),rz(g)
#	i = np.where(np.isnan(delta_angles) == False)
#	for j in i:
#		#print(Rotm[j])
#		points = (np.squeeze(Rotm[j]) @ points.T).T
	
	#print("points")
	#print(points)
	#v1 = points[1] - points[0]
	#v2 = points[2] - points[0]
	#print("vs")
	#print(v1)
	#print(v2)
#	triangle_normal = np.cross(points[1] - points[0], points[2] - points[0])
#	delta_angles = get_vector_angles(plane_normal) - get_vector_angles(triangle_normal)
	
#points = points + centroid
#print("orig222222222222222222222222222222222222222222222")
#print(points)
#print(np.cross(points[1] - points[0], points[2] - points[0]))
		#jos centroid
	




#points = np.array([[1,1,0],[1,0,0],[0,0,0],[0,1,0]])
#plane_normal = np.array([1,0,0])	
#points = rotate_points_to_fit_plane_normal(points, plane_normal)
#print(points)

#	points_normal = np.cross(points[1] - points[0], points[2] - points[0])		#pretpostavka je da su svi pointovi u istom planeu ili da ako nisu to  nije vazno
#	centroid = np.sum(points, axis = 0) / points.shape[0]
#	points = points - centroid
#	delta_angles = get_vector_angles(plane_normal) - get_vector_angles(points_normal)
#	
#	while np.all(delta_angles != 0):
#		Rotm = np.array([Rx(delta_angles[0]), Ry(delta_angles[1]), Rz(delta_angles[2])]) #matrice rotacije rx(a),ry(b),rz(g)
#		i = np.where(np.isnan(delta_angles) == False)
#		for j in i:
#			print(Rotm[j])
#			points = (np.squeeze(Rotm[j]) @ points.T).T
#	
#		print("points")
#		print(points)
#		points_normal = np.cross(points[1] - points[0], points[2] - points[0])
#		delta_angles = get_vector_angles(plane_normal) - get_vector_angles(points_normal)
#	
#	points = points + centroid
#	return points

	
		
#	def calc_a(self, vector):				#atan(z/y)
#		return np.arctan(vector[2]/vector[1])		
#	def calc_b(self, vector):				#atan(x/z)
#		return np.arctan(vector[0]/vector[2])
#	def calc_g(self, vector):				#atan(y/x)
#		return np.arctan(vector[1]/vector[0])
		
#problemi sa dijeljenem sa 0 u angle
#v1 = np.array([1,1,1])

#print(get_vector_angles(v1)*180/np.pi)


def Rm(t, u):	#t je theta(opci kut rotacije u 3d) u je jedinicni vektor oko kojeg rotiramo point(cross produkt dviju normala)
	ux = u[0]		#a1
	uy = u[1]		#a2
	uz = u[2]		#a3
	
	I = np.identity(3)
	uxx = np.array([[0, -uz, uy],[uz, 0, -ux],[-uy, ux, 0]])						# skew-symmetric matrix
	u_out_u = np.outer(u,u)															#outer product
	
	return np.cos(t)*I + np.sin(t)*uxx + (1-np.cos(t))*u_out_u 
	
def	rotate_point(point, centroid, rot_axis_vector, t):			#t is theta
	point = point - centroid		#translacija u centar
	u = rot_axis_vector / np.sum(rot_axis_vector**2)**0.5	#jedinicni vektor
	point = Rm(t, u) @ point
	point = point + centroid
	return point

def rotate_points_to_fit_plane_normal(points, plane_normal):
	centroid = np.sum(points, axis = 0) / points.shape[0]
	for i in range(points.shape[0]):
		points[i] = rotate_point(point, centroid, rot_axis_vector)
	
	
	
	
	points = points + centroid

#point = np.array([-1,0,0])
#centroid = np.full((3), 0)
#rot_axis_vector = np.array([0,1,0])
#t = np.pi
#print(rotate_point(point, centroid, rot_axis_vector, t))





from math import pi ,sin, cos

def R(theta, u):
    return [[cos(theta) + u[0]**2 * (1-cos(theta)), 
             u[0] * u[1] * (1-cos(theta)) - u[2] * sin(theta), 
             u[0] * u[2] * (1 - cos(theta)) + u[1] * sin(theta)],
            [u[0] * u[1] * (1-cos(theta)) + u[2] * sin(theta),
             cos(theta) + u[1]**2 * (1-cos(theta)),
             u[1] * u[2] * (1 - cos(theta)) - u[0] * sin(theta)],
            [u[0] * u[2] * (1-cos(theta)) - u[1] * sin(theta),
             u[1] * u[2] * (1-cos(theta)) + u[0] * sin(theta),
             cos(theta) + u[2]**2 * (1-cos(theta))]]

def Rotate(pointToRotate, point1, point2, theta):


    u= []
    squaredSum = 0
    for i,f in zip(point1, point2):
        u.append(f-i)
        squaredSum += (f-i) **2

    u = [i/squaredSum for i in u]

    r = R(theta, u)
    rotated = []

    for i in range(3):
        rotated.append(round(sum([r[j][i] * pointToRotate[j] for j in range(3)])))

    return rotated


#point = [1,0,0]
#p1 = [0,0,0]
#p2 = [0,1,0]

#print (Rotate(point, p1, p2, pi)) # [-1.0, 0.0, 0.0]











