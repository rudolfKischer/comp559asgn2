# TODO: YOUR NAME AND STUDENT NUMBER HERE

import numpy as np
import polyscope as ps
import igl
import json

class RigidBody:
	def __init__(self, body_desc=None):
		if body_desc == None:
			self.name = "inertial frame"
			self.mass = 0
			self.mass_inv = 0
			self.J0 = np.zeros((3,3))
			self.Jinv0 = np.zeros((3,3))
			self.x0 = np.zeros(3)
			self.R0 = np.eye(3)
			self.v0 = np.zeros(3)
			self.omega0 = np.zeros(3)
			self.reset()
		else:
			self.name = body_desc['name']
			data = json.load(open(body_desc['file_root']+".json"))
			self.V, _, _, self.F, _, _ = igl.read_obj(body_desc['file_root']+".obj")
			self.mass = data['mass']	# scalar
			self.J0 = np.array(data['J'], dtype=float)	
			self.mass_inv = 1.0 / self.mass
			self.Jinv0 = np.linalg.inv(self.J0)
			self.x0 = np.array(body_desc['x0'], dtype=float)
			self.R0 = np.array(body_desc.get('R0',np.eye(3)), dtype=float)
			self.v0 = np.array(body_desc.get('v0',(0,0,0)), dtype=float)
			self.omega0 = np.array(body_desc.get('omega0',(0,0,0)), dtype=float)
			self.reset()
			### Register the mesh
			# `verts` is a Nx3 numpy array of vertex positions
			# `faces` is a Fx3 array of indices, or a nested list
			self.ps_mesh = ps.register_surface_mesh(self.name, self.V, self.F, smooth_shade=False)
			self.update_display()

	def reset(self):
		self.x = self.x0.copy()
		self.R = self.R0.copy()
		self.v = self.v0.copy()
		self.omega = self.omega0.copy()
		self.force = np.zeros(3)
		self.torque = np.zeros(3)
		# TODO: keep track of rotational inertia in the world aligned frame!

	def update_display(self):
		# Construct and set the homogeneous transformation matrix
		T = np.eye(4)
		T[0:3,0:3] = self.R
		T[0:3,3] = self.x
		self.ps_mesh.set_transform( T )
		
	def zero_force_and_torque(self):
		self.force = np.zeros(3)
		self.torque = np.zeros(3)	

	def add_force(self, f):
		self.force += f

	def step_vel(self, h):
		# TODO: implement this function
		return

	def step_pos(self, h):	
		# TODO: implement this function
		return