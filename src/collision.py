# TODO: YOUR NAME AND STUDENT NUMBER HERE

import igl
import numpy as np
import polyscope as ps
from rigidbody import *
from contact import *

class Collision:
	def __init__(self):
		self.contacts = []
		self.ground_body = RigidBody()

	def reset(self):
		self.contacts = []

	def update_display(self, show_contacts):
		if len(self.contacts) == 0 or not show_contacts:
			self.ps_contacts = ps.register_point_cloud("contacts", np.zeros((0,3)))
		else:
			# can only update the points if they have the same number :(
			pos = np.array([c.p for c in self.contacts])
			depth = np.array([c.d for c in self.contacts])
			normal = np.array([c.n for c in self.contacts])
			t1 = np.array([c.t1 for c in self.contacts])
			t2 = np.array([c.t2 for c in self.contacts])
			force_viz_scale = 2
			force = np.array([force_viz_scale*(c.lamb[0] * c.n + c.lamb[1] * c.t1 + c.lamb[2] * c.t2) for c in self.contacts])
			self.ps_contacts = ps.register_point_cloud("contacts", pos)
			self.ps_contacts.add_scalar_quantity("contact depth", depth, enabled=True)	
			# self.ps_contacts.add_vector_quantity("contact normal", normal, enabled=True, radius=0.01, color=(0,0,1), vectortype='ambient')
			# self.ps_contacts.add_vector_quantity("contact t1", t1, enabled=True, radius=0.01, color=(1,0,0), vectortype='ambient')
			# self.ps_contacts.add_vector_quantity("contact t2", t2, enabled=True, radius=0.01, color=(0,1,0), vectortype='ambient')
			self.ps_contacts.add_vector_quantity("contact force", force, enabled=True, radius=0.01, color=(1,1,0), vectortype='ambient')	
				
	# collision check with ground plane
	def check_ground( self, body ):
		# check if any vertex is below the ground plane
		# if so, compute penetration depth and normal and add to contact list
		vt = body.V @ body.R.T + body.x
		for v in vt:
			if v[1] < 0: # y is up
				self.contacts.append( Contact( self.ground_body, body, v, np.array([0,1,0]), v[1] ) )
		
	def check_body_pair( self, body1, body2 ):
		# check if any vertex of one body is inside the other
		# NOTE: this is super gross because the signed distance function is expensive
		# thus we check the larger number of vertices with the smaller body
		# but WATCH OUT becaues this appears to be buggy in general.
		# For vericies inside the other body, compute penetration depth and normal
		v1t = body1.V @ body1.R.T + body1.x
		v2t = body2.V @ body2.R.T + body2.x
		if ( v1t.shape[0] > v2t.shape[0] ):
			S,I,C,N = igl.signed_distance( v1t, v2t, body2.F, return_normals=True )
			for i in range(len(S)):
				if S[i] < 0:
					self.contacts.append( Contact( body1, body2, C[i], -N[i], -S[i] ) )
		else:
			S,I,C,N = igl.signed_distance( v2t, v1t, body1.F, return_normals=True )
			for i in range(len(S)):
				if S[i] < 0:
					self.contacts.append( Contact( body2, body1, C[i], -N[i], -S[i] ) )

	def check( self, rigid_body_list ):
		self.contacts = []
		for i in range(len(rigid_body_list)):
			self.check_ground(rigid_body_list[i])
			for j in range(i+1, len(rigid_body_list)):
				self.check_body_pair(rigid_body_list[i], rigid_body_list[j])
		return len(self.contacts) > 0

	def process(self, rigid_body_list, mu):
		#TODO: implement this function
		return