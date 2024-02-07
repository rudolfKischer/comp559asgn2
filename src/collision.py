# Rudolf C. Kischer
# 260956107

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
      force = np.array([force_viz_scale*(c._lambda[0] * c.n + c._lambda[1] * c.t1 + c._lambda[2] * c.t2) for c in self.contacts])
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

    # when we do this check, lets build build big J


    self.contacts = []
    for i in range(len(rigid_body_list)):
      self.check_ground(rigid_body_list[i])
      for j in range(i+1, len(rigid_body_list)):
        self.check_body_pair(rigid_body_list[i], rigid_body_list[j])
    return len(self.contacts) > 0


  def process(self, rigid_body_list, mu):


    #TODO: implement this function
    # mu = friction coefficient

    # u_{t+1} = u_t + M^-1 * J.T (delta t * lambda) + delta t * M^-1 * f

    # PGS = Projected Gauss Seidel


    # A * lambda + b = w
    # A = effective mass matrix
    # lambda = delta t * lamda
    # b = J(v + delta t * M^-1 * f) 
    # v = generalized velocity
    # f = force
    # J = constraint jacobian

    # bounce
    # restitution = epsilon
    # 0 <= epsilon <= 1
    # J_{k, row} * u_k^{t+1} = epsilon * J_{k, row} * u_k^t = b_k
    # b_bounce = [b_1 0 0 0... b_k 0 0 0...]

    # we want to create a matrix A
    # and we want to create a vector b
    # and we want to create a vector lambda

    # A = J * M^-1 * J.T
    # J is the jacobian of all the contacts
    # J is in R^{3K x 6N}
    # for the kth contact point between i and j
    # all blocks of the jth row of J is zero except for the sub blocks
    # J_k,2i = J
    for rb in rigid_body_list:
      rb.delta_V = np.zeros(6)

    for c in self.contacts:
      # compute delta_V
      delta_V = c.compute_delta_V(100, mu)
      # delta V is a 12 x 1 row vector
      u1_delta = delta_V[0:6]
      u2_delta = delta_V[6:12]
      # reshape to row vector
      u1_delta = u1_delta.reshape(1,-1)
      u2_delta = u2_delta.reshape(1,-1)
      # convert to a 1d row vector
      u1_delta = np.squeeze(u1_delta)
      u2_delta = np.squeeze(u2_delta)
      c.body1.delta_V += u1_delta
      c.body2.delta_V += u2_delta
    
    for rb in rigid_body_list:
      rb.v += rb.delta_V[0:3]
      rb.omega += rb.delta_V[3:6]
      # rb.delta_V = np.zeros(6)

    return