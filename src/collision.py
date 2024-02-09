# Rudolf C. Kischer
# 260956107

import igl
import numpy as np
import polyscope as ps
from rigidbody import *
from contact import *

total_its = 0
total_collisions = 0

class Collision:
  def __init__(self):
    self.contacts = []
    self.prev_contacts = []
    self.ground_body = RigidBody()

  def reset(self):
    self.contacts = []
    self.prev_contacts = []

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
      self.ps_contacts.add_vector_quantity("contact normal", normal, enabled=True, radius=0.01, color=(0,0,1), vectortype='ambient')
      self.ps_contacts.add_vector_quantity("contact t1", t1, enabled=True, radius=0.01, color=(1,0,0), vectortype='ambient')
      self.ps_contacts.add_vector_quantity("contact t2", t2, enabled=True, radius=0.01, color=(0,1,0), vectortype='ambient')
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


  def pair_contacts(self, rigid_body_list):
    # we want to find the contacts that are the same from the previous frame
    # check that the p of the contacts are within some distance of each other
    # check that the normals are within some distance of each other

    unpaired_prev_contacts = self.prev_contacts.copy()

    contact_pairs = []

    epsilon = 0.001



    for c in self.contacts:
        closest_prev_contact = None
        min_distance = float('inf')

        for prev_c in unpaired_prev_contacts:
            # check that the bodies are the same
            if c.body1 != prev_c.body1 or c.body2 != prev_c.body2:
                continue
            
            # calculate distance
            distance = np.linalg.norm(c.p - prev_c.p) + np.linalg.norm(c.n - prev_c.n)
            
            # check if this is the closest so far
            if distance < min_distance and distance < epsilon:
                min_distance = distance
                closest_prev_contact = prev_c

        # Pair with the closest contact (regardless of the distance)
        if closest_prev_contact is not None:
            contact_pairs.append([c, closest_prev_contact])
            # remove the contact from the previous contacts
            unpaired_prev_contacts.remove(closest_prev_contact)
    return contact_pairs




  def process(self, rigid_body_list, mu, pgs_iterations):
    global total_its, total_collisions
    
    for rb in rigid_body_list:
      rb.delta_V = np.zeros(6)

    
    # # warm start the contacts
    contact_pairs = self.pair_contacts(rigid_body_list)
    for c, prev_c in contact_pairs:
      c._lambda = prev_c._lambda.copy()

    

    for c in self.contacts:
      # compute delta_V
      delta_V, k = c.compute_delta_V(pgs_iterations, mu)
      total_its += k
      total_collisions += 1
      # print(f'average its: {total_its / total_collisions}')
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

    
    self.prev_contacts = self.contacts.copy()
    
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

    # A = J * M^-1 * J.T
    # J is the jacobian of all the contacts
    # J is in R^{3K x 6N}
    # for the kth contact point between i and j
    # all blocks of the jth row of J is zero except for the sub blocks
    # J_k,2i = J

    return