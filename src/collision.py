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
    # if len(self.contacts) > 0:
    #     print("contact:")
    #     print(self.contacts[-1].p, self.contacts[-1].n, self.contacts[-1].d)
    #     print("jacobian:")
    #     print(self.contacts[-1].compute_jacobian())
    #     print()
    return len(self.contacts) > 0
  


  def GSStep(self, A, b, _lambda, clamp_fn=lambda x: x):
    # this performs a single Gauss-Seidel step
    # goes through all rows of A and updates lambda
    # and then returns lambda
    # and so is the iteration
    # A * lambda + b = 0
    # are assuming that we are solving for 0
    # make sure to loop in w in b if you are solving for w
    
    # copy lambda

    # 1. Decompose A
    # 2. loop over all rows of A
    # 3. sum up all the linear combination of lambdas for that row, except for the diagonal
    # 4. divide by the diagonal element
    # 5. update lambda

    for i in range(A.shape[0]):
      # sum of all the linear combinations of lambdas for that row
      # except for the diagonal
      _sum = 0
      for j in range(A.shape[1]):
        if j != i:
          _sum += A[i][j] * _lambda[j]
      _lambda_i = (-_sum - b[i]) / A[i][i]
      _lambda[i] = clamp_fn(_lambda_i, i, _lambda)
    return _lambda
  
  def GSStep2(self, A, b, _lambda, clamp_fn=lambda x: x):
    # in this one we take advantage of our knowledge of the structure of A
    # we know that A is a block diagonal matrix
    # we will use a different formula

    # lambda_i ^{k+1} = lambda_i^k - ( b_i - A_{i,i} * lambda_i^k) / A_{i,i}
    # where k is the iteration number

    # we can now also notice that if we are solving for a velocity update we have that
    # lambda_i ^{k+1} = lambda_i^k - b_i' - J_rowi' * delta V
    # where J_rowi' is the ith row of J divided by A_{i,i}
    # and b_i' is b_i divided by A_{i,i}

    # T = M^-1 * J.T
    # delta V^{k+1} = delta V^k - T_col_i * lambda_i^k
    # where T_col_i is the ith column of T

    # lambda^k+1 = T lambda^k + c
    # where c = -(D + L)^-1 * b
    # where T = (D + L)^-1 * U
    pass



  
  def PGSSolve(self, A, b, _lambda, mu, iters=100):
    # A * lambda + b >= 0
    # lambda_lo <= lambda <= lambda_hi
    # Compute normal forces due to contact. You may have computed all 3 rows of the Jacobian of each
    # contact, but here you'll only use the normal direction. Write a PGS solver that does 100 iterations where each
    # iteration loops over all contacts. Compute the new Lagrange multiplier for the normal direction in the Gauss-
    # Seidel manner, then clamp it to zero if it is negative (i.e., acting like glue). Use Erleben's method of maintaining a
    # delta velocity of the bodies, that simplifies the computation of the new Lagrange multiplier, and is cheap to
    # update with any change of the Lagrange multiplier. Note that the easiest way to store this information is in each
    # body, rather than forming a vector for all bodies. Likewise, store the Lagrange multipliers for each contact in the

    # A = L + D + U

    # we want to be able to solve A * lambda + b = 0

    # we need to loop over all the contactis i in (1 , .. K) where K is the number of contacts

    # at each time step we update lambda_i

    # lambda_i ^ {k+1} = ( -Sum_{j=1}^{i-1} L_{i,j} lambda_j^{k+1} - Sum_{j=i+1}^{n-1} U_{i,j} lambda_j^k - b_i ) / D_{i,i}
    # where k is the iteration number
    # n is the number of contacts
    # L is the lower triangular matrix
    # D is the diagonal matrix
    # U is the upper triangular matrix

    # in summary, on iteration k, at row i, we update lambda_i
    # on the numerator we have:
    # -Sum_{j=1}^{i-1} L_{i,j} lambda_j^{k+1} : 
    #    here j is the column index and we are summing over all the columns from 1 to i-1
    #    we are multiplying the lambda_j^{k+1} by the lower triangular matrix L_{i,j}
    # - Sum_{j=i+1}^{n-1} U_{i,j} lambda_j^k - b_i 
    # this is the same, but for all the columns to the right of i, in row i
    # - b_i is the bias term
    # D_{i,i} is the diagonal term


    # after each iteration, we then need to clamp the corresponding lambdas

    # note that for every collision, we have a 3 x 12 jacobian matrix
    # composed of the normal and the two tangent vectors
    # row1 = normal force
    # row2 = frictional tangent force
    # row3 = frictional tangent force

    # each of these have a constraint equation
    # lambda_lo_i <= lambda_i <= lambda_hi_i
    # lambda_lo_i = [ 0, -mu * lambda_i_1, -mu * lambda_i_1]
    # lambda_hi_i = [ inf, mu * lambda_i_1, mu * lambda_i_1]

    # 0 <= lambda_i_1 <= inf
    # -mu * lambda_i_1 <= lambda_i_2 <= mu * lambda_i_1
    # -mu * lambda_i_1 <= lambda_i_3 <= mu * lambda_i_1

    # recall that A = J * M^-1 * J.T in our case
    # w = J * u
    # u = [v1, omega1, v2, omega2]
    # v = generalized velocity
    # omega = angular velocity

    # we want to solve for when w = 0
    # A * lambda + b = w

    # recall 
    # b = J(u.T + delta t * M^-1 * f)


    # also note that 
    # u_{t+1} = u_t + M^-1 * J.T (delta t * lambda) + delta t * M^-1 * f
    # u_{t+1} -  u_t = M^-1 * J.T (delta t * lambda) + delta t * M^-1 * f
    # delta u = M^-1 * J.T (delta t * lambda) + delta t * M^-1 * f

    # now that we have gone through all that, lets take a crack at it

    # we will implement a general PGS solver
    pass


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