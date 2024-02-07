# Rudolf C. Kischer
# 260956107

import numpy as np
import polyscope as ps
import igl
import json
import scipy as sp

from utils import hat

class RigidBody:
  def __init__(self, body_desc=None):
    if body_desc == None:
      self.name = "inertial frame"
      self.mass = 0
      self.mass_inv = 0
      self.J0 = np.zeros((3,3)) # inertia tensor
      self.Jinv0 = np.zeros((3,3)) # inverse of inertia tensor
      self.x0 = np.zeros(3) # position
      self.R0 = np.eye(3) # rotation matrix
      self.v0 = np.zeros(3) # linear velocity
      self.omega0 = np.zeros(3) # angular velocity
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
    # self.v = np.ones(3) # linear velocity
    # self.omega = np.ones(3) * 30
    self.omega = self.omega0.copy()
    self.J = self.J0
    self.Jinv = self.Jinv0
    self.delta_V = np.zeros(6)

    # TODO
    # self.J = self.J0


  def reset(self):
    self.x = self.x0.copy()
    self.R = self.R0.copy()
    self.v = self.v0.copy()
    self.omega = self.omega0.copy()
    self.force = np.zeros(3)
    self.torque = np.zeros(3)
    # TODO: keep track of rotational inertia in the world aligned frame!
    self.J = self.J0
    self.Jinv = self.Jinv0

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
    # we need to update:
    # - linear velocity
    # - angular velocity
    # use the force and torque
    # torque = J * omega_dot + (omega_hat - J * omega) * omega
    # self.v += self.delta_V[0:3]
    # self.omega += self.delta_V[3:6]
    self.J = self.R @ self.J0 @ self.R.T
    self.Jinv = self.R @ self.Jinv0 @ self.R.T


    # update linear velocity
    v_dot = h * self.mass_inv * self.force
    self.v += v_dot



    # update angular velocity
    # omega_dot = h * self.Jinv0 @ (self.torque - omega_hat @ self.J0 @ self.omega)

    # We will apply the effect of inertial forces seperately
    # this way er can account for errror , because we know that the angular velocity is not changing
    # so we can restore the magnitude of the angular velocity after the force is applied
    omega_hat = hat(self.omega)

    omega_dot = h* self.Jinv @ (self.torque - (omega_hat @ self.J @ self.omega))

    # omega_dot = h* self.Jinv @ (-(omega_hat @ self.J @ self.omega))
    # self.omega_prev = self.omega.copy()
    self.omega += omega_dot  

    # normalize omega and multiply by omega prev magnitude
    # if np.linalg.norm(self.omega) != 0:
        # self.omega = self.omega / np.linalg.norm(self.omega) * np.linalg.norm(self.omega_prev) 

    # apply torque after restoring distortions in the angular velocity
    # self.omega += h * self.Jinv @ self.torque


    return
  
  def gram_schmidt(self, A):
    B = np.zeros_like(A)
    for i in range(A.shape[1]):
      v = A[:,i]
      for j in range(i):
        v -= np.dot(B[:,j], A[:,i]) * B[:,j]
      B[:,i] = v / np.linalg.norm(v)
    return B


  def step_pos(self, h):	
    # TODO: implement this function
    # we need to update:
    # - position
    # - rotation matrix
    # - inertia tensor in world frame
    # - inverse inertia tensor in world frame
    self.x += h * self.v

    # update rotation matrix
    # r = self.omega * h
    # cap velocity above 0.0001 to avoid division by zero
    # if np.linalg.norm(r) > 0.00001:
    # theta = np.linalg.norm(r)
    # if theta != 0:
    #   r_hat = hat(r)
      # rodrigues = np.eye(3) + np.sin(theta) * r_hat  + (1 - np.cos(theta)) * r_hat @ r_hat
      # self.R = self.R @ rodrigues
    self.R = sp.linalg.expm(hat(self.omega) * h) @ self.R
    # self.R = self.gram_schmidt(self.R)


    

    # rotation vector
    # R(omega, h) = exp(omega_hat * h)
    # I + omega_hat * sin(h) + omega_hat^2 * (1 - cos(h))
    # omega_hat = log(R) / h

    # J_{t+1} = R * J_t * R^T

    
    # update inverse inertia tensor in world frame
    # Jinv_{t+1} = R * Jinv_t * R^T
    # TODO: use j_inv for M_inv

    return
  
  def generalized_mass_matrix(self):
    A1 = self.mass * np.eye(3)
    A2 = np.zeros((3,3))
    A3 = np.zeros((3,3))
    A4 = self.J

    M = np.block([[A1, A2], [A3, A4]])
    return M
  
  def inv_generalized_mass_matrix(self):
    # if the mass is zero, we can't invert the matrix
    if self.mass == 0:
      return np.zeros((6,6))
    M_inv = np.block([[np.eye(3) / self.mass, np.zeros((3,3))], [np.zeros((3,3)), self.Jinv]])
    return M_inv

  def kinetic_energy(self):
    return 0.5 * (self.v.T @ (self.mass * np.eye(3)) @ self.v) + 0.5 * self.omega.T @ self.J @ self.omega
  
  def potential_energy(self, gravity):
    return self.mass * -gravity[1] * self.x[1]
  
  def linear_momentum(self):
    return self.mass * self.v
  
  def angular_momentum(self):
    return self.J @ self.omega