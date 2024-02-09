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


  def reset(self):
    self.x = self.x0.copy()
    self.R = self.R0.copy()
    self.v = self.v0.copy()
    self.omega = self.omega0.copy()
    self.force = np.zeros(3)
    self.torque = np.zeros(3)
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
    self.J = self.R @ self.J0 @ self.R.T
    self.Jinv = self.R @ self.Jinv0 @ self.R.T


    # update linear velocity
    v_dot = h * self.mass_inv * self.force
    self.v += v_dot

    # update angular velocity
    omega_hat = hat(self.omega)
    omega_dot = h* self.Jinv @ (self.torque - (omega_hat @ self.J @ self.omega))
    self.omega += omega_dot  
    return


  def step_pos(self, h):	
    self.x += h * self.v
    self.R = sp.linalg.expm(hat(self.omega) * h) @ self.R
    return
  
  def generalized_mass_matrix(self):
    return np.block([[self.mass * np.eye(3), np.zeros((3,3))], [np.zeros((3,3)), self.J]])
  
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