# Rudolf C. Kischer
# 260956107

import numpy as np
import polyscope as ps
import igl
import json

def hat(u):
  return np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])

class RigidBody:
  def __init__(self, body_desc=None):
    if body_desc == None:
      self.name = "inertial frame"
      self.mass = 0 # total mass as a scalar
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
    self.omega_dot = np.zeros(3)

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
    # we need to update:
    # - linear velocity
    # - angular velocity
    # use the force and torque
    # torque = J * omega_dot + (omega_hat - J * omega) * omega

    omega_hat = hat(self.omega)
    # update linear velocity
    v_dot = h * self.mass_inv * self.force
    self.v += v_dot


    # update angular velocity
    # omega_dot = h * self.Jinv0 @ (self.torque - omega_hat @ self.J0 @ self.omega)

    # We will apply the effect of inertial forces seperately
    # this way er can account for errror , because we know that the angular velocity is not changing
    # so we can restore the magnitude of the angular velocity after the force is applied
    omega_dot = h* self.Jinv0 @ (-np.cross(self.omega, self.J0 @ self.omega))
    self.omega_dot = omega_dot
    self.omega_prev = self.omega.copy()
    self.omega += omega_dot  
    # normalize omega and multiply by omega prev magnitude
    self.omega = self.omega / np.linalg.norm(self.omega) * np.linalg.norm(self.omega_prev) 

    # apply torque after restoring distortions in the angular velocity
    self.omega += h * self.Jinv0 @ self.torque


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
    r = self.omega * h
    # cap velocity above 0.0001 to avoid division by zero
    if np.linalg.norm(r) > 0.0000001:
      theta = np.linalg.norm(r)
      r_hat = hat(r)
      rodrigues = np.eye(3) + np.sin(theta) * r_hat  + (1 - np.cos(theta)) * r_hat @ r_hat
      self.R = self.R @ rodrigues
    self.R = self.gram_schmidt(self.R)

    # U, _, Vt = np.linalg.svd(self.R)
    # # Ensure a right-handed coordinate system
    # if np.linalg.det(U @ Vt) < 0:
    #     U[:, -1] *= -1  # Change the sign of the last column of U
    # self.R = U @ Vt

    # print(self.R @ self.R.T)


    

    # rotation vector
    # R(omega, h) = exp(omega_hat * h)
    # I + omega_hat * sin(h) + omega_hat^2 * (1 - cos(h))
    # omega_hat = log(R) / h

    # J_{t+1} = R * J_t * R^T
    self.J0 = self.R @ self.J0 @ self.R.T
    
    # update inverse inertia tensor in world frame
    # Jinv_{t+1} = R * Jinv_t * R^T
    self.Jinv0 = self.R @ self.Jinv0 @ self.R.T

    # print(self.J0)
    # print(self.Jinv0)
    




    # correct using svd


    
    # update inertia tensor in world frame
    return


  def kinetic_energy(self):

    return 0.5 * (self.v.T @ (self.mass * np.eye(3)) @ self.v) + 0.5 * self.omega.T @ self.J0 @ self.omega
  
  def potential_energy(self, gravity):
    return self.mass * gravity[1] * self.x[1]
  
  def linear_momentum(self):
    return self.mass * self.v
  
  def angular_momentum(self):
    return self.J0 @ self.omega