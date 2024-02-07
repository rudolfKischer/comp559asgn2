# Rudolf C. Kischer
# 260956107

import numpy as np
from utils import hat

class Contact:
  # p is point of contact
  # n is normal
  # d is penetration depth
  def __init__(self, body1, body2, p, n, d):
    self.body1 = body1
    self.body2 = body2
    self.p = p
    self.n = n
    if (abs(n[0]) < abs(n[1] and abs(n[0]) < abs(n[2]))):
      tmp = np.array([1,0,0])
    elif (abs(n[1]) < abs(n[2])):
      tmp = np.array([0,1,0])	
    else:
      tmp = np.array([0,0,1])
    self.t2 = np.cross(self.n, tmp)
    self.t1 = np.cross(self.t2, self.n)
    self.d = d
    self._lambda = np.zeros(3)

    self.delta_lambda = np.zeros(3)



  def compute_jacobian(self):
    # TODO: implement this function
    
    # J = [-n.T, - (ri_hat * n).T, n.T, (rj_hat * n).T]
    # where r_i = contact point from body i
    # where r_j = contact point from body j
    # where ri_hat = the skew symmetric matrix of r_i
    # where rj_hat = the skew symmetric matrix of r_j
    # where n = the normal vector
    # r_i = self.p - self.body1.x
    # r_j = self.p - self.body2.x
    r_i = self.p - self.body1.x
    r_j = self.p - self.body2.x
    rj_hat = hat(r_j)
    ri_hat = hat(r_i)

    J_row1 = np.array([-self.n.T, -(ri_hat @ self.n).T, self.n.T, (rj_hat @ self.n).T])

    # [ [] [] [] [] ]
    # [              ]

    J_row1 = np.hstack(J_row1)

    # J_row2 = [ -t_1.T, - (ri_hat * t_1).T, t_1.T, (rj_hat * t_1).T]
    # J_row3 = [ -t_2.T, - (ri_hat * t_2).T, t_2.T, (rj_hat * t_2).T]

    # t1 = tangent plane vectors
    # t2 = tangent plane vectors
    J_row2 = np.array([-self.t1.T, -(ri_hat @ self.t1).T, self.t1.T, (rj_hat @ self.t1).T])
    J_row3 = np.array([-self.t2.T, -(ri_hat @ self.t2).T, self.t2.T, (rj_hat @ self.t2).T])
    J_row2 = np.hstack(J_row2)
    J_row3 = np.hstack(J_row3)
    J = np.array([J_row1.T, J_row2.T, J_row3.T])


    return J
  
  def get_b(self, J=None):
    # b = J(u.T + M^-1 * f)
    if J is None:
      J = self.compute_jacobian()
    u = self.compute_u()
    return J @ u
  
  
  
  def compute_delta_V(self, iters, mu):
    # delta V = M^-1 * J.T * lambda_k-1
    # where delta V is the change in velocity
    # lambda_k = lambda_k-1 - b_i' + J'_row_i  * delta V_k-1
    # delta lambda = -b_i' + J'_row_i  * delta V_k-1
    # T = M^-1 * J.T
    # delta V_k = delta V_k-1 + T_col_i * delta lambda_i

    # delta V_k = delta V_k-1 + T_col_i * ( -b_i' + J'_row_i  * delta V_k-1)

    J = self.compute_jacobian()

    M_inv = self.inv_generalized_mass_matrix()
    b = self.get_b(J=J)
    b_prime = b.copy()
    J_prime = J.copy()
    A = self.compute_effective_mass_matrix(J=J)

    # print(f'A: \n{A}')
    for i in range(J.shape[0]):
      J_prime[i,:] = J_prime[i,:] / A[i,i]
      b_prime[i] = b_prime[i] / A[i,i]

    delta_V = np.zeros((12,1))

    T = M_inv @ J.T


    T_cols = [T[:,i].reshape(-1,1) for i in range(3)]
    
    for k in range(iters):
      # print(f'delta_V: \n{delta_V}')
      for i in range(3):
        lambda_i_k = self._lambda[i].copy()
        lambda_i_k_1 = lambda_i_k - b_prime[i] - J_prime[i,:] @ delta_V

        # 0 <= lambda_i_1 <= inf
        # -mu * lambda_i_1 <= lambda_i_2 <= mu * lambda_i_1
        # -mu * lambda_i_1 <= lambda_i_3 <= mu * lambda_i_1

        if i == 0:
          lambda_i_k_1 = max(lambda_i_k_1, 0)
        else:
          abs_bound = mu * self._lambda[0]
          lambda_i_k_1 = min(max(lambda_i_k_1, -abs_bound), abs_bound)


        T_col_i = T_cols[i]
        self._lambda[i] = lambda_i_k_1
    
        delta_V += T_col_i * (lambda_i_k_1 - lambda_i_k)
    
    return delta_V
  
  def generalized_mass_matrix(self):
    M1 = self.body1.generalized_mass_matrix()
    M2 = self.body2.generalized_mass_matrix()
    M = np.block([[M1, np.zeros((6,6))], [np.zeros((6,6)), M2]])
    return M
  
  def inv_generalized_mass_matrix(self):
    M1_inv = self.body1.inv_generalized_mass_matrix()
    M2_inv = self.body2.inv_generalized_mass_matrix()
    M_inv = np.block([[M1_inv, np.zeros((6,6))], [np.zeros((6,6)), M2_inv]])
    return M_inv

  def compute_effective_mass_matrix(self, J=None):
    # TODO: implement this function

    # effective mass matrix
    # Mi = [m_i * I, 0;
    #       0, J_i^-1]

    # M = [M1, 0;
    #      0, M2]

    # M_eff = J * M^-1 * J.T
    # where J is the jacobian
    if J is None:
      J = self.compute_jacobian()
    M_inv = self.inv_generalized_mass_matrix()
    M_eff = J @ M_inv @ J.T
    return M_eff
  
  def compute_u(self):
    v1 = self.body1.v
    omega1 = self.body1.omega
    v2 = self.body2.v
    omega2 = self.body2.omega
    # make vertical vectors
    v1 = v1.reshape(-1,1)
    omega1 = omega1.reshape(-1,1)
    v2 = v2.reshape(-1,1)
    omega2 = omega2.reshape(-1,1)
    u = np.vstack([v1, omega1, v2, omega2])
    return u
  
