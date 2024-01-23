# TODO: YOUR NAME AND STUDENT NUMBER HERE

import numpy as np

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
		self.lamb = np.zeros(3)

	def compute_jacobian(self):
		# TODO: implement this function
		return

	def compute_inv_effective_mass(self):
		# TODO: implement this function
		return