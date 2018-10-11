import numpy as np
class Ekl():
	def __init__(self, strain_key, st):
		if strain_key == 1:
			self.ekl = np.zeros((3,3))
			self.ekl[0,0] = 1 + st
			self.ekl[1,1] = 1
			self.ekl[2,2] = 1
		elif strain_key == 2:
			self.ekl = np.zeros((3,3))
			self.ekl[0,0] = 1
			self.ekl[1,1] = 1 + st
			self.ekl[2,2] = 1
		elif strain_key == 3:
			self.ekl = np.zeros((3,3))
			self.ekl[0,0] = 1
			self.ekl[1,1] = 1
			self.ekl[2,2] = 1 + st
		elif strain_key == 4:
			self.ekl = np.zeros((3,3))
			self.ekl[0,0] = 1
			self.ekl[1,1] = 1
			self.ekl[2,2] = 1
			self.ekl[1,2] = st
			self.ekl[2,1] = st
		elif strain_key == 5:
			self.ekl = np.zeros((3,3))
			self.ekl[0,0] = 1
			self.ekl[1,1] = 1
			self.ekl[2,2] = 1
			self.ekl[0,2] = st
			self.ekl[2,0] = st
		elif strain_key == 6:
			self.ekl = np.zeros((3,3))
			self.ekl[0,0] = 1
			self.ekl[1,1] = 1
			self.ekl[2,2] = 1
			self.ekl[0,1] = st
			self.ekl[1,0] = st

class apply_strain(object):
	def e(self, vectors, dist, i):
		e = Ekl(i, dist)
		return np.dot(vectors, e.ekl)

def define_cij(cclass):
	if cclass == '1' or cclass == '-1':
		print('You have a triclinic Bravais lattice, 6 strains are necessary for the calculation of the 21 elastic coefficients.')
		C11 = 'C11'
		C12 = 'C12'
		C13 = 'C13'
		C14 = 'C14'
		C15 = 'C15'
		C16 = 'C16'
		C22 = 'C22'
		C23 = 'C23'
		C24 = 'C24'
		C25 = 'C25'
		C26 = 'C26'
		C33 = 'C33'
		C34 = 'C34'
		C35 = 'C35'
		C36 = 'C36'
		C44 = 'C44'
		C45 = 'C45'
		C46 = 'C46'
		C55 = 'C55'
		C56 = 'C56'
		C66 = 'C66'
		print('Your Cij tensor has the form:')
		print_cij(C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56, C66)
		nstrain = 6
		return  nstrain
	elif cclass == '2' or cclass == 'm' or cclass == '2/m':
		print('You have a monoclinic Bravais lattice, 6 strains are necessary for the calculation of the 13 elastic coefficients.')
		C11 = 'C11'
		C12 = 'C12'
		C13 = 'C13'
		C14 = '0'
		C15 = 'C15'
		C16 = '0'
		C22 = 'C22'
		C23 = 'C23'
		C24 = '0'
		C25 = 'C25'
		C26 = '0'
		C33 = 'C33'
		C34 = '0'
		C35 = 'C35'
		C36 = '0'
		C44 = 'C44'
		C45 = '0'
		C46 = 'C46'
		C55 = 'C55'
		C56 = '0'
		C66 = 'C66'
		print_cij(C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56, C66)
		nstrain = 6
		return nstrain
	elif cclass == '222' or cclass == 'mm2' or cclass == '2mm' or cclass == 'mm' or cclass == 'mmm':
		print('You have an orthorrombic Bravais lattice, 6 strains are necessary for the calculation of the 9 elastic coefficients.')
		C11 = 'C11'
		C12 = 'C12'
		C13 = 'C13'
		C14 = '0'
		C15 = '0'
		C16 = '0'
		C22 = 'C22'
		C23 = 'C23'
		C24 = '0'
		C25 = '0'
		C26 = '0'
		C33 = 'C33'
		C34 = '0'
		C35 = '0'
		C36 = '0'
		C44 = 'C44'
		C45 = '0'
		C46 = '0'
		C55 = 'C55'
		C56 = '0'
		C66 = 'C66'
		print_cij(C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56, C66)
		nstrain = 6
		return nstrain
	elif cclass == '4' or cclass == '-4' or cclass == '4/m':
		print('You have a tetragonal Bravais lattice, 4 strains are necessary for the calculation of the 7 elastic coefficients.')
		C11 = 'C11'
		C12 = 'C12'
		C13 = 'C13'
		C14 = '0'
		C15 = '0'
		C16 = 'C16'
		C22 = 'C22'
		C23 = 'C13'
		C24 = '0'
		C25 = '0'
		C26 = '-C16'
		C33 = 'C33'
		C34 = '0'
		C35 = '0'
		C36 = '0'
		C44 = 'C44'
		C45 = '0'
		C46 = '0'
		C55 = 'C44'
		C56 = '0'
		C66 = 'C66'
		print_cij(C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56, C66)
		nstrain = 4
		return nstrain
	elif cclass == '4mm' or cclass == '422' or cclass == '42' or cclass == '-42m' or cclass == '4/mmm':
		print('You have a tetragonal Bravais lattice, 4 strains are necessary for the calculation of the 6 elastic coefficients.')
		C11 = 'C11'
		C12 = 'C12'
		C13 = 'C13'
		C14 = '0'
		C15 = '0'
		C16 = '0'
		C22 = 'C22'
		C23 = 'C13'
		C24 = '0'
		C25 = '0'
		C26 = '0'
		C33 = 'C33'
		C34 = '0'
		C35 = '0'
		C36 = '0'
		C44 = 'C44'
		C45 = '0'
		C46 = '0'
		C55 = 'C44'
		C56 = '0'
		C66 = 'C66'
		print_cij(C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56, C66)
		nstrain = 4
		return nstrain
	elif cclass == '23' or cclass == '-43m' or cclass == '432' or cclass =='43' or cclass == 'm3' or cclass == '2/m-3' or cclass == 'm3m' or cclass == 'm-3m':
		print('You have a cubic Bravais lattice, 2 strains are necessary for the calculation of the 3 elastic coefficients.')
		C11 = 'C11'
		C12 = 'C12'
		C13 = 'C12'
		C14 = '0'
		C15 = '0'
		C16 = '0'
		C22 = 'C11'
		C23 = 'C12'
		C24 = '0'
		C25 = '0'
		C26 = '0'
		C33 = 'C33'
		C34 = '0'
		C35 = '0'
		C36 = '0'
		C44 = 'C44'
		C45 = '0'
		C46 = '0'
		C55 = 'C44'
		C56 = '0'
		C66 = 'C44'
		print_cij(C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56, C66)
		nstrain = 2
		return nstrain
	elif cclass == '3' or cclass == '-3':
		print('You have a trigonal Bravais lattice, 3 strains are necessary for the calculation of the 7 elastic coefficients.')
		C11 = 'C11'
		C12 = 'C12'
		C13 = 'C13'
		C14 = 'C14'
		C15 = '-C15'
		C16 = '0'
		C22 = 'C11'
		C23 = 'C13'
		C24 = '-C14'
		C25 = 'C15'
		C26 = '0'
		C33 = 'C33'
		C34 = '0'
		C35 = '0'
		C36 = '0'
		C44 = 'C44'
		C45 = '0'
		C46 = '2*C15'
		C55 = 'C44'
		C56 = '2*C14'
		C66 = '0.5*(C11 - C12)'
		print_cij(C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56, C66)
		nstrain = 3
		return nstrain
	elif cclass == '32' or cclass == '3m' or cclass == '-3m' or cclass == '-32/m':
		print('You have a trigonal Bravais lattice, 3 strains are necessary for the calculation of the 6 elastic coefficients.')
		C11 = 'C11'
		C12 = 'C12'
		C13 = 'C13'
		C14 = '-C14'
		C15 = '0'
		C16 = '0'
		C22 = 'C11'
		C23 = 'C13'
		C24 = 'C14'
		C25 = '0'
		C26 = '0'
		C33 = 'C33'
		C34 = '0'
		C35 = '0'
		C36 = '0'
		C44 = 'C44'
		C45 = '0'
		C46 = '0'
		C55 = 'C44'
		C56 = '2*C14'
		C66 = '0.5*(C11 - C12)'
		print_cij(C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56, C66)
		nstrain = 3
		return nstrain
	elif cclass == '6' or cclass == '-6' or cclass == '622' or cclass == '62' or cclass == '6mm' or cclass == '6/m' or cclass == '6/mmm':
		print('You have a hexagonal Bravais lattice, 3 strains are necessary for the calculation of the 5 elastic coefficients.')
		C11 = 'C11'
		C12 = 'C12'
		C13 = 'C13'
		C14 = '0'
		C15 = '0'
		C16 = '0'
		C22 = 'C11'
		C23 = 'C13'
		C24 = '0'
		C25 = '0'
		C26 = '0'
		C33 = 'C33'
		C34 = '0'
		C35 = '0'
		C36 = '0'
		C44 = 'C44'
		C45 = '0'
		C46 = '0'
		C55 = 'C44'
		C56 = '0'
		C66 = '0.5*(C11 - C12)'
		print_cij(C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56, C66)
		nstrain = 3
		return nstrain
	elif cclass == '-6m2':
		print('You have a hexagonal Bravais lattice, 3 strains are necessary for the calculation of the 6 elastic coefficients.')
		C11 = 'C11'
		C12 = 'C12'
		C13 = 'C13'
		C14 = '0'
		C15 = '0'
		C16 = '0'
		C22 = 'C11'
		C23 = 'C23'
		C24 = '0'
		C25 = '0'
		C26 = '0'
		C33 = 'C33'
		C34 = '0'
		C35 = '0'
		C36 = '0'
		C44 = 'C33'
		C45 = '0'
		C46 = '0'
		C55 = 'C55'
		C56 = '0'
		C66 = '0.5*(C11 - C12)'
		print_cij(C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56, C66)
		nstrain = 3
		return nstrain
	else:
		print('Crystal class not recognized. Make sure you are using the symbols as listed in the International Tables for Crystallography.')
		quit()

def print_cij(C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56, C66):
	print('%4s %4s %4s %4s %4s %4s' %  (C11, C12, C13, C14, C15, C16))
	print('%4s %4s %4s %4s %4s %4s' %  (C12, C22, C23, C24, C25, C26))
	print('%4s %4s %4s %4s %4s %4s' %  (C13, C23, C33, C34, C35, C36))
	print('%4s %4s %4s %4s %4s %4s' %  (C14, C24, C34, C44, C45, C46))
	print('%4s %4s %4s %4s %4s %4s' %  (C15, C25, C35, C45, C55, C56))
	print('%4s %4s %4s %4s %4s %4s' %  (C16, C26, C36, C46, C56, C66))

