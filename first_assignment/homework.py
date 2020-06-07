epsilon = 10**-9

class Matrix: 
	def __init__(self,filePath = None):
		self.rows = 0
		self.collums = 0
		self.matrix = self.initizalizeMatrix(self.rows,self.collums)

		if(filePath != None):
			self.readFromFile(filePath)

	def __add__(self,matrix):
		newMatrix = self.initizalizeMatrix(self.rows,self.collums)
		newObject = Matrix()
		if(matrix.rows == self.rows and matrix.collums == self.collums):
			for i in range(self.rows):
				for j in range(self.collums):
					newMatrix[i][j] = self.matrix[i][j] + matrix.matrix[i][j]
			newObject.setMatrix(newMatrix)
			return newObject
		else:
			exit("Not a valid operation for adding!")

	def __sub__(self,matrix):
		newMatrix = self.initizalizeMatrix(self.rows,self.collums)
		newObject = Matrix()
		if(matrix.rows == self.rows and matrix.collums == self.collums):
			for i in range(self.rows):
				for j in range(self.collums):
					newMatrix[i][j] = self.matrix[i][j] - matrix.matrix[i][j]
			newObject.setMatrix(newMatrix)
			return newObject
		else:
			exit("Not a valid operation for substraction!")

	def __mul__(self,matrix):
		if( isinstance(matrix,int) or isinstance(matrix,float)):
			newMatrix = self.initizalizeMatrix(self.rows,self.collums)
			newObject = Matrix()
			for i in range(self.rows):
				for j in range(self.collums):
					newMatrix[i][j] = self.matrix[i][j]*matrix
			newObject.setMatrix(newMatrix)
			return newObject
		else:
			if(self.collums == matrix.rows):
				newMatrix = self.initizalizeMatrix(self.rows,matrix.collums)
				newObject = Matrix()
				for i in range(self.rows):
					for k in range(matrix.collums):
						newMatrix[i][k] = 0
						for j in range(self.collums):
							newMatrix[i][k] = newMatrix[i][k] + self.matrix[i][j]*matrix.matrix[j][k]
				newObject.setMatrix(newMatrix)
				return newObject
			else:
				exit("Not a valid operation for multiplication!")

	def __eq__(self,matrix):
		if(self.rows == matrix.rows and self.collums == matrix.collums):
			for i in range(self.rows):
				for j in range(self.collums):
					if(self.matrix[i][j] != matrix.matrix[i][j]):
						return False
			return True
		else:
			return False

	def __getitem__(self,coords):
		i,j = coords
		if(i <= self.rows and j <= self.collums):
			return self.matrix[i][j]
		else:
			exit("Not a valid operation!")

	def __setitem__(self,coords,data):
		i,j = coords
		if(i <= self.rows and j <= self.collums):
			self.matrix[i][j] = data
		else:
			exit("Not a valid operation!")

	def __str__(self):
		output = "Output:"
		output = output + "\n"
		output = output + "\n"
		for i in range(self.rows):
			for j in range(self.collums):
				if(j == 0):
					output = output + " "*10
				output = output + str(self.matrix[i][j]) + " "
			if(i < self.rows-1): 
				output = output + "\n"
		output = output + "\n"
		return output

	def traspone(self):
		trasponedMatrix = self.initizalizeMatrix(self.collums,self.rows)
		for i in range(self.rows):
			for j in range(self.collums):
				trasponedMatrix[j][i] = self.matrix[i][j]
		self.setMatrix(trasponedMatrix)

	def inverse(self):
		P_matrix = self.lup_decomposition()
		L_matrix = get_L_matrix(self)
		U_matrix = get_U_matrix(self)
		X_matrix = Matrix()
		X_matrix.setMatrix(self.initizalizeMatrix(self.rows,self.collums))
		for i in range(self.rows):
			e_single_vector = create_single_vector(self.rows,i)
			y_vector = supstitution_foward(L_matrix,P_matrix*e_single_vector)
			x_vector = supstitution_backward(U_matrix,y_vector)
			replaceCollumInMatrix(X_matrix,x_vector,i)
		self.setMatrix(X_matrix.matrix)
			
	def initizalizeMatrix(self,rows,collums):
		return [[0 for i in range(collums)] for j in range(rows)]

	def readFromFile(self,filePath):
		file = open(filePath,"r")
		lines = []
		for line in file:
			lines.append(line[:-1])
		file.close()

		if( not lines ):
			exit("Not a valid matrix file!")
		else:
			self.rows = len(lines)
			self.collums = len(lines[0].split(","))
			self.matrix = self.initizalizeMatrix(self.rows,self.collums)
		
			for i in range(len(lines)):
				values = lines[i].split(",")
				for j in range(len(values)):
					self.matrix[i][j] = float(values[j])

	def writeInFile(self,filePath):
		file = open(filePath,"w")
		for i in range(self.rows):
			line = ""
			for j in range(self.collums):
				line = line + str(self.matrix[i][j]) + str(",")
			line = line[:-1]
			line = line+"\n"
			file.write(line)
		file.close()

	def setMatrix(self,matrix):
		self.matrix = self.initizalizeMatrix(len(matrix),len(matrix[0]))
		self.rows = len(matrix)
		self.collums = len(matrix[0])
		for i in range(self.rows):
			for j in range(self.collums):
				self.matrix[i][j] = matrix[i][j]

	def roundMatrix(self):
		for i in range(self.rows):
			for j in range(self.collums):
				self.matrix[i][j] = round(self.matrix[i][j],2)

	def lu_decomposition(self):
		for i in range(self.rows-1):
			for j in range(i+1,self.rows):
				if(abs(self.matrix[i][i]) < epsilon):
					exit("Divison by zero!")
				self.matrix[j][i] = self.matrix[j][i]/self.matrix[i][i]
				for k in range(i+1,self.rows):
					self.matrix[j][k] = self.matrix[j][k] - self.matrix[j][i]*self.matrix[i][k]

	def lup_decomposition(self):
		P = [ 0 for i in range(self.rows) ]
		for i in range(self.rows):
			P[i] = i
		for i in range(self.rows-1):
			pivot = i
			for j in range(i+1,self.rows):
				if( abs(self.matrix[P[j]][i]) > abs(self.matrix[P[pivot]][i])):
					pivot = j

			temp = P[i]
			P[i] = P[pivot]
			P[pivot] = temp

			for j in range(i+1,self.rows):
				if( abs(self.matrix[P[i]][i]) < epsilon):
					exit("Divison by zero!")
				self.matrix[P[j]][i] = self.matrix[P[j]][i]/self.matrix[P[i]][i]
				for k in range(i+1,self.rows):
					self.matrix[P[j]][k] = self.matrix[P[j]][k] - self.matrix[P[j]][i]*self.matrix[P[i]][k]
		permutated_matrix,permutated_A = permutate_matrix(P,self)
		self.setMatrix(permutated_A)
		return permutated_matrix

def replaceCollumInMatrix(X_matrix,x_vector,position):
	for i in range(X_matrix.rows):
		X_matrix[i,position] = x_vector[i,0]
	
def supstitution_foward(L_matrix,b_vector):
	y_vector = b_vector
	for i in range(L_matrix.rows-1):
		for j in range(i+1,L_matrix.rows):
			y_vector[j,0] = y_vector[j,0]-L_matrix[j,i]*y_vector[i,0]
	return y_vector
	
def supstitution_backward(U_matrix,y_vector):
	x_vector = y_vector
	for i in range(U_matrix.rows-1,-1,-1):
		if(abs(U_matrix[i,i]) < epsilon):
			exit("Divison by zero in supstitution_backward!")
		x_vector[i,0] = x_vector[i,0]/U_matrix[i,i]
		for j in range(i):
			x_vector[j,0] = x_vector[j,0] - U_matrix[j,i] * x_vector[i,0]
	return x_vector

def permutate_matrix(P,A):
	permutation_matrix = Matrix()
	permutation_matrix.setMatrix(A.initizalizeMatrix(len(P),len(P)))
	for i in range(len(P)):
		permutation_matrix[i,P[i]] = 1
	A = permutation_matrix*A
	return permutation_matrix,A.matrix

def determinat_calculator(A_matrix):
	determinat = 1
	S = 0
	P_matrix = A_matrix.lup_decomposition()

	for i in range(P_matrix.rows):
		for j in range(P_matrix.collums):
			if(i == j and P_matrix[i,j] != 1):
				S = S+1
	S = S-1
	L = get_L_matrix(A_matrix)
	U = get_U_matrix(A_matrix)
	for i in range(A_matrix.rows):
		for j in range(A_matrix.collums):
			if(i == j):
				determinat = determinat*L[i,j]*U[i,j]
	p_coef = (-1)**S
	return determinat*p_coef

def get_L_matrix(A_matrix):
	L_matrix = Matrix()
	L_matrix.setMatrix(A_matrix.matrix)
	for i in range(A_matrix.rows):
		for j in range(A_matrix.collums):
			if( i == j ):
				L_matrix[i,j] = 1.0
			elif (i < j):
				L_matrix[i,j] = 0.0
	return L_matrix

def get_U_matrix(A_matrix):
	U_matrix = Matrix()
	U_matrix.setMatrix(A_matrix.matrix)
	for i in range(A_matrix.rows):
		for j in range(A_matrix.collums):
			if( i > j ):
				U_matrix[i,j] = 0.0	
	return U_matrix

def create_single_vector(n,position):
	if(position > n-1):
		exit("Not a valid operation!")

	single_vector = Matrix()
	single_vector.setMatrix(single_vector.initizalizeMatrix(n,1))
	single_vector[position,0] = 1
	return single_vector

def print_inline():
	print("")
	print("")





def main():
	'''
	Prvi zadatak
	'''
	print("PRVI ZADATAK")
	a = 0.5
	b = a/3
	b = a*3
	print(a == b)
	print_inline()



	'''
	Drugi zadatak
	'''

	print("DRUGI ZADATAK:")
	print_inline()
	A2 = Matrix("drugi_A.txt")
	b2 = Matrix("drugi_b.txt")
	#LU decomposition 
	#A2.lu_decomposition()
	#LUP decompostion
	P = A2.lup_decomposition()
	b2 = P*b2

	L = get_L_matrix(A2)
	U = get_U_matrix(A2)

	y2 = supstitution_foward(L,b2)
	x2 = supstitution_backward(U,y2)
	print(x2)
	print_inline()

	'''
	Treci zadatak
	'''
	print("TRECI ZADATAK:")
	print_inline()
	A3 = Matrix("treci_A.txt")
	b3 = Matrix("drugi_b.txt")
	#LU decomposition 
	#A3.lu_decomposition()
	#LUP decompostion 
	P = A3.lup_decomposition()
	b3 = P*b3
	print(A3)

	L = get_L_matrix(A3)
	U = get_U_matrix(A3)

	y3 = supstitution_foward(L,b3)
	#x3 = supstitution_backward(U,y3)
	#print(x3)
	print_inline()

	'''
	Cetvrti zadatak
	'''
	print("CETVRTI ZADATAK:")
	print_inline()
	A4 = Matrix("cetvrti_A.txt")
	b4 = Matrix("cetvrti_b.txt")
	#LU decomposition 
	#A4.lu_decomposition()
	#LUP decompostion
	P = A4.lup_decomposition()
	b4 = P*b4
	print(A4)


	L = get_L_matrix(A4)
	U = get_U_matrix(A4)

	y4 = supstitution_foward(L,b4)
	x4 = supstitution_backward(U,y4)
	print(x4)
	print_inline()

	'''
	Peti zadatak
	'''
	print("PETI ZADATAK:")
	print_inline()
	A5 = Matrix("peti_A.txt")
	b5 = Matrix("peti_b.txt")
	#LU decomposition
	#A5.lu_decomposition()
	#LUP decompostion
	P = A5.lup_decomposition()
	b5 = P*b5

	L = get_L_matrix(A5)
	U = get_U_matrix(A5)

	y5 = supstitution_foward(L,b5)
	x5 = supstitution_backward(U,y5)
	print(x5)
	print_inline()


	'''
	Šesti zadatak
	'''
	print("ŠESTI ZADATAK:")
	print_inline()
	A6 = Matrix("sesti_A.txt")
	b6 = Matrix("sesti_b.txt")
	#LU decomposition
	#A5.lu_decomposition()
	#LUP decompostion
	P = A6.lup_decomposition()
	b6 = P*b6

	L = get_L_matrix(A6)
	U = get_U_matrix(A6)

	y6 = supstitution_foward(L,b6)
	x6 = supstitution_backward(U,y6)
	print(x6)
	print_inline()

	'''
	Sedmi zadatak
	'''
	print("SEDMI ZADATAK:")
	print_inline()
	A7 = Matrix("sedmi_A.txt")
	
	#A7.inverse()
	print(A7)
	
	print_inline()


	'''
	Osmi zadatak
	'''
	print("OSMI ZADATAK:")
	print_inline()
	A8 = Matrix("osmi_A.txt")
	
	A8.inverse()
	print(A8)
	
	print_inline()

	'''
	Deveti zadatak
	'''
	print("DEVETI ZADATAK:")
	print_inline()
	A9 = Matrix("deveti_A.txt")
	
	print("Determinata iznosi: " + str(determinat_calculator(A9)))
	
	print_inline()


	'''
	Deseti zadatak
	'''
	print("DESETI ZADATAK:")
	print_inline()
	A10 = Matrix("deseti_A.txt")
	
	print("Determinata iznosi: " + str(determinat_calculator(A10)))
	
	print_inline()


main()