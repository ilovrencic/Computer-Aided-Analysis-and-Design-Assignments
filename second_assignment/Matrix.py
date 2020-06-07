class Matrix: 
	def __init__(self,rows = 0,collums = 0,filePath = None):
		self.rows = rows
		self.collums = collums
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
		output = ""
		for i in range(self.rows):
			for j in range(self.collums):
				if(j == 0):
					output = output + " "*10
				output = output + str(self.matrix[i][j]) + " "
			if(i < self.rows-1): 
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