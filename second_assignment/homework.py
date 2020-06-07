import math
import random
from Matrix import Matrix
k = 0.5*(math.sqrt(5)-1)
epsilon = 1
max_number_of_calls = 10**3
alpha = 1
beta = 0.5
theta = 2
sigma = 0.5
h = 1

class TwoVariableFunction:
	def __init__(self):
		self.a = 0
		self.b = 0
		self.c = 0
		self.d = 0
		self.n = 2
		self.number_of_calls = 0
		self.dict = dict()

	def checkCache(self,key,value):
		if(key in self.dict):
			return self.dict[key]
		else:
			self.dict[key] = value 
			return None

class RosenbrockFunction(TwoVariableFunction):
	def __init__(self):
		super().__init__()
		self.minimum = Matrix()
		self.minimum.setMatrix([[1,1]])
		self.start = Matrix()
		self.start.setMatrix([[-1.9,2]])

	def f2(self,x):
		self.number_of_calls += 1
		return 100*(x[0,1]-x[0,0])**2 + (1-x[0,0])**2

	def f1(self,x):
		x1 = self.a * x + self.b
		x2 = self.c * x + self.d

		cache = self.checkCache((x1,x2),100*(x2-x1)**2 + (1-x1)**2)
		if(cache == None):
			self.number_of_calls += 1
			return 100*(x2-x1)**2 + (1-x1)**2
		else:
			return cache

class SecondFunction(TwoVariableFunction):
	def __init__(self):
		super().__init__()
		self.minimum = Matrix()
		self.minimum.setMatrix([[4,2]])
		self.start = Matrix()
		self.start.setMatrix([[0.1,0.3]])


	def f2(self,x):
		self.number_of_calls += 1
		return (x[0,0]-4)**2+4*(x[0,1]-2)**2

	def f1(self,x):
		x1 = self.a * x + self.b
		x2 = self.c * x + self.d

		cache = self.checkCache((x1,x2),(x1-4)**2+4*(x2-2)**2)
		if(cache == None):
			self.number_of_calls += 1
			return (x1-4)**2+4*(x2-2)**2
		else:
			return cache

class JakobovicFunction(TwoVariableFunction):
	def __init__(self):
		super().__init__()
		self.minimum = Matrix()
		self.minimum.setMatrix([[0,0]])
		self.start = Matrix()
		self.start.setMatrix([[5.1,1.1]])

	def f2(self,x):
		self.number_of_calls += 1
		return abs((x[0,0]-x[0,1])*(x[0,0]+x[0,1]))+math.sqrt(x[0,0]**2 + x[0,1]**2)

	def f1(self,x):

		x1 = self.a * x + self.b
		x2 = self.c * x + self.d

		cache = self.checkCache((x1,x2),abs((x1-x2)*(x1+x2))+math.sqrt(x1**2 + x2**2))
		if(cache == None):
			self.number_of_calls += 1
			return abs((x1-x2)*(x1+x2))+math.sqrt(x1**2 + x2**2)
		else:
			return cache

class TestTwoVariableFunction(TwoVariableFunction):
	def __init__(self):
		super().__init__()
		self.minimum = Matrix()
		self.minimum.setMatrix([[0,0]])
		self.start = Matrix()
		self.start.setMatrix([[5,5]])

	def f2(self,x):
		self.number_of_calls += 1
		return (x[0,0]-2)**2 + (x[0,1]+1)**2

	def f1(self,x):
		self.number_of_calls += 1
		x1 = self.a * x + self.b
		x2 = self.c * x + self.d
		return (x1-2)**2 + (x2+1)**2

class MultiVariableFunction:
	def __init__(self,number_of_variables):
		self.n = number_of_variables
		self.number_of_calls = 0
		self.alphas = Matrix(1,number_of_variables)
		self.betas = Matrix(1,number_of_variables)
		self.start = Matrix(1,number_of_variables)
		self.dict = dict()

	def checkCache(key,value):
		if(key in self.dict):
			return self.dict[key]
		else:
			self.dict[key] = value 

class ThirdFunction(MultiVariableFunction):
	def __init__(self,number_of_variables):
		super().__init__(number_of_variables)

	def f1(self,x):
		self.number_of_calls += 1
		output = 0
		for i in range(1,self.n+1):
			xi = self.alphas[0,i-1]*x+self.betas[0,i-1]
			output = output + (xi-i)**2
		return output

	def f2(self,x):
		self.number_of_calls += 1
		output = 0
		for i in range(1,self.n+1):
			output = output + (x[0,i-1]-i)**2
		return output

class SchafferFunction(MultiVariableFunction):
	def __init__(self,number_of_variables):
		super().__init__(number_of_variables)

	def f1(self,x):
		self.number_of_calls += 1
		suma = 0
		for i in range(1,self.n+1):
			xi = self.alphas[0,i-1]*x+self.betas[0,i-1]
			suma = suma + xi**2
		return 0.5+(math.sin(math.sqrt(suma))**2 - 0.5)/(1+0.001*suma)**2

	def f2(self,x):
		self.number_of_calls += 1

		

		suma = 0
		for i in range(1,self.n+1):
			suma = suma + x[0,i-1]**2
		return 0.5+(math.sin(math.sqrt(suma))**2 - 0.5)/(1+0.001*suma)**2

class TestMultiVariableFunction(MultiVariableFunction):
	def __init__(self,number_of_variables):
		super().__init__(number_of_variables)
		self.start.setMatrix([[2,3,4]])

	def f1(self,x):
		self.number_of_calls += 1
		output = 0
		for i in range(1,self.n+1):
			xi = self.alphas[0,i-1]*x+self.betas[0,i-1]
			output = output + (xi)**2
		return output

	def f2(self,x):
		self.number_of_calls += 1
		output = 0
		for i in range(self.n):
			output = output + (x[0,i])**2
		return output
		
def unimodal_interval(function,x0,h):
	step = 1
	l = x0-h
	r = x0+h
	m = x0
	
	fm = function.f1(x0)
	fl = function.f1(l)
	fr = function.f1(r)

	if(fm < fl and fm < fr):
		return (l,r)

	elif(fm > fr):
		while(fm > fr):
			l = m
			m = r
			fm = fr
			r = x0+h*(2**step)
			fr = function.f1(r)
			step = step+1
		return (l,r)
	else:
		while(fm > fl):
			r = m
			m = l
			fm = fl
			l = x0-h*(2**step)
			fl = function.f1(l)
			step = step+1
		return (l,r)

def golden_ration(function,arguments):
	interval = ()
	if(isinstance(arguments,int)):
		interval = unimodal_interval(function,arguments,h)
	else:
		interval = arguments

	print(interval)
	a = interval[0]
	b = interval[1]
	c = b - k*(b-a)
	d = a + k*(b-a)

	fc = function.f1(c)
	fd = function.f1(d)
	while((b-a) > epsilon):
		if( fc < fd ):
			b = d
			d = c
			c = b - k*(b-a)
			fd = fc
			fc = function.f1(c)
		else:
			a = c
			c = d
			d = a + k*(b-a)
			fc = fd
			fd = function.f1(d)

	return(a,b)

def minimum(function,x,e):
	function.b = x[0,0]
	function.d = x[0,1]
	function.a = e[0,0]
	function.c = e[0,1]
	gamma = golden_ration(function,0)

	gamma_min = (gamma[0]+gamma[1])/2
	return x + e*gamma_min

def multi_minimum(function,x,e):
	for i in range(x.collums):
		function.alphas[0,i] = e[0,i]
		function.betas[0,i] = x[0,i]
	gamma = golden_ration(function,0)
	gamma_min = (gamma[0]+gamma[1])/2
	return x + e*gamma_min

def checkCondition(xs,x0,n):
	for i in range(n):
		if( abs(xs[0,i]-x0[0,i]) > epsilon ):
			return True
	return False 

def single_vector(i,n):
	e = Matrix(1,n)
	e[0,i] = 1
	return e

def discover(function,xP,delta):
	x = Matrix()
	x.setMatrix(xP.matrix)
	for i in range(x.collums):
		P = function.f2(x)
		x[0,i] = x[0,i] + delta
		N = function.f2(x)
		if(N > P):
			x[0,i] = x[0,i] - 2*delta
			N = function.f2(x)
			if(N > P):
				x[0,i] = x[0,i] + delta
	return x

def create_simplex_object(x0,delta = 1):
	simplex_object = list()
	for i in range(x0.collums):
		xi = Matrix()
		xi.setMatrix(x0.matrix)
		xi[0,i] = xi[0,i]+delta
		simplex_object.append(xi)
	xi = Matrix()
	xi.setMatrix(x0.matrix)
	simplex_object.append(xi)
	return simplex_object

def find_worst_point(function,simplex_object):
	index = 0
	worst = function.f2(simplex_object[index])
	for i in range(len(simplex_object)):
		if( function.f2(simplex_object[i]) > worst ):
			worst = function.f2(simplex_object[i])
			index = i
	return index

def find_best_point(function,simplex_object):
	index = 0
	best = function.f2(simplex_object[index])
	for i in range(len(simplex_object)):
		if( function.f2(simplex_object[i]) < best ):
			best = function.f2(simplex_object[i])
			index = i
	return index

def make_centroid(simplex_object,worst_index,n):
	xC = Matrix(1,n)
	for i in range(len(simplex_object)):
		if(i != worst_index):
			for j in range(n):
				xC[0,j] = xC[0,j]+simplex_object[i][0,j]
	for i in range(n):
		xC[0,i] = xC[0,i]*(1/(n))
	return xC

def checkSimplexCondition(function,xR,simplex_object,worst_index):
	for i in range(len(simplex_object)):
		if(i != worst_index):
			if(function.f2(xR) <= function.f2(simplex_object[i])):
				return False
	return True

def checkSimplexEnd(function,simplex_object,xC):
	value = 0
	for point in simplex_object:
		value = value + (function.f2(point) - function.f2(xC))**2
	value = value*(1/len(simplex_object))
	value = math.sqrt(value)
	if(value <= epsilon):
		return True
	else:
		return False

def translate_points(simplex_object,best_index):
	for point in simplex_object:
		direction = simplex_object[best_index]-point
		point = point + direction*sigma

def hooke_jeeves(function,trace = False,shouldPrint = True,delta = 1):
	xP = function.start
	xB = function.start
	while(delta > epsilon):
		xN = discover(function,xP,delta)
		if(trace):
			print(xB,xP,xN)
		if(function.f2(xN) < function.f2(xB)):
			xP = xN*2 - xB
			xB = xN
		else:
			delta = delta/2
			xP = xB
	if(shouldPrint):
		print("HOOKE JEEVES:")
	#xB.roundMatrix()
	return xB

def simplex(function,delta = 1,trace = False):
	simplex_object = create_simplex_object(function.start,delta)

	while(True):
		worst_index = find_worst_point(function,simplex_object)
		best_index = find_best_point(function,simplex_object)
		xC = make_centroid(simplex_object,worst_index,function.n)
		xR = xC*(1+alpha) - simplex_object[worst_index]*alpha


		if(function.f2(xR) < function.f2(simplex_object[best_index])):
			xE = xC*(1-theta) + xR*theta
			if(function.f2(xE) < function.f2(simplex_object[best_index])):
				simplex_object[worst_index] = xE
			else:
				simplex_object[worst_index] = xR
		else:
			if( checkSimplexCondition(function,xR,simplex_object,worst_index)):
				if(function.f2(xR) < function.f2(simplex_object[worst_index])):
					simplex_object[worst_index] = xR
				xK = xC*(1-beta) + simplex_object[worst_index]*beta
				
				if(function.f2(xK) < function.f2(simplex_object[worst_index])):
					simplex_object[worst_index] = xK
				else:
					translate_points(simplex_object,best_index)
			else:
				simplex_object[worst_index] = xR

		if(checkSimplexEnd(function,simplex_object,xC)):
			break

		
	print("SIMPLEX:")
	#xC.roundMatrix()
	return xC
	
def coordinate_descent(function):
	x = function.start
	xs = x
	n = function.n
	while( True ):
		xs = x
		for i in range(n):
			print(x)
			if(n == 2):
				x = minimum(function,x,single_vector(i,n))
			else:
				x = multi_minimum(function,x,single_vector(i,n))
		if(not checkCondition(xs,x,n)):
			break
	#x.roundMatrix()
	print("COORDINATE DESCENT:")
	return x

def print_inline():
	print("")
	print("")

def print_number_of_calls(function):
	print("        Number of calls: "+str(function.number_of_calls))
	function.number_of_calls = 0

def print_function_minimum(function,minimum,shouldPrint = True):
	function_minimum = function.f2(minimum)
	if(shouldPrint):
		print(minimum)
		print("Function minimum: "+str(function_minimum),end = "")
	else:
		if(abs(function_minimum) <= 10**-4):
			print(minimum)
			print("Function minimum: "+str(function_minimum),end = "")
	return function_minimum

def main():
	#PRVI ZADATAK
	print("PRVI ZADATAK")
	print("#################################")
	print_inline()
	first_function = TestTwoVariableFunction()
	print(coordinate_descent(first_function))
	print_inline()
	

main()