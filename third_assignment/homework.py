import math
import random
from Matrix import Matrix
k = 0.5*(math.sqrt(5)-1)
epsilon = 10**-6
eta = 1
max_number_of_calls = 10**3
alpha = 1.3
beta = 0.5
theta = 2
sigma = 0.5
h = 1

class Function:
	def __init__(self):
		self.a = 0
		self.b = 0
		self.c = 0
		self.d = 0
		self.n = 2
		self.number_of_calls = 0

	def implicit_equations(self,x):
		return True

	def implicit_inqualities(self,x):
		return True

	def explicit_constraints(self,x):
		return True

class RosenbrockFunction(Function):
	def __init__(self):
		super().__init__()
		self.minimum = Matrix()
		self.minimum.setMatrix([[1,1]])
		self.start = Matrix()
		self.start.setMatrix([[-1.9,2]])
		self.xd = Matrix()
		self.xd.setMatrix([[-100,-100]])
		self.xg = Matrix()
		self.xg.setMatrix([[100,100]])

	def first_x1_derivation(self,x):
		return -400*x[0,0]*(x[0,1]-x[0,0]**2)-2*(1-x[0,0])

	def first_x2_derivation(self,x):
		return 200*(x[0,1]-x[0,0]**2)

	def second_x1_derivation(self,x):
		return 1200*x[0,0]**2-400*x[0,1]+2

	def second_x1x2_derivation(self,x):
		return -400*x[0,0]

	def second_x2_derivation(self,x):
		return 200

	def explicit_constarints(self,x):
		for i in range(x.collums):
			if(x[0,i] < self.xd[0,i] or x[0,i] > self.xg[0,i]):
				return False
		return True

	def implicit_inqualities(self,x):
		if(x[0,1] - x[0,0] < 0):
			return False
		elif(2-x[0,0] < 0):
			return False
		else:
			return True

	def f2(self,x):
		self.number_of_calls += 1
		return 100*(x[0,1]-x[0,0]**2)**2+(1-x[0,0])**2

	def f1(self,x):
		x1 = self.a * x + self.b
		x2 = self.c * x + self.d

		self.number_of_calls += 1
		return 100*(x2-x1**2)**2+(1-x1)**2

class SecondFunction(Function):
	def __init__(self):
		super().__init__()
		self.minimum = Matrix()
		self.minimum.setMatrix([[4,2]])
		self.start = Matrix()
		self.start.setMatrix([[0.1,0.3]])
		self.xd = Matrix()
		self.xd.setMatrix([[-100,-100]])
		self.xg = Matrix()
		self.xg.setMatrix([[100,100]])

	def first_x1_derivation(self,x):
		return 2*(x[0,0]-4)

	def first_x2_derivation(self,x):
		return 8*(x[0,1]-2)

	def second_x1_derivation(self,x):
		return 2

	def second_x1x2_derivation(self,x):
		return 0

	def second_x2_derivation(self,x):
		return 8

	def explicit_constarints(self,x):
		for i in range(x.collums):
			if(x[0,i] < self.xd[0,i] or x[0,i] > self.xg[0,i]):
				return False
		return True

	def implicit_inqualities(self,x):
		if(x[0,1] - x[0,0] < 0):
			return False
		elif(2-x[0,0] < 0):
			return False
		else:
			return True

	def f2(self,x):
		self.number_of_calls += 1
		return (x[0,0]-4)**2 + 4*(x[0,1]-2)**2

	def f1(self,x):
		x1 = self.a * x + self.b
		x2 = self.c * x + self.d

		self.number_of_calls += 1
		return (x1-4)**2 + 4*(x2-2)**2

class ThridFunction(Function):
	def __init__(self):
		super().__init__()
		self.minimum = Matrix()
		self.minimum.setMatrix([[2,-3]])
		self.start = Matrix()
		self.start.setMatrix([[0,0]])
		self.xd = Matrix()
		self.xd.setMatrix([[1,1]])
		self.xg = Matrix()
		self.xg.setMatrix([[100,100]])

	def first_x1_derivation(self,x):
		return 2*(x[0,0]-2)

	def first_x2_derivation(self,x):
		return 2*(x[0,1]+3)

	def second_x1_derivation(self,x):
		return 2

	def second_x1x2_derivation(self,x):
		return 0

	def second_x2_derivation(self,x):
		return 2

	def f2(self,x):
		self.number_of_calls += 1
		return (x[0,0]-2)**2 + (x[0,1]+3)**2

	def f1(self,x):
		x1 = self.a * x + self.b
		x2 = self.c * x + self.d

		self.number_of_calls += 1
		return (x1-2)**2 + (x2+3)**2

class FourthFunction(Function):
	def __init__(self):
		super().__init__()
		self.minimum = Matrix()
		self.minimum.setMatrix([[3,0]])
		self.start = Matrix()
		self.start.setMatrix([[0,0]])
		self.xd = Matrix()
		self.xd.setMatrix([[1,1]])
		self.xg = Matrix()
		self.xg.setMatrix([[100,100]])

	def first_x1_derivation(self,x):
		return 2*(x[0,0]-3)

	def first_x2_derivation(self,x):
		return 2*x[0,1]

	def second_x1_derivation(self,x):
		return 2

	def second_x1x2_derivation(self,x):
		return 0

	def second_x2_derivation(self,x):
		return 2

	def f2(self,x):
		self.number_of_calls += 1
		return (x[0,0]-3)**2 + x[0,1]**2

	def f1(self,x):
		x1 = self.a * x + self.b
		x2 = self.c * x + self.d

		self.number_of_calls += 1
		return (x1-3)**2 + x2**2

class RosenbrockConstraintFunction(Function):
	def __init__(self):
		super().__init__()
		self.minimum = Matrix()
		self.minimum.setMatrix([[1,1]])
		self.start = Matrix()
		self.start.setMatrix([[-1.9,2]])

	def f2(self,x,t):
		self.number_of_calls += 1

		if(x[0,1]-x[0,0] <= 0 and 2-x[0,0] <= 0):
			return 100*(x[0,1]-x[0,0]**2)**2+(1-x[0,0])**2 - 10**100-10**100
		elif(2-x[0,0] <= 0):
			return 100*(x[0,1]-x[0,0]**2)**2+(1-x[0,0])**2 - (1/t)*math.log(x[0,1]-x[0,0])-10**100
		elif(x[0,1]-x[0,0] <= 0):
			return 100*(x[0,1]-x[0,0]**2)**2+(1-x[0,0])**2 - 10**100-(1/t)*math.log(2-x[0,0])

		return 100*(x[0,1]-x[0,0]**2)**2+(1-x[0,0])**2 - (1/t)*math.log(x[0,1]-x[0,0])-(1/t)*math.log(2-x[0,0])

class SecondConstraintFunction(Function):
	def __init__(self):
		super().__init__()
		self.minimum = Matrix()
		self.minimum.setMatrix([[4,2]])
		self.start = Matrix()
		self.start.setMatrix([[0.1,0.3]])

	def f2(self,x,t):
		self.number_of_calls += 1

		if(x[0,1]-x[0,0] <= 0 and 2-x[0,0] <= 0):
			return (x[0,0]-4)**2 + 4*(x[0,1]-2)**2 - (10**100)-(10**100)
		elif(2-x[0,0] <= 0):
			return (x[0,0]-4)**2 + 4*(x[0,1]-2)**2 - (1/t)*math.log(x[0,1]-x[0,0])-(10**100)
		elif(x[0,1]-x[0,0] <= 0):
			return (x[0,0]-4)**2 + 4*(x[0,1]-2)**2 - (10**100)-(1/t)*math.log(2-x[0,0])

		return (x[0,0]-4)**2 + 4*(x[0,1]-2)**2 - (1/t)*math.log(x[0,1]-x[0,0])-(1/t)*math.log(2-x[0,0])

class FourthConstraintFunction(Function):
	def __init__(self):
		super().__init__()
		self.minimum = Matrix()
		self.minimum.setMatrix([[3,0]])
		self.start = Matrix()
		self.start.setMatrix([[0,0]])

	def f2(self,x,t):
		self.number_of_calls += 1
		if(3-x[0,0]-x[0,1] <= 0 and 3+1.5*x[0,0]-x[0,1] <= 0):
			return (x[0,0]-3)**2 + x[0,1]**2 - 10*100 - 10*100 + t*(x[0,1]-1)**2
		elif(3-x[0,0]-x[0,1] <= 0):
			return (x[0,0]-3)**2 + x[0,1]**2 - 10*100 - (1/t)*math.log(3+1.5*x[0,0]-x[0,1]) + t*(x[0,1]-1)**2
		elif(3+1.5*x[0,0]-x[0,1] <= 0):
			return (x[0,0]-3)**2 + x[0,1]**2 - (1/t)*math.log(3-x[0,0]-x[0,1]) - 10*100 + t*(x[0,1]-1)**2
		return (x[0,0]-3)**2 + x[0,1]**2 - (1/t)*math.log(3-x[0,0]-x[0,1]) - (1/t)*math.log(3+1.5*x[0,0]-x[0,1]) + t*(x[0,1]-1)**2

class InsidePointFinder(Function):
	def __init__(self):
		super().__init__()
		self.start = Matrix()
		self.start.setMatrix([[0,0]])

	def f2(self,x,t):
		return -1*(2*(3-x[0,0]-x[0,1]) + 4*(3+1.5*x[0,0]-x[0,1]))

def mixed_way_transformation(function):
	x_curr = Matrix()
	x_curr.setMatrix(function.start.matrix)
	t = 1
	while(True):
		x_next = hooke_jeeves(function,t,trace = False)
		t *= 10
		function.start.setMatrix(x_next.matrix)

		if(checkEnding(x_next,x_curr)):
			print("Mixed way transformation:")
			print("")
			print("Minimum:",x_next)
			return x_next

		x_curr.setMatrix(x_next.matrix)

def checkEnding(x1,x2):
	distance = 0
	for i in range(x1.collums):
		distance += (x1[0,i]-x2[0,i])**2

	distance = math.sqrt(distance)
	if(distance < epsilon):
		return True
	return False

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

def gradient_descent(function, useGoldenRation = True,printIteration = True):
	x = function.start
	number_of_iterations = 0
	old_x = Matrix()
	old_x.setMatrix(x.matrix)

	count_divergation = 0
	while( True ):
		number_of_iterations += 1
		if(printIteration):
			print(x)

		v = calculate_gradient(function,x)
		if(v.norm(changeMatrix = False) < epsilon):
			print("Gradient descent:")
			print("")
			x.roundMatrix()
			print(x)
			print("Number of iterations: ",number_of_iterations)
			print("")
			break

		if(useGoldenRation):
			x = minimum_on_line(function,x,v)
		else:
			x = x - v*eta

		if(function.f2(old_x) <= function.f2(x)):
			count_divergation += 1
			old_x.setMatrix(x.matrix)

		if(count_divergation == 100):
			print("Gradient is diverging!")
			break

def newton_raphson(function,useGoldenRation = True,printIteration = True):
	x = function.start
	number_of_iterations = 0

	old_x = Matrix()
	old_x.setMatrix(x.matrix)

	count_divergation = 0
	while(True):
		number_of_iterations += 1
		if(printIteration):
			print(x)

		v = calculate_direction_vector(function,x)
		if(v.norm(changeMatrix = False) < epsilon):
			print("Newton-raphson:")
			print("")
			x.roundMatrix()
			print(x)
			print("Number of iterations: ",number_of_iterations)
			print("")
			break


		if(useGoldenRation):
			x = minimum_on_line(function,x,v)
		else:
			x = x - v*eta

		if(function.f2(old_x) <= function.f2(x)):
			count_divergation += 1
			old_x.setMatrix(x.matrix)

		if(count_divergation == 100):
			print("Gradient is diverging!")
			break

def box(function,x0):
	if(not checkConstraints(function,x0)):
		print("Not a valid point!")
		return None
	xC = Matrix()
	xC.setMatrix(x0.matrix)

	simplex_object = create_simplex_object(function,xC)
	old_xC = Matrix()
	old_xC.setMatrix(xC.matrix)
	count_divergation = 0
	while(True):
		worst_index = find_worst_point(function,simplex_object)
		second_worst_index = find_worst_point(function,simplex_object,second_worst = True)
		xC = make_centroid(simplex_object,worst_index,x0.collums)
		xR = xC*(1+alpha) - simplex_object[worst_index]*alpha

		for i in range(xR.collums):
			if(xR[0,i] < function.xd[0,i]):
				xR[0,i] = function.xd[0,i]
			elif(xR[0,i] > function.xg[0,i]):
				xR[0,i] = function.xg[0,i]

		xR = fitPointToConstraints(function,xR,xC)
		xR = checkFunctionCostConstraint(function,simplex_object[second_worst_index],xC,xR)

		simplex_object[worst_index] = xR

		if(checkSimplexEnd(function,simplex_object,xC)):
			print("Box:")
			print("")
			xC.roundMatrix()
			print("Minimum:",xC)
			break

		if(function.f2(old_xC) <= function.f2(xC)):
			count_divergation += 1

		if(count_divergation == 1000):
			print("Box:")
			print("")
			print("Can't find minimum!")
			print("Last known xC:",xC)
			break

		old_xC.setMatrix(xC.matrix)

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

def checkFunctionCostConstraint(function,xH,xC,xR):
	if(function.f2(xH) > function.f2(xR)):
		return xR
	else:
		xR = (xR + xC)*0.5
		return xR

def make_centroid(simplex_object,worst_index,n):
	xC = Matrix(1,n)
	for i in range(len(simplex_object)):
		if(i != worst_index):
			for j in range(n):
				xC[0,j] = xC[0,j]+simplex_object[i][0,j]
	for i in range(n):
		xC[0,i] = xC[0,i]*(1/(len(simplex_object)-1))
	return xC

def find_worst_point(function,simplex_object,second_worst = False):
	index = 0
	worst = function.f2(simplex_object[index])
	for i in range(len(simplex_object)):
		if( function.f2(simplex_object[i]) > worst ):
			worst = function.f2(simplex_object[i])
			index = i

	second_index = index
	if(index == 0):
		second_index += 1
	else:
		second_index -= 1

	worst = function.f2(simplex_object[second_index])
	for i in range(len(simplex_object)):
		if(i != index):
			if( function.f2(simplex_object[i]) > worst):
				worst = function.f2(simplex_object[i])
				second_index = i

	if(second_worst):
		return second_index
	return index

def find_best_point(function,simplex_object):
	index = 0
	best = function.f2(simplex_object[index])
	for i in range(len(simplex_object)):
		if( function.f2(simplex_object[i]) < best ):
			best = function.f2(simplex_object[i])
			index = i
	return index

def create_simplex_object(function,xC,delta = 1):
	simplex_object = list()
	for i in range(2**xC.collums):
		xi = Matrix()
		xi.setMatrix(xC.matrix)
		for j in range(xC.collums):
			random_number = random.random()
			xi[0,j] = function.xd[0,j] + random_number*(function.xg[0,j]-function.xd[0,j])

		xi = fitPointToConstraints(function,xi,xC)
		simplex_object.append(xi)
		xC = calculate_centroid_build(simplex_object,xC.collums)
	return simplex_object

def calculate_centroid_build(simplex_object,n):
	xC = Matrix(1,n)
	for i in range(len(simplex_object)):
		for j in range(n):
			xC[0,j] = xC[0,j]+simplex_object[i][0,j]
	for i in range(n):
		xC[0,i] = xC[0,i]*(1/(len(simplex_object)))
	return xC

def fitPointToConstraints(function,xi,xC):
	if(checkConstraints(function,xi)):
		return xi
	x_new = Matrix()
	x_new.setMatrix(xi.matrix)
	while(not checkConstraints(function,x_new)):
		x_new = (x_new + xC)*0.5
	return x_new

def hooke_jeeves(function,t,trace = False,shouldPrint = True,delta = 1,):
	xP = function.start
	xB = function.start
	while(delta > epsilon):
		xN = discover(function,xP,delta,t)
		if(trace):
			print(xB,xP,xN)
		if(abs(function.f2(xN,t)) < abs(function.f2(xB,t))):
			xP = xN*2 - xB
			xB = xN
		else:
			delta = delta/2
			xP = xB
	#xB.roundMatrix()
	return xB

def discover(function,xP,delta,t):
	x = Matrix()
	x.setMatrix(xP.matrix)
	for i in range(x.collums):
		P = function.f2(x,t)
		x[0,i] = x[0,i] + delta
		N = function.f2(x,t)
		if(abs(N) > abs(P)):
			x[0,i] = x[0,i] - 2*delta
			N = function.f2(x,t)
			if(abs(N) > abs(P)):
				x[0,i] = x[0,i] + delta
	return x

def checkConstraints(function,x):
	return function.implicit_inqualities(x) and function.explicit_constarints(x) and function.implicit_equations(x)

def calculate_gradient(function,x):
	v = Matrix()
	v.setMatrix([[function.first_x1_derivation(x),function.first_x2_derivation(x)]])
	return v

def calculate_direction_vector(function,x):
	hess = Matrix()
	hess.setMatrix([[function.second_x1_derivation(x),function.second_x1x2_derivation(x)],[function.second_x1x2_derivation(x),function.second_x2_derivation(x)]])
	#print(hess)
	hess.inverse()


	gradient = calculate_gradient(function,x)
	gradient.traspone()

	v = hess*gradient
	v.traspone()
	return v

def minimum_on_line(function,x,v):
	function.b = x[0,0]
	function.d = x[0,1]
	function.a = v[0,0]
	function.c = v[0,1]
	gamma = golden_ration(function,0)

	gamma_min = (gamma[0]+gamma[1])/2
	return x + v*gamma_min

def main():
	#PRVI ZADATAK
	print("PRVI ZADATAK:")
	first_function = ThridFunction()
	gradient_descent(first_function,useGoldenRation = False,printIteration = False)
	print("")
	print("")

	#DRUGI ZADATAK
	print("DRUGI ZADATAK:")
	second_function = RosenbrockFunction()
	third_function = SecondFunction()
	print("-------Rosenbrock Function--------")
	gradient_descent(second_function,printIteration = False)
	newton_raphson(second_function,useGoldenRation = False,printIteration = True)
	print("-------Second function---------")
	gradient_descent(third_function,printIteration = False)
	newton_raphson(third_function,printIteration = False)
	print("")
	print("")

	#TRECI ZADATAK
	print("TRECI ZADATAK:")
	second_function = RosenbrockFunction()
	third_function = SecondFunction()
	box(second_function,second_function.start)
	box(third_function,third_function.start)

	#CETVRTI ZADATAK
	print("CETVRTI ZADATAK:")
	print("------------Rosenbrock function start (-1.9,2)--------------")
	rosenbrock_constrained = RosenbrockConstraintFunction()
	mixed_way_transformation(rosenbrock_constrained)
	print("")

	print("------------Second function start (0.1,0.3)--------------")
	second_function_constrained = SecondConstraintFunction()
	mixed_way_transformation(second_function_constrained)
	print("")

	print("------------Rosenbrock function start (1.5,10)--------------")
	rosenbrock_constrained = RosenbrockConstraintFunction()
	rosenbrock_constrained.start.setMatrix([[1.5,10]])
	mixed_way_transformation(rosenbrock_constrained)
	print("")


	#PETI ZADATAK
	print("PETI ZADATAK:")
	print("------------Fourth function start (5,5)--------------")
	fourth_function_constrained = FourthConstraintFunction()
	fourth_function_constrained.start.setMatrix([[5,5]])
	mixed_way_transformation(fourth_function_constrained)
	print("")

	find_inside_point = mixed_way_transformation(InsidePointFinder())
	print("------------Fourth function start finded by function minimization--------------")
	fourth_function_constrained = FourthConstraintFunction()
	fourth_function_constrained.start.setMatrix(find_inside_point.matrix)
	mixed_way_transformation(fourth_function_constrained)
	print("")


main()
