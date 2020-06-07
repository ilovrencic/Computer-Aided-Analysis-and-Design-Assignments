from Matrix import Matrix
import numpy as np
import math
import random

class TwoVariableFunction:
	def __init__(self):
		self.a = 0
		self.b = 0
		self.c = 0
		self.d = 0
		self.n = 2
		self.number_of_calls = 0
		self.dict = dict()

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

class MultiVariableFunction:
	def __init__(self,number_of_variables):
		self.n = number_of_variables
		self.number_of_calls = 0
		self.alphas = Matrix(1,number_of_variables)
		self.betas = Matrix(1,number_of_variables)
		self.start = Matrix(1,number_of_variables)

class ThirdFunction(MultiVariableFunction):
	def __init__(self,number_of_variables):
		super().__init__(number_of_variables)

	def f2(self,x):
		self.number_of_calls += 1
		output = 0
		for i in range(1,self.n+1):
			output = output + (x[0,i-1]-i)**2
		return output

class SchafferFunction(MultiVariableFunction):
	def __init__(self,number_of_variables):
		super().__init__(number_of_variables)

	def f2(self,x):
		self.number_of_calls += 1
		suma = 0
		for i in range(1,self.n+1):
			suma = suma + x[0,i-1]**2
		return 0.5+(math.sin(math.sqrt(suma))**2 - 0.5)/(1+0.001*suma)**2

class AlmostSchafferFunction(MultiVariableFunction):
    def __init__(self,number_of_variables):
        super().__init__(number_of_variables)

    def f2(self,x):
        self.number_of_calls += 1
        suma = 0
        for i in range(1,self.n+1):
            suma = suma + x[0,i-1]**2
        return ((suma)**0.25) * (1 + math.sin(50*((suma)**0.1))**2)


class GeneticAlgorithm:
    def __init__ (self,
                    function,
                    population_size = 15,
                    mutation_probability = 0.2,
                    crossing_probabilty = 0.5,
                    number_of_evaluations = 5000,
                    solution_type = "binary",
                    crossing_type = "singular",
                    turnir_size = 3,
                    elitisim = 0):
        self.population_size = population_size
        self.solution_type = solution_type
        self.m_prob = mutation_probability
        self.c_prob = crossing_probabilty
        self.n_evaluations = number_of_evaluations
        self.keep = elitisim
        self.population = []
        self.dimension = 0
        self.function = function
        self.crossing_type = crossing_type
        self.k = turnir_size
        self.n = 0
        self.dd = 0
        self.gg = 0
        self.precision = 0
        self.number_of_hits = 0

        #Iteration counter
        self.i = 0

    def calculate_fitness(self,chromoson):
        input = Matrix()
        if(self.solution_type == "binary"):
            matrix = []
            for feature in chromoson:
                b = int(feature,2)
                x = self.dd + (b/(2**self.n))*(self.gg-self.dd)
                matrix.append(x)
            input.setMatrix([matrix])
        else:
            input.setMatrix([chromoson])
        return abs(1/(self.function.f2(input)))

    def crossing_binary(self,picked_chromosoms):
        if(self.crossing_type == "uniform"):
            chromosonA = picked_chromosoms[0][0]
            chromosonB = picked_chromosoms[1][0]
            chromosonR = ""
            for bit in range(self.n*self.dimension):
                chromosonR += str(random.randint(0,1))

            A = ""
            B = ""
            for i in range(self.dimension):
                A += chromosonA[i]
                B += chromosonB[i]

            chromosonA = A
            chromosonB = B

            childX = ""
            child = []
            for i in range(len(chromosonR)):
                suma = int(chromosonA[i])*int(chromosonB[i])+int(chromosonR[i])*int(chromosonA[i])+int(chromosonR[i])*int(chromosonB[i])

                if(suma > 0):
                    childX += "1"
                else:
                    childX += "0"

                if((i+1)%self.n == 0):
                    child.append(childX)
                    childX = ""

            return child
        else:
            chromosonA = picked_chromosoms[0][0]
            chromosonB = picked_chromosoms[1][0]

            A = ""
            B = ""
            for i in range(self.dimension):
                A += chromosonA[i]
                B += chromosonB[i]

            chromosonA = A
            chromosonB = B

            point_of_crossover = round(self.n*self.c_prob)
            A = chromosonA[:point_of_crossover]+chromosonB[point_of_crossover:]
            B = chromosonB[:point_of_crossover]+chromosonA[point_of_crossover:]
            chromosonA = A
            chromosonB = B

            child = []
            childX = ""
            for i in range(len(chromosonA)):
                childX += chromosonA[i]

                if((i+1)%self.n == 0):
                    child.append(childX)
                    childX = ""

            return child

    def mutation_binary(self,chromoson):
        for dimension in range(self.dimension):
            for bit in range(self.n):
                p = random.random()
                mutated_gene = ""
                if(p < self.m_prob):
                    #print(chromoson)
                    if(chromoson[dimension][bit] == "0"):
                        mutated_gene = chromoson[dimension][:bit]+"1"+chromoson[dimension][bit+1:]
                    else:
                        mutated_gene = chromoson[dimension][:bit]+"0"+chromoson[dimension][bit+1:]
                    chromoson[dimension] = mutated_gene
        return chromoson

    def get_random_index(self,n):
        indexs = []
        while(len(indexs) != self.k):
            number = np.random.randint(0,n-1)
            if(number not in indexs):
                indexs.append(number)
        return indexs

    def remove_worst_chromosom(self,picked_chromosoms):
        repicked_chromosoms = []
        worst_fitness = picked_chromosoms[0][1]
        index = 0
        worst_chromosom = []
        for i in range(len(picked_chromosoms)):
            if(picked_chromosoms[i][1] < worst_fitness):
                index = i
                worst_fitness = picked_chromosoms[i][1]

        for i in range(len(picked_chromosoms)):
            if(index != i):
                repicked_chromosoms.append(picked_chromosoms[i])
            else:
                worst_chromosom = picked_chromosoms[i]
        self.population.remove(worst_chromosom)
        return repicked_chromosoms

    def return_best_chromosom(self):
        best_fitness = self.population[0][1]
        index = 0
        for i in range(self.population_size):
            if(self.population[i][1] > best_fitness):
                index = i
                best_fitness = self.population[i][1]

        return self.population[index]

    def init(self,dd,gg,precision,dimension):
        self.dimension = dimension
        self.dd = dd
        self.gg = gg
        self.precision = precision
        if(self.solution_type == "binary"):
            self.init_population_binary()
        else:
            self.init_population_real()

    def init_population_binary(self):
        n = math.ceil(math.log((self.gg-self.dd)*(10**self.precision)+1)/(math.log(2)))
        self.n = n
        for i in range(self.population_size):
            chromoson = []
            for d in range(self.dimension):
                feature = ""
                for bit in range(n):
                    feature += str(random.randint(0,1))
                chromoson.append(feature)
            fitness = self.calculate_fitness(chromoson)
            self.population.append((chromoson,fitness))
        #print(self.population)

    def init_population_real(self):
        for i in range(self.population_size):
            chromoson = []
            for d in range(self.dimension):
                number = random.uniform(self.dd,self.gg)
                chromoson.append(number)
            fitness = self.calculate_fitness(chromoson)
            self.population.append((chromoson,fitness))

    def binary_cromosom_value(self,chromoson):
        input = Matrix()
        matrix = []
        for feature in chromoson[0]:
            b = int(feature,2)
            x = self.dd + (b/(2**self.n))*(self.gg-self.dd)
            matrix.append(x)
        input.setMatrix([matrix])
        return input

    def print_info_about_current_population_binary(self,iteration):
        best_chromosom = self.return_best_chromosom()
        chromoson_value = self.binary_cromosom_value(best_chromosom)
        if(iteration == self.n_evaluations-1):
            print("Iteracija: "+str(iteration)+" ",end="")
            print(chromoson_value,end="")
            print(" Function value: "+str(1/best_chromosom[1]))

        #print(" ",end="")
        #print(best_chromosom[0])

    def crossing_real(self,chromosoms):
        if(self.crossing_type == "aritmetic"):
            X1 = Matrix()
            X2 = Matrix()
            X1.setMatrix([chromosoms[0][0]])
            X2.setMatrix([chromosoms[1][0]])
            alpha = random.random()
            C = X1*alpha + X2*(1-alpha)
            return C.matrix[0]
        else:
            X1 = Matrix()
            X2 = Matrix()
            if(chromosoms[0][1] > chromosoms[1][1]):
                X1.setMatrix([chromosoms[1][0]])
                X2.setMatrix([chromosoms[0][0]])
            else:
                X1.setMatrix([chromosoms[0][0]])
                X2.setMatrix([chromosoms[1][0]])
            alpha = random.random()
            C = (X2-X1)*alpha + X2
            return C.matrix[0]

    def mutation_real(self,chromoson):
        for i in range(len(chromoson)):
            p = random.random()
            if(p < self.m_prob):
                chromoson[i] = random.uniform(self.dd,self.gg)
        return chromoson


    def print_info_about_current_population_real(self,iteration):
        best_chromosom = self.return_best_chromosom()
        chromoson_value = Matrix()
        chromoson_value.setMatrix([best_chromosom[0]])
        if(iteration == self.n_evaluations-1):
            print("Iteracija: "+str(iteration)+" ",end="")
            print(chromoson_value,end="")
            print(" Function value: "+str(1/best_chromosom[1]))

    def start(self,shouldPrint = True):
        while(self.i < self.n_evaluations):
            #Choosing of three random chromosons and removing the worst one
            random_indexs = self.get_random_index(self.population_size)
            picked_chromosoms = []
            for index in random_indexs:
                picked_chromosoms.append(self.population[index])

            repicked_chromosoms = self.remove_worst_chromosom(picked_chromosoms)

            #Crossing of two remaining chromosons
            if(self.solution_type == "binary"):
                new_chromosom = self.crossing_binary(repicked_chromosoms)
                new_chromosom = self.mutation_binary(new_chromosom)
                fitness = self.calculate_fitness(new_chromosom)
                self.population.append((new_chromosom,fitness))
                if(shouldPrint):
                    self.print_info_about_current_population_binary(self.i)
            else:
                new_chromosom = self.crossing_real(repicked_chromosoms)
                new_chromosom = self.mutation_real(new_chromosom)
                fitness = self.calculate_fitness(new_chromosom)
                self.population.append((new_chromosom,fitness))
                if(shouldPrint):
                    self.print_info_about_current_population_real(self.i)
            self.i += 1
        return self.return_best_chromosom()

def main():
    #PRVI ZADATAK
    print("----------------PRVI ZADATAK------------------")
    rosenbrock = RosenbrockFunction()

    print("Rosenbrock binary:")
    GA = GeneticAlgorithm(function = rosenbrock)
    GA.init(dd = -50,gg = 150,precision = 5,dimension = 2)
    GA.start()

    print("Rosenbrock real:")
    GA = GeneticAlgorithm(solution_type = "real",function = rosenbrock)
    GA.init(dd = -50,gg = 150,precision = 5,dimension = 2)
    GA.start()

    print("")
    print("")
    multifunction = ThirdFunction(5)

    print("Multivariablefunction binary:")
    GA = GeneticAlgorithm(function = multifunction)
    GA.init(dd = -50,gg = 150,precision = 5,dimension = 5)
    GA.start()

    print("Multivariablefunction real:")
    GA = GeneticAlgorithm(solution_type = "real",function = multifunction)
    GA.init(dd = -50,gg = 150,precision = 5,dimension = 5)
    GA.start()

    print("")
    print("")
    schaffer = SchafferFunction(3)

    print("Schafferfunction binary:")
    GA = GeneticAlgorithm(function = schaffer)
    GA.init(dd = -50,gg = 150,precision = 5,dimension = 3)
    GA.start()

    print("Schafferfunction real:")
    GA = GeneticAlgorithm(solution_type = "real",function = schaffer)
    GA.init(dd = -50,gg = 150,precision = 5,dimension = 3)
    GA.start()

    print("")
    print("")
    almostSchaffer = AlmostSchafferFunction(3)

    print("AlmostSchafferFunction binary:")
    GA = GeneticAlgorithm(function = almostSchaffer)
    GA.init(dd = -50,gg = 150,precision = 5,dimension = 3)
    GA.start()

    print("AlmostSchafferFunction real:")
    GA = GeneticAlgorithm(solution_type = "real",function = almostSchaffer)
    GA.init(dd = -50,gg = 150,precision = 5,dimension = 3)
    GA.start()

    print("")
    print("")
    print("")
    print("")

    #DRUGI ZADATAK
    print("------------DRUGI ZADATAK--------------")

    N = [1,3,5,7]

    for n in N:
        schaffer = SchafferFunction(n)

        print("Schaffer with "+str(n)+" dimension:")
        GA = GeneticAlgorithm(function = schaffer)
        GA.init(dd = -50,gg = 150,precision = 5,dimension = n)
        GA.start()

        print("")
        print("")

        almostSchaffer = AlmostSchafferFunction(n)

        print("AlmostSchaffer with "+str(n)+" dimension:")
        GA = GeneticAlgorithm(function = almostSchaffer)
        GA.init(dd = -50,gg = 150,precision = 5,dimension = n)
        GA.start()

        print("")
        print("")

    #TRECI ZADATAK
    print("-------------TRECI ZADATAK------------")

    schaffer = SchafferFunction(3)

    print("Schaffer with 3 dimensions:")
    GA = GeneticAlgorithm(function = schaffer,number_of_evaluations =10**5)
    GA.init(dd = -50,gg = 150,precision = 4,dimension = 3)
    GA.start()

    print("")

    schaffer = SchafferFunction(6)

    print("Schaffer with 6 dimenisions:")
    GA = GeneticAlgorithm(function = schaffer,number_of_evaluations = 10**5)
    GA.init(dd = -50,gg = 150,precision = 4,dimension = 6)
    GA.start()

    print("")

    schaffer = SchafferFunction(3)

    print("Schaffer with 3 dimensions real:")
    GA = GeneticAlgorithm(function = schaffer,number_of_evaluations =10**5,solution_type = "real")
    GA.init(dd = -50,gg = 150,precision = 4,dimension = 3)
    GA.start()

    print("")

    schaffer = SchafferFunction(6)

    print("Schaffer with 6 dimenisions real:")
    GA = GeneticAlgorithm(function = schaffer,number_of_evaluations = 10**5,solution_type="real")
    GA.init(dd = -50,gg = 150,precision = 4,dimension = 6)
    GA.start()

    print("")

    schaffer = AlmostSchafferFunction(3)

    print("AlmostSchaffer with 3 dimensions:")
    GA = GeneticAlgorithm(function = schaffer,number_of_evaluations = 10**5)
    GA.init(dd = -50,gg = 150,precision = 4,dimension = 3)
    GA.start()

    print("")

    schaffer = AlmostSchafferFunction(6)

    print("AlmostSchaffer with 6 dimenisions:")
    GA = GeneticAlgorithm(function = schaffer,number_of_evaluations = 10**5)
    GA.init(dd = -50,gg = 150,precision = 4,dimension = 6)
    GA.start()

    print("")

    schaffer = AlmostSchafferFunction(3)

    print("AlmostSchaffer with 3 dimensions:")
    GA = GeneticAlgorithm(function = schaffer,number_of_evaluations = 10**5,solution_type = "real")
    GA.init(dd = -50,gg = 150,precision = 4,dimension = 3)
    GA.start()

    print("")

    schaffer = AlmostSchafferFunction(6)

    print("AlmostSchaffer with 6 dimenisions:")
    GA = GeneticAlgorithm(function = schaffer,number_of_evaluations = 10**5,solution_type = "real")
    GA.init(dd = -50,gg = 150,precision = 4,dimension = 6)
    GA.start()

    print("")

    values = []
    results = []
    for type in ["binary","real"]:
        for i in range(10):
            schaffer = AlmostSchafferFunction(3)

            GA = GeneticAlgorithm(function = schaffer,number_of_evaluations = 10**3,solution_type = "real")
            GA.init(dd = -50,gg = 150,precision = 4,dimension = 3)
            results.append(GA.start(shouldPrint = False)[1])

        values.append(results)

    print("Medians:")
    for value in values:
        array = np.asarray(value)
        print(np.median(array))

    print("")


    #CETVRTI ZADATAK
    print("------------CETVRTI ZADATAK------------")

    N = [30,50,100,200]
    results = []

    for n in N:
        schaffer = SchafferFunction(3)

        print("Schaffer with "+str(n)+" population:")
        GA = GeneticAlgorithm(population_size = n,function = schaffer)
        GA.init(dd = -50,gg = 150,precision = 5,dimension = 3)
        results.append(GA.start()[1])


    print("")
    index1 = results.index(max(results))
    print("Best population size: "+str(N[index1]))

    P = [0.1,0.3,0.6,0.9]
    results = []

    for p in P:
        schaffer = SchafferFunction(3)

        print("Schaffer with "+str(p)+" mutation probability:")
        GA = GeneticAlgorithm(population_size = N[index1],function = schaffer,mutation_probability = p)
        GA.init(dd = -50,gg = 150,precision = 5,dimension = 3)
        results.append(GA.start()[1])

    print("")
    index2 = results.index(max(results))
    print(index2)
    print("Best mutation probability: "+str(P[index2]))

    C = [0.1,0.3,0.6,0.9]
    results = []

    for p in C:
        schaffer = SchafferFunction(3)

        print("Schaffer with "+str(p)+" crossing probability:")
        GA = GeneticAlgorithm(population_size = N[index1],function = schaffer,mutation_probability = P[index2],crossing_probabilty = p)
        GA.init(dd = -50,gg = 150,precision = 5,dimension = 3)
        results.append(GA.start()[1])

    print("")
    index3 = results.index(max(results))
    print("Best crossing probability: "+str(C[index3]))

    print("")
    print("")

    N = [30,50,100,200]
    results = []

    values = []
    for n in N:
        results = []
        for i in range(20):
            schaffer = SchafferFunction(3)

            print("Schaffer with "+str(n)+" population:")
            GA = GeneticAlgorithm(population_size = n,function = schaffer)
            GA.init(dd = -50,gg = 150,precision = 5,dimension = 3)
            results.append(GA.start(shouldPrint = False)[1])

        values.append(results)

    print("Medians:")
    file = open('bok.txt','a')
    for value in values:
        array = np.asarray(value)
        for number in value:
            file.write(str(number)+"\n")
        file.write("\n")
        print(np.median(array))
    file.close()



    #PETI ZADATAK
    print("-----------------PETI ZADATAK----------------------")

    K = [3,5,7,11,13,15]

    for k in K:
        almostSchaffer = AlmostSchafferFunction(3)

        print("AlmostSchaffer with "+str(k)+" turn tables:")
        GA = GeneticAlgorithm(population_size = 20,function = almostSchaffer,turnir_size = k)
        GA.init(dd = -50,gg = 150,precision = 5,dimension = 3)
        GA.start()
        print("")

main()
