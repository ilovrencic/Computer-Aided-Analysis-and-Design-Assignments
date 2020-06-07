from Matrix import Matrix
import matplotlib.pyplot as plt
import math
iter = 100

class FirstRecall:
    def __init__(self):
        pass

    def get_b1(self,t):
        return 1

    def get_b2(self,t):
        return 1
class SecondRecall:
    def __init__(self):
        pass

    def get_b1(self,t):
        return t

    def get_b2(self,t):
        return t

def euler_method(A,x,T,t_max,B = None,r = None,analysis = False):
    print("Euler method:")
    current_t = 0
    I = Matrix()
    I.setMatrix([[1,0],[0,1]])
    M = I + A*T
    X = Matrix()
    X.setMatrix(x.matrix)
    apporx_points = []

    total_error = 0
    x1_next_0 = X[0,0]
    x2_next_0 = X[0,1]

    iteration = 0
    while(current_t <= t_max):
        iteration += 1
        current_t += T
        x1_next = M[0,0]*X[0,0]+M[0,1]*X[0,1]
        x2_next = M[1,0]*X[0,0]+M[1,1]*X[0,1]

        if(r != None):
            N = B*T
            x1_next += N[0,0]*r.get_b1(current_t)+N[0,1]*r.get_b2(current_t)
            x2_next += N[1,0]*r.get_b1(current_t)+N[1,1]*r.get_b2(current_t)

        if(iteration%iter == 0):
            print(x1_next,x2_next)

        if(analysis):
            x1_next_analysis = x1_next_0*math.cos(current_t)+x2_next_0*math.sin(current_t)
            x2_next_analysis = x2_next_0*math.cos(current_t)-x1_next_0*math.sin(current_t)
            total_error += abs(x1_next-x1_next_analysis)+abs(x2_next-x2_next_analysis)

        apporx_points.append([x1_next,x2_next])
        X[0,0] = x1_next
        X[0,1] = x2_next

    X.roundMatrix()
    print(X)
    print("")
    return total_error,apporx_points

def reverse_euler_method(A,x,T,t_max,B = None,r = None,analysis = False):
    print("Reverse euler method:")
    current_t = 0
    I = Matrix()
    I.setMatrix([[1,0],[0,1]])
    P = I - A*T
    P.inverse()
    X = Matrix()
    X.setMatrix(x.matrix)
    apporx_points = []

    total_error = 0
    x1_next_0 = X[0,0]
    x2_next_0 = X[0,1]

    iteration = 0
    while(current_t <= t_max):
        iteration += 1
        current_t += T
        x1_next = P[0,0]*X[0,0]+P[0,1]*X[0,1]
        x2_next = P[1,0]*X[0,0]+P[1,1]*X[0,1]

        if(r != None):
            Q = P*B*T
            x1_next += Q[0,0]*r.get_b1(current_t+T)+Q[0,1]*r.get_b2(current_t+T)
            x2_next += Q[1,0]*r.get_b1(current_t+T)+Q[1,1]*r.get_b2(current_t+T)

        if(iteration%iter == 0):
            print(x1_next,x2_next)

        if(analysis):
            x1_next_analysis = x1_next_0*math.cos(current_t)+x2_next_0*math.sin(current_t)
            x2_next_analysis = x2_next_0*math.cos(current_t)-x1_next_0*math.sin(current_t)
            total_error += abs(x1_next-x1_next_analysis)+abs(x2_next-x2_next_analysis)

        apporx_points.append([x1_next,x2_next])
        X[0,0] = x1_next
        X[0,1] = x2_next

    X.roundMatrix()
    print(X)
    print("")
    return total_error,apporx_points

def trapez_method(A,x,T,t_max,B = None,r = None,analysis = False):
    print("Trapez method:")
    current_t = 0
    I = Matrix()
    I.setMatrix([[1,0],[0,1]])
    S = (I - A*(T/2))
    S.inverse()
    R = S*(I+A*(T/2))
    iteration = 0
    X = Matrix()
    X.setMatrix(x.matrix)
    apporx_points = []

    total_error = 0
    x1_next_0 = X[0,0]
    x2_next_0 = X[0,1]

    while(current_t <= t_max):
        iteration += 1
        current_t += T
        x1_next = R[0,0]*X[0,0]+R[0,1]*X[0,1]
        x2_next = R[1,0]*X[0,0]+R[1,1]*X[0,1]

        if(r != None):
            V = S*B*(T/2)
            x1_next += V[0,0]*(r.get_b1(current_t)+r.get_b1(current_t+T))+V[0,1]*(r.get_b2(current_t)+r.get_b2(current_t+T))
            x2_next += V[1,0]*(r.get_b1(current_t)+r.get_b1(current_t+T))+V[1,1]*(r.get_b2(current_t)+r.get_b2(current_t+T))

        if(iteration%iter == 0):
            print(x1_next,x2_next)

        if(analysis):
            x1_next_analysis = x1_next_0*math.cos(current_t)+x2_next_0*math.sin(current_t)
            x2_next_analysis = x2_next_0*math.cos(current_t)-x1_next_0*math.sin(current_t)
            total_error += abs(x1_next-x1_next_analysis)+abs(x2_next-x2_next_analysis)

        apporx_points.append([x1_next,x2_next])
        X[0,0] = x1_next
        X[0,1] = x2_next

    X.roundMatrix()
    print(X)
    print("")
    return total_error,apporx_points

def runge_kutta(A,x,T,t_max,B = None,r = None,analysis = False):
    print("Runge kutta:")
    current_t = 0
    I = Matrix()
    I.setMatrix([[1,0],[0,1]])
    A = A*T
    M = I+A+(A**2)*0.5+(A**3)*(1/6)+(A**4)*(1/24)
    iteration = 0
    X = Matrix()
    X.setMatrix(x.matrix)
    apporx_points = []

    total_error = 0
    x1_next_0 = X[0,0]
    x2_next_0 = X[0,1]

    while(current_t <= t_max):
        iteration += 1
        current_t += T
        x1_next = M[0,0]*X[0,0]+M[0,1]*X[0,1]
        x2_next = M[1,0]*X[0,0]+M[1,1]*X[0,1]

        if( r != None ):
            C = B*(T/6)
            N = (I+A+(A**2)*0.5+(A**3)*(1/4))*C
            R = (I*4+A*2+(A**2)*0.5)*C

            x1_next += N[0,0]*r.get_b1(current_t)+N[0,1]*r.get_b2(current_t)
            x2_next += N[1,0]*r.get_b1(current_t)+N[1,1]*r.get_b2(current_t)

            x1_next += R[0,0]*r.get_b1(current_t+(T/2))+R[0,1]*r.get_b2(current_t+(T/2))
            x2_next += R[1,0]*r.get_b1(current_t+(T/2))+R[1,1]*r.get_b2(current_t+(T/2))

            x1_next += C[0,0]*r.get_b1(current_t+T)+C[0,1]*r.get_b2(current_t+T)
            x2_next += C[1,0]*r.get_b1(current_t+T)+C[1,1]*r.get_b2(current_t+T)

        if(iteration%iter == 0):
            print(x1_next,x2_next)

        if(analysis):
            x1_next_analysis = x1_next_0*math.cos(current_t)+x2_next_0*math.sin(current_t)
            x2_next_analysis = x2_next_0*math.cos(current_t)-x1_next_0*math.sin(current_t)
            total_error += abs(x1_next-x1_next_analysis)+abs(x2_next-x2_next_analysis)

        apporx_points.append([x1_next,x2_next])
        X[0,0] = x1_next
        X[0,1] = x2_next

    X.roundMatrix()
    print(X)
    print("")
    return total_error,apporx_points

def pece2(A,x,T,t_max,B = None,r = None,analysis = False,k = 2):
    print("PECE2:")
    current_t = 0
    X = Matrix()
    X.setMatrix(x.matrix)
    iteration = 0
    apporx_points = []

    total_error = 0
    x1_next_0 = X[0,0]
    x2_next_0 = X[0,1]

    x1_dot_curr = e_get_dot_curr1(A,X,T,B,r)
    x2_dot_curr = e_get_dot_curr2(A,X,T,B,r)
    while(current_t <= t_max):
        iteration += 1
        current_t += T

        #EULEROV
        x1_next = X[0,0]+T*x1_dot_curr
        x2_next = X[0,1]+T*x2_dot_curr

        x0 = X[0,0]
        x1 = X[0,1]

        X[0,0] = x1_next
        X[0,1] = x2_next

        x1_dot_next = e_get_dot_curr1(A,X,current_t+T,B,r)
        x2_dot_next = e_get_dot_curr2(A,X,current_t+T,B,r)

        for i in range(k):
            #REVERSE EULER
            x1_next = x0+T*x1_dot_next
            x2_next = x1+T*x2_dot_next

            X[0,0] = x1_next
            X[0,1] = x2_next

            x1_dot_next = e_get_dot_curr1(A,X,current_t+T,B,r)
            x2_dot_next = e_get_dot_curr2(A,X,current_t+T,B,r)

        x1_dot_curr = e_get_dot_curr1(A,X,current_t+T,B,r)
        x2_dot_curr = e_get_dot_curr2(A,X,current_t+T,B,r)

        if(iteration%iter == 0):
            print(x1_next,x2_next)

        if(analysis):
            x1_next_analysis = x1_next_0*math.cos(current_t)+x2_next_0*math.sin(current_t)
            x2_next_analysis = x2_next_0*math.cos(current_t)-x1_next_0*math.sin(current_t)
            total_error += abs(x1_next-x1_next_analysis)+abs(x2_next-x2_next_analysis)

        apporx_points.append([x1_next,x2_next])

    X.roundMatrix()
    print(X)
    print("")
    return total_error,apporx_points

def pece(A,x,T,t_max,B = None,r = None,analysis = False):
    print("PECE")
    current_t = 0
    X = Matrix()
    X.setMatrix(x.matrix)
    iteration = 0
    apporx_points = []

    total_error = 0
    x1_next_0 = X[0,0]
    x2_next_0 = X[0,1]

    x1_dot_curr = e_get_dot_curr1(A,X,T,B,r)
    x2_dot_curr = e_get_dot_curr2(A,X,T,B,r)
    while(current_t <= t_max):
        iteration += 1
        current_t += T

        #EULEROV
        x1_next = X[0,0]+T*x1_dot_curr
        x2_next = X[0,1]+T*x2_dot_curr

        x0 = X[0,0]
        x1 = X[0,1]

        X[0,0] = x1_next
        X[0,1] = x2_next

        x1_dot_next = e_get_dot_curr1(A,X,current_t+T,B,r)
        x2_dot_next = e_get_dot_curr2(A,X,current_t+T,B,r)

        #TRAPEZ
        x1_next = x0+(T/2)*(x1_dot_curr+x1_dot_next)
        x2_next = x1+(T/2)*(x2_dot_curr+x2_dot_next)

        X[0,0] = x1_next
        X[0,1] = x2_next

        x1_dot_curr = e_get_dot_curr1(A,X,current_t+T,B,r)
        x2_dot_curr = e_get_dot_curr2(A,X,current_t+T,B,r)

        if(iteration%iter == 0):
            print(x1_next,x2_next)

        if(analysis):
            x1_next_analysis = x1_next_0*math.cos(current_t)+x2_next_0*math.sin(current_t)
            x2_next_analysis = x2_next_0*math.cos(current_t)-x1_next_0*math.sin(current_t)
            total_error += abs(x1_next-x1_next_analysis)+abs(x2_next-x2_next_analysis)

        apporx_points.append([x1_next,x2_next])

    X.roundMatrix()
    print(X)
    print("")
    return total_error,apporx_points

def analitic_solution(X,T,t_max):
    current_t = 0
    points = []
    x1_next_0 = X[0,0]
    x2_next_0 = X[0,1]
    while(current_t <= t_max):
        current_t += T
        x1_next_analysis = x1_next_0*math.cos(current_t)+x2_next_0*math.sin(current_t)
        x2_next_analysis = x2_next_0*math.cos(current_t)-x1_next_0*math.sin(current_t)
        points.append([x1_next_analysis,x2_next_analysis])
    return points

def e_get_dot_curr1(A,X,t,B = None,r = None):
    if( r == None ):
        return A[0,0]*X[0,0]+A[0,1]*X[0,1]
    else:
        return A[0,0]*X[0,0]+A[0,1]*X[0,1]+B[0,0]*r.get_b1(t)+B[0,1]*r.get_b2(t)

def e_get_dot_curr2(A,X,t,B = None,r = None):
    if( r == None ):
        return A[1,0]*X[0,0]+A[1,1]*X[0,1]
    else:
        return A[1,0]*X[0,0]+A[1,1]*X[0,1]+B[1,0]*r.get_b1(t)+B[1,1]*r.get_b2(t)

def main():
    print("------------PRVI ZADATAK----------------")
    X = Matrix()
    X.setMatrix([[1,1]])
    A = Matrix()
    A.setMatrix([[0,1],[-1,0]])
    analysis_points = analitic_solution(X,T = 0.01,t_max = 10)
    euler,e_points = euler_method(A,X,T = 0.01,t_max = 10,analysis = True)
    reverse,r_points = reverse_euler_method(A,X,T = 0.01,t_max = 10,analysis = True)
    trapez,t_points = trapez_method(A,X,T = 0.01,t_max = 10,analysis = True)
    runge,k_points = runge_kutta(A,X,T = 0.01,t_max = 10,analysis = True)
    pece2_error,p2_points = pece2(A,X,T = 0.01,t_max = 10,analysis = True,k = 2)
    pece_error,p_points = pece(A,X,T = 0.01,t_max = 10,analysis = True)
    print("")
    print("Euler error:",euler)
    print("Reverse euler error:",reverse)
    print("Trapez error:",trapez)
    print("Runge error:",runge)
    print("PECE2 error:",pece2_error)
    print("PECE error:",pece_error)
    print("")

    print("------------DRUGI ZADATAK----------------")
    X = Matrix()
    X.setMatrix([[1,-2]])
    A = Matrix()
    A.setMatrix([[0,1],[-200,-102]])
    euler_method(A,X,T = 0.01,t_max = 10)
    reverse_euler_method(A,X,T = 0.1,t_max = 10)
    trapez_method(A,X,T = 0.1,t_max = 10)
    runge_kutta(A,X,T = 0.1,t_max = 10)
    pece2(A,X,T = 0.01,t_max = 10)
    pece(A,X,T = 0.01,t_max = 10)
    print("")
    print("")

    print("------------TRECI ZADATAK----------------")
    X = Matrix()
    X.setMatrix([[1,3]])
    A = Matrix()
    A.setMatrix([[0,-2],[1,-3]])
    B = Matrix()
    B.setMatrix([[2,0],[0,3]])
    r = FirstRecall()
    euler_method(A,X,T = 0.01,t_max = 10,B = B,r = r)
    reverse_euler_method(A,X,T = 0.01,t_max = 10,B = B,r = r)
    trapez_method(A,X,T = 0.01,t_max = 10,B = B,r = r)
    runge_kutta(A,X,T = 0.01,t_max = 10,B = B,r = r)
    pece2(A,X,T = 0.01,t_max = 10,B = B,r = r)
    pece(A,X,T = 0.01,t_max = 10,B = B,r = r)
    print("")
    print("")

    print("----------CETVRTI ZADATAK---------------")
    X = Matrix()
    X.setMatrix([[-1,3]])
    A = Matrix()
    A.setMatrix([[1,-5],[1,-7]])
    B = Matrix()
    B.setMatrix([[5,0],[0,3]])
    r = SecondRecall()
    euler_method(A,X,T = 0.01,t_max = 1,B = B,r = r)
    reverse_euler_method(A,X,T = 0.01,t_max = 1,B = B,r = r)
    trapez_method(A,X,T = 0.01,t_max = 1,B = B,r = r)
    runge_kutta(A,X,T = 0.01,t_max = 1,B = B,r = r)
    pece2(A,X,T = 0.01,t_max = 1,B = B,r = r)
    pece(A,X,T = 0.01,t_max = 1,B = B,r = r)
    print("")
    print("")

    plt.figure(figsize=(15,10))
    plt.subplot(3,3,1)
    plt.title("Euler:")
    plt.plot(e_points,label = "Approximation")
    plt.plot(analysis_points,label = "Real")
    plt.legend()

    plt.subplot(3,3,2)
    plt.title("Reverse euler:")
    plt.plot(r_points,label = "Approximation")
    plt.plot(analysis_points,label = "Real")
    plt.legend()

    plt.subplot(3,3,3)
    plt.title("Trapez:")
    plt.plot(t_points,label = "Approximation")
    plt.plot(analysis_points,label = "Real")
    plt.legend()

    plt.subplot(3,3,4)
    plt.title("Runge Kutta:")
    plt.plot(k_points,label = "Approximation")
    plt.plot(analysis_points,label = "Real")
    plt.legend()

    plt.subplot(3,3,5)
    plt.title("PECE2:")
    plt.plot(p2_points,label = "Approximation")
    plt.plot(analysis_points,label = "Real")
    plt.legend()

    plt.subplot(3,3,6)
    plt.title("PECE:")
    plt.plot(p_points,label = "Approximation")
    plt.plot(analysis_points,label = "Real")
    plt.legend()
    plt.show()



main()
