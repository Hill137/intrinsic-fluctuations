import numpy as np
import math
import random as rnd
import matplotlib.pyplot as plt


#condiciones iniciales
t0=0.0
X0=0
Y0=0


#parametros
alpha = 1
beta = 1.0
gamma = 1
delta = 1
Omega= 100
Steps = 50

S =[[1.,0,-1.,0],[0,1,0,-1]]



def dist_exp(a):
    r = rnd.random()
    return -(1./a)*math.log(r)
    
def dist_reaction(ni,A):
    r = rnd.random() 
    if (r < ni[0]/A):
        return 0
    elif (r < (ni[0] + ni[1])/A):
        return 1
    elif (r < (ni[0] + ni[1] + ni[2])/A):
        return 2
    elif (r <= (ni[0] + ni[1] + ni[2] + ni[3])/A):
        return 3
    


def ev(p,j):
    Y = np.zeros([2,p+1])
    t = np.zeros(p+1)
    Y[0]=X0
    Y[1]=Y0
    for i in range(p):
        ni = np.array([gamma, gamma,beta*Y[0][i]/j,beta*Y[1][i]/j ])
        a = sum(ni)
        tau = dist_exp(a)
        mu = dist_reaction(ni,a)
        Y[0][i+1] = Y[0][i] + S[0][mu]
        Y[1][i+1] = Y[1][i] + S[1][mu]
        t[i+1] = t[i] + tau
    return t , Y[0]/j, Y[1]/j


def val(p,a,b):
    plt.figure(3)
    for j in range(a,b+1):
        T0, Y01, Y1 =ev(p,j)
        plt.plot(T0,Y01,Y1)
    plt.xlabel("t")
    plt.ylabel("x")
    plt.title("Proceso de vida y muerte")
    plt.show()
        

val(5000,10,15)
