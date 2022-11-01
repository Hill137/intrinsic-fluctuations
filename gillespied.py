'Algoritmo de Guillespie' 'Metodo directo'

#En este programa se implementa el agoritmo de Gillespie para reacciones quimicas comunes y convencionales.

# se importan las librerias que seran usadas
import numpy as np
import math 
import random as rnd
import matplotlib.pyplot as plt

#para implementar el algoritmo son necesarias algunas cantidades y definir algunos parametros
#condiciones iniciales de las variables involucradas y el tiempo
t0=0.0 #tiempo inicial
tf= 500 # tiempo final
Omega=60 #tamano del sistema
x=[0,Omega, 0] #vector de condiciones iniciales [X0,Y0,Z0]
N=3  #coeficiente de Hill
P=750 #numero de pasos

#parametros
k1=[1,1]  #vector de k+
k2=[0.1,0.1]  #vector de k-

#matrices de coeficientes estequiometricos
alpha= [[0,0,0],
        [N,1,0]]

beta= [[1,0,0],
       [0,0,1]]

#contruimos la matriz estequiometrica 'evitar editar a partir de esta linea'
S= np.subtract( beta,alpha)
S=S.T
S= np.append(S, -S, axis = 1)

# miscelanea de funciones usadas

def fac(n): #x!
    if n<0:
        return 0
    else:
        l1=0
        l2=1
        while l1 < n:
          l1=l1+1
          l2=l1*l2
        return l2

def prod(x, n): # x!/(x-n)!
    l2=x-n
    if n==0:
        return 1
    if l2==0 or l2==1:
        return fac(x)
    else:
        xprod= l2+1
        if n!=0 and l2>0:
            for i in range(1,n):
                xprod= xprod * (l2+1+i)
            return xprod
        elif l2<0:
             return 0
        else:
            return 1

def elevar(x,n): #x**n
    if n==0:
        return 1
    else:
        return x**n

# elementos para realizar el proceso estocastico

def tazas(x,vector,k,Omega):  # taza de reaccion de una variable
    l1=len(x)
    l2=0
    l3=[]
    while l2<l1:
        l4=(prod(x[l2],vector[l2]))/(elevar(Omega,vector[l2]))
        l2=l2+1
        l3.append(l4)
    l4=1
    l2=0
    while l2<l1:
        l4=l4*l3[l2]
        l2=l2+1
    return k*l4

def ni(x,A,B,k1,k2,Omega): #construimos el vector de propension del sistema
    l1=len(k1)
    l2=[]
    l3=[]
    l4=0
    while l4<l1:
        y1=tazas(x,A[l4],k1[l4],Omega)
        l2.append(y1)
        y2=tazas(x,B[l4],k2[l4],Omega)
        l3.append(y2)
        l4=l4+1
    return l2+l3

def dist_exp(S1):# funcion para calcular tau 
    a=sum(S1)
    r = rnd.random()
    return -(1./a)*math.log(r)
    
def numero(S1): # funcion que devuelve mu
    a=sum(S1)
    l1=0
    l2= S1[0]/a
    l3=0
    r=np.random.rand()
    while l1<len(S1):
        if l3==0:
            if r<=l2:
                if S1[l1]!=0:
                    l3=l3+1
                    return (l1)
            else:
               l1=l1+1
               l2= l2 + S1[l1]/a
 
# la funcion que esta a continuacion nos permite realizar  el proceso estocastico hasta que hayan ocurrido p 
#reacciones continuas para un determinado tamano del sistema.

def ev(x,p,j): #x=condicion inicial, p= numero de pasos, j=tamano del sistema
    l1=len(x)
    Y = np.zeros([l1,p+1])
    t = np.zeros(p+1)
    for j1 in range(l1):
            Y[j1][0] = x[j1]
    t[0]=t0
    for i in range(p):
        ni1 = ni(x,alpha,beta,k1,k2,j) #vector de propension
        tau = dist_exp(ni1)
        mu = numero(ni1)
        t[i+1] = t[i] + tau
        l2=[]
        for j1 in range(l1):
            Y[j1][i+1] = Y[j1][i] + S[j1][mu]
            l2.append(Y[j1][i+1])
        x=l2
    return t ,  Y/j

#con la funcion que ya hemos programado, podemos repetirla algun numero determinado de veces,
# para ello claculamos la funcion anterior algun numero de veces y la graficamos 

def varios(x,p,j,q): #p= pasos  j= tamano del sistema, q= numero de veces que se realiza la misma simulacion 
    l1=0
    l2=len(x)
    l3=np.zeros(l2)
    while l1<q:
        T0, l3 =ev(x,p,j)
        for i in range(l2):
            plt.figure(i)
            plt.plot(T0,l3[i])
            plt.xlabel("t (tiempo)")
            plt.ylabel("x (concentracion de particulas)")
            plt.title("x" )
        l1=l1+1
    plt.legend()
    # Guardar grafico como imagen PNG.
    plt.savefig("hill.png")
    plt.show()

# como alternativa podemos realizar las simulaciones en un intervalo de tamano
# omega={a,b}
def val(x,p,a,b):
    plt.figure(3)
    for j in range(a,b+1):
        T0, Y01, Y1 =ev(x,p,j)
        plt.plot(T0,Y01)
        plt.plot(T0,Y1)
    plt.xlabel("t")
    plt.ylabel("x (concentracion de particulas)")
    plt.title("x")
    plt.show()

#realizar las graficas de las variables
varios(x,P,Omega,10)