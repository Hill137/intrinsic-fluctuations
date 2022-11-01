'Algoritmo de Guillespie' 'Metodo directo'

#En este programa se implementa el agoritmo de Gillespie para reacciones quimicas
#comunes y convencionales.

# se importan las librerias que seran usadas
import numpy as np
import math 
import random as rnd
import matplotlib.pyplot as plt

#para implementar el algoritmo son necesarias algunas cantidades y definir algunos parametros
#condiciones iniciales de las variables involucradas y el tiempo

t0=0.0 #tiempo inicial
tf= 8 # tiempo final
Omega=60 #tamano del sistema
x=[0,Omega, 0] #vector de condiciones iniciales [X0,Y0,Z0]
N=2  #coeficiente de Hill
h= 4 #numero de veces para realizar la misma simulacion

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

def prod(x, n): # x!/(x-n)!
    l2=x-n
    if n==0:
        return 1
    if l2==0 or l2==1:
        return math.factorial(x)
    else:
        xprod= l2+1
        if n!=0 and l2>0:
            for i in range(1,n):
                xprod= xprod * (l2+1+i)
            return xprod
        elif l2<0:
             return 0

# elementos para realizar el proceso estocastico

def tazas(x,vector,k,j):  # taza de reaccion de una variable
    l1=len(x)
    l2=0
    l3=[]
    while l2<l1:
        l4=(prod(x[l2],vector[l2]))/(j**vector[l2])
        l2=l2+1
        l3.append(l4)
    l4=1
    l2=0
    while l2<l1:
        l4=l4*l3[l2]
        l2=l2+1
    return k*l4*j

def ni(x,A,B,k1,k2,j): #construimos el vector de propension del sistema
    l1=len(k1)
    l2=[]
    l3=[]
    l4=0
    while l4<l1:
        y1=tazas(x,A[l4],k1[l4],j)
        l2.append(y1)
        y2=tazas(x,B[l4],k2[l4],j)
        l3.append(y2)
        l4=l4+1
    return l2+l3

def dist_exp(a):# funcion para calcular tau 
    r = rnd.random()
    return -(1./(a))*math.log(r)
    
def numero(S1,a): # funcion que devuelve mu
    l1=0
    l2= S1[0]/a
    l3=0
    r=np.random.rand()
    while l3<1:
        if r<=l2 and S1[l1]!=0:
            l3=l3+1
            return l1
        else:
            l1=l1+1
            l2= l2 + S1[l1]/a
 
# la funcion que esta a continuacion nos permite realizar  el proceso estocastico hasta que haya 
#pasado T tiempo para un determinado tamano del sistema.

def ev(x,j): #x=condicion inicial, p= numero de pasos, j=tamano del sistema
    l1=len(x)
    Y = np.zeros([l1,1])
    t = []
    l2 = np.zeros([l1,1])
    l3=t0
    l4=0
    Y[:,0] = x[:]
    t.append(t0)
    while l3<tf:
        ni1 = ni(Y[:,l4],alpha,beta,k1,k2,j) #vector de propension
        a=sum(ni1)
        tau = dist_exp(a)
        mu = numero(ni1,a)
        l3=l3+tau
        t.append(l3)
        Y= np.append(Y, l2, axis = 1)
        Y[:,l4+1] = Y[:,l4] + S[:,mu]
        l4=l4+1
    return t ,  Y/j

#con la funcion que ya hemos programado, podemos repetirla algun numero determinado de veces,
# para ello claculamos la funcion anterior algun numero de veces y la graficamos 

def varios(x,j,q): #p= pasos  j= tamano del sistema, 
   # q= numero de veces que se realiza la misma simulacion 
    l1=0
    l2=len(x)
    l3=np.zeros(l2)
    while l1<q:
        T0, l3 =ev(x,j)
        for i in range(l2):
            plt.figure(i)
            plt.plot(T0,l3[i])
            plt.xlabel("t (tiempo)")
            plt.ylabel("x (concentracion de particulas)")
            plt.title("x" )
            plt.xlim(0, tf)
            plt.ylim(-0.005)
        l1=l1+1
    plt.legend()
    plt.show()

#como alternativa podemos realizar las simulacion para varios tamanos de sistema
# omega={a,b}
q=1000  #tamano de paso

def val(x,a,b):
    l1=a
    l2=len(x)
    while l1<b:
        x=[0,l1, 0] # condicion inicial
        T0, l3 =ev(x,l1)
        for i in range(l2):
            plt.figure(i)
            plt.plot(T0,l3[i])
            plt.xlabel("t (tiempo)")
            plt.ylabel("x (concentracion de particulas)")
            plt.title("x" )
            plt.xlim(0, tf)
            plt.ylim(-0.005)
        l1=l1+q
    plt.show()
    
#realizar las graficas de las variables
varios(x,Omega,h)
#val(x,60,8000)