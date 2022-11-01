'Redes geneticas estocasticas' 'represilador'

#En este programa se implementa el agoritmo de Gillespie para reacciones bioquimicas,
#especificamente para una red de transcripcion/traduccion.

# se importan las librerias que seran usadas
import numpy as np
import math 
import random as rnd
import matplotlib.pyplot as plt

#para implementar el algoritmo son necesarias algunas cantidades y definir algunos parametros
#condiciones iniciales de las variables involucradas y el tiempo

A=np.array([ [0,0,-1],
             [-1,0,0],
             [0,-1,0]]) # matriz de conexiones

t0=0.0 #tiempo inicial
tf= 150 # tiempo final
t1=0.1 # tamano de paso
N=[2,2,2]  #coeficientes de de Hill para cada uno de los nodos
K=[1,1,1] # constantes de Michelson para cada nodo
h= 1 #numero de veces para realizar la misma simulacion
n= len(A[0]) #numero de nodos
n2= 2 #numero de elementos en cada nodo
Omega=80 # tamano del sistema

#condiciones iniciales para cada nodo
x=[[Omega,0], #vector de condiciones iniciales [ARm0, P0]
    [0,0],                              #  [ARm1, P1]
    [0,0]]                              #  [ARm2, P2]


#parametros, los vectores para cada uno de los nodos 
p1=1
p2=10
k1=[p2,p1,1,p1]  #[sintesis o degradacion de ARm a travez de Hill, 
                 #degradacion ARm, sintesis proteina, degradacion proteina]
k2=[p2,p1,1,p1]  
k3=[p2,p1,1,p1]  

k=[k1,k2,k3] # vector de parametros para todo el sistema

#contruimos la matriz estequiometrica para cada nodo 'evitar editar a partir de esta linea'
S= [[1,-1,0, 0],
    [ 0,0 ,1,-1],]

#construimos la matriz estequimetrica del sistema
def ST(S,n):
    l1=np.zeros([2,4])
    S1=S
    l2=l1
    l3=S
    for i in range(n-1):
        S1= np.append(S1, l2, axis = 1)
        l3=np.append(l1,l3, axis=1)
        S1= np.r_[S1,l3]
        l2= np.r_[l2,l1]
    return S1

S=ST(S,n)

#construimos lista nodos conectados entre si

def vec(A):
    l2=[]
    for i in range(n):
        for j in range(n):
            if A[i][j]!=0:
                l2.append(j)
    return l2

vecco= vec(A)

# miscelanea de funciones usadas

def H_d(x,n,j): #funcion de hill determinista
    return 1/(1 + (x/j)**n)

def H_ap(x,n,j): # funcion de hill con aproximaciones
    return (1 )/(1 + (x/j)**n + (n/(2*j))*(n-1)*(x/j)**(n-1))

# construimos el vector de propension para cada nodo

def ninodo(y1,x1,k1,j,N): #funcion de propension para cada nodo
    l1= [k1[0]*(H_ap(y1,N,j)), k1[1]*x1[0]/j, k1[2]*x1[0]/j, k1[3]*x1[1]/j]
    return l1

def ni(x,k,j): #funcion de porpension para todo el sistema
    l1=[]
    for i in range(n):
        l2=ninodo(x[vecco[i]][1],x[i],k[i],j,N[i])
        l1=l1+l2
    return l1

def dist_exp(a):# funcion para calcular tau 
    r = rnd.random()
    return -(1./(a))*math.log(r)
    
def numero(S1,a): # funcion que devuelve mu
    l1=0
    l2= S1[0]/a
    l3=0
    r=np.random.rand()
    while l1<len(S1) and l3==0:
        if r<=l2 and S1[l1]!=0:
            l3=l3+1
            return l1
        else:
            l1=l1+1
            l2= l2 + S1[l1]/a
 
# la funcion que esta a continuacion nos permite realizar  el proceso estocastico hasta que haya 
#pasado T tiempo para un determinado tamano del sistema.

def ev(x,k,j): #x=condicion inicial, j=tamano del sistema
    t = np.arange(0,tf+t1,t1)
    l3=0
    Z=np.zeros([n,n2,len(t)])
    Y = x
    Z[:,:,0] = x
    l6=1
    while l6<len(t) or l3<tf:
        if l3<(l6)*t1:
            ni1 = ni(Y,k,j) #vector de propension
            a=sum(ni1)
            tau = dist_exp(a)/j
            mu = numero(ni1,a)
            l5=0
            for j1 in range(n):
                for j2 in range(n2):
                   Y[j1][j2] = Y[j1][j2] + S[l5][mu]
                   l5=l5+1
            l3=l3+tau
        else:
            Z[:,:,l6]=Y
            l6=l6+1
    return t ,  Z/j

#con la funcion que ya hemos programado, podemos repetirla algun numero determinado de veces,
# para ello claculamos la funcion anterior algun numero de veces y la graficamos 

def varios(x,k,j,q): #p= pasos  j= tamano del sistema, 
   # q= numero de veces que se realiza la misma simulacion 
    l1=0
    while l1<q:
        T0, l3 =ev(x,k,j)
        for i in range(n):
            for j1 in range(n2):
                if j1==0:
                    plt.figure(1)
                    plt.plot(T0,l3[i][j1], label= i)
                    plt.xlabel("t (tiempo)")
                    plt.ylabel("x (concentracion de ARm)")
                    plt.title("concentracion de las especies de ARm" )
                    plt.xlim(0, tf)
                    plt.ylim(-0.0005)
                    plt.legend()
                if j1==1:
                    plt.figure(2)
                    plt.plot(T0,l3[i][j1], label= i)
                    plt.xlabel("t (tiempo)")
                    plt.ylabel("x (concentracion de proteinas)")
                    plt.title("concentracion de proteinas" )
                    plt.xlim(0, tf)
                    plt.ylim(-0.0001)
                    plt.legend()
        l1=l1+1
    plt.show()

#realizar las graficas de las variables
varios(x,k,Omega,h)