'Programa para modelar la distribucion de probabilidad'
'del proceso de vida y muerte'

#Se importan las librerias para usar
import numpy as np
from numpy import exp,arange
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from scipy import special
import math

#funciones para evaluar la funcion de probabilidad
def A(x,b,d):
    l1= min(x,d)
    l2= abs(d-x)
    return ((special.hyp1f1(-l1,l2+1, b))*((-1/b)**(l1)))/(math.factorial(l1)* math.factorial(l2))

def P(x,y,t,a,b): # funcion de probabilidad
    l2=0
    for d in range(50):
        l1= math.factorial(d)*math.factorial(y)*exp(-b*d*t)*exp(- a/b)* ((a/b)**(x+d))*A(x,a/b,d)*A(y,a/b,d)
        l2=l2+l1
    return l2

def evaluar(x,y,t,a,b): #funcion para crear tuplas [t,x,P]
    l1=[]
    for i in t:
        for j in x:
            l2= P(j,y,i,a,b)
            l1.append([i,j,l2])
    return l1

# parametros
a=3
b=1
y=4  #pocicion inicial
t0=0   #tiempo inicial
tf= 2.5 #tiempo final
xf=10 # x final

# se crean los puntos para evaluar la funcion
x = range(xf+1)
t = arange(t0,tf,0.004)

#Se evalua la funcion segun los valores de x y t
puntos = evaluar(x,y,t,a,b)

# GRAFICA
figura = plt.figure(1, figsize=[8,8])
grafica = figura.add_subplot(111,projection = '3d')

# xi = [:,0] ; yi = [:,1], zi = [:,2]
# selecciona columnas, use la transpuesta de puntos
[xi, yi , zi] = np.transpose(puntos)

grafica.scatter(xi,yi,zi,
                c = 'blue',
                marker='o',
                )
grafica.set_title('Distribucion de probabilidad')
grafica.set_xlabel('t (tiempo en segundos)')
grafica.set_ylabel('x (número de elementos)')
grafica.set_zlabel('P(x,t)')

#limitar el rango de valores que aparecen
plt.xlim(-0.0001, tf)
plt.ylim(-0.0001)

plt.show()