'''
Sistema_de_Lorenz.py
'''

import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import axes3d

from scipy.integrate import odeint

def Hill(x,n):
	return 1 / (1+x**n)

def H_ap(x,n,j): # funcion de hill con aproximaciones
    return 1/(1 + (x)**n + (n/(2*j))*(n-1)*(x)**(n-1))


def Lorenz(State,t,a1,b1,c1,a2,b2,c2):
	x1, y1, z1, x2, y2, z2 = State 

	x1_dot = a1*(H_ap(z2,2,j))- a2*x1
	x2_dot = x1 - a2*x2
	y1_dot = b1*(H_ap(x2,2,j)) - b2*y1
	y2_dot = y1 - b2*y2
	z1_dot = c1*(H_ap(y2,2,j)) - c2*z1
	z2_dot = z1 - c2*z2

	return x1_dot, y1_dot, z1_dot, x2_dot, y2_dot, z2_dot


# Parámetros
a1=b1= c1 = 10
a2= b2= c2 = 1
j=60

# Condición inicial 
x1o, y1o, z1o = 1, 0., 0.
x2o, y2o, z2o = 0., 0., 0.

t = np.arange(0,80,0.001)

# A resolver

State   = odeint(Lorenz,(x1o,y1o,z1o,x2o,y2o,z2o),t,args=(a1,b1,c1,a2, b2, c2))
x1, y1, z1, x2, y2, z2 = State.T


# A graficar
fig = plt.figure(1,figsize=(30,10))
plt.subplot(311)
plt.plot(t,x1,color='red', label='x1')
plt.plot(t,y1,color='blue', label='y1')
plt.plot(t,z1,color='green', label='z1')
plt.xlabel("t (tiempo)")
plt.ylabel("x (concentracion de ARm)")
plt.title("concentracion de las especies de ARm" )
plt.legend()
plt.subplot(313)
plt.plot(t,x2,color='red', label='x2')
plt.plot(t,y2,color='blue', label='y2')
plt.plot(t,z2,color='green', label='z2')
plt.xlabel("t (tiempo)")
plt.ylabel("x (concentracion de la proteina)")
plt.title("concentracion de las especies de proteina" )

#fig2 = plt.figure(2,figsize=(30,10))
#ax=axes3d.Axes3D(fig2)
#ax.plot(x1,y1,z1, color='green')
#ax.plot(x2,y2,z2, color='red')
plt.legend()
plt.show()