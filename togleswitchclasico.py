'''
Togle switch clasico
'''

import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import axes3d

from scipy.integrate import odeint

def Hill(x,n):
	return 1 / (1+x**n)


def Lorenz(State,t,a1,b1,a2,b2):
	x1, y1, x2, y2 = State 

	x1_dot = a1*(Hill(y2,2))- a2*x1
	x2_dot = x1 - a2*x2
	y1_dot = b1*(Hill(x2,2)) - b2*y1
	y2_dot = y1 - b2*y2


	return x1_dot, y1_dot, x2_dot, y2_dot


# Parámetros
a1=b1 = 5
a2= b2 = 1

# Condición inicial 
x1o, y1o = 1, 0.
x2o, y2o = 0., 0.

t = np.arange(0,60,0.001)

# A resolver

State   = odeint(Lorenz,(x1o,y1o,x2o,y2o),t,args=(a1,b1,a2, b2))
x1, y1, x2, y2 = State.T


# A graficar
fig = plt.figure(1,figsize=(20,10))
plt.subplot(311)
plt.plot(t,x1,color='red', label='x1')
plt.plot(t,y1,color='blue', label='y1')
plt.legend()
plt.subplot(313)
plt.plot(t,x2,color='red', label='x2')
plt.plot(t,y2,color='blue', label='y2')

#fig2 = plt.figure(2,figsize=(30,10))
#ax=axes3d.Axes3D(fig2)
#ax.plot(x1,y1,z1, color='green')
#ax.plot(x2,y2,z2, color='red')
plt.legend()
plt.show()