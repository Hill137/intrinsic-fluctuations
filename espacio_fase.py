# Copyright Jorge Velazquez Castro 2018
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
# 
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import matplotlib
matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)
matplotlib.rc('axes', labelsize=18)
#parameteres
param = np.zeros(3)
param[0] = alpha = .183
param[1] = f =13
param[2] = n = 8
#El sistema 
def sistema(y, t, param):
    X1 = y[0]
    X2 = y[1]
    alpha = param[0]
    f = param[1]
    n = param[2]
    dX1 = alpha*(1. + f * np.power(X2,n))/(1.+ np.power(X2,n)) - X1
    dX2 = alpha*(1.+ f * np.power(X1,n))/(1.+ np.power(X1,n)) - X2
    return np.array([ dX1, dX2])
#Condiciones iniciales
X1i=0
X2i=0
tspan = np.linspace(0,600,10000)
ys = odeint(sistema, [X1i,X2i] , tspan, args=(param,))
#plt.plot(ys[:,0], ys[:,1], 'b-') # path
#plt.plot([ys[0,0]], [ys[0,1]], 'o') # start
#plt.plot([ys[-1,0]], [ys[-1,1]], 's') # end
#plt.show()
############
x_range=3
y_range=3
y1 = np.linspace(0, x_range, 30)
y2 = np.linspace(0, y_range, 30)
Y1, Y2 = np.meshgrid(y1, y2)
t = 0
u, v = np.zeros(Y1.shape), np.zeros(Y2.shape)
u,v = sistema([Y1, Y2], t, param)
M = (np.hypot(u,v))
u /= M
v /= M
Q = plt.quiver(Y1, Y2, u, v, M, angles='xy')
plt.xlabel('X1')
plt.ylabel('X2')
plt.xlim([0, x_range])
plt.ylim([0, y_range])
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0), useMathText='True')
plt.show()
#####Formato artistico
from matplotlib.pyplot import cm
plt.streamplot(Y1, Y2, u, v,          # data
               color=M,         # array that determines the colour
               cmap=cm.cool,        # colour map
               linewidth=2,         # line thickness
               arrowstyle='->',     # arrow style
               arrowsize=1.5)       # arrow size
plt.colorbar()                      # add colour bar on the right
plt.show()