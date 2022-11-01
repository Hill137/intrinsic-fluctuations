from pylab import *
import numpy as np
import matplotlib.pyplot as plt

a =  1
b = 1/100. 
c = 1.01
d = 1 / 1000. 
Dt = 0.001

def initialize():
    global x, xresult, y, yresult, t, timesteps
    x = 1
    y = 1
    xresult = [x]
    yresult = [y]
    t = 0.
    timesteps = [t]
    
def observe():
    global x, xresult, y, yresult, t, timesteps
    xresult.append(x)
    yresult.append(y)
    timesteps.append(t)

alpha = .183
f =13
n = 8

def update():
    global x, xresult, y, yresult, t, timesteps
    nextx = x +  ( alpha*(1. + f * y**n)/(1.+ y**n) - x) * Dt
    nexty = y + ( alpha*(1.+ f * x**n)/(1.+ x**n) - y) * Dt
    x, y = nextx, nexty
    t = t + Dt

initialize()
while t < 30:
    update()
    observe()

plt.plot(timesteps, xresult, 'b')
plt.plot(timesteps, yresult, 'g:')
plt.xlabel("t")
plt.ylabel("x")
plt.title("Togle swich")
plt.show()