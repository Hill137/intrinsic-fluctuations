'funciones de Hill'

#programa para graficar las difernetes funciones de Hill

# se importan las librerias
import numpy as np
import math 
import random as rnd
from matplotlib import pyplot

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

def exponencial(y,n):
    if n==0:
        return 1
    else:
        return y**n

def series(x,K,n):
    l1=0
    l2=0
    while l1<n+1:
       l3= exponencial(K,n-l1)* exponencial(x,l1)
       l2=l2+l3
       l1=l1+1
    return l2

def seriei(x,K,n):
    l1=0
    l2=0
    while l1<n+1:
       l3= exponencial(K,n-l1)* exponencial(x,l1) * prod(n,l1)/fac(l1)
       l2=l2+l3
       l1=l1+1
    return l2

def H_d(x,n,K,j): #funcion de hill determinista
    return (x/j)**n/(K + (x/j)**n)

def H_sd(x,n,K,j):# funcion de Hill semideterminista
    if x<n:
        return 0
    if x>=n and x<j/2:
        return prod(x,n)/(K*(j**n) + prod(x,n))
    if x>=n and x>=25:
        return prod(x+n,n)/(K*(j**n) + prod(x+n,n))

def H_ss(x,n,K,j):# funcion de Hill semiestocastica
    if x<n:
        return 0
    if x>=n and x<j/2:
        return prod(x,n)/(K*(j**n)*np.exp((n-1)/j) + prod(x,n))
    if x>=n and x>=j/2:
        return prod(x+n,n)/(K*(j**n)*np.exp((n-1)/j) + prod(x+n,n))

def H_ap(x,n,K,j): # funcion de hill con aproximaciones
    return ((x/j)**n + (n/(2*j))*(n-1)*(x/j)**(n-1) )/(K + (x/j)**n + (n/(2*j))*(n-1)*(x/j)**(n-1))

def H_ds(x,n,K,j): #funcion de hill determinista secuencial
    return (x/j)**n/(series((x/j),K,n))

def H_ds1(x,n,K,j): #funcion de hill determinista secuencial exacto
    return (x/j-K)*(x/j)**n/((x/j)**(n+1)-(K)**(n+1))

def H_di(x,n,K,j): #funcion de hill determinista independiente
    return (x/j)**n/(K+x/j)**n

def H_di1(x,n,K1,K2,K3,j): #funcion de hill determinista independiente
    return (x/j)**n/(1 + K1*x + K1*K2*x**2 + K1*K2*K3*x**3 + x**n )

# Valores del eje X que toma el gráfico.
x = np.arange(0.0000001, 200, 0.0001)
#x=range(0,2000)

y= [math.log(i) for i in x ]
# Graficar ambas funciones.
#pyplot.plot(x, [H_di(i,1,1,1) for i in x], label='n=1')
pyplot.plot(y, [H_ds(i,40,1,1) for i in x], label='n=3.2')
pyplot.plot(y, [H_ds(i,4,1,1) for i in x], label='n=4')
pyplot.plot(y, [H_ds(i,1,1,1) for i in x], label='n=1.7')

# Establecer el color de los ejes.
pyplot.axhline(0, color="black")
pyplot.axvline(0, color="black")

# Limitar los valores de los ejes.
pyplot.xlim(0, 20)
pyplot.ylim(0, 1.05)
pyplot.title("Hill secuencial")
pyplot.xlabel("e (concentración de enzimas)")
pyplot.ylabel("H")
pyplot.legend()

# Guardar gráfico como imágen PNG.
pyplot.savefig("han10.png")

# Mostrarlo.
pyplot.show()