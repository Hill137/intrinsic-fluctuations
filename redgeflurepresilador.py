'Redes geneticas estocasticas' 'togle switch'
#este programa permite cuantificar las fluctuaciones durante
#la evolucion temporal del sistema

#importamos libreria que sera usada
import statistics as st

def varios(x,k,j,q): #p= pasos  j= tamano del sistema, 
   # q= numero de veces que se realiza la misma simulacion 
    l1=0
    while l1<q:
        if l1==0:
            x=[[Omega,0], #vector de condiciones iniciales [ARm0, P0]
               [0,0],                              #  [ARm1, P1]
                [0,0]]
            T0, l2 =ev(x,k,j)
            ky=len(T0)
            l4= np.zeros([n,n2,ky,q])
            l3=np.zeros([n,n2,ky])
            l5=np.zeros([n,n2,ky])
            l6=np.zeros([n,n2,ky])
            l4[:,:,:,l1]=l2
            l1=l1+1
        else:
            x=[[Omega,0], #vector de condiciones iniciales [ARm0, P0]
               [0,0],                              #  [ARm1, P1]
                [0,0]]
            T0, l2 =ev(x,k,j)
            l4[:,:,:,l1]=l2
            l1=l1+1
    for i in range(n):
        for j1 in range(n2):
            for j2 in range(ky):
                l3[i,j1,j2]=st.mean(l4[i,j1,j2,:])
                l5[i,j1,j2]=math.sqrt(st.pvariance(l4[i,j1,j2,:],l3[i,j1,j2]))
                l6[i][j1][j2]= heav(l3[i][j1][j2] - 2*l5[i][j1][j2])
    l7=0
    for i in range(n):
        for j1 in range(n2):
            plt.figure(nom1[l7])
            plt.plot(T0,l3[i][j1]/j,color="green", label= nom3[l7])
            plt.plot(T0,(l3[i][j1]+2*l5[i][j1])/j,'--', color="red", label= '+2*sigma')
            plt.plot(T0,l6[i][j1]/j,':', color="red", label= '-2*sigma')
            plt.xlabel("t (tiempo)")
            plt.ylabel("x (concentracion de ARm)")
            plt.title(nom2[l7] )
            plt.xlim(0, tf)
            plt.ylim(0)
            plt.legend()
            l7= l7 + 1 
    plt.show()

#funciones para poder graficar
def heav(vec): #funcion que devuelve solo valores mayores que cero
    if vec>=0:
        vec=vec                 
    else:
        vec=0  
    return vec

#nombres que usaremos para las graficas y nombrar las variables
nom1=[1,2,3,4]

nom2=['Promedio del numero de moleculas Arm_0', 'Promedio del numero de moleculas P_0',
      'Promedio del numero de moleculas Arm_1', 'Promedio del numero de moleculas P_1']

nom3=['Arm_0', 'P_0', 'Arm_1', 'P_1' ]

#realizar las graficas de las variables
varios(x,k,Omega,10000)