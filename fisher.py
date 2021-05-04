from scipy.integrate import quad
from scipy.integrate import trapz
import numpy as np
import matplotlib.gridspec as gs
from random import seed
from random import random
from copy import deepcopy

import time
import math
# Gráficos
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d



# Velocidade da Luz
c = 3e5 # (Km/s)

#Constante de Hubble
h = 0.7
Ho = 100 * h # (km/s/Mpc)

# Magnitude absoluta das Super Novas (banda B)
M = -19.5

#lambdas
lambv = [0.0,0.1,0.7]
lamb=0.3
#betas
beta=[0.05]
#omega 
OmegaM=0.296

def dl(zv,OmegaM,lamb):
    [integral, err] = quad(lambda x: 1/np.sqrt(OmegaM*(3/(3-lamb**2)*(1+x)**3 + ((1-OmegaM)/OmegaM + lamb**2/(lamb**2-3))*(1+x)**(lamb**2))), 0, zv)

    dist = c/Ho*(1+zv)*integral
    return dist

## moch data z Dl inc_Dl
with open("moch.txt", "r") as f:
    f.readline()

    lista_distancias = []
    lista_incerteza = []
    lista_z = []

    for line in f:
        line = line.split(" ")
        lista_distancias.append(float(line[1]))
        lista_incerteza.append(float(line[2]))
        lista_z.append(float(line[0]))

#plt.figure(1)
#plt.plot(lista_z,lista_distancias,".")
#plt.xlabel("z")
#plt.ylabel("Distância Luminosa (Mpc)")
#plt.show()
#fisher matrix

h=0.001
npar=2
F = np.zeros([npar,npar])
for i in range(npar):
        for j in range(npar):
            if i==j and i==0:
                for k in range(len(lista_incerteza)):
                    #F[i,j] += (lista_incerteza[k]**-2)*((dl(lista_z[k],OmegaM+h,lamb)-2*dl(lista_z[k],OmegaM,lamb)+dl(lista_z[k],OmegaM-h,lamb))/h**2)
                    F[i,j]+=(lista_incerteza[k]**-2)*((dl(lista_z[k],OmegaM+h,lamb)-dl(lista_z[k],OmegaM-h,lamb))/(2*h))**2
                    
            if i!=j:
                #F[i,j]+= (lista_incerteza[k]**-2)* ((dl(lista_z[k],OmegaM+h,lamb+h)-dl(lista_z[k],OmegaM+h,lamb-h)-dl(lista_z[k],OmegaM-h,lamb+h)-dl(lista_z[k],OmegaM-h,lamb-h))/4*h**2)
                F[i,j]+=(lista_incerteza[k]**-2)*((dl(lista_z[k],OmegaM,lamb+h)-dl(lista_z[k],OmegaM,lamb-h))/(2*h))*((dl(lista_z[k],OmegaM+h,lamb)-dl(lista_z[k],OmegaM-h,lamb))/(2*h))
            if i==j and i==1:
                for k in range(len(lista_incerteza)):
                    #F[i,j] += (lista_incerteza[k]**-2)*((dl(lista_z[k],OmegaM,lamb+h)-2*dl(lista_z[k],OmegaM,lamb)+dl(lista_z[k],OmegaM,lamb-h))/h**2) 
                    F[i,j]+=(lista_incerteza[k]**-2)*((dl(lista_z[k],OmegaM,lamb+h)-dl(lista_z[k],OmegaM,lamb-h))/(2*h))**2

k=5
print((lista_incerteza[k]**-2)*((dl(lista_z[k],OmegaM+h,lamb)-dl(lista_z[k],OmegaM-h,lamb))/(2*h))**2,"\n",(lista_incerteza[k]**-2)*((dl(lista_z[k],OmegaM,lamb+h)-dl(lista_z[k],OmegaM,lamb-h))/(2*h))**2)                   
print("FISHER\n", F)                  
FInv=np.mat(F).I
print("covariância\n",FInv)
##ploting the elipse
from matplotlib.patches import Ellipse
a=((FInv[0,0]+FInv[1,1])/2)+np.sqrt(((FInv[0,0]+FInv[1,1]))**2/4+FInv[0,1]**2)

b=((FInv[0,0]+FInv[1,1])/2) - np.sqrt(((FInv[0,0]-FInv[1,1]))**2/4+FInv[0,1]**2)
print("a=",a,"\n b=",b)
angle=(np.arctan(2*FInv[0,1]/(FInv[0,0]-FInv[1,1]))/2)*180/np.pi
print("angulo",((np.arctan(2*FInv[0,1]/(FInv[0,0]-FInv[1,1]))/2)*180/np.pi))
fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
ellipse = Ellipse((OmegaM, lamb), np.sqrt(a)*2.48,np.sqrt(b)*2.48, angle=angle, alpha=0.5)
ax.add_artist(ellipse)
ellipse = Ellipse((OmegaM, lamb), np.sqrt(a)*1.52,np.sqrt(b)*1.52,angle=angle, alpha=0.9)
ax.add_artist(ellipse)
ax.set_xlim(0.2, 0.4)
ax.set_ylim(0.275, 0.330)
#ax.plot(OmegaM,lamb,"o")
plt.show()
