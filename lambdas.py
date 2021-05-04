from scipy.integrate import quad
from scipy.integrate import trapz
import numpy as np

from random import seed
from random import random

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
lambv = [0.1,0.3,0.4]

def Einv(x,OmegaM,lamb):
    return 1/np.sqrt(OmegaM*(3/(3-lamb**2)*(1+x)**3 + ((1-OmegaM)/OmegaM + lamb**2/(lamb**2-3))*(1+x)**(lamb**2)))

def dl_ot(zv,OmegaM,lamb):
    distancias=np.zeros(len(zv)-1)
    f1 = Einv(zv[0],OmegaM,lamb)
    soma = 0
    for i in range(1,len(zv)):
        f0 = f1
        f1 = Einv(zv[i],OmegaM,lamb)
        integral = 0.5*(f0+f1)*(zv[i]-zv[i-1])
        soma+=integral
        dist = (c*(1+zv[i])/Ho)*soma
        distancias[i-1] = dist
    return distancias

def dl(zv,OmegaM,lamb):
    [integral, err] = quad(lambda x: 1/np.sqrt(OmegaM*(3/(3-lamb**2)*(1+x)**3 + ((1-OmegaM)/OmegaM + lamb**2/(lamb**2-3))*(1+x)**(lamb**2))), 0, zv)

    dist = c/Ho*(1+zv)*integral
    return dist
def H(x,OmegaM,lamb):
    return Ho*np.sqrt(OmegaM*(3/(3-lamb**2)*(1+x)**3 + ((1-OmegaM)/OmegaM + lamb**2/(lamb**2-3))*(1+x)**(lamb**2)))
#def H(zv,OmegaM,lamb):
 #   return Ho*np.sqrt(OmegaM*((1+zv)**3)+1-OmegaM)

def fz(zv,OmegaM,lamb):
    if (zv<=1):
        rz=(1+2*zv)
    else:
        rz=(15-3*zv)/4
        
    return (4*np.pi*rz*dl(zv,OmegaM,lamb)**2)/(H(zv,OmegaM,lamb)*(1+zv)**3)

def inc(z,OmegaM,lamb):
    return ((0.1449*z-0.0118*z**2+0.0012*z**3)**2+(0.05*z)**2)**0.5*dl(z,OmegaM,lamb)



fz1=[]
fz2=[]
fz3=[]


zv = np.linspace(0.01,2,500)
for i in (zv):
    fz1.append(fz(i,0.3087,lambv[0]))
    fz2.append(fz(i,0.3087,lambv[1]))
    fz3.append(fz(i,0.3087,lambv[2]))
    
#plt.plot(zv,fz_qe2,'r.',zv,fz_qe3,'b.')
#plt.xlabel('z',size=16)
#plt.ylabel('f(z)',size=16)

#ficheiros para as moch

f= open("moch.txt","w+")

#gerar valores de redshift
z1=[]
z2=[]
z3=[]
seed(3)
i=0;
while (i<1000):
    x=random()*2
    y=random()*max(fz1)
    if (y<=fz(x,0.3087,lambv[0])):
        z1.append(x)
        i+=1
i=0;
while (i<1000):
    x=random()*2
    y=random()*max(fz2)
    if (y<=fz(x,0.3087,lambv[1])):
        z2.append(x)
        i+=1
i=0;
while (i<1000):
    x=random()*2
    y=random()*max(fz3)
    if (y<=fz(x,0.3087,lambv[2])):
        z3.append(x)
        i+=1



zv = np.linspace(0.01, 2, 500)




lista_z=z1

lista_z.sort()

zl=[]
for i in range(len(lista_z)):
    zl.append(lista_z[i])

zl.append(0)
zl.sort()
lista_incerteza=[]
lista_distancias=[]
for i in (lista_z):
    lista_incerteza.append(inc(i,0.3087,lambv[0])/1)
    lista_distancias.append(dl(i,0.3087,lambv[0]))    
for i in range(len(lista_distancias)):
    lista_distancias[i]=np.random.normal(lista_distancias[i],lista_incerteza[i])
    f.write(str(lista_z[i])+" "+str(lista_distancias[i]) +" "+ str(lista_incerteza[i])+"\n")
print("all done")

plt.plot(lista_z,lista_distancias,".")
plt.xlabel("z")
plt.ylabel("Distância Luminosa (Mpc)")
#plt.show()
plt.close()
f.close()
"""
#qui quadrado
lambdaj=np.linspace(0,2,50)
omegaj=np.linspace(0.1,1,50)

qui2=[[] for i in range(len(lambdaj))]
quiq=[]


for i in range(len(lambdaj)):
    for j in range(len(omegaj)):
        soma=0
        dist=dl_ot(zl,omegaj[j],lambdaj[i])
        for k in range(len(lista_z)):
            soma+=(dist[k]-lista_distancias[k])**2/(lista_incerteza[k]**2)
        qui2[i].append(soma)
        quiq.append(soma)
    


minimo=qui2[0][0]
omeg=0
lam=0
for i in range(len(lambdaj)):
    for j in range(len(omegaj)):
        if qui2[i][j]<minimo:
            minimo=qui2[i][j]
            lam=i
            omeg=j
print('lambda=',lambdaj[lam],'\nomega=',omegaj[omeg],'\nMínimo=',minimo)
print('Mínimo do Python',min(quiq))


plt.figure(7)
plt.xlabel(r'$\Omega_m$', size=14)
plt.ylabel(r'$\lambda$', size=14)
plt.contour(omegaj,lambdaj,qui2,[minimo+2.13, minimo+6.18])
plt.plot(omegaj[omeg],lambdaj[lam],".")
new=[]
plt.figure(5)
for i in qui2:
    new.append(i[25])
plt.plot(lambdaj,new)
plt.axvline(x=0.0, color='k', linestyle='--')
plt.figure(6)
plt.plot(omegaj,qui2[0])
plt.axvline(x=0.3087, color='k', linestyle='--')
plt.show()
"""