from scipy.integrate import quad
from scipy.integrate import trapz
import numpy as np
import matplotlib.gridspec as gs
from random import seed
from random import random
from copy import deepcopy
from matplotlib.patches import Ellipse
import time
import math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
#constants
OmegaM=0.3087
lamb=0.2
#covariance matrix 
npar=2
FInv = np.zeros([npar,npar])
FInv[0,0]=0.1
FInv[1,1]=0.6
FInv[0,1]=0.2
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
ax.set_xlim(-2.2, 2.2)
ax.set_ylim(-2.2, 2.2)
#ax.plot(OmegaM,lamb,"o")
plt.show()