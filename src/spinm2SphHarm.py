# Calculate spin-weighted Spherical Harmonics for spin weight s = -2
# function = spinm2_sphharm(l,m,theta,phi)
# Ashok, 2013-07-13 

import numpy as np
import math as mt
from math import factorial as fac

#define spin-weighted Spherical Harmonics with spin weight s = -2 as a function of  l,m, theta(ange with z axis) and phi(angle in x-y plain with x axis ) 
def spinm2_sphharm(l,m,theta,phi):
    k1 = max(0,m-2)	
    k2 = min(l+m,l-2)
    a = sum((((-1)**k)*mt.sqrt(fac(l+m)*fac(l-m)*fac(l+2)*fac(l-2))*((np.cos((theta)/2))**(2*l+m-2*k-2))*((np.sin((theta)/2))**(2*k-m+2)))/(fac(k)*fac(k-m+2)*fac(l+m-k)*fac(l-k-2)) for k in range(k1,k2+1))
    ylm = mt.sqrt((2*l+1)/(4*np.pi))*a*np.exp(m*(phi)*1j)
    return ylm

