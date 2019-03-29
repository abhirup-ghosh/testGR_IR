"""
This module generates plot tickmarks for angular variables
"""
import numpy as np
from numpy import pi
from fractions import Fraction

def ticks(mini, maxi, parts):
  return pi*np.linspace(float(mini), float(maxi), int(parts)+1)

def ticklabels(typ, mini, maxi, parts):
  if typ=='h':
    return [r'%dh'%(ang) for ang in 12*np.linspace(float(mini), float(maxi), int(parts)+1)]
  elif typ=='d':
    return [r'$%d^\circ$'%(ang) for ang in 180*np.linspace(float(mini), float(maxi), int(parts)+1)]
  else:
    ra=[]
    for ii in np.linspace(float(mini), float(maxi), int(parts)+1):
      if int(ii)==ii:
        if ii==0:
          ra.append(r'0')
        elif ii==1:
          ra.append(r'$\pi$')
        elif ii==-1:
          ra.append(r'$-\pi$')
        else:
          ra.append(r'$%d\pi$'%(ii))
      else:
          num = Fraction(ii).limit_denominator(10).numerator
          den = Fraction(ii).limit_denominator(10).denominator
          if num==1:
            ra.append(r'$\frac{\pi}{%d}$'%(den))
          elif num==-1:
            ra.append(r'$-\frac{\pi}{%d}$'%(den))
          else:
            if num>0:
              ra.append(r'$\frac{%d\pi}{%d}$'%(num, den))
            else:
              ra.append(r'$-\frac{%d\pi}{%d}$'%(-num, den))
    return ra
