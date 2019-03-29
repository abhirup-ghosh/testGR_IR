"""
(C) Siddharth Mohite, 2014-12-17
"""

import numpy as np
from scipy import integrate

def snr(h,f,psd):
    return np.sqrt(integrate.simps((np.abs(h)**2)/psd , f, even='avg'))*2
