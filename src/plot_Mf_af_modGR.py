# Script to plot final masses and spins for a variety of different mass ratios and GR modifications (all nonspinning), as well as the allowed region for binaries with the HLZ fit (currently just aligned spins and an approximation to the lower bound)
# NKJ-M, 11.2016, based on earlier code

import numpy as np
import matplotlib.pyplot as plt
import nrutils as nr

from scipy.optimize import minimize

# Setup

npts = 100

massfit = "UIB2016" #"HLZ2014"
spinfit = massfit

# Give data

# modGR values for various a2 values for q = 1 and 2

Mfchifq1a25 = [0.944, 0.67]
Mfchifq1a210 = [0.938, 0.651]
Mfchifq1a220 = [0.928, 0.621]
Mfchifq1a2100 = [0.896, 0.495]
Mfchifq1a2200 = [0.878, 0.410]
Mfchifq1a2300 = [0.866, 0.351]
Mfchifq1a2400 = [0.857, 0.306]

# modGR values for a2 = 20 and various qs (rounded to a single decimal place in the names)

Mfchifq1p5a220 = [0.934, 0.6]
Mfchifq2a220 = [0.943, 0.559]
Mfchifq2p5a220 = [0.95, 0.518]
Mfchifq3a220 = [0.957, 0.482]
Mfchifq3p6a220 = [0.964, 0.446]
Mfchifq4p1a220 = [0.968, 0.417]
Mfchifq4p7a220 = [0.972, 0.393]
Mfchifq5p4a220 = [0.976, 0.363]
Mfchifq6a220 = [0.979, 0.34]
Mfchifq6p8a220 = [0.982, 0.316]

# Define functions for minimization and maximization

def afmin(x):
    return nr.bbh_final_spin_projected_spins(1.,x[0],1.,1.,x[1]*np.pi,x[2]*np.pi, spinfit)

def afmax(x):
    return -nr.bbh_final_spin_projected_spins(1.,x[0],1.,1.,x[1]*np.pi,x[2]*np.pi, spinfit)

# Plot the allowed range of final mass and spin for the fit

Mfmin = nr.bbh_final_mass_projected_spins(1.,1.,1.,1.,0.,0.,massfit)/2.
step = (1. - Mfmin)/npts
qlbmin = np.zeros(npts)
tilt1lbmin = np.zeros(npts)
tilt2lbmin = np.zeros(npts)
aflbmin = np.zeros(npts)
qubmin = np.zeros(npts)
tilt1ubmin = np.zeros(npts)
tilt2ubmin = np.zeros(npts)
afubmin = np.zeros(npts)
Mflubmin = np.linspace(Mfmin,1.-1e-4,npts)
x0 = [0.5, 0.25, 0.25]
for k in range(npts):
    cons = ({'type': 'eq', 'fun': lambda x:  nr.bbh_final_mass_projected_spins(1.,x[0],1.,1.,x[1]*np.pi,x[2]*np.pi,massfit)/(1. + x[0]) - Mflubmin[k]})
    res = minimize(afmin, x0, bounds=((1e-7,1.), (0.,1.), (0.,1.)), constraints=cons, tol=1e-6)
    qlbmin[k] = res.x[0]
    tilt1lbmin[k] = res.x[1]*np.pi
    tilt2lbmin[k] = res.x[2]*np.pi
    aflbmin[k] = afmin(res.x)
    res = minimize(afmax, x0, bounds=((1e-7,1.), (0.,1.), (0.,1.)), constraints=cons, tol=1e-6)
    qubmin[k] = res.x[0]
    tilt1ubmin[k] = res.x[1]*np.pi
    tilt2ubmin[k] = res.x[2]*np.pi
    afubmin[k] = afmin(res.x)
    #x0 = res.x # Doesn't work well, at least with default method                                                                                                                 
    #print("q, tilt1, tilt2, Mf, af: %f %f %f %f %f"%(qlbmin[k], tilt1lbmin[k], tilt2lbmin[k], Mflbmin[k], aflbmin[k]))                                                           
plt.fill_between(Mflubmin, aflbmin, afubmin, alpha=0.5)

# Concatenating values

Mfchifq1a2 = [Mfchifq1a25, Mfchifq1a210, Mfchifq1a220, Mfchifq1a2100, Mfchifq1a2200, Mfchifq1a2300, Mfchifq1a2400]

Mfchifa220 = [Mfchifq1a220, Mfchifq1p5a220, Mfchifq2a220, Mfchifq2p5a220, Mfchifq3a220, Mfchifq3p6a220, Mfchifq4p1a220, Mfchifq4p7a220, Mfchifq5p4a220, Mfchifq6a220, Mfchifq6p8a220]

# Plotting

plt.scatter(*zip(*Mfchifq1a2), color='red', label = '$q = 1$')
plt.scatter(*zip(*Mfchifa220), color='orange', label = '$a_2 = 20$')
plt.xlabel('$M_f$')
plt.ylabel('$\chi_f$')
#plt.xlim(0.88, 1.)
plt.xlim(0.85,1.)
plt.ylim(-1., 1.)
plt.legend(loc='lower left')
plt.show()
