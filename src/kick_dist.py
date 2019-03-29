# Calculate the kick distribution for a uniform distribution of masses and spins

import numpy as np
import matplotlib.pyplot as plt
import nr_fits as nrf

# Set the mass range and the number of samples

mlow = 5.
mhigh = 50.
Nsamp = 1.e3

# Calculate the distribution of masses and spins

m1 = np.random.uniform(mlow,mhigh,Nsamp)
m2 = np.random.uniform(mlow,mhigh,Nsamp)

q = m1/m2

ctheta1 = np.random.uniform(-1.,1.,Nsamp)
stheta1 = (1 - ctheta1**2.)**0.5
phi1 = np.random.uniform(0,2*np.pi,Nsamp)
chimag1 = np.random.uniform(0.,1.,Nsamp)

chi1para = chimag1*ctheta1
chi1perpx = chimag1*stheta1*np.cos(phi1)
chi1perpy = chimag1*stheta1*np.sin(phi1)

ctheta2 = np.random.uniform(-1.,1.,Nsamp)
stheta2 = (1 - ctheta2**2.)**0.5
phi2 = np.random.uniform(0.,2*np.pi,Nsamp)
chimag2 = np.random.uniform(0.,1.,Nsamp)

chi2para = chimag2*ctheta2
chi2perpx = chimag2*stheta2*np.cos(phi2)
chi2perpy = chimag2*stheta2*np.sin(phi2)

# Compute the angle-averaged kicks

kicks = nrf.bbh_recoil_simple_averaged_over_angles(q, chi1perpx, chi1perpy, chi1para, chi2perpx, chi2perpy, chi2para)

# Plot a histogram

plt.hist(kicks,bins=10)
plt.xlabel("Kick (km/s)")
plt.ylabel("Frequency")
plt.show()
