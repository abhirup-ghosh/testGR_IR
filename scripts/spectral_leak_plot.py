import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import template_lalsiminspiral as phhsi_lal
import pycbc
from pycbc import detector

Mc, q, dL, iota, t0, phi0, ra, sin_dec, pol, f_low, dt  = 28.10, 0.81, 500.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 1./2048.

t, hpt, hct = phhsi_lal.lalsiminspiraltd_waveform_SI(Mc,q,dL,iota,t0,(phi0 %(2.*pi)),f_low,dt)

Fp,Fc = detector.overhead_antenna_pattern(ra, np.arcsin(sin_dec), pol)
hoft = Fp*hpt + Fc*hct

phioft = np.unwrap(np.arctan2(hct, hpt))
Foft = np.gradient(phioft)/np.gradient(t)/(2*np.pi)

fig = plt.figure(figsize=(9,9))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

ax1.plot(t, hoft)
plt.savefig('./chapter_5_fig1_motivation.png')
