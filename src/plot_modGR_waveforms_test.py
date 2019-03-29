# Script to plot the modGR waveforms including the higher harmonics
# NKJ-M, 12.2015

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description = 'Plot modGR waveforms to test that they look OK')
parser.add_argument("datfile", help = "waveform data file")
parser.add_argument("tag", help = "tag for output files")
args = parser.parse_args()

data = np.loadtxt(args.datfile, skiprows=1, unpack=True)

# Unpack time and define modulus

t = data[0]

h = (data[1]*data[1] + data[2]*data[2])**0.5

# Plot the modulus and real part of the waveform

plt.plot(t, h, label='$|h_{22}|$',color='k')
plt.plot(t, data[1], label='Re $h_{22}$')
plt.xlabel('$t$')
plt.savefig('%s_plot.png' %(args.tag), dpi=200)
plt.close()

# Plot the same things just for the end of the waveform

# First find the maximum of h

idx_max = np.argmax(h)
idx_zoom = np.arange(2*idx_max-len(t),len(t))

plt.plot(t[idx_zoom], h[idx_zoom], label='$|h_{22}|$',color='k')
plt.plot(t[idx_zoom], data[1][idx_zoom], label='Re $h_{22}$')
plt.xlabel('$t$')
plt.savefig('%s_plot_zoom.png' %(args.tag), dpi=200)
plt.close()
