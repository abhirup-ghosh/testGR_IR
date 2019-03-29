from __future__ import division
from math import ceil, floor, log10
import sys
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib
import matplotlib.ticker as ticker


matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

matplotlib.rcParams['xtick.labelsize'] = 15
matplotlib.rcParams['ytick.labelsize'] = 15


fname = "final-nresults-5e4.dat"

d = np.loadtxt(fname)

match = d[:,0]

inj_mass1 = d[:,9]
inj_mass2 = d[:,10]
inj_spin1z = d[:,13]
inj_spin2z = d[:,16]


mmin = 0.96#min(match)




fig = Figure()
ax = fig.add_subplot(111)
coll = ax.hexbin(inj_mass1, inj_mass2, C=match, cmap=plt.cm.jet, gridsize=60, reduce_C_function = np.min, vmin=mmin, vmax=1.0)
ax.set_xlabel(r"$\boldsymbol{\mathrm{m_{1}}} (\mathrm{M_\odot})$", fontsize=25)
ax.set_ylabel(r"$\boldsymbol{\mathrm{m_{2}}} (\mathrm{M_\odot})$", fontsize=25)
ax.set_xlim([min(inj_mass1), max(inj_mass1)])
ax.set_ylim([min(inj_mass2), max(inj_mass2)])
fig.colorbar(coll, ax=ax).set_label("Minimum  Fitting Factor", fontsize=18)
ax.grid(True)
canvas = FigureCanvas(fig)
fig.savefig("mass1-mass2.png", dpi=300)



