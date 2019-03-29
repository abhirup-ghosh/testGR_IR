import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
import loadline as ll
import matplotlib
import calc_dE_dJ_debug_plots as dp

date = '2015-05-06'
approximant = 'IMRPhenomB'

phase_list = ['inspiral', 'ringdown']
M_list = [11.,22.]
Q_list = [10.]
n_isco_list = [1.0]



# Loading Data

for phase in phase_list:
  for M in M_list:
    for Q in Q_list:
        for n_isco in n_isco_list:

		M1 = M/(1.+Q)
                M2 = M*Q/(1.+Q)
                ETA = (M1*M2)/((M1+M2)**2.)

                MF = M*(1. + (np.sqrt(8./9.)-1.)*ETA - 0.4333*(ETA**2.) - 0.4392*(ETA**3.))
                AF = ETA*np.sqrt(12.) - 3.871*(ETA**2.) + 4.028*(ETA**3.)


                data_file = '../runs/%s/%s/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat'%(date, phase, str(M), str(Q), str(n_isco))
                legend = ll.load_line(data_file, dtype='string')

                logl, m1, m2= np.loadtxt(data_file, skiprows =1, usecols = (list(legend).index('logl'), list(legend).index('m1'), list(legend).index('m2')), unpack=True)

		x1, x2 = dp.debug_plots(m1, m2, logl, phase, M, Q, n_isco)
