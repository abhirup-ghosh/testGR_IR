#!/usr/bin/python

import commands
import numpy as np
import my_plotter as mp
import matplotlib.pyplot as plt


date_I = '2015-02-07'
date_R = '2015-02-08'

M = [100.] #[30., 50., 100.]
Q = [4.] #[1., 2., 4.]
n_isco_I = [0.5] #[2.0, 2.5, 3.0, 3.5, 4.0, 4.3]
n_isco_R = [2.0]

i,j,k=0,0,0

pos_file_I = '../runs/%s/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat'%(date_I, str(M[i]), str(Q[j]), str(n_isco_I[k]))
pos_file_R = '../runs/%s/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat'%(date_R, str(M[i]), str(Q[j]), str(n_isco_R[k]))


def load_line(file_name, n=1, **kwargs):
  """
  Reads the first n lines from a file. Default n=1. If n<0, reads the last n lines.
  """
  num_lines = int(commands.getoutput('cat %s | wc -l'%file_name))+int(commands.getoutput('tail -c 1 %s'%file_name)!='')
  if n<0:
    n = num_lines+1+n
  return np.genfromtxt(file_name, skip_header=n-1, skip_footer=num_lines-n, **kwargs)

legend_I = load_line(pos_file_I, dtype='string')
legend_R = load_line(pos_file_R, dtype='string')

(mc_pos_I, q_pos_I) = np.loadtxt(pos_file_I, usecols=(list(legend_I).index('mc'), list(legend_I).index('q')), skiprows=1, unpack=True)
(mc_pos_R, q_pos_R) = np.loadtxt(pos_file_R, usecols=(list(legend_R).index('mc'), list(legend_R).index('q')), skiprows=1, unpack=True)

plt.figure('Mass estimate', figsize=(6,6))
plt.title('Mass estimate')
plt.subplot2grid((3,3),(0,0),colspan=2)
(mc_av_I, mc_less_I, mc_more_I) = mp.plot1d(mc_pos_I, instance=plt, histtype='stepfilled', bins=20, color='b', limits=True, significance=True, slinestyle='--')
(mc_av_R, mc_less_R, mc_more_R) = mp.plot1d(mc_pos_R, instance=plt, histtype='stepfilled', bins=20, color='b', limits=True, significance=True, slinestyle='--')
#mp.plttext(.2, .85, r'$\mathcal{M}_c$= $%.2f_{\,-%.2f}^{+%.2f}$ $km\,s^{-1}Mpc^{-1}$.'%(mc_av, mc_av-mc_less, mc_more-mc_av), instance=plt)
#plt.xticks([])
#plt.yticks([])
plt.subplot2grid((3,3),(1,2),rowspan=2)
#(q_av, q_less, q_more) = mp.plot1d(q_pos, instance=plt, histtype='stepfilled', bins=20, color='b', limits=True, significance=True, slinestyle='--', orientation='horizontal')
#mp.plttext(.8, .6, r'$q$= $%.2f_{\,-%.2f}^{+%.2f}$.'%(q_av, q_av-q_less, q_more-q_av), instance=plt, rotation=-90)
#plt.xticks([])
#plt.yticks([])
plt.subplot2grid((3,3),(1,0),rowspan=2,colspan=2)
mp.plot2d(mc_pos_I, q_pos_I, instance=plt, histtype='hist2d', bins=100, cmap='Blues', significance=True, sbins=20, slinestyle='--')
mp.plot2d(mc_pos_R, q_pos_R, instance=plt, histtype='hist2d', bins=100, cmap='Blues', significance=True, sbins=20, slinestyle='--')
plt.xlabel(r'$\mathcal{M}_c$ ($M_\odot$)')
plt.ylabel(r'$q$')
#ax = plt.subplot2grid((3,3),(0,2))
#plt.text(0., 0.9, r'Plot of ..')
#plt.xticks([])
#plt.yticks([])
#ax.patch.set_visible(False); ax.axis('off')
plt.tight_layout()
#plt.savefig('../plots/pos.png')
plt.show()
