#!/usr/bin/env python

import os, numpy as np
#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
#import imrtestgr as tgr import tgrplotsettings 

runtag = 'Run2_l2m2'

noise_list = [2, 3]+range(5, 21)
approx = {'S': 'SEOBNRv2_ROM_DoubleSpin', 'I': 'IMRPhenomPv2'}

loc_list = ['/home/archis/Work/WaveformSystematics/imrtgr_runs/2015-12-17/SXS-Ossokine/TGR_frames/%s_noise%d/%sthreePointFivePN_vanilla_seglen8/imrtgr/data'%(runtag, noise, approx['S']) for noise in noise_list]

gr_conf_list = np.array([])

for loc in loc_list:
  gr_conf = np.loadtxt(os.path.join(loc, 'GR_confidence.txt'))
  gr_conf_list = np.append(gr_conf_list, gr_conf)


pp_x = np.linspace(0., 1., 500)
pp_y = np.array([len(np.where(gr_conf_list<p_x)[0]) for p_x in pp_x])/float(len(gr_conf_list))

plt.figure(figsize=(8,8))
plt.plot(pp_x, pp_y, 'b-')
plt.plot(pp_x, pp_x, 'r--')
plt.xlabel('GR confidence level')
plt.ylabel('Fraction of events with GR confidence below a given GR confidence level')

plt.savefig('GR_conf_ppplot_%s_%s.png'%({'Run1': 'Precessing', 'Run2': 'Aligned'}[runtag[:4]], runtag[5:]))
