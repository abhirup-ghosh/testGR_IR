import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import os
import commands
import matplotlib.pyplot as plt



in_dir_root = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations_modGR/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2015-04-08/a2_20/imrtestgr_nonprecspin_Healy2014_bins_401_range_-2.0_2.0'
in_dir_list = commands.getoutput('ls -d %s/injection_*'%in_dir_root).split()
out_dir = '%s/2d_plots_P_dMfbyMf_dchifbychif'%in_dir_root
os.system('mkdir -p %s'%out_dir)

for in_dir in in_dir_list:
    P = np.loadtxt('%s/data/P_dMfbyMf_dchifbychif.dat.gz'%in_dir)
    dMfbyMf_vec = np.loadtxt('%s/data/dMfbyMf_vec.dat.gz'%in_dir)
    dchifbychif_vec = np.loadtxt('%s/data/dchifbychif_vec.dat.gz'%in_dir)
    plt.figure()
    plt.imshow(P)	
    #plt.pcolormesh(dMfbyMf_vec, dchifbychif_vec, P)	
    plt.savefig('%s/%s.png'%(out_dir, in_dir[-2:]))
    plt.close()
