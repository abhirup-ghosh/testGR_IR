import numpy
import matplotlib.pyplot as plt
import os 

m_min, m_max, a_min, a_max = -250., 250., -2., 2.
M = 100.
Q = 1.
n = 1.
debug_plots = 'y' #'n'
#savefig = '/archive/home/abhirup/Documents/Work/testGR_IR/plots/2015-04-10_IMRPhenomB/mf_af_contours_%s_%s_%s.png'%(str(M), str(Q), str(n))

#os.system('nohup python calc_dE_dJ_posterior.py -i /archive/home/abhirup/Documents/Work/testGR_IR/runs/2015-02-25/inspiral/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat -r /archive/home/abhirup/Documents/Work/testGR_IR/runs/2015-04-10/ringdown/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat -f nospin_pan2011 --m_min %f --m_max %f --a_min %f --a_max %f --savefig %s > delta_%s_%s_%s.out 2> delta_%s_%s_%s.err &'%(str(M), str(Q), str(n), str(M), str(Q), str(n), m_min, m_max, a_min, a_max, savefig, str(M), str(Q), str(n), str(M), str(Q), str(n)))

os.system('nohup python calc_dE_dJ_posterior.py -i /archive/home/abhirup/Documents/Work/testGR_IR/runs/2015-02-25/inspiral/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat -r /archive/home/abhirup/Documents/Work/testGR_IR/runs/2015-04-11/ringdown/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat -f nospin_pan2011 -M %f -Q %f -n %f --debug-plots %s --m_min %f --m_max %f --a_min %f --a_max %f > delta_%s_%s_%s.out 2> delta_%s_%s_%s.err &'%(str(M), str(Q), str(n), str(M), str(Q), str(n), M, Q, n, debug_plots, m_min, m_max, a_min, a_max, str(M), str(Q), str(n), str(M), str(Q), str(n)))

#plt.savefig('2dposteriors_mf_af_contours.png', dpi=150)
