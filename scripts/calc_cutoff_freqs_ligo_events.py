import numpy as np
import nr_fits as nr
import lal
from scipy import stats


def maxP(logl, logprior):
        """
        Find the sample with maximum a posteriori probability. Returns value
        of posterior and index of sample .
        """
        max_i=0
        max_pos=logl[0]+logprior[0]
        for i in range(len(logl)):
            if logl[i]+logprior[i] > max_pos:
                max_pos=logl[i]+logprior[i]
                max_i=i

        return max_pos,max_i


event_list = ['GW150914', 'GW170104', 'GW170729','GW170809','GW170814', 'GW170818', 'GW170823']
fcut_map = {'GW150914':132, 'GW170104':143, 'GW170729':91.2,'GW170809':136.4,'GW170814':161, 'GW170818':128.2, 'GW170823':102}


for event in event_list:
	post_loc = '../data/LIGO_events/%s_IMRPPv2_O2catalog.dat'%event
	data = np.genfromtxt(post_loc, dtype=None, names=True)

	max_pos, max_i = maxP(data['logl'], data['logprior'])	

	param_list = ['m1', 'm2', 'a1', 'a2', 'tilt1', 'tilt2', 'phi12', 'mf_evol', 'af_evol']
	median_map = {'m1':[], 'm2':[], 'a1':[], 'a2':[], 'tilt1':[], 'tilt2':[], 'phi12':[], 'mf_evol':[], 'af_evol':[]}	
	maP_map = {'m1':[], 'm2':[], 'a1':[], 'a2':[], 'tilt1':[], 'tilt2':[], 'phi12':[], 'mf_evol':[], 'af_evol':[]}	
	mode_map = {'m1':[], 'm2':[], 'a1':[], 'a2':[], 'tilt1':[], 'tilt2':[], 'phi12':[], 'mf_evol':[], 'af_evol':[]}	

	for param in param_list:
		median_map[param].append(np.median(data[param]))
		maP_map[param].append(data[param][max_i])
		mode_map[param].append(stats.mode(data[param]))

	f_isco_Kerr_median = nr.calc_isco_freq(median_map['af_evol'][0])/(median_map['mf_evol'][0]*lal.MTSUN_SI)
	f_isco_Kerr_maP = nr.calc_isco_freq(maP_map['af_evol'][0])/(maP_map['mf_evol'][0]*lal.MTSUN_SI)
	f_isco_Kerr_mode = nr.calc_isco_freq(mode_map['af_evol'][0][0])/(mode_map['mf_evol'][0][0]*lal.MTSUN_SI)


	print " %s: %.1f %.1f %.1f %.1f"%(event, fcut_map[event],f_isco_Kerr_median, f_isco_Kerr_maP, f_isco_Kerr_mode)	
