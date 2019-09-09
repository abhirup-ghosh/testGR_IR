import numpy as np
import matplotlib.pyplot as plt
import os, sys
np.set_printoptions(threshold=sys.maxsize)

# Module for confidence calculations
class confidence(object):
    def __init__(self, counts):
        # Sort in descending order in frequency
        self.counts_sorted = np.sort(counts.flatten())[::-1]
        # Get a normalized cumulative distribution from the mode
        self.norm_cumsum_counts_sorted = np.cumsum(self.counts_sorted) / np.sum(counts)
        # Set interpolations between heights, bins and levels
        self._set_interp()
    def _set_interp(self):
        self._length = len(self.counts_sorted)
        # height from index
        self._height_from_idx = interpolate.interp1d(np.arange(self._length), self.counts_sorted, bounds_error=False, fill_value=0.)
        # index from height
        self._idx_from_height = interpolate.interp1d(self.counts_sorted[::-1], np.arange(self._length)[::-1], bounds_error=False, fill_value=self._length)
        # level from index
        self._level_from_idx = interpolate.interp1d(np.arange(self._length), self.norm_cumsum_counts_sorted, bounds_error=False, fill_value=1.)
        # index from level
        self._idx_from_level = interpolate.interp1d(self.norm_cumsum_counts_sorted, np.arange(self._length), bounds_error=False, fill_value=self._length)
    def level_from_height(self, height):
        return self._level_from_idx(self._idx_from_height(height))
    def height_from_level(self, level):
        return self._height_from_idx(self._idx_from_level(level))

# gaussian filter with 1-sigma softening
def gf(P):
    return filter.gaussian_filter(P, sigma=1.0)

# compute 1-sigma confidence intervals in 1D
def calc_cred_intervals_in_1d(P, x):

    # find the value of P corresponding to 50% and 9% confidence heights
    conf = confidence(P)
    P_s1 = conf.height_from_level(0.5)
    P_s2 = conf.height_from_level(0.9)

    # calculation of condifence edges (values of x corresponding to the height s1 on the two sides)
    x_s1_l = min(x[np.where(P >= P_s1)[0]])
    x_s1_r = max(x[np.where(P >= P_s1)[0]])

    # calculation of condifence edges (values of x corresponding to the height s2 on the two sides)
    x_s2_l = min(x[np.where(P >= P_s2)[0]])
    x_s2_r = max(x[np.where(P >= P_s2)[0]])

    return P_s1, P_s2, x_s1_l, x_s1_r, x_s2_l, x_s2_r

def matched_filter_snr_module(filename):
    data = np.genfromtxt(filename, dtype=None, names=True, usecols=(0,1,2,3,4,5,6))
    var_names = [d[0] for d in data]
    stat_names = data.dtype.names
    matched_filter_snr = data[var_names.index('matched_filter_snr')][stat_names.index('median')+1]
    return matched_filter_snr

def sample_range(postloc):
    data = np.genfromtxt(postloc, names=True, dtype=None)
    params_list = ['m1','m2','distance']
    range_map = {'m1':[],'m2':[],'distance':[]}

    for param in params_list:
        range_map[param] = [min(data[param]),max(data[param])]

    return range_map

post_loc_root_list = ['imrppv2_injection_o1o2_noise','imrppv2_injection_psdmodel_zeronoise','sxs_injection_o1o2_noise','sxs_injection_psdmodel_zeronoise']
inj_data = np.genfromtxt('/home/abhirup/Documents/Work/testGR_IR/runs/systematics_error_characterisation/SXS_campaign.dat', names=True, dtype=None)   
valid_inj_list_map = {'imrppv2_injection_o1o2_noise':[],'imrppv2_injection_psdmodel_zeronoise':[],'sxs_injection_o1o2_noise':[],'sxs_injection_psdmodel_zeronoise':[]}

for post_loc_root in post_loc_root_list:
    
    post_loc_root = '/home/abhirup/Documents/Work/testGR_IR/runs/systematics_error_characterisation/%s'%post_loc_root

    imrtgr_data_list = []
    valid_inj_list = []

    for idx in range(len(inj_data)):
        tag = inj_data[idx]['tag']
        if os.path.isfile(post_loc_root + '/%s/imrtgr_results/data/P_dMfbyMf_dchifbychif.dat.gz'%(tag)) == True:
            imrtgr_data_list.append(post_loc_root + '/%s/imrtgr_results'%(tag)) 
    
    for (idx, post_loc) in enumerate(imrtgr_data_list):

        M = inj_data[idx]['m1'] + inj_data[idx]['m2']
        snr_i = matched_filter_snr_module(post_loc + '/lalinf_insp/summary_statistics.dat')
        snr_r = matched_filter_snr_module(post_loc + '/lalinf_ring/summary_statistics.dat')
        gr_height = np.loadtxt(post_loc + '/data/GR_confidence.txt')

        if snr_i > 5 and snr_r > 5 and M < 150:
            range_map_i = sample_range(post_loc + '/lalinf_insp/posterior_samples.dat')
            range_map_r = sample_range(post_loc + '/lalinf_ring/posterior_samples.dat')
            if range_map_i['m1'][1] < 198. or range_map_i['m2'][0] > 12. or range_map_i['distance'][1] < 4998 or range_map_r['m1'][1] < 198. or range_map_r['m2'][0] > 12. or range_map_r['distance'][1] < 4998:  
                   valid_inj_list.append(idx)
                    
    valid_inj_list_map[post_loc_root] = valid_inj_list
    print post_loc_root, valid_inj_list

