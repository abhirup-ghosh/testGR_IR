import numpy as np
import matplotlib.pyplot as plt
import glob
import os

def sample_range(postloc):
    data = np.genfromtxt(postloc, names=True, dtype=None)
    params_list = ['m1','m2','distance']
    range_map = {'m1':[],'m2':[],'distance':[]}
    
    for param in params_list:
        range_map[param] = [min(data[param]),max(data[param])]
    
    return range_map

postlocroot_list = glob.glob('/home/abhirup/Documents/Work/testGR_IR/runs/systematics_error_characterisation/*/*')

f = open('./data_wo_prior_railings.txt', 'w')

for postlocroot in postlocroot_list:
    
    try:  
        postloc_imr = glob.glob(postlocroot + '/IMR/lalinferencenest/IMRPhenomPv2pseudoFourPN/*/*/posterior_samples.dat')[0]
        postloc_insp = glob.glob(postlocroot + '/inspiral/lalinferencenest/IMRPhenomPv2pseudoFourPN/*/*/posterior_samples.dat')[0]
        postloc_ring = glob.glob(postlocroot + '/post-inspiral/lalinferencenest/IMRPhenomPv2pseudoFourPN/*/*/posterior_samples.dat')[0]    
        
        if os.path.isfile(postloc_imr) == True and os.path.isfile(postloc_insp) == True and os.path.isfile(postloc_ring) == True:

            print '... data found:', postlocroot  

            range_imr_map = sample_range(postloc_imr)
            range_insp_map = sample_range(postloc_insp)  
            range_ring_map = sample_range(postloc_ring)

            if range_imr_map['m1'][1] < 198. or range_imr_map['m2'][0] > 12. or range_imr_map['distance'][1] < 4998 or range_insp_map['m1'][1] < 198. or range_insp_map['m2'][0] > 12. or range_insp_map['distance'][1] < 4998 or range_ring_map['m1'][1] < 198. or range_ring_map['m2'][0] > 12. or range_ring_map['distance'][1] < 4998:
                print '... data range safe'
                f.write('%s \n'%postlocroot)
                
            else:
                print '... data range not safe'

                
    except IndexError:
        print '... all 3 data not found'
        
f.close()        

