import lal, os, commands, sys, numpy as np
import imrtestgr as tgr
import nr_fits as nr

m1_targ_list = [31.]
m2_targ_list = [19.]
portion_list = ['inspiral', 'post-inspiral', 'IMR']

approx_targ_list = ['SEOBNRv2', 'SEOBNRv2_ROM_DoubleSpin', 'SEOBNRv4', 'SEOBNRv4_ROM']
approx_templ_list = ['SEOBNRv2_ROM_DoubleSpin', 'SEOBNRv4_ROM']
domain_map = {'SEOBNRv2':'time', 'SEOBNRv2_ROM_DoubleSpin':'freq', 'SEOBNRv4':'time', 'SEOBNRv4_ROM':'freq'}

idx = 100

for m1_targ in m1_targ_list:
 for m2_targ in m2_targ_list:
  for portion in portion_list:
   for approx_targ in approx_targ_list:
    for approx_templ in approx_templ_list:

	os.system('nohup python calcmatch_for_waveform_systematics.py --m1-targ %.2f --m2-targ %.2f --portion %s --approx-targ %s --domain-targ %s --approx-templ %s --domain-templ %s > stdout_%d.out 2> stderr_%d.err &'%(m1_targ, m2_targ, portion, approx_targ, domain_map[approx_targ], approx_templ, domain_map[approx_templ], idx, idx))
	idx += 1
