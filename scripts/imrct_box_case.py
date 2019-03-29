"""
A box test case for the IMR consistency test on C02 GW150914 data 
"""

import os

lal_prefix = "/home/abhirup/opt/lalsuite_imrct_kde_implementation_20180827_a7914c9b_updates" 
out_dir = "/home/abhirup/Documents/Work/testGR_IR/runs/IMRCT_code_changes/lalsuite_imrct_kde_implementation_20180827_a7914c9b_updates"

command="python %s/libexec/lalinference/imrtgr_imr_consistency_test.py --insp-post=/home/abhirup/Documents/Work/O2/2017/January/04/1167559936p5991/G268556/lalinference/20180629_prodruns_C02cleaned_IMRPPv2_BWPSDs_dL5000/inspiral/lalinferencenest/IMRPhenomPv2pseudoFourPN/1167559936.6-0/H1L1/posterior_samples.dat --ring-post=/home/abhirup/Documents/Work/O2/2017/January/04/1167559936p5991/G268556/lalinference/20180629_prodruns_C02cleaned_IMRPPv2_BWPSDs_dL5000/post-inspiral/lalinferencenest/IMRPhenomPv2pseudoFourPN/1167559936.6-0/H1L1/posterior_samples.dat --imr-post=/home/abhirup/Documents/Work/O2/2017/January/04/1167559936p5991/G268556/lalinference/20180629_prodruns_C02cleaned_IMRPPv2_BWPSDs_dL5000/IMR/lalinferencenest/IMRPhenomPv2pseudoFourPN/1167559936.6-0/H1L1/posterior_samples.dat --out-dir=%s --insp-fhigh=143 --ring-flow=143 --fit-formula=bbh_average_fits_precessing --waveform=IMRPhenomPv2pseudoFourPN --N_bins=401"%(lal_prefix, out_dir)

print(command)
