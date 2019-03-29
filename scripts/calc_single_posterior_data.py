import os
import glob

post_root = "/home/abhirup/Documents/Work/O2/2017/August/23/1187529256p5179/G298936/lalinference/Reruns/20180629_mcprior_150_dL_10000_IMRPPv2_C02_BWPSD_seglen8_HL_scseg"

imr_post = glob.glob("%s/IMR/lalinferencenest/*/*/*/posterior_samples.dat"%post_root)[0]
insp_post = glob.glob("%s/inspiral/lalinferencenest/*/*/*/posterior_samples.dat"%post_root)[0]
ring_post = glob.glob("%s/post-inspiral/lalinferencenest/*/*/*/posterior_samples.dat"%post_root)[0]
out_dir = post_root + "/imrtgr_results"

fit_formula = "bbh_average_fits_precessing"
insp_fhigh, ring_flow = 102,102
waveform = "IMRPhenomPv2pseudoFourPN" 
N_bins = 401

command = "python /home/abhirup/Documents/Work/testGR_IR/scripts/imrtgr_imr_consistency_test_final.py --insp-post=%s --ring-post=%s --imr-post=%s --out-dir=%s --insp-fhigh=%s --ring-flow=%s --fit-formula=%s --waveform=%s --N_bins=%d"%(insp_post, ring_post, imr_post, out_dir, insp_fhigh, ring_flow, fit_formula, waveform, N_bins)

print command
#os.system(command)
