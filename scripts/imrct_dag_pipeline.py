import numpy as n
import os
import glob
import time

"""
Pre-requisites for running this script:
1. create the run_dir specified here
2. create a universal config.ini (identical for all three runs: IMR, inspiral and post-inspiral, except for the options webdir, baseurl, flow and fhigh) and place it under run_dir.

Workings of this code:
1. the options for webdir, baseurl, flow and fhigh are changed appropriate each of the three LALInference runs: IMR, inspiral and post-inspiral by:
a. creating a subdirectory: 'IMR', 'inspiral', 'post-inspiral' under run_dir
b. copying the universal config.ini to subdir
c. replacing webdir and baseurl lines by paths to subdir
b. replacing flow line with flow appropriate for 'IMR', 'inspiral' and 'post-inspiral'
e. adding an fhigh line for 'inspiral'

2. Once LALInference Runs are done, code performs the IMRCT in the following steps:
a. looks for all 3 posterior_samples.dat files every 6 hours; otherwise sleeps
b. if found, executes IMRCT

"""

# LALInference arguments
run_dir = "/home/abhirup/public_html/O2/2017/August/23/1187529256p5179/G298936/lalinference/Reruns/test"
input_config = run_dir + '/config.ini'

# lines to the replaced in the LALInference config.ini
webdir_line = "webdir = /home/abhirup/public_html/O2/2017/August/23/1187529256p5179/G298936/lalinference/Reruns/test"
baseurl_line = "baseurl = https://dumpty.alice.icts.res.in/~abhirup/O2/2017/August/23/1187529256p5179/G298936/lalinference/Reruns/test"

cbcstage_list = ['IMR', 'inspiral', 'post-inspiral']
f_low = 20.
f_cut = 100.

# IMRCT arguments
fit_formula="bbh_average_fits_precessing"
N_bins=401
insp_fhigh=f_cut; ring_flow=f_cut
waveform="SEOBNRv4_ROMpseudoFourPN"
imrct_out_dir = run_dir + '/imrct_results'
os.system('mkdir -p %s'%imrct_out_dir)

# Manipulating LALInference config.ini's for each run, and launching them

for cbcstage in cbcstage_list:

  # creating subdir
  os.system('mkdir -p %s/%s'%(run_dir, cbcstage))

  # searching and replacing webdir and baseurl lines
  os.system("sed '/^webdir/c\%s\/%s' %s/config.ini > %s/%s/config.ini"%(webdir_line, cbcstage, run_dir, run_dir, cbcstage))
  os.system("sed -i '/^baseurl/c\%s\/%s' %s/%s/config.ini"%(baseurl_line, cbcstage, run_dir, cbcstage))

  # searching and replacing the f_cut lines appropriate to the LALInference run
  if cbcstage == "IMR":
	f_low_line = "flow = {\\'V1\\': %.2f, \\'H1\\': %.2f, \\'L1\\': %.2f}"%(f_low, f_low, f_low)
	os.system("sed -i '/^flow/c\%s' %s/%s/config.ini"%(f_low_line, run_dir, cbcstage))
  if cbcstage == "inspiral":
        f_low_line = "flow = {\\'V1\\': %.2f, \\'H1\\': %.2f, \\'L1\\': %.2f}"%(f_low, f_low, f_low)
        f_high_line = "fhigh = {\\'V1\\': %.2f, \\'H1\\': %.2f, \\'L1\\': %.2f}"%(f_cut, f_cut, f_cut)
	os.system("sed -i '/^flow/c\%s' %s/%s/config.ini"%(f_low_line, run_dir, cbcstage))
	os.system("sed -i '/^flow/a\%s' %s/%s/config.ini"%(f_high_line, run_dir, cbcstage))
  if cbcstage == "post-inspiral":
        f_low_line = "flow = {\\'V1\\': %.2f, \\'H1\\': %.2f, \\'L1\\': %.2f}"%(f_cut, f_cut, f_cut)
        os.system("sed -i '/^flow/c\%s' %s/%s/config.ini"%(f_low_line, run_dir, cbcstage))

  # creating and submitting multidag.dag
  os.system('lalinference_pipe -g %s/gpstime.txt -r %s/%s -p %s/%s %s/%s/config.ini'%(run_dir, run_dir, cbcstage, run_dir, cbcstage, run_dir, cbcstage))
  os.system('condor_submit_dag %s/%s/multidag.dag'%(run_dir, cbcstage))


# Waiting and running the IMRCT

for idx in range(10000):
  
  try:
	imr_post = glob.glob(run_dir + '/IMR/*/*/*/*/posterior_samples.dat')[0]
  	insp_post = glob.glob(run_dir + '/inspiral/*/*/*/*/posterior_samples.dat')[0]
  	ring_post = glob.glob(run_dir + '/post-inspiral/*/*/*/*/posterior_samples.dat')[0]

	if  os.path.isfile(imr_post) == True and os.path.isfile(insp_post) == True and os.path.isfile(ring_post) == True:
		print '... woohoo!'
		command="python /home/abhirup.ghosh/Documents/Work/testGR_IR/scripts/imrtgr_imr_consistency_test_final.py --insp-post=%s --ring-post=%s --imr-post=%s --out-dir=%s --insp-fhigh=%f --ring-flow=%s --fit-formula=%s --waveform=%s --N_bins=%d"%(insp_post,ring_post,imr_post,outdir,insp_fhigh,ring_flow,fit_formula,waveform,N_bins)

  except IndexError:
		print '... posterior files not created; last checked at: %s'%(time.strftime('%d/%m/%Y %H:%M:%S'))
		time.sleep(216000)	

