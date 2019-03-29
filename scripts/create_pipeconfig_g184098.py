#!/usr/bin/env python

"""
This script creates an injection and a .ini file to run LALInference

(C) Archisman Ghosh, Abhirup Ghosh, 2015-09-20
"""

import sys, lal, os, commands, numpy as np
sys.path.insert(1, os.path.join(os.path.expanduser('~'), 'src/lalsuite/lalinference/python'))
import imrtestgr as tgr
import nr_fits as nr
import matplotlib.pyplot as plt
import chirptimes as ct

os.system('. ${HOME}/.profile')
home = os.getenv('HOME')
user = commands.getoutput('whoami')
lal_prefix = os.getenv('LAL_PREFIX')
psd_path = '%s/share/lalsimulation'%lal_prefix
glue_location = os.getenv('GLUE_LOCATION')
pylal_location = os.getenv('PYLAL_LOCATION')

# ### INJECTION ###
# Waveform info
approximant = 'SEOBNRv2_ROM_DoubleSpinpseudoFourPN' # injection waveform 
#approximant = 'IMRPhenomPv2pseudoFourPN' # injection waveform 
amp_order = -1 # amplitude PN order of the injection waveform
f_low = 8. # low-frequency cutoff (Hz)
waveform_info = '--waveform %s --amp-order %d --f-lower %f --taper-injection startend'%(approximant, amp_order, f_low)

# PSD info
ligo_psd = 'LALAdLIGO'
virgo_psd = 'LALAdVirgo'
f_start = 10.
psd_info = '--ifos H1,L1,V1 --ligo-fake-psd %s --virgo-fake-psd %s --ligo-start-freq %f --virgo-start-freq %f'%(ligo_psd, virgo_psd, f_start, f_start)

# Time info
gps_time = 1126259462 # O1 start time 
time_step = 2630./np.pi # time interval between nearby injections 
time_info = '--t-distr fixed --gps-start-time %f --gps-end-time %f --time-step %f'%(gps_time, gps_time, time_step)

# Parameters
m1 = 36.
m2 = 29.
chi1 = 0.3200
chi2 = 0.5798
iota = 0.
phi_c = 0.
psi = 0.
alpha = 0.
delta = 0.
snr_inj = 25.

param_info = '--m-distr fixMasses --fixed-mass1 %f --fixed-mass2 %f --min-spin1 %f --max-spin1 %f --min-spin2 %f --max-spin2 %f --enable-spin --aligned --snr-distr uniform --min-snr %f --max-snr %f --i-distr fixed --fixed-inc %f --coa-phase-distr fixed --fixed-coa-phase %f --polarization %f --l-distr fixed --latitude %f --longitude %f'%(m1, m2, chi1, chi1, chi2, chi2, snr_inj, snr_inj, iota, phi_c, psi, alpha, delta)

# Output info
inj_file = '../injections/gw150914_%s_m1_36_m2_29_snr_25_chi1_0p32_chi2_0p5798_optimallyoriented'%(approximant) # output file name for xml and txt file
output_info = '--output %s.xml'%(inj_file)

# Generate the injection xml file by calling lalapps_inspinj 
run_command = 'lalapps_inspinj %s %s %s %s %s'%(waveform_info, psd_info, time_info, param_info, output_info)
print(run_command)
os.system(run_command)
   
# Print some relevant columns of the xml file to an ASCII table with header before setting spins
os.system('ligolw_print -d " " -c mass1 -c mass2 -c spin1x -c spin2x -c spin1y -c spin2y -c spin1z -c spin2z -c longitude -c latitude -c inclination -c polarization -c distance -c geocent_end_time -c geocent_end_time_ns -c coa_phase %s.xml > %s_all_param.tmp' %(inj_file, inj_file))
os.system('ligolw_print -d " " -c spin1z -c spin2z %s.xml > %s_spin_param.tmp' %(inj_file, inj_file))
os.system('echo \# mass1 mass2 spin1x spin2x spin1y spin2y spin1z spin2z longitude latitude inclination polarization distance geocent_end_time geocent_end_time_ns coa_phase > header.txt')
os.system('cat header.txt %s_all_param.tmp > %s.txt' %(inj_file, inj_file))
os.system('rm %s_spin_param.tmp %s_all_param.tmp header.txt' %(inj_file, inj_file))

# ### RECOVERY ###
ligo_psd = os.path.join(psd_path, 'LIGO-T0900288-v3-ZERO_DET_high_P.txt')
virgo_psd = os.path.join(psd_path, 'LIGO-T0900288-v3-ZERO_DET_high_P.txt')

def write_pipeline(out_folder, seglen, dataseed, randomseed, flow=20., fhigh=1023., nlive=1024, create_dag=True):
    os.system('mkdir -p %s'%(out_folder))

    inj_file = '../injections/gw150914_%s_m1_36_m2_29_snr_25_chi1_0p32_chi2_0p5798_optimallyoriented'%(approximant)
    os.system('cp %s.xml %s/injection.xml'%(inj_file, out_folder))
    os.system('cp %s.txt %s/injection.txt'%(inj_file, out_folder))

    ofile = open(os.path.join(out_folder, 'pipeconfig.ini'), 'w')
    
    ofile.write("[analysis]\n")
    ofile.write("ifos=['H1','L1','V1']\n")
    ofile.write("engine=lalinferencenest\n")
    ofile.write("nparallel=4\n")
    ofile.write("coherence-test=False\n")
    ofile.write("upload-to-gracedb=False\n")
    ofile.write("dataseed=%d\n"%dataseed) # ###FIXME###
    ofile.write("\n")
    
    projectdir = out_folder.replace('%s/Documents/Work/'%home, '')

    ofile.write("[paths]\n")
    ofile.write("webdir=%s/public_html/%s\n"%(home, projectdir))
    ofile.write("basedir=%s/public_html/%s\n"%(home, projectdir))
    #ofile.write("baseurl=https://dogmatix.icts.res.in/~%s/%s\n"%(user, projectdir))
    ofile.write("baseurl=https://humpty.alice.icts.res.in/~%s/%s\n"%(user, projectdir))
    ofile.write("\n")
    
    ofile.write("[input]\n")
    ofile.write("max-psd-length=1024\n")
    ofile.write("padding=16\n")
    ofile.write("timeslides=false\n")
    ofile.write("ignore-science-segments=True\n")
    ofile.write("events=all\n")
    ofile.write("\n")
    
    ofile.write("[datafind]\n")
    ofile.write("types={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo'}\n")
    ofile.write("\n")
    
    ofile.write("[data]\n")
    ofile.write("channels={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo'}\n")
    ofile.write("\n")
    
    ofile.write("[condor]\n")
    ofile.write("mpirun=/bin/true\n")
    ofile.write("datafind=/bin/true\n")
    ofile.write("gracedb=/bin/true\n")
    ofile.write("ligolw_print=%s/bin/ligolw_print\n"%(glue_location))
    ofile.write("lalinferencenest=%s/bin/lalinference_nest\n"%(lal_prefix))
    ofile.write("lalinferencemcmc=%s/bin/lalinference_mcmc\n"%(lal_prefix))
    ofile.write("mergescript=%s/bin/lalapps_nest2pos\n"%(lal_prefix))
    ofile.write("resultspage=%s/bin/cbcBayesPostProc.py\n"%(pylal_location))
    ofile.write("coherencetest=%s/bin/lalapps_coherence_test\n"%(lal_prefix))
    ofile.write("lalinferencedatadump=False\n")
    ofile.write("\n")   
 
    ofile.write("[engine]\n")
    ofile.write("approx=%s\n"%(approximant))
    ofile.write("seglen=%d\n"%(seglen))
    ofile.write("nlive=%d\n"%(nlive))
    ofile.write("srate=2048\n")
    ofile.write("distance-min=1\n") # ###FIXME###
    ofile.write("distance-max=20000\n") # ###FIXME###
    ofile.write("comp-min=10.0\n") 
    ofile.write("comp-max=300.0\n") 
    ofile.write("a_spin1-max=0.99\n") # ###FIXME###
    ofile.write("a_spin2-max=0.99\n") # ###FIXME###
    ofile.write("amporder=-1\n")
    ofile.write("resume=\n")
    ofile.write("progress=\n")
    ofile.write("H1-asd=%s\n"%(ligo_psd))
    ofile.write("L1-asd=%s\n"%(ligo_psd))
    ofile.write("V1-asd=%s\n"%(virgo_psd))
    ofile.write("randomseed=%d\n"%randomseed) # ### FIXME ###
    ofile.write("\n")
    
    ofile.write("[lalinference]\n")
    ofile.write("fake-cache={'H1':'interp:%s','L1':'interp:%s','V1':'interp:%s'}\n"%(ligo_psd, ligo_psd, virgo_psd))
    ofile.write("flow={'H1':%.2f,'L1':%.2f,'V1':%.2f}\n"%(flow, flow, flow))
    ofile.write("fhigh={'H1':%.2f,'L1':%.2f,'V1':%.2f}\n"%(fhigh, fhigh, fhigh))
    ofile.write("\n")
    
    ofile.write("[merge]\n")
    ofile.write("\n")
    
    ofile.write("[resultspage]\n")
    ofile.write("skyres=0.5\n")
    ofile.write("\n")
    
    ofile.close()

    if create_dag:
      os.system('%s/bin/lalinference_pipe -I %s/injection.xml -r %s -p %s %s/pipeconfig.ini'%(lal_prefix, out_folder, out_folder, out_folder, out_folder))
      os.system('/usr/bin/condor_submit_dag %s/multidag.dag'%out_folder)


# ### IMR ###
noise = 'gaussian_noise'
noise_type = 'LALAdLIGO'#'LIGO-P1200087-v18-aLIGO_EARLY_HIGH'
date = '2016-11-11' # ###FIXME###
run_tag = 'SEOB_inj_SEOB_rec_diff_freq_flow20Hz'# ###FIXME###

imr_folder = '%s/Documents/Work/testGR_IR/runs/test/%s/%s/IMR_consistency_robustness_tests/%s/%s/IMR'%(home, noise, noise_type, date, run_tag)
insp_folder = '%s/Documents/Work/testGR_IR/runs/test/%s/%s/IMR_consistency_robustness_tests/%s/%s/inspiral'%(home, noise, noise_type, date, run_tag)
post_insp_folder = '%s/Documents/Work/testGR_IR/runs/test/%s/%s/IMR_consistency_robustness_tests/%s/%s/post-inspiral'%(home, noise, noise_type, date, run_tag)

fit_formula = 'nonprecspin_Healy2014'

# calculate the mass and spin of the final BH 
Mf, af = tgr.calc_final_mass_spin(m1, m2, chi1, chi2, fit_formula)

# calculate the Kerr ISCO freq 
f_isco_Kerr = nr.calc_isco_freq(af)/(Mf*lal.MTSUN_SI)

# calculate the dominant QNM freq 
f_qnm = nr.calc_fqnm_dominant_mode(af)/(Mf*lal.MTSUN_SI)
 
insp_fhigh = f_isco_Kerr

# calculating the appropriate seglen
ctime_r = ct.calc_chirptime(m1, m2, f_min=f_start)
if ctime_r < 6:
        seglen = 8
elif (ctime_r>6) and (ctime_r<14):
        seglen = 16
elif (ctime_r>14) and (ctime_r<30):
        seglen = 32
elif (ctime_r>30) and (ctime_r<60):
        seglen = 64
elif (ctime_r>60) and (ctime_r<120):
        seglen = 128

print seglen

seglen = 8
dataseed = np.random.randint(1,1000)
randomseed = np.random.randint(1,1000)
    
write_pipeline(imr_folder, seglen, dataseed, randomseed)
write_pipeline(insp_folder, seglen, dataseed, randomseed, fhigh=insp_fhigh)
write_pipeline(post_insp_folder, seglen, dataseed, randomseed, flow=insp_fhigh)
