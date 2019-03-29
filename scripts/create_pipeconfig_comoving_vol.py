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
psd_path = os.path.join(home, 'src/lalsuite/lalsimulation/src')
lal_prefix = os.getenv('LAL_PREFIX')
glue_location = os.getenv('GLUE_LOCATION')
pylal_location = os.getenv('PYLAL_LOCATION')

# ### INJECTION ###
no_inj = 10000

# Waveform info
approximant = 'SEOBNRv2_ROM_DoubleSpinthreePointFivePN' # injection waveform 
amp_order = -1 # amplitude PN order of the injection waveform
f_low = 18. # low-frequency cutoff (Hz)
waveform_info = '--waveform %s --amp-order %d --f-lower %f --taper-injection startend'%(approximant, amp_order, f_low)

# PSD info
f_start = 20.
psd_info = '--ligo-start-freq %f --virgo-start-freq %f'%(f_start, f_start)

# Time info
gps_start_time = 1126285216 # O1 start time 
time_step = 2630./np.pi # time interval between nearby injections 
gps_end_time = gps_start_time + time_step*no_inj
time_info = '--t-distr fixed --gps-start-time %f --gps-end-time %f --time-step %f'%(gps_start_time, gps_end_time, time_step)

# Parameters
m1_min = 10.
m2_min = 10.
m1_max = 80.
m2_max = 80.
chi1_min = 0.
chi1_max = 1.
chi2_min = 0.
chi2_max = 1.
z_min = 0.
z_max = 0.25
phi_c = 0. # coalescent phase (degrees)

psi = 0. # polarization (degrees)
param_info = '--m-distr componentMass --min-mass1 %f --max-mass1 %f --min-mass2 %f --max-mass2 %f --enable-spin --aligned --min-spin1 %f --max-spin1 %f --min-spin2 %f --max-spin2 %f --i-distr uniform --coa-phase-distr fixed --fixed-coa-phase %f --polarization %f --l-distr random --d-distr volume --min-z %f --max-z %f'%(m1_min, m1_max, m2_min, m2_max, chi1_min, chi1_max, chi2_min, chi2_max, phi_c, psi, z_min, z_max)

# Output info
inj_file = '../injections/popsynth_injections_%s_distance_reach'%(approximant) # output file name for xml and txt file
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

m1_inj, m2_inj, chi1_inj, chi2_inj = np.loadtxt('%s.txt'%(inj_file), usecols=(0,1,6,7), unpack=True)
exit()
# ### RECOVERY ###

ligo_psd = os.path.join(psd_path, 'LIGO-T0900288-v3-ZERO_DET_high_P.txt')
virgo_psd = os.path.join(psd_path, 'LIGO-T0900288-v3-ZERO_DET_high_P.txt')


def write_pipeline(out_folder, idx, seglen, dataseed, randomseed, flow=10., fhigh=1023., nlive=1024, create_dag=True):
    os.system('mkdir -p %s'%(out_folder))

    os.system('cp %s.xml %s/injection.xml'%(inj_file, out_folder))
    os.system('cp %s.txt %s/injection.txt'%(inj_file, out_folder))

    ofile = open(os.path.join(out_folder, 'pipeconfig.ini'), 'w')
    
    ofile.write("[analysis]\n")
    ofile.write("ifos=['H1','L1','V1']\n")
    ofile.write("engine=lalinferencenest\n")
    ofile.write("nparallel=4\n")
    ofile.write("coherence-test=False\n")
    ofile.write("upload-to-gracedb=False\n")
    #ofile.write("dataseed=%d\n"%dataseed) # ###FIXME###
    ofile.write("%s\n"%dataseed) # ###FIXME###
    #ofile.write("accounting_group=ligo.sim.o1.cbc.testgr.tiger\n")
    ofile.write("\n")
    
    projectdir = out_folder.replace('%s/Documents/Work/'%home, '')

    ofile.write("[paths]\n")
    ofile.write("webdir=%s/public_html/%s\n"%(home, projectdir))
    ofile.write("basedir=%s/public_html/%s\n"%(home, projectdir))
    #ofile.write("baseurl=https://dogmatix.icts.res.in/~%s/%s\n"%(user, projectdir))
    ofile.write("baseurl=https://dumpty.alice.icts.res.in/~%s/%s\n"%(user, projectdir))
    ofile.write("\n")
    
    ofile.write("[input]\n")
    ofile.write("max-psd-length=1024\n")
    ofile.write("padding=16\n")
    ofile.write("timeslides=false\n")
    ofile.write("ignore-science-segments=True\n")
    ofile.write("events=%d\n"%idx)
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
    #ofile.write("dt=0.1\n")
    ofile.write("distance-min=1\n") # ###FIXME###
    ofile.write("distance-max=20000\n") # ###FIXME###
    ofile.write("comp-min=1.0\n") 
    #ofile.write("comp-max=%.2f\n"%comp_mass_max) 
    ofile.write("comp-max=300.0\n") 
    ofile.write("a_spin1-max=0.98\n") # ###FIXME###
    ofile.write("a_spin2-max=0.98\n") # ###FIXME###
    #ofile.write("disable-spin=\n") # ### FIXME ###
    ofile.write("amporder=-1\n")
    ofile.write("margphi=\n") # ###FIXME###
    #ofile.write("marginal-d=\n")
    ofile.write("resume=\n")
    ofile.write("progress=\n")
    ofile.write("H1-psd=%s\n"%(ligo_psd))
    ofile.write("L1-psd=%s\n"%(ligo_psd))
    ofile.write("V1-psd=%s\n"%(virgo_psd))
    #ofile.write("randomseed=%d\n"%randomseed) # ### FIXME ###
    ofile.write("%s\n"%randomseed) # ### FIXME ###
    #ofile.write("0noise=\n") 
    #ofile.write("tolerance=0.05\n")
    #ofile.write("pinparams=[time,phase,costheta_jn,distance,rightascension,declination,polarisation]\n") # ###FIXME####
    ofile.write("\n")
    
    ofile.write("[lalinference]\n")
    ofile.write("fake-cache={'H1':'interp:%s','L1':'interp:%s','V1':'interp:%s'}\n"%(ligo_psd, ligo_psd, virgo_psd))
    ofile.write("flow={'H1':%.2f,'L1':%.2f,'V1':%.2f}\n"%(flow, flow, flow))
    ofile.write("fhigh={'H1':%.2f,'L1':%.2f,'V1':%.2f}\n"%(fhigh, fhigh, fhigh))
    ofile.write("\n")
    
    ofile.write("[merge]\n")
    #ofile.write("npos=5000\n")
    ofile.write("\n")
    
    ofile.write("[resultspage]\n")
    ofile.write("skyres=0.5\n")
    #ofile.write("no2D=\n")
    ofile.write("\n")
    
    ofile.close()

    if create_dag:
      os.system('%s/bin/lalinference_pipe -I %s/injection.xml -r %s -p %s %s/pipeconfig.ini'%(lal_prefix, out_folder, out_folder, out_folder, out_folder))
      os.system('/usr/bin/condor_submit_dag %s/lalinference_*.dag'%out_folder)


#inj_list = np.loadtxt('inj_thresh5_seglen16.txt', dtype='int')
inj_list = np.loadtxt('inj_thresh8.txt')
injection_no_list = inj_list[200:235]

# Main
for idx in injection_no_list:
    # ### IMR ###
    noise = 'gaussian_noise'
    noise_type = 'LALAdLIGO'#'LIGO-P1200087-v18-aLIGO_EARLY_HIGH'
    date = '2016-01-16_test' # ###FIXME###
    run_tag = 'uniform_compmass_spins_comoving_volume'#'realization_2' # ###FIXME###

    imr_folder = '%s/Documents/Work/testGR_IR/runs/simulations/%s/%s/%s/%s/%s/injection_%d/IMR'%(home, noise, noise_type, approximant, date, run_tag, idx+1)
    #insp_folder = '%s/Documents/Work/testGR_IR/runs/simulations/%s/%s/%s/%s/%s/injection_%d/inspiral'%(home, noise, noise_type, approximant, date, run_tag, idx+1)
    #post_insp_folder = '%s/Documents/Work/testGR_IR/runs/simulations/%s/%s/%s/%s/%s/injection_%d/post-inspiral'%(home, noise, noise_type, approximant, date, run_tag, idx+1)
    ring_folder = '%s/Documents/Work/testGR_IR/runs/simulations/%s/%s/%s/%s/%s/injection_%d/ringdown_80percent'%(home, noise, noise_type, approximant, date, run_tag, idx+1)
    
    #fit_formula = 'nospin_Pan2011'
    fit_formula = 'nonprecspin_Healy2014'

    # calculate the mass and spin of the final BH 
    Mf, af = tgr.calc_final_mass_spin(m1_inj[idx], m2_inj[idx], chi1_inj[idx], chi2_inj[idx], fit_formula)

    # calculate the Kerr ISCO freq 
    f_isco_Kerr = nr.calc_isco_freq(af)/(Mf*lal.MTSUN_SI)

    # calculate the dominant QNM freq 
    f_qnm = nr.calc_fqnm_dominant_mode(af)/(Mf*lal.MTSUN_SI)

    # calculating the appropriate seglen
    ctime_r = ct.calc_chirptime(m1_inj[idx], m2_inj[idx], f_min=f_start)
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
 
    print ctime_r, seglen

    insp_fhigh = f_isco_Kerr
    ring_flow = 0.8*f_qnm

    dataseed = commands.getoutput('cat %s/config.ini | grep dataseed'%imr_folder)#np.random.randint(1,1000)
    randomseed = commands.getoutput('cat %s/config.ini | grep randomseed'%imr_folder)#np.random.randint(1,1000)
    
    #write_pipeline(imr_folder, idx, seglen, dataseed, randomseed)
    #write_pipeline(insp_folder, idx, seglen, dataseed, randomseed, fhigh=insp_fhigh)
    #write_pipeline(post_insp_folder, seglen, idx, dataseed, randomseed, flow=insp_fhigh)
    write_pipeline(ring_folder, idx, seglen, dataseed, randomseed, flow=ring_flow)
