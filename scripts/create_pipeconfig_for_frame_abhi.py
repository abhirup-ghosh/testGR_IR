#!/usr/bin/env python

"""
This script creates an injection and a .ini file to run LALInference

(C) Archisman Ghosh, Abhirup Ghosh, 2015-09-20
"""

import sys, lal, os, glob, re, commands, numpy as np
sys.path.insert(1, os.path.join(os.path.expanduser('~'), 'src/lalsuite/lalinference/python'))
import imrtestgr as tgr
import nr_fits as nr
from scipy import random
from write_cache import write_cache
import matplotlib.pyplot as plt
import optparse as op
import chirptimes as ct

os.system('. ${HOME}/.profile')
home = os.getenv('HOME')
user = commands.getoutput('whoami')
lal_prefix = os.getenv('LAL_PREFIX')
psd_path = '%s/share/lalsimulation'%lal_prefix
glue_location = os.getenv('GLUE_LOCATION')
pylal_location = os.getenv('PYLAL_LOCATION')

parser = op.OptionParser()
parser.add_option('--inj-file', dest="inj_file", help="injection file name without any extension")
(options, args) = parser.parse_args()
inj_file = options.inj_file

# Injection Parameters
m1_inj, m2_inj, chi1_inj, chi2_inj = np.loadtxt('%s.txt'%(inj_file), usecols=(0,1,6,7), unpack=True)
fmin = 10.

# ### RECOVERY ###

#approximant = 'SEOBNRv2_ROM_DoubleSpinpseudoFourPN'
approximant = 'IMRPhenomPv2pseudoFourPN'
ligo_psd = os.path.join(psd_path, 'LIGO-T0900288-v3-ZERO_DET_high_P.txt')
virgo_psd = os.path.join(psd_path, 'LIGO-T0900288-v3-ZERO_DET_high_P.txt')



def write_pipeline(out_folder, H1_frame, L1_frame, V1_frame, gpstimefile, seglen, randomseed, flow=10., fhigh=1023., nlive=1024, create_dag=True, submit_dag=True):

    dag_folder = os.path.join(out_folder, 'lalinferencenest', approximant)
    cache_folder = os.path.join(dag_folder, 'caches')
    os.system('mkdir -p %s'%(cache_folder))
  
    (H1_cache,H,H1_channel,start_time,datalen,Hpath) = write_cache(H1_frame, cache_folder).split()
    (L1_cache,L,L1_channel,start_time,datalen,Lpath) = write_cache(L1_frame, cache_folder).split()
    (V1_cache,V,V1_channel,start_time,datalen,Lpath) = write_cache(V1_frame, cache_folder).split()
    
    start_time = int(start_time)
    datalen = int(datalen)
    print "Start time = %d, datalen = %d"%(start_time, datalen)
    ofile = open(os.path.join(out_folder, 'pipeconfig.ini'), 'w')

    ofile.write("[analysis]\n")
    ofile.write("ifos=['H1','L1','V1']\n")
    ofile.write("engine=lalinferencenest\n")
    ofile.write("nparallel=4\n")
    ofile.write("coherence-test=False\n")
    ofile.write("upload-to-gracedb=False\n")
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
    ofile.write("gps-start-time=%d\n"%(start_time))
    ofile.write("gps-end-time=%d\n"%(start_time+datalen+16))
    ofile.write("psd-length=%d\n"%(datalen/2))
    ofile.write("\n")
    
    ofile.write("[datafind]\n")
    ofile.write("types={'H1':'%s','L1':'%s','V1':'%s'}\n"%(H1_channel[3:], L1_channel[3:], V1_channel[3:]))
    ofile.write("\n")
    
    ofile.write("[data]\n")
    ofile.write("channels={'H1':'%s:%s','L1':'%s:%s','V1':'%s:%s'}\n"%(H1_channel[:2], H1_channel[3:], L1_channel[:2], L1_channel[3:], V1_channel[:2], V1_channel[3:]))
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
    #ofile.write("margphi=\n")
    ofile.write("progress=\n")
    ofile.write("H1-asd=%s\n"%(ligo_psd))
    ofile.write("L1-asd=%s\n"%(ligo_psd))
    ofile.write("V1-asd=%s\n"%(virgo_psd))
    ofile.write("randomseed=%d\n"%(randomseed)) # ### FIXME ###
    ofile.write("\n")
    
    ofile.write("[lalinference]\n")
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
      #os.system('echo %.2f > %s'%(start_time+seglen-2, os.path.join(out_folder, 'gpstime.txt')))
      #os.system('echo %.2f > %s'%(start_time+502.05, os.path.join(out_folder, 'gpstime.txt')))
      os.system('cp %s %s/gpstime.txt'%(gpstimefile, out_folder))
      os.system('%s/bin/lalinference_pipe -g %s/gpstime.txt -r %s -p %s %s/pipeconfig.ini'%(lal_prefix, out_folder, out_folder, out_folder, out_folder))
    if submit_dag:
      #dag_folder = os.path.join(out_folder, 'lalinferencenest', approx)
      dag_file = glob.glob(os.path.join(dag_folder, 'lalinference_*.dag'))[0]
      print('condor_submit_dag %s'%(dag_file))
      os.system('condor_submit_dag %s'%(dag_file))


# INJECTION LIST 
inj_file_txt = '/home/abhirup/Documents/Work/testGR_IR/waveforms/IHES_20161129/injection_list.txt'
inj_list = np.loadtxt(inj_file_txt)
inj_list = inj_list.astype(int)


for inj in inj_list[0:1]:
    frame_loc = '/home/abhirup/Documents/Work/testGR_IR/frames/IHES_20161129'
    H1_frame = glob.glob('%s/H-H1*_IHES_%04d_*.gwf'%(frame_loc, inj))[0] 
    L1_frame = glob.glob('%s/L-L1*_IHES_%04d_*.gwf'%(frame_loc, inj))[0] 
    V1_frame = glob.glob('%s/V-V1*_IHES_%04d_*.gwf'%(frame_loc, inj))[0] 
    gpstimefile = glob.glob('%s/*_IHES_%04d_*_gpstime.txt'%(frame_loc, inj))[0]

    print H1_frame, L1_frame, V1_frame

    noise = 'gaussian_noise'
    noise_type = 'LALAdLIGO'#'LIGO-P1200087-v18-aLIGO_EARLY_HIGH'
    date = '2016-11-29' # ###FIXME###
    run_tag = 'flow10Hz'#'realization_2' # ###FIXME###
    out_folder = '%s/Documents/Work/testGR_IR/runs/simulation_modGR/%s/%s/%s/IHES_%s_%s/injection_%04d'%(home, noise, noise_type, date, approximant, run_tag, inj)

    imr_folder = '%s/IMR'%(out_folder)
    insp_folder = '%s/inspiral'%(out_folder)
    post_insp_folder = '%s/post-inspiral'%(out_folder)
   
    #fit_formula = 'nospin_Pan2011'
    fit_formula = 'nonprecspin_Healy2014'

    # calculate the mass and spin of the final BH 
    Mf, af = tgr.calc_final_mass_spin(m1_inj[inj-1], m2_inj[inj-1], 0., 0., fit_formula)

    # calculate the Kerr ISCO freq 
    f_isco_Kerr = nr.calc_isco_freq(af)/(Mf*lal.MTSUN_SI)

    # calculate the dominant QNM freq 
    f_qnm = nr.calc_fqnm_dominant_mode(af)/(Mf*lal.MTSUN_SI)

    # calculating the appropriate seglen
    ctime_r = ct.calc_chirptime(m1_inj[inj-1], m2_inj[inj-1], f_min=fmin)
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

    insp_fhigh = f_isco_Kerr
    ring_flow = 0.8*f_qnm
    
    dataseed = np.random.randint(1,1000)
    randomseed = np.random.randint(1,1000)

    write_pipeline(imr_folder, H1_frame, L1_frame, V1_frame, gpstimefile, seglen, randomseed)
    write_pipeline(insp_folder, H1_frame, L1_frame, V1_frame, gpstimefile, seglen, randomseed, fhigh=insp_fhigh)
    write_pipeline(post_insp_folder, H1_frame, L1_frame, V1_frame, gpstimefile, seglen, randomseed, flow=insp_fhigh)
