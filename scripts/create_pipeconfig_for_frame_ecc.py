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
from optparse import OptionParser


os.system('. ${HOME}/.profile')
home = os.getenv('HOME')
user = commands.getoutput('whoami')
lal_prefix = os.getenv('LAL_PREFIX')
psd_path = '%s/share/lalsimulation'%lal_prefix
glue_location = os.getenv('GLUE_LOCATION')
pylal_location = os.getenv('PYLAL_LOCATION')


approximant = 'EccentricTDIMRv2'
#approximant = 'SEOBNRv2_ROM_DoubleSpinpseudoFourPN'
#approximant = 'IMRPhenomPv2'
seglen = 16

# ### RECOVERY ###
ligo_psd = os.path.join(psd_path, 'LIGO-T0900288-v3-ZERO_DET_high_P.txt')
virgo_psd = os.path.join(psd_path, 'LIGO-T0900288-v3-ZERO_DET_high_P.txt')

def write_pipeline(out_folder, H1_frame, L1_frame, randomseed, flow=15., fhigh=1023., nlive=1024, create_dag=True, submit_dag=True):
  
    cache_folder = os.path.join(out_folder, 'lalinferencenest/%s/caches'%approximant)
    os.system('mkdir -p %s'%(cache_folder))
    
    print cache_folder, H1_frame

    (H1_cache,H,H1_channel,start_time,datalen,Hpath) = write_cache(H1_frame, cache_folder).split()
    (L1_cache,L,L1_channel,start_time,datalen,Lpath) = write_cache(L1_frame, cache_folder).split()
    
    start_time = int(start_time)
    datalen = int(datalen)
    print "Start time = %d, datalen = %d"%(start_time, datalen)
    ofile = open(os.path.join(out_folder, 'pipeconfig.ini'), 'w')

    ofile.write("[analysis]\n")
    ofile.write("ifos=['H1','L1']\n")
    ofile.write("engine=lalinferencenest\n")
    ofile.write("nparallel=4\n")
    ofile.write("coherence-test=False\n")
    ofile.write("upload-to-gracedb=False\n")
    ofile.write("\n")
    
    projectdir = out_folder.replace('%s/Documents/Work/'%home, '')

    ofile.write("[paths]\n")
    ofile.write("webdir=%s/public_html/%s\n"%(home, projectdir))
    ofile.write("basedir=%s/public_html/%s\n"%(home, projectdir))
    ofile.write("baseurl=https://dogmatix.icts.res.in/~%s/%s\n"%(user, projectdir))
    #ofile.write("baseurl=https://dumpty.alice.icts.res.in/~%s/%s\n"%(user, projectdir))
    ofile.write("\n")
    
    ofile.write("[input]\n")
    ofile.write("max-psd-length=1024\n")
    ofile.write("padding=16\n")
    ofile.write("timeslides=false\n")
    ofile.write("ignore-science-segments=True\n")
    ofile.write("gps-start-time=%d\n"%(start_time))
    ofile.write("gps-end-time=%d\n"%(start_time+datalen+16))
    ofile.write("psd-length=%d\n"%(datalen/2))
    #ofile.write("events=all\n")
    ofile.write("\n")
    
    ofile.write("[datafind]\n")
    ofile.write("types={'H1':'%s','L1':'%s'}\n"%(H1_channel[3:], L1_channel[3:]))
    ofile.write("\n")
    
    ofile.write("[data]\n")
    ofile.write("channels={'H1':'%s:%s','L1':'%s:%s'}\n"%(H1_channel[:2], H1_channel[3:], L1_channel[:2], L1_channel[3:]))
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
    #ofile.write("mpirun=/usr/bin/mpirun\n")
    #ofile.write("mpiwrapper=%sbin/lalinference_mpi_wrapper\n"%(lal_prefix))
    ofile.write("\n")
    
    ofile.write("[engine]\n")
    ofile.write("approx=%s\n"%(approximant))
    ofile.write("seglen=%d\n"%(seglen))
    ofile.write("nlive=%d\n"%(nlive))
    ofile.write("srate=2048\n")
    ofile.write("distance-min=1\n") # ###FIXME###
    ofile.write("distance-max=20000\n") # ###FIXME###
    ofile.write("comp-min=10.0\n") 
    ofile.write("comp-max=80.0\n") 
    #ofile.write("chirpmass-min = 5.0\n")
    #ofile.write("chirpmass-max = 20.0\n")
    #ofile.write("q-min = 0.0555555555555556\n")
    #ofile.write("q-max = 1.0\n")
    ofile.write("a_spin1-max=0.99\n") # ###FIXME###
    ofile.write("a_spin2-max=0.99\n") # ###FIXME###
    #ofile.write("disable-spin=\n") # ### FIXME ###
    #ofile.write("aligned-spin=\n") # ### FIXME ###
    ofile.write("amporder=-1\n")
    ofile.write("margphi=\n") # ###FIXME###
    ofile.write("resume=\n")
    ofile.write("progress=\n")
    ofile.write("H1-psd=%s\n"%(ligo_psd))
    ofile.write("L1-psd=%s\n"%(ligo_psd))
    ofile.write("randomseed=%d\n"%(randomseed)) # ### FIXME ###
    ofile.write("fix-lambda2=0\n")
    ofile.write("lambda1-min=0.09\n")
    ofile.write("lambda1-max=0.11\n")
    ofile.write("\n")
    
    ofile.write("[lalinference]\n")
    #ofile.write("cache={'H1':'%s','L1':'%s'}\n"%(H1_cache, L1_cache))
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
      os.system('echo %.2f > %s'%(start_time+seglen-2, os.path.join(out_folder, 'gpstime.txt')))
      os.system('%s/bin/lalinference_pipe -g %s/gpstime.txt -r %s -p %s %s/pipeconfig.ini'%(lal_prefix, out_folder, out_folder, out_folder, out_folder))
    if submit_dag:
      os.system('/usr/bin/condor_submit_dag %s/multidag.dag'%out_folder)



#--- MAIN ----
    
frame_loc = '/home/abhirup/Documents/Work/EccentricTDIMR/frames/20160922_nonoise_EccentricTDIMRv2'
H1_frame = glob.glob('%s/H-H1_modGR_INJECTED-1126285202-16_EccentricTDIMRv2_m1_36_m2_29_f_min_10_emin_0p1_d_410_M65_dist410_incl0p0_ra0p0_dec0p0_psi0p0_flow10_nonoise.gwf'%frame_loc)[0] 
L1_frame = glob.glob('%s/L-L1_modGR_INJECTED-1126285202-16_EccentricTDIMRv2_m1_36_m2_29_f_min_10_emin_0p1_d_410_M65_dist410_incl0p0_ra0p0_dec0p0_psi0p0_flow10_nonoise.gwf'%frame_loc)[0] 

print H1_frame, L1_frame

run_tag = '20160922_nonoise_EccentricTDIMRv2_framefiles'

imr_folder = '%s/Documents/Work/EccentricTDIMR/runs/%s/Ecc_emin_0p1_inj_Ecc_rec_spins_prec_axmodelhack'%(home, run_tag)
    
randomseed = np.random.randint(1,1000)

write_pipeline(imr_folder, H1_frame, L1_frame, randomseed)
