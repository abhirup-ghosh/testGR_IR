#!/usr/bin/env python

"""
This script creates an injection and a .ini file to run LALInference

(C) Archisman Ghosh, Abhirup Ghosh, 2015-09-20
"""

import sys, lal, os, commands, numpy as np
sys.path.insert(1, os.path.join(os.path.expanduser('~'), 'src/lalsuite/lalinference/python'))
import imrtestgr as tgr
import nr_fits as nr
import optparse as op
import chirptimes as ct
import time
from time import sleep
import glob

os.system('. ${HOME}/.profile')
home = os.getenv('HOME')
user = commands.getoutput('whoami')
lal_prefix = os.getenv('LAL_PREFIX')
psd_path = '%s/share/lalsimulation'%lal_prefix
glue_location = os.getenv('GLUE_LOCATION')
pylal_location = os.getenv('PYLAL_LOCATION')

parser = op.OptionParser()
parser.add_option('--inj-no', type='int', dest='inj', help='injection number')
parser.add_option('--fcut', type='float', dest='fcut', help='fcut')
(options, args) = parser.parse_args()
inj = options.inj
fcut = options.fcut

inj_file = '/home/abhirup/Documents/Work/testGR_IR/runs/LIGO_software_injections_HLV/IMRPhenomPv2_inj_SEOBNRv4_ROM_rec/software_inj_IMRPhenomPv2pseudoFourPN_m_10_50_z_0_0p5'
out_folder = '/home/abhirup/Documents/Work/testGR_IR/runs/LIGO_software_injections_HLV/IMRPhenomPv2_inj_SEOBNRv4_ROM_rec/injection_%04d'%inj

# ### RECOVERY ###

approximant = 'SEOBNRv4_ROMpseudoFourPN' # recovery waveform

def write_pipeline(out_folder, flow=20., fhigh=1024., nlive=1024, create_dag=True):

    randomseed = np.random.randint(32768)

    os.system('mkdir -p %s'%(out_folder))
    os.system('cp %s.xml %s/injection.xml'%(inj_file, out_folder))

    ofile = open(os.path.join(out_folder, 'pipeconfig.ini'), 'w')
    
    ofile.write("[analysis]\n")
    ofile.write("ifos=['H1','L1','V1']\n")
    ofile.write("engine=lalinferencenest\n")
    ofile.write("nparallel=4\n")
    ofile.write("coherence-test=True\n")
    ofile.write("accounting_group=ligo.prod.o2.cbc.pe.lalinference\n")
    ofile.write("upload-to-gracedb=False\n")
    ofile.write("add-lvem-tag = False\n")
    ofile.write("\n")
    
    projectdir = out_folder.replace('%s/Documents/Work/'%home, '')

    ofile.write("[paths]\n")
    ofile.write("webdir=%s/public_html/%s\n"%(home, projectdir))
    ofile.write("baseurl=https://dumpty.alice.icts.res.in/~%s/%s\n"%(user, projectdir))
    ofile.write("\n")
   
    ofile.write("[input]\n")
    ofile.write("max-psd-length=1024\n")
    ofile.write("padding=16\n")
    ofile.write("timeslides=false\n")
    ofile.write("ignore-science-segments=True\n")
    ofile.write("events=%d\n"%inj)
    ofile.write("analyse-all-time = False\n")
    ofile.write("ignore-science-segments = True\n")
    ofile.write("ignore-gracedb-psd = True\n")
    ofile.write("\n")
    
    ofile.write("[datafind]\n")
    ofile.write("url-type = file\n")
    ofile.write("types = {'H1':'H1_CLEANED_HOFT_C02','L1':'L1_CLEANED_HOFT_C02','V1':'V1O2Repro2A'}\n")
    ofile.write("\n")
    
    ofile.write("[data]\n")
    ofile.write("channels = {'H1':'H1:DCH-CLEAN_STRAIN_C02','L1':'L1:DCH-CLEAN_STRAIN_C02','V1':'V1:Hrec_hoft_V1O2Repro2A_16384Hz'}\n")
    ofile.write("\n")
    
    ofile.write("[condor]\n")
    ofile.write("mpirun=/bin/true\n")
    ofile.write("datafind=%s/bin/gw_data_find\n"%(glue_location))
    ofile.write("gracedb=/bin/true\n")
    ofile.write("ligolw_print=%s/bin/ligolw_print\n"%(glue_location))
    ofile.write("segfind=%s/bin/ligolw_segment_query\n"%(glue_location))
    ofile.write("lalinferencenest=%s/bin/lalinference_nest\n"%(lal_prefix))
    ofile.write("lalinferencemcmc=%s/bin/lalinference_mcmc\n"%(lal_prefix))
    ofile.write("mergescript=%s/bin/lalapps_nest2pos\n"%(lal_prefix))
    ofile.write("resultspage=%s/bin/cbcBayesPostProc.py\n"%(pylal_location))
    ofile.write("ppanalysis=%s/bin/cbcBayesPPAnalysis.py\n"%(pylal_location))
    ofile.write("pos_to_sim_inspiral=%s/bin/cbcBayesPosToSimInspiral.py\n"%(pylal_location))
    ofile.write("coherencetest=%s/bin/lalapps_coherence_test\n"%(lal_prefix))
    ofile.write("lalinferencedatadump=False\n")
    ofile.write("\n")
    
    ofile.write("[engine]\n")
    ofile.write("approx=%s\n"%(approximant))
    ofile.write("seglen=8\n")
    ofile.write("nlive=%d\n"%(nlive))
    ofile.write("srate=2048\n")
    ofile.write("distance-min=1\n") 
    ofile.write("distance-max=20000\n") 
    ofile.write("comp-min=10.0\n") 
    ofile.write("comp-max=300.0\n") 
    ofile.write("a_spin1-max=0.99\n")
    ofile.write("a_spin2-max=0.99\n")
    ofile.write("amporder=-1\n") ##FIXME##
    ofile.write("resume=\n")
    ofile.write("progress=\n")
    ofile.write("randomseed=%d\n"%(randomseed))
    ofile.write("alignedspin-zprior=\n")
    ofile.write("\n")
    
    ofile.write("[lalinference]\n")
    ofile.write("flow={'H1':%.2f,'L1':%.2f,'V1':%.2f}\n"%(flow, flow, flow))
    ofile.write("fhigh={'H1':%.2f,'L1':%.2f,'V1':%.2f}\n"%(fhigh, fhigh, fhigh))
    ofile.write("\n")
    
    ofile.write("[segfind]\n")
    ofile.write("segment-url = https://segments.ligo.org\n")
    ofile.write("\n")

    ofile.write("[segments]\n")
    ofile.write("h1-analyze = H1:DCH-CLEAN_SCIENCE_C02:1\n")
    ofile.write("l1-analyze = L1:DCH-CLEAN_SCIENCE_C02:1\n")
    ofile.write("v1-analyze = V1:ITF_SCIENCEMODE:7\n")
    ofile.write("\n")

    ofile.write("[resultspage]\n")
    ofile.write("skyres=0.5\n")
    ofile.write("\n")
    
    ofile.close()

    if create_dag:
      os.system('%s/bin/lalinference_pipe -I %s/injection.xml -r %s -p %s %s/pipeconfig.ini'%(lal_prefix, out_folder, out_folder, out_folder, out_folder))
      print('condor_submit_dag %s/multidag.dag'%(out_folder))
      os.system('condor_submit_dag %s/multidag.dag'%(out_folder))


# Main
if __name__ == '__main__':

    imr_folder = '%s/IMR'%(out_folder)
    insp_folder = '%s/inspiral'%(out_folder)
    post_insp_folder = '%s/post-inspiral'%(out_folder)

    insp_fhigh = ring_flow = fcut

    write_pipeline(imr_folder)
    write_pipeline(insp_folder, fhigh=insp_fhigh)
    write_pipeline(post_insp_folder, flow=ring_flow)

