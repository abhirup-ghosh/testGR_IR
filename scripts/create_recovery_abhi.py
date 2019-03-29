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
parser.add_option('--seglen', type='int', dest='seglen', help='seglen')
(options, args) = parser.parse_args()
inj = options.inj
fcut = options.fcut
seglen = options.seglen

inj_file = '../injections/waveform_systematics/waveform_systematics_investigations_SEOBNRv4_ROMpseudoFourPN_m_10_80_z_0_0p5_UCTV'
out_folder = '/home/abhirup/Documents/Work/testGR_IR/runs/waveform_systematics/20170427_SEOBNRv4_ROM_inj_SEOBNRv4_ROM_rec_H1L1V1_flow15Hz/injection_%04d'%inj

# ### RECOVERY ###

approximant = 'SEOBNRv4_ROMpseudoFourPN' # recovery waveform
ligo_psd = os.path.join(psd_path, 'LIGO-T0900288-v3-ZERO_DET_high_P.txt')
virgo_psd = os.path.join(psd_path, 'LIGO-P1200087-v18-AdV_DESIGN.txt')


def write_pipeline(out_folder, dataseed, flow=15, fhigh=1023., nlive=1024, create_dag=True):

    randomseed = np.random.randint(32768)

    os.system('mkdir -p %s'%(out_folder))
    os.system('cp %s.xml %s/injection.xml'%(inj_file, out_folder))

    ofile = open(os.path.join(out_folder, 'pipeconfig.ini'), 'w')
    
    ofile.write("[analysis]\n")
    ofile.write("ifos=['H1','L1','V1']\n")
    ofile.write("engine=lalinferencenest\n")
    ofile.write("nparallel=4\n")
    ofile.write("coherence-test=False\n")
    ofile.write("accounting_group=ligo.prod.o2.cbc.pe.lalinference\n")
    ofile.write("upload-to-gracedb=False\n")
    ofile.write("dataseed=%d\n"%dataseed) # ###FIXME###
    ofile.write("\n")
    
    projectdir = out_folder.replace('%s/Documents/Work/'%home, '')

    ofile.write("[paths]\n")
    ofile.write("webdir=%s/public_html/%s\n"%(home, projectdir))
    ofile.write("basedir=%s/public_html/%s\n"%(home, projectdir))
    ofile.write("baseurl=https://dumpty.alice.icts.res.in/~%s/%s\n"%(user, projectdir))
    ofile.write("\n")
    
    ofile.write("[input]\n")
    ofile.write("max-psd-length=1024\n")
    ofile.write("padding=16\n")
    ofile.write("timeslides=false\n")
    ofile.write("ignore-science-segments=True\n")
    ofile.write("events=%d\n"%(int(inj-1)))
    ofile.write("\n")
    
    ofile.write("[datafind]\n")
    ofile.write("types={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo'}\n")
    ofile.write("\n")
    
    ofile.write("[data]\n")
    ofile.write("channels={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo'}\n")
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
    ofile.write("skyarea=/usr/bin/run_sky_area\n")
    ofile.write("\n")
 
    ofile.write("[engine]\n")
    ofile.write("approx=%s\n"%(approximant))
    ofile.write("seglen=%d\n"%(seglen))
    ofile.write("nlive=%d\n"%(nlive))
    ofile.write("srate=2048\n")
    ofile.write("distance-min=1\n") 
    ofile.write("distance-max=20000\n") 
    ofile.write("comp-min=3.0\n") 
    ofile.write("comp-max=300.0\n") 
    ofile.write("a_spin1-max=0.99\n")
    ofile.write("a_spin2-max=0.99\n")
    ofile.write("amporder=-1\n") ##FIXME##
    ofile.write("resume=\n")
    ofile.write("progress=\n")
    ofile.write("H1-asd=%s\n"%(ligo_psd))
    ofile.write("L1-asd=%s\n"%(ligo_psd))
    ofile.write("V1-asd=%s\n"%(virgo_psd))
    ofile.write("randomseed=%d\n"%(randomseed))
    ofile.write("\n")
    
    ofile.write("[lalinference]\n")
    ofile.write("fake-cache={'H1':'interp:%s','L1':'interp:%s','V1':'interp:%s'}\n"%(ligo_psd, ligo_psd, virgo_psd))
    #ofile.write("fake-cache={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo'}\n")
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
      dag_folder = os.path.join(out_folder, 'lalinferencenest', approximant)
      dag_file = glob.glob(os.path.join(dag_folder, 'lalinference_*.dag'))[0]
      print('condor_submit_dag %s'%(dag_file))
      os.system('condor_submit_dag %s'%(dag_file))


# Main
if __name__ == '__main__':

    imr_folder = '%s/IMR'%(out_folder)
    insp_folder = '%s/inspiral'%(out_folder)
    post_insp_folder = '%s/post-inspiral'%(out_folder)

    insp_fhigh = ring_flow = fcut
    dataseed = np.random.randint(1,1000)

    write_pipeline(imr_folder, dataseed)
    write_pipeline(insp_folder, dataseed, fhigh=insp_fhigh)
    write_pipeline(post_insp_folder, dataseed, flow=ring_flow)


