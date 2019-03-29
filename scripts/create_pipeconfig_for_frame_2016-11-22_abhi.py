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
import optparse as op
import chirptimes as ct
import time
from time import sleep
import glob
from write_cache import write_cache

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

framesuffix = 'IHES_%04d'%inj

# ### FRAME LOCATION FOR Nathan's frames ###
frame_loc = '../frames/IHES_20161130'
out_folder = '/home/abhirup.ghosh/Documents/Work/testGR_IR/runs/paper2_20170417_H1L1/GR_IHES/injection_%04d'%inj
H1_frame = glob.glob(os.path.join(frame_loc, 'H*%s*.gwf'%(framesuffix)))[0]
L1_frame = glob.glob(os.path.join(frame_loc, 'L*%s*.gwf'%(framesuffix)))[0]
V1_frame = glob.glob(os.path.join(frame_loc, 'V*%s*.gwf'%(framesuffix)))[0]
gpstimefile = glob.glob(os.path.join(frame_loc, '*%s*_gpstime.txt'%(framesuffix)))[0]

print "H1_frame:\t%s"%(H1_frame)
print "L1_frame:\t%s"%(L1_frame)
print "V1_frame:\t%s"%(V1_frame)
print "out_folder:\t%s"%(out_folder)

try:
  os.system('ls %s'%(H1_frame))
  os.system('ls %s'%(L1_frame))
  os.system('ls %s'%(V1_frame))
except:
  print "Frame not found .. exiting."
  sys.exit()

os.system('mkdir -p %s'%(out_folder))

# ### RECOVERY ###

approximant = 'SEOBNRv4_ROMpseudoFourPN' # recovery waveform
ligo_psd = os.path.join(psd_path, 'LIGO-T0900288-v3-ZERO_DET_high_P.txt')
virgo_psd = os.path.join(psd_path, 'LIGO-P1200087-v18-AdV_DESIGN.txt')


def write_pipeline(out_folder, flow=20., fhigh=1023., nlive=1024, create_dag=True, submit_dag=True):
  
    randomseed = np.random.randint(32768)
    
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
    ofile.write("ifos=['H1','L1']\n")
    ofile.write("engine=lalinferencenest\n")
    ofile.write("nparallel=4\n")
    ofile.write("coherence-test=False\n")
    ofile.write("accounting_group=ligo.prod.o2.cbc.pe.lalinference\n")
    ofile.write("upload-to-gracedb=False\n")
    ofile.write("\n")
    
    projectdir = out_folder.replace('%s/Documents/Work/'%home, '')

    ofile.write("[paths]\n")
    ofile.write("webdir=%s/public_html/%s\n"%(home, projectdir))
    ofile.write("basedir=%s/public_html/%s\n"%(home, projectdir))
    ofile.write("baseurl=https://www.atlas.aei.uni-hannover.de/~%s/%s\n"%(user, projectdir))
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
    ofile.write("comp-min=3.0\n") 
    ofile.write("comp-max=300.0\n") 
    ofile.write("a_spin1-max=0.98\n") # ###FIXME###
    ofile.write("a_spin2-max=0.98\n") # ###FIXME###
    ofile.write("amporder=-1\n")
    ofile.write("resume=\n")
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
      os.system('cp %s %s/gpstime.txt'%(gpstimefile, out_folder))
      os.system('%s/bin/lalinference_pipe -g %s/gpstime.txt -r %s -p %s %s/pipeconfig.ini'%(lal_prefix, out_folder, out_folder, out_folder, out_folder))
    if submit_dag:
      dag_file = glob.glob(os.path.join(dag_folder, 'lalinference_*.dag'))[0]
      print('condor_submit_dag %s'%(dag_file))
      os.system('condor_submit_dag %s'%(dag_file))

# Main
if __name__ == '__main__':

    imr_folder = '%s/IMR'%(out_folder)
    insp_folder = '%s/inspiral'%(out_folder)
    post_insp_folder = '%s/post-inspiral'%(out_folder)

    insp_fhigh = ring_flow = fcut
    
    write_pipeline(imr_folder)
    write_pipeline(insp_folder, fhigh=insp_fhigh)
    write_pipeline(post_insp_folder, flow=ring_flow)
