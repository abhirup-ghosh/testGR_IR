#!/usr/bin/env python

"""
This script creates an injection and a .ini file to run LALInference

(C) Archisman Ghosh, Abhirup Ghosh, 2015-09-20
"""

import sys, lal, os, glob, re, commands, numpy as np
from scipy import random
from write_cache import write_cache
import matplotlib.pyplot as plt
from optparse import OptionParser

date_tag = '2016-11-29'

parser = OptionParser()
parser.add_option('-f', '--frame', type='string', dest='framesuffix', help='frame suffix')
(options, args) = parser.parse_args()
framesuffix = options.framesuffix

os.system('. ${HOME}/.profile')
home = os.getenv('HOME')
user = commands.getoutput('whoami')
if user in ['abhirup', 'abhirup.ghosh']:
  work_dir = os.path.join(home, 'Documents/Work')
elif user in ['archis', 'archisman.ghosh']:
  work_dir = os.path.join(home, 'Work')

#psd_path = os.path.join(home, 'src/lalsuite/lalsimulation/src')
lal_prefix = os.getenv('LAL_PREFIX')
glue_location = os.getenv('GLUE_LOCATION')
pylal_location = os.getenv('PYLAL_LOCATION')
psd_path = os.path.realpath(os.path.join(lal_prefix, 'share', 'lalsimulation'))

approx = 'SEOBNRv2_ROM_DoubleSpinthreePointFivePN'
#approx = 'IMRPhenomPv2threePointFivePN'
seglen = 32 #4
#seglen = 8

# ### FRAME LOCATION FOR Nathan's frames ###
frame_loc = os.path.join(work_dir, 'testGR_IR', 'NR_frames', date_tag)
print frame_loc
#out_folder = os.path.join(frame_loc.replace('frames', 'runs/simulations_modGR'), date_tag, framesuffix, '%s_seglen%d_nomargphi'%(approx, seglen))
out_folder = os.path.join(os.path.realpath('../runs/simulations_NR'), date_tag, framesuffix, '%s_seglen%d_nomargphi_30Hz'%(approx, seglen))
H1_frame = glob.glob(os.path.join(frame_loc, 'H*%s*.gwf'%(framesuffix)))[0]
L1_frame = glob.glob(os.path.join(frame_loc, 'L*%s*.gwf'%(framesuffix)))[0]
V1_frame = glob.glob(os.path.join(frame_loc, 'V*%s*.gwf'%(framesuffix)))[0]
gpstimefile = glob.glob(os.path.join(frame_loc, '*%s*_gpstime.txt'%(framesuffix)))[0]

print gpstimefile

## ### FRAME LOCATION FOR Patricia's frames ###
#frame_loc = os.path.join(work_dir, 'WaveformSystematics', 'frames')
#out_folder = os.path.join(frame_loc.replace('frames', 'imrtgr_runs'), date_tag, framesuffix, '%s_vanilla_seglen%s'%(approx, seglen))
#if framesuffix[-7:]=='_0noise':
#  foldersuffix = ''
#  filesuffix = ''
#elif framesuffix[-7:]=='_noise1':
#  foldersuffix = ''
#  filesuffix = '_C00'
#elif framesuffix[-7:-1]=='_noise':
#  foldersuffix = framesuffix[-7:]
#  filesuffix = '_C00'
#else:
#  print "Invalid framesuffix!"
#  sys.exit()
#if framesuffix[-12:-7]=='_l2m2':
#  H1_frame = os.path.join(frame_loc, framesuffix[:-12]+foldersuffix, 'H1HWINJ_config%s_l2m2%s.gwf'%(framesuffix[-16:-12], filesuffix))
#  L1_frame = os.path.join(frame_loc, framesuffix[:-12]+foldersuffix, 'L1HWINJ_config%s_l2m2%s.gwf'%(framesuffix[-16:-12], filesuffix))
#elif framesuffix[-10:-7]=='_HM':
#  H1_frame = os.path.join(frame_loc, framesuffix[:-7]+foldersuffix, 'H1HWINJ_config%s%s.gwf'%(framesuffix[-14:-7], filesuffix))
#  L1_frame = os.path.join(frame_loc, framesuffix[:-7]+foldersuffix, 'L1HWINJ_config%s%s.gwf'%(framesuffix[-14:-7], filesuffix))
#else:
#  print "Invalid framesuffix!"
#  sys.exit()
  
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

#ligo_psd = os.path.join(psd_path, 'LIGO-P1200087-v18-aLIGO_EARLY_HIGH.txt')
#virgo_psd = os.path.join(psd_path, 'LIGO-P1200087-v18-AdV_EARLY_HIGH.txt')
ligo_psd = os.path.join(psd_path, 'LIGO-T0900288-v3-ZERO_DET_high_P.txt')
virgo_psd = os.path.join(psd_path, 'LIGO-T0900288-v3-ZERO_DET_high_P.txt')


def write_pipeline(out_folder, flow=30., fhigh=1023., nlive=1024, create_dag=True, submit_dag=True):
  
    randomseed = random.randint(32768)
    
    dag_folder = os.path.join(out_folder, 'lalinferencenest', approx)
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
    #ofile.write("dataseed=1\n") # ###FIXME###
    #ofile.write("accounting_group=ligo.sim.o1.cbc.testgr.tiger\n")
    ofile.write("\n")
    
    projectdir = out_folder.replace('%s/'%work_dir, '')

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
    #ofile.write("events=all\n")
    ofile.write("\n")
    
    ofile.write("[datafind]\n")
    #ofile.write("types={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo','I1':'LALSimAdLIGO','J1':'LALSimAdLIGO'}\n")
    ofile.write("types={'H1':'%s','L1':'%s','V1':'%s'}\n"%(H1_channel[3:], L1_channel[3:], V1_channel[3:]))
    ofile.write("\n")
    
    ofile.write("[data]\n")
    #ofile.write("channels={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo','I1':'LALSimAdLIGO','J1':'LALSimAdLIGO'}\n")
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
    #ofile.write("resultspage=/home/abhirup/opt/lalsuite/pylal/bin/cbcBayesPostProc.py\n")
    ofile.write("coherencetest=%s/bin/lalapps_coherence_test\n"%(lal_prefix))
    ofile.write("lalinferencedatadump=False\n")
    ofile.write("\n")
    
    ofile.write("[engine]\n")
    ofile.write("approx=%s\n"%(approx))
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
    ofile.write("H1-asd=%s\n"%(ligo_psd))
    ofile.write("L1-asd=%s\n"%(ligo_psd))
    ofile.write("V1-asd=%s\n"%(ligo_psd))
    ofile.write("randomseed=%d\n"%(randomseed)) # ### FIXME ###
    #ofile.write("0noise=\n") 
    #ofile.write("tolerance=0.05\n")
    #ofile.write("pinparams=[time,phase,costheta_jn,distance,rightascension,declination,polarisation]\n") # ###FIXME####
    ofile.write("\n")
    
    ofile.write("[lalinference]\n")
    #ofile.write("fake-cache={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo','I1':'LALSimAdLIGO','J1':'LALSimAdLIGO'}\n")
    #ofile.write("cache={'H1':'%s','L1':'%s'}\n"%(H1_cache, L1_cache))
    #ofile.write("fake-cache={'H1':'interp:%s','L1':'interp:%s','V1':'interp:%s'}\n"%(ligo_psd, ligo_psd, virgo_psd))
    ofile.write("flow={'H1':%.2f,'L1':%.2f,'V1':%.2f}\n"%(flow, flow, flow))
    ofile.write("fhigh={'H1':%.2f,'L1':%.2f,'V1':%.2f}\n"%(fhigh, fhigh, fhigh))
    ofile.write("\n")
    
    ofile.write("[merge]\n")
    #ofile.write("npos=5000\n")
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

# Main
if __name__ == '__main__':

    imr_folder = '%s/IMR'%(out_folder)
    insp_folder = '%s/inspiral'%(out_folder)
    post_insp_folder = '%s/post-inspiral'%(out_folder)

    insp_fhigh = ring_flow = 66.960890129648334 #95.7 #134.
    
    write_pipeline(imr_folder)
    write_pipeline(insp_folder, fhigh=insp_fhigh)
    write_pipeline(post_insp_folder, flow=ring_flow)
    #write_pipeline(ring_folder, M, q, chi1, chi2, comp_mass_max=1.5*m2, flow=ring_flow)
