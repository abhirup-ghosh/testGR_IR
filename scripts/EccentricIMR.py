#!/usr/bin/env python

import os, sys
from numpy import random, loadtxt, zeros, log2, ceil
from commands import getoutput
from pycbc.types import TimeSeries

MRSUN_SI = 1476.6250614046494
PC_SI    = 3.085677581491367e+16

def nextpow2( n ): return int(2**ceil(log2(n)))

EXE = {}
#CodeDir = '/projects/ncsa/grav/kuma/src/Eccentric_IMR/Codes/bin/'
CodeDir = '/home/abhirup/Documents/Work/testGR_IR/scripts/bin/'

# inspiral-merger-ringdown executables
#EXE[-1] = os.path.join(CodeDir, 'eccentric_waves_dynPN')
EXE[-1] = os.path.join(CodeDir, 'exc')

## inspiral-only executables # NOT RELATED TO ABOVE AND HAVE BUGS
#EXE[40] = os.path.join(CodeDir, 'eccentric_waves_2PN_finali_intel')
#EXE[50] = os.path.join(CodeDir, 'eccentric_waves_2-5PN_finali_intel')
#EXE[60] = os.path.join(CodeDir, 'eccentric_waves_3PN_finali_intel')
#EXE[70] = os.path.join(CodeDir, 'eccentric_waves_3-5PN_finali_intel')

# Checking that executables exist
for ex in EXE:
    if not os.path.exists(EXE[ex]):
        raise IOError(" Executable %s not found " % EXE[ex])


def generate_eccentric_waveform(\
                      mass1, mass2,\
                      ecc, anomaly, inclination, beta,\
                      tolerance,\
                      r_transition=5.,\
                      phase_order=-1,\
                      distance=1.,\
                      fmin=15,\
                      sample_rate=4096,\
                      length=None,\
                      inspiral_only=False,\
                      verbose=False):
  #{{{
  if phase_order not in EXE.keys():
    print "Only phase orders supported are: ", EXE.keys()
    raise IOError
  #
  tmp_outfile = '/dev/shm/tmp_%07d.dat' % int(1.e7 * random.random())
  if verbose: print "Using tmp file: ", tmp_outfile
  #
  if inspiral_only:
    print "INSPIRAL ONLY"
    tmpO = getoutput("%s %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %s" %\
                    (EXE[phase_order*10],\
                    mass1, mass2,\
                    ecc, anomaly, inclination, beta,\
                    tolerance,\
                    fmin, float(sample_rate),\
                    tmp_outfile))
  elif phase_order == -1:
    CMD = "%s -m %.12e -n %.12e -e %.12e -a %.12e -i %.12e -b %.12e -t %.12e -f %.12e -s %.12e -o %s" %\
                    (EXE[phase_order],\
                    mass1, mass2,\
                    ecc, anomaly, inclination, beta,\
                    tolerance,\
                    fmin, float(sample_rate),\
                    tmp_outfile)
    if verbose: CMD = CMD + " -v"
  elif phase_order >= 6 and phase_order <= 12:
    CMD = "%s -R %d -m %.12e -n %.12e -e %.12e -a %.12e -i %.12e -b %.12e -t %.12e -f %.12e -s %.12e -o %s" %\
                    (EXE[-3], phase_order,\
                    mass1, mass2,\
                    ecc, anomaly, inclination, beta,\
                    tolerance,\
                    fmin, float(sample_rate),\
                    tmp_outfile)
    if verbose: CMD = CMD + " -v"
  else:
    raise IOError("Which code version do you want to use??")
    print "EXPLICIT PHASE PN ORDER SET"
    tmpO = getoutput("%s %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %s" %\
                    (EXE[phase_order],\
                    mass1, mass2,\
                    ecc, anomaly, inclination, beta,\
                    tolerance,\
                    fmin, float(sample_rate),\
                    r_transition, tmp_outfile))
  #
  print >>sys.stdout, "\tCommand to be run: --- \n", CMD
  sys.stdout.flush()
  tmpO = getoutput( CMD )
  #
  #amp0 = (mass1+mass2) * MRSUN_SI / (distance * 1.e6 * PC_SI)
  eta = mass1 * mass2 / (mass1 + mass2)**2.
  DYN_RANGE_FAC = 5.9029581035870565e+20
  print eta, DYN_RANGE_FAC
  amp0 = 1. / DYN_RANGE_FAC
  amp0 *= 100. / distance ## because eccentric code is defaulted to 100Mpc instead of 1Mpc
  #
  d = loadtxt(tmp_outfile)
  times = d[:,0]
  hplus = d[:,1]
  hcross= d[:,2]
  #
  if length is not None and length >= len(times): N = length
  else: N = len(times)
  hp = TimeSeries( zeros(N), delta_t=1./sample_rate )
  hc = TimeSeries( zeros(N), delta_t=1./sample_rate )
  #
  hp.data[:len(hplus)]  = amp0 * hplus
  hc.data[:len(hcross)] = amp0 * hcross
  #
  getoutput('rm -f %s' % tmp_outfile)
  #
  return hp, hc
  #}}}


