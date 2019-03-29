#!/usr/bin/env python

import re, os, glob
from write_cache import write_cache
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--observatory', type='string')
parser.add_option('--url-type', type='string')
parser.add_option('--gps-start-time', type='string')
parser.add_option('--gps-end-time', type='string')
parser.add_option('--output', type='string')
parser.add_option('--lal-cache', action='store_true', default=False)
parser.add_option('--type', type='string', dest='channel')
parser.add_option('--identifier', type='string', help='unique identifier for frame', default='')
parser.add_option('--location', type='string', help='location (directory) of frames')
(options, args) = parser.parse_args()

observatory = options.observatory
identifier = options.identifier
location = options.location
channel = options.channel
output = options.output

frame = glob.glob(os.path.join(location, '%s*%s*.gwf'%(observatory, identifier)))[0]

print "%s frame: \t%s"%(observatory, frame)

try:
  os.system('ls %s'%(frame))
except:
  print "Frame not found .. exiting."
  sys.exit()

cache_folder = os.path.dirname(output)
cachefilename = os.path.basename(output)
datalen = int(re.findall('[0-9]*.gwf', os.path.basename(frame))[0][:-4])

(cache,obs,channel,start_time,datalen,Hpath) = write_cache(frame, cache_folder, cachefilename, channel, datalen, '_').split()

print "Written: \t%s"%(os.path.join(cache_folder, cachefilename))
