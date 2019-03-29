#!/home/haris.k/devel/VirtualPycbc_Oct2017/bin/python

'''
Script to inject a signal to a frame file 
'''
import argparse
import sys
import numpy as np
from pycbc.frame import read_frame, write_frame
from pylal import  Fr
import logging
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Script to inject a signal to a frame file')
parser.add_argument('-i_framefile', '--i_framefile', help='input frame file', required=True)
#parser.add_argument('-o_file_tag', '--o_file_tag', help='tag for output frame file',default='with_inj')
parser.add_argument('-channel', '--channel', help='channel name', required=True)
parser.add_argument('-inj_file', '--inj_file', help='txt file containing inj signal time series', required=True)
parser.add_argument('-inj_start_time', '--inj_start_time', help='injection start time',type=int, required=True)
parser.add_argument('-inj_sample_rate', '--inj_sample_rate', help='injection sample rate', type=int, default=16384)



args = parser.parse_args()
i_framefile = args.i_framefile
channel = args.channel
#o_file_tag = args.o_file_tag
inj_file = args.inj_file
inj_start_time = args.inj_start_time
inj_sample_rate = args.inj_sample_rate

frdata, gps_start, xoffset, xspacing, xunit, yunit = Fr.frgetvect(i_framefile, channel)
gps_end=gps_start+xspacing[0]*len(frdata)
t = np.linspace(gps_start,gps_end,len(frdata))
if (xspacing[0]!=1./inj_sample_rate):
    print 'Time staps of input and output frame files are not matching'
    sys.exit(1)

logging.info('Reading frame file...')
i_data = read_frame(i_framefile, channel, gps_start, gps_end)
time = i_data.get_sample_times().data

logging.info('Reading injection file...')
inj_data=np.loadtxt(inj_file)
start_idx = (np.abs(time - inj_start_time)).argmin()

inj_data_padded = np.pad(inj_data,(start_idx,len(time)-(start_idx+len(inj_data))), 'constant', constant_values=(0, 0))
i_data.data = i_data.data + inj_data_padded
o_file_tag = inj_file.split('.txt')[0]
o_framefile = i_framefile.split('.gwf')[0].split('/')[-1]+'_'+o_file_tag+'.gwf'
write_frame(o_framefile, channel, i_data)
logging.info('Frame file created...')


logging.info('Sanity check...')
frdata_inj, gps_start_inj, xoffset_inj, xspacing_inj, xunit_inj, yunit_inj = Fr.frgetvect(o_framefile, channel)
if (gps_start_inj!=gps_start  or  xspacing_inj[0]!=xspacing[0]):
    print 'Time staps of input and output frame files are not matching'
    sys.exit(1)

t_inj = np.linspace(inj_start_time, inj_start_time+xspacing[0]*len(inj_data), len(inj_data))
plt.figure(figsize=(7,7))
plt.subplot(2,1,1)
plt.plot(t[start_idx:start_idx+len(inj_data)],frdata_inj[start_idx:start_idx+len(inj_data)]-frdata[start_idx:start_idx+len(inj_data)],label=' Residual between the two time series w and w/o an injection')
plt.plot(t_inj,inj_data,alpha=0.4,label='injected signal')
plt.grid()
plt.legend()
plt.xlabel('time')
plt.ylabel('strain')
plt.xlim([t_inj[0]-2,t_inj[-1]+2])

plt.subplot(2,1,2)
plt.plot(t,frdata_inj-frdata-inj_data_padded  ,label=' Residual after substacting signal')
plt.grid()
plt.legend()
plt.xlabel('time')
plt.ylabel('strain')
plt.xlim([t_inj[0]-2,t_inj[-1]+2])
plt.ylim([-5e-22,5e-22])
plt.savefig('sanity_plots/'+o_framefile.replace('.gwf','.png'))

