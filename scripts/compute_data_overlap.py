import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt


h1_times = np.loadtxt('../data/LIGO_real_data_analysis/gpstimes_goodsegments_H1.dat')
l1_times = np.loadtxt('../data/LIGO_real_data_analysis/gpstimes_goodsegments_L1.dat')
h1l1_times = np.loadtxt('../data/LIGO_real_data_analysis/gpstimes_commonsegments_H1L1.dat') 

gps_start_time_global = min(h1_times[0][0], l1_times[0][0]) 
gps_end_time_global = max(h1_times[-1][1], l1_times[-1][1]) 

#gpstimes_list = np.arange(gps_start_time, gps_end_time, 1)
#f = open('../data/LIGO_real_data_analysis/gpstimes_commonsegments_H1L1.dat', 'w')

plt.figure(figsize=(16,12))
#for gpstime in gpstimes_list:
#  if any(lower <= gpstime <= upper for (lower, upper) in h1_times) and any(lower <= gpstime <= upper for (lower, upper) in l1_times):
#    f.write('%d\n'%gpstime)

#f.close()

plt.subplot(411)
for idx in range(len(l1_times)):                                               
     plt.axvspan(l1_times[idx][0], l1_times[idx][1], color='g', alpha=0.5)            
     plt.hold(True)
plt.xlim([gps_start_time_global, gps_end_time_global])

print '... plotted L1 SQ data'
plt.subplot(412)
for idx in range(len(h1_times)):                                               
     plt.axvspan(h1_times[idx][0], h1_times[idx][1], color='r', alpha=0.5)            
     plt.hold(True)
plt.xlim([gps_start_time_global, gps_end_time_global])

print '... plotted H1 SQ data'

f = open('../data/LIGO_real_data_analysis/gpsstartendtimes_commonsegments_H1L1.dat', 'w')
for gps_time in h1l1_times:
     if (gps_time - 1) not in h1l1_times:
        f.write('%d\t'%gps_time)
	print gps_time
     if (gps_time + 1) not in h1l1_times:           
        f.write('%d\n'%gps_time)
	print gps_time
f.close()

h1l1_times = np.genfromtxt('../data/LIGO_real_data_analysis/gpsstartendtimes_commonsegments_H1L1.dat', skip_footer=1)
plt.subplot(413)
for idx in range(len(h1l1_times)):                                               
     plt.axvspan(h1l1_times[idx][0], h1l1_times[idx][1], color='b', alpha=0.5)            
     plt.hold(True)
plt.xlim([gps_start_time_global, gps_end_time_global])

plt.subplot(414)
for idx in range(len(h1l1_times)):
     plt.axvspan(h1l1_times[idx][0], h1l1_times[idx][1], color='b', alpha=0.5)
     plt.hold(True)

print '... plotted H1L1 SQ data'

plt.xlabel('gpstimes')

#leg = plt.legend(['L1 SQ data', 'H1 SQ data', 'H1L1 SQ data'], frameon=False) 
#leg.legendHandles[0].set_color('g')
#leg.legendHandles[1].set_color('r')
#leg.legendHandles[2].set_color('b')
#leg.legendHandles[2].set_linestyle('-')

plt.tight_layout()
plt.grid()
plt.savefig('../data/LIGO_real_data_analysis/SQdata_H1L1.png', dpi=300)
plt.savefig('../data/LIGO_real_data_analysis/SQdata_H1L1.pdf')
plt.close()
