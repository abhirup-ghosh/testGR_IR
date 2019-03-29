import numpy as np
import os
import matplotlib.pyplot as plt


#inj_list = np.loadtxt('../frames/IHES_20161130/injection_list.txt')
inj_list = np.linspace(1, 5000, 5000)
inj_list = inj_list.astype(int)
seglen = 8

data = np.genfromtxt('/home/abhirup/Documents/Work/testGR_IR/injections/popsynth_injections_SEOBNRv2_ROM_DoubleSpinthreePointFivePN.txt', dtype=None, names=True)
mass1, mass2, spin1x, spin2x, spin1y, spin2y, spin1z, spin2z, longitude, latitude, inclination, polarization, distance, geocent_end_time, geocent_end_time_ns, coa_phase = data['mass1'], data['mass2'], data['spin1x'], data['spin2x'], data['spin1y'], data['spin2y'], data['spin1z'], data['spin2z'], data['longitude'], data['latitude'], data['inclination'], data['polarization'], data['distance'], data['geocent_end_time'], data['geocent_end_time_ns'], data['coa_phase']

q = mass1/mass2

out_dir = '/home/abhirup/Documents/Work/testGR_IR/frames/IHES_20170414_GR_NDRP15'
os.system('mkdir -p %s'%out_dir)

for inj in inj_list:
	if q[inj-1] < 1:
		q[inj-1] = 1./q[inj-1]
	print '%.2f'%q[inj-1], mass1[inj-1], mass2[inj-1], inj

	run_cmd = './inject_NR.py -D /home/abhirup/Documents/Work/testGR_IR/waveforms/GR_injections_NDRP15/injection_%d_q_%.2f_NDRP15.dat --data-format IHES_modGR -M %f --trig-time=%d --dist=%f --incl=%f --ra=%f --dec=%f --psi=%f --ifos H1,L1,V1 --noise-models ZERO_DET_high_P --f-low 5 --f-ref 5 --f-start 10 --srate 4096 --seglen %d -o %s --frame-prefix IHES_%04d > %s/injection_IHES_%04d_modGR.log'%(inj, q[inj-1], mass1[inj-1] + mass2[inj-1], geocent_end_time[inj-1], distance[inj-1], inclination[inj-1], longitude[inj-1], latitude[inj-1], polarization[inj-1], seglen, out_dir, inj, out_dir, inj)
	os.system(run_cmd)
	#print run_cmd
