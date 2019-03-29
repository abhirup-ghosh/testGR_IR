import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pycbc
from pycbc.waveform import get_td_waveform
from pycbc.detector import Detector
import standard_gwtransf as gw

apx = 'IMRPhenomPv2'

mc_1, mc_2 = 12., 28.
q, spin1z, spin2z, inclination, coa_phase, delta_t, f_lower, end_time, declination, right_ascension, polarization = 0.6347, 0.7026, 0.2198, 0.4589, 3.7708, 1./4096, 20., 1186741861.54, -0.7085, 1.0623


for mc in [m1_1, mc_2]:

	m1, m2 = gw.comp_from_mcq(mc, q)

	# NOTE: Inclination runs from 0 to pi, with poles at 0 and pi
	#       coa_phase runs from 0 to 2 pi.
	hp, hc = get_td_waveform(approximant=apx,mass1=m1,mass2=m2,spin1z=spin1z,spin2z=spin2z, inclination=inclination,coa_phase=coa_phase,delta_t=delta_t, f_lower=f_lower)

	det_h1 = Detector('H1')
	det_l1 = Detector('L1')
	det_v1 = Detector('V1')

	# Choose a GPS end time, sky location, and polarization phase for the merger
	# NOTE: Right ascension and polarization phase runs from 0 to 2pi
	#       Declination runs from pi/2. to -pi/2 with the poles at pi/2. and -pi/2.
	hp.start_time += end_time
	hc.start_time += end_time

	signal_h1 = det_h1.project_wave(hp, hc,  right_ascension, declination, polarization)
	signal_l1 = det_l1.project_wave(hp, hc,  right_ascension, declination, polarization)
	signal_v1 = det_v1.project_wave(hp, hc,  right_ascension, declination, polarization)

	plt.plot(signal_h1.sample_times, signal_h1, label='H1')
	plt.plot(signal_l1.sample_times, signal_l1, label='L1')
	plt.plot(signal_v1.sample_times, signal_v1, label='V1')

plt.ylabel('Strain')
plt.xlabel('Time (s)')
plt.legend()
plt.show()
