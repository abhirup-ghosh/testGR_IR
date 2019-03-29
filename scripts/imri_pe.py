#import matplotlib as mpl
#mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import os, lal, commands
import emcee
import optparse as op
import corner
import nr_fits as nr
import imrtgrutils_final as tgr
from lal import MTSUN_SI as LAL_MTSUN_SI
import bayesian as ba
import scipy
import time
import noiseutils as nu
import scipy
from scipy.interpolate import interp1d
from pylal import antenna
import td_likelihood_utils as uc
import standard_gwtransf as gw

if __name__ == '__main__':

        start_time = time.time()

	approx = 'IMRPhenomPv2'
	srate = 2048
	flow = 10.
	dL = 100
	m1_inj, m2_inj = 20., 200.
	gpsTime, rightAscension, declination, inclination, polarization, unit, detector =  1126285216, 0., 0., 0, 0, 'radians', 'H1'

	f, hf = uc.wavegen_from_lalsiminspiralFD(approx, m1_inj, m2_inj, dL, flow, srate, gpsTime, rightAscension, declination, inclination, polarization, unit, detector, asd_file=None)

	freqs, psd_adligo = uc.adLIGO_psd(10., 2048.)
	psd_adligo_interp_object = interp1d(freqs, psd_adligo, fill_value=0., bounds_error=False)
	psd_adligo_interp = psd_adligo_interp_object(f)
	asd_adligo = np.sqrt(psd_adligo_interp)

	df = hf + asd_adligo
	

	plt.figure()
	plt.plot(f, hf.real)
	plt.plot(f, asd_adligo)
	plt.plot(f, df.real)
	plt.show()
