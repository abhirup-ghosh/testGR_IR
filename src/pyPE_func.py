import numpy as np
import numpy.linalg as la
import lal
import LoadDetectorData_HLVdetector as ld
from lalinspiral.sbank import psds
import mywavegen as mw

def antenna_pattern_fn(det,ra,dec,psi,gmst):
    D = ld.LoadDetectorData(det)[0].d
    Fp,Fc = lal.ComputeDetAMResponse(D,ra,dec,psi,gmst)
    return Fp,Fc

def apply_antenna_pattern(Fp,Fc,h):
    return np.vectorize(complex)(Fp*np.real(h),Fc*np.imag(h))

def log_fourier_likelihood(x,h,psd):
    l = -la.norm(np.abs(x-h)**2./psd)/2.
    return l

def autocorr(x,lag_max=1):
    return np.array([np.sum(np.corrcoef(np.array([x[0:len(x)-t],x[t:len(x)]]))[0,1]) for t in range(0,lag_max)])

def generate_fourier_colored_psd(f,det_name):
    psd = psds.noise_models[det_name](f)
    return psd

def generate_fourier_colored_noise(f,det_name):
    psd = generate_fourier_colored_psd(f,det_name)
    real_noise = [np.random.normal(0.,1.)*np.sqrt(p/2.) for p in psd]
    imag_noise = [np.random.normal(0.,1.)*np.sqrt(p/2.) for p in psd]
    noise = np.vectorize(complex)(real_noise,imag_noise)
    return noise
    
