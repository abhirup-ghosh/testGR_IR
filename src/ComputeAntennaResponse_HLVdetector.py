
import numpy as np
import scipy
from LoadDetectorData_HLVdetector import *
import matplotlib.pylab as plt


def ComputeAntennaResponse(phi, theta, psi, detector):

    # Local Variables: Fp, phi, psi, d, DetData, Fc, m1, maxLength, theta, detector
    # Function calls: m2, cos, isscalar, m3, ComputeAntennaResponse, max, ischar, isequal, length, ones, LoadDetectorData, error, isnumeric, sin, size
    #% ComputeAntennaResponse - Compute the antenna response factors for a
    #% specfic detector and set of sky positions and polarizations.  
    #%
    #%   [Fp, Fc] = ComputeAntennaResponse(phi,theta,psi,detector)
    #%
    #%  phi      Vector.  Azimuthal angle of the sky position of the source in
    #%           Earth-centered coordinates.  The prime meridian at Greenwich is
    #%           phi=0.
    #%  theta    Vector.  Polar angle of the sky position of the source in
    #%           Earth-centered coordinates.  The north pole is at theta=0; the
    #%           south pole is at theta=pi.
    #%  psi      Vector.  The polarization angle of the gravitational wave,
    #%           defined as the angle counterclockwise about the direction of
    #%           PROPAGATION from a line of constant theta pointing to
    #%           decreasing phi to the positive x axis of the source coordinates
    #%           (the source "+" polarization direction).
    #%  detector String or 3x3 numerical array.  If a string, it specifies the 
    #%           detector site on to which to project the signal, and must be 
    #%           one of the values recognized by the function LoadDetectorData.
    #%           If an array, then it specifies the detector response matrix 
    #%           d^{ab} to be used.
    #%
    #%   Fp      Column vector.  The "plus" antenna response factors for the  
    #%           specified detector.
    #%   Fc      Column vector.  The "cross" antenna response factors for the  
    #%           specified detector.
    #%
    #% The vectors phi, theta, psi must have the same length, except that one or
    #% more may be scalar, in which case the scalars are expanded to vectors of
    #% the same size.
    #%
    #% initial write: Patrick J. Sutton 2004.04.26
    #%
    #% $Id: ComputeAntennaResponse.m 314 2012-04-30 23:36:01Z ajith.parameswaran@LIGO.ORG $
    #%----- Make sure input arguments are column vectors.

	#%----- Make sure input arguments are column vectors.
	if (len(phi.shape)>1):
    	    phi = np.vstack(np.concatenate(phi.T)) #concatenation converts matrix to single array by stiching all rows together
	if (len(theta.shape)>1):
	    theta = np.vstack(np.concatenate(theta.T))
	if (len(psi.shape)>1):
	    psi = np.vstack(np.concatenate(psi.T))

#%----- Verify that input angle vectors are the same size, or are scalars.
	maxLength = max([len(theta), len(phi), len(psi)])
	if (len(phi)!=maxLength):
#    % ---- Make sure it's a scalar, then expand to vector.
 	    if (np.isscalar(phi)):
	        phi = phi*np.ones((maxLength,1))
	    else:
        	print('Sky position, polarization angles phi, theta, psi must have the same length.')

	if (len(theta)!=maxLength):
#    % ---- Make sure it's a scalar, then expand to vector.
 	    if (np.isscalar(theta)):
	        theta = theta*np.ones((maxLength,1))
	    else:
	        print('Sky position, polarization angles phi, theta, psi must have the same length.')
   
	if (len(psi)!=maxLength):
#    % ---- Make sure it's a scalar, then expand to vector.
            if (np.isscalar(psi)):
	        psi = psi*np.ones((maxLength,1))
	    else:
	        print('Sky position, polarization angles phi, theta, psi must have the same length.')


#%----- Get the detector response matrix "d^{ab}_i".If input argument is a detector/site name, retrieve this data by
#%      calling LoadDetectorData.  If the input argument is a 3x3 numerical array, use that for the detector response matrix. 
	if (isinstance(detector,str)):
#   %----- Load data on detector.
	    DetData = LoadDetectorData(detector)
	    d = DetData[0].d
	
	 
	elif (isinstance(detector,np.ndarray) and (np.shape(detector)==(3,3))):
    	    d = detector	
	elif isinstance(detector,list) and (np.shape(detector)==(3,3)):
    	    d = detector
	else:
	    print('Detector not recognized. 4th argument should be a detector/site name or a 3x3 array.')
#%----- Convert to vector.
	d = np.vstack(np.concatenate(d))#d=transpose(d) as d is symmetric matrix
	
#%----- Compute polarization tensors (functions of the sky position and polarization). 
	m1 =  np.sin(phi)*np.cos(psi)-np.cos(phi)*np.cos(theta)*np.sin(psi)
	m2 = -np.cos(phi)*np.cos(psi)-np.sin(phi)*np.cos(theta)*np.sin(psi)
	m3 = np.sin(theta)*np.sin(psi)
	n1 = -np.sin(phi)*np.sin(psi)-np.cos(phi)*np.cos(theta)*np.cos(psi)
	n2 =  np.cos(phi)*np.sin(psi)-np.sin(phi)*np.cos(theta)*np.cos(psi)
	n3 =  np.sin(theta)*np.cos(psi)
	mm = np.array([m1*m1, m1*m2, m1*m3, m2*m1, m2*m2, m2*m3, m3*m1, m3*m2, m3*m3])
	mn = np.array([m1*n1, m1*n2, m1*n3, m2*n1, m2*n2, m2*n3, m3*n1, m3*n2, m3*n3])
	nm = np.array([n1*m1, n1*m2, n1*m3, n2*m1, n2*m2, n2*m3, n3*m1, n3*m2, n3*m3])
	nn = np.array([n1*n1, n1*n2, n1*n3, n2*n1, n2*n2, n2*n3, n3*n1, n3*n2, n3*n3])
	e_plus = mm - nn
	e_cross = mn + nm
	e_plus = np.transpose(e_plus)
	e_cross = np.transpose(e_cross)
	
	
#%----- Compute waveform projected onto antenna pattern.
	Fp=np.zeros(len(psi))
	Fc=np.zeros(len(psi))
	Fp = np.vstack(np.dot(e_plus,d))	
	Fc = np.vstack(np.dot(e_cross,d))

	return Fp,Fc
#% % ---- This direct matrix form is clearer, but only works for a single
#% %      sky position and polarization angle at a time.
#% m = [  sin(phi)*cos(psi)-cos(phi)*cos(theta)*sin(psi) ; ...
#%       -cos(phi)*cos(psi)-sin(phi)*cos(theta)*sin(psi) ; ...
#%        sin(theta)*sin(psi) ];
#% n = [ -sin(phi)*sin(psi)-cos(phi)*cos(theta)*cos(psi) ; ...
#%        cos(phi)*sin(psi)-sin(phi)*cos(theta)*cos(psi) ; ...
#%        sin(theta)*cos(psi) ];
#% e_plus = m*m' - n*n';
#% e_cross = m*n' + n*m';
#% Fp = trace(e_plus*d);
#% Fc = trace(e_cross*d);

#%----- Done


