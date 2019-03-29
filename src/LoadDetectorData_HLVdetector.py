
import numpy as np
import scipy
from RotateVector_HLVdetector import *
import matplotlib.pylab as plt

class detectorclass():
	h=[]
	d=np.array([])
	lambda_ = np.array([])
	phi = np.array([])
	V=np.array([])
	name = []
	n=np.array([])
	psi = []
	X=np.array([])
	Y = np.array([])
	Z=np.array([])

def LoadDetectorData(detector):

    # Local Variables: fs, Omega, GEO_Z, thetaCap, GEO_Y_V, theta, detector, GEO_X, GEO, GEO_X_V, phi, psi, localNorth, LHO, Det, LLO, phiCap, V, Y, X, Z, GEO_Y, d, name, Vi, n, h, lambda
    # Function calls: upper, cos, CartesianPointingVector, error, orderfields, cross, sin, not, recognised, LoadDetectorData, sum, RotateVector, pi, strcmp, norm
    #% LoadDetectorData: Return a struct containing data on the location and 
    #% orientation and orientation of the specified gravitational-wave detector 
    #% or detector site.
    #%
    #%   Det = LoadDetectorData(detector)
    #%
    #%   detector  String. Detector site of interest.  Available detectors are 
    #%             'AIGO' (or 'A1', or 'A'), 'ALLEGRO' (or 'AL'), 'AURIGA' (or
    #%             'AU'), 'EXPLORER' (or 'EX'), 'GEO' (or 'G1', 'G'), 'LCGT',
    #%             'INDIGO' (or I1), 'LHO' (or 'H1', 'H2', 'H'), 'LLO' (or 
    #%             'L1', 'L'), 'NAUTILUS' (or 'NA'), 'NIOBE' (or 'NI'), 'TAMA' 
    #%             (or 'T1', 'T'), and 'VIRGO' (or 'V1', 'V').
    #%
    #%   Det       Struct containing data on the detector site.  The fields are:
    #%
    #%               Det.d       Antenna response matrix.  For IFOs this is
    #%                             Det.d = 1/2*(Det.X*Det.X'-Det.Y*Det.Y');
    #%                           For bars this is
    #%                             Det.d = Det.n*Det.n'-[1 0 0; 0 1 0; 0 0 1]/3;
    #%               Det.h       Detector elevation (m).
    #%               Det.lambda_  Detector longitude (deg E of N).
    #%               Det.n       Bar pointing direction in Cartesian Earth-based
    #%                           coordinates ([] for IFOs).
    #%               Det.name    Name of detector (full, or abbreviated form:
    #%                           one letter for IFOs, two letters for bars).
    #%               Det.phi     Detector latitude (deg N).
    #%               Det.psi     Detector azimuth / bar pointing direction (deg 
    #%                           W of N, x-arm direction or [] for IFOs).
    #%               Det.V       Location of vertex/beamsplitter (IFOs) or bar 
    #%                           in Cartesian Earth-based coordinates (m).
    #%               Det.X       X arm pointing direction ([] for bars).
    #%               Det.Y       Y arm pointing direction ([] for bars).
    #%               Det.Z       Right-handed direction orthogonal to the 
    #%                           detector X-Y plane (IFOs), or zenith sky
    #%                           direction at the detector (bars). 
    #%
    #% Unknown values are reported as [].
    #%
    #% For bar detectors, AIGO, and LCGT we use a spherical Earth model. Data
    #% for all other interferometers follow the WGS-84 model. The antenna
    #% response factors from the spherical Earth model are only accurate to the
    #% 1% level.
    #% $Id: LoadDetectorData.m 314 2012-04-30 23:36:01Z ajith.parameswaran@LIGO.ORG $
    #% Notes:
    #% ------
    #% For a bar detector we have an alternative (equivalent) way to compute F+,
    #% Fx:
    #%
    #%   h(t) = F_+ h_+ + F_x h_x 
    #%
    #% where
    #%
    #%  F_+ = sin^2(theta) * cos(2*phi)
    #%  F_x = sin^2(theta) * sin(2*phi)
    #%
    #%  cos(theta) = \vec{n}\cdot\vec{k}
    #%
    #% where phi is the angle between the polarization vector and the projection
    #% of the long axis of the detector into the transverse plane of the wave, 
    #% and where \vec{n} is the pointing direction of the symmetry axis of the
    #% cylinder and \vec{k} is the direction to to source.
    #%----- Guarantee that the name of the string is all upper case.

   
    detector = detector.upper()
    if (detector== 'LHO') or (detector== 'H1') or (detector =='H2') or (detector== 'H'):
	Det =  detectorclass()      
        #%----- LHO
        #%
        #% Values from Althouse et al, LIGO-P000006-D-E
        #% These give d-matrix values different from those from Allen 
        #% in second or third decimal place.
	LHO =  detectorclass()   
   
        LHO.X = np.vstack((-0.223891216,0.799830697,0.556905359))
	LHO.Y = np.vstack((-0.913978490,0.0260953206,-0.404922650)) 
    	LHO.Z = np.vstack((-0.338402190,-0.599658144,0.725185541)) 

        #%
        LHO.d = np.dot(0.5, np.dot(LHO.X, LHO.X.conj().T)-np.dot(LHO.Y, LHO.Y.conj().T))
        #%----- IFO vertex ("V") position
        #%----------- Cartesian Earth-based coordinates (m)
        LHO.V = np.array(np.vstack((np.hstack([-2.161414928e6]), np.hstack([-3.834695183e6]), np.hstack([4.600350224e6]))))
        #%----------- WGS-84 geodetic coordinates (m/deg/deg)
        #%            Elevation h, longitude lambda, latitude phi
        LHO.h = 142.555
        LHO.lambda_ = -np.array([119.+24./60.+27.565681/3600.])
        LHO.phi = 46.+27./60.+18.527841/3600.
        LHO.fs = 16384.
        #% Hz
        LHO.name = 'H'
        LHO.n = np.array([])
        LHO.psi = np.array([])
        Det = LHO
        
    elif (detector== 'LLO') or (detector== 'L1') or (detector== 'L'):
	LLO =  detectorclass() 
        #%----- LLO 
        #%
        #% Values from Althouse et al, Review of Scientific Instruments 72, 3086 2001,
        #% LIGO-P000006-D-E.
        #% These give d-matrix values different from those from Allen 
        #% in third decimal place.
        LLO.X = np.vstack((-0.954574615,-0.141579994,-0.262187738))
        LLO.Y = np.vstack((0.297740169,-0.487910627,-0.820544948))
        LLO.Z = np.vstack((-0.011751435,-0.861335199,0.507901150))
        #%
        LLO.d = 0.5*( np.dot(LLO.X, LLO.X.conj().T)-np.dot(LLO.Y, LLO.Y.conj().T))
        #%----- IFO vertex ("V") position
        #%----------- Cartesian Earth-based coordinates (m)
        LLO.V = np.vstack((-7.427604192e4,-5.496283721e6,3.224257016e6))
        #%----------- WGS-84 geodetic coordinates (m/deg/deg)
        #%            Elevation h, longitude lambda, latitude phi
        LLO.h = -6.574
        LLO.lambda_ = -np.array(np.hstack([90.+46./60.+27.265294/3600.]))
        LLO.phi = 30.+33./60.+46.419531/3600.
        LLO.fs = 16384.
        #% Hz
        LLO.name = 'L'
        LLO.n = np.array([])
        LLO.psi = np.array([])
        Det = LLO
        
    elif (detector== 'VIRGO') or (detector== 'V1') or (detector== 'V'):
	Det =  detectorclass()      
        #%----- Virgo
        #%
        #% Data are taken from Anderson, Brady, Creighton, 
        #% and Flanagan, PRD 63 042003, 2001.
        #% Since the data are only given to 4 significant digits, 
        #% remormalize the units vectors so that d_ab will be traceless.
        Det.X = np.vstack([-0.70050, 0.20850,0.68260])
	Det.X = Det.X/((np.linalg.norm(Det.X)))
	#print Det.X
	Det.Y = np.vstack((-0.05380,-0.96910,0.24080))
        Det.Y = Det.Y/((np.linalg.norm(Det.Y)))
	#print Det.Y
	Det.Z = np.cross(np.hstack(Det.X),np.hstack(Det.Y))
	
	Det.d = np.dot(0.5, np.dot(Det.X, Det.X.conj().T)-np.dot(Det.Y, Det.Y.conj().T))
	#print Det.d
        #%
        #%----- IFO vertex ("V") position
        #%----------- Cartesian Earth-based coordinates (m)
        Det.V = 1e6* np.vstack((4.546374,0.842990,4.378577))
        #%----------- WGS-84 geodetic coordinates (m/deg/deg)
        #% Elevation h, longitude lambda, latitude phi
        Det.h = 51.8840
        Det.lambda_ = np.array([10.+30./60.+16.1878/3600.])
        Det.phi = np.array([43.+37./60.+53.0921/3600.])
        #% psi_x =  70.5674  % N of E  arm pointing directions
        #% psi_y = 160.5674  % N of E
        Det.name = 'V'
        Det.n = np.array([])
        Det.psi = np.array([])
        
    else:
        print('Not LIGO or VIRGO Detector')
        
    
    #% ---- Sort fields.
    #Det = orderfields(Det)
    #%----- Done
    #return []
    return [Det]

#et = LoadDetectorData('A')
#	print et
 
