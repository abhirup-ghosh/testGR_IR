
import numpy as np
import scipy
#import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def RotateVector(X, Y, psi):

    # Local Variables: psi, S, R, Y, X, Z
    # Function calls: cos, sum, error, isequal, RotateVector, sin, size
    #%
    #% Z = RotateVector(X,Y,psi)
    #%
    #% Rotate vector X around the vector Y by angle psi.  All vectors 
    #% are assumed to be 3-dimensional and in Cartesian coordinates; 
    #% Y must be non-null.  The resulting vector Z is returned as a 
    #% column vector.
    #% -- Patrick J. Sutton 2004/04/23
    #%----- Check that X,Y are 3-vectors and convert X, Y to column vectors
    #%      if necessary.
    if (len(np.hstack(X).shape)==3.):
	X=np.vstack(X)
	print('X is', X)
    else:
	print('Error encountered - X is not a (1,3) vector')
	print X.shape
    
    if (len(np.hstack(Y)) == 3.):
	Y=np.vstack(Y)
	print('Y is',Y)
    else:
	print('Error encountered because Y is {} vector'.format(Y.shape))

    if ((np.array(X.shape) == np.array((3., 1.))).all() and (np.array(Y.shape)== np.array((3., 1.))).all()):
    	print('Input vectors must be 3-vectors.')
    
    
    #%----- Check that Y is non-null.
    if np.linalg.norm(Y)== 0.:
        print('Input vector Y (axis of rotation) must be non-null.')
    
    
    #%----- This script takes the vector X (Cartesian coordinates) and rotates it 
    #%      by an angle psi about the vector Y.  It returns the rotated vector Z. 
    #%      Note: the axis must be specified by a unit vector, therefore normalize Y.
    Y = Y/((np.linalg.norm(Y))**(0.5))
    #%----- Rotate.  (Checked matrix is correct for general angle of rotation 
    #%      psi of general initial vector X about general unit vector Y).
    S = np.array(np.vstack((np.hstack([0., -Y[2], Y[1]]), np.hstack([Y[2], 0., -Y[0]]), np.hstack([-Y[1], Y[0], 0.]))))
    R = np.array(np.vstack((np.hstack((1., 0., 0.)), np.hstack((0., 1., 0.)), np.hstack((0., 0., 1.)))))+np.dot(np.sin(psi), S)+np.dot(1.-np.cos(psi), S**2.)
    Z = np.dot(R, X)
    #%----- Done.
    return []
    return [Z]
