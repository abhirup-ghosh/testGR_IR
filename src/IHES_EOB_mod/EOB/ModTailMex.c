/* 
   Modulus of the tail term in the EOB waveform
   MEX-file to compile in Matlab
*/


#include "mex.h"
#include <stdlib.h>
#include "math.h"

#define Pi 3.1415926535897932384626433832795028
#define SQ(x) ((x)*(x))
#define PR 0

/* factorial */
double f14[] = {1,         1,          2,     
		6,         24,         120,
		720,       5040,       40320, 
		362880,    3628800,    39916800,   
		479001600, 6227020800, 87178291200};

double fact(double n)
{

  if(n < 0.) {
    mexErrMsgTxt(" can not compute a negative factorial.");
    return 0;
  }
  
  if (n<=14.) {
    return(f14[(int)n]);
  } else {
    n = n*fact(n-1);
    return(n);
  }

}


/* main routine */
void ModTail(int kmax, int nx, 
	     double *L, double  *M, double *w, 
	     double *mTlm) 
{
  
  int i,j,k, l;
  double f,hatk,x2, y,prod;  

  double j2[kmax+1]; /* 0^2 1^2 ... kmax^2 */
  double oofact2[kmax];
  
  if (PR)
    printf(" kmax = %d nx = %d\n",kmax,nx);


  /* pre-compute some stuffs */
  for( j = 0; j < kmax; j++ ) {
    oofact2[j] = 1./( SQ(fact(L[j])) );    
    j2[j]      = (double)(SQ(j)); 
    /* test fact: Ok */
    /*
      if (PR)       
      printf(" %d -> %g %g\n",j,L[j],fact(L[j])); 
    */
  }
  j2[kmax] = (double)SQ(kmax);
  

  /* main loop */
  for( i = 0; i < nx; i++ ) {      
    
    if (PR) 
      printf(" i = %d (%d) \n",i,nx); 

    for( k = 0; k < kmax; k++ ) {    
      
      f = oofact2[k];
      l = (int)(L[k]);

      if (PR) 
	printf(" k = %d (%d) l = %d\n",k,kmax, l); 
      
      hatk = M[k]*w[i];      
      
      x2   = 4.*(SQ(hatk));			  
      prod = 1.;
      for( j = 1; j <= l; j++ ) {       
	prod *= (j2[j]+x2); 
	/*
	  if (PR)
	  printf(" %d  -> %g %g -> %g\n",j,j2[j],x2, prod);   
	*/
      }

      /*
	if (PR)
	printf(" %d (%g %g) %d -> %g -> %g\n",k,M[k],L[k],i, hatk, prod);  
      */

      y = 4.*Pi*hatk;
      y = y/( 1. - exp(-y) );
      
      /* put the three pieces together */
      mTlm[(i*kmax) + k] = sqrt( f * y * prod ); 

    }

  }
  
  return;
}


/* GATEWAY ROUTINE */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {

  double *L, *M, *w;  
  double *mTlm;
  int kmax, nx;
  
  /* Check for the proper number of arguments */
  if (nrhs != 3) {
    mexErrMsgTxt("Three input required.");
  } else if (nlhs != 1) {
    mexErrMsgTxt("One output required.");
  }

  /* Check for real double vectors */
  if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) )
    mexErrMsgTxt("1st input argument must be a noncomplex double vector.");
  if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) )
    mexErrMsgTxt("2nd input argument must be a noncomplex double vector.");
  if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
    mexErrMsgTxt("3rd input argument must be a noncomplex double vector.");

  /* Get the number of elements */
  kmax = mxGetNumberOfElements(prhs[0]);
  nx   = mxGetNumberOfElements(prhs[2]);

  /* Create a pointer to the input vectors */
  L = mxGetPr(prhs[0]);
  M = mxGetPr(prhs[1]);
  w = mxGetPr(prhs[2]);

  /* Create the output vectors */
  plhs[0] = mxCreateDoubleMatrix(kmax,nx,mxREAL);
  
  /* Assign a pointer to the output */
  mTlm = mxGetPr(plhs[0]);
  
  /* Call the C subroutine */
  ModTail(kmax,nx, L,M, w, mTlm);
    
}


