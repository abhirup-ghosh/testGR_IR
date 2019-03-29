/* 
   Lagrangian interpolation of 1d data
   MEX-file to compile in Matlab
   Original C routines by P.Galaviz 
*/


#include "mex.h"


int locate(double x[],unsigned int Nx,double xval);
double LagrangePoly(const double xval, double *x,int i,int ni,int nf);
double Lagrange(const double xval,double *f,double *x,const int Nx,const int order); 


/* GATEWAY ROUTINE */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {

  double *ff, *xx, *xv;  
  double *fv;
  int mrows, ncols, k, order, sizex, sizef, size;
  
  /* Check for the proper number of arguments */
  if (nrhs != 4) {
    mexErrMsgTxt("Four input required.");
  } else if (nlhs != 1) {
    mexErrMsgTxt("One output required.");
  }

  /* 1st input must be a noncomplex scalar double */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(mrows == 1 && ncols == 1)) {
    mexErrMsgTxt("Bad interpolation order");
  }
    
  /* Other input args must be real doubles */
  if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) )
    mexErrMsgTxt("1st input argument must be a noncomplex double vector.");
  if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
    mexErrMsgTxt("2nd input argument must be a noncomplex double vector.");
  if ( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) )
    mexErrMsgTxt("3rd input argument must be a noncomplex double vector.");
 
  /* Get the number of elements */
  sizef = mxGetNumberOfElements(prhs[1]);
  sizex = mxGetNumberOfElements(prhs[2]);
  if ( sizex != sizef )
    mexErrMsgTxt(" f and x must be of the same size ");
 
  size = mxGetNumberOfElements(prhs[3]);

  /* Create a pointer to the input vectors */
  order = (int)mxGetScalar(prhs[0]); 
  xx = mxGetPr(prhs[1]);
  ff = mxGetPr(prhs[2]);
  xv = mxGetPr(prhs[3]);

  /* Create the output vectors */
  plhs[0] = mxCreateDoubleMatrix(1,size,mxREAL);
  
  /* Assign a pointer to the output */
  fv = mxGetPr(plhs[0]);
  
  /* Call the C subroutine */
  for( k = 0; k < size; k++ ) 
    fv[k] = Lagrange( xv[k], ff, xx, sizex, order );   
  
  
}


/* Locate NN */
/* NR routine for 1-offset array, invoke it with x-1 */
int locate( double x[], unsigned int Nx, double xval )
{
  
  int ju,jm,jl;
  int j;  
  int ascnd;
  
  jl=0;
  ju=Nx+1;

  ascnd=(x[Nx] >= x[1]);  
  while (ju-jl > 1) {
    jm=(ju+jl) >> 1;
    if (xval >= x[jm] == ascnd)
      jl=jm;
    else
      ju=jm;
  }
  
  if (xval == x[1]) j=1;
  else if(xval == x[Nx]) j=Nx-1;
  else j=jl;
  return j-1; 

}


/* Routines for 1D Lagrangian interpolation */
double LagrangePoly( const double xval, 
		     double *x, 
		     int i, 
		     int ni,  int nf )
{
  
  if( ni >= 0 && i >= ni && i <= nf ) {

    double lp = 1.0;
    const double xi = x[i];
    int j;
    for( j = ni; j < i; j++ )
      lp *= ( xval - x[j] ) / (xi - x[j] );
    for( j = i+1; j <= nf; j++ )
      lp *= ( xval - x[j] ) / (xi - x[j] );
    
    return lp;
    
  } else return 0.0;
  
}

double Lagrange( const double xval,       /* Value to interpolate */
		 double *f, double *x,    /* Tabulated values */
		 const  int Nx,           /* Size of x[0...Nx-1], f[0...Nx-1] */
		 const int order          /* Order of the interpolation */
		 )
{
  
  double fval = 0.0;
  
  /* Locate the index of the nearest point to xval */
  int i = locate( x-1, Nx, xval );
  
  /* Find the interval */
  int desp_i, desp_f;
  
  if( order%2 == 0 ) {
    desp_i = order / 2;
    desp_f = desp_i - 1;
  } else {
    desp_i = (order - 1) / 2;
    desp_f = desp_i;
  }

  int ni_i = i-desp_i;
  int ni_f = i+desp_f;
  int ord = ni_f - ni_i;
  
  if( ni_i < 0 ) {
    ni_i = 0;
    ni_f = ord;
  }
  
  if( ni_f >= Nx ) {
    ni_f = Nx-1;
    ni_i = ni_f - ord;
  }
  
  /* Here interpolation */
  int ii;
  for ( ii = ni_i; ii <= ni_f; ii++ )
    fval += f[ii] * LagrangePoly( xval, x, ii, ni_i, ni_f ); 
  return fval;
  
}
