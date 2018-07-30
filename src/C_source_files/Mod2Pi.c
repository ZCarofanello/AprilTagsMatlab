/*
 * Mod2Pi.c - wraps the value of NumA between -Pi - Pi 
 *
 * The calling syntax is:
 *
 *		out = Mod2Pi(NumA)
 *
 * This is a MEX file for MATLAB.
*/

#include "mex.h"
#include "matrix.h"
#include <nmmintrin.h>
#include <stdlib.h>

#define M_PI   3.14159265358979323846264338327950288


	// Returns a result in [-Pi, Pi]
	double mod2pi(double vin) {
		const double twopi = 2 * (float)M_PI;
		const double twopi_inv = 1.f / (2.f * (double)M_PI);
		double absv = abs(vin);
		double q = absv*twopi_inv + 0.5f;
		int qi = (int) q;
		double r = absv - qi*twopi;
		return (vin<0) ? -r : r;
	}


/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    
//unsigned __int64 NumA;
double NumA;
double* NumB;
    
/* variable declarations here */
if(nrhs != 1) {
	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Only one input");
}

if(nlhs != 1) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "One output required.");
}

/* make sure the first input argument is scalar */
if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
    mxGetNumberOfElements(prhs[0]) != 1) {
		
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar",
                      "First input must be scalar.");
}
  plhs[0] = mxCreateDoubleMatrix(1,1,0);
  NumA = mxGetScalar(prhs[0]);
  NumB = mxGetPr(plhs[0]);
  
  
  
  *NumB = mod2pi(NumA);

return;

}