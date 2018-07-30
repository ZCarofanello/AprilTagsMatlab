/*
 * HammingDistMex.c - Calculates the hamming distance using x86 instructions
 *
 * Calculates the hamming distance between two numbers
 *
 * The calling syntax is:
 *
 *		out = HammingDistMex(NumA, NumB)
 *
 * This is a MEX file for MATLAB.
*/

#include "mex.h"
#include "matrix.h"
#include <nmmintrin.h>
#include <stdint.h>

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    
//unsigned __int64 NumA;
double NumA, NumB;
unsigned __int64 Test;
int PopCount = 0;
double* NumC;
    
/* variable declarations here */
if(nrhs != 2) {
	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Needs two inputs");
}

if(nlhs != 1) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "One output required.");
}

/* make sure the first input argument is scalar */
if( !mxIsUint64(prhs[0]) || mxIsComplex(prhs[0]) ||
    mxGetNumberOfElements(prhs[0]) != 1) {
		
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar",
                      "First input must be scalar.");
}
  plhs[0] = mxCreateDoubleMatrix(1,1,0);
  NumA = mxGetScalar(prhs[0]);
  NumB = mxGetScalar(prhs[1]);
  
  NumA = (int)NumA ^ (int)NumB;
  
  NumC = mxGetPr(plhs[0]);
  PopCount = _mm_popcnt_u64((unsigned __int64)NumA);
  
  *NumC = (double)PopCount;

return;

}