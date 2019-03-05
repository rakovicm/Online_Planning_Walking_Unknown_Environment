#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "matops.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{    
   
    /// ulazne promenljive
    double  *aN,  *aMs, *aDel, *aBet;
  
    /// izlazne promenljive
    double *J,*JA;
    double mass;
   
    /// privremene promenljive
    int rows,j,k;
	 mwSize   s3N[3]={3,10};
    mwSize   s31[2]={3,1};
   
	 if(nrhs!=1)
      mexErrMsgTxt("One input required\n [J A]=km_JCM(a)\na -flier object \n calculates jacobian associated with Center of Mass of the Flier");
    else if(nlhs >2)
      mexErrMsgTxt("Too many output arguments.");
    else if(!mxIsClass(prhs[0],"k_flier"))
      mexErrMsgTxt("First input must be flier object.");       
	

	aN   =mxGetPr(mxGetField(prhs[0],0,"N"));
    aBet =mxGetPr(mxGetField(prhs[0],0,"bet"));
    aMs  =mxGetPr(mxGetField(prhs[0],0,"ms"));
   aDel =mxGetPr(mxGetField(prhs[0],0,"del"));

	
   
	 rows=(int)aN[0];  
    s3N[1]=rows;
	/// izlazne promenljive
   
	plhs[0] = mxCreateNumericArray(2, s3N,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericArray(2, s31,mxDOUBLE_CLASS,mxREAL);    
    
   
    J  = mxGetPr(plhs[0]);
    JA = mxGetPr(plhs[1]);    
    /// nadjemo prethodnika
	mass=0;
   for (k=0;k<rows;k++)
   {
	   for (j=0; j<=k;j++)
			vecaddcoef(J+3*j,aMs[k], aBet+3*k+3*rows*j, J+3*j);
	   vecaddcoef(JA, aMs[k], aDel, JA);
	   mass+=aMs[k];
		aDel+=3;
   }
 
   vecmulscalar(JA,1/mass, JA);
   matmulscalar(J, 1/mass, 3, rows, J);
}
   