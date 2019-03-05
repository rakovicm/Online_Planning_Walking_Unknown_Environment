#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "matops.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{    
   
    /// ulazne promenljive
    double  *aN,  *aMs, *aRc, *aV;//, *aCh,*aM,*aBs,*aCl,*lnk,*gsv,*aAlf,*aBet, *aGam, *aDel,*aGe, *aOmg, *aDq, *cnp, *cvec, *aA; 
  
    /// izlazne promenljive
    double *rcm,*vcm;
    double mass;
   
    /// privremene promenljive
    int rows,j,k;
	 mwSize   s3N[3]={3,10};
    mwSize   s31[2]={3,1};
   
	 if(nrhs!=1)
      mexErrMsgTxt("One input required\n [rcm vcm]=km_CentarMase(a)\na -flier object \n position and velocity of Center of Mass of the Flier");
    else if(nlhs >2)
      mexErrMsgTxt("Too many output arguments.");
    else if(!mxIsClass(prhs[0],"k_flier"))
      mexErrMsgTxt("First input must be flier object.");       
	

	aN   =mxGetPr(mxGetField(prhs[0],0,"N"));
    aMs  =mxGetPr(mxGetField(prhs[0],0,"ms"));
    aRc  =mxGetPr(mxGetField(prhs[0],0,"rc"));
	aV   =mxGetPr(mxGetField(prhs[0],0,"v"));
  // aDel =mxGetPr(mxGetField(prhs[0],0,"del"));
//	aGe  =mxGetPr(mxGetField(prhs[0],0,"ge"));    
 //   aOmg =mxGetPr(mxGetField(prhs[0],0,"omg"));
//	aDq  =mxGetPr(mxGetField(prhs[0],0,"dq"));	
//	aA	 =mxGetPr(mxGetField(prhs[0],0,"A"));	
//	tea  =mxGetField(prhs[0],0,"Con");	

	
   
   rows=(int)aN[0];  
    s3N[1]=rows;
	/// izlazne promenljive
   
	plhs[0] = mxCreateNumericArray(2, s31,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericArray(2, s31,mxDOUBLE_CLASS,mxREAL);    
    
   
    rcm  = mxGetPr(plhs[0]);
    vcm = mxGetPr(plhs[1]);    
    /// nadjemo prethodnika
	mass=0;
   for (k=0;k<rows;k++)
   {
	   vecaddcoef(rcm,aMs[k], aRc+3*k,rcm);
	   vecaddcoef(vcm,aMs[k], aV,vcm);
	   mass+=aMs[k];
	   aV+=3;
	   
   }
 
   vecmulscalar(rcm,1/mass, rcm);
   vecmulscalar(vcm,1/mass, vcm);
   
}
   