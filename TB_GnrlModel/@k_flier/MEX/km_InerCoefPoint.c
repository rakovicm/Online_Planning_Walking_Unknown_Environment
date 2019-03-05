#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "matops.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{    
   
    /// ulazne promenljive
    double  *aN,*aA,*aJ,*aAlf, *aBet, *aGam, *aDel, *aMs, *aRc, *aOmg;

	double *ar;
    /// izlazne promenljive
    double *am,*a0m;
	double *bm,*b0m;
   
   
    /// privremene promenljive
	int rows;
	int i,j;

	 double *tmp,*tmp1,*tmp2, *cvec;

    //int rows,nlnk,nCh,c_chain,c_si,i,ip,k,kk,currSeg, cn;
    //double *tmp,*tmp1;
	mwSize   s3N[3]={3,10};
    mwSize   s31[2]={3,1};
   
	 if(nrhs!=2)
      mexErrMsgTxt("Two inputs required\n [am a0m bm b0m]=km_jakPP(a,r)\na -flier object \nr-point coordinate\n calculates coefficients to so that\n Inertial force equals Fi=am*q''+am0\nInertial moment Mi=bm*q''+bm0 ");
    else if(nlhs > 4)
      mexErrMsgTxt("Too many output arguments.");
    else if(!mxIsClass(prhs[0],"k_flier"))
      mexErrMsgTxt("First input must be flier object.");       
	

	aN   =mxGetPr(mxGetField(prhs[0],0,"N"));
	aMs  =mxGetPr(mxGetField(prhs[0],0,"ms"));
    //aM   =mxGetPr(mxGetField(prhs[0],0,"M"));
   // aBs  =mxGetPr(mxGetField(prhs[0],0,"bs"));
 //   aCl  =mxGetPr(mxGetField(prhs[0],0,"cl"));
    aAlf =mxGetPr(mxGetField(prhs[0],0,"alf"));
    aBet =mxGetPr(mxGetField(prhs[0],0,"bet"));
    aGam =mxGetPr(mxGetField(prhs[0],0,"gam"));
    aDel =mxGetPr(mxGetField(prhs[0],0,"del"));
//	aGe  =mxGetPr(mxGetField(prhs[0],0,"ge"));    
    aOmg =mxGetPr(mxGetField(prhs[0],0,"omg"));
//	aDq  =mxGetPr(mxGetField(prhs[0],0,"dq"));	
	aA	 =mxGetPr(mxGetField(prhs[0],0,"A"));	
	aJ	 =mxGetPr(mxGetField(prhs[0],0,"J"));
	aRc  =mxGetPr(mxGetField(prhs[0],0,"rc"));
	ar  =mxGetPr(prhs[1]);
	
	 rows=(int)aN[0];   

    s3N[1]=rows;
	/// izlazne promenljive
   
	plhs[0] = mxCreateNumericArray(2, s3N,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericArray(2, s31,mxDOUBLE_CLASS,mxREAL);
   	plhs[2] = mxCreateNumericArray(2, s3N,mxDOUBLE_CLASS,mxREAL);
    plhs[3] = mxCreateNumericArray(2, s31,mxDOUBLE_CLASS,mxREAL);
   
    am  = mxGetPr(plhs[0]);
    a0m = mxGetPr(plhs[1]); 
	bm  = mxGetPr(plhs[2]);
    b0m = mxGetPr(plhs[3]);    
   
	tmp=mxCalloc(9,sizeof(double));
	tmp1=mxCalloc(9,sizeof(double));
	tmp2=mxCalloc(9,sizeof(double));
	cvec=mxCalloc(3,sizeof(double));
	for ( i=0; i<rows; i++)
	{
		matprod(aA+9*i,aJ+9*i, 3,3,3,tmp);
		tran3x3(aA+9*i, tmp1);
		matprod(tmp,tmp1,3,3,3,tmp2);
		vecsub(aRc+3*i, ar,cvec);
		for ( j=0; j<rows; j++)
		{
			vecaddcoef(am+3*j,-aMs[i],aBet+3*(i)+3*rows*j, am+3*j);
			matprod(tmp2, aAlf+3*i+3*rows*j, 3,3,1, tmp);
			vecsub(bm+3*j,tmp,bm+3*j); 
			veccrossprod(cvec, aBet+3*i+3*rows*j, tmp);
			vecaddcoef(bm+3*j,-aMs[i],tmp,bm+3*j);
		}

		matprod(tmp2,aGam+3*i,3,3,1,tmp);
		vecsub(b0m, tmp, b0m);
		
		matprod(tmp1, aOmg+3*i, 3,3,1, tmp);
		matprod(aJ+9*i, tmp, 3,3,1,tmp+3);
		veccrossprod(tmp,tmp+3, tmp+6);
		matprod(aA+9*i, tmp+6, 3,3,1, tmp);
		vecsub(b0m, tmp, b0m);

		veccrossprod(cvec,  aDel+3*i, tmp);
		vecaddcoef(b0m, -aMs[i], tmp ,b0m);
		vecaddcoef(a0m, -aMs[i], aDel+3*i ,a0m);
	}

    mxFree(tmp); 
    mxFree(tmp1); 
	mxFree(tmp2); 
	mxFree(cvec); 
    
}