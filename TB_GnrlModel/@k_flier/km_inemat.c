#include "mex.h"
#include "matops.h"


void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{    
   
    
    double *aN,*aCh,*aM, *aBs,*aCl,*aMs, *aA0, *aAc,*aB0, *aBc, *aRc, *aGe;

	double *H, *h0;
	double *G,*tmp,*tmp1;
    int c_chain, c_si, i, ip, ii,jj, currSeg,cs2, rows,nCh, rows3; 
	int dc; 
    mwSize   sNN[2]={45,45},sN1[2]={45,1};

	if(nrhs!=1)
      mexErrMsgTxt("One input required\n Calculates the inertia matrix and adjoint vector of the system\n  Usage :  [H, h0] = km_inemat(k_flier_object)");
    else if(nlhs > 2)
      mexErrMsgTxt("Too many output arguments.");
    else if(!mxIsClass(prhs[0],"k_flier"))
      mexErrMsgTxt("First input must be single flier object.");
    
    aN   =mxGetPr(mxGetField(prhs[0],0,"N"));
    aCh  =mxGetPr(mxGetField(prhs[0],0,"ch"));
    aM   =mxGetPr(mxGetField(prhs[0],0,"M"));
    aBs  =mxGetPr(mxGetField(prhs[0],0,"bs"));
    aCl  =mxGetPr(mxGetField(prhs[0],0,"cl"));
	aMs  =mxGetPr(mxGetField(prhs[0],0,"ms"));
	aGe  =mxGetPr(mxGetField(prhs[0],0,"ge"));
	aRc  =mxGetPr(mxGetField(prhs[0],0,"rc"));
	aAc = mxGetPr(mxGetField(prhs[0],0,"ac"));
	aA0 = mxGetPr(mxGetField(prhs[0],0,"ac0"));
	aBc = mxGetPr(mxGetField(prhs[0],0,"bc"));
	aB0 = mxGetPr(mxGetField(prhs[0],0,"bc0"));

    rows=(int)aN[0];   
    nCh=(int)aCh[0];

	G=mxCalloc(3,sizeof(double));
	tmp1=mxCalloc(3,sizeof(double));
	tmp=mxCalloc(3,sizeof(double));
	G[0]=G[1]=0; G[2]=-9.81;
	rows3=3*rows;
	sNN[0]=sNN[1]=sN1[0]=rows;

	plhs[0] = mxCreateNumericArray(2, sNN,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericArray(2, sN1,mxDOUBLE_CLASS,mxREAL);  
    H  = mxGetPr(plhs[0]);
    h0 = mxGetPr(plhs[1]);    
	dc=0;
   for (c_chain=0;c_chain<nCh;c_chain++)        
    {         
		if (c_chain==1)
			dc=1;
         for(c_si=(int)aBs[c_chain]-1+dc;c_si<(int)aCl[c_chain];c_si++)           
         {             
             i=(int)aM[c_chain+nCh*c_si]-1;
             for (ii=0;ii<=c_si;ii++)
                    {
						currSeg=(int)aM[c_chain+nCh*ii]-1;   
						if (currSeg>=3)
						{
							vecaddcoef(aA0+3*i, aMs[i], G, tmp1);
							veccrossprod(aRc+3*i+rows3*currSeg, tmp1, tmp);
							vecadd(tmp, aB0+3*i, tmp);
							h0[currSeg]-=vecdotprod(aGe+3*currSeg, tmp);
							for (jj=0;jj<=c_si;jj++)
							{
								cs2=(int)aM[c_chain+nCh*jj]-1;   
								 veccrossprod(aRc+3*i+rows3*currSeg, aAc+3*i+rows3*cs2,tmp);
								 vecadd(aBc+3*i+rows3*cs2,tmp, tmp);
								 H[currSeg+cs2*rows]-=vecdotprod(aGe+3*currSeg,tmp) ;
							}
							
						 }
						else
						{							
							vecaddcoef(aA0+3*i, aMs[i], G, tmp1);
							h0[currSeg]=h0[currSeg]-vecdotprod(aGe+3*currSeg, tmp1);
							for (jj=0;jj<=c_si;jj++)
							{
								 cs2=(int)aM[c_chain+nCh*jj]-1; 
								 H[currSeg+cs2*rows]-=vecdotprod(aGe+3*currSeg,aAc+3*i+rows3*cs2) ;
							}
						}
                 
                    }
                                        
         }         
        
    }
}