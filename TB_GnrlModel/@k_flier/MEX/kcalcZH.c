#include "mex.h"
#include "matrix.h"
#include "matops.h"


void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{    
   
    /// ulazne promenljive
    double *aN, *aCh, *aM, *aBs, *aCl, *aRc, *aA, *aAlf,*aBet,*aGam,*aDel,*aOmg,*aMs,*aJ,*ra;   
    /// izlazne promenljive
    double *ZH,*Zh0;
    
   
    /// privremene promenljive
    int rows,c_chain,nCh;
    int rows3;
    int i,c_si,kk,ip;
	int dSi;
    //int currSeg;
    double *tmp, *tmp1, *tmp2, *tmp3, *Atran;
	double *tm, *Jmod, *gammod, *rcnov, *G, *Alfmod;
    mwSize   s3N[2]={3,51};
	mwSize   s31[2]={3,1};
    
    aN   =mxGetPr(prhs[0]);
    aCh  =mxGetPr(prhs[1]);
    aM   =mxGetPr(prhs[2]);
    aBs  =mxGetPr(prhs[3]);
    aCl  =mxGetPr(prhs[4]);
    aRc  =mxGetPr(prhs[5]);
	aA   =mxGetPr(prhs[6]);
	aAlf =mxGetPr(prhs[7]);
	aBet =mxGetPr(prhs[8]);
	aGam =mxGetPr(prhs[9]);
	aDel =mxGetPr(prhs[10]);
	aOmg =mxGetPr(prhs[11]);
	aMs  =mxGetPr(prhs[12]);
	aJ   =mxGetPr(prhs[13]);
    ra   =mxGetPr(prhs[14]);
 
    rows=(int)aN[0];   
    nCh=(int)aCh[0];
    rows3=3*rows;
    s3N[1]=rows;
    
    
	/// izlazne promenljive
   
	plhs[0] = mxCreateNumericArray(2, s3N,mxDOUBLE_CLASS,mxREAL);
    plhs[1] =  mxCreateNumericArray(2, s31,mxDOUBLE_CLASS,mxREAL);
    ZH = mxGetPr(plhs[0]);
    Zh0 = mxGetPr(plhs[1]);   
   
 
    tmp=mxCalloc(3, sizeof(double));
    tmp1=mxCalloc(3, sizeof(double));
	tmp2=mxCalloc(3, sizeof(double));
	tmp3=mxCalloc(3, sizeof(double));
	rcnov=mxCalloc(3, sizeof(double));	   
	gammod=mxCalloc(3, sizeof(double));
	G=mxCalloc(3, sizeof(double));	

	Atran=mxCalloc(3*3,sizeof(double));
	Alfmod=mxCalloc(3,sizeof(double));
	tm=mxCalloc(3*3,sizeof(double));
	Jmod=mxCalloc(3*3,sizeof(double));
	
	dSi=-1;
	G[0]=G[1]=0;
	G[2]=-9.81;

   for (c_chain=0;c_chain<nCh;c_chain++)        
    {         
         for(c_si=(int)aBs[c_chain]+dSi;c_si<(int)aCl[c_chain];c_si++)           
         {             
             i=(int)aM[c_chain+nCh*c_si]-1;       
			 
			 tran3x3(aA+9*i,Atran);
 			 matprod(aA+9*i,aJ+9*i,3,3,3,tm);
			 matprod(tm,Atran,3,3,3,Jmod);

			 vecsub(aRc+3*i,ra,rcnov);
			 // racunamo gama modifikovano
			 matprod(Jmod,aGam+3*i,3,3,1,tmp);
			 matprod(Atran, aOmg+3*i, 3,3,1,tmp1);
			 matprod(aJ+9*i,tmp1, 3,3,1,tmp2);
			 veccrossprod(tmp1,tmp2,tmp3);
			 matprod(aA+9*i, tmp3, 3,3,1,tmp2);
			 vecadd(tmp2, tmp,gammod);

			 // racunamo Zh0
			 vecsub(G,aDel+3*i,tmp);
			 veccrossprod(rcnov, tmp, tmp1);
			 vecaddcoef(Zh0,aMs[i],tmp1,tmp2);
			 vecsub(tmp2, gammod,Zh0);		
			 for (kk=0; kk<=c_si; kk++)
			 {
				 // izracunam alf mod
				 
				 ip=(int)aM[c_chain+nCh*kk]-1;
				 matprod(Jmod,aAlf+3*i+rows3*ip,3,3,1,Alfmod);
				 veccrossprod(rcnov, aBet+3*i+rows3*ip,tmp);
				 vecaddcoef(ZH+3*ip,-aMs[i],tmp,tmp1);
				 vecsub(tmp1,Alfmod,ZH+3*ip);
				 

			 }        
         }     
		 dSi=0;
        
    }
   
	mxFree(tmp); 
    mxFree(tmp1);
	mxFree(tmp2);
	mxFree(tmp3);
	mxFree(rcnov);
	mxFree(gammod);
	mxFree(G);
	mxFree(Atran);
	mxFree(Alfmod);
	mxFree(tm);
	mxFree(Jmod);
   
   
}
    