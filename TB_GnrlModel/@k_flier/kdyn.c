#include "mex.h"
#include "matops.h"


void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{    
   
    /// ulazne promenljive
    double *aN, *aCh, *aM, *aBs, *aCl, *aMs, *aJ, *aA, *aAlf,*aBet,*aGam,*aDel,*aOmg,*aV;
        
    /// izlazne promenljive
    double *aAc,*aA0, *aBc, *aB0;
    
    double *Nula;
    /// privremene promenljive
    int rows,c_chain,nCh;
    int rows3;
    int i,ip,c_si,ii;
    int currSeg;
    double *tmp, *T, *aTr, *tmat, *lam;
    mwSize   s3[3]={3,51,51};
    mwSize   s2[2]={3,51};
    
    aN   =mxGetPr(prhs[0]);
     aCh  =mxGetPr(prhs[1]);
    aM   =mxGetPr(prhs[2]);
    aBs  =mxGetPr(prhs[3]);
    aCl  =mxGetPr(prhs[4]);
	aMs  =mxGetPr(prhs[5]);
    aJ   =mxGetPr(prhs[6]);
	aA   =mxGetPr(prhs[7]);
	aAlf =mxGetPr(prhs[8]);
    aBet =mxGetPr(prhs[9]);    
    aGam =mxGetPr(prhs[10]);
    aDel =mxGetPr(prhs[11]);    
    aOmg =mxGetPr(prhs[12]);   
    aV   =mxGetPr(prhs[13]);

    rows=(int)aN[0];   
    nCh=(int)aCh[0];
    rows3=3*rows;
    s3[1]=s3[2]=rows;
    s2[1]=rows;
    
	/// izlazne promenljive
   
	plhs[0] = mxCreateNumericArray(3, s3,mxDOUBLE_CLASS,mxREAL);
    plhs[1] =  mxCreateNumericArray(2, s2,mxDOUBLE_CLASS,mxREAL);
    plhs[2] = mxCreateNumericArray(3, s3,mxDOUBLE_CLASS,mxREAL);
    plhs[3] =  mxCreateNumericArray(2, s2,mxDOUBLE_CLASS,mxREAL);   
   
    
   
 
	aAc = mxGetPr(plhs[0]); 
	aA0 = mxGetPr(plhs[1]); 
	aBc = mxGetPr(plhs[2]); 
	aB0 = mxGetPr(plhs[3]); ;
 
    

    tmp=mxCalloc(3, sizeof(double));
  	lam=mxCalloc(3, sizeof(double));
	T=mxCalloc(3*3,sizeof(double));
	aTr=mxCalloc(3*3,sizeof(double));
	tmat=mxCalloc(3*3,sizeof(double));
	Nula=mxCalloc(3, sizeof(double));
	zeros(Nula,3,1);
	for (c_chain=0;c_chain<nCh;c_chain++)        
    {         
         for(c_si=(int)aBs[c_chain];c_si<(int)aCl[c_chain];c_si++)           
         {             
             i=(int)aM[c_chain+nCh*c_si]-1;
             ip=(int)aM[c_chain+nCh*(c_si-1)]-1;     
             if (i>=3)
             {
				 tran3x3(aA+9*i,aTr);
				 matprod(aJ+9*i,aTr ,3,3,3,tmat);
				 matprod(aA+9*i,tmat,3,3,3,T);           
			     matprod(T,aOmg+3*i,3,3,1,tmp);
				 veccrossprod(aOmg+3*i,tmp,lam);

                 for (ii=0;ii<c_si;ii++)
                    {
                     currSeg=(int)aM[c_chain+nCh*ii]-1;     
                     vecaddcoef(Nula,-aMs[i],aBet+3*i+currSeg*rows3,aAc+3*i+currSeg*rows3);		
					 matprod(T,aAlf+3*ip+currSeg*rows3, 3,3,1,tmp);
					 vecsub(Nula, tmp,aBc+3*i+currSeg*rows3);
                    }
				matprod(T,aAlf+3*i+i*rows3,3,3,1,tmp);
				vecsub(Nula, tmp,aBc+3*i+i*rows3);
				vecaddcoef(Nula,-aMs[i],aBet+3*i+i*rows3,aAc+3*i+i*rows3);		
				
				/// racunanje a0
				vecaddcoef(Nula,-aMs[i],aDel+3*i,aA0+3*i);		
				matprod(T,aGam+3*i,3,3,1,tmp);
				vecadd(tmp,lam,tmp);
				vecsub(Nula,tmp,aB0+3*i);   
             }
                     
         }         
        
    }
   mxFree(tmp); 
  
}
    