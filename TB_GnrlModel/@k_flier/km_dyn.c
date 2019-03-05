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
    
	if(nrhs!=1)
		mexErrMsgTxt("One input required\n calculates the dynamics coeffcients of the system\nUsage : km_dyn(k_flier_object)\n this function calculates a bunch of coefficients needed for later calculation\nof different dynamics-related quantities:\n* inertia matrices H and h0 <-> k_inemat");
	else if(nlhs > 1)
      mexErrMsgTxt("Too many output arguments.");

    else if(!mxIsClass(prhs[0],"k_flier"))
      mexErrMsgTxt("First input must be single flier object.");

	aN   =mxGetPr(mxGetField(prhs[0],0,"N"));
    aCh  =mxGetPr(mxGetField(prhs[0],0,"ch"));
    aM   =mxGetPr(mxGetField(prhs[0],0,"M"));
    aBs  =mxGetPr(mxGetField(prhs[0],0,"bs"));
    aCl  =mxGetPr(mxGetField(prhs[0],0,"cl"));
	aMs  =mxGetPr(mxGetField(prhs[0],0,"ms"));
    aJ   =mxGetPr(mxGetField(prhs[0],0,"J"));
	aA   =mxGetPr(mxGetField(prhs[0],0,"A"));
	aAlf =mxGetPr(mxGetField(prhs[0],0,"alf"));
    aBet =mxGetPr(mxGetField(prhs[0],0,"bet"));
    aGam =mxGetPr(mxGetField(prhs[0],0,"gam"));
    aDel =mxGetPr(mxGetField(prhs[0],0,"del"));
    aOmg =mxGetPr(mxGetField(prhs[0],0,"omg"));
    aV   =mxGetPr(mxGetField(prhs[0],0,"v"));

    rows=(int)aN[0];   
    nCh=(int)aCh[0];
    rows3=3*rows; 
    
	/// izlazne promenljive  
 
	aAc = mxGetPr(mxGetField(prhs[0],0,"ac"));
	aA0 = mxGetPr(mxGetField(prhs[0],0,"ac0"));
	aBc = mxGetPr(mxGetField(prhs[0],0,"bc"));
	aB0 = mxGetPr(mxGetField(prhs[0],0,"bc0"));
 
    

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
    