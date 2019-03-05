#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "matops.h"
void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{    
   
    /// ulazne promenljive
    double *aN, *aCh,*aM,*aBs,*aCl,*lnk,*gsv,*aAlf,*aBet, *aGam, *aDel,*aGe, *aOmg, *aDq; 
        
    /// izlazne promenljive
    double *J,*JA;
   
   
    /// privremene promenljive
    int rows,nlnk,nCh,c_chain,c_si,i,ip,k,kk,currSeg;
    double *tmp,*tmp1;
   
    mwSize   s6N[3]={6,10};
    mwSize   s61[2]={6,1};
    
    
    aN   =mxGetPr(prhs[0]);
    aCh  =mxGetPr(prhs[1]);
    aM   =mxGetPr(prhs[2]);
    aBs  =mxGetPr(prhs[3]);
    aCl  =mxGetPr(prhs[4]);
    aAlf =mxGetPr(prhs[5]);
    aBet =mxGetPr(prhs[6]);	
    aGam =mxGetPr(prhs[7]);	
    aDel =mxGetPr(prhs[8]);	            
    aGe  =mxGetPr(prhs[9]);
	aOmg =mxGetPr(prhs[10]);
	aDq  =mxGetPr(prhs[11]);
    lnk  =mxGetPr(prhs[12]);
    gsv  =mxGetPr(prhs[13]);
   
   
    rows=(int)aN[0];   
    nCh=(int)aCh[0];
    nlnk=(int)lnk[0];
    
    s6N[1]=rows;
	/// izlazne promenljive
   
	plhs[0] = mxCreateNumericArray(2, s6N,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericArray(2, s61,mxDOUBLE_CLASS,mxREAL);
   
   
    
   
    J  = mxGetPr(plhs[0]);
    JA = mxGetPr(plhs[1]);    
    /// nadjemo prethodnika
    ip=-1;
    for (c_chain=0;((c_chain<nCh) && (ip==-1));c_chain++)
	{
       for(c_si=(int)aBs[c_chain];c_si<(int)aCl[c_chain];c_si++)     
           if ((int)aM[c_chain +nCh*((int)c_si) ]==nlnk)
           {
                ip=(int)aM[c_chain+nCh*(c_si-1)]-1;     			
                break;
           }
	}
   for(k=0;k<rows;k++)
   {
        matcopy(J+6*k,aBet+3*(nlnk-1)+3*rows*k,3);
       matcopy(J+3+6*k,aAlf+3*(nlnk-1)+3*rows*k,3);
   }
    tmp=mxCalloc(3,sizeof(double));
    tmp1=mxCalloc(3,sizeof(double));

	c_chain=c_chain-1;
   for (k=0;k<c_si;k++)
   {
	    currSeg=(int)aM[c_chain+nCh*k]-1;      
		veccrossprod(aAlf+3*ip+3*rows*currSeg,gsv,tmp);
        vecadd(J+6*currSeg,tmp,J+6*currSeg);        
   }
   veccrossprod(aGe+3*(nlnk-1),gsv,tmp);
   vecadd(J+6*(nlnk-1),tmp,J+6*(nlnk-1));        
   
   matcopy(JA,aDel+3*(nlnk-1),3);
   veccrossprod(aGam+3*ip,gsv,tmp);
   vecadd(JA,tmp,JA);
   
   veccrossprod(aOmg+3*ip,aGe+3*(nlnk-1),tmp1);
   veccrossprod(tmp1,gsv,tmp);
   vecaddcoef(JA,aDq[nlnk-1],tmp,JA);
   veccrossprod(aOmg+3*(nlnk-1),gsv,tmp1);
   veccrossprod(aOmg+3*(nlnk-1),tmp1,tmp);
   vecadd(JA,tmp,JA);

   matcopy(JA+3,aGam+3*(nlnk-1),3);
   mxFree(tmp); 
   mxFree(tmp1);  
    
}