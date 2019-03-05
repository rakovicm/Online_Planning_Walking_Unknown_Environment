#include "mex.h"
#include "matrix.h"
#include "matops.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{    
   
    /// ulazne promenljive
	double *aN, *aCh, *aM, *aBs, *aCl, *aGe, *aDq, *aGpr,*aGdr,*aJ, *aA, *aMs;
        
    /// izlazne promenljive
    double *aAlf,*aBet,*aGam,*aDel,*omg, *aV, *aAc,*aA0, *aBc, *aB0;
    
	double *Nula;
   
    /// privremene promenljive
    int rows,c_chain,nCh;
    int rows3;
    int i,ip,c_si,ii;
    int currSeg;
    double *tmp, *tmp1, *difR, *T, *aTr, *tmat, *lam;
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
    aGe  =mxGetPr(prhs[8]);
    aDq  =mxGetPr(prhs[9]);
    aGpr =mxGetPr(prhs[10]);
    aGdr =mxGetPr(prhs[11]);
    rows=(int)aN[0];   
    nCh=(int)aCh[0];
    rows3=3*rows;
    s3[1]=rows;s3[2]=rows;
    s2[1]=rows;
    
	/// izlazne promenljive
   
	plhs[0] = mxCreateNumericArray(3, s3,mxDOUBLE_CLASS,mxREAL);
    plhs[1] =  mxCreateNumericArray(3, s3,mxDOUBLE_CLASS,mxREAL);
    plhs[2] = mxCreateNumericArray(2, s2,mxDOUBLE_CLASS,mxREAL);
    plhs[3] =  mxCreateNumericArray(2, s2,mxDOUBLE_CLASS,mxREAL);
    plhs[4] = mxCreateNumericArray(2, s2,mxDOUBLE_CLASS,mxREAL);
    plhs[5] =  mxCreateNumericArray(2, s2,mxDOUBLE_CLASS,mxREAL);
	plhs[6] =  mxCreateNumericArray(3, s3,mxDOUBLE_CLASS,mxREAL);
	plhs[7] =  mxCreateNumericArray(2, s2,mxDOUBLE_CLASS,mxREAL);
	plhs[8] =  mxCreateNumericArray(3, s3,mxDOUBLE_CLASS,mxREAL);
	plhs[9] =  mxCreateNumericArray(2, s2,mxDOUBLE_CLASS,mxREAL);
    //plhs[4] = mxCreateDoubleMatrix(rowsa,collb,mxREAL);
    
   
    aAlf = mxGetPr(plhs[0]);
    aBet = mxGetPr(plhs[1]);    
    aGam = mxGetPr(plhs[2]);
    aDel = mxGetPr(plhs[3]);    
    omg = mxGetPr(plhs[4]);   
    aV= mxGetPr(plhs[5]); 
	aAc = mxGetPr(plhs[6]); 
	aA0 = mxGetPr(plhs[7]); 
	aBc = mxGetPr(plhs[8]); 
	aB0 = mxGetPr(plhs[9]); ;
    
    aV[0]=aV[3]=aV[6]=aDq[0];
    aV[4]=aV[7]=aDq[1];
    aV[8]=aDq[2];
    
    
    aBet[0]=aBet[3]=aBet[6]=1;
    aBet[4+rows3]=aBet[7+rows3]=1;
    aBet[8+6*rows]=1;
    tmp=mxCalloc(3, sizeof(double));
    tmp1=mxCalloc(3, sizeof(double));
    difR=mxCalloc(3, sizeof(double));
	lam=mxCalloc(3, sizeof(double));
	Nula=mxCalloc(3, sizeof(double));
	zeros(Nula,3,1);
	T=mxCalloc(3*3,sizeof(double));
	aTr=mxCalloc(3*3,sizeof(double));
	tmat=mxCalloc(3*3,sizeof(double));
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
                 vecsub(aGdr+3*i,aGpr+3*ip+rows3*c_chain,difR);
				 vecaddcoef(omg+3*ip,aDq[i],aGe+3*i,omg+3*i);          
				
				 matprod(T,omg+3*i,3,3,1,tmp);
				 veccrossprod(omg+3*i,tmp,lam);

                 for (ii=0;ii<c_si;ii++)
                    {
                     currSeg=(int)aM[c_chain+nCh*ii]-1;     
                     matcopy(aAlf+3*i+currSeg*rows3,aAlf+3*ip+currSeg*rows3,3);
                     veccrossprod(aAlf+3*ip+currSeg*rows3, difR,tmp);
                     vecadd(aBet+3*ip+currSeg*rows3,tmp,aBet+3*i+currSeg*rows3);
					 vecaddcoef(Nula,-aMs[i],aBet+3*i+currSeg*rows3,aAc+3*i+currSeg*rows3);		
					 matprod(T,aAlf+3*ip+currSeg*rows3, 3,3,1,tmp);
					 vecsub(Nula, tmp,aBc+3*i+currSeg*rows3);
                    }
			    matcopy(aAlf+3*i+i*3*rows,aGe+3*i,3);
                veccrossprod(aGe+3*i,aGdr+3*i,aBet+ 3*i+i*3*rows);
               
				matprod(T,aAlf+3*i+i*rows3,3,3,1,tmp);
				vecsub(Nula, tmp,aBc+3*i+i*rows3);
				vecaddcoef(Nula,-aMs[i],aBet+3*i+i*rows3,aAc+3*i+i*rows3);		
					
                               
                veccrossprod(omg+3*ip,aGe+3*i,tmp);  
                vecaddcoef(aGam+3*ip,aDq[i],tmp,aGam+3*i);         

                veccrossprod(omg+3*ip, aGpr+3*ip+rows3*c_chain,tmp);
                veccrossprod(omg+3*i , aGdr+3*i                ,tmp1);                
                vecsub(aV+3*ip,tmp,aV+3*i);
                vecadd(aV+3*i,tmp1,aV+3*i);
                
                veccrossprod(aGam+3*ip,difR,tmp);
                vecadd(tmp,aDel+3*ip,aDel+3*i);
                
                veccrossprod(omg+3*ip,aGe+3*i,tmp);
                veccrossprod(tmp,aGdr+3*i, tmp1);
                vecaddcoef(aDel+3*i,aDq[i],tmp1,aDel+3*i);
                
                veccrossprod(omg+3*ip,aGpr+3*ip+rows3*c_chain,tmp);
                veccrossprod(omg+3*ip,tmp,tmp1);
                vecsub(aDel+3*i,tmp1,aDel+ 3*i);
                veccrossprod(omg+3*i,aGdr+3*i,tmp);
                veccrossprod(omg+3*i,tmp,tmp1);
                vecadd(aDel+3*i,tmp1,aDel+3*i);
				/// racunanje a0
				vecaddcoef(aA0+3*i,-aMs[i],aDel+3*i,aA0+3*i);		
				matprod(T,aGam+3*i,3,3,1,tmp);
				vecadd(tmp,lam,tmp);
				vecsub(aB0+3*i,tmp,aB0+3*i);                    
             }
                     
         }         
        
    }
   mxFree(tmp); 
   mxFree(tmp1); 
   mxFree(difR);
   
}
    