#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "matops.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{    
   
    /// ulazne promenljive
    double *aN,*aCh,*aM, *aBs,*aCl, *aQ, *aLne, *aAs0, *aLe, *aLdr, *aLpr;
    char *aOrs;
        
    /// izlazne promenljive
    double *aA,*aRc, *aGe, *aGd,*aGdr, *aGpr;
   
   
    /// privremene promenljive
    int rows,rows3,nCh,c_chain,c_si,i,ip,k,kk;
    double *TT,*tmp,*tvc;
   
    mwSize   s33N[3]={3,3,51};
    mwSize   s3NN[3]={3,51,51};
	mwSize   s3NC[3]={3,51,1};
    mwSize   s2[2]={3,51};
    
    aN   =mxGetPr(prhs[0]);
    aCh  =mxGetPr(prhs[1]);
    aM   =mxGetPr(prhs[2]);
    aBs  =mxGetPr(prhs[3]);
    aCl  =mxGetPr(prhs[4]);
	aOrs=mxCalloc(4,sizeof(char));
    mxGetString(prhs[5],aOrs,3);  
   
    //aOrs =mxGetPr(prhs[5]);
    aQ   =mxGetPr(prhs[6]);
    aLne =mxGetPr(prhs[7]);
    aAs0 =mxGetPr(prhs[8]);
    aLe  =mxGetPr(prhs[9]);
    aLdr =mxGetPr(prhs[10]);
    aLpr =mxGetPr(prhs[11]);

    rows=(int)aN[0];   
    nCh=(int)aCh[0];
    
    s3NC[1]=s2[1]=s3NN[1]=s3NN[2]=s33N[2]=rows;
	s3NC[2]=nCh;
    
	/// izlazne promenljive
   
	plhs[0] = mxCreateNumericArray(3, s33N,mxDOUBLE_CLASS,mxREAL);
    plhs[1] =  mxCreateNumericArray(3, s3NN,mxDOUBLE_CLASS,mxREAL);
    plhs[2] = mxCreateNumericArray(2, s2,mxDOUBLE_CLASS,mxREAL);
    plhs[3] =  mxCreateNumericArray(2, s2,mxDOUBLE_CLASS,mxREAL);
    plhs[4] =  mxCreateNumericArray(3, s3NC,mxDOUBLE_CLASS,mxREAL);
  
    
   
    aA  = mxGetPr(plhs[0]);
    aRc = mxGetPr(plhs[1]);    
    aGe = mxGetPr(plhs[2]);    
    aGdr= mxGetPr(plhs[3]);   
    aGpr= mxGetPr(plhs[4]);   
   
    eye(aA,3);
    TT=mxCalloc(3*3,sizeof(double));
    tmp=mxCalloc(3*3,sizeof(double));
    tvc=mxCalloc(3,sizeof(double));
    rows3=3*rows;
    for (c_chain=0;c_chain<nCh;c_chain++)        
    {         
        // printf("%d",rows);
        zeros(TT,3,3);
        i=(int)aM[c_chain +nCh*((int)aBs[c_chain]-1) ]-1;
        matprod(aA+9*i,aLpr+3*i+rows3*c_chain,3,3,1,aGpr+3*i+rows3*c_chain);
        for(c_si=(int)aBs[c_chain];c_si<(int)aCl[c_chain];c_si++)           
         {             
              
             i=(int)aM[c_chain+nCh*c_si]-1;
             ip=(int)aM[c_chain+nCh*(c_si-1)]-1;     
        
             if (i>=6)
             {                
                 //// rodrigova formula
               //  printf("podrigo");
                 for (k=0;k<3;k++)
                 {
                     zeros(TT+3*k,3,1);
                     veccrossprod(aLne+3*ip+rows3*c_chain,aAs0+3*k+9*i,tmp);
                     vecaddcoef(TT+3*k,sin(aQ[i]),tmp,TT+3*k);
                //     printf("%d",c_chain*rows3);
                 //  +3*ip+rows3*c_chain   veccrossprod(aLne+3*ip+rows3*c_chain,aAs0+3*k+9*i,tmp);
                     vecaddcoef(TT+k*3,(1-cos(aQ[i]))*vecdotprod(aLne+3*ip+rows3*c_chain,aAs0+3*k+9*i),aLne+3*ip+rows3*c_chain,TT+k*3);
                     vecaddcoef(TT+k*3,cos(aQ[i]),aAs0+3*k+9*i,TT+k*3);
                 }
                 matprod(aA+9*ip,TT,3,3,3,aA+9*i);
                 matprod(aA+9*i,aLe+3*i,3,3,1,aGe+3*i);
                 matprod(aA+9*i,aLdr+3*i,3,3,1,aGdr+3*i);
             }
             else
             {
               
                 if (i==1)
                 {
                     aGdr[0]=aQ[0];
                     eye(aA    ,3);
                     aRc[0]=aQ[0];
                     matcopy(aGe,aLe,3);
                     
                     aGdr[4]=aQ[1];
                     eye(aA+9*1,3);
                 }
                 else if (i==2)             
                 {
                      aGdr[8]=aQ[2];
                      eye(aA+9*2,3);                 
                 }
                 else
                 {
                   
                    if (aOrs[i-3]=='x')   
                         RotX(TT,aQ[i]);            
                    else if (aOrs[i-3]=='y')                      
                         RotY(TT,aQ[i]); 
                    else
                         RotZ(TT,aQ[i]);
                    matprod(aA+9*ip,TT,3,3,3,aA+9*i);
                 }
                 matprod(aA+9*i,aLe+3*i,3,3,1,aGe+3*i);
             }
             matprod(aA+9*i,aLpr+3*i+rows3*c_chain,3,3,1,aGpr+3*i+rows3*c_chain);   
             matcopy(aRc+3*i+rows3*i,aGdr+3*i,3);
             
             vecaddcoef(aGdr+3*i,-1,aGpr+3*ip+rows3*c_chain,tvc);
             for(k=0;k<=c_si-1;k++)
             {
                 kk=(int)aM[c_chain+nCh*k]-1;
                 vecadd(aRc+3*ip+rows3*kk,tvc,aRc+3*i+rows3*kk);
             }
                     
         }     
        for(i=5;i<(int)aBs[c_chain]-1;i++)
        {
            k=(int)aM[c_chain+nCh*i]-1;              
            matprod(aA+9*k,aLpr+3*k+rows3*c_chain,3,3,1,aGpr+3*k+rows3*c_chain);    
        }
        
    }
  mxFree(TT);   
  mxFree(aOrs);
  mxFree(tmp); 
  mxFree(tvc); 

}