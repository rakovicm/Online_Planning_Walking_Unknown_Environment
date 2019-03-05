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

	if(nrhs!=2)
      mexErrMsgTxt("Two inputs required\n Calculates the geometry and configuration of the mechanism\n Usage : km_geo( k_flier_object, coordinate vector q )");
    else if(nlhs > 0)
      mexErrMsgTxt("Too many output arguments.");

    else if(!mxIsClass(prhs[0],"k_flier"))
      mexErrMsgTxt("First input must be single flier object.");
    
    aN   =mxGetPr(mxGetField(prhs[0],0,"N"));
    aCh  =mxGetPr(mxGetField(prhs[0],0,"ch"));
    aM   =mxGetPr(mxGetField(prhs[0],0,"M"));
    aBs  =mxGetPr(mxGetField(prhs[0],0,"bs"));
    aCl  =mxGetPr(mxGetField(prhs[0],0,"cl"));
	aOrs=mxCalloc(4,sizeof(char));
    mxGetString(mxGetField(prhs[0],0,"ors"),aOrs,3);  
   
    //aOrs =mxGetPr(prhs[5]);
    aQ   =mxGetPr(mxGetField(prhs[0],0,"q"));
    aLne =mxGetPr(mxGetField(prhs[0],0,"lne"));
    aAs0 =mxGetPr(mxGetField(prhs[0],0,"As0"));
    aLe  =mxGetPr(mxGetField(prhs[0],0,"le"));
    aLdr =mxGetPr(mxGetField(prhs[0],0,"ldr"));
    aLpr =mxGetPr(mxGetField(prhs[0],0,"lpr"));

    rows=(int)aN[0];   
    nCh=(int)aCh[0];
    
	if (mxGetM(prhs[1])!= rows)
		mexErrMsgTxt("Coordinate vector length illegal error!");

	matcopy(aQ,mxGetPr(prhs[1]), rows);
	/// izlazne promenljive    
   
    aA  = mxGetPr(mxGetField(prhs[0],0,"A"));
    aRc = mxGetPr(mxGetField(prhs[0],0,"rc"));
    aGe = mxGetPr(mxGetField(prhs[0],0,"ge"));
    aGdr= mxGetPr(mxGetField(prhs[0],0,"gdr"));
    aGpr= mxGetPr(mxGetField(prhs[0],0,"gpr"));

   
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