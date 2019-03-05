#include "mex.h"
#include "matrix.h"
#include "matops.h"
void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{    
   
    /// ulazne promenljive
    double *aN, *aCh, *aM, *aBs, *aCl, *aGe, *aDq, *aGpr,*aGdr;
        
    /// izlazne promenljive
    double *aAlf,*aBet,*aGam,*aDel,*aOmg, *aV;
    
   
    /// privremene promenljive
    int rows,c_chain,nCh;
    int rows3;
    int i,ip,c_si,ii;
    int currSeg;
    double *tmp, *tmp1, *difR;   
	if(nrhs!=2)
      mexErrMsgTxt("Two inputs required\n calculates the kinematics of the system ...\n Usage :  km_kine(k_kine_object, dq)");
    else if(nlhs > 0)
      mexErrMsgTxt("Too many output arguments.");
    else if(!mxIsClass(prhs[0],"k_flier"))
      mexErrMsgTxt("First input must be single flier object.");
    
	aN   =mxGetPr(mxGetField(prhs[0],0,"N"));
    aCh  =mxGetPr(mxGetField(prhs[0],0,"ch"));
    aM   =mxGetPr(mxGetField(prhs[0],0,"M"));
    aBs  =mxGetPr(mxGetField(prhs[0],0,"bs"));
    aCl  =mxGetPr(mxGetField(prhs[0],0,"cl"));
    aGe  =mxGetPr(mxGetField(prhs[0],0,"ge"));
    aDq  =mxGetPr(mxGetField(prhs[0],0,"dq"));
    aGpr =mxGetPr(mxGetField(prhs[0],0,"gpr"));
    aGdr =mxGetPr(mxGetField(prhs[0],0,"gdr"));
    rows=(int)aN[0];   
    nCh=(int)aCh[0];
    rows3=3*rows;  
    
	if (mxGetM(prhs[1])!= rows)
		mexErrMsgTxt("Coordinate vector length illegal error!");
	matcopy(aDq,mxGetPr(prhs[1]), rows);
	/// izlazne promenljive   
   
    aAlf = mxGetPr(mxGetField(prhs[0],0,"alf"));
    aBet = mxGetPr(mxGetField(prhs[0],0,"bet"));
    aGam = mxGetPr(mxGetField(prhs[0],0,"gam"));
    aDel = mxGetPr(mxGetField(prhs[0],0,"del"));
    aOmg = mxGetPr(mxGetField(prhs[0],0,"omg"));
    aV   = mxGetPr(mxGetField(prhs[0],0,"v"));
    
    aV[0]=aV[3]=aV[6]=aDq[0];
    aV[4]=aV[7]=aDq[1];
    aV[8]=aDq[2];
    
    
    aBet[0]=aBet[3]=aBet[6]=1;
    aBet[4+rows3]=aBet[7+rows3]=1;
    aBet[8+6*rows]=1;
    tmp=mxCalloc(3, sizeof(double));
    tmp1=mxCalloc(3, sizeof(double));
    difR=mxCalloc(3, sizeof(double));
    for (c_chain=0;c_chain<nCh;c_chain++)        
    {         
         for(c_si=(int)aBs[c_chain];c_si<(int)aCl[c_chain];c_si++)           
         {             
             i=(int)aM[c_chain+nCh*c_si]-1;
             ip=(int)aM[c_chain+nCh*(c_si-1)]-1;     
             if (i>=3)
             {
                 vecsub(aGdr+3*i,aGpr+3*ip+rows3*c_chain,difR);
                 for (ii=0;ii<c_si;ii++)
                    {
                     currSeg=(int)aM[c_chain+nCh*ii]-1;     
                     matcopy(aAlf+3*i+currSeg*rows3,aAlf+3*ip+currSeg*rows3,3);
                     veccrossprod(aAlf+3*ip+currSeg*rows3, difR,tmp);
                     vecadd(aBet+3*ip+currSeg*rows3,tmp,aBet+3*i+currSeg*rows3);
                    }
                matcopy(aAlf+3*i+i*rows3,aGe+3*i,3);
                veccrossprod(aGe+3*i,aGdr+3*i,aBet+ 3*i+i*rows3);     
                vecaddcoef(aOmg+3*ip,aDq[i],aGe+3*i,aOmg+3*i);          
                               
                veccrossprod(aOmg+3*ip,aGe+3*i,tmp);  
                vecaddcoef(aGam+3*ip,aDq[i],tmp,aGam+3*i);         

                veccrossprod(aOmg+3*ip, aGpr+3*ip+rows3*c_chain,tmp);
                veccrossprod(aOmg+3*i , aGdr+3*i                ,tmp1);                
                vecsub(aV+3*ip,tmp,aV+3*i);
                vecadd(aV+3*i,tmp1,aV+3*i);
                
                veccrossprod(aGam+3*ip,difR,tmp);
                vecadd(tmp,aDel+3*ip,aDel+3*i);
                
                veccrossprod(aOmg+3*ip,aGe+3*i,tmp);
                veccrossprod(tmp,aGdr+3*i, tmp1);
                vecaddcoef(aDel+3*i,aDq[i],tmp1,aDel+3*i);
                
                veccrossprod(aOmg+3*ip,aGpr+3*ip+rows3*c_chain,tmp);
                veccrossprod(aOmg+3*ip,tmp,tmp1);
                vecsub(aDel+3*i,tmp1,aDel+ 3*i);
                veccrossprod(aOmg+3*i,aGdr+3*i,tmp);
                veccrossprod(aOmg+3*i,tmp,tmp1);
                vecadd(aDel+3*i,tmp1,aDel+ 3*i);
                               
             }
                     
         }         
        
    }
   mxFree(tmp); 
   mxFree(tmp1); 
   mxFree(difR);
   
}
    