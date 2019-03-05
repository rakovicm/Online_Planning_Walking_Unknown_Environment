#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "matops.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{    
   
    /// ulazne promenljive
    double  *aN, *aCh,*aM,*aBs,*aCl,*lnk,*gsv,*aAlf,*aBet, *aGam, *aDel,*aGe, *aOmg, *aDq, *cnp, *cvec, *aA; 
    mxArray *tea, *t2;
    /// izlazne promenljive
    double *J,*JA;
   
   
    /// privremene promenljive
    int rows,nlnk,nCh,c_chain,c_si,i,ip,k,kk,currSeg, cn;
    double *tmp,*tmp1;
	 mwSize   s6N[3]={6,10};
    mwSize   s61[2]={6,1};
   
	 if(nrhs!=2)
      mexErrMsgTxt("Two inputs required\n [J A]=kjakPP(a,cn)\na -flier object \ncn-contact point number");
    else if(nlhs > 2)
      mexErrMsgTxt("Too many output arguments.");
    else if(!mxIsClass(prhs[0],"k_flier"))
      mexErrMsgTxt("First input must be flier object.");       
	

	aN   =mxGetPr(mxGetField(prhs[0],0,"N"));
    aCh  =mxGetPr(mxGetField(prhs[0],0,"ch"));
    aM   =mxGetPr(mxGetField(prhs[0],0,"M"));
    aBs  =mxGetPr(mxGetField(prhs[0],0,"bs"));
    aCl  =mxGetPr(mxGetField(prhs[0],0,"cl"));
    aAlf =mxGetPr(mxGetField(prhs[0],0,"alf"));
    aBet =mxGetPr(mxGetField(prhs[0],0,"bet"));
    aGam =mxGetPr(mxGetField(prhs[0],0,"gam"));
    aDel =mxGetPr(mxGetField(prhs[0],0,"del"));
	aGe  =mxGetPr(mxGetField(prhs[0],0,"ge"));    
    aOmg =mxGetPr(mxGetField(prhs[0],0,"omg"));
	aDq  =mxGetPr(mxGetField(prhs[0],0,"dq"));	
	aA	 =mxGetPr(mxGetField(prhs[0],0,"A"));	
	tea  =mxGetField(prhs[0],0,"Con");	

	cnp  =mxGetPr(prhs[1]);
	cn	 =(int) cnp[0];
	gsv  =mxCalloc(3,sizeof(double));
	if (mxIsStruct(tea))
	{	
		t2	=mxGetField(tea,0,"lnr");
		if (mxGetN(t2)<cn)
			 mexErrMsgTxt("Invalid contact number");
		lnk=mxGetPr(t2);
		nlnk=(int)lnk[cn-1];
		t2 =mxGetField(tea,0,"cvec");
		cvec=mxGetPr(t2);
		matprod(aA+9*(nlnk-1),cvec+3*(cn-1),3,3,1,gsv);
		//mexPrintf("%d %f %f %f",nlnk,cvec[0+3*(cn-1)],cvec[1+3*(cn-1)],cvec[2+3*(cn-1)]);
		if (nlnk<4)
			 mexErrMsgTxt(" 'lnk' is not allowed to be 1,2 or 3.");
	
	}
    
   
	 rows=(int)aN[0];   
	//mexPrintf("%d \n ",rows);
   
   
    nCh=(int)aCh[0];    
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