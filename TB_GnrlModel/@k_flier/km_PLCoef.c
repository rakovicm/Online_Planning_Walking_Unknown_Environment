#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "matops.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{   	 
    /// ulazne promenljive
    double *aN, *aCh,*aM, *aBs, *aCl, *aRc, *aV, *aAc, *aBc,*r0;
        
    /// izlazne promenljive
    double *PC,*LC;
   
   
    /// privremene promenljive
    int rows,nCh,chain, c_si, c_pi, i, ip,ip3,i3,irow, row3,ind;
    double *tmp,*tmp2;
   
    mwSize   s3N[2]={3,10};

	if(nrhs!=2)
      mexErrMsgTxt("Two inputs required\n [Pc Lc]=kk_PLCoef(a,r)\n  Racuna koeficiente Pc Lc koji sluze za racunanje kolicine kretanja\n  i momenta kolicine kretanja za tacku r\n Kolicina kretanja:\n     P=Pc*q\n Moment kolicine kretanja:\n L=Lc*q");
    else if(nlhs > 2)
      mexErrMsgTxt("Too many output arguments.");
    else if(!mxIsClass(prhs[0],"k_flier"))
      mexErrMsgTxt("First input must be single flier object.");
    
    
    aN   =mxGetPr(mxGetField(prhs[0],0,"N"));
    aCh  =mxGetPr(mxGetField(prhs[0],0,"ch"));
    aM   =mxGetPr(mxGetField(prhs[0],0,"M"));
    aBs  =mxGetPr(mxGetField(prhs[0],0,"bs"));
    aCl  =mxGetPr(mxGetField(prhs[0],0,"cl"));
	aRc	 =mxGetPr(mxGetField(prhs[0],0,"rc"));
	aAc  =mxGetPr(mxGetField(prhs[0],0,"ac"));
    aBc  =mxGetPr(mxGetField(prhs[0],0,"bc"));
    r0   =mxGetPr(prhs[1]);
   
   
    rows=(int)aN[0];   
    nCh =(int)aCh[0];
    s3N[1]=rows;
	/// izlazne promenljive
   
	plhs[0] = mxCreateNumericArray(2, s3N,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericArray(2, s3N,mxDOUBLE_CLASS,mxREAL);
   
   
    
   
    PC  = mxGetPr(plhs[0]);
	zeros(PC,3,rows);
    LC  = mxGetPr(plhs[1]);  
	zeros(LC,3,rows);
	tmp=mxCalloc(3,sizeof(double));
	tmp2=mxCalloc(3*3,sizeof(double));
	row3=3*rows;
	for (chain=0; chain<nCh; chain++)
	{	
	    for(c_si=(int)aBs[chain];c_si<(int)aCl[chain];c_si++)       
	
		{
			i=(int)aM[chain+nCh*c_si]-1;
			i3=3*i;
			vecsub(aRc+i3,r0,aRc+i3);			
			CrossMat(aRc+i3,tmp2);			
		//	mexPrintf("\n %d",i);
			for (c_pi=0;c_pi<=c_si; c_pi++)
			{
				ip=(int)aM[chain+nCh*c_pi]-1;		
				ip3=ip*3;
				ind=row3*ip+i3;
			    vecsub(PC+ip3,aAc+ind,PC+ip3);
				vecsub(LC+ip3,aBc+ind,LC+ip3);
		//		zeros(tmp2,3,3);
				matprod(tmp2  ,aAc+ind,3,3,1,tmp);				
		//		mexPrintf("%2d %8.5f %8.5f %8.5f\n",ip,tmp[0],tmp[1],tmp[2]);
				
				vecsub(LC+ip3,tmp,LC+ip3);
				
			}
		}
	} 
      mxFree(tmp);    
}