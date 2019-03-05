#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "matops.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{    
	double *aA,*aRc,*cnp, *lnk;
	char* ors;
	double *cvec,*cQ, *gsv, *Res;
	double *tmat;
	mxArray *tea,*t2;
	int cn, nlnk;
	mwSize   s61[2]={6,1};

	if(nrhs!=3)
      mexErrMsgTxt("Three inputs required\n Calculates the world coordinates of one point of the link\n Usage : X = km_linkX(k_flier, cn, rorder)\nX - position of the contact point of the link\nk_flier - an instance of the k_flier object\ncn - inedex number of the flier contact\nrorder - a string describing the rotation order - p.e. 'xyz'");
    else if(nlhs > 1)
      mexErrMsgTxt("Too many output arguments.");
    else if(!mxIsClass(prhs[0],"k_flier"))
      mexErrMsgTxt("First input must be single flier object.");

	aA  = mxGetPr(mxGetField(prhs[0],0,"A"));
    aRc = mxGetPr(mxGetField(prhs[0],0,"rc"));

	plhs[0] = mxCreateNumericArray(2, s61,mxDOUBLE_CLASS,mxREAL);
	Res		= mxGetPr(plhs[0]);

	tea  =mxGetField(prhs[0],0,"Con");	
	cnp  =mxGetPr(prhs[1]);
	ors =mxCalloc(3,sizeof(char));
    mxGetString(prhs[2],ors,4);  
	cn	 =(int) cnp[0];
	gsv =mxCalloc(3,sizeof(double));
	tmat=mxCalloc(3*3,sizeof(double));
	if (mxIsStruct(tea))
	{	
		t2	=mxGetField(tea,0,"lnr");
		if (mxGetN(t2)<cn)
			 mexErrMsgTxt("Invalid contact number");
		lnk=mxGetPr(t2);
		nlnk=(int)lnk[cn-1];
		if (nlnk<4)
			 mexErrMsgTxt(" 'lnk' is not allowed to be 1,2 or 3.");

		t2 =mxGetField(tea,0,"cvec");
		cvec=mxGetPr(t2);
		matprod(aA+9*(nlnk-1),cvec+3*(cn-1),3,3,1,gsv);
		vecadd(aRc+3*(nlnk-1),gsv,Res);
		cQ =mxGetPr(mxGetField(tea,0,"cQ"));
		matprod(aA+9*(nlnk-1),cQ+9*(cn-1),3,3,3,tmat);
		
		if (!strcmp(ors,"xyz"))
			TRotXYZ(tmat,Res+3);
		else if (!strcmp(ors,"xzy"))		
			TRotXZY(tmat,Res+3);		
		else if (!strcmp(ors,"yxz"))		
			TRotYXZ(tmat,Res+3);		
		else if (!strcmp(ors,"yzx"))		
			TRotYZX(tmat,Res+3);		
		else if (!strcmp(ors,"zxy"))		
			TRotZXY(tmat,Res+3);		
		else if (!strcmp(ors,"zyx"))		
			TRotZYX(tmat,Res+3);		
		else
			mexErrMsgTxt("The requested order ''%s'' is not supported!",ors);
		//mexPrintf("%d %f %f %f",nlnk,cvec[0+3*(cn-1)],cvec[1+3*(cn-1)],cvec[2+3*(cn-1)]);
		
	
	}
}
