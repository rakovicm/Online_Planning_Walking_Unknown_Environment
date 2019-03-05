#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "matops.h"

double sqr(double a) {return a*a;}

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{    
   
    /// ulazne promenljive
//	, *aCh, *aM, *aBs, *aCl, *aGe, *aDq, *aGpr,*aGdr,*aJ, *aA, *aMs;      
 //   double *aAlf,*aBet,*aGam,*aDel,*aOmg, *aV, *aAc,*aA0, *aBc, *aB0;
	double *aN, *aA, *aRc, *aLnk, *aCvec, *aOmg, *aV,*aAlf,*aBet;
	int nCP, nCF;
	double *COP,*del;
	double Cn,Kn,mu;
	// izlazne promenljive 
	double *delP, *tau, *F, *lam;    

    /// privremene promenljive
    int rows,c_cp,nCh, nlnk;
	double *tmpBet,*tmpAlf;
	
	 mwSize   sN1[2]={6,1};
	 mwSize	  s3NCP[2]={3,1};
	 mwSize   sNCP1[2]={3,1};
	 mxArray  *tmpA, *aCon;

    int cp,i;
    int currSeg;
	double rc[3], vc[3] ,gsv[3],tmpV[3], M[3];
    double *tmp; //, *tmp1, *difR, *T, *aTr, *tmat, *lam;


	if(nrhs!=6)
      mexErrMsgTxt("Six inputs required\n calculates the kinematics and dynamics of the system ...\n Usage :  km_CalcContacts(flier, CONPOINTS,del,Kn,Cn,mi)");
    else if(nlhs !=4)
      mexErrMsgTxt("Need four outputs.");
    else if(!mxIsClass(prhs[0],"k_flier"))
      mexErrMsgTxt("First input must be single flier object.");
    
    
    aN		=	mxGetPr(mxGetField(prhs[0],0,"N"));
	aCon	=	mxGetField(prhs[0],0,"Con");
	tmpA	=	mxGetField(aCon,0,"lnr");
    nCF		=   (int)mxGetN(tmpA);
	aLnk	=	mxGetPr(tmpA);
	aCvec   =   mxGetPr(mxGetField(aCon,0,"cvec"));
	
 //   aCh  =mxGetPr(mxGetField(prhs[0],0,"ch"));
  //  aM   =mxGetPr(mxGetField(prhs[0],0,"M"));
  //  aBs  =mxGetPr(mxGetField(prhs[0],0,"bs"));
 //   aCl  =mxGetPr(mxGetField(prhs[0],0,"cl"));
//	aMs  =mxGetPr(mxGetField(prhs[0],0,"ms"));
  //  aJ   =mxGetPr(mxGetField(prhs[0],0,"J"));
	aA   =mxGetPr(mxGetField(prhs[0],0,"A"));
	 aRc = mxGetPr(mxGetField(prhs[0],0,"rc"));
	 aOmg =mxGetPr(mxGetField(prhs[0],0,"omg"));
	 aV   = mxGetPr(mxGetField(prhs[0],0,"v"));
 //   aGe  =mxGetPr(mxGetField(prhs[0],0,"ge"));
  //  aDq  =mxGetPr(mxGetField(prhs[0],0,"dq"));
//    aGpr =mxGetPr(mxGetField(prhs[0],0,"gpr"));
 //   aGdr =mxGetPr(mxGetField(prhs[0],0,"gdr"));

    rows=(int)aN[0];   
	sN1[0]=rows;
  
  
    aAlf = mxGetPr(mxGetField(prhs[0],0,"alf"));
    aBet = mxGetPr(mxGetField(prhs[0],0,"bet"));
 //   aGam = mxGetPr(mxGetField(prhs[0],0,"gam"));
//    aDel = mxGetPr(mxGetField(prhs[0],0,"del"));
 //   aOmg = mxGetPr(mxGetField(prhs[0],0,"omg"));
   
//	aAc = mxGetPr(mxGetField(prhs[0],0,"ac"));
//	aA0 = mxGetPr(mxGetField(prhs[0],0,"ac0"));
//	aBc = mxGetPr(mxGetField(prhs[0],0,"bc"));
//	aB0 = mxGetPr(mxGetField(prhs[0],0,"bc0"));	  

	nCP=(int)mxGetM(prhs[1]);
	if (nCP==1)
		nCP=(int)mxGetN(prhs[1]);	 
	sNCP1[0]=s3NCP[1]=nCP;
	COP=mxGetPr(prhs[1]);
	del=mxGetPr(prhs[2]);
	tmp=mxGetPr(prhs[3]);
	Kn=(double)(tmp[0]);
	tmp=mxGetPr(prhs[4]);
	Cn=(double)(tmp[0]);
	tmp=mxGetPr(prhs[5]);
	mu=(double)(tmp[0]);
    
	plhs[0] = mxCreateNumericArray(1, sN1,mxDOUBLE_CLASS,mxREAL); 
	plhs[1] = mxCreateNumericArray(2, s3NCP,mxDOUBLE_CLASS,mxREAL); 
	plhs[2] = mxCreateNumericArray(2, s3NCP,mxDOUBLE_CLASS,mxREAL); 
	plhs[3] = mxCreateNumericArray(2, sNCP1,mxDOUBLE_CLASS,mxREAL);


	tau	= mxGetPr(plhs[0]); 
    F	= mxGetPr(plhs[1]); 
	delP	= mxGetPr(plhs[2]);
	lam	= mxGetPr(plhs[3]);

    for (c_cp=0;c_cp<nCP;c_cp++, del+=3, delP+=3, F+=3, lam++)        
    {         
		cp=(int)COP[c_cp];
		nlnk=(int)aLnk[cp-1];
		matprod(aA+9*(nlnk-1),aCvec+3*(cp-1),3,3,1,gsv);
		vecadd(aRc+3*(nlnk-1),gsv,rc);   // pozicija tacke
		veccrossprod(aOmg+3*(nlnk-1),gsv,tmpV);
		vecadd(aV+3*(nlnk-1),tmpV, vc);	
		if (rc[2]<0)
		{  //contact
			vecmulscalar(vc,-1,delP);
			vecmulscalar(del,Kn,tmpV);
			vecaddcoef(tmpV, Cn,delP, F);
            if (F[2]<0)
            {
                F[0]=F[1]=F[2]=0;
                vecmulscalar(del, -Kn/Cn, delP);
            }
            else if (F[0]*F[0]+F[1]*F[1]>mu*mu*F[2]*F[2]) //slip
			{
				vecmulscalar(del, Kn/Cn, tmpV);
				vecsub(tmpV, vc,M);
				M[2]=0;
				*lam=1/Cn-sqrt( sqr(vecmod(M))/sqr(mu*F[2]));
				vecmulscalar(M,1.0/(1/Cn-*lam),tmpV);
				F[0]=tmpV[0]; F[1]=tmpV[1];
				delP[0]=*lam*F[0]-vc[0];
				delP[1]=*lam*F[1]-vc[1];
			}
		}
		else
		{  //no contact
			vecmulscalar(del, -Kn/Cn, delP);
		}
		vecadd(gsv, del, tmpV);
		veccrossprod(tmpV,F,M);
		for (i=0; i<rows; i++)
		{
			tmpBet=aBet+3*(nlnk-1)+3*rows*i;
			tmpAlf=aAlf+3*(nlnk-1)+3*rows*i;
			tau[i]+=tmpBet[0]*F[0]+tmpBet[1]*F[1]+tmpBet[2]*F[2];
			tau[i]+=tmpAlf[0]*M[0]+tmpAlf[1]*M[1]+tmpAlf[2]*M[2];
		}
		
	} 
}
    
