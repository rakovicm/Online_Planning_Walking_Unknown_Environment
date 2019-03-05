#include "math.h"

void matprod(double* ar1, double* ar2,int n1, int m1,  int m2,double *res)
{
    int brn1,brm1,brm2;
    for(brn1=0;brn1<n1;brn1++)
        for(brm2=0;brm2<m2;brm2++)
        {
            res[brn1+n1*brm2]=0;
            for (brm1=0;brm1<m1;brm1++)         
                res[brn1+n1*brm2]=res[brn1+n1*brm2]+ar1[brn1+n1*brm1]*ar2[brm1+m1*brm2];    
        }
   
}

void matmulscalar(double *ar1, double c, int n, int m, double* res)
{
	int i, j;
	for (i=0; i<n; i++)
		for (j=0; j<m; j++)
			*(res++)=c* (*(ar1++));
}


double vecmod(double* vec)
{
	return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}

void vecmulscalar(double* ar1,double c,  double* res)
{
    res[0]=c*ar1[0];
    res[1]=c*ar1[1];
    res[2]=c*ar1[2];
}

void vecadd(double* ar1, double* ar2, double* res)
{
    res[0]=ar1[0]+ar2[0];
    res[1]=ar1[1]+ar2[1];
    res[2]=ar1[2]+ar2[2];
}

void vecaddcoef(double* ar1,double c, double* ar2, double* res)
{
    res[0]=ar1[0]+c*ar2[0];
    res[1]=ar1[1]+c*ar2[1];
    res[2]=ar1[2]+c*ar2[2];
}

void vecsub(double* ar1, double* ar2, double* res)
{
    res[0]=ar1[0]-ar2[0];
    res[1]=ar1[1]-ar2[1];
    res[2]=ar1[2]-ar2[2];
}
void matcopy(double* dest, double* src,int n)
{
    int i;
    for(i=0; i<n;i++) 
         dest[i]=src[i];
    
}

void veccrossprod(double* ar1,double* ar2,double *result)
{
    result[0]=ar1[1]*ar2[2]-ar1[2]*ar2[1];
    result[1]=ar1[2]*ar2[0]-ar1[0]*ar2[2];
    result[2]=ar1[0]*ar2[1]-ar1[1]*ar2[0];
}

double vecdotprod(double* ar1,double* ar2)
{
    double temp;
    temp=ar1[1]*ar2[1];
    temp=temp+ar1[2]*ar2[2];
    temp=temp+ar1[0]*ar2[0];
    return temp;
}

void tran3x3(double* ar1,double *result)
{
    result[0]=ar1[0];
    result[1]=ar1[3];
    result[2]=ar1[6];
    result[3]=ar1[1];
    result[4]=ar1[4];
    result[5]=ar1[7];
    result[6]=ar1[2];
    result[7]=ar1[5];
    result[8]=ar1[8];
}

void zeros(double *res, int n, int m)
{
    int i;
    for(i=0;i<n*m;i++)
        res[i]=0;
    
}
void eye(double *res, int n)
{
    int i;
    for(i=0;i<n;i++)
        res[i+3*i]=1;
    
}

void RotX(double* res,double q)
{ 
    res[0]=1;
    res[1]=res[2]=res[3]=res[6]=0;
    res[4]=res[8]=cos(q);
    res[5]=sin(q);
    res[7]=-res[5];
 
    
}

void RotY(double* res,double q)
{ 
    res[4]=1;
    res[1]=res[3]=res[5]=res[7]=0;
    res[0]=res[8]=cos(q);
    res[6]=sin(q);
    res[2]=-res[6];
}

void RotZ(double* res,double q)
{ 
    res[8]=1;
    res[2]=res[5]=res[6]=res[7]=0;
    res[0]=res[4]=cos(q);
    res[3]=-sin(q);
    res[1]=-res[3];
}


void CrossMat(double* vec, double* mat)
{	
	mat[0]=mat[4]=mat[8]=0;
	mat[1]=vec[2];
	mat[2]=-vec[1];
	mat[3]=-vec[2];
	mat[5]=vec[0];
	mat[6]=vec[1];
	mat[7]=-vec[0];
}

void TRotXYZ(double *mat,double *res)
{
	res[0]=atan(-mat[7]/mat[8]);
	res[1]=atan(mat[6]/(mat[8]*cos(res[0])- mat[7]*sin(res[0])));
	res[2]=atan(-mat[3]/mat[0]);
}

void TRotXZY(double *mat,double *res)
{
	res[0]=atan(mat[5]/mat[4]);
	res[1]=atan(-mat[3]/(mat[4]*cos(res[0])+mat[5]*sin(res[0])));
	res[2]=atan(mat[6]/mat[0]);
}

void TRotYXZ(double *T,double *res)
{
	res[0]=atan(T[6]/T[8]);
	res[1]=atan(-T[7]/(T[8]*cos(res[0]) +  T[6]*sin(res[0])));
	res[2]=atan(T[1]/T[4]);
}

void TRotYZX(double *T,double *res)
{
	res[0]=atan(-T[2]/T[0]);
	res[1]= atan(T[1]/(T[0]*cos(res[0]) - T[2]*sin(res[0])));
	res[2]=atan(-T[7]/T[4]);
}

void TRotZXY(double *T,double *res)
{
	res[0]=atan(-T[3]/T[4]);
	res[1]= atan(T[5]/(T[4]*cos(res[0]) - T[3]*sin(res[0])));
	res[2]=atan(-T[2]/T[8]);
}

void TRotZYX(double *T,double *res)
{
	res[0]=atan(T[1]/T[0]);
	res[1]=atan(-T[2]/(T[0]*cos(res[0]) + T[1]*sin(res[0])));
	res[2]=atan(T[5]/T[8]);
}

