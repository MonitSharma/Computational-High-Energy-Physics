
#include <math.h>
#define NRANSI
#include "../include/nrutilD.h"
#define EPS 1.0e-4

void fdjac(int n, double x[], double fvec[], double **df,
	void (*vecfunc)(int, double [], double []))
{
	int i,j;
	double h,temp,*f;

	f=vector(1,n);
	for (j=1;j<=n;j++) {
		temp=x[j];
		h=EPS*fabs(temp);
		if (h == 0.0) h=EPS;
		x[j]=temp+h;
		h=x[j]-temp;
		(*vecfunc)(n,x,f);
		x[j]=temp;
		for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h;
	}
	free_vector(f,1,n);
}
#undef EPS
#undef NRANSI
