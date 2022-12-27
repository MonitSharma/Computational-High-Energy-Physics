
#define NRANSI
#include "../include/nrutilD.h"

extern int nn;
extern double *fvec;
extern void (*nrfuncv)(int n, double v[], double f[]);

double fminimum(double x[])
{
	int i;
	double sum;

	(*nrfuncv)(nn,x,fvec);
	for (sum=0.0,i=1;i<=nn;i++) sum += SQR(fvec[i]);
	return 0.5*sum;
}
#undef NRANSI
