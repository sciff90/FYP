#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
int main(int argc, char *argv[])
{
	const gsl_rng_type * T;
	gsl_rng * r;

	//select random number generator
	r = gsl_rng_alloc(gsl_rng_mt19937);
	FILE *fp;

	const int N = 1000,Iters = 100000;
	int ii;

	double x_0,alpha,candidate,a;

	fp = fopen("./data/test.dat","w");
	//generate random starting point
	x_0 = gsl_ran_flat(r,-20,20);
	alpha = 0.2;

	for (ii = 0; ii < Iters; ii++)
	{
		candidate = x_0+gsl_ran_gaussian(r,alpha);
		a = ((gsl_ran_gaussian_pdf(candidate,1)+gsl_ran_gaussian_pdf(candidate-10,2)))/((gsl_ran_gaussian_pdf(x_0,1)+gsl_ran_gaussian_pdf(x_0-10,2)));
		//printf("a = %f\n",a);
		if(gsl_rng_uniform(r)<a)
			x_0 = candidate;
		fprintf(fp,"%f\n",x_0);

	}
	fclose(fp);
	gsl_rng_free(r);
	system("gnuplot mcmc_serial.plt");
	

	return 0;
}
