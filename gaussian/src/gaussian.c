#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
int main(int argc, char *argv[])
{
	const gsl_rng_type * T;
	gsl_rng * r;

	//select random number generator
	r = gsl_rng_alloc(gsl_rng_mt19937);

	FILE *fp1,*fp2;
	

	const int N = 10000;
	double u[N],v[N];
	double sigma = 1;
	int ii;

	fp1 = fopen("./data/uniform.dat","w");
	fp2 = fopen("./data/normal.dat","w");


	for (ii = 0; ii < N; ii++)
	{
		u[ii] = gsl_rng_uniform(r);
		fprintf(fp1,"%f\n",u[ii]);
		v[ii] = gsl_ran_gaussian(r,sigma);
		fprintf(fp2,"%f\n",v[ii]);
	}

	gsl_rng_free(r);
	fclose(fp1);
	fclose(fp2);
	system("gnuplot gaussian.plt");

	return 0;
}
