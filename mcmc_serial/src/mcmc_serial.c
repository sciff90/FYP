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
	
	int accepted = 0;

	int ii;

	double x_0,alpha,candidate,a,a_rate;

	fp = fopen("./data/mcmc_serial.dat","w");
	//generate random starting point
	x_0 = gsl_ran_flat(r,-20,20);
	alpha = 1.0;

	for (ii = 0; ii < Iters; ii++)
	{
		candidate = x_0+gsl_ran_gaussian(r,alpha);
		a = ((gsl_ran_gaussian_pdf(candidate,1)+gsl_ran_gaussian_pdf(candidate-10,2)))/((gsl_ran_gaussian_pdf(x_0,1)+gsl_ran_gaussian_pdf(x_0-10,2)));
		//printf("a = %f\n",a);
		if(gsl_rng_uniform(r)<a)
		{
			x_0 = candidate;
			accepted++;
		}
		fprintf(fp,"%f\n",x_0);
		if(ii%100==0&&ii!=0)
		{
			a_rate = accepted/(double)(ii);
			printf("acceptance rate = %f\n",a_rate);
			printf("alpha = %f\n",alpha);
			if(a_rate>0.45)
			{
				alpha = alpha*1.5;
			}
			else if(a_rate<0.25)
			{
				alpha = alpha/1.5;
			}
			

		}

	}
	fclose(fp);
	gsl_rng_free(r);
	system("gnuplot mcmc_serial.plt");
	printf("Number accpeted = %d",accepted);
	

	return 0;
}
