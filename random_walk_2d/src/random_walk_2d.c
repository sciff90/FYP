#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int main(int argc, char *argv[])
{
	//Make stack size unlimited
	system("ulimit -s unlimited");

	int ii;
	const long N = 1000;
	double x[N],y[N],dx ,dy;

	const gsl_rng_type * T;
	gsl_rng *r;
	FILE *fp1;

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	x[0] = 0;
	y[0] = 0;
	
	fp1 = fopen("./data/random_walk.dat","w");
	
	for (ii = 1; ii < N; ii++)
	{
		gsl_ran_dir_2d(r,&dx,&dy);	
		x[ii] = x[ii-1]+dx;
		y[ii] = y[ii-1]+dy;
		//printf("%f\n",dx);
		//printf("%f\n",dy);
		fprintf(fp1,"%f %f\n",x[ii],y[ii]);		
	}

	gsl_rng_free(r);
	fclose(fp1);

	//Call gnuplot
	system("gnuplot random_walk.plt");
	return 0;
}
