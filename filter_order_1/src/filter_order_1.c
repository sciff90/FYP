#include<stdio.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_statistics.h>
//Function Declarations
void filter_out(double a[],double b[], double y[],double u[],int N,int n_order);
void add_noise(double y[],double sigma,int N);
double var(double y1[],double y2[],int N);
int main()
{
	printf("Butter Worth nth-Order Filter MCMC Approximation");
	
	//Random Generator Setup
	const gsl_rng_type * T;
	gsl_rng * r;
	r = gsl_rng_alloc(gsl_rng_mt19937);
	
	//Filter Design
	int const n_order = 1;

	//Time constraints
	double const fs = 400,fc = 20,f_norm = fc*2/fs;

	double const t0 = 0.0,t1 =1.0,dt = 1/fs;

	int const N = (int)((t1-t0)/dt);

	//Output arrays
	double u[N],y_proposed[N],y_out[N],t[N],y_candidate[N];

	//Coefficient Arrays
	double a[] = {1.0,-0.7265};	
	double b[] = {0.1367,0.1367};
	double a_proposed[n_order+1];
	double b_proposed[n_order+1];
	double a_candidate[n_order+1];
	double b_candidate[n_order+1];

	//Counters
	int ii,jj,kk;


	//populate time array and unit step array
	for (ii = 0; ii < N; ii++)
	{
		if (ii==0)
		{
			u[ii] = 0;
			t[ii] = 0;
		}
		else 
		{
			u[ii] = 1;
			t[ii] = ii*dt;
		}
		y_out[ii] = 0;
		y_proposed[ii] = 0;
		y_candidate[ii] = 0;
	}

	//Get output for known filter
	filter_out(a,b,y_out,u,N,n_order);
	//Add noise to output
	add_noise(y_out,0.05,N);

	//Generate initial proposed values
	for (ii = 0; ii <=n_order; ii++)
	{
		a_proposed[ii] = gsl_ran_flat(r,-1,1);
		b_proposed[ii] = gsl_ran_flat(r,-1,1);
	}

	a_proposed[0] = 1.0;
	a_candidate[0] = 1.0;

	//Begin MCMC Algorithm
	double accept;
	double alpha = 1.0;
	//Generate Output for proposal
	filter_out(a_proposed,b_proposed,y_proposed,u,N,n_order);
	add_noise(y_proposed,0.05,N);
	printf("Testing\n");
	for (ii = 0; ii < 100; ii++)
	{
		//Generate Random Coefficients
		for (jj = 0; jj <=(n_order+1); jj++)
		{
			if(jj==0) b_candidate[jj] += gsl_ran_gaussian(r,alpha);
			else
			{
				a_candidate[jj] += gsl_ran_gaussian(r,alpha);
				b_candidate[jj] += gsl_ran_gaussian(r,alpha);
			}

		}
		//Generate Candidate
		filter_out(a_candidate,b_candidate,y_candidate,u,N,n_order);
		add_noise(y_candidate,0.05,N);
		//Calculate acceptance
		accept = var(y_out,y_candidate,N)/var(y_out,y_proposed,N);
		//Check a value
		if(gsl_rng_uniform(r)<accept)
		{
			for (jj = 0; jj <= n_order+1; jj++)
			{
				if(jj=0) b_proposed[jj] = b_candidate[jj];
				else
				{	
					a_proposed[jj] = a_candidate[jj];
					b_proposed[jj] = b_candidate[jj];		
				}
			}
			filter_out(a_proposed,b_proposed,y_proposed,u,N,n_order);
			add_noise(y_proposed,0.05,N);
			printf("accepted\n");
		}
		printf("Not accepted\n");

		
	}
	//Output To Files
	FILE *fp,*fp1;
	fp = fopen("data/filter.dat","w");
	fp1 = fopen("data/difference.dat","w");
		
	for (ii = 0; ii < N; ii++)
	{
		fprintf(fp,"%E	%E	%E	%E\n",t[ii],u[ii],y_out[ii],y_proposed[ii]);	
		
	}

	fclose(fp);
	return 0;
}

void filter_out(double a[],double b[],double y[],double u[],int N,int n_order)
{
	int ii,jj;

	for (ii = n_order; ii < (N); ii++)
	{
		for (jj = 1; jj <= n_order; jj++)
		{
			y[ii] = y[ii] - a[jj]*y[ii-jj];
		}
		for (jj = 0; jj <= n_order; jj++)
		{
			y[ii] = y[ii]+b[jj]*u[ii-jj];
		}
		
		y[ii] = y[ii]/a[0];
	}
}
void add_noise(double y[],double sigma,int N)
{
	double x_rand;
	int ii;
	const gsl_rng_type * T;
	gsl_rng * r;
	r = gsl_rng_alloc(gsl_rng_mt19937);

	for (ii = 0; ii < N; ii++)
	{
		x_rand = gsl_ran_gaussian(r,sigma);
		y[ii] += x_rand;
	}	

}
double var(double y1[],double y2[],int N)
{
	int ii;
	double total = 0;

	for (ii = 0; ii < N; ii++)
	{
		total += pow((y1[ii]-y2[ii]),2);
	}

	return total;
}
