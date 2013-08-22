#include<stdio.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_statistics.h>
//Function Declarations
double RSS(double y1[],double y2[],int Npts);
void add_noise(double y[],double sigma,int N);
void filter_out(double a[],double b[],double y[],double u[],int N,int n_order);

int main()
{
	printf("Butterworth Nth order MCMC Approximation\n");

	//Random Generator Setup
	const gsl_rng_type *T;
	gsl_rng *r;
	r = gsl_rng_alloc(gsl_rng_mt19937);

	//Filter Design
	int const n_order = 1;

	//Time constraints
	double const fs = 400,fc = 20,f_norm = fc*2/fs;
	double const t0 = 0.0,t1 = 1.0,dt = 1/fs;
	int const Npts = (int)((t1-t0)/dt);

	//Coeffecient Arrays
	double a_mdl[] = {1.0,-0.7265};
	double b_mdl[] = {0.1367,0.1367};
	double a_curr[n_order+1];
	double b_curr[n_order+1];
	double a_cand[n_order+1];
	double b_cand[n_order+1];
	
	//Output Arrays
	double u[Npts],y_mdl[Npts],t[Npts],y_curr[Npts],y_cand[Npts];

	//Counters
	int ii,jj,kk;

	//MCMC Variables
	int const N_Iters = 100000000;
	int flg = 0;
	int burnin;
	int accepted = 0;
	double alpha = 1.0;
	double rss_cand,rss_curr;
	double a_ratio;
	FILE *fp1,*fp2;

	fp1 = fopen("data/test1.dat","w");
	fp2 = fopen("data/test2.dat","w");

	//Populate t,u and zero outputs
	for (ii = 0; ii < Npts; ii++)
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
		y_mdl[ii] = 0;
		y_cand[ii] = 0;
		y_curr[ii] = 0;
	}
	//Get output for known filter
	filter_out(a_mdl,b_mdl,y_mdl,u,Npts,n_order);
	//Add Noise to output
	//add_noise(y_mdl,sigma,Npts);
	
	//Generate first guess at parameters
	for (ii = 0; ii <(n_order+1); ii++)
	{
		a_curr[ii] = gsl_ran_flat(r,0,1);
		b_curr[ii] = gsl_ran_flat(r,0,1);
	}
	a_curr[0] = 1.0;
	//b_curr[0] = 0.136;
	b_curr[1] = 0.136;
	
	
	//Get output for current guess
	filter_out(a_curr,b_curr,y_curr,u,Npts,n_order);
	//Add noise to current output
	//add_noise(y_curr,sigma,Npts);
		
	//MCMC
	//Set the kk value to zero this is used to check burnin
	kk = 0;

	for (ii = 0; ii < N_Iters; ii++)
	{
		//Generate Random Coefficients
		for (jj = 0; jj <(n_order+1); jj++)
		{
			if(jj ==0) b_cand[jj] = alpha*gsl_ran_gaussian(r,1)+b_curr[jj];
			else
			{
				//printf("rand step = %E\n",alpha*gsl_ran_gaussian(r,1));
				a_cand[jj] = alpha*gsl_ran_gaussian(r,1)+a_curr[jj];
				b_cand[jj] = alpha*gsl_ran_gaussian(r,1)+b_curr[jj];
			}			
		}
		a_cand[0] = 1.0;
		//b_cand[0] = 0.136;
		b_cand[1] = 0.136;
		//printf("a_cand[1] = %E\n",a_cand[1]);
		
		//Generate Output from Candidate
		filter_out(a_cand,b_cand,y_cand,u,Npts,n_order);
		for (jj = 0; jj < Npts; jj++)
		{
			//fprintf(fp1,"%E	%E	%E\n",y_mdl[jj],y_cand[jj],y_curr[jj]);
		
		}
		//Add noise
		//add_noise(y_cand,sigma,Npts);

		//printf("a_cand[1] = %E\n",a_cand[1]);
		//printf("b_cand = %E   %E\n",b_cand[0],b_cand[1]);
		//Generate RSS for current parameters
		rss_curr = RSS(y_mdl,y_curr,Npts);
		//printf("rss_curr = %E\n",rss_curr);
		//Generate RSS for candidate parameters
		rss_cand = RSS(y_mdl,y_cand,Npts);
		//printf("rss_cand = %E\n\n\n",rss_cand);
		a_ratio = exp(-rss_cand+rss_curr);
		//printf("a_ratio = %E\n",a_ratio);
		if(gsl_rng_uniform(r)<a_ratio)
		{
			for (jj = 0; jj < n_order+1; jj++)
			{
				a_curr[jj] = a_cand[jj];
				b_curr[jj] = b_cand[jj];
			}
			filter_out(a_curr,b_curr,y_curr,u,Npts,n_order);
			accepted++;
			//printf("a_curr[1] = %E\n",a_curr[1]);
		}

		if(ii%10000==0 && ii!=0 && flg == 0)
		{	
			//printf("A_Ratio = %E\n",(double)(accepted)/(double)(kk));
			//printf("alpha = %E\n\n",alpha);
			if ((double)(accepted)/kk<0.3)
			{
				alpha = alpha/1.2;
				kk = 0;
				accepted = 0;
			}
			else if((double)(accepted)/kk>0.35)
			{
				alpha = alpha*1.2;
				kk = 0;
				accepted = 0;
			}
			else
			{
				burnin = ii;
				flg = 1;
			}
		}
		if(flg==1)
		{
			fprintf(fp1,"%E\n",a_curr[1]);
			fprintf(fp2,"%E\n",b_curr[0]);
		}
		kk++;

	}
	
	printf("alpha = %E\n",alpha);
	printf("accepted = %d\n",accepted);
	printf("flg = %d\n",flg);
	printf("burnin = %d\n",burnin);
	fclose(fp1);
	fclose(fp2);
	return 0;
}

void filter_out(double a[],double b[],double y[],double u[],int N,int n_order)
{
	int ii,jj;

	for(ii=0;ii<N;ii++)y[ii] = 0;
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
		if(abs(y[ii])>1000)y[ii] = 0;
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

double RSS(double y1[],double y2[],int Npts)
{
	int ii;
	double total = 0;
	for (ii = 0; ii < Npts; ii++)
	{
		total = total + pow((y1[ii]-y2[ii]),2);
	}
	return total;
}
