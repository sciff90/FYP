#include<stdio.h>
#include<math.h>
int main()
{
	printf("First Order Filter MCMC Approximation");
	
	//Filter Design
	int const P = 4,Q = 4,N = 1000;

	//Time constraints
	double const t1 = 0.0,t2 =2,dt = (t2-t1)/N;

	//Output arrays
	double u[N],y[N],t[N];

	//Coefficient Arrays
	;

	//Counters
	int ii,jj,kk;

	double a[] = {1.0,-0.5772,0.4218,-0.0563};
	double b[] = {0.09840,0.2956,0.2956,0.985};

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
		y[ii] = 0;
	}

	//populate output array
	FILE *fp;
	fp = fopen("data/filter.dat","w");

	for (ii = 0; ii < (N-Q); ii++)
	{
		for (jj = 1; jj <= Q; jj++)
		{
			y[ii] = y[ii] - a[jj]*y[ii-jj];
		}
		for (jj = 0; jj <= P; jj++)
		{
			y[ii] = y[ii]+b[jj]*u[ii-jj];
		}
		y[ii] = y[ii]/a[0];
		fprintf(fp,"%E	%E	%E\n",t[ii],u[ii],y[ii]);
	}

	fclose(fp);
	return 0;
}
