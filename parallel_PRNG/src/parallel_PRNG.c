#include<omp.h>
#include<stdlib.h>
#include <stdio.h>
int seed[10];
int main(int argc, char *argv[])
{
	int ii,s;
	for (ii = 0; ii < 10; ii++)
	{
		seed[ii] = rand();
	}

	#pragma omp parallel private(s)
	{
		s = seed[omp_get_thread_num()];
		#pragma omp for
		for (ii = 0; ii < 1000; ii++)
		{
			printf("%d %d %d\n",ii,omp_get_thread_num(),s);
			s = (s*17931+7391);
		}
		seed[omp_get_thread_num()] = s;
	}
		
	return 0;
}
