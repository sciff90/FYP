#include<stdio.h>
#include<time.h>

int main(int argc, char *argv[])
{
	clock_t start = clock(),diff;
	int ii,msec;

	for (ii = 0; ii <100000; ii++)
	{
		printf("Testing\n");
	}

	diff = clock()-start;
	msec = diff*1000/CLOCKS_PER_SEC;
	printf("Time taken %d seconds %d milliseconds\n",msec/1000,msec%1000);


	return 0;
}
