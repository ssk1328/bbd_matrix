#include "stdio.h"
#include "stdlib.h"
#include "float.h"
#include "omp.h"
#include "time.h"

void main()
{

long int MAX = 1000000;
double res[MAX];
int i;

clock_t begin, end;
begin = clock();

#pragma omp prallel for
for (i = 0; i<MAX; i++)
{
	res[i] = 2*i+i/4;
}

end = clock();
double time_spent = (double)(end - begin); // CLOCKS_PER_SEC;
printf("%f\n", time_spent);

}

