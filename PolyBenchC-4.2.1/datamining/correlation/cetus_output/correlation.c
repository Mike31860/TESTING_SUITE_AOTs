/*
Copyright (C) 1991-2012 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it andor
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http:www.gnu.org/licenses/>. 
*/
/*
This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it. 
*/
/* We do support the IEC 559 math functionality, real and complex.  */
/*
wchar_t uses ISOIEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0. 
*/
/* We do not support C11 <threads.h>.  */
/*

 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http:polybench.sourceforge.net

*/
/* correlation.c: this file is part of PolyBenchC */
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
/* Include polybench common header. */
#include <polybench.h>
/* Include benchmark-specific header. */
#include "correlation.h"
/* Array initialization. */
static void init_array(int m, int n, double * float_n, double data[((1400*3)+0)][((1200*3)+0)])
{
	int i, j;
	( * float_n)=(((double)1400)*3);
	#pragma cetus private(i, j) 
	#pragma loop name init_array#0 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, j)
	for (i=0; i<(1400*3); i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name init_array#0#0 
		for (j=0; j<(1200*3); j ++ )
		{
			data[i][j]=(((((double)(i*j))/1200)*3)+i);
		}
	}
	return ;
}

/*
DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output.
*/
static void print_array(int m, double corr[((1200*3)+0)][((1200*3)+0)])
{
	int i, j;
	fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
	fprintf(stderr, "begin dump: %s", "corr");
	#pragma cetus private(i, j) 
	#pragma loop name print_array#0 
	for (i=0; i<m; i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name print_array#0#0 
		for (j=0; j<m; j ++ )
		{
			if ((((i*m)+j)%20)==0)
			{
				fprintf(stderr, "\n");
			}
			fprintf(stderr, "%0.2lf ", corr[i][j]);
		}
	}
	fprintf(stderr, "\nend   dump: %s\n", "corr");
	fprintf(stderr, "==END   DUMP_ARRAYS==\n");
	return ;
}

/*
Main computational kernel. The whole function will be timed,
   including the call and return.
*/
static void kernel_correlation(int m, int n, double float_n, double data[((1400*3)+0)][((1200*3)+0)], double corr[((1200*3)+0)][((1200*3)+0)], double mean[((1200*3)+0)], double stddev[((1200*3)+0)])
{
	int i, j, k;
	double eps = 0.1;
	#pragma scop 
	#pragma cetus firstprivate(mean) 
	#pragma cetus private(i, j) 
	#pragma cetus lastprivate(mean) 
	#pragma loop name kernel_correlation#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((1L+(5L*m))+((3L*m)*n)))) private(i, j) firstprivate(mean) lastprivate(mean)
	for (j=0; j<m; j ++ )
	{
		mean[j]=0.0;
		#pragma cetus private(i) 
		#pragma loop name kernel_correlation#0#0 
		/* #pragma cetus reduction(+: mean[j])  */
		for (i=0; i<n; i ++ )
		{
			mean[j]+=data[i][j];
		}
		mean[j]/=float_n;
	}
	#pragma cetus firstprivate(stddev) 
	#pragma cetus private(i, j) 
	#pragma cetus lastprivate(stddev) 
	#pragma loop name kernel_correlation#1 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((1L+(107L*m))+((3L*m)*n)))) private(i, j) firstprivate(stddev) lastprivate(stddev)
	for (j=0; j<m; j ++ )
	{
		stddev[j]=0.0;
		#pragma cetus private(i) 
		#pragma loop name kernel_correlation#1#0 
		/* #pragma cetus reduction(+: stddev[j])  */
		for (i=0; i<n; i ++ )
		{
			stddev[j]+=((data[i][j]-mean[j])*(data[i][j]-mean[j]));
		}
		stddev[j]/=float_n;
		stddev[j]=sqrt(stddev[j]);
		/*
		The following in an inelegant but usual way to handle
		         near-zero std. dev. values, which below would cause a zero-
		         divide.
		*/
		stddev[j]=((stddev[j]<=eps) ? 1.0 : stddev[j]);
	}
	/* Center and reduce the column vectors. */
	#pragma cetus private(i, j) 
	#pragma loop name kernel_correlation#2 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((1L+(3L*n))+((104L*m)*n)))) private(i, j)
	for (i=0; i<n; i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name kernel_correlation#2#0 
		for (j=0; j<m; j ++ )
		{
			data[i][j]-=mean[j];
			data[i][j]/=(sqrt(float_n)*stddev[j]);
		}
	}
	/* Calculate the m m correlation matrix. */
	#pragma cetus firstprivate(corr) 
	#pragma cetus private(i, j, k) 
	#pragma cetus lastprivate(corr) 
	#pragma loop name kernel_correlation#3 
	#pragma cetus parallel 
	#pragma omp parallel for private(i, j, k) firstprivate(corr) lastprivate(corr)
	for (i=0; i<(m-1); i ++ )
	{
		corr[i][i]=1.0;
		#pragma cetus firstprivate(corr) 
		#pragma cetus private(j, k) 
		#pragma cetus lastprivate(corr) 
		#pragma loop name kernel_correlation#3#0 
		for (j=(i+1); j<m; j ++ )
		{
			corr[i][j]=0.0;
			#pragma cetus private(k) 
			#pragma loop name kernel_correlation#3#0#0 
			/* #pragma cetus reduction(+: corr[i][j])  */
			for (k=0; k<n; k ++ )
			{
				corr[i][j]+=(data[k][i]*data[k][j]);
			}
			corr[j][i]=corr[i][j];
		}
	}
	corr[m-1][m-1]=1.0;
	#pragma endscop 
	return ;
}

int main(int argc, char * * argv)
{
	/* Retrieve problem size. */
	int n = 1400*3;
	int m = 1200*3;
	/* Variable declarationallocation. */
	double float_n;
	double (* data)[((1400*3)+0)][((1200*3)+0)];
	double (* corr)[((1200*3)+0)][((1200*3)+0)];
	double (* mean)[((1200*3)+0)];
	double (* stddev)[((1200*3)+0)];
	int _ret_val_0;
	data=((double (* )[((1400*3)+0)][((1200*3)+0)])polybench_alloc_data(((1400*3)+0)*((1200*3)+0), sizeof (double)));
	;
	corr=((double (* )[((1200*3)+0)][((1200*3)+0)])polybench_alloc_data(((1200*3)+0)*((1200*3)+0), sizeof (double)));
	;
	mean=((double (* )[((1200*3)+0)])polybench_alloc_data((1200*3)+0, sizeof (double)));
	;
	stddev=((double (* )[((1200*3)+0)])polybench_alloc_data((1200*3)+0, sizeof (double)));
	;
	/* Initialize array(s). */
	init_array(m, n,  & float_n,  * data);
	/* Start timer. */
	;
	/* Run kernel. */
        struct timeval start, end; 
        gettimeofday(&start, NULL);
	kernel_correlation(m, n, float_n,  * data,  * corr,  * mean,  * stddev);
        gettimeofday(&end, NULL);
        long seconds = (end.tv_sec - start.tv_sec);
        long micros = ((seconds * 1000000) + end.tv_usec) - (start.tv_usec);
        printf("The elapsed time is %ld seconds and %ld micros\n", seconds, micros);
	/* Stop and print timer. */
	;
	;
	/*
	Prevent dead-code elimination. All live-out data must be printed
	     by the function call in argument.
	*/
	if ((argc>42)&&( ! strcmp(argv[0], "")))
	{
		print_array(m,  * corr);
	}
	/* Be clean. */
	free((void * )data);
	;
	free((void * )corr);
	;
	free((void * )mean);
	;
	free((void * )stddev);
	;
	_ret_val_0=0;
	return _ret_val_0;
}
