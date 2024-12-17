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
/* 3mm.c: this file is part of PolyBenchC */
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
/* Include polybench common header. */
#include "polybench.h"
#include <sys/time.h>
/* Include benchmark-specific header. */
#include "3mm.h"
/* Array initialization. */
static void init_array(int ni, int nj, int nk, int nl, int nm, double A[((800*4)+0)][((1000*4)+0)], double B[((1000*4)+0)][((900*4)+0)], double C[((900*4)+0)][((1200*4)+0)], double D[((1200*4)+0)][((1100*4)+0)])
{
	int i, j;
	#pragma cetus private(i, j) 
	#pragma loop name init_array#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((1L+(3L*ni))+((3L*ni)*nk)))) private(i, j)
	for (i=0; i<ni; i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name init_array#0#0 
		for (j=0; j<nk; j ++ )
		{
			A[i][j]=(((double)(((i*j)+1)%ni))/(5*ni));
		}
	}
	#pragma cetus private(i, j) 
	#pragma loop name init_array#1 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((1L+(3L*nk))+((3L*nj)*nk)))) private(i, j)
	for (i=0; i<nk; i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name init_array#1#0 
		for (j=0; j<nj; j ++ )
		{
			B[i][j]=(((double)(((i*(j+1))+2)%nj))/(5*nj));
		}
	}
	#pragma cetus private(i, j) 
	#pragma loop name init_array#2 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((1L+(3L*nj))+((3L*nj)*nm)))) private(i, j)
	for (i=0; i<nj; i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name init_array#2#0 
		for (j=0; j<nm; j ++ )
		{
			C[i][j]=(((double)((i*(j+3))%nl))/(5*nl));
		}
	}
	#pragma cetus private(i, j) 
	#pragma loop name init_array#3 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<((1L+(3L*nm))+((3L*nl)*nm)))) private(i, j)
	for (i=0; i<nm; i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name init_array#3#0 
		for (j=0; j<nl; j ++ )
		{
			D[i][j]=(((double)(((i*(j+2))+2)%nk))/(5*nk));
		}
	}
	return ;
}

/*
DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output.
*/
static void print_array(int ni, int nl, double G[((800*4)+0)][((1100*4)+0)])
{
	int i, j;
	fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
	fprintf(stderr, "begin dump: %s", "G");
	#pragma cetus private(i, j) 
	#pragma loop name print_array#0 
	for (i=0; i<ni; i ++ )
	{
		#pragma cetus private(j) 
		#pragma loop name print_array#0#0 
		for (j=0; j<nl; j ++ )
		{
			if ((((i*ni)+j)%20)==0)
			{
				fprintf(stderr, "\n");
			}
			fprintf(stderr, "%0.2lf ", G[i][j]);
		}
	}
	fprintf(stderr, "\nend   dump: %s\n", "G");
	fprintf(stderr, "==END   DUMP_ARRAYS==\n");
	return ;
}

/*
Main computational kernel. The whole function will be timed,
   including the call and return.
*/
static void kernel_3mm(int ni, int nj, int nk, int nl, int nm, double E[((800*4)+0)][((900*4)+0)], double A[((800*4)+0)][((1000*4)+0)], double B[((1000*4)+0)][((900*4)+0)], double F[((900*4)+0)][((1100*4)+0)], double C[((900*4)+0)][((1200*4)+0)], double D[((1200*4)+0)][((1100*4)+0)], double G[((800*4)+0)][((1100*4)+0)])
{
	int i, j, k;
	/* #pragma scop */
	/* E := AB */
	#pragma cetus firstprivate(E) 
	#pragma cetus private(i, j, k) 
	#pragma cetus lastprivate(E) 
	#pragma loop name kernel_3mm#0 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(((1L+(3L*ni))+((4L*ni)*nj))+(((3L*ni)*nj)*nk)))) private(i, j, k) firstprivate(E) lastprivate(E)
	for (i=0; i<ni; i ++ )
	{
		#pragma cetus firstprivate(E) 
		#pragma cetus private(j, k) 
		#pragma cetus lastprivate(E) 
		#pragma loop name kernel_3mm#0#0 
		for (j=0; j<nj; j ++ )
		{
			E[i][j]=0.0;
			#pragma cetus private(k) 
			#pragma loop name kernel_3mm#0#0#0 
			/* #pragma cetus reduction(+: E[i][j])  */
			for (k=0; k<nk;  ++ k)
			{
				E[i][j]+=(A[i][k]*B[k][j]);
			}
		}
	}
	/* F := CD */
	#pragma cetus firstprivate(F) 
	#pragma cetus private(i, j, k) 
	#pragma cetus lastprivate(F) 
	#pragma loop name kernel_3mm#1 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(((1L+(3L*nj))+((4L*nj)*nl))+(((3L*nj)*nl)*nm)))) private(i, j, k) firstprivate(F) lastprivate(F)
	for (i=0; i<nj; i ++ )
	{
		#pragma cetus firstprivate(F) 
		#pragma cetus private(j, k) 
		#pragma cetus lastprivate(F) 
		#pragma loop name kernel_3mm#1#0 
		for (j=0; j<nl; j ++ )
		{
			F[i][j]=0.0;
			#pragma cetus private(k) 
			#pragma loop name kernel_3mm#1#0#0 
			/* #pragma cetus reduction(+: F[i][j])  */
			for (k=0; k<nm;  ++ k)
			{
				F[i][j]+=(C[i][k]*D[k][j]);
			}
		}
	}
	/* G := EF */
	#pragma cetus firstprivate(G) 
	#pragma cetus private(i, j, k) 
	#pragma cetus lastprivate(G) 
	#pragma loop name kernel_3mm#2 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(((1L+(3L*ni))+((4L*ni)*nl))+(((3L*ni)*nj)*nl)))) private(i, j, k) firstprivate(G) lastprivate(G)
	for (i=0; i<ni; i ++ )
	{
		#pragma cetus firstprivate(G) 
		#pragma cetus private(j, k) 
		#pragma cetus lastprivate(G) 
		#pragma loop name kernel_3mm#2#0 
		for (j=0; j<nl; j ++ )
		{
			G[i][j]=0.0;
			#pragma cetus private(k) 
			#pragma loop name kernel_3mm#2#0#0 
			/* #pragma cetus reduction(+: G[i][j])  */
			for (k=0; k<nj;  ++ k)
			{
				G[i][j]+=(E[i][k]*F[k][j]);
			}
		}
	}
	/*
	#pragma endscop
	
	*/
	return ;
}

float time_diff(struct timeval * start, struct timeval * end)
{
	float _ret_val_0;
	_ret_val_0=((end->tv_sec-start->tv_sec)+(1.0E-6*(end->tv_usec-start->tv_usec)));
	return _ret_val_0;
}

int main(int argc, char * * argv)
{
	/* Retrieve problem size. */
	int ni = 800*4;
	int nj = 900*4;
	int nk = 1000*4;
	int nl = 1100*4;
	int nm = 1200*4;
	/* Variable declarationallocation. */
	double (* E)[((800*4)+0)][((900*4)+0)];
	double (* A)[((800*4)+0)][((1000*4)+0)];
	double (* B)[((1000*4)+0)][((900*4)+0)];
	double (* F)[((900*4)+0)][((1100*4)+0)];
	double (* C)[((900*4)+0)][((1200*4)+0)];
	double (* D)[((1200*4)+0)][((1100*4)+0)];
	double (* G)[((800*4)+0)][((1100*4)+0)];
	struct timeval start, end;
	int _ret_val_0;
	E=((double (* )[((800*4)+0)][((900*4)+0)])polybench_alloc_data(((800*4)+0)*((900*4)+0), sizeof (double)));
	;
	A=((double (* )[((800*4)+0)][((1000*4)+0)])polybench_alloc_data(((800*4)+0)*((1000*4)+0), sizeof (double)));
	;
	B=((double (* )[((1000*4)+0)][((900*4)+0)])polybench_alloc_data(((1000*4)+0)*((900*4)+0), sizeof (double)));
	;
	F=((double (* )[((900*4)+0)][((1100*4)+0)])polybench_alloc_data(((900*4)+0)*((1100*4)+0), sizeof (double)));
	;
	C=((double (* )[((900*4)+0)][((1200*4)+0)])polybench_alloc_data(((900*4)+0)*((1200*4)+0), sizeof (double)));
	;
	D=((double (* )[((1200*4)+0)][((1100*4)+0)])polybench_alloc_data(((1200*4)+0)*((1100*4)+0), sizeof (double)));
	;
	G=((double (* )[((800*4)+0)][((1100*4)+0)])polybench_alloc_data(((800*4)+0)*((1100*4)+0), sizeof (double)));
	;
	/* Initialize array(s). */
	init_array(ni, nj, nk, nl, nm,  * A,  * B,  * C,  * D);
	/* Start timer. */
	;
	/* Run kernel. */
	gettimeofday( & start, (void * )0);
	kernel_3mm(ni, nj, nk, nl, nm,  * E,  * A,  * B,  * F,  * C,  * D,  * G);
	gettimeofday( & end, (void * )0);
	printf("Time in seconds %0.8f \n", time_diff( & start,  & end));
	/* Stop and print timer. */
	;
	;
	/*
	Prevent dead-code elimination. All live-out data must be printed
	     by the function call in argument.
	*/
	if ((argc>42)&&( ! strcmp(argv[0], "")))
	{
		print_array(ni, nl,  * G);
	}
	/* Be clean. */
	free((void * )E);
	;
	free((void * )A);
	;
	free((void * )B);
	;
	free((void * )F);
	;
	free((void * )C);
	;
	free((void * )D);
	;
	free((void * )G);
	;
	_ret_val_0=0;
	return _ret_val_0;
}
