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
/* polybench.c: this file is part of PolyBenchC */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sched.h>
#include <math.h>
/*

 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http:polybench.sourceforge.net

*/
/*

 polybench.h: this file is part of PolyBench/C
 *
 * Polybench header for instrumentation.
 *
 * Programs must be compiled with `-I utilities utilities/polybench.c'
 *
 * Optionally, one can define:
 *
 * -DPOLYBENCH_TIME, to report the execution time,
 *   OR (exclusive):
 * -DPOLYBENCH_PAPI, to use PAPI H/W counters (defined in polybench.c)
 *
 *
 * See README or utilities/polybench.c for additional options.


*/
/*
Copyright (C) 1991-2007, 2009-2011, 2012 Free Software Foundation, Inc.
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

	ISO C99 Standard: 7.20 General utilities	<stdlib.h>

*/
/* Array padding. By default, none is used. */
/* default: */
/* Inter-array padding, for use with . By default, none is used. */
/* default: */
/* C99 arrays in function prototype. By default, do not use. */
/* default: */
/* Scalar loop bounds in SCoPs. By default, use parametric loop bounds. */
/* default: */
/*
Use the 'restrict' keyword to declare that the different arrays do not
 alias. By default, we do not use it as it is only supported in C99 and
 * even here several compilers do not properly get it.

*/
/* default: */
/*
Macros to reference an array. Generic for heap and stack arrays
   (C99).  Each array dimensionality has his own macro, to be used at
   declaration or as a function argument.
   Example:
   int b[x] => POLYBENCH_1D_ARRAY(b, x)
   int A[N][N] => POLYBENCH_2D_ARRAY(A, N, N)

*/
/* Macros for using arrays in the function prototypes. */
/* Macros for using arrays within the functions. */
/*
Macros to allocate heap arrays.
   Example:
   polybench_alloc_2d_array(N, M, double) => allocates N x M x sizeof(double)
					  and returns a pointer to the 2d array

*/
/* Macros for array declaration. */
/* Dead-code elimination macros. Use argcargv for the run-time check. */
/* Performance-related instrumentation. See polybench.c */
/* PAPI support. */
/* Timing support. */
/* PAPI support. */
/* Function prototypes. */
extern void *polybench_alloc_data(unsigned long long int n, int elt_size);
extern void polybench_free_data(void * ptr);
/* PolyBench internal functions that should not be directly called by */
/* the user, unless when designing customized execution profiling */
/* approaches. */
extern void polybench_flush_cache();
extern void polybench_prepare_instruments();
/* By default, collect PAPI counters on thread 0. */
/* Total LLC cache size. By default 32+MB.. */
int polybench_papi_counters_threadid = 0;
double polybench_program_total_flops = 0;
/*

 Allocation table, to enable inter-array padding. All data allocated
 * with polybench_alloc_data should be freed with polybench_free_data.


*/
struct polybench_data_ptrs
{
	void * * user_view;
	void * * real_ptr;
	int nb_entries;
	int nb_avail_entries;
};

static struct polybench_data_ptrs * _polybench_alloc_table = (void * )0;
static size_t polybench_inter_array_padding_sz = 0;
/* Timer code (gettimeofday). */
double polybench_t_start, polybench_t_end;
/* Timer code (RDTSC). */
unsigned long long int polybench_c_start, polybench_c_end;
static double rtclock()
{
	double _ret_val_0;
	_ret_val_0=0;
	return _ret_val_0;
}

void polybench_flush_cache()
{
	int cs = (32770*1024)/sizeof (double);
	double * flush = (double * )calloc(cs, sizeof (double));
	int i;
	double tmp = 0.0;
	#pragma cetus private(i) 
	#pragma loop name polybench_flush_cache#0 
	#pragma cetus reduction(+: tmp) 
	#pragma cetus parallel 
	#pragma omp parallel for if((10000<(1L+(3L*cs)))) private(i) reduction(+: tmp)
	for (i=0; i<cs; i ++ )
	{
		tmp+=flush[i];
	}
	((tmp<=10.0) ? ((void)0) : __assert_fail("tmp <= 10.0", "polybench.c", 123, __PRETTY_FUNCTION__));
	free(flush);
	return ;
}

/* ! POLYBENCH_PAPI */
void polybench_prepare_instruments()
{
	polybench_flush_cache();
	return ;
}

void polybench_timer_start()
{
	polybench_prepare_instruments();
	polybench_t_start=rtclock();
	return ;
}

void polybench_timer_stop()
{
	polybench_t_end=rtclock();
	return ;
}

void polybench_timer_print()
{
	printf("%0.6f\n", polybench_t_end-polybench_t_start);
	return ;
}

/*

 These functions are used only if the user defines a specific
 * inter-array padding. It grows a global structure,
 * _polybench_alloc_table, which keeps track of the data allocated via
 * polybench_alloc_data (on which inter-array padding is applied), so
 * that the original, non-shifted pointer can be recovered when
 * calling polybench_free_data.


*/
static void *xmalloc(size_t alloc_sz)
{
	void * ret = (void * )0;
	/* By default, post-pad the arrays. Safe behavior, but likely useless. */
	size_t padded_sz = alloc_sz+polybench_inter_array_padding_sz;
	int err = posix_memalign( & ret, 4096, padded_sz);
	polybench_inter_array_padding_sz+=0;
	if (( ! ret)||err)
	{
		fprintf(stderr, "[PolyBench] posix_memalign: cannot allocate memory");
		exit(1);
	}
	/*
	Safeguard: this is invoked only if polybench.c has been compiled
	     with inter-array padding support from polybench.h. If so, move
	     the starting address of the allocation and return it to the
	     user. The original pointer is registered in an allocation table
	     internal to polybench.c. Data must then be freed using
	     polybench_free_data, which will inspect the allocation table to
	     free the original pointer.
	*/
	return ret;
}

void polybench_free_data(void * ptr)
{
	free(ptr);
	return ;
}

void *polybench_alloc_data(unsigned long long int n, int elt_size)
{
	/* FIXME: detect overflow! */
	size_t val = n;
	void * ret = xmalloc(val);
	val*=elt_size;
	return ret;
}
