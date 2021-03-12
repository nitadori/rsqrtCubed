#include <x86intrin.h>

#define RSQRTCUBED_MINLAT
__attribute__((noinline))
void rsqrtCubed_avx512(
		const float * __restrict xs, 
		float       * __restrict zs, 
		const float * __restrict ms, 
		const int                N)
{
	/*
	 * x : r^2
	 * y : r^-1
	 * z : m * r^-3
	 * z1 := z0 * (1 + ah + bh^2), where
	 *   a = 3/2
	 *   b = 15/8
	 *   h = 1 - x * y0^2
	 * z1 = z0 + (z0 * h) * (a + b * h)
	 * a+b*h = (a+b) - (b*x)*y0^2, for minlat mode
	 */

	const __m512 one  = _mm512_set1_ps(1.0);
	const __m512 b    = _mm512_set1_ps(15./8.);
#ifdef RSQRTCUBED_MINLAT
	const __m512 apb  = _mm512_set1_ps(27./8.);
#else
	const __m512 a    = _mm512_set1_ps( 3./2.);
#endif

	for(int i=0; i<N; i+=16){
		__m512 x   = _mm512_loadu_ps(xs + i);
		__m512 m   = _mm512_loadu_ps(ms + i);

#ifdef RSQRTCUBED_MINLAT
		__m512 bx  = _mm512_mul_ps(b, x);
#endif

		__m512 y   = _mm512_rsqrt14_ps(x);
		
		__m512 y2  = _mm512_mul_ps(y, y);
		__m512 my  = _mm512_mul_ps(m, y);

		__m512 z   = _mm512_mul_ps(my, y2);
		__m512 h   = _mm512_fnmadd_ps(x,  y2, one); // c - a * b
#ifdef RSQRTCUBED_MINLAT
		__m512 abh = _mm512_fnmadd_ps(bx, y2, apb); // c - a * b
#else
		__m512 abh = _mm512_fmadd_ps(b, h, a); // a + b*h
#endif

		__m512 zh  = _mm512_mul_ps(z, h);

		// abh    = _mm512_set1_ps( 3./2.);  // force 2nd order
		
		__m512 z1  = _mm512_fmadd_ps(zh, abh, z);

		_mm512_storeu_ps(zs + i, z1);
	}
}

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<algorithm>

#include "timer.hpp"

int main(){
	enum{
		N = 2048+32,
		NLOOP = 10000,
	};

	static float x[N], y0[N], y1[N], mass[N], err[N] __attribute__((aligned(256)));

	srand48(20210309);
	for(int i=0; i<N; i++){
		x[i]    = std::exp(10.*(drand48() - 0.5));
		mass[i] = drand48() + 0.125;
	}
	for(int i=0; i<N; i++){
		double ri = 1.0 / std::sqrt( (double)(x[i]) );
		y0[i] = (mass[i] * ri) * (ri * ri);
	}
	rsqrtCubed_avx512(x, y1, mass, N);

	auto rel_err = [](auto val, auto ref){
		return (val - ref) / ref;
	};

	double err_max = 0.0, err_min = 0.0;
	for(int i=0; i<N; i++){
		err[i] = rel_err(y0[i], y1[i]);
		// printf("%e %e %e %e\n", x[i], y0[i], y1[i], err[i]);
		double e = err[i];
		err_max = std::max(err_max, e);
		err_min = std::min(err_min, e);
	}
	printf("err in [%e, %e]\n", err_min, err_max);

#if 0
	std::sort(err, err+N);
	for(int i=0; i<N; i++){
		printf("%f %e\n", i * (1.0/N), err[i]);
	}
#endif


	// time measurement
	//
	
	for(int j=0; j<10; j++){
		auto tick0 = get_utime();
		for(int k=0; k<NLOOP; k++){
			rsqrtCubed_avx512(x, y1, mass, N);
		}
		auto tick1 = get_utime();
		double dt = tick2second(tick1 - tick0);
		double iter = N/16.0 * NLOOP;
		double nsec = dt/iter * 1.e9;
		printf("%f nsec/loop\n", nsec);
	}

	return 0;
}
