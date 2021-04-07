#include <cmath>
#include <arm_sve.h>

__attribute__((noinline))
void rsqrtCubed_autovec(
		const float * __restrict xs, 
		float       * __restrict zs, 
		const float * __restrict ms, 
		const int                N)
{
	for(int i=0; i<N; i++){
		float x = xs[i];
		float m = ms[i];
		float y = 1.0f / sqrtf(x);
		float z = (m*y)*(y*y);
		zs[i] = z;
	}
}

__attribute__((noinline))
void rsqrtCubed_sve(
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
	/*
	 * msb(p,a,b,c) : c - a*b
	 * mad(p,a,b,c) : c + a*b
	 * mls(p,a,b,c) : a - b*c
	 * mla(p,a,b,c) : a + b*c
	 */

	const svfloat32_t one  = svdup_f32(1.0);
	const svfloat32_t b    = svdup_f32(15./8.);
	const svfloat32_t a    = svdup_f32( 3./2.);

	svbool_t p0 = svptrue_b32();

	for(int i=0; i<N; i+=16){
		svfloat32_t x   = svld1_f32(p0, xs + i);
		svfloat32_t m   = svld1_f32(p0, ms + i);

		svfloat32_t y   = svrsqrte_f32(x);
		
		svfloat32_t y2  = svmul_f32_x(p0, y, y);
		svfloat32_t my  = svmul_f32_x(p0, m, y);

		svfloat32_t z   = svmul_f32_x(p0, my, y2);
#if 0
		svfloat32_t h   = svmls_f32_x(p0, one, x,  y2);
#else
		svfloat32_t h   = svmsb_f32_x(p0, x,  y2, one);
#endif

		svfloat32_t zh  = svmul_f32_x(p0, z, h);
#if 1
		svfloat32_t abh = svmad_f32_x(p0, b, h, a); // a + b*h
#else
		svfloat32_t abh = svmla_f32_x(p0, a, b, h); // a + b*h
#endif


		// abh    = _mm512_set1_ps( 3./2.);  // force 2nd order
		
		svfloat32_t z1  = svmad_f32_x(p0, zh, abh, z);

		svst1_f32(p0, zs + i, z1);
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
		NLOOP = 1000,
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
	rsqrtCubed_autovec(x, y1, mass, N);

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
			rsqrtCubed_autovec(x, y1, mass, N);
		}
		auto tick1 = get_utime();
		double dt = tick2second(tick1 - tick0);
		double iter = N/16.0 * NLOOP;
		double nsec = dt/iter * 1.e9;
		printf("%f nsec/loop\n", nsec);
	}

	return 0;
}
