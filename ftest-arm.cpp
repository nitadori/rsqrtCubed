#include <cmath>
#include <arm_sve.h>

__attribute__((noinline))
void rsqrtCubed_autovec(
		const float * __restrict xs, 
		float       * __restrict zs, 
		const float * __restrict ms, 
		const int                N)
{
#pragma loop unroll 4
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

#pragma loop unroll 4
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

__attribute__((noinline))
void rsqrtCubed_swp(
		const float * __restrict xs, 
		float       * __restrict zs, 
		const float * __restrict ms, 
		const int                N)
{
	const svfloat32_t one  = svdup_f32(1.0);
	const svfloat32_t b    = svdup_f32(15./8.);
	const svfloat32_t a    = svdup_f32( 3./2.);

	svbool_t p0 = svptrue_b32();

	svfloat32_t s1_x;
	svfloat32_t s2_x, s2_y, s2_m;
	svfloat32_t s3_x, s3_y2, s3_my;
	svfloat32_t s4_h, s4_z;
	svfloat32_t s5_abh, s5_zh, s5_z;
	svfloat32_t s6_z1;

#pragma loop unroll 4
	for(int i=-6*16; i<N; i+=16){
		svst1_f32(p0, zs + i, s6_z1);

		s6_z1  = svmad_f32_x(p0, s5_zh, s5_abh, s5_z);

		s5_zh  = svmul_f32_x(p0, s4_z, s4_h);
		s5_abh = svmad_f32_x(p0, b, s4_h, a);
		s5_z   = s4_z;

		s4_z   = svmul_f32_x(p0, s3_my, s3_y2);
		s4_h   = svmsb_f32_x (p0, s3_x, s3_y2, one);

		s3_y2  = svmul_f32_x(p0, s2_y, s2_y);
		s3_my  = svmul_f32_x(p0, s2_m, s2_y);
		s3_x   = s2_x;

		s2_m   = svld1_f32(p0, ms + (i + 5*16));
		s2_y   = svrsqrte_f32(s1_x);
		s2_x   = s1_x;

		s1_x   = svld1_f32(p0, xs + (i + 6*16));
	}
}

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<cassert>
#include<algorithm>

#include "timer.hpp"

int main(){
	enum{
		N = 2048+32,
		NLOOP = 1000,
	};

#if 0
	static float x[N], y0[N], pad0[128], y1[N], pad1[128], mass[N], err[N] __attribute__((aligned(256)));
	memcpy(pad0, pad1, 128*4);
#else
	static float x[N], y0[N], yy[N+256], *y1=yy+128, mass[N], err[N] __attribute__((aligned(256)));
#endif

	srand48(20210309);
	for(int i=0; i<N; i++){
		x[i]    = std::exp(10.*(drand48() - 0.5));
		mass[i] = drand48() + 0.125;
	}
	for(int i=0; i<N; i++){
		double ri = 1.0 / std::sqrt( (double)(x[i]) );
		y0[i] = (mass[i] * ri) * (ri * ri);
	}
	// rsqrtCubed_autovec(x, y1, mass, N);

	auto rel_err = [](double val, double ref){
		assert(ref > 0.0);
		return (val - ref) / ref;
	};

	auto verify = [=](auto kernel){
		for(int i=0; i<N; i++) y1[i] = 0.0/0.0;
		kernel(x, y1, mass, N);

		double err_max = 0.0, err_min = 0.0;
		for(int i=0; i<N; i++){
			double e = rel_err(y1[i], y0[i]);
			err[i] = e;
			// printf("%e %e %e %e\n", x[i], y0[i], y1[i], err[i]);
			err_max = std::max(err_max, e);
			err_min = std::min(err_min, e);

			if(!std::isfinite(e) || fabs(e) > 1.e-6){
				printf("%4d %+e %e %e\n", i, e, y1[i], y0[i]);
				(void)err[i];
			}
		}
		printf("err in [%e, %e]\n", err_min, err_max);
	};
#if 0
	verify(rsqrtCubed_swp);
	for(int i=0; i<16; i++){
		printf("%4d %+e %e %e\n", i, err[i], y1[i], y0[i]);
	}
#endif
	verify(rsqrtCubed_autovec);
	verify(rsqrtCubed_sve);
	verify(rsqrtCubed_swp);

#if 0
	std::sort(err, err+N);
	for(int i=0; i<N; i++){
		printf("%f %e\n", i * (1.0/N), err[i]);
	}
#endif


	// time measurement
	//
	
	auto benchmark = [=](auto kernel, int ntimes=10){
		double nsecs[ntimes];
		for(int j=0; j<ntimes; j++){
			auto tick0 = get_utime();
			for(int k=0; k<NLOOP; k++){
				kernel(x, y1, mass, N);
			}
			auto tick1 = get_utime();
			double dt = tick2second(tick1 - tick0);
			double iter = N/16.0 * NLOOP;
			nsecs[j] = dt/iter * 1.e9;
		}
		for(int j=0; j<ntimes; j++){
#ifdef __aarch64__ 
			// Just assume 2.0 GHz of Fugaku
			printf("%f nsec/loop, %f cycles\n", nsecs[j], 2.0*nsecs[j]);
#else
			printf("%f nsec/loop\n", nsecs[j]);
#endif
		}
		puts("");
		fflush(stdout);
	};
	benchmark(rsqrtCubed_autovec);
	benchmark(rsqrtCubed_sve);
	benchmark(rsqrtCubed_swp);

	return 0;
}
