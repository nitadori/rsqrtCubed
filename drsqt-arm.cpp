#include <cmath>
#include <arm_sve.h>

__attribute__((noinline))
void drsqrt_autovec(
		const double * __restrict xs, 
		double       * __restrict ys, 
		const int                N)
{
	for(int i=0; i<N; i++){
		ys[i] = 1.0 / sqrt(xs[i]);
	}
}

__attribute__((noinline))
void drsqrt_nr(
		const double * __restrict xs, 
		double       * __restrict ys, 
		const int                N)
{
	svbool_t p0 = svptrue_b64();

	const svfloat64_t half   = svdup_f64(1./2.);

	svfloat64_t x   = svld1_f64(p0, xs + 0);
#pragma loop unroll 4
	for(int i=0; i<N; i+=8){
		// svfloat64_t x   = svld1_f64(p0, xs + i);
		svfloat64_t x2  = svmul_f64_x(p0, x, half);
		svfloat64_t y   = svrsqrte_f64(x);
		x   = svld1_f64(p0, &xs[i+8]);

		svfloat64_t y2  = svmul_f64_x(p0, y, y);
		svfloat64_t h2  = svmsb_f64_x(p0, x2,  y2, half);
		y  = svmad_f64_x(p0, y, h2, y);

		y2  = svmul_f64_x(p0, y, y);
		h2  = svmsb_f64_x(p0, x2,  y2, half);
		y  = svmad_f64_x(p0, y, h2, y);

		y2  = svmul_f64_x(p0, y, y);
		h2  = svmsb_f64_x(p0, x2,  y2, half);
		y  = svmad_f64_x(p0, y, h2, y);

		svst1_f64(p0, ys + i, y);
	}
}

__attribute__((noinline))
void drsqrt_sve(
		const double * __restrict xs, 
		double       * __restrict ys, 
		const int                N)
{
	svbool_t p0 = svptrue_b64();

	const svfloat64_t one  = svdup_f64(1.0);
	const svfloat64_t c1   = svdup_f64(1./2.);
	const svfloat64_t c2   = svdup_f64(3./8.);
	const svfloat64_t c3   = svdup_f64(5./16.);
	const svfloat64_t c4   = svdup_f64(35./128.);
	const svfloat64_t c5   = svdup_f64(63./256.);
	const svfloat64_t c6   = svdup_f64(231./1024.);

	svfloat64_t x   = svld1_f64(p0, xs + 0);
#pragma loop unroll 4
	for(int i=0; i<N; i+=8){
		// svfloat64_t x   = svld1_f64(p0, xs + i);
		svfloat64_t y   = svrsqrte_f64(x);
		
		svfloat64_t y2  = svmul_f64_x(p0, y, y);
		svfloat64_t h   = svmsb_f64_x(p0, x,  y2, one);
		x   = svld1_f64(p0, &xs[i+8]);
		svfloat64_t yh  = svmul_f64_x(p0, y, h);

		svfloat64_t d5 = svmad_f64_x(p0, c6, h, c5); // c5 + h*c6
		svfloat64_t d4 = svmad_f64_x(p0, d5, h, c4); // c4 + h*d5
		svfloat64_t d3 = svmad_f64_x(p0, d4, h, c3); // c3 + h*d4
		svfloat64_t d2 = svmad_f64_x(p0, d3, h, c2); // c2 + h*d3
		svfloat64_t d1 = svmad_f64_x(p0, d2, h, c1); // c1 + h*d2

		svfloat64_t y1  = svmad_f64_x(p0, yh, d1, y);

		svst1_f64(p0, ys + i, y1);
	}
}

__attribute__((noinline))
void drsqrt_sve_thr( // tree-height-reduction
		const double * __restrict xs, 
		double       * __restrict ys, 
		const int                N)
{
	svbool_t p0 = svptrue_b64();

	const svfloat64_t one  = svdup_f64(1.0);
	const svfloat64_t c1   = svdup_f64(1./2.);
	const svfloat64_t c2   = svdup_f64(3./8.);
	const svfloat64_t c3   = svdup_f64(5./16.);
	const svfloat64_t c4   = svdup_f64(35./128.);
	const svfloat64_t c5   = svdup_f64(63./256.);
	const svfloat64_t c6   = svdup_f64(231./1024.);

	svfloat64_t x   = svld1_f64(p0, xs + 0);
#pragma loop unroll 4
	for(int i=0; i<N; i+=8){
		// svfloat64_t x   = svld1_f64(p0, xs + i);
		svfloat64_t y   = svrsqrte_f64(x);
		
		svfloat64_t y2  = svmul_f64_x(p0, y, y);
		svfloat64_t h   = svmsb_f64_x(p0, x,  y2, one);
		x   = svld1_f64(p0, &xs[i+8]);
		svfloat64_t yh  = svmul_f64_x(p0, y, h);

		svfloat64_t d5 = svmad_f64_x(p0, c6, h, c5); // c5 + h*c6
		svfloat64_t d3 = svmad_f64_x(p0, c4, h, c3); // c3 + h*c4
		svfloat64_t d1 = svmad_f64_x(p0, c2, h, c1); // c1 + h*c2
		svfloat64_t h2  = svmul_f64_x(p0, h, h);

		svfloat64_t e3 = svmad_f64_x(p0, d5, h2, d3); // d3 + h^2*d5
		svfloat64_t e1 = svmad_f64_x(p0, e3, h2, d1); // d1 + h^2*e5

		svfloat64_t y1  = svmad_f64_x(p0, yh, e1, y);

		svst1_f64(p0, ys + i, y1);
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
		N = 1024+16,
		NLOOP = 1000,
	};

	auto err_rsqrt = [](double x, double y){
		double y2 = y*y;
		double yc = fma(y, y, -y2);
		double h  = fma(-x, y2, 1.0);
		double hc = fma(-x, yc, h);

		return hc;
	};

	static double x[N], y1[N], err[N] __attribute__((aligned(256)));

	srand48(20210309);
	for(int i=0; i<N; i++){
		x[i]    = std::exp(10.*(drand48() - 0.5));
	}

	auto verify = [=](auto kernel){
		for(int i=0; i<N; i++) y1[i] = 0.0/0.0;

		kernel(x, y1,N);

		double err_max = -1.0/0.0, err_min = +1.0/0.0;
		for(int i=0; i<N; i++){
			double e = err_rsqrt(x[i], y1[i]);

			err_max = std::max(err_max, e);
			err_min = std::min(err_min, e);

			if(!std::isfinite(e) || fabs(e) > 1.e-6){
				printf("%4d %+e %e\n", i, e, y1[i]);
				(void)err[i];
			}
		}
		printf("err in [%e, %e]\n", err_min, err_max);
	};

	verify(drsqrt_autovec);
	verify(drsqrt_nr);
	verify(drsqrt_sve);
	verify(drsqrt_sve_thr);

	auto benchmark = [=](auto kernel, int ntimes=10){
		double nsecs[ntimes];
		for(int j=0; j<ntimes; j++){
			auto tick0 = get_utime();
			for(int k=0; k<NLOOP; k++){
				kernel(x, y1, N);
			}
			auto tick1 = get_utime();
			double dt = tick2second(tick1 - tick0);
			double iter = N/8.0 * NLOOP;
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

	benchmark(drsqrt_autovec);
	benchmark(drsqrt_nr);
	benchmark(drsqrt_sve);
	benchmark(drsqrt_sve_thr);

	return 0;
}
