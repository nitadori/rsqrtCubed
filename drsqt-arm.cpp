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

	return 0;
}
