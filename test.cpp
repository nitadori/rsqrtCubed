#include <x86intrin.h>

#if 0
void rsqrtCubed_avx512(const double * __restrict x, double * __restrict y, const int N){
	const __m512d c1 = _mm512_set1_pd(3./2.);
	for(int i=0; i<N; i+=8){
		__m512d r2  = _mm512_loadu_pd(x+i);
		__m512d r2c = _mm512_mul_pd(r2, c1);
		__m512d ri1 = _mm512_rsqrt14_pd(r2);
		__m512d ri2 = _mm512_mul_pd(ri1, ri1);
		__m512d ri3 = _mm512_mul_pd(ri2, ri1);
		__m512d hc  = _mm512_fnmadd_pd(r2c, ri2, c1); // c - a * b
#if 0
		const __m512d cc = _mm512_set1_pd(5./6.);
		__m512d hc2 = _mm512_mul_pd(hc, cc);
#else
		const __m512d c2 = _mm512_set1_pd(5./4.);
		__m512d r2c2 = _mm512_mul_pd(r2, c2);
		__m512d hc2  = _mm512_fnmadd_pd(r2c2, ri2, c2); // c - a * b
#endif
		__m512d poly  = _mm512_fmadd_pd(hc, hc2, hc);
		__m512d ri3p  = _mm512_fmadd_pd(poly, ri3, ri3);

		_mm512_storeu_pd(y+i, ri3p);
	}
}
#else
// #define RSQRTCUBED_MINLAT
void rsqrtCubed_avx512(const double * __restrict xs, double * __restrict ys, const int N){
	/*
	 * x : r^2
	 * y : r^-1
	 * z : r^-3
	 * z1 := z0 * (1 + ah + bh^2), where
	 *   a = 3/2
	 *   b = 15/8
	 *   h = 1 - x * y0^2
	 * z1 = z0 + (z0 * h) * (a + b * h)
	 * a+b*h = (a+b) - (b*x)*y0^2, for minlat mode
	 */

	const __m512d one  = _mm512_set1_pd(1.0);
	const __m512d b    = _mm512_set1_pd(15./8.);
#ifdef RSQRTCUBED_MINLAT
	const __m512d apb  = _mm512_set1_pd(27./8.);
#else
	const __m512d a    = _mm512_set1_pd( 3./2.);
#endif

	for(int i=0; i<N; i+=8){
		__m512d x   = _mm512_loadu_pd(xs + i);

#ifdef RSQRTCUBED_MINLAT
		__m512d bx  = _mm512_mul_pd(b, x);
#endif

		__m512d y   = _mm512_rsqrt14_pd(x);
		
		__m512d y2  = _mm512_mul_pd(y, y);

		__m512d z   = _mm512_mul_pd(y, y2);
		__m512d h   = _mm512_fnmadd_pd(x,  y2, one); // c - a * b
#ifdef RSQRTCUBED_MINLAT
		__m512d abh = _mm512_fnmadd_pd(bx, y2, apb); // c - a * b
#else
		__m512d abh = _mm512_fmadd_pd(b, h, a); // a + b*h
#endif

		__m512d zh  = _mm512_mul_pd(z, h);

		__m512d zc  = _mm512_fmadd_pd(zh, abh, z);

		_mm512_storeu_pd(ys + i, zc);
	}
}
#endif

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<algorithm>

int main(){
	enum{
		N = 1000,
	};

	static double x[N], y0[N], y1[N], err[N] __attribute__((aligned(256)));

	srand48(20210309);
	for(int i=0; i<N; i++){
		x[i] = std::exp(10.*(drand48() - 0.5));
	}
	for(int i=0; i<N; i++){
		double ri = 1.0 / std::sqrt(x[i]);
		y0[i] = ri * ri * ri;
	}
	rsqrtCubed_avx512(x, y1, N);

	auto rel_err = [](auto val, auto ref){
		return (val - ref) / ref;
	};

	for(int i=0; i<N; i++){
		err[i] = rel_err(y0[i], y1[i]);
		printf("%e %e %e %e\n", x[i], y0[i], y1[i], err[i]);
	}

#if 0
	std::sort(err, err+N);
	for(int i=0; i<N; i++){
		printf("%f %e\n", i * (1.0/N), err[i]);
	}
#endif

	return 0;
}
