#include <x86intrin.h>

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

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<algorithm>

int main(){
	enum{
		N = 1000,
	};

	static double x[N], y0[N], y1[N], err[N];

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

	std::sort(err, err+N);
	for(int i=0; i<N; i++){
		printf("%f %e\n", i * (1.0/N), err[i]);
	}

	return 0;
}
