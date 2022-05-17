#include <x86intrin.h>

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<algorithm>

#include "timer.hpp"

struct Body{
	double x, y, z, m;
};

struct Acceleration{
	double ax, ay, az, pot;
};

__attribute__((noinline))
void nbody_ref(
	const int n,
	const double eps2,
	const Body body[],
	Acceleration acc[])
{
	for(int i=0; i<n; i++){ 
		const double xi=body[i].x, yi=body[i].y, zi=body[i].z;
		double ax=0, ay=0, az=0, pot=0;

		for(int j=0; j<n; j++){
			double dx = body[j].x - xi;
			double dy = body[j].y - yi;
			double dz = body[j].z - zi;

			double r2 = eps2 + dx*dx;
			r2 += dy*dy;
			r2 += dz*dz;

			double ri = 1.0 / sqrt(r2);

			double mri = body[j].m * ri;
			double ri2 = ri * ri;

			double mri3 = mri * ri2;

			ax += mri3 * dx;
			ay += mri3 * dy;
			az += mri3 * dz;
			pot -= mri;
		}
		acc[i] = {ax, ay, az, pot};
	}
	return ;
}

inline void transpose_4ymm_pd(__m256d &r0, __m256d &r1, __m256d &r2, __m256d &r3){
	__m256d tmp0 = _mm256_unpacklo_pd(r0, r1); // | r12 | r02 | r10 | r00 |
	__m256d tmp1 = _mm256_unpackhi_pd(r0, r1); // | r13 | r03 | r11 | r01 |
	__m256d tmp2 = _mm256_unpacklo_pd(r2, r3); // | r32 | r22 | r30 | r20 |
	__m256d tmp3 = _mm256_unpackhi_pd(r2, r3); // | r33 | r23 | r31 | r21 |
	r0 = _mm256_permute2f128_pd(tmp0, tmp2, (0)+(2<<4));
	r1 = _mm256_permute2f128_pd(tmp1, tmp3, (0)+(2<<4));
	r2 = _mm256_permute2f128_pd(tmp0, tmp2, (1)+(3<<4));
	r3 = _mm256_permute2f128_pd(tmp1, tmp3, (1)+(3<<4));
}

inline void rsqrt_nr_2ymm_pd(__m256d &x0, __m256d &x1){
	const __m256 half = _mm256_set1_ps(0.5f);
	const __m256 c1p5 = _mm256_set1_ps(1.5f);

	__m128 lo = _mm256_cvtpd_ps(x0);
	__m128 hi = _mm256_cvtpd_ps(x1);

	__m256 x = _mm256_set_m128(hi, lo);
	__m256 y = _mm256_rsqrt_ps(x);
	__m256 xh = _mm256_mul_ps(x, half); 
	__m256 y2 = _mm256_mul_ps(y, y); 
	__m256 p  = _mm256_fnmadd_ps(xh, y2, c1p5);

	y = _mm256_mul_ps(y, p); 

	lo = _mm256_extractf128_ps(y, 0);
	hi = _mm256_extractf128_ps(y, 1);
	
	x0 = _mm256_cvtps_pd(lo);
	x1 = _mm256_cvtps_pd(hi);
}

inline __m256d rsqrtCubed(
		const __m256d x,
		const __m256d y,
		const __m256d m)
{
	const __m256d one  = _mm256_set1_pd(1.0);
	const __m256d a    = _mm256_set1_pd(3./2.);
	const __m256d b    = _mm256_set1_pd(15./8.);

	__m256d y2  = _mm256_mul_pd(y, y);
	__m256d my  = _mm256_mul_pd(m, y);

	__m256d z   = _mm256_mul_pd(my, y2);
	__m256d h   = _mm256_fnmadd_pd(x,  y2, one); // c - a * b
	__m256d abh = _mm256_fmadd_pd(b, h, a); // a + b*h

	__m256d zh  = _mm256_mul_pd(z, h);

	// abh = a;  // force 2nd order

	__m256d z1  = _mm256_fmadd_pd(zh, abh, z);

	return z1;
}

__attribute__((noinline))
void nbody_m256d(
	const int n,
	const double eps2_sd,
	const Body body[],
	Acceleration acc[])
{
	const __m256d eps2 = _mm256_set1_pd(eps2_sd);

	for(int i=0; i<n; i+=4){
		__m256d xi = _mm256_load_pd(&body[i+0].x);
		__m256d yi = _mm256_load_pd(&body[i+1].x);
		__m256d zi = _mm256_load_pd(&body[i+2].x);
		__m256d mi = _mm256_load_pd(&body[i+3].x);
		transpose_4ymm_pd(xi, yi, zi, mi);

		__m256d ax, ay, az, pot;
		ax = ay = az = pot = _mm256_set1_pd(0);

		for(int j=0; j<n; j+=2){
			__m256d dx_0 = _mm256_sub_pd(xi, _mm256_set1_pd(body[j+0].x));
			__m256d dy_0 = _mm256_sub_pd(yi, _mm256_set1_pd(body[j+0].y));
			__m256d dz_0 = _mm256_sub_pd(zi, _mm256_set1_pd(body[j+0].z));
			__m256d mj_0 = _mm256_set1_pd(body[j+0].m);

			__m256d r2_0  = 
				_mm256_fmadd_pd(dz_0, dz_0, 
					_mm256_fmadd_pd(dy_0, dy_0,
						_mm256_fmadd_pd(dx_0, dx_0, eps2)));

			__m256d dx_1 = _mm256_sub_pd(xi, _mm256_set1_pd(body[j+1].x));
			__m256d dy_1 = _mm256_sub_pd(yi, _mm256_set1_pd(body[j+1].y));
			__m256d dz_1 = _mm256_sub_pd(zi, _mm256_set1_pd(body[j+1].z));
			__m256d mj_1 = _mm256_set1_pd(body[j+1].m);

			__m256d r2_1  = 
				_mm256_fmadd_pd(dz_1, dz_1, 
					_mm256_fmadd_pd(dy_1, dy_1,
						_mm256_fmadd_pd(dx_1, dx_1, eps2)));

			__m256d rinv_0 = r2_0;
			__m256d rinv_1 = r2_1;
			rsqrt_nr_2ymm_pd(rinv_0, rinv_1);

		       __m256d mri3_0 = rsqrtCubed(r2_0, rinv_0, mj_0);

		       ax  = _mm256_fnmadd_pd(mri3_0, dx_0, ax);
		       ay  = _mm256_fnmadd_pd(mri3_0, dy_0, ay);
		       az  = _mm256_fnmadd_pd(mri3_0, dz_0, az);
		       pot = _mm256_fnmadd_pd(mri3_0, r2_0, pot);

		       __m256d mri3_1 = rsqrtCubed(r2_1, rinv_1, mj_1);

		       ax  = _mm256_fnmadd_pd(mri3_1, dx_1, ax);
		       ay  = _mm256_fnmadd_pd(mri3_1, dy_1, ay);
		       az  = _mm256_fnmadd_pd(mri3_1, dz_1, az);
		       pot = _mm256_fnmadd_pd(mri3_1, r2_1, pot);
		}

		transpose_4ymm_pd(ax, ay, az, pot);
		_mm256_store_pd(&acc[i+0].ax, ax);
		_mm256_store_pd(&acc[i+1].ax, ay);
		_mm256_store_pd(&acc[i+2].ax, az);
		_mm256_store_pd(&acc[i+3].ax, pot);
	}
}

inline __m256d rsqrtCubed_x5(
		const __m256d x,
		const __m256d m)
{
	const __m256d one  = _mm256_set1_pd(1.0);
	const __m256d a  = _mm256_set1_pd(  3./ 2.);
	const __m256d b  = _mm256_set1_pd( 15./ 8.);
	const __m256d c  = _mm256_set1_pd( 35./ 16.);
	const __m256d d  = _mm256_set1_pd(315./128.);

	__m256d y = _mm256_cvtps_pd(
			_mm_rsqrt_ps(
				_mm256_cvtpd_ps(x)));

	__m256d y2  = _mm256_mul_pd(y, y);
	__m256d my  = _mm256_mul_pd(m, y);

	__m256d z   = _mm256_mul_pd(my, y2);
	__m256d h   = _mm256_fnmadd_pd(x,  y2, one); // c - a * b

	__m256d poly = _mm256_fmadd_pd(d, h, c); // d*h + c
	poly = _mm256_fmadd_pd(poly, h, b); // (d*h + c)*h + b
	poly = _mm256_fmadd_pd(poly, h, a); // ((d*h + c)*h + b)*h + a

	__m256d zh  = _mm256_mul_pd(z, h);

	__m256d z1  = _mm256_fmadd_pd(zh, poly, z);

	return z1;
}

__attribute__((noinline))
void nbody_m256d_nj1(
	const int n,
	const double eps2_sd,
	const Body body[],
	Acceleration acc[])
{
	const __m256d eps2 = _mm256_set1_pd(eps2_sd);

	for(int i=0; i<n; i+=4){
		__m256d xi = _mm256_load_pd(&body[i+0].x);
		__m256d yi = _mm256_load_pd(&body[i+1].x);
		__m256d zi = _mm256_load_pd(&body[i+2].x);
		__m256d mi = _mm256_load_pd(&body[i+3].x);
		transpose_4ymm_pd(xi, yi, zi, mi);

		__m256d ax, ay, az, pot;
		ax = ay = az = pot = _mm256_set1_pd(0);

		for(int j=0; j<n; j+=1){
			__m256d dx_0 = _mm256_sub_pd(xi, _mm256_set1_pd(body[j+0].x));
			__m256d dy_0 = _mm256_sub_pd(yi, _mm256_set1_pd(body[j+0].y));
			__m256d dz_0 = _mm256_sub_pd(zi, _mm256_set1_pd(body[j+0].z));
			__m256d mj_0 = _mm256_set1_pd(body[j+0].m);

			__m256d r2_0  = 
				_mm256_fmadd_pd(dz_0, dz_0, 
					_mm256_fmadd_pd(dy_0, dy_0,
						_mm256_fmadd_pd(dx_0, dx_0, eps2)));

		       __m256d mri3_0 = rsqrtCubed_x5(r2_0, mj_0);

		       ax  = _mm256_fnmadd_pd(mri3_0, dx_0, ax);
		       ay  = _mm256_fnmadd_pd(mri3_0, dy_0, ay);
		       az  = _mm256_fnmadd_pd(mri3_0, dz_0, az);
		       pot = _mm256_fnmadd_pd(mri3_0, r2_0, pot);
		}

		transpose_4ymm_pd(ax, ay, az, pot);
		_mm256_store_pd(&acc[i+0].ax, ax);
		_mm256_store_pd(&acc[i+1].ax, ay);
		_mm256_store_pd(&acc[i+2].ax, az);
		_mm256_store_pd(&acc[i+3].ax, pot);
	}
}

#if 0
int main(){
	__attribute__((aligned(64)))
	double a[4][4] = {
		{1, 2, 3, 4},
		{5, 6, 7, 8},
		{9, 10, 11, 12},
		{13, 14, 15, 16},
	};
	for(auto &b : a){ 
		for(auto c : b)  printf("%f, ", c);
		puts("");
	}

	transpose_4ymm_pd(*(__m256d*)a[0], *(__m256d*)a[1], *(__m256d*)a[2], *(__m256d*)a[3]); 

	puts("");
	for(auto &b : a){ 
		for(auto c : b)  printf("%f, ", c);
		puts("");
	}
}
#endif


#if 0
int main(){
	__attribute__((aligned(64)))
	double x[8], y[8];

	srand48(20220517);
	for(int i=0; i<8; i++){
		x[i] = y[i] = drand48() * drand48();
	}
	rsqrt_nr_2ymm_pd(*(__m256d *)&y[0], *(__m256d *)&y[4]);

	for(int i=0; i<8; i++){
		auto err = 1.0 - x[i]*(y[i]*y[i]);
		printf("%e\n", err);
	}
}
#endif

#if 1
int main(){
	enum{
		N = 2048,
	};

	const double eps = 1./256.;
	const double eps2 = eps*eps;

	static Body          body    [N];
	static Acceleration  acc     [N];
	static Body          dbl_body[N];
	static Acceleration  dbl_acc [N];

	srand48(20220517);
	for(int i=0; i<N; i++){
		body[i].x = drand48() - 0.5;
		body[i].y = drand48() - 0.5;
		body[i].z = drand48() - 0.5;
		body[i].m = (1./N) * (drand48() + 0.5);

		const Body &b = body[i];
		dbl_body[i] = {b.x, b.y, b.z, b.m};
	}

	nbody_ref(N, eps2, dbl_body, dbl_acc);

	auto verify = [=](auto kernel){
		kernel(N, eps2, body, acc);

		double fx=0, fy=0, fz=0, ff=0;
		for(int i=0; i<N; i++){
			auto norm = [](auto x, auto y, auto z){
				return std::sqrt(x*x + y*y + z*z);
			};
			fx += (double)body[i].m * acc[i].ax;
			fy += (double)body[i].m * acc[i].ay;
			fz += (double)body[i].m * acc[i].az;
			ff += (double)body[i].m * norm(acc[i].ax, acc[i].ay, acc[i].az);
		}

		printf("(%e, %e, %e) : %e\n", fx, fy, fz, ff);

		auto rel_err = [](auto val, auto ref){
			double dx = val.ax - ref.ax;
			double dy = val.ay - ref.ay;
			double dz = val.az - ref.az;

			double num = dx*dx + dy*dy + dz*dz;
			double den = ref.ax*ref.ax + ref.ay*ref.ay + ref.az*ref.az;

			return sqrt(num / den);
		};

		double err_max = 0.0, err_min = 1.0;
		double dp_max = 0.0, dp_min = 1.0;
		for(int i=0; i<N; i++){
			double e = rel_err(acc[i], dbl_acc[i]);
			double dp = fabs(acc[i].pot - dbl_acc[i].pot);
			// err[i] = e;
			// printf("%e %e %e %e\n", x[i], y0[i], y1[i], err[i]);
			err_max = std::max(err_max, e);
			err_min = std::min(err_min, e);
			dp_max = std::max(dp_max, dp);
			dp_min = std::min(dp_min, dp);

			if(!std::isfinite(e) || fabs(e) > 1.e-5){
				printf("%4d %+e (%e %e %e)\n", i, e, acc[i].ax, acc[i].ay, acc[i].az);
			}
			if(!std::isfinite(dp) || fabs(dp) > 1.e-5){
				printf("%4d %+e (%e %e)\n", i, dp, acc[i].pot, dbl_acc[i].pot);
			}
		}
		printf("err in [%e, %e], perr in [%e, %e]\n", err_min, err_max, dp_min, dp_max);
		puts("");
	};

	verify(nbody_m256d);
	verify(nbody_m256d_nj1);

	auto warmup = [=](auto kernel, int ntimes=100){
		for(int j=0; j<ntimes; j++){
			kernel(N, eps2, body, acc);
		}
	};

	warmup(nbody_m256d);

	auto benchmark = [=](auto kernel, int ntimes=10){
		double nsecs[ntimes];
		for(int j=0; j<ntimes; j++){
			for(int i=0; i<N; i++){
				acc[i] = {0,0,0};
			}

			auto tick0 = get_utime();
			kernel(N, eps2, body, acc);
			auto tick1 = get_utime();

			double dt = tick2second(tick1 - tick0);
			double iter = N/4.0 * N;
			nsecs[j] = dt/iter * 1.e9;
		}
		for(int j=0; j<ntimes; j++){
#ifdef __aarch64__ 
			// Just assume 2.0 GHz of Fugaku
			printf("%f nsec/loop, %f cycles\n", nsecs[j], 2.0*nsecs[j]);
#else
			double cycle = nsecs[j] * 4.3; // 4.3 GHz?
			double eff = 100.0 * 11.25 / cycle; // 22.5 fp-inst/loop?
			double Gflops = 38.*4. / nsecs[j];
			printf("%f nsec/loop, %f cycles, %f%%, %f Gflops\n", nsecs[j], cycle, eff, Gflops);
#endif
		}
		puts("");
		fflush(stdout);
	};
	
	benchmark(nbody_ref);
	benchmark(nbody_m256d);
	benchmark(nbody_m256d_nj1);

	return 0;
}
#endif

