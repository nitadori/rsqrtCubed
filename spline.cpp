#include <cstdio>
#include <cstdlib>
#include <cmath>
// #include <cstring>
#include <cassert>
#include <algorithm>

#include <x86intrin.h>

#include "timer.hpp"

// __attribute__((noinline))
// __attribute__((always_inline))
inline double pcut(double u0){
	using F = decltype(u0);
#if 0
	using std::min;
	using std::max;
#elif 1
	auto min = fmin;
	auto max = fmax;
#else
	auto min = [](auto a, auto b) { return a<b ? a : b; };
	auto max = [](auto a, auto b) { return a>b ? a : b; };
#endif

	double u  = min(F(2), u0);
	double u2 = u*u;

	double pl = fma(u, -3, 9);
	pl = fma(pl, u2, -20);
	pl = fma(pl, u2, 42);
	pl *= u;
	
	double pr = fma(u, 4, 2);

	double usub = max(F(0), u-1);
	double usub2 = usub  * usub;
	double usub4 = usub2 * usub2;
	double usub5 = usub  * usub4;

	double sum = fma(usub5, pr, pl);

	return F(1./30.) * sum;
} 

// __attribute__((noinline))
// __attribute__((always_inline))
inline double acut(double u0){
	using F = decltype(u0);
#if 0
	using std::min;
	using std::max;
#elif 1
	auto min = fmin;
	auto max = fmax;
#else
	auto min = [](auto a, auto b) { return a<b ? a : b; };
	auto max = [](auto a, auto b) { return a>b ? a : b; };
#endif
	double u  = min(F(2), u0);
	double u2 = u*u;
	double u3 = u*u2;

	double pl = fma(u, 15, -36);
	pl = fma(pl, u2, 40);
	pl *= u3;

	double pr = fma(u, 20, 8);
	pr = fma(pr, u, 2);

	double usub = max(F(0), u-1);
	double usub2 = usub  * usub;
	double usub4 = usub2 * usub2;

	double sum = fma(usub4, -pr, pl);

	return F(1./30.) * sum;
}

struct Body{
	double x, y, z, m;
};

struct Acceleration{
	double ax, ay, az, pot;
};

__attribute__((noinline))
void nbody_spline_ref(
	const int n,
	const double /* rcut2 */, // = 4 eps^2
	const double epsinv,
	const Body body[],
	Acceleration acc[])
{
	for(int i=0; i<n; i++){ 
		const double xi=body[i].x, yi=body[i].y, zi=body[i].z;
		double ax=0, ay=0, az=0, pot=0;

		for(int j=0; j<n; j++){
			if(j == i) continue;

			double dx = body[j].x - xi;
			double dy = body[j].y - yi;
			double dz = body[j].z - zi;

			double r2 = dx*dx + dy*dy + dz*dz;

			double ri = 1.0 / sqrt(r2);

			double mri = body[j].m * ri;
			double ri2 = ri * ri;

			double mri3 = mri * ri2;

			double u = r2 * ri * epsinv;

			mri3 *= acut(u);
			mri  *= pcut(u);

			ax += mri3 * dx;
			ay += mri3 * dy;
			az += mri3 * dz;
			pot -= mri;
		}
		assert(std::isfinite(pot));
		acc[i] = {ax, ay, az, pot};
	}
	return ;
}

__attribute__((noinline))
void nbody_spline_list(
	const int n,
	const double rcut2, // = 4 eps^2
	const double epsinv,
	const Body body[],
	Acceleration acc[])
{
	enum{ LEN = 64, };
	int list[LEN];
	for(int i=0; i<n; i++){ 
		int num = 0;
		const double xi=body[i].x, yi=body[i].y, zi=body[i].z;
		double ax=0, ay=0, az=0, pot=0;

		for(int j=0; j<n; j++){
			double dx = body[j].x - xi;
			double dy = body[j].y - yi;
			double dz = body[j].z - zi;

			double r2 = dx*dx + dy*dy + dz*dz;

			double ri = 1.0 / sqrt(r2);

			double mri = body[j].m * ri;
			double ri2 = ri * ri;

			double mri3 = mri * ri2;

			if(r2 < rcut2){
				list[num++%LEN] = j;
				mri3 = 0.0;
			}

			ax += mri3 * dx;
			ay += mri3 * dy;
			az += mri3 * dz;
			pot -= mri3 * r2;
		}
		assert(std::isfinite(pot));
		// printf("%d : %d\n", i, num);
		assert(num <= LEN);
		for(int jj=0; jj<num; jj++){
			int j = list[jj];

			double dx = body[j].x - xi;
			double dy = body[j].y - yi;
			double dz = body[j].z - zi;

			double r2 = dx*dx + dy*dy + dz*dz;

			if(0.0 == r2) continue;

			double ri = 1.0 / sqrt(r2);

			double mri = body[j].m * ri;
			double ri2 = ri * ri;

			double mri3 = mri * ri2;
			double u = r2 * ri * epsinv;

			mri3 *= acut(u);
			mri  *= pcut(u);

			ax += mri3 * dx;
			ay += mri3 * dy;
			az += mri3 * dz;
			pot -= mri;

		}
		assert(std::isfinite(pot));
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
void nbody_spline_m256d(
	const int n,
	const double rcut2_sd, // = 4 eps^2
	const double epsinv,
	const Body body[],
	Acceleration acc[])
{
	enum{ LEN = 64, };
	int list[4][LEN];

	const __m256d rcut2 = _mm256_set1_pd(rcut2_sd);

	for(int i=0; i<n; i+=4){
		__m256d xi = _mm256_load_pd(&body[i+0].x);
		__m256d yi = _mm256_load_pd(&body[i+1].x);
		__m256d zi = _mm256_load_pd(&body[i+2].x);
		__m256d mi = _mm256_load_pd(&body[i+3].x);
		transpose_4ymm_pd(xi, yi, zi, mi);

		__m256d ax, ay, az, pot;
		ax = ay = az = pot = _mm256_set1_pd(0);

		int n0, n1, n2, n3;
		n0 = n1 = n2 = n3 = 0;;

		for(int j=0; j<n; j+=1){
			__m256d dx = _mm256_sub_pd(xi, _mm256_set1_pd(body[j+0].x));
			__m256d dy = _mm256_sub_pd(yi, _mm256_set1_pd(body[j+0].y));
			__m256d dz = _mm256_sub_pd(zi, _mm256_set1_pd(body[j+0].z));
			__m256d mj = _mm256_set1_pd(body[j+0].m);

			__m256d r2  = 
				_mm256_fmadd_pd(dz, dz, 
					_mm256_fmadd_pd(dy, dy,
						_mm256_mul_pd(dx, dx)));

		       __m256d mask = _mm256_cmp_pd(r2, rcut2, _CMP_LT_OQ);

		       __m256d mri3 = rsqrtCubed_x5(r2, mj);

		       mri3 = _mm256_andnot_pd(mask, mri3);

		       ax  = _mm256_fnmadd_pd(mri3, dx, ax);
		       ay  = _mm256_fnmadd_pd(mri3, dy, ay);
		       az  = _mm256_fnmadd_pd(mri3, dz, az);
		       pot = _mm256_fnmadd_pd(mri3, r2, pot);

		       int m = _mm256_movemask_pd(mask);
		       if(m){
			       if(m&1){
				       list[0][n0++%LEN] = j;
			       }
			       if(m&2){
				       list[1][n1++%LEN] = j;
			       }
			       if(m&4){
				       list[2][n2++%LEN] = j;
			       }
			       if(m&8){
				       list[3][n3++%LEN] = j;
			       }
		       }
		}

		transpose_4ymm_pd(ax, ay, az, pot);
		_mm256_store_pd(&acc[i+0].ax, ax);
		_mm256_store_pd(&acc[i+1].ax, ay);
		_mm256_store_pd(&acc[i+2].ax, az);
		_mm256_store_pd(&acc[i+3].ax, pot);

		int ns[4] = {n0, n1, n2, n3};
		for(int ii=0; ii<4; ii++){
			int num = ns[ii];
			assert(num <= LEN);

			double xi = body[i+ii].x;
			double yi = body[i+ii].y;
			double zi = body[i+ii].z;

			double ax=0, ay=0, az=0, pot=0;

			for(int jj=0; jj<num; jj++){
				int j = list[ii][jj];

				double dx = body[j].x - xi;
				double dy = body[j].y - yi;
				double dz = body[j].z - zi;

				double r2 = dx*dx + dy*dy + dz*dz;

				if(0.0 == r2) continue;

				double ri = 1.0 / sqrt(r2);

				double mri = body[j].m * ri;
				double ri2 = ri * ri;

				double mri3 = mri * ri2;
				double u = r2 * ri * epsinv;

				mri3 *= acut(u);
				mri  *= pcut(u);

				ax += mri3 * dx;
				ay += mri3 * dy;
				az += mri3 * dz;
				pot -= mri;
			}
			acc[i+ii].ax  += ax;
			acc[i+ii].ay  += ay;
			acc[i+ii].az  += az;
			acc[i+ii].pot += pot;
		}
	}
}



#if 0
void prlong(long l){
	l <<= 32;
	double d;
	memcpy(&d, &l, 8);
	printf("%lx : %f\n", l>>32, d);
}

int main(){
	enum { N = 20, };

	prlong(1072693248);
	prlong(1076101120);
	prlong(1077149696);
	prlong(1073741824);
	prlong(-1067941888);
	prlong(1080033280);

	for(int i=0; i<N+3; i++){
		double u = i * (2.0/N);
		double f = pcut(u);
		double g = acut(u);

		printf("%f\t%f\t%f\t%A\n", u, f, f*30, f*30);
		printf("%f\t%f\t%f\t%A\n", u, g, g*30, g*30);
	}
}
#endif

#if 1
int main(){
	enum{
		N = 2048,
	};

	const double eps = 0.05;
	const double epsinv = 1.0 / eps;
	const double rcut2 = 4.0 * eps * eps;

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

	nbody_spline_ref(N, rcut2, epsinv, dbl_body, dbl_acc);

	auto verify = [=](auto kernel){
		kernel(N, rcut2, epsinv, body, acc);

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

	verify(nbody_spline_list);
	verify(nbody_spline_m256d);

	auto warmup = [=](auto kernel, int ntimes=100){
		for(int j=0; j<ntimes; j++){
			kernel(N, rcut2, epsinv, body, acc);
		}
	};

	warmup(nbody_spline_list, 5);

	auto benchmark = [=](auto kernel, int ntimes=10){
		double nsecs[ntimes];
		for(int j=0; j<ntimes; j++){
			for(int i=0; i<N; i++){
				acc[i] = {0,0,0};
			}

			auto tick0 = get_utime();
			kernel(N, rcut2, epsinv, body, acc);
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
	
	benchmark(nbody_spline_ref);
	benchmark(nbody_spline_list);
	benchmark(nbody_spline_m256d);

	return 0;
}
#endif
