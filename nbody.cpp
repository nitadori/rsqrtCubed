#include <x86intrin.h>

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<algorithm>

#include "timer.hpp"

struct Body{
	float x, y, z, m;
};

struct Acceleration{
	float ax, ay, az;
};

struct BodyD{
	double x, y, z, m;
};

struct AccelerationD{
	double ax, ay, az;
};



__attribute__((noinline))
void nbody_ref(
	const int n,
	const double eps2,
	const BodyD body[],
	AccelerationD acc[])
{
	for(int i=0; i<n; i++){ 
		const double xi=body[i].x, yi=body[i].y, zi=body[i].z;
		double ax=0, ay=0, az=0;

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
		}
		acc[i] = {ax, ay, az};
	}
	return ;
}

__attribute__((noinline))
void nbody_ipar(
	const int n,
	const float eps2,
	const Body body[],
	Acceleration acc[])
{
#pragma omp simd
	for(int i=0; i<n; i++){ 
		const float xi=body[i].x, yi=body[i].y, zi=body[i].z;
		float ax=0, ay=0, az=0;

		for(int j=0; j<n; j++){
			float dx = xi - body[j].x;
			float dy = yi - body[j].y;
			float dz = zi - body[j].z;

			float r2 = eps2 + dx*dx;
			r2 += dy*dy;
			r2 += dz*dz;

			float ri = 1.f / sqrtf(r2);

			float mri = body[j].m * ri;
			float ri2 = ri * ri;

			float mri3 = mri * ri2;

			ax -= mri3 * dx;
			ay -= mri3 * dy;
			az -= mri3 * dz;
		}
		acc[i] = {ax, ay, az};
	}
	return ;
}

__attribute__((noinline))
void nbody_jpar(
	const int n,
	const float eps2,
	const Body body[],
	Acceleration acc[])
{
	for(int i=0; i<n; i++){ 
		const float xi=body[i].x, yi=body[i].y, zi=body[i].z;
		float ax=0, ay=0, az=0;

#pragma omp simd
		for(int j=0; j<n; j++){
			float dx = xi - body[j].x;
			float dy = yi - body[j].y;
			float dz = zi - body[j].z;

			float r2 = eps2 + dx*dx;
			r2 += dy*dy;
			r2 += dz*dz;

			float ri = 1.f / sqrtf(r2);

			float mri = body[j].m * ri;
			float ri2 = ri * ri;

			float mri3 = mri * ri2;

			ax -= mri3 * dx;
			ay -= mri3 * dy;
			az -= mri3 * dz;
		}
		acc[i] = {ax, ay, az};
	}
	return ;
}

static void transpose_4AoStoSoA(
		const void *base,
		__m512 &x,
		__m512 &y,
		__m512 &z,
		__m512 &w)
{
	const float *p = (const float *)base;

	__m512 xyzw_00_03 = _mm512_loadu_ps(p +  0);
	__m512 xyzw_04_07 = _mm512_loadu_ps(p + 16);
	__m512 xyzw_08_11 = _mm512_loadu_ps(p + 32);
	__m512 xyzw_12_15 = _mm512_loadu_ps(p + 48);

	__m512 x0x4y0y4_ = _mm512_unpacklo_ps(xyzw_00_03, xyzw_04_07);
	__m512 x8xcy8yc_ = _mm512_unpacklo_ps(xyzw_08_11, xyzw_12_15);
	__m512 z0z4w0w4_ = _mm512_unpackhi_ps(xyzw_00_03, xyzw_04_07);
	__m512 z8zcw8wc_ = _mm512_unpackhi_ps(xyzw_08_11, xyzw_12_15);

	__m512i ilo = _mm512_setr_epi32(0,4,8,12,   1,5,9,13,   16,20,24,28,  17,21,25,29);
	__m512i ihi = _mm512_setr_epi32(2,6,10,14,  3,7,11,15,  18,22,26,30,  19,23,27,31);

	x = _mm512_permutex2var_ps(x0x4y0y4_, ilo, x8xcy8yc_);
	y = _mm512_permutex2var_ps(x0x4y0y4_, ihi, x8xcy8yc_);
	z = _mm512_permutex2var_ps(z0z4w0w4_, ilo, z8zcw8wc_);
	w = _mm512_permutex2var_ps(z0z4w0w4_, ihi, z8zcw8wc_);
}

static void transpose_3SoAtoAoS(
		const __m512 xv,
		const __m512 yv,
		const __m512 zv,
		void *base)
{
	float *p = (float *)base;

	__m512 xylow = _mm512_shuffle_f32x4 (xv, yv, 0b01000100); // 1010(4)
	__m512 xymid = _mm512_shuffle_f32x4 (xv, yv, 0b10011001); // 2121(4)
	__m512 xyhig = _mm512_shuffle_f32x4 (xv, yv, 0b11101110); // 3232(4)

	//                              | x, y, z | x, y, z | x, y, z | x, y, z | x, y, z | x, y, z |
	__m512i ilow = _mm512_setr_epi32( 0, 8,16,  1, 9,17,  2,10,18,  3,11,19,  4,12,20,  5      );
	__m512i imid = _mm512_setr_epi32(    9,21,  2,10,22,  3,11,23,  4,12,24,  5,13,25,  6,14   );
	__m512i ihig = _mm512_setr_epi32(      26,  3,11,27,  4,12,28,  5,13,29,  6,14,30,  7,15,31);

	__m512 xyzlow = _mm512_permutex2var_ps(xylow, ilow, zv);
	__m512 xyzmid = _mm512_permutex2var_ps(xymid, imid, zv);
	__m512 xyzhig = _mm512_permutex2var_ps(xyhig, ihig, zv);

	_mm512_storeu_ps(p+ 0, xyzlow);
	_mm512_storeu_ps(p+16, xyzmid);
	_mm512_storeu_ps(p+32, xyzhig);
}

static inline __m512 rsqrtCubed(
		const __m512 x,
		const __m512 m)
{
	const __m512 one  = _mm512_set1_ps(1.0);
	const __m512 a    = _mm512_set1_ps(3./2.);
	const __m512 b    = _mm512_set1_ps(15./8.);

	__m512 y   = _mm512_rsqrt14_ps(x);

	__m512 y2  = _mm512_mul_ps(y, y);
	__m512 my  = _mm512_mul_ps(m, y);

	__m512 z   = _mm512_mul_ps(my, y2);
	__m512 h   = _mm512_fnmadd_ps(x,  y2, one); // c - a * b
	__m512 abh = _mm512_fmadd_ps(b, h, a); // a + b*h

	__m512 zh  = _mm512_mul_ps(z, h);

	abh = a;  // force 2nd order

	__m512 z1  = _mm512_fmadd_ps(zh, abh, z);

	return z1;
}

__attribute__((noinline))
void nbody_zmm(
	const int n,
	const float eps2_ss,
	const Body body[],
	Acceleration acc[])
{
	const __m512 eps2 = _mm512_set1_ps(eps2_ss);

	for(int i=0; i<n; i+=16){
		__m512 xi, yi, zi, mi;
		transpose_4AoStoSoA(body+i, xi, yi, zi, mi);

		__m512 ax, ay, az;
		ax = ay = az = _mm512_set1_ps(0);

		for(int j=0; j<n; j++){
			__m512 xj = _mm512_set1_ps(body[j].x);
			__m512 yj = _mm512_set1_ps(body[j].y);
			__m512 zj = _mm512_set1_ps(body[j].z);
			__m512 mj = _mm512_set1_ps(body[j].m);

			__m512 dx = _mm512_sub_ps(xi, xj);
			__m512 dy = _mm512_sub_ps(yi, yj);
			__m512 dz = _mm512_sub_ps(zi, zj);

			__m512 r2 = _mm512_fmadd_ps(dx, dx, eps2);
			       r2 = _mm512_fmadd_ps(dy, dy, r2);
			       r2 = _mm512_fmadd_ps(dz, dz, r2);

		       __m512 mri3 = rsqrtCubed(r2, mj);

		       ax = _mm512_fnmadd_ps(mri3, dx, ax);
		       ay = _mm512_fnmadd_ps(mri3, dy, ay);
		       az = _mm512_fnmadd_ps(mri3, dz, az);

		}

		transpose_3SoAtoAoS(ax, ay, az, acc+i);
	}
}

__attribute__((noinline))
void nbody_zmmomp(
	const int n,
	const float eps2_ss,
	const Body body[],
	Acceleration acc[])
{
	const __m512 eps2 = _mm512_set1_ps(eps2_ss);

// KMP_AFFINITY=granularity=fine,compact
// KMP_NUM_THREADS=2
#pragma omp parallel for
	for(int i=0; i<n; i+=16){
		__m512 xi, yi, zi, mi;
		transpose_4AoStoSoA(body+i, xi, yi, zi, mi);

		__m512 ax, ay, az;
		ax = ay = az = _mm512_set1_ps(0);

		for(int j=0; j<n; j++){
			__m512 xj = _mm512_set1_ps(body[j].x);
			__m512 yj = _mm512_set1_ps(body[j].y);
			__m512 zj = _mm512_set1_ps(body[j].z);
			__m512 mj = _mm512_set1_ps(body[j].m);

			__m512 dx = _mm512_sub_ps(xi, xj);
			__m512 dy = _mm512_sub_ps(yi, yj);
			__m512 dz = _mm512_sub_ps(zi, zj);

			__m512 r2 = _mm512_fmadd_ps(dx, dx, eps2);
			       r2 = _mm512_fmadd_ps(dy, dy, r2);
			       r2 = _mm512_fmadd_ps(dz, dz, r2);

		       __m512 mri3 = rsqrtCubed(r2, mj);

		       ax = _mm512_fnmadd_ps(mri3, dx, ax);
		       ay = _mm512_fnmadd_ps(mri3, dy, ay);
		       az = _mm512_fnmadd_ps(mri3, dz, az);

		}

		transpose_3SoAtoAoS(ax, ay, az, acc+i);
	}
}

int main(){
	enum{
		N = 2048,
	};

	const float eps = 1./256.;
	const float eps2 = eps*eps;

	static Body          body    [N];
	static Acceleration  acc     [N];
	static BodyD         dbl_body[N];
	static AccelerationD dbl_acc [N];

	srand48(20210309);
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
		for(int i=0; i<N; i++){
			double e = rel_err(acc[i], dbl_acc[i]);
			// err[i] = e;
			// printf("%e %e %e %e\n", x[i], y0[i], y1[i], err[i]);
			err_max = std::max(err_max, e);
			err_min = std::min(err_min, e);

			if(!std::isfinite(e) || fabs(e) > 1.e-5){
				printf("%4d %+e (%e %e %e)\n", i, e, acc[i].ax, acc[i].ay, acc[i].az);
			}
		}
		printf("err in [%e, %e]\n", err_min, err_max);
		puts("");
	};

	verify(nbody_ipar);
	verify(nbody_jpar);
	verify(nbody_zmm);
	verify(nbody_zmmomp);

	auto benchmark = [=](auto kernel, int ntimes=32){
		double nsecs[ntimes];
		for(int j=0; j<ntimes; j++){
			for(int i=0; i<N; i++){
				acc[i] = {0,0,0};
			}

			auto tick0 = get_utime();
			kernel(N, eps2, body, acc);
			auto tick1 = get_utime();

			double dt = tick2second(tick1 - tick0);
			double iter = N/16.0 * N;
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
	
	benchmark(nbody_ipar);
	benchmark(nbody_ipar);
	benchmark(nbody_zmm);
	benchmark(nbody_zmmomp);

	return 0;
}
