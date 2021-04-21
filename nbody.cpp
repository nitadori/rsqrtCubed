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
void nbody_base(
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
			float dx = body[j].x - xi;
			float dy = body[j].y - yi;
			float dz = body[j].z - zi;

			float r2 = eps2 + dx*dx;
			r2 += dy*dy;
			r2 += dz*dz;

			float ri = 1.f / sqrtf(r2);

			float mri = body[j].m * ri;
			float ri2 = ri * ri;

			float mri3 = mri * ri2;

			ax += mri3 * dx;
			ay += mri3 * dy;
			az += mri3 * dz;
		}
		acc[i] = {ax, ay, az};
	}
	return ;
}

static void transpose_4AoStoSoA(
		void *base,
		__m512 &x,
		__m512 &y,
		__m512 &z,
		__m512 &w)
{
	float *p = (float *)base;

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

__attribute__((noinline))
void nbody_zmm(
	const int n,
	const float eps2,
	const Body body[],
	Acceleration acc[])
{
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
	};

	verify(nbody_base);

	auto benchmark = [=](auto kernel, int ntimes=20){
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
	
	benchmark(nbody_base);

	return 0;
}
