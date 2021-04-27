#include <arm_sve.h>
#include <cmath>

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

static inline svfloat32_t rsqrtCubed(
		const svfloat32_t x,
		const svfloat32_t m,
		const svbool_t p0,
		const svfloat32_t one,
		const svfloat32_t a,
		const svfloat32_t b)
{
	svfloat32_t y   = svrsqrte_f32(x);

	svfloat32_t y2  = svmul_f32_x(p0, y, y);
	svfloat32_t my  = svmul_f32_x(p0, m, y);

	svfloat32_t z   = svmul_f32_x(p0, my, y2);
	svfloat32_t h   = svmsb_f32_x(p0, x,  y2, one);

	svfloat32_t zh  = svmul_f32_x(p0, z, h);
	svfloat32_t abh = svmad_f32_x(p0, b, h, a); // a + b*h

	svfloat32_t z1 = svmad_f32_x(p0, zh, abh, z);

	// 8 ops

	return z1;
}

__attribute__((noinline))
void nbody_sve(
	const int n,
	const float eps2_ss,
	const Body body[],
	Acceleration acc[])
{
	const svfloat32_t eps2 = svdup_f32(eps2_ss);

	const svfloat32_t one  = svdup_f32(1.0);
	const svfloat32_t b    = svdup_f32(15./8.);
	const svfloat32_t a    = svdup_f32( 3./2.);

	const svbool_t p0 = svptrue_b32();

	for(int i=0; i<n; i+=16){
		svfloat32_t xi, yi, zi;
		// transpose_4AoStoSoA(body+i, xi, yi, zi, mi);
		svfloat32x4_t ibody = svld4_f32(p0, (const float *)(body+i));
#if 0
		xi = svget4_f32(ibody, 0);
		yi = svget4_f32(ibody, 1);
		zi = svget4_f32(ibody, 2);
#else
		xi = ibody.v0;
		yi = ibody.v1;
		zi = ibody.v2;
#endif

		svfloat32_t ax, ay, az;
		ax = ay = az = svdup_f32(0);

		for(int j=0; j<n; j++){
			svfloat32_t xj = svdup_f32(body[j].x);
			svfloat32_t yj = svdup_f32(body[j].y);
			svfloat32_t zj = svdup_f32(body[j].z);
			svfloat32_t mj = svdup_f32(body[j].m);

			svfloat32_t dx = svsub_f32_x(p0, xj, xi);
			svfloat32_t dy = svsub_f32_x(p0, yj, yi);
			svfloat32_t dz = svsub_f32_x(p0, zj, zi);

			svfloat32_t r2 = svmad_f32_x(p0, dx, dx, eps2);
			r2 = svmad_f32_x(p0, dy, dy, r2);
			r2 = svmad_f32_x(p0, dz, dz, r2);


#if 0
			svfloat32_t y   = svrsqrte_f32(r2);

			svfloat32_t y2  = svmul_f32_x(p0, y, y);
			svfloat32_t my  = svmul_f32_x(p0, mj, y);

			svfloat32_t z   = svmul_f32_x(p0, my, y2);
			svfloat32_t h   = svmsb_f32_x(p0, r2,  y2, one);

			svfloat32_t zh  = svmul_f32_x(p0, z, h);
			svfloat32_t abh = svmad_f32_x(p0, b, h, a); // a + b*h

			svfloat32_t mri3 = svmad_f32_x(p0, zh, abh, z);
#else
			svfloat32_t mri3 = rsqrtCubed(r2, mj, p0, one, a, b);
#endif

			ax = svmla_f32_x(p0, ax, mri3, dx);
			ay = svmla_f32_x(p0, ay, mri3, dy);
			az = svmla_f32_x(p0, az, mri3, dz);

			// 17-ops, 8.5-cycle in the best case
		}

		// transpose_3SoAtoAoS(ax, ay, az, acc+i);
#if 0
		svfloat32x3_t acci = svcreate3_f32(ax, ay, az);
#else
		svfloat32x3_t acci = {ax, ay, az};
#endif
		svst3_f32(p0, (float *)(acc+i), acci);
	}
}

__attribute__((noinline))
void nbody_ni32(
	const int n,
	const float eps2_ss,
	const Body body[],
	Acceleration acc[])
{
	const svfloat32_t eps2 = svdup_f32(eps2_ss);

	const svfloat32_t one  = svdup_f32(1.0);
	const svfloat32_t b    = svdup_f32(15./8.);
	const svfloat32_t a    = svdup_f32( 3./2.);

	const svbool_t p0 = svptrue_b32();

	for(int i=0; i<n; i+=32){
		svfloat32_t xi_0, yi_0, zi_0;
		svfloat32_t xi_1, yi_1, zi_1;
		svfloat32x4_t ibody_0 = svld4_f32(p0, (const float *)(body+i));
		svfloat32x4_t ibody_1 = svld4_f32(p0, (const float *)(body+(i+16)));
#if 0
		xi = svget4_f32(ibody, 0);
		yi = svget4_f32(ibody, 1);
		zi = svget4_f32(ibody, 2);
#else
		xi_0 = ibody_0.v0;
		yi_0 = ibody_0.v1;
		zi_0 = ibody_0.v2;
		xi_1 = ibody_1.v0;
		yi_1 = ibody_1.v1;
		zi_1 = ibody_1.v2;
#endif

		svfloat32_t ax_0, ay_0, az_0;
		ax_0 = ay_0 = az_0 = svdup_f32(0);
		svfloat32_t ax_1, ay_1, az_1;
		ax_1 = ay_1 = az_1 = svdup_f32(0);

		for(int j=0; j<n; j++){
			svfloat32_t xj = svdup_f32(body[j].x);
			svfloat32_t yj = svdup_f32(body[j].y);
			svfloat32_t zj = svdup_f32(body[j].z);
			svfloat32_t mj = svdup_f32(body[j].m);

			svfloat32_t dx_0 = svsub_f32_x(p0, xj, xi_0);
			svfloat32_t dy_0 = svsub_f32_x(p0, yj, yi_0);
			svfloat32_t dz_0 = svsub_f32_x(p0, zj, zi_0);

			svfloat32_t dx_1 = svsub_f32_x(p0, xj, xi_1);
			svfloat32_t dy_1 = svsub_f32_x(p0, yj, yi_1);
			svfloat32_t dz_1 = svsub_f32_x(p0, zj, zi_1);

			svfloat32_t r2_0 = svmad_f32_x(p0, dx_0, dx_0, eps2);
			r2_0 = svmad_f32_x(p0, dy_0, dy_0, r2_0);
			r2_0 = svmad_f32_x(p0, dz_0, dz_0, r2_0);

			svfloat32_t r2_1 = svmad_f32_x(p0, dx_1, dx_1, eps2);
			r2_1 = svmad_f32_x(p0, dy_1, dy_1, r2_1);
			r2_1 = svmad_f32_x(p0, dz_1, dz_1, r2_1);


			svfloat32_t mri3_0 = rsqrtCubed(r2_0, mj, p0, one, a, b);
			svfloat32_t mri3_1 = rsqrtCubed(r2_1, mj, p0, one, a, b);

			ax_0 = svmla_f32_x(p0, ax_0, mri3_0, dx_0);
			ay_0 = svmla_f32_x(p0, ay_0, mri3_0, dy_0);
			az_0 = svmla_f32_x(p0, az_0, mri3_0, dz_0);

			ax_1 = svmla_f32_x(p0, ax_1, mri3_1, dx_1);
			ay_1 = svmla_f32_x(p0, ay_1, mri3_1, dy_1);
			az_1 = svmla_f32_x(p0, az_1, mri3_1, dz_1);

			// 17-ops, 8.5-cycle in the best case
		}
#if 0
		svfloat32x3_t acci = svcreate3_f32(ax, ay, az);
#else
		svfloat32x3_t acci_0 = {ax_0, ay_0, az_0};
		svfloat32x3_t acci_1 = {ax_1, ay_1, az_1};
#endif
		svst3_f32(p0, (float *)(acc+i), acci_0);
		svst3_f32(p0, (float *)(acc+(i+16)), acci_1);
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
	verify(nbody_sve);
	verify(nbody_ni32);

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
			double iter = N/16.0 * N;
			nsecs[j] = dt/iter * 1.e9;
		}
		for(int j=0; j<ntimes; j++){
#ifdef __aarch64__ 
			// Just assume 2.0 GHz of Fugaku
			double cycle = nsecs[j] * 2.0;
			double eff = 100.0 * 8.5 / cycle;
			printf("%f nsec/loop, %f cycles, %f%%\n", nsecs[j], cycle, eff);
#else
			printf("%f nsec/loop\n", nsecs[j]);
#endif
		}
		puts("");
		fflush(stdout);
	};
	
	benchmark(nbody_ipar);
	benchmark(nbody_jpar);
	benchmark(nbody_sve);
	benchmark(nbody_ni32);

	return 0;
}
