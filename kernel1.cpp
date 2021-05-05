#include <arm_sve.h>

struct Body{
	float x, y, z, m;
};

struct Acceleration{
	float ax, ay, az;
};


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
	// svfloat32_t my  = svmul_f32_x(p0, y, m);

	svfloat32_t z   = svmul_f32_x(p0, my, y2);
	svfloat32_t h   = svmsb_f32_x(p0, x,  y2, one);

	svfloat32_t zh  = svmul_f32_x(p0, z, h);
	svfloat32_t abh = svmad_f32_x(p0, b, h, a); // a + b*h

	svfloat32_t z1 = svmad_f32_x(p0, zh, abh, z);

	// 8 ops (3 fma, 3 mul, 1 rsqr)

	return z1;
}

extern void nbody_ext_inner(
		const int n,
		const svfloat32_t eps2,
		const Body body[],
		Acceleration acc[],
		const svfloat32_t one,
		const svfloat32_t a,
		const svfloat32_t b,
		const Body body2[])
{
	const svbool_t p0 = svptrue_b32();

	for(int i=0; i<n; i+=32){
		svfloat32_t xi_0, yi_0, zi_0;
		svfloat32_t xi_1, yi_1, zi_1;
		svfloat32x4_t ibody_0 = svld4_f32(p0, (const float *)(body+i));
		svfloat32x4_t ibody_1 = svld4_f32(p0, (const float *)(body+(i+16)));
#if !defined(__FUJITSU)
		xi_0 = svget4_f32(ibody_0, 0);
		yi_0 = svget4_f32(ibody_0, 1);
		zi_0 = svget4_f32(ibody_0, 2);
		xi_1 = svget4_f32(ibody_1, 0);
		yi_1 = svget4_f32(ibody_1, 1);
		zi_1 = svget4_f32(ibody_1, 2);
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

// #pragma loop unroll 4
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

			// load again
			xj = svdup_f32(body2[j].x);
			yj = svdup_f32(body2[j].y);
			zj = svdup_f32(body2[j].z);
			mj = svdup_f32(body2[j].m);

			// subtract again
			dx_0 = svsub_f32_x(p0, xj, xi_0);
			dy_0 = svsub_f32_x(p0, yj, yi_0);
			dz_0 = svsub_f32_x(p0, zj, zi_0);

			dx_1 = svsub_f32_x(p0, xj, xi_1);
			dy_1 = svsub_f32_x(p0, yj, yi_1);
			dz_1 = svsub_f32_x(p0, zj, zi_1);

			ax_0 = svmla_f32_x(p0, ax_0, mri3_0, dx_0);
			ay_0 = svmla_f32_x(p0, ay_0, mri3_0, dy_0);
			az_0 = svmla_f32_x(p0, az_0, mri3_0, dz_0);

			ax_1 = svmla_f32_x(p0, ax_1, mri3_1, dx_1);
			ay_1 = svmla_f32_x(p0, ay_1, mri3_1, dy_1);
			az_1 = svmla_f32_x(p0, az_1, mri3_1, dz_1);
		}
#if !defined(__FUJITSU)
		svfloat32x3_t acci_0 = svcreate3_f32(ax_0, ay_0, az_0);
		svfloat32x3_t acci_1 = svcreate3_f32(ax_1, ay_1, az_1);
#else
		svfloat32x3_t acci_0 = {ax_0, ay_0, az_0};
		svfloat32x3_t acci_1 = {ax_1, ay_1, az_1};
#endif
		svst3_f32(p0, (float *)(acc+i), acci_0);
		svst3_f32(p0, (float *)(acc+(i+16)), acci_1);
	}
}
