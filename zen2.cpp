#include <x86intrin.h>

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<algorithm>

// #include "timer.hpp"

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

void transpose_4ymm_pd(__m256d &r0, __m256d &r1, __m256d &r2, __m256d &r3){
	__m256d tmp0 = _mm256_unpacklo_pd(r0, r1); // | r12 | r02 | r10 | r00 |
	__m256d tmp1 = _mm256_unpackhi_pd(r0, r1); // | r13 | r03 | r11 | r01 |
	__m256d tmp2 = _mm256_unpacklo_pd(r2, r3); // | r32 | r22 | r30 | r20 |
	__m256d tmp3 = _mm256_unpackhi_pd(r2, r3); // | r33 | r23 | r31 | r21 |
	r0 = _mm256_permute2f128_pd(tmp0, tmp2, (0)+(2<<4));
	r1 = _mm256_permute2f128_pd(tmp1, tmp3, (0)+(2<<4));
	r2 = _mm256_permute2f128_pd(tmp0, tmp2, (1)+(3<<4));
	r3 = _mm256_permute2f128_pd(tmp1, tmp3, (1)+(3<<4));
}

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
