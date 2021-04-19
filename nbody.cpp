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


__attribute__((noinline))
void nbody_base(
	const int n,
	const float eps2,
	const Body body[],
	Acceleration acc[])
{
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

int main(){
	enum{
		N = 2048,
	};

	const float eps = 1./256.;
	const float eps2 = eps*eps;

	static Body body[N];
	static Acceleration acc[N];

	srand48(20210309);
	for(int i=0; i<N; i++){
		body[i].x = drand48() - 0.5;
		body[i].y = drand48() - 0.5;
		body[i].z = drand48() - 0.5;
		body[i].m = (1./N) * (drand48() + 0.5);
	}

	nbody_base(N, eps2, body, acc);

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

	return 0;
}
