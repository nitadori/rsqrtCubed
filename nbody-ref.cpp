#include <cmath>

struct Body{
	float x, y, z, m;
};

struct Acceleration{
	float ax, ay, az;
};

void nbody(
	const int n,
	const float eps2,
	const Body body[],
	Acceleration acc[])
{
	for(int i=0; i<n; i++){ 
		float xi=body[i].x, yi=body[i].y, zi=body[i].z;
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
