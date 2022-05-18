#include <cstdio>
#include <cmath>
#include <cstring>
#include <algorithm>

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
	const double epsinv,
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
		acc[i] = {ax, ay, az, pot};
	}
	return ;
}


#if 1
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
