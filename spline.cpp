#include <cstdio>
#include <cstdlib>
#include <cmath>
// #include <cstring>
#include <cassert>
#include <algorithm>

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

	return 0;
}
#endif
