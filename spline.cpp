#include <cstdio>
#include <cmath>
#include <cstring>
#include <algorithm>

__attribute__((noinline))
static double pcut(double u0){
	using F = decltype(u0);
	double u  = std::min(F(2), u0);
	double u2 = u*u;

	double pl = fma(u, -3, 9);
	pl = fma(pl, u2, -20);
	pl = fma(pl, u2, 42);
	pl *= u;
	
	double pr = fma(u, 4, 2);

	double usub = std::max(F(0), u-1);
	double usub2 = usub  * usub;
	double usub4 = usub2 * usub2;
	double usub5 = usub  * usub4;

	double sum = fma(usub5, pr, pl);

	return F(1./30.) * sum;
} 

__attribute__((noinline))
static double acut(double u0){
	using F = decltype(u0);
	double u  = std::min(F(2), u0);
	double u2 = u*u;
	double u3 = u*u2;

	double pl = fma(u, 15, -36);
	pl = fma(pl, u2, 40);
	pl *= u3;

	double pr = fma(u, 20, 8);
	pr = fma(pr, u, 2);

	double usub = std::max(F(0), u-1);
	double usub2 = usub  * usub;
	double usub4 = usub2 * usub2;

	double sum = fma(usub4, -pr, pl);

	return F(1./30.) * sum;
}

void prlong(long l){
	l <<= 32;
	double d;
	memcpy(&d, &l, 8);
	printf("%x : %f\n", l, d);
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

		printf("%f %f %f %A\n", u, f, f*30, f*30);
		printf("%f %f %f %A\n", u, g, g*30, g*30);
	}
}
