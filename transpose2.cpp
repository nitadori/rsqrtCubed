#include <cstdio>
#include <cassert>
#include <x86intrin.h>

int main(){
	float AoS[16][ 4];
	float SoA[ 4][16];

	for(int i=0; i<16; i++){
		for(int c=0; c<4; c++){
			SoA[c][i] = -1.0;
			AoS[i][c] = (100*(c+1)) + (i+1);
		}
	}

	asm volatile ("#BEG");

	__m512 xyzw_00_03 = _mm512_loadu_ps(AoS[0] +  0);
	__m512 xyzw_04_07 = _mm512_loadu_ps(AoS[0] + 16);
	__m512 xyzw_08_11 = _mm512_loadu_ps(AoS[0] + 32);
	__m512 xyzw_12_15 = _mm512_loadu_ps(AoS[0] + 48);

	__m512 x0x4y0y4_ = _mm512_unpacklo_ps(xyzw_00_03, xyzw_04_07);
	__m512 x8xcy8yc_ = _mm512_unpacklo_ps(xyzw_08_11, xyzw_12_15);
	__m512 z0z4w0w4_ = _mm512_unpackhi_ps(xyzw_00_03, xyzw_04_07);
	__m512 z8zcw8wc_ = _mm512_unpackhi_ps(xyzw_08_11, xyzw_12_15);

	__m512i ilo = _mm512_setr_epi32(0,4,8,12,   1,5,9,13,   16,20,24,28,  17,21,25,29);
	__m512i ihi = _mm512_setr_epi32(2,6,10,14,  3,7,11,15,  18,22,26,30,  19,23,27,31);

	__m512 xs = _mm512_permutex2var_ps(x0x4y0y4_, ilo, x8xcy8yc_);
	__m512 ys = _mm512_permutex2var_ps(x0x4y0y4_, ihi, x8xcy8yc_);
	__m512 zs = _mm512_permutex2var_ps(z0z4w0w4_, ilo, z8zcw8wc_);
	__m512 ws = _mm512_permutex2var_ps(z0z4w0w4_, ihi, z8zcw8wc_);

	_mm512_storeu_ps(SoA[0], xs);
	_mm512_storeu_ps(SoA[1], ys);
	_mm512_storeu_ps(SoA[2], zs);
	_mm512_storeu_ps(SoA[3], ws);

	asm volatile ("#END");

	for(int i=0; i<16; i++){
		for(int c=0; c<4; c++){
			printf("%3.0f ", SoA[c][i]);
		}
		puts("");
	}
	fflush(stdout);
	for(int i=0; i<16; i++){
		for(int c=0; c<4; c++){
			assert(AoS[i][c] == SoA[c][i]);
		}
	}
}
