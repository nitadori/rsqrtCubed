#include <cstdio>
#include <cassert>
#include <x86intrin.h>

int main(){
	float SoA[ 3][16];
	float AoS[16][ 3];

	for(int i=0; i<16; i++){
		for(int c=0; c<3; c++){
			SoA[c][i] = (100*(c+1)) + (i+1);
			AoS[i][c] = -1.0;
		}
	}

	asm volatile ("#BEG");

	__m512i xv = _mm512_loadu_epi32(SoA[0]);
	__m512i yv = _mm512_loadu_epi32(SoA[1]);
	__m512i zv = _mm512_loadu_epi32(SoA[2]);

#if 0
	// make [x0, y0, x1, y1,...,x7, y7]
	int ixylo[16] = { 0, 16,  1, 17,  2, 18,  3, 19,  4, 20,  5, 21,  6, 22,  7, 23};

	// make [x8, y8, x1, y1,...,x15, y15]
	int ixyhi[16] = { 8, 24,  9, 25, 10, 26, 11, 27, 12, 28, 13, 29, 14, 30, 15, 31};

	// align to make [x5,y5,...,x12,y12]
	
	// make [x0,y0,z0,...,x5]
	int ixyzlow[16] = {0, 1, 16,    2,3,17,    4,5,18,    6,7,19,    8,9,20,    10};

	// make [y5,z5, ... , x10,y10]
	int ixyzmid[16] = {3,21,    4,5,22,    6,7,23,    8,9,24,    10,11,25,    12,13};

	// make [z10, ... , x15,y15,a15]
	int ixyzhig[16] = {26,    6,7,27,    8,9,28,    10,11,29,    12,13,30,    14,15,31};

	__m512i mxylo = _mm512_loadu_epi32(ixylo);
	__m512i mxyhi = _mm512_loadu_epi32(ixyhi);
	__m512i mxyzlow = _mm512_loadu_epi32(ixyzlow);
	__m512i mxyzmid = _mm512_loadu_epi32(ixyzmid);
	__m512i mxyzhig = _mm512_loadu_epi32(ixyzhig);

	__m512i xylow = _mm512_permutex2var_epi32(xv, mxylo, yv);
	__m512i xyhig = _mm512_permutex2var_epi32(xv, mxyhi, yv);
	__m512i xymid =_mm512_alignr_epi32 (xyhig, xylow, 8);

	__m512i xyzlow = _mm512_permutex2var_epi32(xylow, mxyzlow, zv);
	__m512i xyzmid = _mm512_permutex2var_epi32(xymid, mxyzmid, zv);
	__m512i xyzhig = _mm512_permutex2var_epi32(xyhig, mxyzhig, zv);
#else
	__m512i xylow = _mm512_shuffle_i32x4 (xv, yv, 0b01000100); // 1010(4)
	__m512i xymid = _mm512_shuffle_i32x4 (xv, yv, 0b10011001); // 2121(4)
	__m512i xyhig = _mm512_shuffle_i32x4 (xv, yv, 0b11101110); // 3232(4)

	__m512i ilow = _mm512_setr_epi32( 0, 8,16,  1, 9,17,  2,10,18,  3,11,19,  4,12,20,  5      );
	__m512i imid = _mm512_setr_epi32(    9,21,  2,10,22,  3,11,23,  4,12,24,  5,13,25,  6,14   );
	__m512i ihig = _mm512_setr_epi32(      26,  3,11,27,  4,12,28,  5,13,29,  6,14,30,  7,15,31);

	__m512i xyzlow = _mm512_permutex2var_epi32(xylow, ilow, zv);
	__m512i xyzmid = _mm512_permutex2var_epi32(xymid, imid, zv);
	__m512i xyzhig = _mm512_permutex2var_epi32(xyhig, ihig, zv);
#endif

	_mm512_storeu_epi32(AoS[0]+ 0, xyzlow);
	_mm512_storeu_epi32(AoS[0]+16, xyzmid);
	_mm512_storeu_epi32(AoS[0]+32, xyzhig);

	asm volatile ("#END");

	for(int i=0; i<16; i++){
		for(int c=0; c<3; c++){
			printf("%3.0f ", AoS[i][c]);
		}
		puts("");
	}
	for(int c=0; c<3; c++){
		for(int i=0; i<16; i++){
			printf("%3.0f ", ((float *)AoS[0])[16*c+i]);
			// assert(AoS[i][c] == SoA[i][c]);
		}
		puts("");
	}
	fflush(stdout);
	for(int i=0; i<16; i++){
		for(int c=0; c<3; c++){
			assert(AoS[i][c] == SoA[c][i]);
		}
	}

	return 0;
}
