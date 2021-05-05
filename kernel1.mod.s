..text.b:
	.ident	"$Options: Fujitsu C/C++ Compiler 4.5.0 (Mar  4 2021 13:29:06) --preinclude /opt/FJSVxtclanga/tcsds-1.2.31/bin/../lib/FCC.pre --g++ -Dunix -Dlinux -D__FUJITSU -D__FCC_major__=4 -D__FCC_minor__=5 -D__FCC_patchlevel__=0 -D__FCC_version__=\"4.5.0\" -D__aarch64__ -D__unix -D__PRAGMA_REDEFINE_EXTNAME -D__FCC_VERSION=800 -D__USER_LABEL_PREFIX__= -D__OPTIMIZE__ -D__ARM_ARCH=8 -D__ARM_FEATURE_SVE -D__FP_FAST_FMA -D__ELF__ -D__linux -Asystem(unix) -D__LIBC_6B -D_LP64 -D__LP64__ --K=noocl --zmode=64 --sys_include=/opt/FJSVxtclanga/tcsds-1.2.31/bin/../include/libc++/v371 --sys_include=/opt/FJSVxtclanga/tcsds-1.2.31/bin/../include --sys_include=/opt/FJSVxos/devkit/aarch64/rfs/usr/include --K=opt --exceptions kernel1.cpp -- -ncmdname=FCCpx -zobe=cplus -zcfc=target_sve -O3 -x- -Komitfp,mfunc,eval,fp_relaxed,fz,fast_matmul,fp_contract,ilfunc,simd_packed_promotion -Klargepage kernel1.s $"
	.file	"kernel1.cpp"
	.ident	"$Compiler: Fujitsu C/C++ Compiler 4.5.0 (Mar  4 2021 13:29:06) kernel1.cpp _Z15nbody_ext_inneriDvfPK4BodyP12AccelerationS_S_S_S2_ $"
	.text
	.align	2
	.global	_Z15nbody_ext_inneriDvfPK4BodyP12AccelerationS_S_S_S2_
	.type	_Z15nbody_ext_inneriDvfPK4BodyP12AccelerationS_S_S_S2_, %function
_Z15nbody_ext_inneriDvfPK4BodyP12AccelerationS_S_S_S2_:
	.file 1 "kernel1.cpp"
	.loc 1 39 0
..LDL1:
.LFB0:
	.cfi_startproc
/*    ??? */	addvl	sp, sp, -32
/*    ??? */	addvl	sp, sp, -32
/*    ??? */	addvl	sp, sp, -3
/*    ??? */	stp	x29, x30, [sp]	//  (*)
/*     39 */	add	x29, sp, 0
	.cfi_escape 0xf,0xb,0x92,0x1d,0x0,0x11,0x98,0x4,0x92,0x2e,0x0,0x1e,0x22
	.cfi_escape 0x10,0x1d,0x8,0x11,0xe8,0x7b,0x92,0x2e,0x0,0x1e,0x22
	.cfi_escape 0x10,0x1e,0xb,0x11,0x8,0x22,0x11,0xe8,0x7b,0x92,0x2e,0x0,0x1e,0x22
/*     39 */	sub	sp, sp, 16
/*    ??? */	str	x19, [x29, -8]	//  (*)
	.cfi_escape 0x10,0x13,0xb,0x11,0x78,0x22,0x11,0xe8,0x7b,0x92,0x2e,0x0,0x1e,0x22
/*     39 */	add	x19, sp, 0
/*     39 */	and	sp, x19, -64
/*    ??? */	str	z8, [x29, 54, mul vl]	//  (*)
	.cfi_escape 0x10,0x68,0x8,0x11,0x98,0x7f,0x92,0x2e,0x0,0x1e,0x22
/*    ??? */	str	z9, [x29, 55, mul vl]	//  (*)
	.cfi_escape 0x10,0x69,0x8,0x11,0xa0,0x7f,0x92,0x2e,0x0,0x1e,0x22
/*    ??? */	str	z10, [x29, 56, mul vl]	//  (*)
	.cfi_escape 0x10,0x6a,0x8,0x11,0xa8,0x7f,0x92,0x2e,0x0,0x1e,0x22
/*    ??? */	str	z11, [x29, 57, mul vl]	//  (*)
	.cfi_escape 0x10,0x6b,0x8,0x11,0xb0,0x7f,0x92,0x2e,0x0,0x1e,0x22
/*    ??? */	str	z12, [x29, 58, mul vl]	//  (*)
	.cfi_escape 0x10,0x6c,0x8,0x11,0xb8,0x7f,0x92,0x2e,0x0,0x1e,0x22
/*    ??? */	str	z13, [x29, 59, mul vl]	//  (*)
	.cfi_escape 0x10,0x6d,0x7,0x11,0x40,0x92,0x2e,0x0,0x1e,0x22
/*    ??? */	str	z14, [x29, 60, mul vl]	//  (*)
	.cfi_escape 0x10,0x6e,0x7,0x11,0x48,0x92,0x2e,0x0,0x1e,0x22
/*    ??? */	str	z15, [x29, 61, mul vl]	//  (*)
	.cfi_escape 0x10,0x6f,0x7,0x11,0x50,0x92,0x2e,0x0,0x1e,0x22
/*    ??? */	str	z16, [x29, 62, mul vl]	//  (*)
	.cfi_escape 0x10,0x70,0x7,0x11,0x58,0x92,0x2e,0x0,0x1e,0x22
/*    ??? */	str	z17, [x29, 63, mul vl]	//  (*)
	.cfi_escape 0x10,0x71,0x7,0x11,0x60,0x92,0x2e,0x0,0x1e,0x22
/*    ??? */	str	z18, [x29, 64, mul vl]	//  (*)
	.cfi_escape 0x10,0x72,0x7,0x11,0x68,0x92,0x2e,0x0,0x1e,0x22
/*    ??? */	str	z19, [x29, 65, mul vl]	//  (*)
	.cfi_escape 0x10,0x73,0x7,0x11,0x70,0x92,0x2e,0x0,0x1e,0x22
/*    ??? */	str	z20, [x29, 66, mul vl]	//  (*)
	.cfi_escape 0x10,0x74,0x7,0x11,0x78,0x92,0x2e,0x0,0x1e,0x22
/*    ??? */	str	z3, [x29, 15, mul vl]	//  (*)
/*    ??? */	str	z2, [x29, 16, mul vl]	//  (*)
/*    ??? */	str	z1, [x29, 14, mul vl]	//  (*)
/*    ??? */	str	z0, [x29, 13, mul vl]	//  (*)
	.loc 1 49 0
..LDL2:
/*     49 */	ptrue	p0.s, ALL
	.loc 1 51 0
..LDL3:
/*     51 */	sub	w5, w0, 1
/*     51 */	cmp	w0, 0
/*     51 */	asr	w4, w5, 4
/*     51 */	add	w4, w5, w4, lsr #27
/*     51 */	asr	w4, w4, 5
/*     51 */	add	w13, w4, 1
/*     51 */	ble	.L81
	.loc 1 54 0 is_stmt 0
..LDL4:
/*     54 */	mov	x8, 0
	.loc 1 55 0
..LDL5:
/*     55 */	mov	x9, 256
	.loc 1 131 0
..LDL6:
/*    131 */	addvl	x11, x29, 4
	.loc 1 135 0
..LDL7:
/*    135 */	mov	x12, 192
	.loc 1 79 0
..LDL8:
/*     79 */	add	x10, x1, 12
	.loc 1 132 0
..LDL9:
/*    132 */	addvl	x7, x29, 1
	.loc 1 105 0
..LDL10:
/*    105 */	add	x15, x3, 8
	.loc 1 132 0
..LDL11:
/*    132 */	mov	x14, 0
.L73:					// :entr
	.loc 1 54 0 is_stmt 1
..LDL12:
/*     54 */	mov	x4, x8
	.loc 1 55 0
..LDL13:
/*     55 */	mov	x3, x9
..D1.pchi:
	.loc 1 73 0 is_stmt 0
..LDL14:
	.loc 1 73 0 is_stmt 1
..LDL15:
/*     73 */	fmov	z27.s, 0.000000e+00
	.loc 1 54 0
..LDL16:
/*     54 */	add	x4, x4, x1
	.loc 1 55 0
..LDL17:
/*     55 */	add	x3, x3, x1
	.loc 1 54 0
..LDL18:
/*     54 */	ld4w	{z4.s, z5.s, z6.s, z7.s}, p0/z, [x4, 0, mul vl]	//  (*)
	.loc 1 55 0
..LDL19:
/*     55 */	ld4w	{z0.s, z1.s, z2.s, z3.s}, p0/z, [x3, 0, mul vl]	//  (*)
	.loc 1 54 0
..LDL20:
/*    ??? */	str	z7, [x29, 33, mul vl]	//  (*)
/*    ??? */	str	z6, [x29, 12, mul vl]	//  (*)
/*    ??? */	str	z5, [x29, 11, mul vl]	//  (*)
/*    ??? */	str	z4, [x29, 10, mul vl]	//  (*)
	.loc 1 55 0
..LDL21:
/*    ??? */	str	z3, [x29, 53, mul vl]	//  (*)
/*    ??? */	str	z2, [x29, 36, mul vl]	//  (*)
/*    ??? */	str	z1, [x29, 35, mul vl]	//  (*)
/*    ??? */	str	z0, [x29, 34, mul vl]	//  (*)
	.loc 1 78 0
..LDL22:
/*     78 */	cmp	w0, 0
/*     78 */	mov	z1.d, z27.d
/*     78 */	mov	z5.d, z27.d
/*     78 */	mov	z16.d, z27.d
/*     78 */	mov	z3.d, z27.d
/*     78 */	mov	z13.d, z27.d
/*     78 */	ble	.L78
	.loc 1 105 0 is_stmt 0
..LDL23:
/*    105 */	mov	z1.d, z27.d
/*    105 */	mov	z5.d, z27.d
/*    105 */	mov	w6, w0
/*    105 */	mov	x4, x10
/*    105 */	mov	z16.d, z27.d
/*    105 */	mov	z3.d, z27.d
/*    105 */	mov	x3, x15
	.loc 1 78 0
..LDL24:
/*     78 */	cmp	w6, 6
	.loc 1 105 0
..LDL25:
/*    105 */	mov	z13.d, z27.d
	.loc 1 78 0
..LDL26:
/*     78 */	blt	.L92
	.loc 1 79 0
..LDL27:
/*     79 */	ldr	s0, [x4, -12]	//  (*)
/*    ??? */	ldr	z4, [x29, 34, mul vl]	//  (*)
/*    126 */	add	x5, x4, 16
/*    126 */	sub	w6, w6, 2
/*    ??? */	ldr	z2, [x29, 10, mul vl]	//  (*)
/*    126 */	sub	w6, w6, 2
/*     81 */	ldr	s7, [x4, -4]	//  (*)
/*     79 */	ldr	s30, [x4, 4]	//  (*)
/*     82 */	ld1rw	{z25.s}, p0/z, [x4]	//  (*)
/*     81 */	ldr	s28, [x4, 12]	//  (*)
/*     79 */	dup	z0.s, z0.s[0]
/*     84 */	fsub	z2.s, z0.s, z2.s
/*     88 */	fsub	z31.s, z0.s, z4.s
/*    ??? */	ldr	z4, [x29, 13, mul vl]	//  (*)
/*     80 */	ldr	s0, [x4, -8]	//  (*)
/*     80 */	dup	z6.s, z0.s[0]
/*    ??? */	ldr	z0, [x29, 35, mul vl]	//  (*)
/*     96 */	fmad	z31.s, p0/m, z31.s, z4.s
/*     92 */	movprfx	z14.s, p0/z, z2.s
/*     92 */	fmad	z14.s, p0/m, z2.s, z4.s
/*    ??? */	ldr	z4, [x29, 11, mul vl]	//  (*)
/*     81 */	dup	z2.s, z7.s[0]
/*     89 */	fsub	z0.s, z6.s, z0.s
/*     85 */	fsub	z26.s, z6.s, z4.s
/*    ??? */	ldr	z4, [x29, 12, mul vl]	//  (*)
/*    ??? */	ldr	z6, [x29, 13, mul vl]	//  (*)
/*     97 */	fmla	z31.s, p0/m, z0.s, z0.s
/*     79 */	dup	z0.s, z30.s[0]
/*     93 */	fmad	z26.s, p0/m, z26.s, z14.s
/*     81 */	dup	z30.s, z28.s[0]
/*     86 */	fsub	z7.s, z2.s, z4.s
/*    ??? */	ldr	z4, [x29, 36, mul vl]	//  (*)
/*     94 */	fmla	z26.s, p0/m, z7.s, z7.s
/*     90 */	fsub	z10.s, z2.s, z4.s
/*    ??? */	ldr	z2, [x29, 10, mul vl]	//  (*)
/*    ??? */	ldr	z4, [x29, 34, mul vl]	//  (*)
/*     98 */	fmad	z10.s, p0/m, z10.s, z31.s
/*    ??? */	ldr	z31, [x29, 36, mul vl]	//  (*)
/*     84 */	fsub	z2.s, z0.s, z2.s
/*     88 */	fsub	z4.s, z0.s, z4.s
/*     80 */	ldr	s0, [x4, 8]	//  (*)
/*    126 */	add	x4, x4, 32
/*     81 */	ldr	s14, [x4, -4]	//  (*)
/*     20 */	frsqrte	z24.s, z10.s
/*     80 */	dup	z7.s, z0.s[0]
/*    ??? */	ldr	z0, [x29, 35, mul vl]	//  (*)
/*     92 */	fmad	z2.s, p0/m, z2.s, z6.s
/*     96 */	fmad	z4.s, p0/m, z4.s, z6.s
/*    ??? */	ldr	z6, [x29, 11, mul vl]	//  (*)
/*     90 */	fsub	z11.s, z30.s, z31.s
/*     22 */	fmul	z29.s, z24.s, z24.s
/*     23 */	fmul	z28.s, z25.s, z24.s
/*    ??? */	ldr	z24, [x29, 12, mul vl]	//  (*)
/*     82 */	ld1rw	{z31.s}, p0/z, [x5]	//  (*)
/*     81 */	dup	z18.s, z14.s[0]
/*    126 */	add	x5, x4, 16
/*     89 */	fsub	z0.s, z7.s, z0.s
/*     85 */	fsub	z7.s, z7.s, z6.s
/*     20 */	frsqrte	z6.s, z26.s
/*     86 */	fsub	z24.s, z30.s, z24.s
/*     79 */	ldr	s30, [x4, -12]	//  (*)
/*     23 */	fmul	z25.s, z25.s, z6.s
/*     22 */	fmul	z6.s, z6.s, z6.s
/*     97 */	fmad	z0.s, p0/m, z0.s, z4.s
/*    ??? */	ldr	z4, [x29, 14, mul vl]	//  (*)
/*     93 */	fmad	z7.s, p0/m, z7.s, z2.s
/*    ??? */	ldr	z2, [x29, 14, mul vl]	//  (*)
/*     79 */	dup	z30.s, z30.s[0]
/*     98 */	fmad	z11.s, p0/m, z11.s, z0.s
/*    105 */	ldr	s0, [x3, -8]	//  (*)
/*     94 */	fmad	z24.s, p0/m, z24.s, z7.s
/*    ??? */	ldr	z7, [x29, 15, mul vl]	//  (*)
/*     27 */	fmsb	z10.s, p0/m, z29.s, z4.s
/*     26 */	fmul	z4.s, z28.s, z29.s
/*    107 */	ld1rw	{z29.s}, p0/z, [x3]	//  (*)
/*     27 */	fmsb	z26.s, p0/m, z6.s, z2.s
/*     26 */	fmul	z2.s, z25.s, z6.s
/*    ??? */	ldr	z6, [x29, 15, mul vl]	//  (*)
/*    ??? */	ldr	z25, [x29, 16, mul vl]	//  (*)
/*    105 */	dup	z28.s, z0.s[0]
/*    ??? */	ldr	z0, [x29, 16, mul vl]	//  (*)
/*     30 */	fmad	z7.s, p0/m, z26.s, z25.s
/*    106 */	ldr	s25, [x3, -4]	//  (*)
/*     29 */	fmul	z26.s, z26.s, z2.s
/*     30 */	fmad	z6.s, p0/m, z10.s, z0.s
/*     29 */	fmul	z0.s, z10.s, z4.s
/*     20 */	frsqrte	z10.s, z24.s
/*    106 */	dup	z9.s, z25.s[0]
/*    ??? */	ldr	z25, [x29, 34, mul vl]	//  (*)
/*     32 */	fmla	z2.s, p0/m, z7.s, z26.s
/*     23 */	fmul	z15.s, z31.s, z10.s
/*    ??? */	ldr	z26, [x29, 11, mul vl]	//  (*)
/*    ??? */	ldr	z7, [x29, 12, mul vl]	//  (*)
/*     32 */	fmla	z4.s, p0/m, z6.s, z0.s
/*    ??? */	ldr	z0, [x29, 10, mul vl]	//  (*)
/*    115 */	fsub	z12.s, z28.s, z25.s
/*     20 */	frsqrte	z25.s, z11.s
/*     86 */	fsub	z7.s, z18.s, z7.s
/*     22 */	fmul	z8.s, z25.s, z25.s
/*     23 */	fmul	z17.s, z31.s, z25.s
/*    111 */	fsub	z6.s, z28.s, z0.s
/*    ??? */	ldr	z0, [x29, 12, mul vl]	//  (*)
/*    123 */	fmad	z12.s, p0/m, z4.s, z13.s
/*    113 */	fsub	z28.s, z29.s, z0.s
/*    ??? */	ldr	z0, [x29, 36, mul vl]	//  (*)
/*    117 */	fsub	z29.s, z29.s, z0.s
/*    ??? */	ldr	z0, [x29, 35, mul vl]	//  (*)
/*    125 */	fmad	z29.s, p0/m, z4.s, z16.s
/*    ??? */	ldr	z16, [x29, 14, mul vl]	//  (*)
/*    116 */	fsub	z25.s, z9.s, z0.s
/*    ??? */	ldr	z0, [x29, 11, mul vl]	//  (*)
/*     27 */	fmsb	z11.s, p0/m, z8.s, z16.s
/*     26 */	fmul	z8.s, z17.s, z8.s
/*    112 */	fsub	z0.s, z9.s, z0.s
/*     80 */	ldr	s9, [x4, -8]	//  (*)
/*     80 */	dup	z31.s, z9.s[0]
/*    ??? */	ldr	z9, [x29, 35, mul vl]	//  (*)
/*     85 */	fsub	z26.s, z31.s, z26.s
/*     89 */	fsub	z31.s, z31.s, z9.s
/*    ??? */	ldr	z9, [x29, 34, mul vl]	//  (*)
/*     88 */	fsub	z14.s, z30.s, z9.s
/*    ??? */	ldr	z9, [x29, 13, mul vl]	//  (*)
/*     96 */	fmad	z14.s, p0/m, z14.s, z9.s
/*     22 */	fmul	z9.s, z10.s, z10.s
/*    ??? */	ldr	z10, [x29, 10, mul vl]	//  (*)
/*     97 */	fmad	z31.s, p0/m, z31.s, z14.s
/*    ??? */	ldr	z14, [x29, 13, mul vl]	//  (*)
/*     84 */	fsub	z30.s, z30.s, z10.s
/*    ??? */	ldr	z10, [x29, 36, mul vl]	//  (*)
/*     92 */	fmla	z14.s, p0/m, z30.s, z30.s
/*     90 */	fsub	z10.s, z18.s, z10.s
	.p2align 5
.L76:					// :entr:term:swpl
// X5 == x4+16
/*     93 */	fmad	z26.s, p0/m, z26.s, z14.s
/*    126 */	sub	w6, w6, 2
/*    ??? */	ld1rw	{z13.s}, p0/z, [x4, 4]	//  (*)
/* #00001 */	ldr	z30, [x29, 14, mul vl]	//  (*)
/*    101 */	cmp	w6, 2
/*    124 */	fmad	z4.s, p0/m, z25.s, z3.s
/*    121 */	fmad	z28.s, p0/m, z2.s, z27.s
/*    120 */	fmad	z0.s, p0/m, z2.s, z1.s
/*     27 */	fmls	z30.s, p0/m, z9.s, z24.s
/*     26 */	fmul	z15.s, z15.s, z9.s
/*    119 */	fmad	z6.s, p0/m, z2.s, z5.s
/* #00001 */	ldr	z2, [x29, 10, mul vl]	//  (*)
/*    126 */	add	x16, x3, 16
/* #00001 */	ldr	z5, [x29, 16, mul vl]	//  (*)
/*     98 */	fmad	z10.s, p0/m, z10.s, z31.s
/*     84 */	fsub	z9.s, z13.s, z2.s
/* #00001 */	ldr	z2, [x29, 34, mul vl]	//  (*)
/*    107 */	ld1rw	{z25.s}, p0/z, [x16]	//  (*)
/*     88 */	fsub	z14.s, z13.s, z2.s
/*    ??? */	ld1rw	{z3.s}, p0/z, [x3, 8]	//  (*)
/* #00001 */	ldr	z1, [x29, 15, mul vl]	//  (*)
/*     82 */	ld1rw	{z16.s}, p0/z, [x4]	//  (*)
/*     29 */	fmul	z2.s, z11.s, z8.s
/*     30 */	fmad	z1.s, p0/m, z11.s, z5.s
/* #00001 */	ldr	z5, [x29, 15, mul vl]	//  (*)
/*     94 */	fmla	z26.s, p0/m, z7.s, z7.s
/* #00001 */	ldr	z11, [x29, 16, mul vl]	//  (*)
/*    ??? */	ld1rw	{z27.s}, p0/z, [x4, 8]	//  (*)
/*    ??? */	ld1rw	{z7.s}, p0/z, [x3, 12]	//  (*)
/* #00001 */	ldr	z24, [x29, 34, mul vl]	//  (*)
/*     30 */	fmla	z11.s, p0/m, z30.s, z5.s
/*     29 */	fmul	z5.s, z30.s, z15.s
/*    115 */	fsub	z13.s, z3.s, z24.s
/* #00001 */	ldr	z24, [x29, 35, mul vl]	//  (*)
/*     20 */	frsqrte	z18.s, z10.s
/*     89 */	fsub	z31.s, z27.s, z24.s
/* #00001 */	ldr	z24, [x29, 13, mul vl]	//  (*)
/*     96 */	fmad	z14.s, p0/m, z14.s, z24.s
/*     92 */	fmad	z9.s, p0/m, z9.s, z24.s
/* #00001 */	ldr	z24, [x29, 36, mul vl]	//  (*)
/*    ??? */	ld1rw	{z19.s}, p0/z, [x4, 12]	//  (*)
/*    117 */	fsub	z30.s, z25.s, z24.s
/*     32 */	fmad	z2.s, p0/m, z1.s, z8.s
/* #00001 */	ldr	z1, [x29, 11, mul vl]	//  (*)
/*     22 */	fmul	z17.s, z18.s, z18.s
/*     23 */	fmul	z18.s, z16.s, z18.s
/*     85 */	fsub	z24.s, z27.s, z1.s
/* #00001 */	ldr	z1, [x29, 12, mul vl]	//  (*)
/*     20 */	frsqrte	z20.s, z26.s
/*    113 */	fsub	z27.s, z25.s, z1.s
/* #00001 */	ldr	z1, [x29, 10, mul vl]	//  (*)
/*    111 */	fsub	z25.s, z3.s, z1.s
/* #00001 */	ldr	z1, [x29, 35, mul vl]	//  (*)
/*    116 */	fsub	z3.s, z7.s, z1.s
/* #00001 */	ldr	z1, [x29, 11, mul vl]	//  (*)
/*    112 */	fsub	z1.s, z7.s, z1.s
/* #00001 */	ldr	z7, [x29, 12, mul vl]	//  (*)
/*     86 */	fsub	z7.s, z19.s, z7.s
/*     32 */	fmad	z5.s, p0/m, z11.s, z15.s
/* #00001 */	ldr	z11, [x29, 36, mul vl]	//  (*)
/*     23 */	fmul	z8.s, z16.s, z20.s
/*     22 */	fmul	z16.s, z20.s, z20.s
/*     90 */	fsub	z11.s, z19.s, z11.s
/*     97 */	fmad	z31.s, p0/m, z31.s, z14.s
/*    125 */	fmad	z30.s, p0/m, z2.s, z29.s
/* #00001 */	ldr	z29, [x29, 14, mul vl]	//  (*)
/*    123 */	fmad	z13.s, p0/m, z2.s, z12.s
/*    126 */	add	x4, x5, 16
// x4 == x5+16
/*     27 */	fmsb	z10.s, p0/m, z17.s, z29.s
/*     26 */	fmul	z15.s, z18.s, z17.s
/*     93 */	fmad	z24.s, p0/m, z24.s, z9.s
/*    ??? */	ld1rw	{z29.s}, p0/z, [x5, 4]	//  (*)
/*    124 */	fmad	z3.s, p0/m, z2.s, z4.s
/* #00001 */	ldr	z2, [x29, 14, mul vl]	//  (*)
/*    121 */	fmad	z27.s, p0/m, z5.s, z28.s
/* #00001 */	ldr	z4, [x29, 10, mul vl]	//  (*)
/*    120 */	fmad	z1.s, p0/m, z5.s, z0.s
/*     27 */	fmls	z2.s, p0/m, z16.s, z26.s
/* #00001 */	ldr	z26, [x29, 16, mul vl]	//  (*)
/*     26 */	fmul	z9.s, z8.s, z16.s
/*    119 */	fmad	z5.s, p0/m, z25.s, z6.s
/*     98 */	fmad	z11.s, p0/m, z11.s, z31.s
/*     84 */	fsub	z14.s, z29.s, z4.s
/* #00001 */	ldr	z4, [x29, 34, mul vl]	//  (*)
/*    107 */	ld1rw	{z25.s}, p0/z, [x3, 32]	//  (*)
/*     88 */	fsub	z8.s, z29.s, z4.s
/*    106 */	ldr	s29, [x3, 28]	//  (*)
/*    ??? */	ld1rw	{z6.s}, p0/z, [x3, 24]	//  (*)
/* #00001 */	ldr	z0, [x29, 15, mul vl]	//  (*)
/*     82 */	ld1rw	{z18.s}, p0/z, [x5]	//  (*)
/*     80 */	ldr	s28, [x5, 8]	//  (*)
/*     29 */	fmul	z4.s, z10.s, z15.s
/*     30 */	fmad	z0.s, p0/m, z10.s, z26.s
/* #00001 */	ldr	z26, [x29, 15, mul vl]	//  (*)
/*     94 */	fmla	z24.s, p0/m, z7.s, z7.s
/* #00001 */	ldr	z10, [x29, 16, mul vl]	//  (*)
/*     80 */	dup	z28.s, z28.s[0]
/*    106 */	dup	z7.s, z29.s[0]
/*    126 */	add	x3, x3, 32
/*     30 */	fmla	z10.s, p0/m, z2.s, z26.s
/* #00001 */	ldr	z26, [x29, 34, mul vl]	//  (*)
/*     29 */	fmul	z2.s, z2.s, z9.s
/*    115 */	fsub	z12.s, z6.s, z26.s
/* #00001 */	ldr	z26, [x29, 35, mul vl]	//  (*)
/*     20 */	frsqrte	z17.s, z11.s
/*     81 */	ldr	s29, [x5, 12]	//  (*)
/*     89 */	fsub	z31.s, z28.s, z26.s
/* #00001 */	ldr	z26, [x29, 13, mul vl]	//  (*)
/*     96 */	fmad	z8.s, p0/m, z8.s, z26.s
/*     92 */	fmad	z14.s, p0/m, z14.s, z26.s
/* #00001 */	ldr	z26, [x29, 36, mul vl]	//  (*)
/*     81 */	dup	z19.s, z29.s[0]
/*    117 */	fsub	z29.s, z25.s, z26.s
/*     32 */	fmad	z4.s, p0/m, z0.s, z15.s
/* #00001 */	ldr	z0, [x29, 11, mul vl]	//  (*)
/*     22 */	fmul	z16.s, z17.s, z17.s
/*     23 */	fmul	z17.s, z18.s, z17.s
/*     85 */	fsub	z26.s, z28.s, z0.s
/* #00001 */	ldr	z0, [x29, 12, mul vl]	//  (*)
/*     20 */	frsqrte	z20.s, z24.s
/*    113 */	fsub	z28.s, z25.s, z0.s
/* #00001 */	ldr	z0, [x29, 10, mul vl]	//  (*)
/*    111 */	fsub	z6.s, z6.s, z0.s
/* #00001 */	ldr	z0, [x29, 35, mul vl]	//  (*)
/*    116 */	fsub	z25.s, z7.s, z0.s
/* #00001 */	ldr	z0, [x29, 11, mul vl]	//  (*)
/*    112 */	fsub	z0.s, z7.s, z0.s
/* #00001 */	ldr	z7, [x29, 12, mul vl]	//  (*)
/*     86 */	fsub	z7.s, z19.s, z7.s
/*     32 */	fmad	z2.s, p0/m, z10.s, z9.s
/* #00001 */	ldr	z10, [x29, 36, mul vl]	//  (*)
/*     23 */	fmul	z15.s, z18.s, z20.s
/*     22 */	fmul	z9.s, z20.s, z20.s
/*     90 */	fsub	z10.s, z19.s, z10.s
/*     97 */	fmad	z31.s, p0/m, z31.s, z8.s
/* #00001 */	ldr	z8, [x29, 14, mul vl]	//  (*)
/*    125 */	fmad	z29.s, p0/m, z4.s, z30.s
/*    123 */	fmad	z12.s, p0/m, z4.s, z13.s
/*    126 */	add	x5, x5, 32
/*     79 */	ldr	s30, [x5, -12]	//  (*)
/*     27 */	fmsb	z11.s, p0/m, z16.s, z8.s
/*     26 */	fmul	z8.s, z17.s, z16.s
/*    101 */	bge	.L76 // end of swp_body
..LDL287:
/*    121 */	fmad	z28.s, p0/m, z2.s, z27.s
/*    119 */	fmad	z6.s, p0/m, z2.s, z5.s
/*    126 */	add	x16, x3, 16
/*    ??? */	ldr	z5, [x29, 16, mul vl]	//  (*)
/*     79 */	dup	z27.s, z30.s[0]
/*    ??? */	ldr	z30, [x29, 14, mul vl]	//  (*)
/*    120 */	fmla	z1.s, p0/m, z2.s, z0.s
/*    124 */	fmad	z25.s, p0/m, z4.s, z3.s
/*     26 */	fmul	z4.s, z15.s, z9.s
/*    ??? */	ldr	z2, [x29, 10, mul vl]	//  (*)
/*     93 */	fmad	z26.s, p0/m, z26.s, z14.s
/*     98 */	fmad	z10.s, p0/m, z10.s, z31.s
/*    105 */	ldr	s0, [x3, 8]	//  (*)
/*    107 */	ld1rw	{z31.s}, p0/z, [x16]	//  (*)
/*     29 */	fmul	z3.s, z11.s, z8.s
/*    105 */	dup	z17.s, z0.s[0]
/*    ??? */	ldr	z0, [x29, 16, mul vl]	//  (*)
/*     27 */	fmls	z30.s, p0/m, z9.s, z24.s
/*     94 */	fmla	z26.s, p0/m, z7.s, z7.s
/*    ??? */	ldr	z24, [x29, 34, mul vl]	//  (*)
/*     84 */	fsub	z9.s, z27.s, z2.s
/*    ??? */	ldr	z2, [x29, 15, mul vl]	//  (*)
/*     29 */	fmul	z7.s, z30.s, z4.s
/*    115 */	fsub	z13.s, z17.s, z24.s
/*    ??? */	ldr	z24, [x29, 36, mul vl]	//  (*)
/*     30 */	fmad	z2.s, p0/m, z11.s, z0.s
/*    ??? */	ldr	z0, [x29, 15, mul vl]	//  (*)
/*     82 */	ld1rw	{z11.s}, p0/z, [x4]	//  (*)
/*    126 */	add	x4, x5, 16
/*     32 */	fmad	z3.s, p0/m, z2.s, z8.s
/*    ??? */	ldr	z2, [x29, 11, mul vl]	//  (*)
/*    117 */	fsub	z16.s, z31.s, z24.s
/*     20 */	frsqrte	z24.s, z26.s
/*     30 */	fmad	z0.s, p0/m, z30.s, z5.s
/*    106 */	ldr	s5, [x3, 12]	//  (*)
/*    126 */	add	x3, x3, 32
/*    ??? */	ldr	z30, [x29, 12, mul vl]	//  (*)
/*     22 */	fmul	z8.s, z24.s, z24.s
/*    123 */	fmad	z13.s, p0/m, z3.s, z12.s
/*     23 */	fmul	z12.s, z11.s, z24.s
/*    ??? */	ldr	z24, [x29, 36, mul vl]	//  (*)
/*    106 */	dup	z18.s, z5.s[0]
/*    ??? */	ldr	z5, [x29, 13, mul vl]	//  (*)
/*    125 */	fmad	z16.s, p0/m, z3.s, z29.s
/*     32 */	fmad	z7.s, p0/m, z0.s, z4.s
/*     81 */	ldp	s4, s0, [x5, -8]	//  (*)
/*    ??? */	ldr	z29, [x29, 14, mul vl]	//  (*)
/*    113 */	fsub	z30.s, z31.s, z30.s
/*    ??? */	ldr	z31, [x29, 35, mul vl]	//  (*)
/*    112 */	fsub	z2.s, z18.s, z2.s
/*     81 */	dup	z0.s, z0.s[0]
/*     92 */	fmad	z9.s, p0/m, z9.s, z5.s
/*     20 */	frsqrte	z5.s, z10.s
/*     22 */	fmul	z15.s, z5.s, z5.s
/*     23 */	fmul	z14.s, z11.s, z5.s
/*    ??? */	ldr	z5, [x29, 10, mul vl]	//  (*)
/*    116 */	fsub	z31.s, z18.s, z31.s
/*     90 */	fsub	z11.s, z0.s, z24.s
/*    ??? */	ldr	z24, [x29, 12, mul vl]	//  (*)
/*    120 */	fmla	z1.s, p0/m, z7.s, z2.s
/*    ??? */	ldr	z2, [x29, 14, mul vl]	//  (*)
/*     27 */	fmls	z29.s, p0/m, z15.s, z10.s
/*     82 */	ld1rw	{z10.s}, p0/z, [x5]	//  (*)
/*    126 */	add	x5, x3, 16
/*    124 */	fmad	z3.s, p0/m, z31.s, z25.s
/*    ??? */	ldr	z31, [x29, 34, mul vl]	//  (*)
/*    111 */	fsub	z5.s, z17.s, z5.s
/*    ??? */	ldr	z25, [x29, 35, mul vl]	//  (*)
/*     86 */	fsub	z0.s, z0.s, z24.s
/*     80 */	dup	z17.s, z4.s[0]
/*    ??? */	ldr	z4, [x29, 11, mul vl]	//  (*)
/*     27 */	fmsb	z26.s, p0/m, z8.s, z2.s
/*     26 */	fmul	z2.s, z12.s, z8.s
/*    119 */	fmad	z5.s, p0/m, z7.s, z6.s
/*    105 */	ldr	s6, [x3, -8]	//  (*)
/*     88 */	fsub	z27.s, z27.s, z31.s
/*    ??? */	ldr	z31, [x29, 13, mul vl]	//  (*)
/*     89 */	fsub	z25.s, z17.s, z25.s
/*     85 */	fsub	z24.s, z17.s, z4.s
/*     26 */	fmul	z4.s, z14.s, z15.s
/*     93 */	fmad	z24.s, p0/m, z24.s, z9.s
/*     96 */	fmad	z27.s, p0/m, z27.s, z31.s
/*    105 */	dup	z9.s, z6.s[0]
/*    ??? */	ldr	z6, [x29, 16, mul vl]	//  (*)
/*     97 */	fmad	z25.s, p0/m, z25.s, z27.s
/*     94 */	fmla	z24.s, p0/m, z0.s, z0.s
/*    106 */	ldr	s0, [x3, -4]	//  (*)
/*    121 */	movprfx	z27.s, p0/z, z28.s
/*    121 */	fmla	z27.s, p0/m, z7.s, z30.s
/*    ??? */	ldr	z7, [x29, 15, mul vl]	//  (*)
/*    ??? */	ldr	z30, [x29, 15, mul vl]	//  (*)
/*    107 */	ld1rw	{z28.s}, p0/z, [x3]	//  (*)
/*    106 */	dup	z15.s, z0.s[0]
/*    ??? */	ldr	z0, [x29, 34, mul vl]	//  (*)
/*     98 */	fmad	z11.s, p0/m, z11.s, z25.s
/*     29 */	fmul	z25.s, z29.s, z4.s
/*     20 */	frsqrte	z14.s, z24.s
/*     30 */	fmad	z7.s, p0/m, z29.s, z6.s
/*     30 */	fmad	z30.s, p0/m, z26.s, z6.s
/*     29 */	fmul	z26.s, z26.s, z2.s
/*    ??? */	ldr	z6, [x29, 36, mul vl]	//  (*)
/*    115 */	fsub	z12.s, z9.s, z0.s
/*     20 */	frsqrte	z0.s, z11.s
/*     32 */	fmla	z4.s, p0/m, z7.s, z25.s
/*    ??? */	ldr	z7, [x29, 35, mul vl]	//  (*)
/*     32 */	fmla	z2.s, p0/m, z30.s, z26.s
/*    107 */	ld1rw	{z26.s}, p0/z, [x5]	//  (*)
/*     22 */	fmul	z31.s, z0.s, z0.s
/*     23 */	fmul	z8.s, z10.s, z0.s
/*    ??? */	ldr	z0, [x29, 10, mul vl]	//  (*)
/*    117 */	fsub	z29.s, z28.s, z6.s
/*    105 */	ldr	s30, [x3, 8]	//  (*)
/*    123 */	fmad	z12.s, p0/m, z4.s, z13.s
/*    116 */	fsub	z25.s, z15.s, z7.s
/*    ??? */	ldr	z7, [x29, 14, mul vl]	//  (*)
/*     26 */	fmul	z8.s, z8.s, z31.s
/*    125 */	fmad	z29.s, p0/m, z4.s, z16.s
/*    111 */	fsub	z6.s, z9.s, z0.s
/*    ??? */	ldr	z0, [x29, 12, mul vl]	//  (*)
/*     22 */	fmul	z9.s, z14.s, z14.s
/*    105 */	dup	z30.s, z30.s[0]
/*    124 */	fmla	z3.s, p0/m, z4.s, z25.s
/*    ??? */	ldr	z4, [x29, 16, mul vl]	//  (*)
/*     27 */	fmsb	z11.s, p0/m, z31.s, z7.s
/*    119 */	fmla	z5.s, p0/m, z2.s, z6.s
/*    ??? */	ldr	z6, [x29, 15, mul vl]	//  (*)
/*    113 */	fsub	z28.s, z28.s, z0.s
/*    ??? */	ldr	z0, [x29, 11, mul vl]	//  (*)
/*    121 */	fmla	z27.s, p0/m, z2.s, z28.s
/*    106 */	ldr	s28, [x3, 12]	//  (*)
/*    126 */	add	x3, x3, 32
/*    112 */	fsub	z0.s, z15.s, z0.s
/*     23 */	fmul	z15.s, z10.s, z14.s
/*    120 */	fmla	z1.s, p0/m, z2.s, z0.s
/*     27 */	movprfx	z0.s, p0/z, z24.s
/*     27 */	fmsb	z0.s, p0/m, z9.s, z7.s
/*    ??? */	ldr	z2, [x29, 16, mul vl]	//  (*)
/*     26 */	fmul	z25.s, z15.s, z9.s
/*     29 */	fmul	z7.s, z11.s, z8.s
/*     30 */	fmad	z6.s, p0/m, z11.s, z2.s
/*    ??? */	ldr	z2, [x29, 15, mul vl]	//  (*)
/*     32 */	fmad	z7.s, p0/m, z6.s, z8.s
/*    ??? */	ldr	z6, [x29, 10, mul vl]	//  (*)
/*     30 */	fmad	z2.s, p0/m, z0.s, z4.s
/*     29 */	fmul	z4.s, z0.s, z25.s
/*    ??? */	ldr	z0, [x29, 36, mul vl]	//  (*)
/*    111 */	fsub	z6.s, z30.s, z6.s
/*     32 */	fmad	z4.s, p0/m, z2.s, z25.s
/*    ??? */	ldr	z2, [x29, 35, mul vl]	//  (*)
/*    ??? */	ldr	z25, [x29, 11, mul vl]	//  (*)
/*    117 */	fsub	z16.s, z26.s, z0.s
/*    ??? */	ldr	z0, [x29, 12, mul vl]	//  (*)
/*    119 */	fmla	z5.s, p0/m, z4.s, z6.s
/*    125 */	fmad	z16.s, p0/m, z7.s, z29.s
/*    113 */	fsub	z24.s, z26.s, z0.s
/*    ??? */	ldr	z26, [x29, 34, mul vl]	//  (*)
/*    106 */	dup	z0.s, z28.s[0]
/*    116 */	fsub	z2.s, z0.s, z2.s
/*    112 */	fsub	z0.s, z0.s, z25.s
/*    121 */	fmla	z27.s, p0/m, z4.s, z24.s
/*    115 */	fsub	z13.s, z30.s, z26.s
/*    124 */	fmla	z3.s, p0/m, z7.s, z2.s
/*    120 */	fmla	z1.s, p0/m, z4.s, z0.s
/*    123 */	fmad	z13.s, p0/m, z7.s, z12.s
/*    126 */	cbz	w6, .L89
.L92:
	.p2align 5
.L95:					// :entr:term:mod:swpl
/*     79 */	ldr	s2, [x4, -12]	//  (*)
/*    126 */	subs	w6, w6, 1
/*    105 */	ldr	s4, [x3, -8]	//  (*)
/*     81 */	ldr	s7, [x4, -4]	//  (*)
/*    106 */	ldr	s6, [x3, -4]	//  (*)
/*     80 */	ldr	s0, [x4, -8]	//  (*)
/* #00002 */	ldr	z28, [x29, 36, mul vl]	//  (*)
/* #00002 */	ldr	z29, [x29, 11, mul vl]	//  (*)
/*     82 */	ld1rw	{z25.s}, p0/z, [x4]	//  (*)
/*    126 */	add	x4, x4, 16
/*    107 */	ld1rw	{z31.s}, p0/z, [x3]	//  (*)
/*    126 */	add	x3, x3, 16
/*    105 */	dup	z24.s, z4.s[0]
/* #00002 */	ldr	z4, [x29, 10, mul vl]	//  (*)
/*     79 */	dup	z2.s, z2.s[0]
/*     81 */	dup	z26.s, z7.s[0]
/*    106 */	dup	z30.s, z6.s[0]
/*     80 */	dup	z0.s, z0.s[0]
/*    112 */	fsub	z9.s, z30.s, z29.s
/* #00002 */	ldr	z29, [x29, 34, mul vl]	//  (*)
/*     84 */	fsub	z7.s, z2.s, z4.s
/* #00002 */	ldr	z4, [x29, 34, mul vl]	//  (*)
/*    115 */	fsub	z29.s, z24.s, z29.s
/*     88 */	fsub	z6.s, z2.s, z4.s
/* #00002 */	ldr	z2, [x29, 11, mul vl]	//  (*)
/* #00002 */	ldr	z4, [x29, 35, mul vl]	//  (*)
/*     85 */	fsub	z2.s, z0.s, z2.s
/*     89 */	fsub	z0.s, z0.s, z4.s
/* #00002 */	ldr	z4, [x29, 12, mul vl]	//  (*)
/*     86 */	fsub	z4.s, z26.s, z4.s
/*     90 */	fsub	z26.s, z26.s, z28.s
/* #00002 */	ldr	z28, [x29, 10, mul vl]	//  (*)
/*    111 */	fsub	z8.s, z24.s, z28.s
/* #00002 */	ldr	z24, [x29, 35, mul vl]	//  (*)
/* #00002 */	ldr	z28, [x29, 12, mul vl]	//  (*)
/*    116 */	fsub	z30.s, z30.s, z24.s
/* #00002 */	ldr	z24, [x29, 36, mul vl]	//  (*)
/*    113 */	fsub	z28.s, z31.s, z28.s
/*    117 */	fsub	z31.s, z31.s, z24.s
/* #00002 */	ldr	z24, [x29, 13, mul vl]	//  (*)
/*     92 */	fmad	z7.s, p0/m, z7.s, z24.s
/*     96 */	fmad	z6.s, p0/m, z6.s, z24.s
/*     93 */	fmad	z2.s, p0/m, z2.s, z7.s
/*     97 */	fmad	z0.s, p0/m, z0.s, z6.s
/*     94 */	fmad	z4.s, p0/m, z4.s, z2.s
/*     98 */	fmad	z26.s, p0/m, z26.s, z0.s
	.loc 1 20 0
..LDL487:
/*     20 */	frsqrte	z6.s, z4.s
/*     20 */	frsqrte	z0.s, z26.s
/*     22 */	fmul	z2.s, z6.s, z6.s
/*     23 */	fmul	z7.s, z25.s, z6.s
/*     22 */	fmul	z6.s, z0.s, z0.s
/*     23 */	fmul	z24.s, z25.s, z0.s
/* #00002 */	ldr	z0, [x29, 14, mul vl]	//  (*)
/* #00002 */	ldr	z25, [x29, 16, mul vl]	//  (*)
/*     26 */	fmul	z7.s, z7.s, z2.s
/*     27 */	fmsb	z4.s, p0/m, z2.s, z0.s
/*     26 */	fmul	z0.s, z24.s, z6.s
/* #00002 */	ldr	z2, [x29, 14, mul vl]	//  (*)
/* #00002 */	ldr	z24, [x29, 15, mul vl]	//  (*)
/*     27 */	fmsb	z26.s, p0/m, z6.s, z2.s
/*     29 */	fmul	z6.s, z4.s, z7.s
/* #00002 */	ldr	z2, [x29, 16, mul vl]	//  (*)
/*     30 */	fmad	z24.s, p0/m, z4.s, z2.s
/*     29 */	fmul	z2.s, z26.s, z0.s
/* #00002 */	ldr	z4, [x29, 15, mul vl]	//  (*)
/*     32 */	fmad	z6.s, p0/m, z24.s, z7.s
/*     30 */	fmad	z4.s, p0/m, z26.s, z25.s
	.loc 1 119 0
..LDL506:
/*    119 */	fmla	z5.s, p0/m, z6.s, z8.s
/*    120 */	fmla	z1.s, p0/m, z6.s, z9.s
/*    121 */	fmla	z27.s, p0/m, z6.s, z28.s
	.loc 1 32 0
..LDL509:
/*     32 */	fmad	z2.s, p0/m, z4.s, z0.s
	.loc 1 123 0
..LDL510:
/*    123 */	fmla	z13.s, p0/m, z2.s, z29.s
/*    124 */	fmla	z3.s, p0/m, z2.s, z30.s
/*    125 */	fmla	z16.s, p0/m, z2.s, z31.s
/*    126 */	bne	.L95
.L89:
.L78:					// :term
	.loc 1 131 0 is_stmt 1
..LDL514:
/*    131 */	st1w	{z5.s}, p0, [x11, 0, mul vl]	//  "acci_0"
	.loc 1 134 0
..LDL515:
/*    134 */	mov	x3, x14
	.loc 1 136 0
..LDL516:
/*    136 */	add	x9, x9, 512
	.loc 1 131 0
..LDL517:
/*    131 */	st1w	{z1.s}, p0, [x11, 1, mul vl]	//  "acci_0"
	.loc 1 134 0
..LDL518:
/*    134 */	add	x3, x3, x2
	.loc 1 136 0
..LDL519:
/*    136 */	add	x14, x14, 384
	.loc 1 131 0
..LDL520:
/*    131 */	st1w	{z27.s}, p0, [x11, 2, mul vl]	//  "acci_0"
	.loc 1 136 0
..LDL521:
/*    136 */	add	x8, x8, 512
/*    136 */	subs	w13, w13, 1
	.loc 1 132 0
..LDL522:
/*    132 */	st1w	{z13.s}, p0, [x7, 0, mul vl]	//  "acci_1"
/*    132 */	st1w	{z3.s}, p0, [x7, 1, mul vl]	//  "acci_1"
/*    132 */	st1w	{z16.s}, p0, [x7, 2, mul vl]	//  "acci_1"
	.loc 1 134 0
..LDL523:
/*    134 */	ld1w	{z0.s}, p0/z, [x11, 0, mul vl]	//  "acci_0"
/*    134 */	ld1w	{z1.s}, p0/z, [x11, 1, mul vl]	//  "acci_0"
/*    134 */	ld1w	{z2.s}, p0/z, [x11, 2, mul vl]	//  "acci_0"
/*    134 */	st3w	{z0.s, z1.s, z2.s}, p0, [x3, 0, mul vl]	//  (*)
	.loc 1 135 0
..LDL524:
/*    135 */	mov	x3, x12
	.loc 1 136 0
..LDL525:
/*    136 */	add	x12, x12, 384
	.loc 1 135 0
..LDL526:
/*    135 */	add	x3, x3, x2
/*    135 */	ld1w	{z0.s}, p0/z, [x7, 0, mul vl]	//  "acci_1"
/*    135 */	ld1w	{z1.s}, p0/z, [x7, 1, mul vl]	//  "acci_1"
/*    135 */	ld1w	{z2.s}, p0/z, [x7, 2, mul vl]	//  "acci_1"
/*    135 */	st3w	{z0.s, z1.s, z2.s}, p0, [x3, 0, mul vl]	//  (*)
	.loc 1 136 0
..LDL527:
/*    136 */	bne	.L73
.L81:					// :epi:term
	.loc 1 137 0
..LDL528:
/*    ??? */	ldr	z8, [x29, 54, mul vl]	//  (*)
/*    ??? */	ldr	z9, [x29, 55, mul vl]	//  (*)
/*    ??? */	ldr	z10, [x29, 56, mul vl]	//  (*)
/*    ??? */	ldr	z11, [x29, 57, mul vl]	//  (*)
/*    ??? */	ldr	z12, [x29, 58, mul vl]	//  (*)
/*    ??? */	ldr	z13, [x29, 59, mul vl]	//  (*)
/*    ??? */	ldr	z14, [x29, 60, mul vl]	//  (*)
/*    ??? */	ldr	z15, [x29, 61, mul vl]	//  (*)
/*    ??? */	ldr	z16, [x29, 62, mul vl]	//  (*)
/*    ??? */	ldr	z17, [x29, 63, mul vl]	//  (*)
/*    ??? */	ldr	z18, [x29, 64, mul vl]	//  (*)
/*    ??? */	ldr	z19, [x29, 65, mul vl]	//  (*)
/*    ??? */	ldr	z20, [x29, 66, mul vl]	//  (*)
/*    ??? */	ldr	x19, [x29, -8]	//  (*)
	.cfi_restore 19
/*    ??? */	add	sp, x29, 0
/*    ??? */	ldp	x29, x30, [sp]	//  (*)
	.cfi_restore 29
	.cfi_restore 30
/*    ??? */	addvl	sp, sp, 31
/*    ??? */	addvl	sp, sp, 31
/*    ??? */	addvl	sp, sp, 5
	.cfi_def_cfa_offset 0
/*    137 */	ret	
..D2.pchi:
	.cfi_endproc
.LFE0:
	.size	_Z15nbody_ext_inneriDvfPK4BodyP12AccelerationS_S_S_S2_, .-_Z15nbody_ext_inneriDvfPK4BodyP12AccelerationS_S_S_S2_
	.file 2 "/opt/FJSVxtclanga/tcsds-1.2.31/bin/../include/arm_sve.h"
	.file 3 "/opt/FJSVxos/devkit/aarch64/rfs/usr/include/bits/stdint-intn.h"
	.file 4 "/opt/FJSVxos/devkit/aarch64/rfs/usr/include/bits/types.h"
	.file 5 "/opt/FJSVxos/devkit/aarch64/rfs/usr/include/bits/stdint-uintn.h"
	.file 6 "/opt/FJSVxos/devkit/aarch64/rfs/usr/include/stdint.h"
	.pushsection	.text
..text.e:
	.popsection
	.section	.debug_info
	.4byte	.LSEdebug_info-.LSBdebug_info	// Length of .debug_info section
.LSBdebug_info:
	.2byte	0x4	// Version of DWARF information
	.4byte	.Ldebug_abbrev	// Offset into .debug_abbrev section
	.byte	0x8	// Address size
	.uleb128	0x1	// DW_TAG_compile_unit (0xb)
	.ascii	"kernel1.cpp\0"	// DW_AT_name
	.4byte	.Ldebug_line	// DW_AT_stmt_list
	.byte	0x4	// DW_AT_language
	.ascii	"/vol0004/rccs-atd/a01005/rsqrtCubed\0"	// DW_AT_comp_dir
	.ascii	"ccpcompx: Fujitsu C/C++ Compiler 4.5.0 (Mar  4 2021 13:29:06)\0"	// DW_AT_producer
	.4byte	.Ldebug_ranges2	// DW_AT_ranges
	.uleb128	0x2	// DW_TAG_subprogram (0x83)
	.4byte	0x11c	// DW_AT_sibling
	.ascii	"nbody_ext_inner\0"	// DW_AT_name
	.8byte	_Z15nbody_ext_inneriDvfPK4BodyP12AccelerationS_S_S_S2_	// DW_AT_low_pc
	.8byte	..D2.pchi-_Z15nbody_ext_inneriDvfPK4BodyP12AccelerationS_S_S_S2_	// DW_AT_high_pc
	.byte	0x1	// DW_AT_decl_file
	.byte	0x27	// DW_AT_decl_line
			// DW_AT_external
	.uleb128	0x1	// DW_AT_frame_base
	.byte	0x9c	// DW_OP_call_frame_cfa
	.ascii	"_Z15nbody_ext_inneriDvfPK4BodyP12AccelerationS_S_S_S2_\0"	// DW_AT_linkage_name
	.uleb128	0x3	// DW_TAG_FJ_loop (0xe3)
	.byte	0x1	// DW_AT_decl_file
	.byte	0x4e	// DW_AT_FJ_loop_start_line
	.byte	0x7e	// DW_AT_FJ_loop_end_line
	.byte	0x2	// DW_AT_FJ_loop_nest_level
	.byte	0x5	// DW_AT_FJ_loop_type
	.uleb128	0x3	// DW_TAG_FJ_loop (0xe9)
	.byte	0x1	// DW_AT_decl_file
	.byte	0x33	// DW_AT_FJ_loop_start_line
	.byte	0x88	// DW_AT_FJ_loop_end_line
	.byte	0x1	// DW_AT_FJ_loop_nest_level
	.byte	0x5	// DW_AT_FJ_loop_type
	.uleb128	0x4	// DW_TAG_inlined_subroutine (0xef)
	.4byte	0x11c	// DW_AT_abstract_origin
	.4byte	.Ldebug_ranges1	// DW_AT_ranges
	.byte	0x1	// DW_AT_call_file
	.byte	0x65	// DW_AT_call_line
	.uleb128	0x4	// DW_TAG_inlined_subroutine (0xfa)
	.4byte	0x11c	// DW_AT_abstract_origin
	.4byte	.Ldebug_ranges1	// DW_AT_ranges
	.byte	0x1	// DW_AT_call_file
	.byte	0x66	// DW_AT_call_line
	.uleb128	0x4	// DW_TAG_inlined_subroutine (0x105)
	.4byte	0x11c	// DW_AT_abstract_origin
	.4byte	.Ldebug_ranges1	// DW_AT_ranges
	.byte	0x1	// DW_AT_call_file
	.byte	0x65	// DW_AT_call_line
	.uleb128	0x4	// DW_TAG_inlined_subroutine (0x110)
	.4byte	0x11c	// DW_AT_abstract_origin
	.4byte	.Ldebug_ranges1	// DW_AT_ranges
	.byte	0x1	// DW_AT_call_file
	.byte	0x66	// DW_AT_call_line
	.byte	0x0	// End of children (0x83)
	.uleb128	0x5	// DW_TAG_subprogram (0x11c)
	.ascii	"_ZN33_INTERNAL_11_kernel1_cpp_1b15b59910rsqrtCubedEDvfS0_DvbS0_S0_S0_\0"	// DW_AT_name
	.byte	0x1	// DW_AT_inline
			// DW_AT_declaration
	.byte	0x0	// End of children (0xb)
.LSEdebug_info:
	.section	.debug_abbrev
.Ldebug_abbrev:
	.uleb128	0x1	// Abbreviation code
	.uleb128	0x11	// DW_TAG_compile_unit
	.byte	0x1	// DW_CHILDREN_yes
	.uleb128	0x3	// DW_AT_name
	.uleb128	0x8	// DW_FORM_string
	.uleb128	0x10	// DW_AT_stmt_list
	.uleb128	0x17	// DW_FORM_sec_offset
	.uleb128	0x13	// DW_AT_language
	.uleb128	0xb	// DW_FORM_data1
	.uleb128	0x1b	// DW_AT_comp_dir
	.uleb128	0x8	// DW_FORM_string
	.uleb128	0x25	// DW_AT_producer
	.uleb128	0x8	// DW_FORM_string
	.uleb128	0x55	// DW_AT_ranges
	.uleb128	0x17	// DW_FORM_sec_offset
	.byte	0x0
	.byte	0x0
	.uleb128	0x2	// Abbreviation code
	.uleb128	0x2e	// DW_TAG_subprogram
	.byte	0x1	// DW_CHILDREN_yes
	.uleb128	0x1	// DW_AT_sibling
	.uleb128	0x13	// DW_FORM_ref4
	.uleb128	0x3	// DW_AT_name
	.uleb128	0x8	// DW_FORM_string
	.uleb128	0x11	// DW_AT_low_pc
	.uleb128	0x1	// DW_FORM_addr
	.uleb128	0x12	// DW_AT_high_pc
	.uleb128	0x7	// DW_FORM_data8
	.uleb128	0x3a	// DW_AT_decl_file
	.uleb128	0xb	// DW_FORM_data1
	.uleb128	0x3b	// DW_AT_decl_line
	.uleb128	0xb	// DW_FORM_data1
	.uleb128	0x3f	// DW_AT_external
	.uleb128	0x19	// DW_FORM_flag_present
	.uleb128	0x40	// DW_AT_frame_base
	.uleb128	0x18	// DW_FORM_exprloc
	.uleb128	0x6e	// DW_AT_linkage_name
	.uleb128	0x8	// DW_FORM_string
	.byte	0x0
	.byte	0x0
	.uleb128	0x3	// Abbreviation code
	.uleb128	0xf000	// DW_TAG_FJ_loop
	.byte	0x0	// DW_CHILDREN_no
	.uleb128	0x3a	// DW_AT_decl_file
	.uleb128	0xb	// DW_FORM_data1
	.uleb128	0x3300	// DW_AT_FJ_loop_start_line
	.uleb128	0xb	// DW_FORM_data1
	.uleb128	0x3301	// DW_AT_FJ_loop_end_line
	.uleb128	0xb	// DW_FORM_data1
	.uleb128	0x3302	// DW_AT_FJ_loop_nest_level
	.uleb128	0xb	// DW_FORM_data1
	.uleb128	0x3303	// DW_AT_FJ_loop_type
	.uleb128	0xb	// DW_FORM_data1
	.byte	0x0
	.byte	0x0
	.uleb128	0x4	// Abbreviation code
	.uleb128	0x1d	// DW_TAG_inlined_subroutine
	.byte	0x0	// DW_CHILDREN_no
	.uleb128	0x31	// DW_AT_abstract_origin
	.uleb128	0x13	// DW_FORM_ref4
	.uleb128	0x55	// DW_AT_ranges
	.uleb128	0x17	// DW_FORM_sec_offset
	.uleb128	0x58	// DW_AT_call_file
	.uleb128	0xb	// DW_FORM_data1
	.uleb128	0x59	// DW_AT_call_line
	.uleb128	0xb	// DW_FORM_data1
	.byte	0x0
	.byte	0x0
	.uleb128	0x5	// Abbreviation code
	.uleb128	0x2e	// DW_TAG_subprogram
	.byte	0x0	// DW_CHILDREN_no
	.uleb128	0x3	// DW_AT_name
	.uleb128	0x8	// DW_FORM_string
	.uleb128	0x20	// DW_AT_inline
	.uleb128	0xb	// DW_FORM_data1
	.uleb128	0x3c	// DW_AT_declaration
	.uleb128	0x19	// DW_FORM_flag_present
	.byte	0x0
	.byte	0x0
	.byte	0x0
	.section	.debug_line
.Ldebug_line:
	.section	.debug_ranges
.Ldebug_ranges1:
	.8byte	0xffffffffffffffff	// Base addr selection entry ID
	.8byte	0x0
	.8byte	..LDL487
	.8byte	..LDL506
	.8byte	..LDL509
	.8byte	..LDL510
	.8byte	0x0
	.8byte	0x0
.Ldebug_ranges2:
	.8byte	0xffffffffffffffff	// Base addr selection entry ID
	.8byte	0x0
	.8byte	_Z15nbody_ext_inneriDvfPK4BodyP12AccelerationS_S_S_S2_
	.8byte	..D2.pchi
	.8byte	0x0
	.8byte	0x0
	.section	.note.GNU-stack,"",%progbits
	.section	.fj.compile_info, "e"
	.ascii	"C++::trad-libc++"
