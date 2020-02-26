/*********
 * MIT License
 * 
 * Copyright (c) 2018 NECSTLab, Politecnico di Milano
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *********/

/*
 *
 * utils.h
 *
 * authors:
 * Marco Rabozzi (marco.rabozzi@polimi.it)
 * Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)
 * Lorenzo Di Tucci (lorenzo.ditucci@polimi.it)
 *
 */

#ifndef UTILS_H
#define UTILS_H

#include "ap_int.h"
#include "hls_stream.h"

/*
 * This file contains a series of common parameters and functions that will
 * be used within the design
 */

#define PTS_SIZE 5
#define R_SIZE_0 PTS_SIZE
#define R_SIZE_1 9
#define E_SIZE_0 4
#define E_SIZE_1 R_SIZE_1
#define E_SIZE (E_SIZE_0 * E_SIZE_1)
#define A_SIZE 10
#define B_SIZE_0 3
#define B_SIZE_1 13
#define B_SIZE (B_SIZE_0 * B_SIZE_1)
#define C_SIZE 11
#define COEFFS_NUM C_SIZE
#define OUT_AP 512
#define OUT_LIMIT 287
#define OUT_SIZE Q_SIZE_1
#define ROOT_NUM 10
#define SOL_NUM 9
#define OUT_MODEL_SIZE (ROOT_NUM*10)
#define MODEL_SIZE 3

typedef ap_uint<64> point;
typedef ap_uint<288> rType;
typedef ap_uint<OUT_AP> outType;
typedef float my_type;

const int my_depth = 10000;

/*
 * stream duplication
 */
template<
	int SIZE,
	typename Type>
void dup_stream_inner(hls::stream<Type> &in, hls::stream<Type> &out_0, hls::stream<Type> &out_1){

	for(int i = 0; i < SIZE; i++){
#pragma HLS PIPELINE
		Type tmp = in.read();
		out_0.write(tmp);
		out_1.write(tmp);
	}
}

/*
 * wrapper for stream duplication
 */
template<
	int SIZE,
	typename Type>
void dup_stream(hls::stream<Type> &in, hls::stream<Type> &out_0, hls::stream<Type> &out_1, const int iter){

	for(int i = 0; i < iter; i++){
		dup_stream_inner<SIZE, Type>(in, out_0, out_1);
	}
}

/*
 * stream triplication
 */
template<
	int SIZE,
	typename Type>
void tri_stream_inner(hls::stream<Type> &in, hls::stream<Type> &out_0, hls::stream<Type> &out_1, hls::stream<Type> &out_2){

	for(int i = 0; i < SIZE; i++){
#pragma HLS PIPELINE
		Type tmp = in.read();
		out_0.write(tmp);
		out_1.write(tmp);
		out_2.write(tmp);
	}
}

/*
 * wrapper for stream triplication
 */
template<
	int SIZE,
	typename Type>
void tri_stream(hls::stream<Type> &in, hls::stream<Type> &out_0, hls::stream<Type> &out_1, hls::stream<Type> &out_2, const int iter){

	for(int i = 0; i < iter; i++){
		tri_stream_inner<SIZE, Type>(in, out_0, out_1, out_2);
	}
}

/*
 * stream to matrix
 */
template<
	int SIZE_0,
	int SIZE_1,
	typename Type>
void stream2mat(hls::stream<Type> &in, Type out[SIZE_0][SIZE_1]){
	for(int i = 0; i < SIZE_0; i++){
		for(int j = 0; j < SIZE_1; j++){
#pragma HLS PIPELINE
			Type val = in.read();
			out[i][j] = val;
		}
	}
}

/*
 * matrix to stream
 */
template<
	int SIZE_0,
	int SIZE_1,
	typename Type>
void mat2stream(Type in[SIZE_0][SIZE_1], hls::stream<Type> &out){
	for(int i = 0; i < SIZE_0; i++){
		for(int j = 0; j < SIZE_1; j++){
#pragma HLS PIPELINE
			Type val = in[i][j];
			out.write(val);
		}
	}
}

/*
 * stream to vector
 */
template<
	int SIZE_0,
	typename Type>
void stream2vector(hls::stream<Type> &in, Type out[SIZE_0]){
	for(int i = 0; i < SIZE_0; i++){
#pragma HLS PIPELINE
		Type val = in.read();
		out[i] = val;
	}
}

/*
 * vector to stream
 */
template<
	int SIZE_0,
	typename Type>
void vector2stream(Type in[SIZE_0], hls::stream<Type> &out){
	for(int i = 0; i < SIZE_0; i++){
#pragma HLS PIPELINE
		Type val = in[i];
		out.write(val);
	}
}

/*
 * optimized matrix multiplication between two 10x10 matrices
 */
template<
	typename Type>
void myMatMul10_10(hls::stream<Type> &in_streamA, hls::stream<Type> &in_streamB, hls::stream<Type> &out_streamC){

	Type A[A_SIZE][A_SIZE];
#pragma HLS ARRAY_PARTITION variable=A complete dim=2
	Type B[A_SIZE][A_SIZE];
#pragma HLS ARRAY_PARTITION variable=B complete dim=2
	Type C[A_SIZE][A_SIZE];

	for(int i = 0; i < A_SIZE; i++){
		for(int j = 0; j < A_SIZE; j++){
#pragma HLS PIPELINE
			Type tmpA = in_streamA.read();
			Type tmpB = in_streamB.read();
			A[i][j] = tmpA;
			B[i][j] = tmpB;
		}
	}

	for(int i = 0; i < A_SIZE; i++){
		for(int j = 0; j < A_SIZE; j++){
#pragma HLS PIPELINE
			Type tmp0 = A[i][0] * B[0][j];
			Type tmp1 = A[i][1] * B[1][j];
			Type tmp2 = A[i][2] * B[2][j];
			Type tmp3 = A[i][3] * B[3][j];
			Type tmp4 = A[i][4] * B[4][j];
			Type tmp5 = A[i][5] * B[5][j];
			Type tmp6 = A[i][6] * B[6][j];
			Type tmp7 = A[i][7] * B[7][j];
			Type tmp8 = A[i][8] * B[8][j];
			Type tmp9 = A[i][9] * B[9][j];

			Type tmp10 = tmp0 + tmp1;
			Type tmp11 = tmp2 + tmp3;
			Type tmp12 = tmp4 + tmp5;
			Type tmp13 = tmp6 + tmp7;
			Type tmp14 = tmp8 + tmp9;

			Type tmp15 = tmp10 + tmp11;
			Type tmp16 = tmp12 + tmp13;
			Type tmp17 = tmp14;

			Type tmp18 = tmp15 + tmp16;
			Type tmp19 = tmp17;

			C[i][j] = tmp18 + tmp19;

		}
	}

	for(int i = 0; i < A_SIZE; i++){
		for(int j = 0; j < A_SIZE; j++){
#pragma HLS PIPELINE
			Type tmpC = C[i][j];
			out_streamC.write(tmpC);
		}
	}

}

bool uchar_greater_3(unsigned char ch) {
	return ch & 0xFC;
}

bool is_float_zero(float in) {
#pragma HLS INLINE
	return ((*(unsigned long *) &in) & 0x7FFFFFFF) == 0x00000000;
}

float float_negate(float in) {
#pragma HLS INLINE
	union un_t {
		float f;
		unsigned long u;
	};

	union un_t *un = ((union un_t*) &in);
	un->u ^= 0x80000000;

	return un->f;
}

float float_div2(float in) {
#pragma HLS INLINE
	union un_t {
		float f;
		unsigned long u;
	};

	union un_t *un = ((union un_t*) &in);

	unsigned char exp = (un->u >> 23) & 0xFF;
	exp = exp == 0 ? 0 : exp - 1;
	un->u = ((unsigned long) exp) << 23 | (un->u & 0x807FFFFF);

	return un->f;
}

float float_mul2(float in) {
#pragma HLS INLINE

	union un_t {
		float f;
		unsigned long u;
	};

	union un_t *un = ((union un_t*) &in);

	unsigned char exp = (un->u >> 23) & 0xFF;
	exp = exp == 255 ? 255 : exp + 1;
	un->u = ((unsigned long) exp) << 23 | (un->u & 0x807FFFFF);

	return un->f;
}

float rsqrt(float x)
{
	#pragma HLS INLINE
	// int ihalf = *(int *)&x - 0x00800000; // Alternative to next line,
	// float xhalf = *(float *)&ihalf;      // for sufficiently large nos.
	float xhalf = 0.5f*x;
	int i = *(int *)&x;          // View x as an int.
	// i = 0x5f3759df - (i >> 1);   // Initial guess (traditional).
	i = 0x5f375a82 - (i >> 1);   // Initial guess (slightly better).
	x = *(float *)&i;            // View i as float.
	x = x*(1.5f - xhalf*x*x);    // Newton step.
	// x = x*(1.5008908 - xhalf*x*x);  // Newton step for a balanced error.
	return x;
}

/* This is rsqrt with an additional step of the Newton iteration, for
increased accuracy. The constant 0x5f37599e makes the relative error
range from 0 to -0.00000463.
   You can't balance the error by adjusting the constant. */
float rsqrt1(float x)
{
	#pragma HLS INLINE
	float xhalf = 0.5f*x;
	int i = *(int *)&x;          // View x as an int.
	i = 0x5f37599e - (i >> 1);   // Initial guess.
	x = *(float *)&i;            // View i as float.
	x = x*(1.5f - xhalf*x*x);    // Newton step.
	x = x*(1.5f - xhalf*x*x);    // Newton step again.
	return x;
}


bool is_float_negative(float in) {
#pragma HLS INLINE
	return ((*(unsigned long *) &in) & 0x80000000) != 0x00000000;
}



#endif
