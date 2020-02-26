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
 * step_4.h
 *
 * authors:
 * Marco Rabozzi (marco.rabozzi@polimi.it)
 * Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)
 * Lorenzo Di Tucci (lorenzo.ditucci@polimi.it)
 *
 * Perform step 4 of the five-point algorithm:
 * computes the polynomial system B in unknown z and derives
 * the coefficient of a corresponding 10-th degree polynomial
 *
 */

#ifndef STEP_4_H
#define STEP_4_H

#include "utils.h"
#include "core_coeffs_C.h"

/*
 * Matrix B generation
 */
template<
	int IN_SIZE,
	int OUT_SIZE_0,
	int OUT_SIZE_1,
	typename Type>
void genB_inner(hls::stream<Type> &in_stream, hls::stream<Type> &out_stream){

	Type A[IN_SIZE][IN_SIZE];
#pragma HLS ARRAY_PARTITION variable=A complete dim=2
	Type B[OUT_SIZE_0*OUT_SIZE_1];
#pragma HLS ARRAY_PARTITION variable=B cyclic factor=OUT_SIZE_1 dim=1


	for(int i = 0; i < IN_SIZE; i++){
		for(int j = 0; j < IN_SIZE; j++){
#pragma HLS PIPELINE
			Type out = in_stream.read();
			A[i][j] = out;
		}
	}


	b_loop:for(int i = 0, j = 0; i < OUT_SIZE_0; i++, j+=OUT_SIZE_1){
	#pragma HLS PIPELINE
		B[j] = -A[i*2 + 5][0];
		B[j + 1] = A[i*2 + 4][0] - A[i*2 + 5][1];
		B[j + 2] = A[i*2 + 4][1] - A[i*2 + 5][2];
		B[j + 3] = A[i*2 + 4][2];
		B[j + 4] = -A[i*2 + 5][3];
		B[j + 5] = A[i*2 + 4][3] - A[i*2 + 5][4];
		B[j + 6] = A[i*2 + 4][4] - A[i*2 + 5][5];
		B[j + 7] = A[i*2 + 4][5];
		B[j + 8] = -A[i*2 + 5][6];
		B[j + 9] = A[i*2 + 4][6] - A[i*2 + 5][7];
		B[j + 10] = A[i*2 + 4][7] - A[i*2 + 5][8];
		B[j + 11] = A[i*2 + 4][8] - A[i*2 + 5][9];
		B[j + 12] = A[i*2 + 4][9];

	}

	for(int i = 0; i < OUT_SIZE_0*OUT_SIZE_1; i++){
#pragma HLS PIPELINE
		Type out = B[i];
		out_stream.write(out);
	}
}

/*
 * Wrapper for the generation of B matrix
 */
template<
	int IN_SIZE,
	int OUT_SIZE_0,
	int OUT_SIZE_1,
	typename Type>
void genB(hls::stream<Type> &in_stream, hls::stream<Type> &out_stream, const int iter){
	for(int i = 0; i < iter; i++){
		genB_inner<IN_SIZE, OUT_SIZE_0, OUT_SIZE_1, Type>(in_stream, out_stream);
	}
}

/*
 * Wrapper for the computation of the 10-th degree polynomial coefficients
 */
void coeffs_C(hls::stream<my_type> &b_stream, hls::stream<my_type> &c_stream, const int iter){
	for(int i = 0; i < iter; i++){
		CFC_multiCore(b_stream, c_stream);
	}
}

#endif
