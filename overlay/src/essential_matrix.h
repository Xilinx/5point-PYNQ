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
 * essential_matrix.h
 *
 * authors:
 * Marco Rabozzi (marco.rabozzi@polimi.it)
 * Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)
 * Lorenzo Di Tucci (lorenzo.ditucci@polimi.it)
 *
 * utility functions used by step_56.h to compute the essential matrix.
 */

#ifndef ESSENTIAL_MATRIX_H
#define ESSENTIAL_MATRIX_H

#include "hls_stream.h"
#include "ap_int.h"
#include "polysolve.h"
#include <math.h>
#include "qrUtils.h"

/*****************************************************************************
 * Functions
 ******************************************************************************/

const int t6_input_A2_s_0_depth=9*ROOT_NUM*2;
const int t6_input_A2_s_1_depth=9*ROOT_NUM*18;
const int t6_eigen_s_depth=6*ROOT_NUM*2;
const int t6_eigen_p_depth=2*ROOT_NUM*2;
const int t6_eigen_qV_depth=4*ROOT_NUM*2;
const int t6_V_s_depth=9*ROOT_NUM*2;

/*
 * Apply a 3x3 QR decomposition to solve a 3x3 underdetermined linear system.
 */
void solvez_streaming(hls::stream<my_type> &in, hls::stream<my_type> &out,
		const int iterations) {

#pragma HLS INLINE

	hls::my_qrf_alt_top_33<3, 3, my_type>(in, out, iterations);

}

/*
 * Utility function used to generate the underdetermined 3x3 linear system.
 */
void getSystemsMatrix_inner_read(hls::stream<float> &in_b,
		hls::stream<my_type> &in_r, my_type b[39], my_type z1[10],
		my_type z2[10], my_type z3[10], my_type z4[10])
{
	#pragma HLS INLINE off

	for (int i = 0; i < 39; i++) {
		#pragma HLS PIPELINE
		b[i] = in_b.read();
	}

	for (int i = 0; i < 10; i++) {
		#pragma HLS PIPELINE
		my_type tmp = in_r.read();
		my_type tmp2 = tmp*tmp;
		my_type tmp3 = tmp2*tmp;
		my_type tmp4 = tmp2*tmp2;

		z1[i] = tmp;
		z2[i] = tmp2;
		z3[i] = tmp3;
		z4[i] = tmp4;
	}
}

/*
 * Computes the matrix associated to a 3x3 linear system.
 */
void getSystemsMatrix_inner(hls::stream<float> &in_b,
		hls::stream<my_type> &in_r, hls::stream<my_type> &out)
{
	my_type b[39];
	my_type z1[10];
	my_type z2[10];
	my_type z3[10];
	my_type z4[10];
	#pragma HLS ARRAY_PARTITION variable=b dim=1

	getSystemsMatrix_inner_read(in_b, in_r, b, z1, z2, z3, z4);

	for (int i = 0; i < MAX_ORDER; i++) {
		for (int j = 0; j < 3; j++) {
			#pragma HLS PIPELINE II=3
			my_type tmp, tmp2;
			tmp = b[j*13 + 0] * z3[i] + b[j*13 + 1] * z2[i] + b[j*13 + 2] * z1[i] + b[j*13 + 3];
			tmp2 = z1[i] != -INFINITY ? tmp : -INFINITY;
			out.write(tmp2);
			tmp = b[j*13 + 4] * z3[i] + b[j*13 + 5] * z2[i] + b[j*13 + 6] * z1[i] + b[j*13 + 7];
			tmp2 = z1[i] != -INFINITY ? tmp : -INFINITY;
			out.write(tmp2);
			tmp = b[j*13 + 8] * z4[i] + b[j*13 + 9] * z3[i] + b[j*13 + 10] * z2[i] + b[j*13 + 11] * z1[i] + b[j*13 + 12];
			tmp2 = z1[i] != -INFINITY ? tmp : -INFINITY;
			out.write(tmp2);
		}
	}
}

/*
 * Wrapper function for computing the matrix associated to a 3x3 linear system.
 */
void getSystemsMatrix(hls::stream<my_type> &in_b,
		hls::stream<my_type> &in_r, hls::stream<my_type> &out,
		const int iterations) {

	for(int i = 0; i < iterations; i++) {
		getSystemsMatrix_inner(in_b, in_r, out);
	}
}

/*
 * Recover the unknown X and Y from the 3x3 system solution.
 * Such unknown are used together with the root of the polynomial
 * to recover the essential matrix.
 */
void recoverMatrixFromXY_inner(hls::stream<my_type> &in_ee_s,
		hls::stream<my_type> &in_sol, hls::stream<my_type> &in_roots,
		hls::stream<my_type> &out) {

	my_type sol[3];
	my_type in_ee[4][9];

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 9; j++) {
			#pragma HLS PIPELINE
			in_ee[i][j] = in_ee_s.read();
		}
	}

	for (int s = 0; s < MAX_ORDER; s++) {
		#pragma HLS PIPELINE

		for(int i = 0; i < 3; i++) {
			sol[i] = in_sol.read();
		}

		float scale = sol[2];
		float ee_out[9];
		float w[4];

		w[0] = sol[0] / scale;
		w[1] = sol[1] / scale;
		w[2] = in_roots.read();
		w[3] = 1.0f;

		// compute EE
		for (int i = 0; i < 9; i++) {
			#pragma HLS UNROLL
			ee_out[i] = 0;
		}

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 9; j++) {
				#pragma HLS UNROLL
				ee_out[j] += in_ee[i][j] * w[i];
			}
		}

		// normalize EE
		float normScale = 0.0f;

		for (int i = 0; i < 9; i++) {
			#pragma HLS UNROLL
			normScale += ee_out[i] * ee_out[i];
		}

		normScale = 1.0f / sqrt(normScale);

		for (int i = 0; i < 9; i++) {
			#pragma HLS UNROLL
			my_type tmp = sol[0] != -INFINITY ? ee_out[i]*normScale : -INFINITY;
			out.write(tmp);
		}
	}
}

/*
 * Wrapper function for recovering the X and Y values from the 3x3 system.
 */
void recoverMatrixFromXY(hls::stream<my_type> &in_ee,
		hls::stream<my_type> &in_sol, hls::stream<my_type> &in_roots,
		hls::stream<my_type> &out,
		const int iterations) {
	for(int i = 0; i < iterations; i++){
		recoverMatrixFromXY_inner(in_ee, in_sol, in_roots, out);
	}
}


#endif
