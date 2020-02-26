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
#include "qr_utils.h"

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
 * Generate of the initial Q and R matrices for a 3x3 QR decomposition.
 */
template <
	int RowsA,
	int ColsA,
	typename Type>
 void genQrData33_inner(hls::stream<Type> &in_stream, hls::stream<Type> &out_streamQ, hls::stream<Type> &out_streamR, const int iter) {
	// Copy input data to local R memory and initialize Q
	for(int it = 0; it < ROOT_NUM; it++) {
		#pragma HLS PIPELINE
		row_copy : for(int r=0; r<RowsA; r++){
		  // Merge loops to parallelize the A input read and the Q matrix prime.
		  #pragma HLS LOOP_MERGE force
		  col_copy_q_i : for(int c=0; c<RowsA; c++) {
			#pragma HLS PIPELINE
			Type q_val;
			if (r == c) {
			  q_val = 1.0;
			} else {
			  q_val = 0.0;
			}
			out_streamQ << q_val;
		  }
		  col_copy_r_i : for(int c=0; c<ColsA; c++) {
			#pragma HLS PIPELINE
			 Type tmp = in_stream.read();
			 out_streamR.write(tmp);
		  }
		}
	}
}

/*
 * Wrapper function for the generation of the Q and R matrices for the 3x3 QR decomposition.
 */
template <
	int RowsA,
	int ColsA,
	typename Type>
 void genQrData33(hls::stream<Type> &in_stream, hls::stream<Type> &out_streamQ, hls::stream<Type> &out_streamR, const int iter) {
	for(int i = 0; i < iter; i++){
		genQrData33_inner<RowsA, ColsA, Type>(in_stream, out_streamQ, out_streamR, iter);
	}

}

/*
 * Computes the batches for the 3x3 QR decomposition.
 * NOTE: this function returns only the last column of the Q matrix
 * since the others are not used for the recovery of the essential matrix
 * procedure. Indeed the 3x3 QR is used here to find the solution of
 * an underdetermined 3x3 system.
 */
template <
    int RowsA,
    int ColsA,
    typename Type>
  void my_qrf_alt_top_33(hls::stream<Type> &in_stream,
          hls::stream<Type> &out_streamQ, const int iter) {

#pragma HLS INLINE

    int sequence_0[1][3] = {
    	{0, 1, 0}
    };

    int sequence_1[1][3] = {
    	{0, 2, 0}
    };

    int sequence_2[1][3] = {
       {1, 2, 1}
    };

    hls::stream<Type> q_stream[3];
#pragma HLS STREAM variable=q_stream depth=10*RowsA*ColsA*2 dim=1
    hls::stream<Type> r_stream[3];
#pragma HLS STREAM variable=r_stream depth=10*RowsA*ColsA*2 dim=1

    genQrData33<RowsA, ColsA, Type>(in_stream, q_stream[0], r_stream[0], iter);
	hls::batch_trR<1, RowsA, ColsA, false, Type>(q_stream[0], r_stream[0], q_stream[1], r_stream[1], sequence_0, iter);
	hls::batch_loop<1, RowsA, ColsA, false, Type>(q_stream[1], r_stream[1], q_stream[2], r_stream[2], sequence_1, iter);
	hls::batch_last_Q_last_col<1, RowsA, ColsA, RowsA, true, Type>(q_stream[2], r_stream[2], out_streamQ, sequence_2, iter);

} // end qrf_alt

/*
 * Apply a 3x3 QR decomposition to solve a 3x3 underdetermined linear system.
 */
void solvez_streaming(hls::stream<my_type> &in, hls::stream<my_type> &out,
		const int iterations) {

#pragma HLS INLINE

	my_qrf_alt_top_33<3, 3, my_type>(in, out, iterations);
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
