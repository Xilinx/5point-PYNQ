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
 * step_3.h
 *
 * authors:
 * Marco Rabozzi (marco.rabozzi@polimi.it)
 * Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)
 * Lorenzo Di Tucci (lorenzo.ditucci@polimi.it)
 *
 * Perform step 3 of the five-point algorithm:
 * 1. generates matrices A1 and A2 starting from the null space,
 * 2. computes inv(A1) using the QR decomposition and back substitution
 * 3. computes A* = inv(A1)*A2
 *
 */

#ifndef STEP_3_H
#define STEP_3_H

#include "utils.h"
#include "core_coeffs_A1.h"
#include "core_coeffs_A2.h"
#include "qrUtils.h"

/*---------- part 1 ----------*/

/*
 * Generates matrix A1 starting from the null space
 */
void coeffs_A1(hls::stream<my_type> &e_stream, hls::stream<my_type> &A1_stream, const int iter){
	for(int i = 0; i < iter; i++){
		CFA1_multiCore(e_stream, A1_stream);
	}
}

/*
 * Generates matrix A2 starting from the null space
 */
void coeffs_A2(hls::stream<my_type> &e_stream, hls::stream<my_type> &A1_stream, const int iter){
	for(int i = 0; i < iter; i++){
		CFA2_multiCore(e_stream, A1_stream);
	}
}

/*---------- part 2 ----------*/

namespace hls{

/*
 * Generates initial Q and R matrices for the QR decomposition of A1
 */
template <
	int RowsA,
	int ColsA,
	typename Type>
 void genQrData_inner(hls::stream<Type> &in_stream, hls::stream<Type> &out_streamQ, hls::stream<Type> &out_streamR, const int iter) {
	// Copy input data to local R memory and initialize Q
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

/*
 * Wrapper for the generation of initial Q and R matrices
 */
template <
	int RowsA,
	int ColsA,
	typename Type>
 void genQrData(hls::stream<Type> &in_stream, hls::stream<Type> &out_streamQ, hls::stream<Type> &out_streamR, const int iter) {
	for(int i = 0; i < iter; i++){
		genQrData_inner<RowsA, ColsA, Type>(in_stream, out_streamQ, out_streamR, iter);
	}

}

/*
 * Custom QR factorization for 10x10 input matrix A1
 */
template <
    int RowsA,
    int ColsA,
    typename Type>
  void my_qrf_alt_top_10(hls::stream<Type> &in_stream,
          hls::stream<Type> &out_streamQ, hls::stream<Type> &out_streamR, const int iter) {

#pragma HLS INLINE
#pragma HLS DATAFLOW

    int sequence_0[5][3] = {
  		  {0, 1, 0},
			  {2, 3, 0},
			  {4, 5, 0},
			  {6, 7, 0},
			  {8, 9, 0}
    };

    int sequence_1[4][3] = {
        		  {1, 3, 1},
    			  {5, 7, 1},
    			  {0, 2, 0},
    			  {4, 6, 0}
          };

    int sequence_2[4][3] = {
        		  {3, 7, 2},
    			  {1, 9, 1},
    			  {2, 5, 1},
    			  {0, 8, 0}
          };

    int sequence_3[4][3] = {
        		  {3, 9, 2},
    			  {1, 6, 1},
    			  {2, 8, 1},
    			  {0, 4, 0}
          };

    int sequence_4[4][3] = {
        		  {7, 9, 3},
    			  {3, 5, 2},
    			  {6, 8, 2},
    			  {1, 2, 1}
          };

    int sequence_5[3][3] = {
        		  {5, 7, 3},
    			  {3, 6, 2},
    			  {1, 4, 1}
          };

    int sequence_6[3][3] = {
        		  {7, 9, 4},
    			  {5, 8, 3},
    			  {2, 3, 2}
          };

    int sequence_7[3][3] = {
        		  {7, 8, 4},
    			  {5, 6, 3},
    			  {2, 4, 2}
          };

    int sequence_8[3][3] = {
        		  {8, 9, 5},
    			  {6, 7, 4},
    			  {3, 5, 3}
          };

    int sequence_9[3][3] = {
        		  {7, 8, 5},
    			  {5, 6, 4},
    			  {3, 4, 3}
          };

    int sequence_10[3][3] = {
        		  {8, 9, 6},
    			  {6, 7, 5},
    			  {4, 5, 4}
          };

    int sequence_11[2][3] = {
        		  {7, 8, 6},
    			  {5, 6, 5}
          };

    int sequence_12[2][3] = {
        		  {8, 9, 7},
    			  {6, 7, 6}
          };

    int sequence_13[1][3] = {
        		  {7, 8, 7}
          };

    int sequence_14[1][3] = {
        		  {8, 9, 8}
    		};

    hls::stream<Type> q_stream[16];
#pragma HLS STREAM variable=q_stream depth=RowsA*ColsA*2 dim=1
    hls::stream<Type> r_stream[16];
#pragma HLS STREAM variable=r_stream depth=RowsA*ColsA*2 dim=1


      genQrData<RowsA, ColsA, Type>(in_stream, q_stream[0], r_stream[0], iter);
	  batch<5, RowsA, ColsA, false, Type>(q_stream[0], r_stream[0], q_stream[1], r_stream[1], sequence_0, iter);
	  batch<4, RowsA, ColsA, false, Type>(q_stream[1], r_stream[1], q_stream[2], r_stream[2], sequence_1, iter);
	  batch<4, RowsA, ColsA, false, Type>(q_stream[2], r_stream[2], q_stream[3], r_stream[3], sequence_2, iter);
	  batch<4, RowsA, ColsA, false, Type>(q_stream[3], r_stream[3], q_stream[4], r_stream[4], sequence_3, iter);
	  batch<4, RowsA, ColsA, false, Type>(q_stream[4], r_stream[4], q_stream[5], r_stream[5], sequence_4, iter);
	  batch<3, RowsA, ColsA, false, Type>(q_stream[5], r_stream[5], q_stream[6], r_stream[6], sequence_5, iter);
	  batch<3, RowsA, ColsA, false, Type>(q_stream[6], r_stream[6], q_stream[7], r_stream[7], sequence_6, iter);
	  batch<3, RowsA, ColsA, false, Type>(q_stream[7], r_stream[7], q_stream[8], r_stream[8], sequence_7, iter);
	  batch<3, RowsA, ColsA, false, Type>(q_stream[8], r_stream[8], q_stream[9], r_stream[9], sequence_8, iter);
	  batch<3, RowsA, ColsA, false, Type>(q_stream[9], r_stream[9], q_stream[10], r_stream[10], sequence_9, iter);
	  batch<3, RowsA, ColsA, false, Type>(q_stream[10], r_stream[10], q_stream[11], r_stream[11], sequence_10, iter);
	  batch<2, RowsA, ColsA, false, Type>(q_stream[11], r_stream[11], q_stream[12], r_stream[12], sequence_11, iter);
	  batch<2, RowsA, ColsA, false, Type>(q_stream[12], r_stream[12], q_stream[13], r_stream[13], sequence_12, iter);
	  batch<1, RowsA, ColsA, false, Type>(q_stream[13], r_stream[13], q_stream[14], r_stream[14], sequence_13, iter);
	  batch<1, RowsA, ColsA, true, Type>(q_stream[14], r_stream[14], out_streamQ, out_streamR, sequence_14, iter);

  }


}

/*---------- part 3 ----------*/

/*
 * Computation of reciprocal
 */
template <typename T> T my_back_substitute_recip(T x) {
    const T ONE = 1.0;
    return ONE/x;
  }

/*
 * First step of back substitute
 */
template<
	int RowsColsA,
	typename Type>
	void back_substitute_first_step_inner(hls::stream<Type> &A_in, hls::stream<Type> &A_out, hls::stream<Type> &B_out, hls::stream<Type> &rowSum_out, hls::stream<Type> &diag_in, hls::stream<Type> &diag_out, int i){
	Type diag[RowsColsA];
	Type diag_recip;
	Type diag_recip_low;
	Type subst_prod;
	Type subst_prod_sum;
	Type final_sum;
	Type subst_sum;
	Type subst_prod_m1;
	Type column_multiplier[RowsColsA];
	Type select_column_multiplier;
	Type neg_diag_prod;


	for(int k = 0; k < RowsColsA; k++){
#pragma HLS PIPELINE
		Type tmp = diag_in.read();
		diag[k] = tmp;
	}

	diag_recip = diag[i];
	diag_recip_low = diag_recip;

	Type A_local[RowsColsA];

	a_row_loop: for (int j=0;j<RowsColsA;j++) {
		 if (j>=i) {
			 b_col_if_loop: for (int k=0;k<RowsColsA;k++) {
#pragma HLS PIPELINE II = 1
				Type A_tmp = A_in.read();
				if(j == i){
					 A_local[k] = A_tmp;
				}
				if (k<=i) {
					if (i==j) {
						if (k==i) {
						  select_column_multiplier = diag_recip;
						}
						column_multiplier[k] = select_column_multiplier;
						A_out.write(A_tmp);
						B_out.write(select_column_multiplier);
						rowSum_out.write(0.0);
					} else {
						subst_prod_m1 = A_local[j]; // (A[j][i]) Working with a upper triangular matrix
						subst_prod = subst_prod_m1 * column_multiplier[k];
						subst_prod_sum = subst_prod; // Resize
						A_out.write(A_tmp);
						B_out.write(0.0);
						subst_sum = subst_prod_sum;
						rowSum_out.write(subst_sum);
					  }
				} else {
	              A_out.write(A_tmp);
	              B_out.write(0.0);
				  rowSum_out.write(0.0);
	            }
			 }
		 } else {
			 b_col_else_loop: for (int k=0;k<RowsColsA;k++) {
#pragma HLS PIPELINE II = 1
				Type A_tmp = A_in.read();
				A_out.write(A_tmp);
				B_out.write(0.0);
				rowSum_out.write(0.0);
			 }
		 }
	}

	for(int k = 0; k < RowsColsA; k++){
#pragma HLS PIPELINE
		Type tmp = diag[k];
		diag_out.write(tmp);
	}

}

/*
 * Wrapper for the first step of back substitute
 */
template<
	int RowsColsA,
	typename Type>
	void back_substitute_first_step(hls::stream<Type> &A_in, hls::stream<Type> &A_out, hls::stream<Type> &B_out, hls::stream<Type> &rowSum_out, hls::stream<Type> &diag_in, hls::stream<Type> &diag_out, int i, const int iter){
	for(int k = 0; k < iter; k++){
		back_substitute_first_step_inner<RowsColsA, Type>(A_in, A_out, B_out, rowSum_out, diag_in, diag_out, i);
	}
}

/*
 * Generic step of back substitute
 */
template<
	int RowsColsA,
	typename Type>
	void back_substitute_step_inner(hls::stream<Type> &A_in, hls::stream<Type> &A_out, hls::stream<Type> &B_in, hls::stream<Type> &B_out, hls::stream<Type> &rowSum_in, hls::stream<Type> &rowSum_out, hls::stream<Type> &diag_in, hls::stream<Type> &diag_out, int i){

	Type diag[RowsColsA];
	Type diag_recip;
	Type diag_recip_low;
	Type subst_prod;
	Type subst_prod_sum;
	Type final_sum;
	Type subst_sum;
	Type subst_prod_m1;
	Type column_multiplier[RowsColsA];
	Type select_column_multiplier;
	Type neg_diag_prod;


	for(int k = 0; k < RowsColsA; k++){
#pragma HLS PIPELINE
		Type tmp = diag_in.read();
		diag[k] = tmp;
	}

	diag_recip = diag[i];
	diag_recip_low = diag_recip;

	Type A_local[RowsColsA];

	a_row_loop: for (int j=0;j<RowsColsA;j++) {
		 if (j>=i) {
			 b_col_if_loop: for (int k=0;k<RowsColsA;k++) {
#pragma HLS PIPELINE II = 1
				 Type A_tmp = A_in.read();
				 Type B_tmp = B_in.read();
				 Type rowSum_tmp = rowSum_in.read();
				 if(j == i){
					 A_local[k] = A_tmp;
				 }
				 if (k<=i) {
					if (i==j) {
						if (k==i) {
						  select_column_multiplier = diag_recip;
						} else {
						  final_sum                = rowSum_tmp;
						  neg_diag_prod            = -final_sum * diag_recip_low;
						  select_column_multiplier = neg_diag_prod;
						}
						column_multiplier[k] = select_column_multiplier;
						A_out.write(A_tmp);
						B_out.write(select_column_multiplier);
						rowSum_out.write(rowSum_tmp);
					} else {
						subst_prod_m1 = A_local[j]; // (A[j][i]) Working with a upper triangular matrix
						subst_prod = subst_prod_m1 * column_multiplier[k];
						subst_prod_sum = subst_prod; // Resize
						A_out.write(A_tmp);
						B_out.write(B_tmp);
						if (k==i) {
						  // First accumulation in the row sum
						  subst_sum = subst_prod_sum;
						} else {
						  subst_sum = rowSum_tmp + subst_prod_sum;
						}
						rowSum_out.write(subst_sum);
					  }
				} else {
	              B_tmp = 0; // Zero lower triangle
	              A_out.write(A_tmp);
	              B_out.write(B_tmp);
	              rowSum_out.write(rowSum_tmp);
	            }
			 }
		 } else {
			 b_col_else_loop: for (int k=0;k<RowsColsA;k++) {
#pragma HLS PIPELINE II = 1
				Type A_tmp = A_in.read();
				Type B_tmp = B_in.read();
				Type rowSum_tmp = rowSum_in.read();
				A_out.write(A_tmp);
				B_out.write(B_tmp);
				rowSum_out.write(rowSum_tmp);
			 }
		 }
	}

	for(int k = 0; k < RowsColsA; k++){
#pragma HLS PIPELINE
		Type tmp = diag[k];
		diag_out.write(tmp);
	}
}

/*
 * Wrapper for a generic step of back substitute
 */
template<
	int RowsColsA,
	typename Type>
	void back_substitute_step(hls::stream<Type> &A_in, hls::stream<Type> &A_out, hls::stream<Type> &B_in, hls::stream<Type> &B_out, hls::stream<Type> &rowSum_in, hls::stream<Type> &rowSum_out, hls::stream<Type> &diag_in, hls::stream<Type> &diag_out, int i, const int iter){
	for(int k = 0; k < iter; k++){
		back_substitute_step_inner<RowsColsA, Type>(A_in, A_out, B_in, B_out, rowSum_in, rowSum_out, diag_in, diag_out, i);
	}
}

/*
 * Last step of back substitute
 */
template<
	int RowsColsA,
	typename Type>
	void back_substitute_last_step_inner(hls::stream<Type> &A_in, hls::stream<Type> &B_in, hls::stream<Type> &B_out, hls::stream<Type> &rowSum_in, hls::stream<Type> &diag_in, int i){

	Type diag[RowsColsA];
	Type diag_recip;
	Type diag_recip_low;
	Type subst_prod;
	Type subst_prod_sum;
	Type final_sum;
	Type subst_sum;
	Type subst_prod_m1;
	Type column_multiplier[RowsColsA];
	Type select_column_multiplier;
	Type neg_diag_prod;

	for(int k = 0; k < RowsColsA; k++){
#pragma HLS PIPELINE
		Type tmp = diag_in.read();
		diag[k] = tmp;
	}

	diag_recip = diag[i];
	diag_recip_low = diag_recip;

	Type A_local[RowsColsA];

	a_row_loop: for (int j=0;j<RowsColsA;j++) {
		 if (j>=i) {
			 b_col_if_loop: for (int k=0;k<RowsColsA;k++) {
#pragma HLS PIPELINE II = 1
				Type A_tmp = A_in.read();
				Type B_tmp = B_in.read();
				Type rowSum_tmp = rowSum_in.read();
				if(j == i){
					 A_local[k] = A_tmp;
				 }
				if (k<=i) {
					if (i==j) {
						if (k==i) {
						  select_column_multiplier = diag_recip;
						} else {
						  final_sum                = rowSum_tmp;
						  neg_diag_prod            = -final_sum * diag_recip_low;
						  select_column_multiplier = neg_diag_prod;
						}
						column_multiplier[k] = select_column_multiplier;
						B_out.write(select_column_multiplier);
					  } else {
						subst_prod_m1 = A_local[j]; // (A[j][i]) Working with a upper triangular matrix
						subst_prod = subst_prod_m1 * column_multiplier[k];
						subst_prod_sum = subst_prod; // Resize
						B_out.write(B_tmp);
					  }
				} else {
	              B_tmp = 0; // Zero lower triangle
	              B_out.write(B_tmp);
	            }
			 }
		 } else {
			 b_col_else_loop: for (int k=0;k<RowsColsA;k++) {
#pragma HLS PIPELINE II = 1
				Type A_tmp = A_in.read();
				Type B_tmp = B_in.read();
				Type rowSum_tmp = rowSum_in.read();
				B_out.write(B_tmp);
			 }
		 }
	}
}

/*
 * Wrapper for the last step of back substitute
 */
template<
	int RowsColsA,
	typename Type>
	void back_substitute_last_step(hls::stream<Type> &A_in, hls::stream<Type> &B_in, hls::stream<Type> &B_out, hls::stream<Type> &rowSum_in, hls::stream<Type> &diag_in, int i, const int iter){
	for(int k = 0; k < iter; k++){
		back_substitute_last_step_inner<RowsColsA, Type>(A_in, B_in, B_out, rowSum_in, diag_in, i);
	}
}

/*
 * Extraction of reciprocal diagonal matrix
 */
template<
    int RowsColsA,
    typename Type>
	void genDiag_inner(hls::stream<Type> &in_stream,
			hls::stream<Type> &out_stream, hls::stream<Type> &out_streamD){

	Type A[RowsColsA][RowsColsA];

		for(int i = 0; i < RowsColsA; i++){
			for(int j = 0; j < RowsColsA; j++){
#pragma HLS PIPELINE
				Type tmp = in_stream.read();
				A[i][j] = tmp;
			}
		}

		for(int i = 0; i < RowsColsA; i++){
#pragma HLS PIPELINE
			Type diag_recip_calc = my_back_substitute_recip(A[i][i]);
			out_streamD.write(diag_recip_calc);
		}

		for(int i = 0; i < RowsColsA; i++){
			for(int j = 0; j < RowsColsA; j++){
#pragma HLS PIPELINE
				Type tmp = A[i][j];
				out_stream.write(tmp);
			}
		}

}

/*
 * Wrapper for the extraction of reciprocal diagonal matrix
 */
template<
    int RowsColsA,
    typename Type>
	void genDiag(hls::stream<Type> &in_stream,
			hls::stream<Type> &out_stream, hls::stream<Type> &out_streamD, const int iter){
	for(int i = 0; i < iter; i++){
		genDiag_inner<RowsColsA, Type>(in_stream, out_stream, out_streamD);
	}
}

/*
 * Computation of matrix transpose
 */
template<
	int RowsColsA,
    typename Type>
void transpose_B_inner(hls::stream<Type> &in_stream, hls::stream<Type> &out_stream){

	Type B[RowsColsA][RowsColsA];
	for(int i = 0; i < RowsColsA; i++){
		for(int j = 0; j < RowsColsA; j++){
#pragma HLS PIPELINE
			Type tmp = in_stream.read();
			B[j][i] = tmp;
		}
	}

	for(int i = 0; i < RowsColsA; i++){
		for(int j = 0; j < RowsColsA; j++){
#pragma HLS PIPELINE
			Type tmp = B[i][j];
			out_stream.write(tmp);
		}
	}


}

/*
 * Wrapper for the computation of matrix transpose
 */
template<
	int RowsColsA,
    typename Type>
void transpose_B(hls::stream<Type> &in_stream, hls::stream<Type> &out_stream, const int iter){
	for(int i = 0; i < iter; i++){
		transpose_B_inner<RowsColsA, Type>(in_stream, out_stream);
	}
}

/*
 * Custom back substitute from the result of the QR factorization
 */
template<
    int RowsColsA,
    typename Type>
void my_back_substitute(hls::stream<Type> &in_stream,
		hls::stream<Type> &out_stream, const int iter) {
#pragma HLS INLINE

	hls::stream<Type> streamA[RowsColsA];
#pragma HLS STREAM variable=streamA depth=my_depth dim=1
	hls::stream<Type> streamB[RowsColsA];
#pragma HLS STREAM variable=streamB depth=my_depth dim=1
	hls::stream<Type> streamRowSum[RowsColsA];
#pragma HLS STREAM variable=streamRowSum depth=my_depth dim=1
	hls::stream<Type> d_stream[RowsColsA];
#pragma HLS STREAM variable=d_stream depth=my_depth dim=1

	genDiag<RowsColsA, Type>(in_stream, streamA[0], d_stream[0], iter);
	back_substitute_first_step<RowsColsA, Type>(streamA[0], streamA[1], streamB[0], streamRowSum[0], d_stream[0], d_stream[1], 0, iter);

	for(int i = 1; i < RowsColsA-1; i++){
#pragma HLS UNROLL
		back_substitute_step<RowsColsA, Type>(streamA[i], streamA[i+1], streamB[i-1], streamB[i], streamRowSum[i-1], streamRowSum[i], d_stream[i], d_stream[i+1], i, iter);
	}

   back_substitute_last_step<RowsColsA, Type>(streamA[RowsColsA-1], streamB[RowsColsA-2], streamB[RowsColsA-1], streamRowSum[RowsColsA-2], d_stream[RowsColsA-1], 9, iter);

   transpose_B<RowsColsA, Type>(streamB[RowsColsA-1], out_stream, iter);

}

/*---------- part 4 ----------*/

/*
 * Wrapper for the custom 10x10 matrix multiplication to obtain A* = inv(A1) * A2
 */
template<
	typename Type>
void matMulStep(hls::stream<Type> &in_streamA, hls::stream<Type> &in_streamB, hls::stream<Type> &out_streamC, const int iter){
	for(int i = 0; i < iter; i++){
		myMatMul10_10<Type>(in_streamA, in_streamB, out_streamC);
	}
}

#endif
