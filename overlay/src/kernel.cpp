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
 * kernel.cpp
 *
 * Authors:
 * Marco Rabozzi (marco.rabozzi@polimi.it)
 * Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)
 * Lorenzo Di Tucci (lorenzo.ditucci@polimi.it)
 *
 *
 * This file contains the top function of our kernel
 *
 *********
 * Modifications:
 *
 * 01/31/2020 gn - changed kernel function name to "fivept_nister"
 *
 *********/

#include "step_1.h"
#include "step_2.h"
#include "step_3.h"
#include "step_4.h"
#include "step_5.h"
#include "step_6.h"
#include "step_7.h"
#include "memInterface.h"

/*
 * Streams depth
 * We set the depth of each stream according to the number of stages it crosses, plus a slack
 */
const int pts_depth=PTS_SIZE*4;
const int pts_1_depth=PTS_SIZE*105;
const int r_stream_depth=R_SIZE_0*2;
const int e_stream_depth=E_SIZE*2;
const int e_stream_dup_2_depth=E_SIZE*90;
const int A1_stream_depth=A_SIZE*A_SIZE*2;
const int A2_stream_depth=A_SIZE*A_SIZE*32;
const int Q_stream_depth=A_SIZE*A_SIZE*14;
const int B_stream_depth=B_SIZE_0*B_SIZE_1*2;
const int B_stream_dup_1_depth=B_SIZE_0*B_SIZE_1*38;
const int C_stream_depth=COEFFS_NUM*2;
const int t5_root_s_depth=ROOT_NUM*2;
const int sol_stream_depth=ROOT_NUM*SOL_NUM*2;
const int out_model_stream_depth=OUT_MODEL_SIZE*4;

extern "C" {

void fivept_nister(point* pts1_in, point* pts2_in, my_type* out, int iterations, my_type thresh){

#pragma HLS INTERFACE m_axi depth=10 port=pts1_in bundle=gmem0
#pragma HLS INTERFACE m_axi depth=10 port=pts2_in bundle=gmem1
#pragma HLS INTERFACE m_axi depth=200 port=out bundle=gmem2

#pragma HLS INTERFACE s_axilite register port=pts1_in bundle=control
#pragma HLS INTERFACE s_axilite register port=pts2_in bundle=control
#pragma HLS INTERFACE s_axilite register port=out bundle=control
#pragma HLS INTERFACE s_axilite register port=iterations bundle=control
#pragma HLS INTERFACE s_axilite register port=thresh bundle=control

#pragma HLS INTERFACE s_axilite register port=return bundle=control

#pragma HLS DATAFLOW

	// 4 times the proper depth to make sure that we hide the time for memory transfers
	hls::stream<point> pts1("pts1_stream"), pts2("pts2_stream");
#pragma HLS STREAM variable=pts1 depth=pts_depth dim=1
#pragma HLS STREAM variable=pts2 depth=pts_depth dim=1

	hls::stream<point> pts1_0("pts1_0_stream"), pts2_0("pts2_0_stream");
#pragma HLS STREAM variable=pts1_0 depth=pts_depth dim=1
#pragma HLS STREAM variable=pts2_0 depth=pts_depth dim=1

	hls::stream<point> pts1_1("pts1_1_stream"), pts2_1("pts2_1_stream");
#pragma HLS STREAM variable=pts1_1 depth=pts_1_depth dim=1
#pragma HLS STREAM variable=pts2_1 depth=pts_1_depth dim=1

	hls::stream<rType> r_stream("r_stream");
#pragma HLS STREAM variable=r_stream depth=r_stream_depth dim=1

	hls::stream<my_type> e_stream("e_stream");
#pragma HLS STREAM variable=e_stream depth=e_stream_depth dim=1

	hls::stream<my_type> e_stream_dup_0("e_stream_dup_0");
#pragma HLS STREAM variable=e_stream_dup_0 depth=e_stream_depth dim=1
	hls::stream<my_type> e_stream_dup_1("e_stream_dup_1");
#pragma HLS STREAM variable=e_stream_dup_1 depth=e_stream_depth dim=1

	// this stream crosses 84 stages
	hls::stream<my_type> e_stream_dup_2("e_stream_dup_2");
#pragma HLS STREAM variable=e_stream_dup_2 depth=e_stream_dup_2_depth dim=1

	hls::stream<my_type> A1_stream("A1_stream");
#pragma HLS STREAM variable=A1_stream depth=A1_stream_depth dim=1

	// this stream crosses 29 stages
	hls::stream<my_type> A2_stream("A2_stream");
#pragma HLS STREAM variable=A2_stream depth=A2_stream_depth dim=1

	// this stream crosses 12 stages
	hls::stream<my_type> Q_stream("Q_stream");
#pragma HLS STREAM variable=Q_stream depth=Q_stream_depth dim=1

	hls::stream<my_type> R_stream("R_stream");
#pragma HLS STREAM variable=R_stream depth=A1_stream_depth dim=1

	hls::stream<my_type> R_inv_stream("R_inv_stream");
#pragma HLS STREAM variable=R_inv_stream depth=A1_stream_depth dim=1

	hls::stream<my_type> A1_inv_stream("A1_inv_stream");
#pragma HLS STREAM variable=A1_inv_stream depth=A1_stream_depth dim=1

	hls::stream<my_type> A_stream("A_stream");
#pragma HLS STREAM variable=A_stream depth=A1_stream_depth dim=1

	hls::stream<my_type> B_stream("B_stream");
#pragma HLS STREAM variable=B_stream depth=B_stream_depth dim=1

	hls::stream<my_type> B_stream_dup_0("B_stream_dup_0");
#pragma HLS STREAM variable=B_stream_dup_0 depth=B_stream_depth dim=1

	// this stream crosses 34 stages
	hls::stream<my_type> B_stream_dup_1("B_stream_dup_1");
#pragma HLS STREAM variable=B_stream_dup_1 depth=B_stream_dup_1_depth dim=1

	hls::stream<my_type> C_stream("C_stream");
#pragma HLS STREAM variable=C_stream depth=C_stream_depth dim=1

	hls::stream<my_type> t5_root_s("t5_root_s");
#pragma HLS STREAM variable=t5_root_s depth=t5_root_s_depth dim=1

	// 4 times the proper depth to make sure that we hide the time for memory transfers
	hls::stream<my_type> sol_stream("sol_stream");
#pragma HLS STREAM variable=sol_stream depth=sol_stream_depth dim=1

	hls::stream<my_type> out_model_stream("out_model_stream");
#pragma HLS STREAM variable=out_model_stream depth=out_model_stream_depth dim=1

	int iter = iterations;
	my_type threshold = thresh;

	/*-----Read from memory-----*/

	axi2stream<point, PTS_SIZE>(pts1_in, pts1, iter);
	axi2stream<point, PTS_SIZE>(pts2_in, pts2, iter);

	dup_stream<PTS_SIZE, point>(pts1, pts1_0, pts1_1, iter);
	dup_stream<PTS_SIZE, point>(pts2, pts2_0, pts2_1, iter);

	/*-----Step 1-----*/

	genR(pts1_0, pts2_0, r_stream, iter);

	/*-----Step 2-----*/

	hls::my_qrf_alt_top_9<R_SIZE_1, R_SIZE_0, E_SIZE_0, rType, my_type>(r_stream, e_stream, iter);

	/*-----Step 3-----*/

	tri_stream<E_SIZE, my_type>(e_stream, e_stream_dup_0, e_stream_dup_1, e_stream_dup_2, iter);

	coeffs_A1(e_stream_dup_0, A1_stream, iter);
	coeffs_A2(e_stream_dup_1, A2_stream, iter);

	/*-----Step 4-----*/

	hls::my_qrf_alt_top_10<A_SIZE, A_SIZE, my_type>(A1_stream, Q_stream, R_stream, iter);

	my_back_substitute<A_SIZE, my_type>(R_stream, R_inv_stream, iter);

	matMulStep<my_type>(R_inv_stream, Q_stream, A1_inv_stream, iter);

	matMulStep<my_type>(A1_inv_stream, A2_stream, A_stream, iter);

	genB<A_SIZE, B_SIZE_0, B_SIZE_1, my_type>(A_stream, B_stream, iter);

	dup_stream<B_SIZE_0*B_SIZE_1, my_type>(B_stream, B_stream_dup_0, B_stream_dup_1, iter);

	coeffs_C(B_stream_dup_0, C_stream, iter);

	/*-----Step 5-----*/

	step_5_streaming(C_stream, t5_root_s, iter);

	/*-----Step 6-----*/

	step_6_streaming(t5_root_s, B_stream_dup_1, e_stream_dup_2, sol_stream, iter);

	/*-----Step 7-----*/

	findInliers(pts1_1, pts2_1, sol_stream, out_model_stream, threshold, iter);

	/*-----Write to memory-----*/

	stream2axi<my_type, OUT_MODEL_SIZE>(out_model_stream, out, iter);

}

}
