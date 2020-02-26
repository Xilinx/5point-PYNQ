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
 * step_6.h
 *
 * authors:
 * Marco Rabozzi (marco.rabozzi@polimi.it)
 * Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)
 * Lorenzo Di Tucci (lorenzo.ditucci@polimi.it)
 *
 * Step 6 of the five-point algorithm:
 * recover the essential matrices, one for each root of the 10-th degree polynomial
 */

#ifndef STEP_6_H
#define STEP_6_H

#include "utils.h"
#include "essential_matrix.h"

const int t6_ls_v_s_depth=9*ROOT_NUM*2;
const int t6_sol_ls_s_depth=3*ROOT_NUM*2;
const int t6_root_s_depth_dup1=ROOT_NUM*2;
const int t6_root_s_depth_dup2=ROOT_NUM*22;

/*
 * Recover the essential matrix starting form the roots of the 10-th degree polynomial
 * the system B and the null space stored in EE_stream.
 */
void step_6_streaming(
		hls::stream<my_type> &root_s,
		hls::stream<my_type> &B_in,
		hls::stream<my_type> &EE_stream,
		hls::stream<my_type> &out,
		const int iterations)
{
	#pragma HLS INLINE
	#pragma HLS DATAFLOW

	hls::stream<my_type> t6_ls_v_s("t6_ls_v_s");
#pragma HLS STREAM variable=t6_ls_v_s depth=t6_ls_v_s_depth dim=1

	hls::stream<my_type> t6_sol_ls_s("t6_sol_ls_s");
#pragma HLS STREAM variable=t6_sol_ls_s depth=t6_sol_ls_s_depth dim=1

	hls::stream<my_type> t6_root_s_dup1("t6_root_s_dup1");
#pragma HLS STREAM variable=t6_root_s_dup1 depth=t6_root_s_depth_dup1 dim=1

	hls::stream<my_type> t6_root_s_dup2("t6_root_s_dup2");
#pragma HLS STREAM variable=t6_root_s_dup2 depth=t6_root_s_depth_dup2 dim=1

	dup_stream<ROOT_NUM, my_type>(root_s, t6_root_s_dup1, t6_root_s_dup2, iterations);

	getSystemsMatrix(B_in, t6_root_s_dup1, t6_ls_v_s, iterations);

	solvez_streaming(t6_ls_v_s, t6_sol_ls_s, iterations);

	recoverMatrixFromXY(EE_stream, t6_sol_ls_s, t6_root_s_dup2, out, iterations);
}

#endif
