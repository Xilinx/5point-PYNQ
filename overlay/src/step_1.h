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
 * step_1.h
 *
 * authors:
 * Marco Rabozzi (marco.rabozzi@polimi.it)
 * Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)
 * Lorenzo Di Tucci (lorenzo.ditucci@polimi.it)
 *
 * This file contains the utilities to generate the 5x9 epipolar constraints matrix
 */

#ifndef STEP_1_H
#define STEP_1_H

#include "utils.h"

/*
 * Generation of the epipolar constraints matrix
 */
void genR_inner(hls::stream<point> &pts1, hls::stream<point> &pts2, hls::stream<rType> &out_stream){

	pts_loop:for(int i = 0; i < R_SIZE_0; i++){
	#pragma HLS PIPELINE

			point pts1_val = pts1.read();
			unsigned int tmp_pts1_x = pts1_val.range(31, 0);
			unsigned int tmp_pts1_y = pts1_val.range(63, 32);
			my_type pts1_x = *((my_type *)&tmp_pts1_x);
			my_type pts1_y = *((my_type *)&tmp_pts1_y);

			point pts2_val = pts2.read();
			unsigned int tmp_pts2_x = pts2_val.range(31, 0);
			unsigned int tmp_pts2_y = pts2_val.range(63, 32);
			my_type pts2_x = *((my_type *)&tmp_pts2_x);
			my_type pts2_y = *((my_type *)&tmp_pts2_y);

			my_type r0 = pts1_x * pts2_x;
			my_type r1 = pts1_y * pts2_x;
			my_type r2 = pts2_x;
			my_type r3 = pts1_x * pts2_y;
			my_type r4 = pts1_y * pts2_y;
			my_type r5 = pts2_y;
			my_type r6 = pts1_x;
			my_type r7 = pts1_y;
			my_type r8 = 1.0;

			rType r_val;
			r_val.range(31, 0) = *((unsigned int *)&r0);
			r_val.range(63, 32) = *((unsigned int *)&r1);
			r_val.range(95, 64) = *((unsigned int *)&r2);
			r_val.range(127, 96) = *((unsigned int *)&r3);
			r_val.range(159, 128) = *((unsigned int *)&r4);
			r_val.range(191, 160) = *((unsigned int *)&r5);
			r_val.range(223, 192) = *((unsigned int *)&r6);
			r_val.range(255, 224) = *((unsigned int *)&r7);
			r_val.range(287, 256) = *((unsigned int *)&r8);

			out_stream << r_val;

		}
}

/*
 * Wrapper of the generation of the epipolar constraints matrix
 */
void genR(hls::stream<point> &pts1, hls::stream<point> &pts2, hls::stream<rType> &out_stream, const int iter){

	for(int i = 0; i < iter; i++){
		genR_inner(pts1, pts2, out_stream);
	}
}

#endif
