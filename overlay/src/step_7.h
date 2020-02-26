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
 * step_7.h
 *
 * authors:
 * Marco Rabozzi (marco.rabozzi@polimi.it)
 * Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)
 * Lorenzo Di Tucci (lorenzo.ditucci@polimi.it)
 *
 * This file contains the computation of the forward error
 *
 */

#ifndef STEP_7_H
#define STEP_7_H

#include "utils.h"


typedef float DATA_TYPE;
typedef ap_uint<64> point;

/*
 * Computation of the forward error
 */
void computeError(DATA_TYPE m1[PTS_SIZE][2], DATA_TYPE m2[PTS_SIZE][2], DATA_TYPE model[ROOT_NUM][MODEL_SIZE*MODEL_SIZE+1], DATA_TYPE err[ROOT_NUM][PTS_SIZE]){

	computeError_loop_0:for(int k = 0; k < ROOT_NUM; k++){
		computeError_loop_1:for(int i = 0; i < PTS_SIZE; i++){
	#pragma HLS PIPELINE
			DATA_TYPE Ex1[3];
			DATA_TYPE Etx2[3];
			DATA_TYPE x2tEx1;
			DATA_TYPE x1[3] = {m1[i][0], m1[i][1], 1.0f};
			DATA_TYPE x2[3] = {m2[i][0], m2[i][1], 1.0f};

			build_vect:for(int j = 0; j < 3; j++){
				Ex1[j] = 0;
				Etx2[j] = 0;
			}
			x2tEx1 = 0;

			DATA_TYPE red_0 = model[k][1] * x1[0];
			DATA_TYPE red_1 = model[k][2] * x1[1];
			DATA_TYPE red_2 = model[k][3] * x1[2];

			Ex1[0] = red_0 + red_1 + red_2;

			DATA_TYPE red_3 = model[k][4] * x1[0];
			DATA_TYPE red_4 = model[k][5] * x1[1];
			DATA_TYPE red_5 = model[k][6] * x1[2];

			Ex1[1] = red_3 + red_4 + red_5;

			DATA_TYPE red_6 = model[k][7] * x1[0];
			DATA_TYPE red_7 = model[k][8] * x1[1];
			DATA_TYPE red_8 = model[k][9] * x1[2];

			Ex1[2] = red_6 + red_7 + red_8;

			DATA_TYPE red_9 = model[k][1] * x2[0];
			DATA_TYPE red_10 = model[k][2] * x2[1];
			DATA_TYPE red_11 = model[k][3] * x2[2];

			Etx2[0] =red_9 + red_10 + red_11;

			DATA_TYPE red_12 = model[k][4] * x2[0];
			DATA_TYPE red_13 = model[k][5] * x2[1];
			DATA_TYPE red_14 = model[k][6] * x2[2];

			Etx2[1] =red_12 + red_13 + red_14;

			DATA_TYPE red_15 = model[k][7] * x2[0];
			DATA_TYPE red_16 = model[k][8] * x2[1];
			DATA_TYPE red_17 = model[k][9] * x2[2];

			Etx2[1] =red_15 + red_16 + red_17;

	// Maronn how much is bell rabotz

			DATA_TYPE red_18 = x2[0] * Ex1[0];
			DATA_TYPE red_19 = x2[1] * Ex1[1];
			DATA_TYPE red_20 = x2[2] * Ex1[2];

			x2tEx1 = red_18 + red_19 + red_20;


			DATA_TYPE a = Ex1[0] * Ex1[0];
			DATA_TYPE b = Ex1[1] * Ex1[1];
			DATA_TYPE c = Etx2[0] * Etx2[0];
			DATA_TYPE d = Etx2[1] * Etx2[1];

			err[k][i] =(float)(x2tEx1 * x2tEx1 / (a + b + c + d));

		}
	}

}

/*
 * Calculation of the number of errors below a provided threshold
 */
void findInliersInner(hls::stream<point> &m1, hls::stream<point> &m2, hls::stream<DATA_TYPE> &model, hls::stream<DATA_TYPE> &model_out, DATA_TYPE thresh/*model out ha prima nz_out e poi il modello*/){
	DATA_TYPE err[ROOT_NUM][PTS_SIZE];
	DATA_TYPE m1_l[PTS_SIZE][2], m2_l[PTS_SIZE][2], model_l[ROOT_NUM][MODEL_SIZE*MODEL_SIZE+1];
#pragma HLS ARRAY_PARTITION variable=model_l complete dim=2

	read_points:for(int i = 0; i < PTS_SIZE; i++){
#pragma HLS PIPELINE
		point p1, p2;
		p1 = m1.read();
		p2 = m2.read();

		unsigned int tmp_pts1_x = p1.range(31, 0);
		unsigned int tmp_pts1_y = p1.range(63, 32);
		m1_l[i][0] = *((DATA_TYPE *)&tmp_pts1_x);
		m1_l[i][1] = *((DATA_TYPE *)&tmp_pts1_y);

		unsigned int tmp_pts2_x = p2.range(31, 0);
		unsigned int tmp_pts2_y = p2.range(63, 32);
		m2_l[i][0] = *((DATA_TYPE *)&tmp_pts2_x);
		m2_l[i][1] = *((DATA_TYPE *)&tmp_pts2_y);
	}

	read_model_0:for(int k = 0; k < ROOT_NUM; k++){
		read_model_1:for(int i = 1; i < MODEL_SIZE*MODEL_SIZE + 1; i++){
			//read_model_2:for(int j = 0; j < MODEL_SIZE; j++){
	#pragma HLS PIPELINE
			model_l[k][i] = model.read();
		}
	}

	computeError(m1_l, m2_l, model_l, err);

	DATA_TYPE t = (float)(thresh * thresh);

	int nz[ROOT_NUM] = {0};
#pragma HLS ARRAY_PARTITION variable=nz complete dim=1
	compute_nz_0:for(int k = 0; k < ROOT_NUM; k++){
		compute_nz_1:for(int i = 0; i < PTS_SIZE; i++){
	#pragma HLS PIPELINE
			int f = err[k][i] <= t;
			nz[k] += f;
		}
	}

	write_nz:for(int k = 0; k < ROOT_NUM; k++){
#pragma HLS PIPELINE
		model_l[k][0] = nz[k];
	}

	write_out_0:for(int k = 0; k < ROOT_NUM; k++){
		write_out_1:for(int i = 0; i < MODEL_SIZE*MODEL_SIZE+1; i++){
		#pragma HLS PIPELINE
				model_out.write((DATA_TYPE)model_l[k][i] );
		}
	}

}

/*
 * Wrapper for the calculation of the number of errors below a provided threshold
 */
void findInliers(hls::stream<point> &pts1, hls::stream<point> &pts2, hls::stream<DATA_TYPE> &model, hls::stream<DATA_TYPE> &out_model, DATA_TYPE threshold, const int iter){

	for(int i = 0; i < iter; i++){
		findInliersInner(pts1, pts2, model, out_model, threshold);
	}
}

#endif



