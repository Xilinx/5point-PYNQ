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
 * step_2.h
 *
 * authors:
 * Marco Rabozzi (marco.rabozzi@polimi.it)
 * Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)
 * Lorenzo Di Tucci (lorenzo.ditucci@polimi.it)
 *
 * This file contains a custom implementation of QR factorization tailored
 * to a 5x9 input matrix
 *
 */

#ifndef STEP_2_H
#define STEP_2_H

#include "utils.h"
#include "qrUtils.h"

namespace hls {

/*
 * Computes the null space of a 5x9 matrix leveraging a custom QR factorization.
 * The output of the method is a 9x4 matrix in which each column is a column
 * vector of the null space of the matrix.
 */
  template <
          int RowsA,
          int ColsA,
    	  int RowsQ,
          typename InputStreamType,
          typename OutputType>
        void my_qrf_alt_top_9(hls::stream<InputStreamType> &in_stream,
                           hls::stream<OutputType> &out_stream, const int iter) {
#pragma HLS INLINE
    #pragma HLS DATAFLOW

          int sequence_0_3[4][3] = {
        		  {0, 1, 0},
    			  {2, 3, 0},
    			  {4, 5, 0},
    			  {6, 7, 0}
          };

          int sequence_4_7[4][3] = {
    			  {1, 3, 1},
    			  {5, 7, 1},
    			  {0, 8, 0},
    			  {2, 4, 0}
    	  };

          int sequence_8_11[4][3] = {
    			  {3, 7, 2},
    			  {1, 5, 1},
    			  {4, 8, 1},
    			  {0, 6, 0}
    	  };

          int sequence_12_14[3][3] = {
    			  {3, 5, 2},
    			  {1, 4, 1},
    			  {0, 2, 0}
    	  };

          int sequence_15_17[3][3] = {
    			  {5, 7, 3},
    			  {3, 8, 2},
    			  {1, 6, 1}
    	  };

          int sequence_18_20[3][3] = {
    			  {5, 8, 3},
    			  {3, 4, 2},
    			  {1, 2, 1}
    	  };

          int sequence_21_23[3][3] = {
    			  {7, 8, 4},
    			  {4, 5, 3},
    			  {3, 6, 2}
    	  };

          int sequence_24_26[3][3] = {
    			  {5, 7, 4},
    			  {4, 6, 3},
    			  {2, 3, 2}
    	  };

          int sequence_27_28[2][3] = {
    			  {5, 6, 4},
    			  {3, 4, 3}
    	  };

          int sequence_29[1][3] = {
    			  {4, 5, 4}
    	  };

          hls::stream<OutputType> q_stream[9];
#pragma HLS STREAM variable=q_stream depth=RowsA*RowsA*2 dim=1

          hls::stream<OutputType> r_stream[9];
#pragma HLS STREAM variable=r_stream depth=RowsA*ColsA*2 dim=1

          batch_first_ap<4, RowsA, ColsA, false, rType, OutputType>(in_stream, q_stream[0], r_stream[0], sequence_0_3, iter);
    	  batch<4, RowsA, ColsA, false, OutputType>(q_stream[0], r_stream[0], q_stream[1], r_stream[1], sequence_4_7, iter);
    	  batch<4, RowsA, ColsA, false, OutputType>(q_stream[1], r_stream[1], q_stream[2], r_stream[2], sequence_8_11, iter);
    	  batch<3, RowsA, ColsA, false, OutputType>(q_stream[2], r_stream[2], q_stream[3], r_stream[3], sequence_12_14, iter);
    	  batch<3, RowsA, ColsA, false, OutputType>(q_stream[3], r_stream[3], q_stream[4], r_stream[4], sequence_15_17, iter);
    	  batch<3, RowsA, ColsA, false, OutputType>(q_stream[4], r_stream[4], q_stream[5], r_stream[5], sequence_18_20, iter);
    	  batch<3, RowsA, ColsA, false, OutputType>(q_stream[5], r_stream[5], q_stream[6], r_stream[6], sequence_21_23, iter);
    	  batch<3, RowsA, ColsA, false, OutputType>(q_stream[6], r_stream[6], q_stream[7], r_stream[7], sequence_24_26, iter);
    	  batch<2, RowsA, ColsA, false, OutputType>(q_stream[7], r_stream[7], q_stream[8], r_stream[8], sequence_27_28, iter);
    	  batch_last_Q<1, RowsA, ColsA, RowsQ, true, OutputType>(q_stream[8], r_stream[8], out_stream, sequence_29, iter);

        }

}

#endif
