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
 * qrUtils.h
 *
 * authors:
 * Marco Rabozzi (marco.rabozzi@polimi.it)
 * Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)
 * Lorenzo Di Tucci (lorenzo.ditucci@polimi.it)
 *
 * This file contains utilities to implement custom dataflow versions
 * of QR factorization
 */

#ifndef QR_UTILS_H
#define QR_UTILS_H

#include "utils.h"
#include "hls/linear_algebra/hls_qrf.h"

namespace hls {

/*
 * First QR batch
 * This batch supports as input a stream of ap_uint values
 */
  template<
  	  int iter,
	  int RowsA,
	  int ColsA,
	  bool last_batch,
	  typename streamType,
	  typename type>
  void batch_first_ap_inner(hls::stream<streamType> &r_stream_in, hls::stream<type> &q_stream_out, hls::stream<type> &r_stream_out, const int SEQUENCE[iter][3]){

	  	hls::stream<int> to_rot[3];
#pragma HLS STREAM variable=to_rot depth=RowsA/2
	  	int seq_cnt   = 0;
	  	int extra_pass  = 0;
		int extra_pass2 = 0;
		int use_mag     = 0;
		int px_row1, px_row2, px_col, rot_row1, rot_row2, rot_col;
		type G[2][2];
		type mag = 0;
		hls::stream<type> rotations[5];
		#pragma HLS STREAM variable=rotations depth=RowsA/2
		type G_delay[2][2];
		type mag_delay;

		type q_i[RowsA][RowsA] = {0.0};
		type r_i[RowsA][ColsA];
#pragma HLS ARRAY_PARTITION variable=r_i complete dim=1

		in_row_copy : for(int c=0; c<ColsA; c++){
#pragma HLS PIPELINE
			streamType r_val = r_stream_in.read();
			unsigned int tmp_r_0 = r_val.range(31, 0);
			unsigned int tmp_r_1 = r_val.range(63, 32);
			unsigned int tmp_r_2 = r_val.range(95, 64);
			unsigned int tmp_r_3 = r_val.range(127, 96);
			unsigned int tmp_r_4 = r_val.range(159, 128);
			unsigned int tmp_r_5 = r_val.range(191, 160);
			unsigned int tmp_r_6 = r_val.range(223, 192);
			unsigned int tmp_r_7 = r_val.range(255, 224);
			unsigned int tmp_r_8 = r_val.range(287, 256);

			r_i[0][c] = *((my_type *)&tmp_r_0);
			r_i[1][c] = *((my_type *)&tmp_r_1);
			r_i[2][c] = *((my_type *)&tmp_r_2);
			r_i[3][c] = *((my_type *)&tmp_r_3);
			r_i[4][c] = *((my_type *)&tmp_r_4);
			r_i[5][c] = *((my_type *)&tmp_r_5);
			r_i[6][c] = *((my_type *)&tmp_r_6);
			r_i[7][c] = *((my_type *)&tmp_r_7);
			r_i[8][c] = *((my_type *)&tmp_r_8);

		  }
		for(int r = 0; r < RowsA; r++){
#pragma HLS PIPELINE
			q_i[r][r] = 1.0;
		}

		  calc_rotations: for(int px_cnt = 0; px_cnt < iter; px_cnt++) {
			#pragma HLS PIPELINE II=1
			px_row1 = SEQUENCE[seq_cnt][0];
			px_row2 = SEQUENCE[seq_cnt][1];
			px_col  = SEQUENCE[seq_cnt][2];
			seq_cnt++;
			extra_pass = 0;
			qrf_givens(extra_pass, r_i[px_row1][px_col], r_i[px_row2][px_col], G[0][0], G[0][1], G[1][0], G[1][1], mag);
			// Pass on rotation to next block to apply rotations
			rotations[0].write(G[0][0]);
			rotations[1].write(G[0][1]);
			rotations[2].write(G[1][0]);
			rotations[3].write(G[1][1]);
			rotations[4].write(mag);
			to_rot[0].write(px_row1);
			to_rot[1].write(px_row2);
			to_rot[2].write(px_col);
		  }

		  rotate: for(int px_cnt = 0; px_cnt < iter; px_cnt++) {
			G_delay[0][0] = rotations[0].read();
			G_delay[0][1] = rotations[1].read();
			G_delay[1][0] = rotations[2].read();
			G_delay[1][1] = rotations[3].read();
			mag_delay     = rotations[4].read();
			rot_row1      = to_rot[0].read();
			rot_row2      = to_rot[1].read();
			rot_col       = to_rot[2].read();

			extra_pass2 = 0;

			// Merge the loops to maximize throughput, otherwise HLS will execute them sequentially and
			// share hardware.
			#pragma HLS LOOP_MERGE force
			update_r : for(int k=0; k<ColsA; k++) {
			  #pragma HLS PIPELINE II=4
			  #pragma HLS UNROLL FACTOR=1
			  use_mag = 0;
			  if (k==rot_col) {
				use_mag = 1;
			  }
			  qrf_mm_or_mag(G_delay, r_i[rot_row1][k], r_i[rot_row2][k], mag_delay, use_mag, extra_pass2);
			}
			update_q : for(int k=0; k<RowsA; k++) {
	#pragma HLS PIPELINE II=1
			  #pragma HLS PIPELINE II=4
			  #pragma HLS UNROLL FACTOR=1
			  qrf_mm(G_delay, q_i[rot_row1][k], q_i[rot_row2][k]);
			}
		  }

		  out_row_copy : for(int r=0; r<RowsA; r++){
			// Merge loops to parallelize the A input read and the Q matrix prime.
			#pragma HLS LOOP_MERGE force
			out_col_copy_q_i : for(int c=0; c<RowsA; c++) {
			  #pragma HLS PIPELINE
			  q_stream_out << q_i[r][c];
			}
			out_col_copy_r_i : for(int c=0; c<ColsA; c++) {
			  #pragma HLS PIPELINE
				r_stream_out << r_i[r][c];
			}
		  }

  }

/*
 * First QR batch wrapper
 * This wrapper supports as input a stream of ap_uint values
 */
  template<
      int batchIter,
  	  int RowsA,
  	  int ColsA,
  	  bool last_batch,
	  typename streamType,
  	  typename type>
    void batch_first_ap(hls::stream<streamType> &r_stream_in, hls::stream<type> &q_stream_out, hls::stream<type> &r_stream_out, const int SEQUENCE[batchIter][3], const int iter){
	  	  for(int i = 0; i < iter; i++){
	  		batch_first_ap_inner<batchIter, RowsA, ColsA, last_batch, streamType, type>(r_stream_in, q_stream_out, r_stream_out, SEQUENCE);
	  	  }
  }

  /*
   * QR batch
   */
  template<
  	  int iter,
	  int RowsA,
	  int ColsA,
	  bool last_batch,
	  typename type>
  void batch_inner(hls::stream<type> &q_stream_in, hls::stream<type> &r_stream_in, hls::stream<type> &q_stream_out, hls::stream<type> &r_stream_out, const int SEQUENCE[iter][3]){

	  	hls::stream<int> to_rot[3];
#pragma HLS STREAM variable=to_rot depth=RowsA/2
	  	int seq_cnt   = 0;
	  	int extra_pass  = 0;
		int extra_pass2 = 0;
		int use_mag     = 0;
		int px_row1, px_row2, px_col, rot_row1, rot_row2, rot_col;
		type G[2][2];
		type mag = 0;
		hls::stream<type> rotations[5];
		#pragma HLS STREAM variable=rotations depth=RowsA/2
		type G_delay[2][2];
		type mag_delay;

		type q_i[RowsA][RowsA];
		type r_i[RowsA][ColsA];

		in_row_copy : for(int r=0; r<RowsA; r++){
			// Merge loops to parallelize the A input read and the Q matrix prime.
			#pragma HLS LOOP_MERGE force
			in_col_copy_q_i : for(int c=0; c<RowsA; c++) {
			  #pragma HLS PIPELINE
			  q_i[r][c] = q_stream_in.read();
			}
			in_col_copy_r_i : for(int c=0; c<ColsA; c++) {
			  #pragma HLS PIPELINE
			  r_i[r][c] = r_stream_in.read();
			}
		  }

		  calc_rotations: for(int px_cnt = 0; px_cnt < iter; px_cnt++) {
			#pragma HLS PIPELINE II=1
			px_row1 = SEQUENCE[seq_cnt][0];
			px_row2 = SEQUENCE[seq_cnt][1];
			px_col  = SEQUENCE[seq_cnt][2];
			seq_cnt++;
			extra_pass = 0;
			qrf_givens(extra_pass, r_i[px_row1][px_col], r_i[px_row2][px_col], G[0][0], G[0][1], G[1][0], G[1][1], mag);
			// Pass on rotation to next block to apply rotations
			rotations[0].write(G[0][0]);
			rotations[1].write(G[0][1]);
			rotations[2].write(G[1][0]);
			rotations[3].write(G[1][1]);
			rotations[4].write(mag);
			to_rot[0].write(px_row1);
			to_rot[1].write(px_row2);
			to_rot[2].write(px_col);
		  }

		  rotate: for(int px_cnt = 0; px_cnt < iter; px_cnt++) {
			G_delay[0][0] = rotations[0].read();
			G_delay[0][1] = rotations[1].read();
			G_delay[1][0] = rotations[2].read();
			G_delay[1][1] = rotations[3].read();
			mag_delay     = rotations[4].read();
			rot_row1      = to_rot[0].read();
			rot_row2      = to_rot[1].read();
			rot_col       = to_rot[2].read();

			extra_pass2 = 0;

			// Merge the loops to maximize throughput, otherwise HLS will execute them sequentially and
			// share hardware.
			#pragma HLS LOOP_MERGE force
			update_r : for(int k=0; k<ColsA; k++) {
			  #pragma HLS PIPELINE II=4
			  #pragma HLS UNROLL FACTOR=1
			  use_mag = 0;
			  if (k==rot_col) {
				use_mag = 1;
			  }
			  qrf_mm_or_mag(G_delay, r_i[rot_row1][k], r_i[rot_row2][k], mag_delay, use_mag, extra_pass2);
			}
			update_q : for(int k=0; k<RowsA; k++) {
	#pragma HLS PIPELINE II=1
			  #pragma HLS PIPELINE II=4
			  #pragma HLS UNROLL FACTOR=1
			  qrf_mm(G_delay, q_i[rot_row1][k], q_i[rot_row2][k]);
			}
		  }

		  out_row_copy : for(int r=0; r<RowsA; r++){
			// Merge loops to parallelize the A input read and the Q matrix prime.
			#pragma HLS LOOP_MERGE force
			out_col_copy_q_i : for(int c=0; c<RowsA; c++) {
			  #pragma HLS PIPELINE
			  q_stream_out << q_i[r][c];
			}
			out_col_copy_r_i : for(int c=0; c<ColsA; c++) {
			  #pragma HLS PIPELINE
				r_stream_out << r_i[r][c];
			}
		  }

  }

  /*
   * QR batch wrapper
   */
  template<
    	  int batchIter,
  	  int RowsA,
  	  int ColsA,
  	  bool last_batch,
  	  typename type>
   void batch(hls::stream<type> &q_stream_in, hls::stream<type> &r_stream_in, hls::stream<type> &q_stream_out, hls::stream<type> &r_stream_out, const int SEQUENCE[batchIter][3], const int iter){
	  	  for(int i = 0; i < iter; i++){
	  		batch_inner<batchIter, RowsA, ColsA, last_batch, type>(q_stream_in, r_stream_in, q_stream_out, r_stream_out, SEQUENCE);
	  	  }
  }

  /*
   * QR last batch
   */
  template<
	  int iter,
	  int RowsA,
	  int ColsA,
  int RowsQ,
	  bool last_batch,
	  typename type>
  void batch_last_Q_inner(hls::stream<type> &q_stream_in, hls::stream<type> &r_stream_in, hls::stream<type> &q_stream_out, const int SEQUENCE[iter][3]){

		hls::stream<int> to_rot[3];
#pragma HLS STREAM variable=to_rot depth=RowsA/2
		int seq_cnt   = 0;
		int extra_pass  = 0;
		int extra_pass2 = 0;
		int use_mag     = 0;
		int px_row1, px_row2, px_col, rot_row1, rot_row2, rot_col;
		type G[2][2];
		type mag = 0;
		hls::stream<type> rotations[5];
		#pragma HLS STREAM variable=rotations depth=RowsA/2
		type G_delay[2][2];
		type mag_delay;

		type q_i[RowsA][RowsA];
		type r_i[RowsA][ColsA];

		in_row_copy : for(int r=0; r<RowsA; r++){
			// Merge loops to parallelize the A input read and the Q matrix prime.
			#pragma HLS LOOP_MERGE force
			in_col_copy_q_i : for(int c=0; c<RowsA; c++) {
			  #pragma HLS PIPELINE
			  q_i[r][c] = q_stream_in.read();
			}
			in_col_copy_r_i : for(int c=0; c<ColsA; c++) {
			  #pragma HLS PIPELINE
			  r_i[r][c] = r_stream_in.read();
			}
		  }

		  calc_rotations: for(int px_cnt = 0; px_cnt < iter; px_cnt++) {
			#pragma HLS PIPELINE II=1
			px_row1 = SEQUENCE[seq_cnt][0];
			px_row2 = SEQUENCE[seq_cnt][1];
			px_col  = SEQUENCE[seq_cnt][2];
			seq_cnt++;
			extra_pass = 0;
			qrf_givens(extra_pass, r_i[px_row1][px_col], r_i[px_row2][px_col], G[0][0], G[0][1], G[1][0], G[1][1], mag);
			// Pass on rotation to next block to apply rotations
			rotations[0].write(G[0][0]);
			rotations[1].write(G[0][1]);
			rotations[2].write(G[1][0]);
			rotations[3].write(G[1][1]);
			rotations[4].write(mag);
			to_rot[0].write(px_row1);
			to_rot[1].write(px_row2);
			to_rot[2].write(px_col);
		  }

		  rotate: for(int px_cnt = 0; px_cnt < iter; px_cnt++) {
			G_delay[0][0] = rotations[0].read();
			G_delay[0][1] = rotations[1].read();
			G_delay[1][0] = rotations[2].read();
			G_delay[1][1] = rotations[3].read();
			mag_delay     = rotations[4].read();
			rot_row1      = to_rot[0].read();
			rot_row2      = to_rot[1].read();
			rot_col       = to_rot[2].read();

			extra_pass2 = 0;

			// Merge the loops to maximize throughput, otherwise HLS will execute them sequentially and
			// share hardware.
			#pragma HLS LOOP_MERGE force
			update_r : for(int k=0; k<ColsA; k++) {
			  #pragma HLS PIPELINE II=4
			  #pragma HLS UNROLL FACTOR=1
			  use_mag = 0;
			  if (k==rot_col) {
				use_mag = 1;
			  }
			  qrf_mm_or_mag(G_delay, r_i[rot_row1][k], r_i[rot_row2][k], mag_delay, use_mag, extra_pass2);
			}
			update_q : for(int k=0; k<RowsA; k++) {
	#pragma HLS PIPELINE II=1
			  #pragma HLS PIPELINE II=4
			  #pragma HLS UNROLL FACTOR=1
			  qrf_mm(G_delay, q_i[rot_row1][k], q_i[rot_row2][k]);
			}
		  }

		  out_row_copy : for(int r=RowsA-RowsQ; r<RowsA; r++){
			// Merge loops to parallelize the A input read and the Q matrix prime.
			out_col_copy_q_i : for(int c=0; c<RowsA; c++) {
			  #pragma HLS PIPELINE
			  q_stream_out << q_i[r][c];
			}
		  }

  }

  /*
   * QR last batch wrapper
   */
  template<
      	  int batchIter,
    	  int RowsA,
    	  int ColsA,
		  int RowsQ,
    	  bool last_batch,
    	  typename type>
      void batch_last_Q(hls::stream<type> &q_stream_in, hls::stream<type> &r_stream_in, hls::stream<type> &q_stream_out, const int SEQUENCE[batchIter][3], const int iter){
  	  	  for(int i = 0; i < iter; i++){
  	  		batch_last_Q_inner<batchIter, RowsA, ColsA, RowsQ, last_batch, type>(q_stream_in, r_stream_in, q_stream_out, SEQUENCE);
  	  	  }
    }

  /*
   * Loop of batches
   */
  template<
  	  int iter,
  	  int RowsA,
  	  int ColsA,
  	  bool last_batch,
  	  typename type>
  void batch_loop_inner(hls::stream<type> &q_stream_in, hls::stream<type> &r_stream_in, hls::stream<type> &q_stream_out, hls::stream<type> &r_stream_out, const int SEQUENCE[iter][3]){

  	for(int it = 0; it < ROOT_NUM; it++) {
  		#pragma HLS PIPELINE

  	  	hls::stream<int> to_rot[3];
  #pragma HLS STREAM variable=to_rot depth=RowsA/2
  	  	int seq_cnt   = 0;
  	  	int extra_pass  = 0;
  		int extra_pass2 = 0;
  		int use_mag     = 0;
  		int px_row1, px_row2, px_col, rot_row1, rot_row2, rot_col;
  		type G[2][2];
  		type mag = 0;
  		hls::stream<type> rotations[5];
  		#pragma HLS STREAM variable=rotations depth=RowsA/2
  		type G_delay[2][2];
  		type mag_delay;

  		type q_i[RowsA][RowsA];
  		type r_i[RowsA][ColsA];

  		in_row_copy : for(int r=0; r<RowsA; r++){
  			// Merge loops to parallelize the A input read and the Q matrix prime.
  			#pragma HLS LOOP_MERGE force
  			in_col_copy_q_i : for(int c=0; c<RowsA; c++) {
  			  #pragma HLS PIPELINE
  			  q_i[r][c] = q_stream_in.read();
  			}
  			in_col_copy_r_i : for(int c=0; c<ColsA; c++) {
  			  #pragma HLS PIPELINE
  			  r_i[r][c] = r_stream_in.read();
  			}
  		  }

  		  calc_rotations: for(int px_cnt = 0; px_cnt < iter; px_cnt++) {
  			#pragma HLS PIPELINE II=1
  			px_row1 = SEQUENCE[seq_cnt][0];
  			px_row2 = SEQUENCE[seq_cnt][1];
  			px_col  = SEQUENCE[seq_cnt][2];
  			seq_cnt++;
  			extra_pass = 0;
  			qrf_givens(extra_pass, r_i[px_row1][px_col], r_i[px_row2][px_col], G[0][0], G[0][1], G[1][0], G[1][1], mag);
  			// Pass on rotation to next block to apply rotations
  			rotations[0].write(G[0][0]);
  			rotations[1].write(G[0][1]);
  			rotations[2].write(G[1][0]);
  			rotations[3].write(G[1][1]);
  			rotations[4].write(mag);
  			to_rot[0].write(px_row1);
  			to_rot[1].write(px_row2);
  			to_rot[2].write(px_col);
  		  }

  		  rotate: for(int px_cnt = 0; px_cnt < iter; px_cnt++) {
  			G_delay[0][0] = rotations[0].read();
  			G_delay[0][1] = rotations[1].read();
  			G_delay[1][0] = rotations[2].read();
  			G_delay[1][1] = rotations[3].read();
  			mag_delay     = rotations[4].read();
  			rot_row1      = to_rot[0].read();
  			rot_row2      = to_rot[1].read();
  			rot_col       = to_rot[2].read();

  			extra_pass2 = 0;

  			// Merge the loops to maximize throughput, otherwise HLS will execute them sequentially and
  			// share hardware.
  			#pragma HLS LOOP_MERGE force
  			update_r : for(int k=0; k<ColsA; k++) {
  			  #pragma HLS PIPELINE II=4
  			  #pragma HLS UNROLL FACTOR=1
  			  use_mag = 0;
  			  if (k==rot_col) {
  				use_mag = 1;
  			  }
  			  qrf_mm_or_mag(G_delay, r_i[rot_row1][k], r_i[rot_row2][k], mag_delay, use_mag, extra_pass2);
  			}
  			update_q : for(int k=0; k<RowsA; k++) {
  	#pragma HLS PIPELINE II=1
  			  #pragma HLS PIPELINE II=4
  			  #pragma HLS UNROLL FACTOR=1
  			  qrf_mm(G_delay, q_i[rot_row1][k], q_i[rot_row2][k]);
  			}
  		  }

  		  out_row_copy : for(int r=0; r<RowsA; r++){
  			// Merge loops to parallelize the A input read and the Q matrix prime.
  			#pragma HLS LOOP_MERGE force
  			out_col_copy_q_i : for(int c=0; c<RowsA; c++) {
  			  #pragma HLS PIPELINE
  			  q_stream_out << q_i[r][c];
  			}
  			out_col_copy_r_i : for(int c=0; c<ColsA; c++) {
  			  #pragma HLS PIPELINE
  				r_stream_out << r_i[r][c];
  			}
  		  }
  	}
  }

  /*
   * Wrapper for a loop of batches
   */
  template<
    	  int batchIter,
  	  int RowsA,
  	  int ColsA,
  	  bool last_batch,
  	  typename type>
   void batch_loop(hls::stream<type> &q_stream_in, hls::stream<type> &r_stream_in, hls::stream<type> &q_stream_out, hls::stream<type> &r_stream_out, const int SEQUENCE[batchIter][3], const int iter){
  	  	  for(int i = 0; i < iter; i++){
  	  		batch_loop_inner<batchIter, RowsA, ColsA, last_batch, type>(q_stream_in, r_stream_in, q_stream_out, r_stream_out, SEQUENCE);
  	  	  }
  }

  /*
   * Last QR batch
   * This version of the last QR batch returns the last column of Q matrix
   */
  template<
  	  int iter,
  	  int RowsA,
  	  int ColsA,
  int RowsQ,
  	  bool last_batch,
  	  typename type>
  void batch_last_Q_last_col_inner(hls::stream<type> &q_stream_in, hls::stream<type> &r_stream_in, hls::stream<type> &q_stream_out, const int SEQUENCE[iter][3]){

  	for(int it = 0; it < ROOT_NUM; it++) {
  		#pragma HLS PIPELINE

  		hls::stream<int> to_rot[3];
  #pragma HLS STREAM variable=to_rot depth=RowsA/2
  		int seq_cnt   = 0;
  		int extra_pass  = 0;
  		int extra_pass2 = 0;
  		int use_mag     = 0;
  		int px_row1, px_row2, px_col, rot_row1, rot_row2, rot_col;
  		type G[2][2];
  		type mag = 0;
  		hls::stream<type> rotations[5];
  		#pragma HLS STREAM variable=rotations depth=RowsA/2
  		type G_delay[2][2];
  		type mag_delay;

  		type q_i[RowsA][RowsA];
  		type r_i[RowsA][ColsA];

  		in_row_copy : for(int r=0; r<RowsA; r++){
  			// Merge loops to parallelize the A input read and the Q matrix prime.
  			#pragma HLS LOOP_MERGE force
  			in_col_copy_q_i : for(int c=0; c<RowsA; c++) {
  			  #pragma HLS PIPELINE
  			  q_i[r][c] = q_stream_in.read();
  			}
  			in_col_copy_r_i : for(int c=0; c<ColsA; c++) {
  			  #pragma HLS PIPELINE
  			  r_i[r][c] = r_stream_in.read();
  			}
  		  }

  		  calc_rotations: for(int px_cnt = 0; px_cnt < iter; px_cnt++) {
  			#pragma HLS PIPELINE II=1
  			px_row1 = SEQUENCE[seq_cnt][0];
  			px_row2 = SEQUENCE[seq_cnt][1];
  			px_col  = SEQUENCE[seq_cnt][2];
  			seq_cnt++;
  			extra_pass = 0;
  			qrf_givens(extra_pass, r_i[px_row1][px_col], r_i[px_row2][px_col], G[0][0], G[0][1], G[1][0], G[1][1], mag);
  			// Pass on rotation to next block to apply rotations
  			rotations[0].write(G[0][0]);
  			rotations[1].write(G[0][1]);
  			rotations[2].write(G[1][0]);
  			rotations[3].write(G[1][1]);
  			rotations[4].write(mag);
  			to_rot[0].write(px_row1);
  			to_rot[1].write(px_row2);
  			to_rot[2].write(px_col);
  		  }

  		  rotate: for(int px_cnt = 0; px_cnt < iter; px_cnt++) {
  			G_delay[0][0] = rotations[0].read();
  			G_delay[0][1] = rotations[1].read();
  			G_delay[1][0] = rotations[2].read();
  			G_delay[1][1] = rotations[3].read();
  			mag_delay     = rotations[4].read();
  			rot_row1      = to_rot[0].read();
  			rot_row2      = to_rot[1].read();
  			rot_col       = to_rot[2].read();

  			extra_pass2 = 0;

  			// Merge the loops to maximize throughput, otherwise HLS will execute them sequentially and
  			// share hardware.
  			#pragma HLS LOOP_MERGE force
  			update_r : for(int k=0; k<ColsA; k++) {
  			  #pragma HLS PIPELINE II=4
  			  #pragma HLS UNROLL FACTOR=1
  			  use_mag = 0;
  			  if (k==rot_col) {
  				use_mag = 1;
  			  }
  			  qrf_mm_or_mag(G_delay, r_i[rot_row1][k], r_i[rot_row2][k], mag_delay, use_mag, extra_pass2);
  			}
  			update_q : for(int k=0; k<RowsA; k++) {
  	#pragma HLS PIPELINE II=1
  			  #pragma HLS PIPELINE II=4
  			  #pragma HLS UNROLL FACTOR=1
  			  qrf_mm(G_delay, q_i[rot_row1][k], q_i[rot_row2][k]);
  			}
  		  }

  		  out_row_copy : for(int r=0; r<RowsA; r++){
  			// Merge loops to parallelize the A input read and the Q matrix prime.
  			  #pragma HLS PIPELINE
  			  q_stream_out << q_i[RowsA-1][r];
  		  }
  	}
  }

  /*
   * Last QR batch wrapper
   */
  template<
      	  int batchIter,
    	  int RowsA,
    	  int ColsA,
  		  int RowsQ,
    	  bool last_batch,
    	  typename type>
      void batch_last_Q_last_col(hls::stream<type> &q_stream_in, hls::stream<type> &r_stream_in, hls::stream<type> &q_stream_out, const int SEQUENCE[batchIter][3], const int iter){
  	  	  for(int i = 0; i < iter; i++){
  	  		batch_last_Q_last_col_inner<batchIter, RowsA, ColsA, RowsQ, last_batch, type>(q_stream_in, r_stream_in, q_stream_out, SEQUENCE);
  	  	  }
    }

  /*
   * QR batch trR
   * This version of the batch transposes matrix R
   */
  template<
  	  int iter,
  	  int RowsA,
  	  int ColsA,
  	  bool last_batch,
  	  typename type>
  void batch_trR_inner(hls::stream<type> &q_stream_in, hls::stream<type> &r_stream_in, hls::stream<type> &q_stream_out, hls::stream<type> &r_stream_out, const int SEQUENCE[iter][3]){

  	for(int it = 0; it < ROOT_NUM; it++) {
  		#pragma HLS PIPELINE

  	  	hls::stream<int> to_rot[3];
  #pragma HLS STREAM variable=to_rot depth=RowsA/2
  	  	int seq_cnt   = 0;
  	  	int extra_pass  = 0;
  		int extra_pass2 = 0;
  		int use_mag     = 0;
  		int px_row1, px_row2, px_col, rot_row1, rot_row2, rot_col;
  		type G[2][2];
  		type mag = 0;
  		hls::stream<type> rotations[5];
  		#pragma HLS STREAM variable=rotations depth=RowsA/2
  		type G_delay[2][2];
  		type mag_delay;

  		type q_i[RowsA][RowsA];
  		type r_i[RowsA][ColsA];

  		in_row_copy : for(int r=0; r<RowsA; r++){
  			// Merge loops to parallelize the A input read and the Q matrix prime.
  			#pragma HLS LOOP_MERGE force
  			in_col_copy_q_i : for(int c=0; c<RowsA; c++) {
  			  #pragma HLS PIPELINE
  			  q_i[r][c] = q_stream_in.read();
  			}
  			in_col_copy_r_i : for(int c=0; c<ColsA; c++) {
  			  #pragma HLS PIPELINE
  			  r_i[c][r] = r_stream_in.read();
  			}
  		  }

  		  calc_rotations: for(int px_cnt = 0; px_cnt < iter; px_cnt++) {
  			#pragma HLS PIPELINE II=1
  			px_row1 = SEQUENCE[seq_cnt][0];
  			px_row2 = SEQUENCE[seq_cnt][1];
  			px_col  = SEQUENCE[seq_cnt][2];
  			seq_cnt++;
  			extra_pass = 0;
  			qrf_givens(extra_pass, r_i[px_row1][px_col], r_i[px_row2][px_col], G[0][0], G[0][1], G[1][0], G[1][1], mag);
  			// Pass on rotation to next block to apply rotations
  			rotations[0].write(G[0][0]);
  			rotations[1].write(G[0][1]);
  			rotations[2].write(G[1][0]);
  			rotations[3].write(G[1][1]);
  			rotations[4].write(mag);
  			to_rot[0].write(px_row1);
  			to_rot[1].write(px_row2);
  			to_rot[2].write(px_col);
  		  }

  		  rotate: for(int px_cnt = 0; px_cnt < iter; px_cnt++) {
  			G_delay[0][0] = rotations[0].read();
  			G_delay[0][1] = rotations[1].read();
  			G_delay[1][0] = rotations[2].read();
  			G_delay[1][1] = rotations[3].read();
  			mag_delay     = rotations[4].read();
  			rot_row1      = to_rot[0].read();
  			rot_row2      = to_rot[1].read();
  			rot_col       = to_rot[2].read();

  			extra_pass2 = 0;

  			// Merge the loops to maximize throughput, otherwise HLS will execute them sequentially and
  			// share hardware.
  			#pragma HLS LOOP_MERGE force
  			update_r : for(int k=0; k<ColsA; k++) {
  			  #pragma HLS PIPELINE II=4
  			  #pragma HLS UNROLL FACTOR=1
  			  use_mag = 0;
  			  if (k==rot_col) {
  				use_mag = 1;
  			  }
  			  qrf_mm_or_mag(G_delay, r_i[rot_row1][k], r_i[rot_row2][k], mag_delay, use_mag, extra_pass2);
  			}
  			update_q : for(int k=0; k<RowsA; k++) {
  	#pragma HLS PIPELINE II=1
  			  #pragma HLS PIPELINE II=4
  			  #pragma HLS UNROLL FACTOR=1
  			  qrf_mm(G_delay, q_i[rot_row1][k], q_i[rot_row2][k]);
  			}
  		  }

  		  out_row_copy : for(int r=0; r<RowsA; r++){
  			// Merge loops to parallelize the A input read and the Q matrix prime.
  			#pragma HLS LOOP_MERGE force
  			out_col_copy_q_i : for(int c=0; c<RowsA; c++) {
  			  #pragma HLS PIPELINE
  			  q_stream_out << q_i[r][c];
  			}
  			out_col_copy_r_i : for(int c=0; c<ColsA; c++) {
  			  #pragma HLS PIPELINE
  				r_stream_out << r_i[r][c];
  			}
  		  }
  	}
  }

  /*
   * QR batch trR wrapper
   */
  template<
    	  int batchIter,
  	  int RowsA,
  	  int ColsA,
  	  bool last_batch,
  	  typename type>
   void batch_trR(hls::stream<type> &q_stream_in, hls::stream<type> &r_stream_in, hls::stream<type> &q_stream_out, hls::stream<type> &r_stream_out, const int SEQUENCE[batchIter][3], const int iter){
  	  	  for(int i = 0; i < iter; i++){
  	  		batch_trR_inner<batchIter, RowsA, ColsA, last_batch, type>(q_stream_in, r_stream_in, q_stream_out, r_stream_out, SEQUENCE);
  	  	  }
  }

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
  	  batch_trR<1, RowsA, ColsA, false, Type>(q_stream[0], r_stream[0], q_stream[1], r_stream[1], sequence_0, iter);
  	  batch_loop<1, RowsA, ColsA, false, Type>(q_stream[1], r_stream[1], q_stream[2], r_stream[2], sequence_1, iter);
  	  batch_last_Q_last_col<1, RowsA, ColsA, RowsA, true, Type>(q_stream[2], r_stream[2], out_streamQ, sequence_2, iter);

    }
}




#endif
