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
 * polysolve.h
 *
 * authors:
 * Marco Rabozzi (marco.rabozzi@polimi.it)
 * Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)
 * Lorenzo Di Tucci (lorenzo.ditucci@polimi.it)
 *
 *
 * This library provides functionalities to find the roots of a Nth-degree
 * polynomial. The main usage of this library is through the polysolve_streaming
 * function while the other functions are used in order to fulfill this task.
 *
 * The method is based on the recursive evaluation of the Sturm Sequence associated
 * with the polynomial (see: https://en.wikipedia.org/wiki/Sturm%27s_theorem) efficiently
 * computed as described in "An Efficient Solution to the Five-Point Relative Pose Problem"
 * from David Nister in 2003.
 */

#ifndef POLYSOLVE_H
#define POLYSOLVE_H

#include "hls_stream.h"
#include "ap_int.h"
#include "utils.h"
#include <stdint.h>
#include <math.h>

/*****************************************************************************
 * Constants
 ******************************************************************************/

// interval to search for roots: [-MAX_DISTANCE, +MAX_DISTANCE]
#define MAX_DISTANCE 4e3f

/* Max number of iterations. */

// the number of bisection iteration given by the product
// INTERNAL_ITERATIONS*MAX_ITERATIONS
// MAX_ITERATION dictates how many bisector module to use on each chain
// for the computation of a root, while INTERNAL_ITERAITONS is the number
// of bisections performed within each module of the chain.
//
// to reduce the initiation interval reduce INTERNAL_ITERATIONS and
// increase MAX_ITERATIONS.
// to reduce area at the cost of the initiation interval, increase
// INTERNAL_ITERATIONS and decrease MAX_ITERATIONS
#define INTERNAL_ITERATIONS 2
#define MAX_ITERATIONS 16

// order of the polynomial to solve
#define MAX_ORDER 10
#define MAX_ORDER_LOG2 4

// A coefficient smaller than SMALL_ENOUGH is considered to be zero (0.0)
#define SMALL_ENOUGH 1.0e-10

// Length of the Sturm sequence for counting number of roots in an interval
#define STURM_LEN (MAX_ORDER*2+2)
#define STURM_LEN_LOG2 (MAX_ORDER_LOG2 * 2 + ((MAX_ORDER + 2) >> MAX_ORDER_LOG2))

/*****************************************************************************
 * Typedefs
 ******************************************************************************/

// Data type used for streaming across bisection modules
// defines a Sturm Sequence
typedef struct {
	my_type val[STURM_LEN];
	ap_uint<STURM_LEN> mark;
} STURM_T;

// the state of the polynomial solver propagated between
// bisection stages.
typedef struct {
	my_type sturm_val[STURM_LEN];
	my_type min_value[MAX_ORDER];
	my_type max_value[MAX_ORDER];
	ap_uint<MAX_ORDER_LOG2> offset[MAX_ORDER];
	ap_uint<STURM_LEN> sturm_mark;
	ap_uint<MAX_ORDER_LOG2> atmin[MAX_ORDER];
	ap_uint<MAX_ORDER_LOG2> atmax[MAX_ORDER];
} BISECT_MULTI_T;

// utility union for float to raw bits conversion and vice-versa
typedef union {
	my_type f;
	uint32_t i;
} DATA_TYPE_U;

/*****************************************************************************
 * Functions
 ******************************************************************************/

/*
 * converts a my_type stream to a single STURM_T data type
 */
void data_stream_to_sturm_t(hls::stream<my_type> &in, STURM_T &out) {
#pragma HLS INLINE
	DATA_TYPE_U tmp;

	for (int i = 0; i < STURM_LEN; i++) {
#pragma HLS PIPELINE
		out.val[i] = in.read();
	}

	tmp.f = in.read();
	out.mark = ap_uint<STURM_LEN>(tmp.i);
}

/*
 * converts a STURM_T data type to a my_type stream
 */
void sturm_t_to_data_stream(STURM_T &in, hls::stream<my_type> &out) {
#pragma HLS INLINE
	DATA_TYPE_U tmp;

	for (int i = 0; i < STURM_LEN; i++) {
#pragma HLS PIPELINE
		out.write(in.val[i]);
	}
	tmp.i = in.mark.to_uint();
	out.write(tmp.f);
}

/*
 * converts a STURM_T data type to a my_type stream
 */
void data_stream_to_bisect_multi_t(hls::stream<my_type> &in,
		BISECT_MULTI_T &out) {
#pragma HLS INLINE
	DATA_TYPE_U tmp;

	for (int i = 0; i < STURM_LEN; i++) {
#pragma HLS PIPELINE
		out.sturm_val[i] = in.read();
	}

	tmp.f = in.read();
	out.sturm_mark = ap_uint<STURM_LEN>(tmp.i);
	for (int i = 0; i < MAX_ORDER; i++) {
#pragma HLS PIPELINE II=5
		out.min_value[i] = in.read();
		out.max_value[i] = in.read();
		tmp.f = in.read();
		out.offset[i] = ap_uint<MAX_ORDER_LOG2>(tmp.i);
		tmp.f = in.read();
		out.atmin[i] = ap_uint<MAX_ORDER_LOG2>(tmp.i);
		tmp.f = in.read();
		out.atmax[i] = ap_uint<MAX_ORDER_LOG2>(tmp.i);
	}
}

/*
 * converts a BISECT_MULTI_T data type to a my_type stream
 */
void bisect_multi_t_to_data_stream(BISECT_MULTI_T &in,
		hls::stream<my_type> &out) {
#pragma HLS INLINE
	DATA_TYPE_U tmp;

	for (int i = 0; i < STURM_LEN; i++) {
#pragma HLS PIPELINE
		out.write(in.sturm_val[i]);
	}
	tmp.i = in.sturm_mark.to_uint();
	out.write(tmp.f);
	for (int i = 0; i < MAX_ORDER; i++) {
#pragma HLS PIPELINE II=5
		out.write(in.min_value[i]);
		out.write(in.max_value[i]);
		tmp.i = in.offset[i].to_uint();
		out.write(tmp.f);
		tmp.i = in.atmin[i].to_uint();
		out.write(tmp.f);
		tmp.i = in.atmax[i].to_uint();
		out.write(tmp.f);
	}
}

/*
 * main component for building the Sturm Sequence.
 * It computes the values associated to the Sturm Sequence in the range start / end.
 */
void buildsturm_body(int ord, double num[MAX_ORDER + 1],
		double denom[MAX_ORDER + 1], int& ord_num, int& ord_denom,
		ap_uint<32>& sturm_mark_bits, my_type sturm_val[STURM_LEN], int& s,
		const int start, const int end) {

	#pragma HLS INLINE

	int tmp_i;
	double tmp;

	for (int i = start; i < end; i++) {

		if (ord_num > 0 && i < 2 * ord + 1) {

			if (fabs(num[ord_num]) <= SMALL_ENOUGH) {
				// put a zero in the sequence
				s--;
				ord_num--;
			} else if (ord > 1) {
				if (ord_denom > ord_num) {
					// swap and change sign of the computed remainder
					for (int j = 0; j < MAX_ORDER; j++) {
						#pragma HLS UNROLL
						tmp = num[j];
						num[j] = denom[j];
						denom[j] = -tmp;
					}

					tmp_i = ord_num;
					ord_num = ord_denom;
					ord_denom = tmp_i;

					sturm_mark_bits.set(s, 1);
				}

				double q = (double) num[ord_num] / (double) denom[ord_denom];

				// division step
				for (int j = 0; j <= MAX_ORDER; j++) {
					#pragma HLS UNROLL
					if (ord_num - ord_denom <= j && j <= ord_num)
						num[j] -= q * denom[j - (ord_num - ord_denom)];
				}

				// mark the end of a new polynomial of the sequence
				sturm_val[s] = q;
				s--;
				ord_num--;
			}
		}
	}
}

/*
 * The buildsturm operation is split in two stages within the pipeline
 * to reduce the dataflow initiation interval.
 *
 * Initialize the computation of the Sturm Sequence and computes the first
 * half of the values.
 */
void buildsturm_part1_inner(hls::stream<my_type> &coeffs_in,
		hls::stream<double> &num_denom_out,
		hls::stream<my_type> &partials_out) {

	int i, j;
	int ord, ord_num, ord_denom;
	int tmp_i;
	double tmp;

	my_type sturm_val[STURM_LEN];
	ap_uint<32> sturm_mark_bits = 0;

	double num[MAX_ORDER + 1];
#pragma HLS ARRAY_PARTITION variable=num complete dim=1
	double denom[MAX_ORDER + 1];
#pragma HLS ARRAY_PARTITION variable=denom complete dim=1
	my_type Coeffs[MAX_ORDER + 1];
#pragma HLS ARRAY_PARTITION variable=Coeffs complete dim=1

	for(int i = 0; i <= MAX_ORDER; i++) {
		#pragma HLS PIPELINE
		Coeffs[i] = coeffs_in.read();
	}

	// compute the order
	ord = 0;
	for (int i = MAX_ORDER; i >= 0; i--) {
		#pragma HLS UNROLL
		if (fabs(Coeffs[i]) > SMALL_ENOUGH) {
			ord = MAX_ORDER - i;
		}
	}

	// reset data structure
	for (int i = 0; i < STURM_LEN; i++) {
		#pragma HLS UNROLL
		sturm_val[i] = 0.0;
	}

	/* Put the coefficients into the top of the stack. */
	// prepare the first polynomial (numerator)
	for (int i = 0; i <= MAX_ORDER; i++) {
		#pragma HLS PIPELINE
		num[i] =
				(i <= ord) ?
						((double) (Coeffs[MAX_ORDER - i])
								/ (double) (Coeffs[MAX_ORDER - ord])) :
						0.0;
	}
	ord_num = ord;

	/* calculate the derivative for the second polynomial. */
	for (int i = 1; i <= MAX_ORDER; i++) {
		#pragma HLS PIPELINE
		denom[i - 1] = num[i] * i;
	}
	ord_denom = ord - 1;

	/* construct the rest of the Sturm sequence */

	int s = STURM_LEN - 1;

	sturm_mark_bits.set(s, 1); // end of the last polynomial in the sequence

	buildsturm_body(ord, num, denom, ord_num, ord_denom, sturm_mark_bits,
				sturm_val, s, 0, MAX_ORDER);

	// stream out num, denom, ord, ord_num, ord_denum, sturm_val, sturm_mark_bits
	DATA_TYPE_U tmp_u;

	for(int i = 0; i < MAX_ORDER + 1; i++) {
		#pragma HLS PIPELINE
		num_denom_out.write(num[i]);
	}
	for(int i = 0; i < MAX_ORDER + 1; i++) {
		#pragma HLS PIPELINE
		num_denom_out.write(denom[i]);
	}
	for(int i = 0; i < STURM_LEN; i++) {
		#pragma HLS PIPELINE
		partials_out.write(sturm_val[i]);
	}

	tmp_u.i = ord;
	partials_out.write(tmp_u.f);
	tmp_u.i = ord_num;
	partials_out.write(tmp_u.f);
	tmp_u.i = ord_denom;
	partials_out.write(tmp_u.f);
	tmp_u.i = s;
	partials_out.write(tmp_u.f);
	tmp_u.i = sturm_mark_bits.to_uint();
	partials_out.write(tmp_u.f);
}

/*
 * The buildsturm operation is split in two stages within the pipeline
 * to reduce the dataflow initiation interval.
 *
 * Computes the second half the Sturm Sequence values and finilize the sequence
 */
void buildsturm_part2_inner(hls::stream<double> &num_denom_in,
		hls::stream<my_type> &partials_in, hls::stream<my_type> &sturm_out) {

	int ord, ord_num, ord_denom;
	int tmp_i, s;
	double tmp;
	DATA_TYPE_U tmp_u;

	my_type sturm_val[STURM_LEN];

	ap_uint<32> sturm_mark_bits = 0;

	double num[MAX_ORDER + 1];
#pragma HLS ARRAY_PARTITION variable=num complete dim=1
	double denom[MAX_ORDER + 1];
#pragma HLS ARRAY_PARTITION variable=denom complete dim=1
	my_type Coeffs[MAX_ORDER + 1];
#pragma HLS ARRAY_PARTITION variable=Coeffs complete dim=1

	// read data from previous core
	for(int i = 0; i < MAX_ORDER + 1; i++) {
		#pragma HLS PIPELINE
		num[i] = num_denom_in.read();
	}
	for(int i = 0; i < MAX_ORDER + 1; i++) {
		#pragma HLS PIPELINE
		denom[i] = num_denom_in.read();
	}
	for(int i = 0; i < STURM_LEN; i++) {
		#pragma HLS PIPELINE
		sturm_val[i] = partials_in.read();
	}

	tmp_u.f = partials_in.read();
	ord = tmp_u.i;
	tmp_u.f = partials_in.read();
	ord_num = tmp_u.i;
	tmp_u.f = partials_in.read();
	ord_denom = tmp_u.i;
	tmp_u.f = partials_in.read();
	s = tmp_u.i;
	tmp_u.f = partials_in.read();
	sturm_mark_bits = ap_uint<32>(tmp_u.i);

	buildsturm_body(ord, num, denom, ord_num, ord_denom, sturm_mark_bits,
			sturm_val, s, MAX_ORDER, MAX_ORDER * 2 + 1);

	if (ord <= 1) {
		for (int i = 0; i <= 1; i++) {
			#pragma HLS UNROLL
			tmp = num[i];
			num[i] = denom[i];
			denom[i] = tmp;
		}

		tmp_i = ord_denom;
		ord_denom = ord_num;
		ord_num = tmp_i;
	}

	// add left over denominator and remainder
	sturm_mark_bits.set(s, 1);
	for (int i = ord_denom; i >= 0; i--) {
		#pragma HLS LOOP_TRIPCOUNT max=11
		#pragma HLS PIPELINE
		sturm_val[s] = denom[i];
		s--;
	}

	// add what's left in the last remainder of the sequence
	if (ord_num >= 0) {
		sturm_mark_bits.set(s, 1);
		for (int i = ord_num; i >= 0; i--) {
			#pragma HLS LOOP_TRIPCOUNT max=11
			#pragma HLS PIPELINE
			sturm_val[s] = -num[i];
			s--;
		}
	}

	// mark the end of initial empty sequence of 0s
	if (s >= 0) {
		sturm_mark_bits.set(s, 1);
	}

	for (int i = 0; i < STURM_LEN; i++) {
		#pragma HLS PIPELINE
		sturm_out.write(sturm_val[i]);
	}
	tmp_u.i = sturm_mark_bits.to_uint();
	sturm_out.write(tmp_u.f);
}

/*
 * Wrapper for the computation of the first part of the Sturm Sequence
 */
void buildsturm_part1(hls::stream<my_type> &coeffs_in,
		hls::stream<double> &num_denom_out,
		hls::stream<my_type> &partials_out,
		const int iterations) {
	for (int i = 0; i < iterations; i++) {
		buildsturm_part1_inner(coeffs_in, num_denom_out, partials_out);
	}
}

/*
 * Wrapper for the computation of the second part of the Sturm Sequence
 */
void buildsturm_part2(hls::stream<double> &num_denom_in,
		hls::stream<my_type> &partials_in, hls::stream<my_type> &sturm_out,
		const int iterations)
{
	for (int i = 0; i < iterations; i++) {
		buildsturm_part2_inner(num_denom_in, partials_in, sturm_out);
	}
}

/*
 * Computes the number of sign changes of the Sturm Sequence at a certain value
 */
ap_uint<MAX_ORDER_LOG2> numchanges(my_type sturm_val[STURM_LEN],
		ap_uint<STURM_LEN> sturm_mark, my_type a) {

#pragma HLS ALLOCATION instances=fmul limit=1 operation
#pragma HLS ALLOCATION instances=fadd limit=1 operation
#pragma HLS ALLOCATION instances=fsub limit=0 operation

	ap_uint<MAX_ORDER_LOG2> num_changes = -2; // it is guarantee to reach at least zero
	ap_uint<STURM_LEN_LOG2> polys = 1;
	my_type polyVal = 0.0f, fullVal;
	my_type curVal = 1.0f;
	my_type p0, p1, res;
	my_type val;
	ap_uint<2> f, lf = 2;

	for (int i = 0; i < STURM_LEN; i++) {
#pragma HLS PIPELINE
		val = sturm_val[i];

		polyVal += curVal * val;
		fullVal = polyVal * p1 + p0;
		curVal *= a;

		if (sturm_mark.get_bit(i)) {
			res = uchar_greater_3(polys) ? fullVal : polyVal;
			f = is_float_zero(res) ?
					ap_uint<2>(2) :
					(is_float_negative(res) ? ap_uint<2>(1) : ap_uint<2>(0));

			num_changes += (lf.get_bit(1) | (lf.get_bit(0) ^ f.get_bit(0)));

			lf = f;
			p0 = float_negate(p1);
			p1 = res;

			polyVal = 0.0f;
			curVal = 1.0f;
		}

		polys += sturm_mark.get_bit(i);
	}

	return num_changes;
}

/*
 * Initialize the the N intervals bisector.
 */
void bisect_init(my_type mid[MAX_ORDER], BISECT_MULTI_T& data,
		ap_uint<MAX_ORDER_LOG2> atmid[MAX_ORDER],
		ap_uint<STURM_LEN_LOG2> polys[MAX_ORDER], my_type polyVal[MAX_ORDER],
		my_type curVal[MAX_ORDER], my_type p0[MAX_ORDER], my_type p1[MAX_ORDER],
		ap_uint<2> lf[MAX_ORDER]) {
#pragma HLS INLINE
	for (int of = 0; of < MAX_ORDER; of++) {
#pragma HLS PIPELINE
		mid[of] = float_div2(data.min_value[of] + data.max_value[of]);

		atmid[of] = -2; // it is guarantee to reach at least zero
		polys[of] = 1;
		polyVal[of] = 0.0f;
		curVal[of] = 1.0f;
		p0[of] = 0.0f;
		p1[of] = 0.0f;
		lf[of] = 2;
	}
}

/*
 * Writes the result of the N intervals bisector.
 */
void bisect_write(ap_uint<MAX_ORDER_LOG2> diff, BISECT_MULTI_T& data,
		ap_uint<MAX_ORDER_LOG2> atmid[MAX_ORDER],
		bool check, my_type mid[MAX_ORDER]) {
#pragma HLS INLINE
	for (int of = 0; of < MAX_ORDER; of++) {
#pragma HLS PIPELINE

		diff = data.atmin[of] - atmid[of];
		check = data.offset[of] >= diff || diff == 0;
		data.min_value[of] = check ? mid[of] : data.min_value[of];
		data.atmin[of] = check ? atmid[of] : data.atmin[of];
		data.max_value[of] = check ? data.max_value[of] : mid[of];
		data.atmax[of] = check ? data.atmax[of] : atmid[of];
		data.offset[of] =
				check ? ap_uint<MAX_ORDER_LOG2>(data.offset[of] - diff) : data.offset[of];
	}
}

/*
 * Core computation of the N intervals bisector.
 */
void bisect_main_body(my_type val, BISECT_MULTI_T& data,
		my_type polyVal[MAX_ORDER], my_type curVal[MAX_ORDER], my_type fullVal,
		my_type p1[MAX_ORDER], my_type p0[MAX_ORDER], my_type mid[MAX_ORDER],
		my_type res, ap_uint<STURM_LEN_LOG2> polys[MAX_ORDER], ap_uint<2> f,
		ap_uint<MAX_ORDER_LOG2> atmid[MAX_ORDER], ap_uint<2> lf[MAX_ORDER]) {
#pragma HLS INLINE
	for (int i = 0; i < STURM_LEN; i++) {
		for (int of = 0; of < MAX_ORDER; of++) {
#pragma HLS PIPELINE

			val = data.sturm_val[i];

			polyVal[of] += curVal[of] * val;
			fullVal = polyVal[of] * p1[of] + p0[of];
			curVal[of] *= mid[of];

			if (data.sturm_mark.get_bit(i)) {
				res = uchar_greater_3(polys[of]) ? fullVal : polyVal[of];
				f =
						is_float_zero(res) ?
								ap_uint<2>(2) :
								(is_float_negative(res) ?
										ap_uint<2>(1) : ap_uint<2>(0));

				atmid[of] += (lf[of].get_bit(1)
						| (lf[of].get_bit(0) ^ f.get_bit(0)));

				lf[of] = f;
				p0[of] = float_negate(p1[of]);
				p1[of] = res;

				polyVal[of] = 0.0f;
				curVal[of] = 1.0f;
			}

			polys[of] += data.sturm_mark.get_bit(i);
		}
	}
}

/*
 * Bisect all the N intervals enclosing the N-roots of the polynomial.
 */
void bracket_multiple_root_stream_inner(hls::stream<my_type> &in,
		hls::stream<my_type> &out) {

	BISECT_MULTI_T data;

	ap_uint<MAX_ORDER_LOG2> atmid[MAX_ORDER];
	my_type mid[MAX_ORDER];

	ap_uint<MAX_ORDER_LOG2> num_changes[MAX_ORDER]; // it is guarantee to reach at least zero
	ap_uint<STURM_LEN_LOG2> polys[MAX_ORDER];
	my_type polyVal[MAX_ORDER];
	my_type curVal[MAX_ORDER];
	my_type p0[MAX_ORDER], p1[MAX_ORDER];
	my_type res, val, fullVal;
	ap_uint<2> f;
	ap_uint<2> lf[MAX_ORDER];

	ap_uint<MAX_ORDER_LOG2> diff;
	bool check;

	data_stream_to_bisect_multi_t(in, data);

	bisect_init(mid, data, atmid, polys, polyVal, curVal, p0, p1, lf);

	bisect_main_body(val, data, polyVal, curVal, fullVal, p1, p0, mid, res,
			polys, f, atmid, lf);

	bisect_write(diff, data, atmid, check, mid);

	bisect_multi_t_to_data_stream(data, out);
}

/*
 * Wrapper function for the N intervals bisector.
 */
void bracket_multiple_root_stream(hls::stream<my_type> &in,
		hls::stream<my_type> &out, const int iterations) {
	for (int i = 0; i < iterations; i++) {
		bracket_multiple_root_stream_inner(in, out);
	}
}

/*
 * Extracts the current estimation of the roots from the polynomial bisector state.
 * NOTE: if MAX_ORDER roots are searched but k < MAX_ORDER roots are found, the
 * remaining MAX_ORDER - k roots are set to -INFINITY to indicate that they are not
 * valid roots
 */
void extract_roots_inner(hls::stream<my_type> &in, hls::stream<my_type> &out) {
	BISECT_MULTI_T data;

	data_stream_to_bisect_multi_t(in, data);

	for (int of = 0; of < MAX_ORDER; of++) {
#pragma HLS PIPELINE
		my_type r =
				data.atmin[of] - data.atmax[of] > 0 ?
						float_div2(data.max_value[of] + data.min_value[of]) :
						-INFINITY;
		out.write(r);
	}
}

/*
 * Wrapper function for the root extractor.
 */
void extract_roots(hls::stream<my_type> &in, hls::stream<my_type> &out,
		const int iterations) {
	for (int i = 0; i < iterations; i++) {
		extract_roots_inner(in, out);
	}
}

/*
 * Computes the number of roots of a polynomial within an interval
 * and creates the bisector state.
 */
void roots_in_interval_inner(hls::stream<my_type> &in,
		hls::stream<my_type> &out) {
#pragma INLINE RECURSIVE

	STURM_T sturm;
	int i;

	my_type min_value = -MAX_DISTANCE;
	my_type max_value = MAX_DISTANCE;

	data_stream_to_sturm_t(in, sturm);

	ap_uint<MAX_ORDER_LOG2> atmin = numchanges(sturm.val, sturm.mark,
			min_value);
	ap_uint<MAX_ORDER_LOG2> atmax = numchanges(sturm.val, sturm.mark,
			max_value);

	DATA_TYPE_U tmp;

	for (i = 0; i < STURM_LEN; i++) {
#pragma HLS PIPELINE
		out.write(sturm.val[i]);
	}
	tmp.i = sturm.mark.to_uint();
	out.write(tmp.f);
	for (i = 0; i < MAX_ORDER; i++) {
#pragma HLS PIPELINE II=5
		out.write(min_value);
		out.write(max_value);
		tmp.i = i;
		out.write(tmp.f);
		tmp.i = atmin.to_uint();
		out.write(tmp.f);
		tmp.i = atmax.to_uint();
		out.write(tmp.f);
	}
}

/*
 * Wrapper function for finding the number of roots in a polynomial within an interval.
 */
void roots_in_interval(hls::stream<my_type> &in, hls::stream<my_type> &out,
		const int iterations) {
	for (int i = 0; i < iterations; i++) {
		roots_in_interval_inner(in, out);
	}
}

const int T5_STURM_S_DEPTH = (STURM_LEN + 1) * 2;
const int T5_BISECT_MULTI_S_DEPTH = (STURM_LEN + 1 + MAX_ORDER * 5) * 2;
const int T5_NUM_DENOM_S_DEPTH = (MAX_ORDER*2 + 2)*2;
const int T5_PARTIALS_S_DEPTH = (STURM_LEN + 5)*2;

/*
 * computes the roots of a polynomial of order MAX_ORDER given its coefficients.
 */
void polysolve_streaming(hls::stream<my_type> &in, hls::stream<my_type> &out,
		const int iterations) {
#pragma HLS INLINE
#pragma HLS DATAFLOW

	hls::stream<my_type> t5_sturm_s("t5_sturm_s");
#pragma HLS STREAM variable=t5_sturm_s depth=T5_STURM_S_DEPTH dim=1
	hls::stream<my_type> t5_bisect_multi_s[MAX_ITERATIONS * INTERNAL_ITERATIONS + 1];
#pragma HLS STREAM variable=t5_bisect_multi_s depth=T5_BISECT_MULTI_S_DEPTH dim=1
	hls::stream<double> num_denom_s("num_denom_s");
#pragma HLS STREAM variable=num_denom_s depth=T5_NUM_DENOM_S_DEPTH dim=1
	hls::stream<my_type> partials_s("partials_s");
#pragma HLS STREAM variable=partials_s depth=T5_PARTIALS_S_DEPTH dim=1

	// 5. generate sturm sequence
	buildsturm_part1(in, num_denom_s, partials_s, iterations);

	buildsturm_part2(num_denom_s, partials_s, t5_sturm_s, iterations);

	// 5. compute number of roots in the considered interval
	roots_in_interval(t5_sturm_s, t5_bisect_multi_s[0], iterations);

	// 5. chain of roots bisectors
	for (int i = 0; i < MAX_ITERATIONS * INTERNAL_ITERATIONS; i++) {
#pragma HLS UNROLL
		bracket_multiple_root_stream(t5_bisect_multi_s[i],
				t5_bisect_multi_s[i + 1], iterations);
	}

	// 5. extract roots from bisectors
	extract_roots(t5_bisect_multi_s[MAX_ITERATIONS * INTERNAL_ITERATIONS], out,
			iterations);
}

#endif
