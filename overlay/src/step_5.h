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
 * step_5.h
 *
 * authors:
 * Marco Rabozzi (marco.rabozzi@polimi.it)
 * Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)
 * Lorenzo Di Tucci (lorenzo.ditucci@polimi.it)
 *
 * Perform step 5 of five-point algorithm:
 * solve the 10-degree polynomial and recover its roots.
 */

#ifndef STEP_5_H
#define STEP_5_H

#include "utils.h"
#include "polysolve.h"

/*
 * Solve the 10-th degree polynomial and find up to 10 roots.
 * NOTE: 10 values are returned and invalid roots are set to -INFINITY.
 */
void step_5_streaming(
		hls::stream<my_type> &coeff_stream,
		hls::stream<my_type> &out,
		const int iterations)
{
	#pragma HLS INLINE
	#pragma HLS DATAFLOW

	polysolve_streaming(coeff_stream, out, iterations);

}

#endif
