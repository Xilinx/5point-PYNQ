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
 * memInterface.h
 *
 * authors:
 * Marco Rabozzi (marco.rabozzi@polimi.it)
 * Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)
 * Lorenzo Di Tucci (lorenzo.ditucci@polimi.it)
 *
 *
 * This file contains a set of utility to read/write data from/to memory,
 * and push/pop to/from hls::stream
 * The purpose of wrapper functions is to invoke N times the inner function,
 * in order to enable overlapping computations among different couples of
 * 5 points
 *
 */

#ifndef MEM_INTERFACE_H
#define MEM_INTERFACE_H

#include "utils.h"

/*
 * hls::stream to master axi
 */
template<
	typename Type,
	int size>
void stream2axi_inner(hls::stream<Type> &in, Type* out){
	s2a_loop_0:for(int i = 0; i < size; i++){
#pragma HLS PIPELINE
		Type val = in.read();
		out[i] = val;
	}
}

/*
 * hls::stream to master axi wrapper
 */
template<
	typename Type,
	int size>
void stream2axi(hls::stream<Type> &in, Type* out, const int iter){
	s2a_loop_0:for(int i = 0; i < size*iter; i++){
#pragma HLS PIPELINE
		Type val = in.read();
		out[i] = val;
	}
}

/*
 * Master axi to hls::stream
 */
template<
	typename Type,
	int size>
void axi2stream_inner(Type* in, hls::stream<Type> &out){
	a2s_loop_0:for(int i = 0; i < size; i++){
#pragma HLS PIPELINE
		out.write(in[i]);
	}
}

/*
 * Master axi to hls::stream wrapper
 */
template<
	typename Type,
	int size>
void axi2stream(Type* in, hls::stream<Type> &out, const int iter){
	a2s_loop_0:for(int i = 0; i < size*iter; i++){
#pragma HLS PIPELINE
		out.write(in[i]);
	}

}

#endif
