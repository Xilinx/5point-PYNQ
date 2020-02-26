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
 * testbench.cpp
 *
 * authors:
 * Marco Rabozzi (marco.rabozzi@polimi.it)
 * Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)
 * Lorenzo Di Tucci (lorenzo.ditucci@polimi.it)
 *
 *
 *
 * testbench file to evaluate the correctness of the computation
 *
 *********
 * Modifications:
 *
 * 01/31/2020 gn - changed kernel function name to "fivept_nister"
 *
 *********/

#include "ap_int.h"
typedef ap_uint<64> point;
#define OUT_AP 512
typedef ap_uint<OUT_AP> outType;
typedef float my_type;
typedef float DATA_TYPE;
#define NUM_ITER 1
#define PTS_SIZE 10
#define OUT_SIZE_0 10
#define OUT_SIZE_1 10

void fivept_nister(point* pts1_in, point* pts2_in, my_type* out, int iter, my_type thresh);

int test_all(){

	my_type pts1[PTS_SIZE*NUM_ITER];
	my_type pts2[PTS_SIZE*NUM_ITER];
	my_type out[OUT_SIZE_0*NUM_ITER][OUT_SIZE_1];

/*
 * golden inputs
 */
	pts1[0] = 0.067;
	pts1[1] = 0.287;
	pts1[2] = 0.254;
	pts1[3] = 0.0646;
	pts1[4] = 0.239;
	pts1[5] = -0.213;
	pts1[6] = -0.710;
	pts1[7] = -0.693;
	pts1[8] = 0.661;
	pts1[9] = -0.307;

	pts2[0] = 0.329;
	pts2[1] = 1.297;
	pts2[2] = 0.523;
	pts2[3] = 1.0807;
	pts2[4] = 0.517;
	pts2[5] = 0.645;
	pts2[6] = -0.141;
	pts2[7] = 0.157;
	pts2[8] = 0.950;
	pts2[9] = 0.773;


	fivept_nister((point*)pts1, (point*)pts2, (my_type*) out, NUM_ITER, 1e-6);

	for(int i = 0; i < OUT_SIZE_0*NUM_ITER; i++){
		for(int j = 0; j < OUT_SIZE_1; j++){
			printf("%f ", out[i][j]);
		}
		printf("\n");
	}

/*
 * golden outputs
 */
	my_type golden[OUT_SIZE_0][OUT_SIZE_1] = {
			{5.0, -0.077525, -0.149317, 0.656514, 0.109315, -0.006998, -0.194794, -0.599610, 0.361774, -0.018087},
			{5.0, -0.182699, -0.360083, 0.563042, 0.225435, -0.005875, -0.096941, -0.616093, 0.272700, -0.075871},
			{5.0, -0.306013, -0.401334, 0.353560, -0.027458, -0.105090, 0.411128, -0.368289, -0.324397, -0.445632},
			{5.0, -0.216725, -0.632374, 0.149477, 0.461331, -0.027739, -0.071133, -0.476819, 0.289955, 0.026535},
			{5.0, -0.281905, -0.420364, 0.272932, -0.054526, -0.163980, 0.459042, -0.330547, -0.335147, -0.455154},
			{5.0, -0.204233, -0.634220, 0.063670, 0.508168, -0.002469, -0.110494, -0.421275, 0.311961, 0.082203},
			{5.0, -0.128642, -0.493354, -0.263300, 0.608425, 0.144048, -0.271032, -0.113786, 0.326008, 0.295142},
			{5.0, -0.218241, -0.445984, 0.050055, -0.099937, -0.276238, 0.543620, -0.207262, -0.355111, -0.447307},
			{5.0, -0.111726, -0.492821, -0.260939, -0.085134, -0.401232, 0.556796, -0.062614, -0.269288, -0.349086},
			{5.0, -0.010308, -0.494998, -0.433113, -0.068668, -0.484313, 0.503188, 0.005036, -0.132471, -0.239245}
			};

	int errors = 0;

	for(int i = 0; i < OUT_SIZE_0; i++){
		for(int j = 0; j < OUT_SIZE_1; j++){
			if(fabs(golden[i][j] - out[i][j]) > 3e-3 && golden[i][j] != out[i][j]) {
				errors++;
			}
		}
	}

	printf("Errors: %d\n", errors);

	return errors;
}



int main()
{
	int status = 0;

	status |= test_all();

	return status;
}
