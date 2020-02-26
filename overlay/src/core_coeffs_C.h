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
 * authors:
 * Marco Rabozzi (marco.rabozzi@polimi.it)
 * Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)
 * Lorenzo Di Tucci (lorenzo.ditucci@polimi.it)
 *
 * This code was auto generated using the load-store architecture compiler.
 * Original source code file: coeffC_in.h
 */
#include <stdint.h>
#include "hls_stream.h"

typedef float CFC_DATA_TYPE_IN;
typedef double CFC_DATA_TYPE_INTERNAL;
typedef float CFC_DATA_TYPE_OUT;

#define CFC_NUM_CORES 1
#define CFC_ISTR 418
#define CFC_NUM_REGS 468
#define CFC_IN_SIZE 39
#define CFC_CONSTS 0
#define CFC_OUT_SIZE 11

const uint32_t CFC_cores_instr[CFC_NUM_CORES][CFC_ISTR] = {
		{0x01100C27, 0x01A01828, 0x01207429, 0x01D0102A, 0x00207C2B, 0x0000442C, 0x00D0782D, 0x0010502E, 0x00F0782F, 0x00E07C30, 0x01C01831, 0x01E00C32, 0x00D08033, 0x01A01C34, 0x00D07C35, 0x02000C36, 0x01100837, 0x01200C38, 0x01001C39, 0x00D01C3A, 0x0010803B, 0x01D0143C, 0x00D0143D, 0x0140743E, 0x01B0143F, 0x01407040, 0x01008441, 0x01306842, 0x00004843, 0x01A01444, 0x01C01445, 0x01D01C46, 0x00E08047, 0x01106848, 0x01107049, 0x0130704A, 0x0100184B, 0x01F00C4C, 0x0010844D, 0x00E0184E, 0x01A0104F, 0x01300850, 0x00E01451, 0x00F08452, 0x01007C53, 0x01D01854, 0x00F01455, 0x01400C56, 0x00208457, 0x02100C58, 0x00E08459, 0x01B01C5A, 0x0130745B, 0x01C01C5C, 0x00007C5D, 0x0130045E, 0x0100805F, 0x01406C60, 0x01106C61, 0x01001462, 0x01200463, 0x01200864, 0x01300C65, 0x00208066, 0x00107867, 0x01306C68, 0x01107469, 0x00F07C6A, 0x00D0846B, 0x0000786C, 0x0100106D, 0x01C0106E, 0x00F0106F, 0x00F01C70, 0x01B01071, 0x00004C72, 0x00D01873, 0x00F01874, 0x00F08075, 0x00008076, 0x00E01077, 0x00E07878, 0x01406879, 0x0020787A, 0x0120707B, 0x0120687C, 0x00D0107D, 0x00107C7E, 0x01B0187F, 0x00005080, 0x01206C81, 0x01400882, 0x01007883, 0x00008484, 0x01100485, 0x00E01C86, 0x42A0D087, 0x8641B488, 0x86A1EC89, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0x4811248A, 0x4601288B, 0x8381888C, 0x85C15C8D, 0x47A1D88E, 0x8301088F, 0x8721BC90, 0x8410F891, 0x87F0EC92, 0x46B11C93, 0x4780D494, 0x4861D095, 0x8560E496, 0x82C1F497, 0x44D13098, 0x8590A499, 0x8281F89A, 0x8452109B, 0x8271549C, 0x45A0C49D, 0x47514C9E, 0x8521009F, 0x4611F0A0, 0x42E140A1, 0x46E0FCA2, 0x85E0E8A3, 0x42F0CCA4, 0x8850F4A5, 0x837144A6, 0x871174A7, 0x8540D8A8, 0x85F16CA9, 0x880138AA, 0x82D120AB, 0x482194AC, 0x846160AD, 0x84F1B0AE, 0x47012CAF, 0x8631CCB0, 0x83C198B1, 0x4320ACB2, 0x4791A0B3, 0x8431DCB4, 0x84419CB5, 0x8831A4B6, 0xC0000000, 0xC0000000, 0x00A2ACB7, 0x02425CB8, 0x0162B4B9, 0x00B2ACBA, 0x0172B8BB, 0x009244BC, 0x0152B8BD, 0x026258BE, 0x89B2C8BF, 0xC0000000, 0x4B529CC0, 0x499278C1, 0x4A3270C2, 0x49F2A4C3, 0x48D2A0C4, 0x8AC2BCC5, 0x8932CCC6, 0x0092ACC7, 0x0172B4C8, 0x00B244C9, 0x022258CA, 0x02225CCB, 0x0162B8CC, 0x0152B4CD, 0x00C244CE, 0x0192B4CF, 0x02625CD0, 0x894280D1, 0x4892D8D2, 0x88F228D3, 0x4AA220D4, 0x8A1254D5, 0x02525CD6, 0x0182B8D7, 0x024258D8, 0x008244D9, 0x00A244DA, 0x0182B4DB, 0x025258DC, 0x4902C0DD, 0x023258DE, 0x0082ACDF, 0x49221CE0, 0x00C2ACE1, 0x8A2238E2, 0x49D2C4E3, 0x4B4294E4, 0x02325CE5, 0x0192B8E6, 0x022314E7, 0x024314E8, 0x4D62E8E9, 0x015310EA, 0x018300EB, 0x00A344EC, 0x026314ED, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0x022390EE, 0x00B30CEF, 0x026390F0, 0x016310F1, 0x023314F2, 0x025390F3, 0x015300F4, 0x4D3290F5, 0x4D5230F6, 0x4C6348F7, 0x49A388F8, 0x8E3260F9, 0x8C122CFA, 0x4C2350FB, 0x4BF380FC, 0x4DD298FD, 0x023390FE, 0x025314FF, 0x4DF32D00, 0x4DE2F101, 0x01930102, 0x00C30D03, 0x00930D04, 0x01831105, 0x00B34506, 0x00830D07, 0x4CE2F908, 0x01731109, 0x0243910A, 0x0163010B, 0x0083450C, 0x00C3450D, 0x0093450E, 0x0193110F, 0x01730110, 0x00A30D11, 0x0183F112, 0x00C3DD13, 0x0253D914, 0x4DA3FD15, 0x0253F516, 0xC0000000, 0xC0000000, 0xC0000000, 0x4CC43117, 0x4EE31D18, 0x5023A919, 0x5073C11A, 0x4E74351B, 0x0193E11C, 0x00B3D51D, 0x0163E11E, 0x0263ED1F, 0x0093DD20, 0x0263F521, 0x0223F522, 0x0243ED23, 0x0223D924, 0x0243D925, 0x0263D926, 0x0083DD27, 0x0153F128, 0x0253ED29, 0x00C3E92A, 0x4BD401C9, 0x0093D52B, 0x4CF421D3, 0x5094452C, 0x4B93A12D, 0x4BB4392E, 0x4B83F92F, 0x4B742D30, 0x4E141931, 0x4F439532, 0x4D742933, 0x0173E534, 0x5103B135, 0x00A3E936, 0x0163F137, 0x0083D538, 0x00B3DD39, 0x0153E13A, 0x00A3D53B, 0x4D93C53C, 0x0083E93D, 0x0193E53E, 0x4C93B53F, 0x0233F540, 0x0093E941, 0x00C3D542, 0x0173E143, 0x4C841544, 0x0183E145, 0x4DC43D46, 0x50336D47, 0x4D83BD48, 0x4CD41149, 0x4F23294A, 0x4E63AD4B, 0x4F33414C, 0x0163E54D, 0x0173F14E, 0x00A3DD4F, 0x00B3E950, 0x0183E551, 0x0223ED52, 0x0243F553, 0x0193F154, 0x0233ED55, 0x0233D956, 0x0153E557, 0x54452158, 0x5154F959, 0x5234695A, 0x5454755B, 0xC0000000, 0x54F4655C, 0x51B5055D, 0x5273A55E, 0x51347D5F, 0xC0000000, 0xC0000000, 0xC0000000, 0x54D55960, 0x5474FD61, 0x51845D62, 0x5404A163, 0x5264A964, 0x55352D65, 0x53B49166, 0x51E4D567, 0x52B4CD68, 0x5364A569, 0x5214716A, 0x5125256B, 0x5254E56C, 0x5424F16D, 0x54A4D16E, 0x5014B16F, 0x55154170, 0x54C4DD71, 0x55550D72, 0x53D55D73, 0x5204C574, 0x5224E175, 0x52E4C176, 0x5144B577, 0x52F4E978, 0x54E45979, 0x532589CA, 0x546585D2, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0x57757D7A, 0x55A5757B, 0x57956D7C, 0x5585917D, 0x5545C17E, 0x55C5817F, 0x56857980, 0x5765E181, 0x57259982, 0x5715D183, 0x57359584, 0x56359D85, 0x56B5B586, 0x56C5B987, 0x56A5A588, 0xC0000000, 0xC0000000, 0x5595F5D1, 0xC0000000, 0xC0000000, 0xC0000000, 0x575605CB, 0x56F5E989, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0x58661D8A, 0x57C5ED8B, 0x5846098C, 0x5526018D, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0xC0000000, 0x588629CF, 0x57E625D0, 0xC0000000, 0x585635CC, 0x583631CD, 0x57F62DCE}
};

const uint32_t CFC_cores_num_outputs[CFC_NUM_CORES] = {
11
};

const uint32_t CFC_cores_output_offset[CFC_NUM_CORES] = {
0
};

void CFC_decode_instr(uint32_t instr, uint16_t &addr1, uint16_t &addr2, 
	uint16_t &addr3, uint8_t &op_code)
{
	op_code = instr >> 30;
	addr1 = (instr >> 20) & 0x3FF;
	addr2 = (instr >> 10) & 0x3FF;
	addr3 = instr & 0x3FF;
}

void CFC_core(CFC_DATA_TYPE_INTERNAL tmp1[CFC_NUM_REGS], CFC_DATA_TYPE_INTERNAL tmp2[CFC_NUM_REGS], const uint32_t instr[CFC_ISTR], const int num_outputs) 
{	
	// initialized used constants



	for(int i = 0; i < CFC_ISTR; i++) {
		#pragma HLS PIPELINE
		#pragma HLS DEPENDENCE variable=tmp1 inter false
		#pragma HLS DEPENDENCE variable=tmp2 inter false
		uint16_t addr1, addr2, addr3;
		uint8_t op;
		CFC_decode_instr(instr[i], addr1, addr2, addr3, op);

		CFC_DATA_TYPE_INTERNAL v;
		switch(op) {
			case 0:
				v = tmp1[addr1]*tmp2[addr2];
				tmp1[addr3] = v;
				tmp2[addr3] = v;
				break;
			case 1:
				v = tmp1[addr1]+tmp2[addr2];
				tmp1[addr3] = v;
				tmp2[addr3] = v;
				break;
			case 2:
				v = tmp1[addr1]-tmp2[addr2];
				tmp1[addr3] = v;
				tmp2[addr3] = v;
				break;
		}
	}
}

void CFC_write_back(CFC_DATA_TYPE_INTERNAL tmp1[CFC_NUM_CORES][CFC_NUM_REGS], hls::stream<CFC_DATA_TYPE_OUT> &out)
{
#pragma HLS INLINE off
	for(int j = 0; j < CFC_NUM_CORES; j++) {
		#pragma HLS UNROLL
		for(int i = 0; i < CFC_cores_num_outputs[j]; i++) {
			#pragma HLS UNROLL
			out.write(tmp1[j][CFC_ISTR + CFC_IN_SIZE + CFC_CONSTS + i]);
		}
	}
}

void dataRead(double tmp1[CFC_NUM_CORES][CFC_NUM_REGS],
		double tmp2[CFC_NUM_CORES][CFC_NUM_REGS],
		hls::stream<CFC_DATA_TYPE_IN>& in) {
#pragma HLS INLINE off
	for (int i = 0; i < CFC_IN_SIZE; i++) {
#pragma HLS PIPELINE
		double in_val = in.read();
		for (int j = 0; j < CFC_NUM_CORES; j++) {
			tmp1[j][i] = in_val;
			tmp2[j][i] = in_val;
		}
	}
}

void CFC_multiCore(hls::stream<CFC_DATA_TYPE_IN> &in, hls::stream<CFC_DATA_TYPE_OUT> &out)
{
	CFC_DATA_TYPE_INTERNAL tmp1[CFC_NUM_CORES][CFC_NUM_REGS];
	CFC_DATA_TYPE_INTERNAL tmp2[CFC_NUM_CORES][CFC_NUM_REGS];

	#pragma HLS ARRAY_PARTITION variable=CFC_cores_instr complete dim=1
	#pragma HLS ARRAY_PARTITION variable=CFC_cores_num_outputs complete dim=1
	#pragma HLS ARRAY_PARTITION variable=tmp1 complete dim=1
	#pragma HLS ARRAY_PARTITION variable=tmp2 complete dim=1

	#pragma HLS RESOURCE variable=tmp1 core=RAM_S2P_BRAM
	#pragma HLS RESOURCE variable=tmp2 core=RAM_S2P_BRAM

	dataRead(tmp1, tmp2, in);

	for(int j = 0; j < CFC_NUM_CORES; j++) {
		#pragma HLS UNROLL
		CFC_core(tmp1[j], tmp2[j], CFC_cores_instr[j], CFC_cores_num_outputs[j]);
	}

	CFC_write_back(tmp1, out);
}
