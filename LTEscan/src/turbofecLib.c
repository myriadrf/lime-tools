/*
 Copyright 2017 Lime Microsystems Ltd.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
// provide an easy to use library for handling viterbi and turbo code data
// LTE max segment length 6144 (6120bit + 24bit CRC)

// note type issues, lte_rate_matcher expects dd and e to be 'signed char'
// note lte_turbo_decode_unpack expects dd and e to be int8_t, and d to be uint8_t

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "turbofecLib.h"
#ifdef __cplusplus // If this is a C++ compiler, use C linkage
extern "C" {
#endif
#include "turbofec/conv.h"
//#include "turbofec/turbo.h"
#include "turbofec/rate_match.h"
#ifdef __cplusplus // If this is a C++ compiler, end C linkage
}
#endif

#ifndef _CONV_H_
	#warning No _CONV_H_, check if turbofec installed
#endif // check if #include "turbofec/conv.h" found ok

#include <errno.h>

void conv_bin2bpsk( signed char e[], unsigned int eLen );
void conv_zero_dec( unsigned char in[],unsigned int dLen );
void rate_match_demo( void );

static signed char *dd[3]={NULL,NULL,NULL};
static int8_t *bsv=NULL;
static uint8_t *bu0=NULL;
static uint8_t *bu1=NULL;
static uint8_t *bu2=NULL;
static uint8_t *buv=NULL;

#define MAX_LEN_BITS 512

static struct lte_conv_code lte_conv_pbch = 
{
	.n = 3, // 2:4
	.k = 7, // 5 or 7
	.len = MAX_LEN_BITS, // was 512
	.rgen = 0,
	.gen = { 0133, 0171, 0165 },
	.punc = NULL, // NULL=>unpunctured
	.term = CONV_TERM_TAIL_BITING, // recursive not supported
}; // was const

static struct lte_rate_matcher pbch_rm =
{
	.E = 0, // automatically generated if E and D are zero
	.D = 0, // automatically generated if E and D are zero
	.V = 0, // automatically generated if E and D are zero
	.rows = 0, // automatically generated if E and D are zero
	.w_null = NULL, // int *, must be set to NULL prior to use, free/calloc
	.w_null_cnt = 0, // int
	.rv = 0, // int - revision number
	.z = {NULL,NULL,NULL}, // signed char [3]
	.v = {NULL,NULL,NULL}, // signed char [3]
	.pi = NULL, // * int
};

static struct lte_rate_matcher_io pbch_rm_io =
{
	.D = MAX_LEN_BITS, // int
	.E = 3*MAX_LEN_BITS, // int was 2*
	.d={NULL,NULL,NULL}, // signed char* [3] needs to be defined before use, need to reformat output of viterbi prior to calling!
	.e = NULL, // signed char* needs to be defined before use
};

void ViterbiInit( void )
{
	unsigned char ci=0;
	buv = (unsigned char *) malloc(sizeof(unsigned char) * 3 * MAX_LEN_BITS);  // viterbi dec [d0,d1,d2]
	bsv = (signed char *) malloc(sizeof(signed char) * 3 * MAX_LEN_BITS);  // viterbi dec [d0,d1,d2]
	for( ci=0; ci<3; ci++ )
	{
		dd[ci]=(signed char *)malloc(sizeof(signed char)*(MAX_LEN_BITS+4)); // viterbi/turbo d0,d1,d2
		pbch_rm_io.d[ci]=dd[ci];
	}
}

void ViterbiFree( void )
{
	unsigned char ci;
	if(buv!=NULL)
		free(buv);
	if(bsv!=NULL)
		free(bsv);
	for( ci=0; ci<3; ci++ )
		if(dd[ci]!=NULL)
			free(dd[ci]);
}

void ViterbiRateEnc( unsigned char d[],unsigned int dLen,signed char e[],unsigned int eLen,unsigned char rv )
{ // note in turbofec code, convention if( ptr ) <==> if( ptr!=NULL )
	int l;
	unsigned int cj=0;
	unsigned char ci=0;
	if(dLen>MAX_LEN_BITS)
	{
		printf("ViterbiRateEnc:ERROR size %i exceeds storage %i\n",dLen,MAX_LEN_BITS);
		return;
	}
	lte_conv_pbch.len=dLen;
	pbch_rm_io.e=e;
	pbch_rm_io.D=dLen;
	pbch_rm_io.E=eLen;
	pbch_rm.rv=rv;
	pbch_rm.D=0; // force memory update, rows and V autocalculated.
	pbch_rm.E=0;
	l = lte_conv_encode(&lte_conv_pbch, (uint8_t *)d, (uint8_t *)buv); 
	if (l != (lte_conv_pbch.n*dLen)) // returns code->n*l
	{
		printf("ViterbiRateEnc: ERROR failed encoding length check! %i\n",l);
		return;
	}
	for( cj=0; cj<dLen; cj++ ) // (unsigned) buv= (signed)[d0_0,d1_0,d2_0,d0_1,...]	
		for( ci=0; ci<3; ci++ )
			dd[ci][cj]=buv[cj*3+ci];
	lte_conv_rate_match_fw(&pbch_rm, &pbch_rm_io );
}

void ViterbiRateDec( signed char e[],unsigned int eLen,unsigned char d[],unsigned int dLen,unsigned char rv )
{
	int err=0;
	unsigned int cj=0;
	unsigned char ci=0;
	if(dLen>MAX_LEN_BITS)
	{
		printf("ViterbiRateDec:ERROR size %i exceeds storage %i\n",dLen,MAX_LEN_BITS);
		return;
	}
	lte_conv_pbch.len=dLen;
	pbch_rm_io.e=e;
	pbch_rm_io.D=dLen;
	pbch_rm_io.E=eLen;
	pbch_rm.rv=(int)rv;
	pbch_rm.D=0; // force memory update, rows and V autocalculated.
	pbch_rm.E=0;
	err=lte_conv_rate_match_rv(&pbch_rm, &pbch_rm_io);
	if( err ) // err != 0
		printf("ViterbiRateDec:RateMatch err=%i\n",err);
	for( cj=0; cj<dLen; cj++ )
		for( ci=0; ci<3; ci++ )
			bsv[cj*3+ci]=dd[ci][cj];  
	lte_conv_decode(&lte_conv_pbch, (int8_t *)bsv, (uint8_t *)d);
}

void conv_bin2bpsk( signed char e[], unsigned int eLen )
{
	int ci=0;
	for( ci=0; ci<eLen; ++ci )  // vectorizes
	{
		e[ci]*=124;
		e[ci]+=-62;
	}
}

void conv_zero_dec( unsigned char in[],unsigned int dLen )
{
	unsigned char ci=0;
	unsigned int cj=0;
	for( cj=0; cj<dLen; cj++ )
	{
		for( ci=0; ci<3; ci++ )
			dd[ci][cj]=7;
		in[cj]=7;
	}
}


