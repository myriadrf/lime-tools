/*
 Copyright 2019 Lime Microsystems Ltd.

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "lime.h"
#ifdef USE_C99_COMPLEX
#include <complex.h> // double _Complex, compliant with double, fftw_complex[2] Re=0, Im=1
#include <fftw3.h>  // should come after complex.h
#else // #ifdef USE_FFTW_COMPLEX
#include <fftw3.h>  // Do before complex.h to force typedef double fftw_complex[2]
#include <complex.h> 
#endif

#include "CellItem.hpp"
#include "CellList.hpp"
#include "LTEscan.hpp"

void LTEscan::LoadDictPSCH( void )
{ // note differences in TDD and FDD on location of PSCH
	unsigned char u[3]={25,29,34}; // N2ID; TS 136.211 13.0.0 Table 6.11.1.1-1
	signed short cn=0;
	unsigned short cj;
	unsigned char c2id=0;
	float _Complex temp;
	//printf("LTEscan::LoadDictPSCH\n");
	for( c2id=0; c2id<3; c2id++ )
	{
		for(cj=0;cj<NFFT;cj++)
			dictPSCHzdcNFFT[c2id][cj]=0.0;
		for(cj=0;cj<64;cj++)
			dictPSCHzdc[c2id][cj]=0.0;
		temp=-I*PI*u[c2id]/63.0;
		for( cn=0;cn<31;cn++ ) // 0:30 (Negative frequencies)
			dictPSCHzdc[c2id][cn]=cexpf(temp*cn*(cn+1));
		for( cn=31;cn<62;cn++ ) // 32:62 (Positive frequenices) 31 DC
			dictPSCHzdc[c2id][cn+1]=cexpf(temp*(cn+1)*(cn+2)); // [c2id][cn+1], implies DC=0
		for( cn=0;cn<31;cn++ ) // 0:30 (Negative frequencies)
			dictPSCHzdcNFFT[c2id][NFFT+cn-31]=dictPSCHzdc[c2id][cn];
		for( cn=31;cn<63;cn++ ) // 31:61 (Positive frequenices), note [0] is DC.
			dictPSCHzdcNFFT[c2id][cn-31]=dictPSCHzdc[c2id][cn];
#ifdef USE_C99_COMPLEX
		for( cn=0;cn<NFFT;cn++ ) // vectorized
			in[cn]=dictPSCHzdc[c2id][cn];
#else // #ifdef USE_FFTW_COMPLEX
		for(cn=0;cn<NFFT;cn++) // copy SDR buffer to FFTW and window (Vectorized)
		{
			in[cn][0]=crealf(dictPSCHzdc[c2id][cn]);
			in[cn][1]=cimagf(dictPSCHzdc[c2id][cn]);
		}
#endif
		fftw_execute(pfft); // note this is an fft of a pattern for ACF.
#ifdef USE_C99_COMPLEX
		for( cn=0;cn<NFFT;cn++ ) // vectorized
			dictPSCHconjFFT[c2id][cn]=conj(out[cn]);
#else // #ifdef USE_FFTW_COMPLEX
		for( cn=0;cn<NFFT;cn++ ) // vectorized	
			dictPSCHconjFFT[c2id][cn]=out[cn][0]-I*out[cn][1]; // conj(x)=re(x)-im(x)
#endif
		for(cj=0;cj<64;cj++)
			dictPSCHzdcPh[c2id][cj]=cargf(dictPSCHzdc[c2id][cj]);
	}	
}

void LTEscan::LoadDictSSCH(void) // TS 136.211 Section 6.11.2.2
{ // mostly depends on N1ID, but note co/c1 depend on N2ID
	signed char states[5];
	signed char statec[5];
	signed char statez[5];

	signed char stilda[31];
	signed char ctilda[31];
	signed char ztilda[31];

	signed char s0[31];
	signed char c0[31];
	signed char z0[31];
	signed char s1[31];
	signed char c1[31];
	signed char z1[31];
	unsigned char c2id=0;
	short c1id=0;
	short idx;
	unsigned char c=0;
//	printf("LTEDL::CalcAllSSCH\n");
	for( c1id=0; c1id<168; c1id++)
	{
		signed short qdash=floor(c1id/30);
		signed short q=floor((c1id+qdash*(qdash+1)/2)/30);
		signed short mdash=c1id+q*(q+1)/2;
		signed short m0=mdash % 31; // mod
		signed short m1=(m0+floor(mdash/31)+1); // Also in TS 136.211 13.0.0 Table 6.11.2.1-1
		m1=m1 % 31; // mod
		for( c=0; c<5; c++ ) // initial state 0 ---> +1
		{
			states[c]=1;
			statec[c]=1; 
			statez[c]=1; 
		} // inital state 1 ---> -1
		states[4]=-1;
		statec[4]=-1;
		statez[4]=-1;
			
		for( c=0; c<31; c++ )
		{
			stilda[c]=states[0];
			ztilda[c]=statez[0];
			ctilda[c]=statec[0];
			LTEDLGold1SSSR( states );
			LTEDLGold2SSSR( statez );
			LTEDLGold3SSSR( statec );
		}
		for( c2id=0; c2id<3; c2id++ )
		{
			idx=c2id+3*c1id;
			Copy31(stilda,s0);
			Copy31(stilda,s1);
			Copy31(ctilda,c0);
			Copy31(ctilda,c1);
			Copy31(ztilda,z0);
			Copy31(ztilda,z1);
			CycRot31(s0, -m0);
			CycRot31(s1, -m1);
			CycRot31(c0, -c2id);
			CycRot31(c1, -(c2id+3));
			CycRot31(z0, -(m0 % 8)); // mod
			CycRot31(z1, -(m1 % 8)); // mod
			for( c=0; c<31; c++ )
			{
				dictSSCH0zdc[idx][2*c]=s0[c]*c0[c]; // nSlot==0 (SF=0-4 TDD?)
				dictSSCH0zdc[idx][2*c+1]=s1[c]*c1[c]*z0[c];
				dictSSCH10zdc[idx][2*c]=s1[c]*c0[c]; // nSlot==10 (SF=5-9 TDD?)
				dictSSCH10zdc[idx][2*c+1]=s0[c]*c1[c]*z1[c];
			}
			for( c=61; c>30; c-- ) // make room for DC
				dictSSCH0zdc[idx][c+1]=dictSSCH0zdc[idx][c];
			for( c=61; c>30; c-- )
				dictSSCH10zdc[idx][c+1]=dictSSCH10zdc[idx][c];
			dictSSCH0zdc[idx][31]=0.0;
			dictSSCH10zdc[idx][31]=0.0;
			dictSSCH0zdc[idx][63]=0.0;
			dictSSCH10zdc[idx][63]=0.0;
			for(c=0;c<64;c++)
				dictSSCH0zdcPh[idx][c]=cargf(dictSSCH0zdc[idx][c]);
			for(c=0;c<64;c++)
				dictSSCH10zdcPh[idx][c]=cargf(dictSSCH10zdc[idx][c]);
		}
	}
}

void LTEscan::LTEDLGold1SSSR( signed char *state )
{
	signed char fb=SXOR(state[0],state[2]);
	for( unsigned char i=0; i<4; ++i )
		state[i]=state[i+1];
	state[4]=fb;
}

void LTEscan::LTEDLGold2SSSR( signed char *state )
{
	signed char fb=SXOR(SXOR(SXOR(state[0],state[1]),state[2]),state[4]);
	for( unsigned char i=0; i<4; ++i )
		state[i]=state[i+1];
	state[4]=fb;
}

void LTEscan::LTEDLGold3SSSR( signed char *state )
{
	signed char fb=SXOR(state[0],state[3]);
	for( unsigned char i=0; i<4; ++i )
		state[i]=state[i+1];
	state[4]=fb;
}

void LTEscan::Copy31( signed char *src, signed char *dest )
{
	for( unsigned char i=0; i<31; ++i ) // vectorizes
		dest[i]=src[i];
}

void LTEscan::CycRot31( signed char *seq, signed short rot ) // same as shift(deq,rot) in octave
{ // rot increased to signed short from signed char as |rot|>127
	signed char e=0;
	if( rot <0 )
		for( signed short c=0; c<-rot; ++c )
		{
			e=seq[0];
			for( unsigned char i=0; i<30; ++i )
				seq[i]=seq[i+1];
			seq[30]=e;
		}

	if( rot >0 )
		for( signed short c=0; c<rot; ++c )
		{
			e=seq[30];
			for( unsigned char i=29; i>=0; --i )
				seq[i+1]=seq[i];
			seq[0]=e;
		}
	// if rot==0 do nothing!
}

signed char LTEscan::SXOR( signed char a, signed char b ) // we are doing signed XOR with inverted BPSK 0=1, 1=-1
{ // XOR a==b => 0 ---> 1 XOR a!=b => 1 ---> -1
	signed char c=-1; // a != b
	if( a==b )
		c=1; // 1 ---> -1
//	else
//		c=-1;
	return(c);
}

void LTEscan::LTESC( unsigned char *psc,unsigned int cinit,unsigned int cnt )
{ // TS 136.211 Section 7.2 gold code calculation, x10 fast as original
	uint32_t state[4]={1,cinit,0,0}; // use vector form to help compiler vectorize
	//uint8_t fb[2];
	unsigned short cj;
	unsigned int ck;
	for(cj=1;cj<1600;++cj) // use parallel form to reduce loop overhead and help vectorize
	{ // does not vectorize on SIMD 128?
		state[2]=(__builtin_popcountll(state[0]&9))&1; // mask=1001b, popcount=hamming distance= bitcount
		state[3]=(__builtin_popcountll(state[1]&15))&1; // mask=1111b, popcount=hamming distance= bitcount
		state[0]>>=1;
		state[1]>>=1;
		state[0]|=(1<<30)*state[2]; // multiply by 2^30
		state[1]|=(1<<30)*state[3]; // multiply by 2^30
	}
	for(ck=0;ck<cnt;++ck)
	{ // does not vectorize on SIMD 128?
		state[2]=(__builtin_popcountll(state[0]&9))&1; // mask=1001b, popcount=hamming distance= bitcount
		state[3]=(__builtin_popcountll(state[1]&15))&1; // mask=1111b, popcount=hamming distance= bitcount
		state[0]>>=1;
		state[1]>>=1;
		state[0]|=(1<<30)*state[2]; // multiply by 2^30
		state[1]|=(1<<30)*state[3]; // multiply by 2^30
		psc[ck]=(state[0]&1)^(state[1]&1);
	}
}

// ref prs is the SAME for ALL MIMO antennas with SAME l, regardless of BW
// but pilot location mapped to UNIQUE subcarrier locations for EACH antenna
// the purpose of the refs (pilots) is to equalise other OFDM symbols via interpolation of amp and phase
// in theory we could use multiple slots to improve SNR, but this involves multiple FFTs, ignore for now.
// PBCH 6 RB wide, each with 2 refs, 12 refs total.
// PBCH in middle, so {0,9,19,44,69,94} is the displacement from beginning of PRS sequence, for a given NDLRB.
// m=0..(2*NDLRB-1)
// PRS(2*m')+PRS(2*m'+1) m'=m+110-NDLRB
// PRS     0..439  = 4*110-1 (2 subcarriers, R+I for each point).  Actually used is 20..419
// Refs    0..219
// 20MHz  10..209, PBCH 94..105 --> 104..115 --> PRS 
// 15MHz  35..184, PBCH 69..80  --> 104..115
// 10MHz  60..159, PBCH 44..55  --> 104..115
//  5MHz  85..134, PBCH 19..31  --> 104..115
//  3MHz  95..124, PBCH  9..21  --> 104..115
//  1MHz 104..115, PBCH  0..11  --> 104..115
// So PBCH PRS is independent of NDLRB
void LTEscan::LoadRefs(unsigned char nSlot,unsigned char nSymbol,unsigned char NID,unsigned char NCP)
{ // TS 136.211 Section 6.10.1.2.  Implement mapping of QPSK pilot sequence.  See see fig 6.10.1.2-1 and -2
	printf("LTEDL::LoadRefs ns%i l%i NID%i NCP%i\n",nSlot,nSymbol,NID,NCP);
	unsigned char cprs[440];
	unsigned char nspri=nSlot;
	unsigned int cinit=1024*(7*(nspri+1)+nSymbol+1)*(2*NID+1)+2*NID+NCP; // 6.10.1.1
	unsigned short ci;
	LTESC(cprs,cinit,440); // init at the start of each symbol
	for(ci=0;ci<12;ci++ )
	{
		pilotVal[ci]=1+I;
		pilotVal[ci]-=2*cprs[2*ci+208]+(I*2)*cprs[2*ci+209]; // 2*(110-6)=208
		pilotVal[ci]*=rsqrt2;
	}
}


