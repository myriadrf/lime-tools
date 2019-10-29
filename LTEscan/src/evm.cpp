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
#include <math.h>
#include <time.h> // tools for telling the time!
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

#ifdef MXD_CMPLR // needed if C compiled with gcc/clang, and cpp compiled g++/clang++
extern "C" {
#endif
void FftShift( float _Complex *buffer );
#ifdef MXD_CMPLR 
}
#endif

void LTEscan::CalcRMS62( short foff,float *p ) // vecm is c * conj(c), not 
{ // Sync blocks consist of 2x3 RBs (2x36 s/c), separated by DC, offset is assumed to lie on DC s/c
	unsigned char cf=0;
	float *ptrm=NULL;
	if( (foff<36) || (NFFTm37<foff) ) // [0] is invalid
		printf("Range Fail CalcRMS62 foff=%i\n",foff);
	else
	{
		p[0]=0.0f;
		ptrm=ffttmpMg+foff-36; // move from DC, to edge of negative 36 s/c, first 5 s/c are not used
		for(cf=5;cf<36;cf++) // vectorized
			p[0]+=ptrm[cf];
		for(cf=37;cf<68;cf++) // vectorized
			p[0]+=ptrm[cf];
		p[0]*=r62;
	}
} // this is independent of c2id

void LTEscan::CalcRMS10( short foff,float *p )
{ // Sync blocks consist of 2x3 RBs (2x36 s/c), separated by DC, offset is assumed to lie on DC s/c
	unsigned char cf=0;
	float *ptrm=NULL;
	if( (foff<36) || (NFFTm37<foff) ) // [0] is invalid
		printf("Range Fail CalcRMS10 foff=%i\n",foff);
	else
	{
		p[1]=0.0f;
		ptrm=ffttmpMg+foff-36; // move from DC, to edge of negative 36 s/c, first 5 s/c are not used
 		for( cf=0;cf<5;cf++) // vectorized
			p[1]+=ptrm[cf]; // these unused s/c contain noise, so can use to make SNR measurement
		for( cf=68;cf<73;cf++)  // vectorized
			p[1]+=ptrm[cf]; // these unused s/c contain noise, so can use to make SNR measurement
		p[1]*=0.1f; // r10
	}
} // this is independent of c2id

void LTEscan::MakeEqvVecFromXSCH1(short cid,short foff,float _Complex **dict,float **dictPh,float *p,float *sd)
{ // FFT X(w)exp(jwTD/NFFT) <--> x(t+TD) ==> we can estimate delay p[4]+p[3]*cf linear phase
	float phd[63];
	float _Complex eqv[64];
	float _Complex evm[64];
	signed char cf=0;
	float *ptrPh=ffttmpPh+foff-31; //'zero DC' version of dictionary to generate 'zero DC' equalisation vector
	float _Complex *ptr;
	float _Complex wph[2]={0.0f}; // 0 wphc, 1 wphm
	if( (foff<36) || (NFFTm37<foff) )
	{
		printf("Range Fail MakeEqvVecFromXSCH foff=%i\n",foff);
		return;
	}

	phd[30]=ptrPh[30]-dictPh[cid][30]; // negative freq error
	phd[32]=ptrPh[32]-dictPh[cid][32]; // positive freq error
	p[3]=(phd[30]-phd[32])*0.5f; // time delay is linear phase shift
	p[4]=(phd[30]+phd[32])*0.5f; // calculate constant phase offset - phErr - averaging usually more accurate.
	sd[0]=(phd[30]-p[4])*(phd[30]-p[4]); // SIMD friendly
	sd[1]=(phd[30]-p[4])*(phd[30]-p[4]);
	sd[2]=(phd[30]-p[4])*(phd[30]-p[4]);
	sd[0]+=(phd[32]-p[4])*(phd[32]-p[4]);
	sd[1]+=(phd[32]-p[4])*(phd[32]-p[4]);
	sd[2]+=(phd[32]-p[4])*(phd[32]-p[4]);
	sd[0]-=float(2*PI2)*(phd[32]-p[4]);
	sd[2]+=float(2*PI2)*(phd[32]-p[4]);
	sd[0]+=float(PI2*PI2); // was +=(phd[32]-p[4]-PI2)*(phd[32]-p[4]-PI2);
	sd[2]+=float(PI2*PI2); // was +==(phd[32]-p[4]+PI2)*(phd[32]-p[4]+PI2)
	sd[0]*=0.5f;
	sd[1]*=0.5f;
	sd[2]*=0.5f;
	if( (sd[1]>sd[0]) && (sd[2]>sd[0]) )
		p[4]-=float(PI); // was =(phd[30]+phd[32]-PI2)*0.5;
	if( (sd[1]>sd[2]) && (sd[2]>sd[0]) )
		p[4]+=float(PI); // was =(phd[30]+phd[32]+PI2)*0.5;
	//printf("foff=%i phm=%.3f phc=%.3f\n",foff,phm,phc);
	//for( unsigned char ci=28;ci<34;ci++)
	//	printf("ci=%i |%.3f| <%.3f |%.3f| <%.3f %.3f %.3f\n",ci,cabsf(ptr[ci]),cargf(ptr[ci]),cabsf(dict[cid][ci]),cargf(dict[cid][ci]),t[0],t[1]);
		
	// expect delay to be within +/-0.5 rad (+/-(10*2*PI)/128) if within CP, and phase err=0
	p[3]+=float(PI)*((p[3]<float(-PI))-(float(PI)<p[3]));
	p[4]+=float(PI)*((p[4]<float(-PI))-(float(PI)<p[4]));
	wph[0]=cexpf(-I*p[4] ); // FAST WAY - reduces calc time from 65ms to 16ms
	wph[1]=cexpf(-I*p[3] );
	wph[0]*=rrms;
	eqv[31]=wph[0]; // DC has zero linear phase shift
	eqv[63]=1.0f; // not used
	//printf("phm=%.3f phc=%.3f rrms=%.3f t[0]=%.3f t[1]=%.3f td=%.3f fdop=%.4f\n",phm,phc,rrms,t[0],t[1],td,fdop);
	eqv[30]=wph[0]*wph[1]; // work from smallest phase shift to minimise cumulative error
	for(cf=29;cf>=0;cf--) // DOES NOT VECTORIZE - history dependent
		eqv[cf]=eqv[cf+1]*wph[1]; // calculate negative frequency linear phase shift
//	wphm=cexpf(I*p[3] );
	wph[1]=1.0f/wph[1]; // comparable to wphm=cexpf(I*phm ); ?
	eqv[32]=wph[0]*wph[1]; // work from smallest phase shift to minimise cumulative error
	for(cf=33;cf<63;cf++) // DOES NOT VECTORIZE - history dependent
		eqv[cf]=eqv[cf-1]*wph[1]; // calculate negative frequency linear phase shift
	p[3]*=-NFFTdPI2; // td - real, not integer.

	p[2]=0.0f; // evm^2 total
	ptr=ffttmp+foff-31; // start at active part of PSCH/SSCH,
	for(cf=0;cf<63;cf++) // vectorized
		evm[cf]=ptr[cf]*eqv[cf]; // copy part of fft data
	for(cf=0;cf<63;cf++) // vectorized
		evm[cf]-=dict[cid][cf]; // calc EVM (zero DC ver of Dict)
	evm[31]=0.0;
	for(cf=0;cf<63;cf++) // vectorized
		evm[cf]*=conjf(evm[cf]); // convert to ||^2
	for(cf=0;cf<63;cf++) // vectorized
		p[2]+=crealf(evm[cf]);
	p[2]*=r62;
} // this is dependent on c2id code

void LTEscan::MakeEqvVecFromXSCH1A(short cid,short foff,float _Complex **dict,float **dictPh,float *p,float *sd)
{ // FFT X(w)exp(jwTD/NFFT) <--> x(t+TD) ==> we can estimate delay p[4]+p[3]*cf linear phase
	float phd[63];
	float _Complex eqv[64];
	float _Complex evm[64];
	signed char cf=0;
	float *ptrPh=ffttmpPh+foff-31; //'zero DC' version of dictionary to generate 'zero DC' equalisation vector
	float _Complex *ptr;
	float _Complex wph[2]={0.0f}; // 0 wphc, 1 wphm
	if( (foff<36) || (NFFTm37<foff) )
	{
		printf("Range Fail MakeEqvVecFromXSCH foff=%i\n",foff);
		return;
	}
	phd[27]=ptrPh[27]-dictPh[cid][27]; // negative freq error [27,28,29,30], use +/-4th s/c instead of 1st s/c
	phd[35]=ptrPh[35]-dictPh[cid][35]; // positive freq error [32,33,34,35], use +/-4th s/c instead of 1st s/c
	p[3]=(phd[27]-phd[35])*0.125f; // time delay is linear phase shift
	p[4]=(phd[27]+phd[35])*0.5f; // calculate constant phase offset - phErr - averaging usually more accurate.
	sd[0]=(phd[27]-p[4])*(phd[27]-p[4]); // SIMD friendly
	sd[1]=(phd[27]-p[4])*(phd[27]-p[4]);
	sd[2]=(phd[27]-p[4])*(phd[27]-p[4]);
	sd[0]+=(phd[35]-p[4])*(phd[35]-p[4]);
	sd[1]+=(phd[35]-p[4])*(phd[35]-p[4]);
	sd[2]+=(phd[35]-p[4])*(phd[35]-p[4]);
	sd[0]-=float(2*PI2)*(phd[35]-p[4]);
	sd[2]+=float(2*PI2)*(phd[35]-p[4]);
	sd[0]+=float(PI2*PI2); // was +=(phd[32]-p[4]-PI2)*(phd[32]-p[4]-PI2);
	sd[2]+=float(PI2*PI2); // was +==(phd[32]-p[4]+PI2)*(phd[32]-p[4]+PI2)
	sd[0]*=0.5f;
	sd[1]*=0.5f;
	sd[2]*=0.5f;
	if( (sd[1]>sd[0]) && (sd[2]>sd[0]) )
		p[4]-=float(PI); // was =(phd[30]+phd[32]-PI2)*0.5;
	if( (sd[1]>sd[2]) && (sd[2]>sd[0]) )
		p[4]+=float(PI); // was =(phd[30]+phd[32]+PI2)*0.5;
	//printf("foff=%i phm=%.3f phc=%.3f\n",foff,phm,phc);
	//for( unsigned char ci=28;ci<34;ci++)
	//	printf("ci=%i |%.3f| <%.3f |%.3f| <%.3f %.3f %.3f\n",ci,cabsf(ptr[ci]),cargf(ptr[ci]),cabsf(dict[cid][ci]),cargf(dict[cid][ci]),t[0],t[1]);
		
	// expect delay to be within +/-0.5 rad (+/-(10*2*PI)/128) if within CP, and phase err=0
	p[3]+=float(PI)*((p[3]<float(-PI))-(float(PI)<p[3]));
	p[4]+=float(PI)*((p[4]<float(-PI))-(float(PI)<p[4]));
	wph[0]=cexpf(-I*p[4] ); // FAST WAY - reduces calc time from 65ms to 16ms
	wph[1]=cexpf(-I*p[3] );
	wph[0]*=rrms;
	eqv[31]=wph[0]; // DC has zero linear phase shift
	eqv[63]=1.0f; // not used
	//printf("phm=%.3f phc=%.3f rrms=%.3f t[0]=%.3f t[1]=%.3f td=%.3f fdop=%.4f\n",phm,phc,rrms,t[0],t[1],td,fdop);
	eqv[30]=wph[0]*wph[1]; // work from smallest phase shift to minimise cumulative error
	for(cf=29;cf>=0;cf--) // DOES NOT VECTORIZE - history dependent
		eqv[cf]=eqv[cf+1]*wph[1]; // calculate negative frequency linear phase shift
//	wphm=cexpf(I*p[3] );
	wph[1]=1.0f/wph[1]; // comparable to wphm=cexpf(I*phm ); ?
	eqv[32]=wph[0]*wph[1]; // work from smallest phase shift to minimise cumulative error
	for(cf=33;cf<63;cf++) // DOES NOT VECTORIZE - history dependent
		eqv[cf]=eqv[cf-1]*wph[1]; // calculate negative frequency linear phase shift
	p[3]*=-NFFTdPI2; // td - real, not integer.

	p[2]=0.0f; // evm^2 total
	ptr=ffttmp+foff-31; // start at active part of PSCH/SSCH,
	for(cf=0;cf<63;cf++) // vectorized
		evm[cf]=ptr[cf]*eqv[cf]; // copy part of fft data
	for(cf=0;cf<63;cf++) // vectorized
		evm[cf]-=dict[cid][cf]; // calc EVM (zero DC ver of Dict)
	evm[31]=0.0;
	for(cf=0;cf<63;cf++) // vectorized
		evm[cf]*=conjf(evm[cf]); // convert to ||^2
	for(cf=0;cf<63;cf++) // vectorized
		p[2]+=crealf(evm[cf]);
	p[2]*=r62;
} // this is dependent on c2id code

void LTEscan::MakeEqvVecFromXSCH1B(short cid,short foff,float _Complex **dict,float **dictPh,float *p,float *sd)
{ // FFT X(w)exp(jwTD/NFFT) <--> x(t+TD) ==> we can estimate delay p[4]+p[3]*cf linear phase
	float phd[63];
	float _Complex eqv[64];
	float _Complex evm[64];
	signed char cf=0;
	float *ptrPh=ffttmpPh+foff-31; //'zero DC' version of dictionary to generate 'zero DC' equalisation vector
	float _Complex *ptr;
	float _Complex wph[2]={0.0f}; // 0 wphc, 1 wphm
	if( (foff<36) || (NFFTm37<foff) )
	{
		printf("Range Fail MakeEqvVecFromXSCH foff=%i\n",foff);
		return;
	} // more low amp detects than +/-6, but less than +/-4, EVM generally better.
	phd[26]=ptrPh[26]-dictPh[cid][26]; // negative freq error [27,28,29,30], use +/-5th s/c instead of 1st s/c
	phd[36]=ptrPh[36]-dictPh[cid][36]; // positive freq error [32,33,34,35], use +/-5th s/c instead of 1st s/c
	p[3]=(phd[26]-phd[36])*0.1f; // time delay is linear phase shift
	p[4]=(phd[26]+phd[36])*0.5f; // calculate constant phase offset - phErr - averaging usually more accurate.
	sd[0]=(phd[26]-p[4])*(phd[26]-p[4]); // SIMD friendly
	sd[1]=(phd[26]-p[4])*(phd[26]-p[4]);
	sd[2]=(phd[26]-p[4])*(phd[26]-p[4]);
	sd[0]+=(phd[36]-p[4])*(phd[36]-p[4]);
	sd[1]+=(phd[36]-p[4])*(phd[36]-p[4]);
	sd[2]+=(phd[36]-p[4])*(phd[36]-p[4]);
	sd[0]-=float(2*PI2)*(phd[36]-p[4]);
	sd[2]+=float(2*PI2)*(phd[36]-p[4]);
	sd[0]+=float(PI2*PI2); // was +=(phd[32]-p[4]-PI2)*(phd[32]-p[4]-PI2);
	sd[2]+=float(PI2*PI2); // was +==(phd[32]-p[4]+PI2)*(phd[32]-p[4]+PI2)
	sd[0]*=0.5f;
	sd[1]*=0.5f;
	sd[2]*=0.5f;
	if( (sd[1]>sd[0]) && (sd[2]>sd[0]) )
		p[4]-=float(PI); // was =(phd[30]+phd[32]-PI2)*0.5;
	if( (sd[1]>sd[2]) && (sd[2]>sd[0]) )
		p[4]+=float(PI); // was =(phd[30]+phd[32]+PI2)*0.5;
	//printf("foff=%i phm=%.3f phc=%.3f\n",foff,phm,phc);
	//for( unsigned char ci=28;ci<34;ci++)
	//	printf("ci=%i |%.3f| <%.3f |%.3f| <%.3f %.3f %.3f\n",ci,cabsf(ptr[ci]),cargf(ptr[ci]),cabsf(dict[cid][ci]),cargf(dict[cid][ci]),t[0],t[1]);
		
	// expect delay to be within +/-0.5 rad (+/-(10*2*PI)/128) if within CP, and phase err=0
	p[3]+=float(PI)*((p[3]<float(-PI))-(float(PI)<p[3]));
	p[4]+=float(PI)*((p[4]<float(-PI))-(float(PI)<p[4]));
	wph[0]=cexpf(-I*p[4] ); // FAST WAY - reduces calc time from 65ms to 16ms
	wph[1]=cexpf(-I*p[3] );
	wph[0]*=rrms;
	eqv[31]=wph[0]; // DC has zero linear phase shift
	eqv[63]=1.0f; // not used
	//printf("phm=%.3f phc=%.3f rrms=%.3f t[0]=%.3f t[1]=%.3f td=%.3f fdop=%.4f\n",phm,phc,rrms,t[0],t[1],td,fdop);
	eqv[30]=wph[0]*wph[1]; // work from smallest phase shift to minimise cumulative error
	for(cf=29;cf>=0;cf--) // DOES NOT VECTORIZE - history dependent
		eqv[cf]=eqv[cf+1]*wph[1]; // calculate negative frequency linear phase shift
//	wphm=cexpf(I*p[3] );
	wph[1]=1.0f/wph[1]; // comparable to wphm=cexpf(I*phm ); ?
	eqv[32]=wph[0]*wph[1]; // work from smallest phase shift to minimise cumulative error
	for(cf=33;cf<63;cf++) // DOES NOT VECTORIZE - history dependent
		eqv[cf]=eqv[cf-1]*wph[1]; // calculate negative frequency linear phase shift
	p[3]*=-NFFTdPI2; // td - real, not integer.

	p[2]=0.0f; // evm^2 total
	ptr=ffttmp+foff-31; // start at active part of PSCH/SSCH,
	for(cf=0;cf<63;cf++) // vectorized
		evm[cf]=ptr[cf]*eqv[cf]; // copy part of fft data
	for(cf=0;cf<63;cf++) // vectorized
		evm[cf]-=dict[cid][cf]; // calc EVM (zero DC ver of Dict)
	evm[31]=0.0;
	for(cf=0;cf<63;cf++) // vectorized
		evm[cf]*=conjf(evm[cf]); // convert to ||^2
	for(cf=0;cf<63;cf++) // vectorized
		p[2]+=crealf(evm[cf]);
	p[2]*=r62;
} // this is dependent on c2id code

// if we have an arbitrary window in time, we might see only part of the Sync, leading to SNR degradation.
// by using all subcarriers we can eliminate effect of noise
void LTEscan::MakeEqvVecFromXSCH2(short cid,short foff,float _Complex **dict,float **dictPh,float *p,float *sd)
{ //  // FFT X(w)exp(jwTD/NFFT) <--> x(t+TD), search delay for lowest err
	signed char cf=0;
	short cdly=0;
	float evm=0.0f;
//	float _Complex eqv[64];
	float phd[64];
	float bestEvm=3.0f;
	short bestDly=0;
	float bestDly2=0.0f;
	float *ptrPh=ffttmpPh+foff-31;;
//	float _Complex wph[2]={0.0}; // use linear phase model, ph = phm * fsc + phc
	if( (foff<36) || (NFFTm37<foff) )
	{
		printf("Range Fail MakeEqvVecFromXSCH2 foff=%i\n",foff);
		return;
	}	
	for(cf=0;cf<63;cf++) // vectorized
		phd[cf]=ptrPh[cf]-dictPh[cid][cf];// note using 'zero DC' version of dictionary
	for(cf=0;cf<63;cf++) // vectorized
		phd[cf]+=float(PI2)*((phd[cf]<float(-PI))-(float(PI)<phd[cf]));	
	p[4]=(phd[30]+phd[32])*0.5f; // normally average is more accurate
	// but if phase flip occurs, we need to check standard deviation to find best fit
	sd[0]=(phd[30]-p[4])*(phd[30]-p[4]); // SIMD friendly
	sd[1]=(phd[30]-p[4])*(phd[30]-p[4]);
	sd[2]=(phd[30]-p[4])*(phd[30]-p[4]);
	sd[0]+=(phd[32]-p[4])*(phd[32]-p[4]);
	sd[1]+=(phd[32]-p[4])*(phd[32]-p[4]);
	sd[2]+=(phd[32]-p[4])*(phd[32]-p[4]);
	sd[0]-=float(2*PI2)*(phd[32]-p[4]);
	sd[2]+=float(2*PI2)*(phd[32]-p[4]);
	sd[0]+=float(PI2*PI2); // was +=(phd[32]-p[4]-PI2)*(phd[32]-p[4]-PI2);
	sd[2]+=float(PI2*PI2); // was +==(phd[32]-p[4]+PI2)*(phd[32]-p[4]+PI2)
	sd[0]*=0.5f;
	sd[1]*=0.5f;
	sd[2]*=0.5f;
	if( (sd[1]>sd[0]) && (sd[2]>sd[0]) )
		p[4]-=float(PI); // was =(phd[30]+phd[32]-PI2)*0.5;
	if( (sd[1]>sd[2]) && (sd[2]>sd[0]) )
		p[4]+=float(PI); // was =(phd[30]+phd[32]+PI2)*0.5;
	
	//printf("[%.3f %.3f] {%.3f %.3f} [%.3f %.3f] %.3f\n",cargf(ptr[30]),cargf(ptr[32]),cargf(dict[cid][30]),cargf(dict[cid][32]),phd[30],phd[32],phc);
	p[4]+=float(PI2)*((p[4]<float(-PI))-(float(PI)<p[4])); // constrain to -PI:PI
	p[4]+=float(PI2)*((p[4]<float(-PI))-(float(PI)<p[4])); // constrain to -PI:PI
	for( cdly=-(NFFTd2);cdly<(NFFTd2);cdly+=3) // /4 = >>2, biggest step we can make is 3
	{
		evm=CalcFastEVM( phd,cdly,p[4] ); // evm is 4*sin2(phd/2) - we use cubic approx
		//printf("cdly=%i evm=%.3f phc=%f\n",cdly,evm,p[4]);
		if(bestEvm>evm)
		{
			bestEvm=evm;
			bestDly=cdly; // td
			//printf("BBB cdly=%i evm=%.5f phc=%.3f\n",bestDly,bestEvm,p[4]);
		}
	}
	evm=CalcFastEVM( phd,(bestDly-1),p[4] ); // evm is 4*sin2(phd/2) - we use cubic approx
	if(evm<bestEvm) // refine last two points
	{
		bestEvm=evm;
		bestDly-=1; // td
		//printf("BBB-cdly=%i evm=%.5f phc=%.3f\n",bestDly,bestEvm,p[4]);
	}
	evm=CalcFastEVM( phd,(bestDly+1),p[4] ); // evm is 4*sin2(phd/2) - we use cubic approx
	if(evm<bestEvm)
	{
		bestEvm=evm;
		bestDly+=1; // td
		//printf("BBB+cdly=%i evm=%.5f phc=%.3f\n",bestDly,bestEvm,p[4]);
	}
	// fractional search - would bisection be faster?
	float frac=0.0f;
	for(frac=-0.9f;frac<1.0f;frac+=0.1f)
	{
		evm=CalcFastEVM( phd,((float)bestDly+frac),p[4] ); // evm is 4*sin2(phd/2) - we use cubic approx
		if(evm<bestEvm)
		{
			bestEvm=evm;
			bestDly2=(float)bestDly+frac; // td
			//printf("BBB+cdly=%i evm=%.5f phc=%.3f\n",bestDly,bestEvm,p[4]);
		}
	}
	// bestDly2=bestDly;
		
	p[3]=-bestDly2*PI2dNFFT; // calculate 'zero dc' eq vector
/*	wph[0]=cexpf(-I*p[4] ); // small doppler looks like a constant phase shift
	wph[1]=cexpf(-I*p[3] );
	wph[0]*=rrms;
	eqv[31]=wph[0]; // DC has zero linear phase shift
	eqv[63]=1.0; // not used
	//printf("phm=%.3f phc=%.3f rrms=%.3f btd=%i berr=%.5f\n",p[3],p[4],rrms,bestDly,bestEvm);
	eqv[30]=wph[0]*wph[1]; // work from smallest phase shift to minimise cumulative error
	for(cf=29;cf>=0;cf--) // DOES NOT VECTORIZE - history dependent
		eqv[cf]=eqv[cf+1]*wph[1]; // calculate negative frequency linear phase shift
	wph[1]=1.0/wph[1]; // faster than //wphm=cexpf(I*phm ); ?
	eqv[32]=wph[0]*wph[1]; // work from smallest phase shift to minimise cumulative error
	for(cf=33;cf<63;cf++) // DOES NOT VECTORIZE - history dependent
		eqv[cf]=eqv[cf-1]*wph[1]; // calculate negative frequency linear phase shift */
	p[2]=bestEvm; // evm^2 ?
	p[3]=(float)bestDly2; // td
	//printf("MakeEqvVecFromXSCH2 cid=%i evm=%.5f %.5f td=%.1f phErr=%.3f\n",cid,100*sqrt(p[2]),p[2],p[3],p[4]);
}

void LTEscan::makeRotWt( float *ptr )
{ // assumes ptr[] is already constrained -PI:PI - not always true!
// move cyclic bracket, to give best chance to m/c calculator
// Whilst 4 versions may seem wasteful, it does improve accuracy in high noise.
// [0] -PI:PI, [1] -1.5*PI:0.5*PI, [2] -0.5*PI:1.5*PI, [3] 0:2*PI
	signed char cf=0;
	for(cf=0;cf<63;cf++)
	{
		rotWt[0][cf]=float(PI2)*((ptr[cf]<float(-PI))-(float(PI)<ptr[cf]));
		rotWt[1][cf]=float(PI2)*((ptr[cf]<float(-PI3d2))-(float(PId2)<ptr[cf])); // -1.5*PI to 0.5*PI
		rotWt[2][cf]=float(PI2)*((ptr[cf]<float(-PId2))-(float(PI3d2)<ptr[cf])); // -0.5*PI to 1.5*PI
		rotWt[3][cf]=float(PI2)*((ptr[cf]<0.0f)-(float(PI2)<ptr[cf])); //  0 to 2*PI
	}
}

void LTEscan::clampPh( float *ptr,unsigned char idx )
{ // [0] -PI:PI, [1] -1.5*PI:0.5*PI, [2] -0.5*PI:1.5*PI, [3] 0:2*PI
	signed char cf=0;
	if( idx==0 ) // [0] -PI:PI
		for(cf=0;cf<63;cf++)
			ptr[cf]+=float(PI2)*((ptr[cf]<float(-PI))-(float(PI)<ptr[cf]));
	if( idx==1 ) // [1] -1.5*PI:0.5*PI
		for(cf=0;cf<63;cf++)
			ptr[cf]+=float(PI2)*((ptr[cf]<float(-PI3d2))-(float(PId2)<ptr[cf]));
	if( idx==2 ) // [2] -0.5*PI:1.5*PI
		for(cf=0;cf<63;cf++)
			ptr[cf]+=float(PI2)*((ptr[cf]<float(-PId2))-(float(PI3d2)<ptr[cf]));
	if( idx==3 ) // [3] 0:2*PI
		for(cf=0;cf<63;cf++)
			ptr[cf]+=float(PI2)*((ptr[cf]<0.0f)-(float(PI2)<ptr[cf]));
}

void LTEscan::makeRotWtF( float *ptr )
{ // assumes ptr[] is already constrained -PI:PI - not always true!
	signed char cfm=(32+PTS_PRS); // (31+PTS_PRS+1)
	signed char cf=0;
	for(cf=(31-PTS_PRS);cf<cfm;cf++) // move cyclic bracket, to give best chance to m/c calculator
	{ // [0] -PI:PI, [1] -1.5*PI:0.5*PI, [2] -0.5*PI:1.5*PI, [3] 0:2*PI, wasteful, but improve accuracy in high noise.
		rotWt[0][cf]=float(PI2)*((ptr[cf]<float(-PI))-(float(PI)<ptr[cf]));
		rotWt[1][cf]=float(PI2)*((ptr[cf]<float(-PI3d2))-(float(PId2)<ptr[cf])); // -1.5*PI to 0.5*PI
		rotWt[2][cf]=float(PI2)*((ptr[cf]<float(-PId2))-(float(PI3d2)<ptr[cf])); // -0.5*PI to 1.5*PI
		rotWt[3][cf]=float(PI2)*((ptr[cf]<0.0f)-(float(PI2)<ptr[cf])); //  0 to 2*PI
	}
}

void LTEscan::clampPhF( float *ptr,unsigned char idx )
{ // [0] -PI:PI, [1] -1.5*PI:0.5*PI, [2] -0.5*PI:1.5*PI, [3] 0:2*PI
	unsigned char cfm=(32+PTS_PRS); // (31+PTS_PRS+1)
	unsigned char cf=0;
	if( (idx&1)==1 )
	{
		if( idx==1 ) // [1] -1.5*PI:0.5*PI
			for(cf=(31-PTS_PRS);cf<cfm;cf++)
				ptr[cf]+=float(PI2)*((ptr[cf]<float(-PI3d2))-(float(PId2)<ptr[cf]));
		else //if( idx==3 ) // [3] 0:2*PI
			for(cf=(31-PTS_PRS);cf<cfm;cf++)
				ptr[cf]+=float(PI2)*((ptr[cf]<0.0f)-(float(PI2)<ptr[cf]));
	}
	else
	{
		if( idx==0 ) // [0] -PI:PI
			for(cf=(31-PTS_PRS);cf<cfm;cf++)
				ptr[cf]+=float(PI2)*((ptr[cf]<float(-PI))-(float(PI)<ptr[cf]));
		else //if( idx==2 ) // [2] -0.5*PI:1.5*PI
			for(cf=(31-PTS_PRS);cf<cfm;cf++)
				ptr[cf]+=float(PI2)*((ptr[cf]<float(-PId2))-(float(PI3d2)<ptr[cf]));
	}
}

void LTEscan::fitF( float *ptr,float *prms,unsigned char idx )
{
	float ph[64];
	float data[64];
	float mod[64];//={0.0};
	unsigned char cfm=(32+PTS_PRS); // (31+PTS_PRS+1)
	unsigned char cf;
	prms[5]=0.0f; // variance (reduced range - model fit restricted to pair points)
	prms[4]=0.0f; // c;
	prms[3]=0.0f; // m;
//	prms[2]=0.0f; // variance (sd*sd);
	for(cf=(31-PTS_PRS);cf<cfm;cf++) // vectorizes
		data[cf]=ptr[cf]+rotWt[idx][cf];
	for(cf=(31-PTS_PRS);cf<cfm;cf++) // vectorizes
	{ // faster to do both as single loop
		prms[4]+=sumWt[cf]*data[cf]; // signed char * float - vectorizes as single loop
		prms[3]+=(rmWt[cf]*data[cf]); // float * float - fails to vectorize as single loop
	}
	prms[4]*=float(r2PTS_PRS);
	prms[3]*=float(r2PTS_PRS);
	if(prms[3]<-0.405f) // constrain m<=0.45, 0.45*31=2.22*PI2, worse case phase shift 3*PI2 => >|m|<=0.6
		prms[3]=-0.405f; // needs less constraint than main version as using a subset of the 62 points.
	if(0.405f<prms[3])
		prms[3]=0.405f;
	for(cf=(31-PTS_PRS);cf<cfm;cf++) // vectorizes
		mod[cf]=prms[4]+mWt[cf]*prms[3];
	clampPhF( mod,idx ); // constrain m<=0.45, 0.45*31=2.22*PI2, worse case phase shift 3*PI2 => >|m|<=0.6
	clampPhF( mod,idx );
//	if(fabsf(prms[3])>0.203)
//		clampPhF( mod,idx );
//	for(cf=0;cf<63;cf++)
	for(cf=(31-PTS_PRS);cf<cfm;cf++) // vectorizes
		ph[cf]=data[cf]-mod[cf];
	ph[31]=0.0f;
	for(cf=(31-PTS_PRS);cf<cfm;cf++) // vectorizes
		prms[5]+=ph[cf]*ph[cf]; // reduced range variance - restricted to pair points
//	for(cf=0;cf<63;cf++) // vectorizes
//		printf("%i\t%.3f\t%.3f\t%.3f\t%.3f\n",cf,ph[cf],data[cf],mod[cf],ptr[cf]);
//	printf("\n");
/*	if( (prms[5]*r2PTS_PRS) <0.25 )
	{
		printf("prms[5]=%.3f\n",(prms[5]*r2PTS_PRS));
		for(cf=(31-PTS_PRS);cf<cfm;cf++) // vectorizes
			printf("%i\t%.3f\t%.3f\t%.3f\t%.3f\n",cf,ph[cf],data[cf],mod[cf],ptr[cf]);
		printf("\n");
	}*/
}

void LTEscan::fit( float *ptr,float *prms,unsigned char idx )
{
	float ph[64];
	float data[64];
	float mod[64];
	signed char cfm=(32+PTS_PRS); // (31+PTS_PRS+1)
	signed char cf=0;
	prms[5]=0.0f; // variance (reduced range - model fit restricted to pair points)
	prms[4]=0.0f; // c;
	prms[3]=0.0f; // m;
	prms[2]=0.0f; // variance (sd*sd);
	for(cf=0;cf<63;cf++) // vectorizes
		data[cf]=ptr[cf]+rotWt[idx][cf];
	for(cf=(31-PTS_PRS);cf<cfm;cf++) // vectorizes
	{
		prms[4]+=sumWt[cf]*data[cf]; // signed char * float - vectorizes as single loop
		prms[3]+=(rmWt[cf]*data[cf]); // float * float - fails to vectorize as single loop
	}
	prms[4]*=float(r2PTS_PRS);
	prms[3]*=float(r2PTS_PRS);
//	printf("c=%.4f m=%.4f ",prms[4],prms[3]);
	if(prms[3]<-0.607f) // constrain m<=0.45, 0.45*31=2.22*PI2, worse case phase shift 3*PI2 => >|m|<=0.6
		prms[3]=-0.607f;
	if(0.607f<prms[3])
		prms[3]=0.607f;
	for(cf=0;cf<63;cf++) // vectorizes
		mod[cf]=prms[4]+mWt[cf]*prms[3];
	clampPh( mod,idx ); // constrain m<=0.45, 0.45*31=2.22*PI2, worse case phase shift 3*PI2 => >|m|<=0.6
	if(fabsf(prms[3])>0.203f)
		clampPh( mod,idx );
	if(fabsf(prms[3])>0.405f)
		clampPh( mod,idx );
	for(cf=0;cf<63;cf++) // vectorizes
		ph[cf]=(data[cf]-mod[cf]);
	ph[31]=0.0f;
	for(cf=0;cf<63;cf++) // vectorizes
		prms[2]+=ph[cf]*ph[cf];
	for(cf=(31-PTS_PRS);cf<cfm;cf++) // vectorizes
		prms[5]+=ph[cf]*ph[cf]; // reduced range variance - restricted to pair points
/*	if( (prms[5]*r2PTS_PRS) <0.25 )
	{
		printf("%.4f %.6f %.6f\n",prms[3],prms[2],(prms[5]*r2PTS_PRS));
		for(cf=0;cf<63;cf++) // vectorizes
			printf("%i\t%.6f\t%.3f\t%.3f\t%.3f\n",cf,(float)ph[cf],(float)data[cf],(float)mod[cf],(float)ptr[cf]);
		printf("\n");
	}*/
}

void LTEscan::MakeEqvVecFromXSCH3F(short cid,short foff,float _Complex **dict,float **dictPh,float *p,float *sd)
{ // 10-13us
	float phd[64];
	float prms[6]; // [0]=c,[1]=m,[2]=variance (sd*sd)
	float *ptrPh=ffttmpPh+foff-31;
	unsigned char cf;
	unsigned char ci;
	unsigned char cfm=(32+PTS_PRS); // (31+PTS_PRS+1) // faster to use single max variable than min/max array.
	if( (foff<36) || (NFFTm37<foff) )
	{
		printf("Range Fail MakeEqvVecFromXSCH3F foff=%i\n",foff);
		return;
	}
	memcpy(prms,p,2*sizeof(float)); // use p as best result, and use prms as scratch pad
	for(cf=(31-PTS_PRS);cf<cfm;cf++) // vectorizes
		phd[cf]=ptrPh[cf]-dictPh[cid][cf];
	for(cf=(31-PTS_PRS);cf<cfm;cf++) // vectorizes
	{ // using individual pointers does not speed up code
		rotWt[0][cf]=float(PI2)*((phd[cf]<float(-PI))-(float(PI)<phd[cf]));
		rotWt[1][cf]=float(PI2)*((phd[cf]<float(-PI3d2))-(float(PId2)<phd[cf])); // -1.5*PI to 0.5*PI
		rotWt[2][cf]=float(PI2)*((phd[cf]<float(-PId2))-(float(PI3d2)<phd[cf])); // -0.5*PI to 1.5*PI
		rotWt[3][cf]=float(PI2)*((phd[cf]<0.0f)-(float(PI2)<phd[cf])); //  0 to 2*PI
	}
	p[5]=100.0f;
	for(ci=0;ci<4;ci++)
	{ // search different brackets, and find best Variance
		fitF( phd,prms,ci ); // need to parse results
		if(prms[5]<p[5])
			memcpy(p,prms,6*sizeof(float)); // bstPrms[0]=prms[0] etc
	}
	p[5]*=float(r2PTS_PRS); // mean of phErr2
	p[4]+=float(PI2)*((p[4]<float(-PI))-(float(PI)<p[4]));
	p[2]=p[5];
	//printf("idx=%i c=%.4f m=%.6f td=%.4f var=%.8f var2=%.8f\n",foff,p[4],p[3],p[3]*(NFFTdPI2),p[2],p[5]);
	p[3]*=NFFTdPI2;
	float sdt=sqrtf(p[2]);
	//p[2]=sdt*(1-sdt*sdt/24); // conver to proper EVM number for test only.
//	printf("m=%.3f c=%.3f var=%.6f var2=%.6f sd=%.3f evm=%.1f%% snr=%.1fdB TD=%.3f\n",(PI2dNFFT*p[3]),p[4],p[2],p[5],sdt,200*sin(0.5*sdt),-10*log10(p[2]),p[3]);
//	printf("\n");
}

void LTEscan::MakeEqvVecFromXSCH3(short cid,short foff,float _Complex **dict,float **dictPh,float *p,float *sd)
{ // 13-18us
	float prms[8]; // [0]=c,[1]=m,[2]=variance (sd*sd)
	float phd[64];
	float *ptrPh=ffttmpPh+foff-31;
	signed char cf;
	signed char ci;
	if( (foff<36) || (NFFTm37<foff) )
	{
		printf("Range Fail MakeEqvVecFromXSCH3 foff=%i\n",foff);
		return;
	}
	memcpy(prms,p,6*sizeof(float)); // use p as best result, and use prms as scratch pad
	for(cf=0;cf<63;cf++) // vectorized
		phd[cf]=ptrPh[cf]-dictPh[cid][cf];
	makeRotWt( phd );
	p[5]=100.0f;
	for(ci=0;ci<4;ci++)
	{ // search different brackets, and find best Variance
		fit( phd,prms,ci ); // need to parse results
		if(prms[5]<p[5])
			memcpy(p,prms,6*sizeof(float)); // bstPrms[0]=prms[0] etc
	}
	p[2]*=r62; // r62
	p[5]*=float(r2PTS_PRS);
	p[4]+=float(PI2)*((p[4]<float(-PI))-(float(PI)<p[4]));
//	printf("idx=%i c=%.4f m=%.6f td=%.4f var=%.8f var2=%.8f\n",idx,prms[4],idx,prms[3],prms[3]*(NFFT*rPI2),idx,prms[2],prms[5]);
	p[3]*=NFFTdPI2;
	float sdt=sqrtf(p[2]);
//	printf("m=%.3f c=%.3f var=%.6f var2=%.6f sd=%.3f evm=%.1f%% snr=%.1fdB TD=%.3f\n",(PI2dNFFT*p[3]),p[4],p[2],p[5],sdt,200*sin(0.5*sdt),-10*log10(p[2]),p[3]);
//	printf("\n");
}

float LTEscan::CalcFastEVM( float *phd,float dly,float phc )
{ // treat signal as pure PSK (ignore Rx amplitude), problem is argf() is expensive
	float evm[64];
	float ph[64];
	unsigned char cf; // initialised in for loop
	float evmt=0.0f;
	float phm=-dly*PI2dNFFT;
	for(cf=0;cf<63;cf++) // vectorizes
		ph[cf]=phd[cf]+(phm*(cf-31)-phc); // typically phd=argf(ffttmp)-argf(dict[])
	for( cf=0;cf<63;cf++) // vectorizes, truncate to 0:(2*PI)
		ph[cf]-=float(PI2)*floorf(ph[cf]*float(rPI2)); // floor(ph[cf]*rPI2) int(ph[cf]*rPI2)
	for( cf=0;cf<63;cf++) // vectorizes, set range to -PI:PI
		ph[cf]-=float(PI2)*(ph[cf]>float(PI));
	ph[31]=0.0f;
	for(cf=0;cf<63;cf++) // vectorized
		evm[cf]=ph[cf]*(1-ph[cf]*ph[cf]*r24);// evm=2*sin(ph/2), apprx ph*(1-ph*ph/24), (acc 0:pi/2), >1.8  pi/2:pi.
	evm[31]=0.0f;
//	evm[63]=0.0f;
	for(cf=0;cf<63;cf++) // vectorized
		evmt+=evm[cf]*evm[cf];
//	for(cf=0;cf<63;cf++)
//		printf("%i %.3f %.3f %.3f %.6f\n",cf,phd[cf],ph[cf],evm[cf],evm[cf]*evm[cf]);
//	printf("%.6f\n",evmt*r62);
	return(evmt*r62);
}

void LTEscan::EVMpsch( int offset,signed char fftFIdxOff,short fL,short fH )
{ // Take FFT of a time block, evmt
	short cf=0;
	int cntCalcs=0;
	int cntDet=0;
	unsigned char c2id=0;
	float p[6]={0.0f}; // p[0]=pwr, p[1]=noise, p[2]=evm^2 p[3]=td p[4]=phErr
	float sd[3]={0.0f,0.0f,0.0f};
	clock_t mytref = clock();
	Blk2FFT( offset,fftFIdxOff,(1+fftFIdxOff) ); // approx 300us, NFFT 2048,  // fftFIdxOff={-1,0,1} *5kHz
	for( cf=fL;cf<=fH;cf+=20 ) // fL[fftFIdxOff+1] fH[fftFIdxOff+1]
	{
		CalcRMS62( cf,p ); // p[0], needed for Eqv
		rrms=1.0f/sqrtf(p[0]);
		for( c2id=0; c2id<3; c2id++ )
		{
			cntCalcs++;
/*			MakeEqvVecFromXSCH1A(c2id,cf,dictPSCHzdc,dictPSCHzdcPh,p,sd); // MakeEqvVecFromXSCH1 poor sens. MakeEqvVecFromXSCH2 too slow.
			if( p[2]<1.22f ) // 1.432 was lowest false alarm during testing.  Was 1.3
			{ // also see some false alarms when |cf|<36 away from valid code
				CalcRMS10( cf,p ); // p[1], needed for SNR
				tmpCellLog.AddLastCellPSCH(offset,c2id,cf,p); // -(int)(p[3]+0.5)
				//printf("cf=%i c2id=%i offset=%i evm=%.3f pwr=%.1f td=%.1f phErr=%.3f\n",cf,c2id,offset,100*sqrt(fabs(p[2])),10*log10(p[0]),p[3],p[4]);
//				printf("cBlk=%i cf=%i c2id=%i offset=%i evm=%.3f snr=%.1f td=%.1f phErr=%.3f\n",(offset/NFFTdKBLK),cf,c2id,offset,100*sqrtf(p[2]),10*log10(p[0]/p[1]),p[3],p[4]);
			}*/

			// uses only a few points to 'fast test' for coherence, this gives false alarms slowing down program
			// by using a second version to double check false alarms, can maintain high speed and reduce false alarms
			// still not as sensitive as the other routine.
			MakeEqvVecFromXSCH3F(c2id,cf,dictPSCHzdc,dictPSCHzdcPh,p,sd); // v3f 11us use to search for likely points
			//MakeEqvVecFromXSCH3(c2id,cf,dictPSCHzdc,dictPSCHzdcPh,p,sd); // v3 15us, use to double check 
			if( p[5]<1.44f ) // 1.432 was lowest false alarm during testing.  Was 1.3
			{ // also see some false alarms when |cf|<36 away from valid code
				cntDet++;
				//printf("cf=%i p[2]=%.8f\n",cf,p[2]);
				// we can use p[] to make this faster.
				MakeEqvVecFromXSCH3(c2id,cf,dictPSCHzdc,dictPSCHzdcPh,p,sd); // v3 15us, use to double check 
				if( p[2]<2.2f )
				{
					CalcRMS10( cf,p ); // p[1], needed for SNR
					tmpCellLog.AddLastCellPSCH(offset,c2id,cf,p); // -(int)(p[3]+0.5)
					printf("cf=%i p[2]=%.8f p[5]=%.8f\n",cf,p[2],p[5]);
//					printf("cf=%i c2id=%i offset=%i var=%.8f evm=%.3f pwr=%.1f td=%.1f phErr=%.3f\n",cf,c2id,offset,p[2],100*sqrtf(p[2]),10*log10(p[0]),p[3],p[4]);
//					printf("cBlk=%i cf=%i c2id=%i offset=%i evm=%.3f snr=%.1f td=%.1f phErr=%.3f\n",(offset/NFFTdKBLK),cf,c2id,offset,100*sqrtf(p[2]),10*log10(p[0]/p[1]),p[3],p[4]);
				}
			}

		}
	}
	//printf("cntDet=%i cntCalcs=%i\n",cntDet,cntCalcs);
	//printf("FFT + Post Proc took %.1fms\n",((double)(clock()-mytref))/CLOCKS_PER_mSEC);
}

void LTEscan::CheckCellPSCHevm( void )
{ // we add small timing error to ensure we are within CP, else EVM rises
	short cf=tmpCellLog.LoopGetFftFIdx(); // location in FFT
	unsigned char c2id=tmpCellLog.LoopGetN2ID();
	int offset=tmpCellLog.LoopGetTimePSCH(0)-(int)(tmpCellLog.LoopGetTdPSCH(0)+0.5f); // retime 
	signed char fftFIdxOff=tmpCellLog.LoopGetFftFIdxOff(); // 5kHz offsets
	unsigned char evmIdx=0;
	int FFTd4=32*NFFTd128;
	int FFT=105*NFFTd128; // was 420 (80% of FFT length), was 220 before that
	signed char fftOff=(1+fftFIdxOff);
	//printf("%i %i %i\n",offset,evmIdx,FrmSizeD2);
	Blk2FFT( offset,fftFIdxOff,fftOff ); // approx 300us, NFFT 2048
	CheckCellPSCHevmFoff( offset,evmIdx,cf,c2id );
	if( tmpCellLog.LoopGetEvmPSCH(0)>50 ) 
	{ // for poor SNR, search +/-300kHz with best timing to improve marginal detection rate
		if(offset>FFT) // check +/-0.8 FFT Length, if not too close to beginning of buffer.  Best EVM wins
			Blk2FFT( offset-FFT,fftFIdxOff,fftOff );
		if( (cf-20)>=36 )
			CheckCellPSCHevmFoff( offset,evmIdx,(cf-20),c2id ); // -300kHz=-20*15kHz
		Blk2FFT( offset+FFT,fftFIdxOff,fftOff );
		if( NFFTm37>=(cf+20) ) 
			CheckCellPSCHevmFoff( offset,evmIdx,(cf+20),c2id ); // +300kHz=+20*15kHz
		//if( cf!=tmpCellLog.LoopGetFftFIdx() )
		{ // if frequency changed, may need to retime
		//	printf("Retiming %i-->%i\n",offset,tmpCellLog.LoopGetTime()-tmpCellLog.LoopGetTdPSCH()-1 );
			offset=tmpCellLog.LoopGetTimePSCH(0)-(int)(tmpCellLog.LoopGetTdPSCH(0)+0.5f);
			cf=tmpCellLog.LoopGetFftFIdx();
			Blk2FFT( offset,fftFIdxOff,0 ); // approx 300us, NFFT 2048
			CheckCellPSCHevmFoff( offset,evmIdx,cf,c2id );
		}
	}
	if( tmpCellLog.LoopGetTdPSCH(0)>0 )
	{
		//printf("pos error %i %i\n",cf,tmpCellLog.LoopGetTdPSCH(0));
		offset=tmpCellLog.LoopGetTimePSCH(0)-tmpCellLog.LoopGetTdPSCH(0);
		Blk2FFT( offset,fftFIdxOff,fftOff ); // approx 300us, NFFT 2048
		CheckCellPSCHevmFoff( offset,evmIdx,cf,c2id );
	}
	if( (tmpCellLog.LoopGetTdPSCH(0)<(-FFTd4)) ) // 32*NFFTd128 =0.25 FFT length
	{
		//printf("neg error %i %i\n",cf,tmpCellLog.LoopGetTdPSCH(0));
		offset=tmpCellLog.LoopGetTimePSCH(0)-tmpCellLog.LoopGetTdPSCH(0);
		Blk2FFT( offset,fftFIdxOff,fftOff ); // approx 300us, NFFT 2048
		CheckCellPSCHevmFoff( offset,evmIdx,cf,c2id );
	}
	evmIdx=1; // check half frame, from best result
	offset+=(evmIdx*FrmSizeD2);
	Blk2FFT( offset,fftFIdxOff,fftOff ); // approx 300us, NFFT 2048
	cf=tmpCellLog.LoopGetFftFIdx();
	CheckCellPSCHevmFoff( offset,evmIdx,cf,c2id );
}

void LTEscan::CheckCellPSCHevmFoff( int offset,unsigned char evmIdx,short cf,unsigned char c2id )
{
	float p[8]={0.0f}; // p[0]=pwr, p[1]=noise, p[2]=evm^2 ph[3]=td ph[4]=phErr
	float sd[3]={0.0f};
	CalcRMS62( cf,p );
	CalcRMS10( cf,p );
	rrms=1.0f/sqrtf(p[0]);
	MakeEqvVecFromXSCH2(c2id,cf,dictPSCHzdc,dictPSCHzdcPh,p,sd);
	tmpCellLog.LoopUpdatePSCH(offset,cf,c2id,evmIdx,p); // best results with foff
	//printf("CheckCellPSCHevmFoff: evmIdx=%i offset=%i cf=%i c2id=%i td=%.2f evm=%.1f phErr=%.3f\n",
	//	evmIdx,offset,cf,c2id,p[3],100.0*sqrtf(p[2]),p[4] );
}

void LTEscan::CellSSCH( void )
{
	float p[8]={0.0f}; // p[0]=pwr, p[1]=noise, p[2]=evm^2 ph[3]=td ph[4]=phErr
	float sd[3]={0.0f};
	unsigned char cncp=0;
	short cf=tmpCellLog.LoopGetFftFIdx();
	unsigned char cns=0;
	short cid=0;
	unsigned char c1id=0;
	unsigned char c2id=tmpCellLog.LoopGetN2ID();
	long offset=tmpCellLog.LoopGetTimePSCH(0)-NFFT-cp1[0]; // PSCH, Assume NCP=0, FDD
	signed char fftFIdxOff=tmpCellLog.LoopGetFftFIdxOff();
	unsigned char bestNCP=0;
	unsigned char bestC1ID=0;
	unsigned char bestNs=0;
	int bestOffset=0;
	float bestP[8]={0.0f};
	bestP[2]=2.0f;
	//printf("CellSSCH\n");
	for( cncp=0;cncp<2;cncp++)
	{
		offset+=cncp*(cp1[0]-cp1[1]);
		Blk2FFT2( offset,fftFIdxOff,0 ); // approx 300us, NFFT 2048
		CalcRMS62( cf,p );
		CalcRMS10( cf,p );
		//printf("cncp=%i offset=%i p[0]=%.3f p[1]=%.3f\n",cncp,offset,p[0],p[1]);
		for( cns=0;cns<11;cns+=10 )
			for( c1id=0;c1id<168;c1id++ ) // 0:167
			{
				cid=(3*c1id+c2id); // SSCH depends on both N1ID and N2ID
				if(cns==0)
					MakeEqvVecFromXSCH2(cid,cf,dictSSCH0zdc,dictSSCH0zdcPh,p,sd);
				else // cns==10
					MakeEqvVecFromXSCH2(cid,cf,dictSSCH10zdc,dictSSCH10zdcPh,p,sd);
				// printf("cid=%i c1id=%i cs=%i cncp=%i time=%i p={%.3f %.3f %.3f %.3f %.3f}\n\n",cid,c1id,cns,cncp,offset,p[0],p[1],p[2],p[3],p[4]);
				if( bestP[2]>p[2] )
				{
					bestOffset=offset;
					bestC1ID=c1id;
					memcpy(bestP,p,sizeof(float)*5);
					bestNCP=cncp;
					bestNs=cns;
					// printf("CellSSCH cf=%i cid=%i c1id=%i cs=%i cncp=%i time=%i p={%.3f %.3f %.3f %.3f %.3f}\n",cf,cid,bestC1ID,bestNs,bestNCP,bestOffset,p[0],p[1],p[2],p[3],p[4]);
				}
			}
	}
	//printf("BestCellSSCH cf=%i cid=%i c1id=%i cs=%i cncp=%i time=%i p={%.3f %.3f %.3f %.3f %.3f}\n",cf,(3*bestC1ID+c2id),bestC1ID,bestNs,bestNCP,bestOffset,bestP[0],bestP[1],bestP[2],bestP[3],bestP[4]);
	tmpCellLog.LoopUpdateSSCH(bestOffset,bestC1ID,0,bestP,bestNs,bestNCP); // evmIdx
}

void LTEscan::FftShift( float _Complex *buffer )
{ // [-ve,0,+ve] <---> [0,ve,-ve], i.e. cyclically rotates spectrum by NFFT/2
	short cf;
	float _Complex tmp[1024];
	float _Complex *buffer2=buffer+NFFTd2;
	memcpy(tmp,buffer,sizeof(float _Complex)*NFFTd2);
	memcpy(buffer,buffer2,sizeof(float _Complex)*NFFTd2);
	memcpy(buffer2,tmp,sizeof(float _Complex)*NFFTd2);
}

void LTEscan::Blk2FFT( int offset,signed char fftFIdxOff,unsigned char sparseIdx )
{ // fftFIdxOff -1,0,1 * 5kHz, includes fftshift so DC is in middle (NFFT/2)
	unsigned short ct=0;
	unsigned short cf=0;
#ifdef USE_F32
	float _Complex *ptr=iqData+offset;
	for( ct=0;ct<NFFT;ct++ )
		ffttmp[ct]=ptr[ct]; // load depends of F32,I16,I12
#else
	int16_t *ptr=iqData+2*offset;
	for( ct=0;ct<NFFT;ct++ )
		ffttmp[ct]=ptr[2*ct]+I*ptr[2*ct+1]; // load depends of F32,I16,I12
#endif
	//printf("Blk2FFT=%i\n",fftFIdxOff);
	for( ct=0;ct<NFFT;ct++ )
		ffttmp[ct]*=fshift[(fftFIdxOff+1)][ct];
#ifdef USE_C99_COMPLEX
	for( ct=0;ct<NFFT;ct++ ) // vectorized
		in[ct]=ffttmp[ct];
#else
	for( ct=0;ct<NFFT;ct++ ) // vectorized
	{
		in[ct][0]=crealf(ffttmp[ct]);
		in[ct][1]=cimagf(ffttmp[ct]);
	}
#endif
	fftw_execute(pfft);
#ifdef USE_C99_COMPLEX // convert result to C99
	for( cf=0;cf<NFFT;cf++ ) // vectorized
		ffttmp[cf]=out[cf]; // convert to ||^2
#else // #ifdef USE_FFTW_COMPLEX
	for( cf=0;cf<NFFT;cf++ ) // vectorized
		ffttmp[cf]=out[cf][0]+I*out[cf][1]; // convert to ||^2
#endif
	FftShift(ffttmp); // fftshift here or in dictionary, but not both, to make ACF lineup with DC.
//	memcpy(fftacc,ffttmp,sizeof(float _Complex)*NFFT);
//	for( cf=0;cf<NFFT;cf++ ) // DOES NOT VECTORIZE
//		fftacc[cf]*=conjf(ffttmp[cf]);
//	for( cf=0;cf<NFFT;cf++ ) // DOES NOT VECTORIZE
//		ffttmpMg[cf]=crealf(fftacc[cf]);

//	for( cf=0;cf<NFFT;cf++ ) // DOES NOT VECTORIZE
//		fftacc[cf]=ffttmp[cf]*conjf(ffttmp[cf]); // convert to ||^2
//	for( cf=1;cf<NFFT;cf++ ) // DOES NOT VECTORIZE
//		ffttmpMg[cf]=crealf(fftacc[cf]); // bulk calc for search as cargf expensive

	for( cf=0;cf<NFFT;cf++ ) // DOES NOT VECTORIZE
		ffttmpMg[cf]=crealf(ffttmp[cf]*conjf(ffttmp[cf])); // convert to ||^2
	for( cf=1;cf<NFFT;cf++ ) // DOES NOT VECTORIZE - expensive
		ffttmpPh[cf]=cargf(ffttmp[cf]); // bulk calc for search as cargf expensive
}


void LTEscan::Blk2FFT2( int offset,signed char fftFIdxOff,unsigned char sparseIdx )
{
	unsigned short ct=0;
	unsigned short cf=0;
	float _Complex **sparsePtr=tmpCellLog.LoopGetSparseCor();
#ifdef USE_F32
	float _Complex *ptr=iqData+offset;
	for( ct=0;ct<NFFT;ct++ )
		ffttmp[ct]=ptr[ct]; // load depends of F32,I16,I12
#else
	int16_t *ptr=iqData+2*offset;
	for( ct=0;ct<NFFT;ct++ )
		ffttmp[ct]=ptr[2*ct]+I*ptr[2*ct+1]; // load depends of F32,I16,I12
#endif
	//printf("Blk2FFT2 sparseIdx=%i\n",sparseIdx);
	for( ct=0;ct<NFFT;ct++ )
		ffttmp[ct]*=sparsePtr[sparseIdx][ct];
#ifdef USE_C99_COMPLEX
	for( ct=0;ct<NFFT;ct++ ) // vectorized
		in[ct]=ffttmp[ct];
#else
	for( ct=0;ct<NFFT;ct++ ) // vectorized
	{
		in[ct][0]=crealf(ffttmp[ct]);
		in[ct][1]=cimagf(ffttmp[ct]);
	}
#endif
	fftw_execute(pfft);
#ifdef USE_C99_COMPLEX // convert result to C99
	for( cf=0;cf<NFFT;cf++ ) // vectorized
		ffttmp[cf]=out[cf]; // convert to ||^2
#else // #ifdef USE_FFTW_COMPLEX
	for( cf=0;cf<NFFT;cf++ ) // vectorized
		ffttmp[cf]=out[cf][0]+I*out[cf][1]; // convert to ||^2
#endif
	FftShift(ffttmp); // fftshift here or in dictionary, but not both, to make ACF lineup with DC.
//	memcpy(fftacc,ffttmp,sizeof(float _Complex)*NFFT);
//	for( cf=0;cf<NFFT;cf++ ) // DOES NOT VECTORIZE
//		fftacc[cf]*=conjf(ffttmp[cf]);
//	for( cf=0;cf<NFFT;cf++ ) // DOES NOT VECTORIZE
//		ffttmpMg[cf]=crealf(fftacc[cf]);

//	for( cf=0;cf<NFFT;cf++ ) // DOES NOT VECTORIZE
//		fftacc[cf]=ffttmp[cf]*conjf(ffttmp[cf]); // convert to ||^2
//	for( cf=1;cf<NFFT;cf++ ) // DOES NOT VECTORIZE
//		ffttmpMg[cf]=crealf(fftacc[cf]); // bulk calc for search as cargf expensive

	for( cf=0;cf<NFFT;cf++ ) // DOES NOT VECTORIZE
		ffttmpMg[cf]=crealf(ffttmp[cf]*conjf(ffttmp[cf])); // convert to ||^2
	for( cf=1;cf<NFFT;cf++ ) // DOES NOT VECTORIZE - expensive
		ffttmpPh[cf]=cargf(ffttmp[cf]); // bulk calc for search as cargf expensive
}

void LTEscan::EqualiseData( int toff,short foff,signed char fftFIdxOff,unsigned char NCP,unsigned char frm,unsigned char ns,unsigned char sym,float td,float phc,float rrms,float _Complex *data,unsigned char prtDiag )
{ // data[73]
	float _Complex eqv[73];
	float _Complex wphc=0.0f;
	float _Complex wphm=0.0f;
	short cf=0;
	float _Complex *ptr=ffttmp+foff-36; // 36<foff<(NFFT-36)
	int offset=toff+cp0[NCP];//+td;
	offset+=sym*(NFFT+cp1[NCP]);
	offset+=ns*(6+NCP)*(NFFT+cp1[NCP]);
	offset+=ns*(cp0[NCP]-cp1[NCP]); // start of actual OFDM symbol
	//printf("foff=%i toff=%i(%i) NCP=%i frm=%i ns=%i sym=%i td=%i ph%.3f rrms=%.3g\n",foff,toff,offset,NCP,frm,ns,sym,td,phc,rrms);
	if(offset>(NFFT*150*SRCH_FRMS))
		printf("WARNING: EqualiseData %i>%i\n",offset,(NFFT*150*SRCH_FRMS));
	Blk2FFT2( (offset),fftFIdxOff,CalcSparseIdx( frm,ns,sym,NCP ) ); // beginning of OFDM symbol	
	wphc=cexpf(-I*phc );
	wphm=cexpf(I*(td*PI2dNFFT) ); // phm=-td*PI2*rNFFT, wphm=cexpf(-I*phm );
	
	eqv[36]=rrms*wphc; // DC has zero linear phase shift
	if(prtDiag!=0)
		printf("phm=%.3f phc=%.3f rrms=%.3g td=%.3f\n",td*(PI2dNFFT),phc,rrms,td);

	eqv[35]=rrms*wphc*wphm; // work from smallest phase shift to minimise cumulative error
	for(cf=34;cf>=0;cf--) // DOES NOT VECTORIZE - history dependent
		eqv[cf]=eqv[cf+1]*wphm; // calculate negative frequency linear phase shift
	wphm=1.0f/wphm; // faster than //wphm=cexpf(I*phm ); ?
	eqv[37]=rrms*wphc*wphm; // work from smallest phase shift to minimise cumulative error
	for(cf=38;cf<73;cf++) // DOES NOT VECTORIZE - history dependent
		eqv[cf]=eqv[cf-1]*wphm; // calculate negative frequency linear phase shift 
	for(cf=0;cf<73;cf++)
		data[cf]=eqv[cf]*ptr[cf];
// convert to EVM
//	for(cf=0;cf<73;cf++)
//		printf("%i |%.3f| <%.3f |%.3f| <%.3f |%.3f| <%.3f\n",cf,cabsf(ptr[cf]),cargf(ptr[cf]),cabsf(eqv[cf]),cargf(eqv[cf]),cabsf(data[cf]),cargf(data[cf])); // calculate negative frequency linear phase shift  
}

void LTEscan::BulkPSCHSSCH( void )
{
	float _Complex data[292]; // PSCH,SSCH 36+DC+36 PBCH =4x73
	short cf=0;
	signed char fftFIdxOff=tmpCellLog.LoopGetFftFIdxOff();
	short foff=tmpCellLog.LoopGetFftFIdx();
	unsigned char NCP=tmpCellLog.LoopGetNCP();
	float td=tmpCellLog.LoopGetTdPSCH(0);
	float phErr=tmpCellLog.LoopGetPhErrPSCH(0);
	unsigned char ns=tmpCellLog.LoopGetNs(0);
	float rrms=tmpCellLog.LoopGetRrmsPSCH(0); // previously calculated
	unsigned char evmIdx=0;
	int toff=tmpCellLog.LoopGetTref();
	unsigned char frm=0;
	float _Complex *ptr=tmpCellLog.LoopGetPSCHptr();
	//printf("LTEscan::BulkPSCHSSCH\n");
	for(evmIdx=0;evmIdx<FERR_COR_HFRMS;evmIdx++) // update equalisation angles first
	{ // calc multiple PSCH and SSCH so that dphErrdt can be averaged to remove noise
		FastPSCHevm(toff,frm,ns,evmIdx); // calc evm,rms,snr etc. for each, but do not search again for NCP,NID
		FastSSCHevm(toff,frm,ns,evmIdx); // calc evm,rms,snr etc. for each, but do not search again for NCP,NID
		ns+=10;
		if(ns>10)
		{
			ns=0;
			frm++;
			toff+=FrmSize;
		}
	}	
	evmIdx=0;
	frm=0;
	ns=tmpCellLog.LoopGetNs(evmIdx);
	toff=tmpCellLog.LoopGetTref();
	ptr=tmpCellLog.LoopGetPSCHptr();
	td=tmpCellLog.LoopGetTdPSCH(evmIdx);
	phErr=tmpCellLog.LoopGetPhErrPSCH(evmIdx);
	rrms=tmpCellLog.LoopGetRrmsPSCH(evmIdx); // previously calculated
	EqualiseData( toff,foff,fftFIdxOff,NCP,frm,ns,(5+NCP),td,phErr,rrms,data,0 ); // fdd PSCH Ns=0
	// need EVM calculation based on phase and data
	//unsigned char N2ID=tmpCellLog.LoopGetN2ID();
	//printf("td=%.1f phErr=%.3f N2ID=%i Tref=%i ns=%i\n",td,phErr,N2ID,toff,ns);
	//CheckEVM(data,dictPSCHzdc,dictPSCHzdcPh,N2ID);
	for(cf=0;cf<63;cf++)
		ptr[cf]=data[5+cf]; // includes DC at ptr[31];
	ptr=tmpCellLog.LoopGetSSCHptr();
	td=tmpCellLog.LoopGetTdSSCH(evmIdx);
	phErr=tmpCellLog.LoopGetPhErrSSCH(evmIdx);
	rrms=tmpCellLog.LoopGetRrmsSSCH(evmIdx); // previously calculated
	EqualiseData( toff,foff,fftFIdxOff,NCP,frm,ns,(4+NCP),td,phErr,rrms,data,0 ); // fdd SSCH Ns=0
	unsigned short NID=tmpCellLog.LoopGetNID();
	if(ns==0) // updated Ns
		CheckEVM(data,dictSSCH0zdc,dictSSCH0zdcPh,NID);
	else // cns==10
		CheckEVM(data,dictSSCH10zdc,dictSSCH10zdcPh,NID);
	for(cf=0;cf<63;cf++)
		ptr[cf]=data[5+cf]; // includes DC at ptr[31];
}

void LTEscan::CheckEVM(float _Complex *ptr,float _Complex **dict,float **dictPh,short cid)
{ // diagnostic routine to close gap between 
// FastPSCHevm --> MakeEqvVecFromXSCH2 --> CalcFastEVM 
// EqualiseData --> CheckEVM
	float evm[64];
	float _Complex evm2[64];
	float phd[64];
	unsigned char cf; // initialised in for loop
	float evmt=0.0f;
	float evmt2=0.0f;
	for(cf=0;cf<63;cf++) // vectorized
		phd[cf]=cargf(ptr[cf+5]);
	for(cf=0;cf<63;cf++)  // vectorized
		phd[cf]-=dictPh[cid][cf]; // note using 'zero DC' version of dictionary
	for(cf=0;cf<63;cf++)  // vectorized
		phd[cf]+=float(PI2)*((phd[cf]<float(-PI))-(phd[cf]>float(PI)));
	phd[31]=0.0f;
	for(cf=0;cf<63;cf++) // vectorized
		evm[cf]=1-phd[cf]*phd[cf]*r24;
	for(cf=0;cf<63;cf++) // vectorized
		evm[cf]*=phd[cf]; // evm=2*sin( ph/2 ), approx ph*(1-ph*ph/24), (accurate to pi/2), remains above 1.8 to pi.
	evm[31]=0.0f;
	evm[63]=0.0f;
	for(cf=0;cf<63;cf++) // vectorized
		evmt+=evm[cf]*evm[cf];
	evmt*=r62;
	for(cf=0;cf<63;cf++) // vectorized
		evm2[cf]=ptr[cf+5]-dict[cid][cf];
	for(cf=0;cf<63;cf++) // vectorized
		evm2[cf]*=conjf(evm2[cf]);
	for(cf=0;cf<63;cf++) // vectorized
		evmt2+=crealf(evm2[cf]);
	evmt2*=r62;
	// ptr2
	// phd2
	// evmt3
	// dly=
	// float phDop=
	//evmt3=CalcFastEVM( float *phd2,float dly,float phDop );
//	for(cf=0;cf<63;cf++)
//		printf("%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",cf,phd[cf],cabsf(ptr[cf+5]),cargf(ptr[cf+5]),dictPh[cid][cf],evm[cf],sqrt(crealf(evm2[cf])));
//	printf("%.2f %.2f\n",sqrt(evmt)*100.0,sqrt(evmt2)*100.0);
}
void LTEscan::FindPBCH( void )
{
	unsigned char ns=tmpCellLog.LoopGetNs(0);
	signed char crc=0;
	unsigned char evmIdx=0;
	int toff=tmpCellLog.LoopGetTref();
//	int toff=tmpCellLog.LoopGetTref()-36; // 36 is for test only
	unsigned char frm=0;
	printf("FindPBCH\n");
	
// test routine to compare 
	short foff=tmpCellLog.LoopGetFftFIdx();	
	signed char fftFIdxOff=tmpCellLog.LoopGetFftFIdxOff();
	unsigned char NID=tmpCellLog.LoopGetNID();
	unsigned char N2ID=tmpCellLog.LoopGetN2ID();
	unsigned char NCP=tmpCellLog.LoopGetNCP();
	unsigned char sym=6; // fdd
	float p[6]={0.0f};
	float sd[3]={0.0f};
	int offset=toff+cp0[NCP];
	offset+=sym*(NFFT+cp1[NCP]);
	offset+=ns*(6+NCP)*(NFFT+cp1[NCP]); // was 6 for PSCH, something else for junk
	offset+=ns*(cp0[NCP]-cp1[NCP]);
	Blk2FFT2( offset,fftFIdxOff,CalcSparseIdx( frm,ns,(5+NCP),NCP ) ); // approx 300us, NFFT 2048
	CalcRMS62( foff,p ); // --> rms
	CalcRMS10( foff,p ); // --> snr
	clock_t mytref2;
	clock_t mytref = clock();
	MakeEqvVecFromXSCH3F(N2ID,foff,dictPSCHzdc,dictPSCHzdcPh,p,sd);
	mytref2=clock();
	printf("v3F p[2]=%.6f p[3]=%.3f p[4]=%.3f p[5]=%.6f %.3fms\n",p[2],p[3],p[4],p[5],((double)(mytref2-mytref))/CLOCKS_PER_mSEC);
	Blk2FFT2( offset,fftFIdxOff,CalcSparseIdx( frm,ns,(5+NCP),NCP ) ); // approx 300us, NFFT 2048
	CalcRMS62( foff,p ); // --> rms
	CalcRMS10( foff,p ); // --> snr
	mytref=clock();
	MakeEqvVecFromXSCH3(N2ID,foff,dictPSCHzdc,dictPSCHzdcPh,p,sd);
	mytref2=clock();
	printf("v3 p[2]=%.6f p[3]=%.3f p[4]=%.3f p[5]=%.6f %.3fms\n",p[2],p[3],p[4],p[5],double(mytref2-mytref)/CLOCKS_PER_mSEC);
	Blk2FFT2( offset,fftFIdxOff,CalcSparseIdx( frm,ns,(5+NCP),NCP ) ); // approx 300us, NFFT 2048
	CalcRMS62( foff,p ); // --> rms
	CalcRMS10( foff,p ); // --> snr
	rrms=tmpCellLog.LoopGetRrmsPSCH(0);
	mytref = clock();
	MakeEqvVecFromXSCH1(N2ID,foff,dictPSCHzdc,dictPSCHzdcPh,p,sd);
	mytref2=clock();
	printf("v1 p[2]=%.6f p[3]=%.3f p[4]=%.3f p[5]=%.6f %.3fms\n",p[2],p[3],p[4],p[5],((double)(mytref2-mytref))/CLOCKS_PER_mSEC);
	Blk2FFT2( offset,fftFIdxOff,CalcSparseIdx( frm,ns,(5+NCP),NCP ) ); // approx 300us, NFFT 2048
	CalcRMS62( foff,p ); // --> rms
	CalcRMS10( foff,p ); // --> snr
	rrms=tmpCellLog.LoopGetRrmsPSCH(0);
	mytref = clock();
	MakeEqvVecFromXSCH1A(N2ID,foff,dictPSCHzdc,dictPSCHzdcPh,p,sd);
	mytref2=clock();
	printf("v1A p[2]=%.6f p[3]=%.3f p[4]=%.3f p[5]=%.6f %.3fms\n",p[2],p[3],p[4],p[5],((double)(mytref2-mytref))/CLOCKS_PER_mSEC);
	Blk2FFT2( offset,fftFIdxOff,CalcSparseIdx( frm,ns,(5+NCP),NCP ) ); // approx 300us, NFFT 2048
	CalcRMS62( foff,p ); // --> rms
	CalcRMS10( foff,p ); // --> snr
	rrms=tmpCellLog.LoopGetRrmsPSCH(0);
	mytref = clock();
	MakeEqvVecFromXSCH1B(N2ID,foff,dictPSCHzdc,dictPSCHzdcPh,p,sd);
	mytref2=clock();
	printf("v1B p[2]=%.6f p[3]=%.3f p[4]=%.3f p[5]=%.6f %.3fms\n",p[2],p[3],p[4],p[5],((double)(mytref2-mytref))/CLOCKS_PER_mSEC);
	Blk2FFT2( offset,fftFIdxOff,CalcSparseIdx( frm,ns,(5+NCP),NCP ) ); // approx 300us, NFFT 2048
	CalcRMS62( foff,p ); // --> rms
	CalcRMS10( foff,p ); // --> snr
	mytref = clock();
	MakeEqvVecFromXSCH2(N2ID,foff,dictPSCHzdc,dictPSCHzdcPh,p,sd);
	mytref2=clock();
	printf("v2 p[2]=%.6f p[3]=%.3f p[4]=%.3f p[5]=%.6f %.3fms\n",p[2],p[3],p[4],p[5],((double)(mytref2-mytref))/CLOCKS_PER_mSEC);

	LoadRefs(1,0,NID,NCP); // PBCH ns=1,l=0 for FDD.
	for(evmIdx=0;evmIdx<SRCH_HFRMS;evmIdx++)
	{ // PBCH can only be decoded in 1/4 PBCHs, so must read up to x4 PBCH to decode
		if( (ns==0)&&(crc!=1) ) // search for first successful PBCH
			crc=FastPBCHevm(toff,frm,evmIdx);
		printf("PBCH toff=%i frm=%i ns=%i crc=%i\n",toff,frm,ns,crc);
		ns+=10;
		if(ns>10)
		{
			ns=0;
			toff+=FrmSize;
			frm++;
		}
	}
}

void LTEscan::FastPSCHevm( int toff,unsigned char frm,unsigned char ns,unsigned char evmIdx )
{
	float p[8]={0.0f};
	float sd[3]={0.0f};
	short foff=tmpCellLog.LoopGetFftFIdx();	
	signed char fftFIdxOff=tmpCellLog.LoopGetFftFIdxOff();
	unsigned char N2ID=tmpCellLog.LoopGetN2ID();
	unsigned char NCP=tmpCellLog.LoopGetNCP();
	unsigned char sym=6; // fdd
	int offset=toff+cp0[NCP];
	offset+=sym*(NFFT+cp1[NCP]);
	offset+=ns*(6+NCP)*(NFFT+cp1[NCP]);
	offset+=ns*(cp0[NCP]-cp1[NCP]);
	//printf("FastPSCHevm\n");
//	printf("FastPSCHevm toff=%i ns=%i evmIdx=%i NCP=%i NID=%i\n",toff,ns,evmIdx,NCP,N2ID);
//	if(evmIdx<1) // already done
//		return;
	Blk2FFT2( offset,fftFIdxOff,CalcSparseIdx( frm,ns,(5+NCP),NCP ) ); // approx 300us, NFFT 2048
	CalcRMS62( foff,p ); // --> rms
	CalcRMS10( foff,p ); // --> snr
	// rrms=tmpCellLog.LoopGetRrmsSSCH(evmIdx);
	// EqualiseData( toff,foff,fftFIdxOff,NCP,ns,(5+NCP),td,phErr,rrms,data ); // fdd SSCH Ns=0
	MakeEqvVecFromXSCH2(N2ID,foff,dictPSCHzdc,dictPSCHzdcPh,p,sd);
	tmpCellLog.LoopUpdatePSCH(offset,foff,N2ID,evmIdx,p); // best results with foff
	//printf("FastPSCHevm toff=%i ns=%i evmIdx=%i NCP=%i N2ID=%i  %f %f %f %f\n",toff,ns,evmIdx,NCP,N2ID,p[4],sd[0],sd[1],sd[2]);
}

void LTEscan::FastSSCHevm( int toff,unsigned char frm,unsigned char ns,unsigned char evmIdx )
{
	float p[8]={0.0f};
	float sd[3]={0.0f};
	short foff=tmpCellLog.LoopGetFftFIdx();	
	signed char fftFIdxOff=tmpCellLog.LoopGetFftFIdxOff();
	unsigned short NID=tmpCellLog.LoopGetNID();
	unsigned char N1ID=tmpCellLog.LoopGetN1ID();
	unsigned char NCP=tmpCellLog.LoopGetNCP();
	unsigned char sym=5; // fdd
	int offset=toff+cp0[NCP];
	offset+=sym*(NFFT+cp1[NCP]);
	offset+=ns*(6+NCP)*(NFFT+cp1[NCP]);
	offset+=ns*(cp0[NCP]-cp1[NCP]);
	//printf("FastSSCHevm\n");
//	printf("FastSSCHevm toff=%i ns=%i evmIdx=%i NCP=%i NID=%i\n",toff,ns,evmIdx,NCP,NID);
//	if(evmIdx<1) // already done
//		return;
	Blk2FFT2( offset,fftFIdxOff,CalcSparseIdx( frm,ns,(4+NCP),NCP ) ); // approx 300us, NFFT 2048
	CalcRMS62( foff,p ); // --> rms
	CalcRMS10( foff,p ); // --> snr
	if(ns==0) // updated Ns
		MakeEqvVecFromXSCH2(NID,foff,dictSSCH0zdc,dictSSCH0zdcPh,p,sd);
	else // cns==10
		MakeEqvVecFromXSCH2(NID,foff,dictSSCH10zdc,dictSSCH10zdcPh,p,sd);
	tmpCellLog.LoopUpdateSSCH(offset,N1ID,evmIdx,p,ns,NCP );
	//printf("FastSSCHevm toff=%i ns=%i evmIdx=%i NCP=%i NID=%i  %f %f %f %f\n",toff,ns,evmIdx,NCP,NID,p[4],sd[0],sd[1],sd[2]);
}

signed char LTEscan::FastPBCHevm( int toff,unsigned char frm,unsigned char evmIdx )
{
	char *pbchMsg=tmpCellLog.LoopGetPBCHmsgPtr();
	short cf=0;
	unsigned char csym=0;
	float _Complex data[74];
	float _Complex data2[292]; // PSCH,SSCH 36+DC+36 PBCH =4x73
	float _Complex *data2ptr=data2;
	signed char crc=0;
	float _Complex *ptrPBCH=tmpCellLog.LoopGetPBCHptr(); // array to store PBCH
	float _Complex *ptrPilotRaw=tmpCellLog.LoopGetPBCHpilotRawPtr(); // array to store PBCH Pilots
	float _Complex *ptrPilotDeSpread=tmpCellLog.LoopGetPBCHpilotDeSpreadPtr(); // array to store despread PBCH Pilots
	signed char fftFIdxOff=tmpCellLog.LoopGetFftFIdxOff();
	short foff=tmpCellLog.LoopGetFftFIdx();	// 0 at most negative frequency
	short fftIdx=tmpCellLog.LoopGetFftFIdx()-NFFTd2; // signed relative to DC
	float rrms=tmpCellLog.LoopGetRrmsPSCH(evmIdx); // Use PSCH to amplitude equalize PBCH
	float phErr=tmpCellLog.LoopGetPhErrPSCH(evmIdx);
	//float dphErrdt=phErr-tmpCellLog.LoopGetPhErrSSCH(evmIdx);
	float td=tmpCellLog.LoopGetTdPSCH(evmIdx);
	unsigned char NCP=tmpCellLog.LoopGetNCP();
	unsigned short NID=tmpCellLog.LoopGetNID();	
	float dphErrdt=tmpCellLog.LoopGetdpErrdtAv(); // average dphErrdt to remove noise
	//printf("FastPBCHevm\n");
	phErr+=dphErrdt*cp0pNFFTdcp1pNFFT[NCP]-0.5f*cp0mcp1dNFFT_2pi[NCP]; // seems to need both?
	phErr=CalcSlotBndryPhShft( phErr,fftIdx,NCP );
	for(csym=0;csym<4;csym++) // always 4 OFDM symbols, regardless of NCP.  But pilots change with P,NCP
	{ // PBCH always in slot 1 for FDD
		EqualiseData( toff,foff,fftFIdxOff,NCP,frm,1,csym,td,phErr,rrms,data,1 );
		//for(cf=0;cf<73;cf++)
		//	data2[csym*73+cf]=data[cf];
		memcpy(data2ptr,data,sizeof(data2ptr)*73);
		data2ptr+=73;
		phErr+=dphErrdt;
		phErr+=float(PI2)*((phErr<float(-PI))-(phErr>float(PI)));	
	}
	crc=DemodPBCH( data2,NCP,NID,pbchMsg ); // --> pbchEVM
	if(tmpCellLog.LoopGetPBCHcrc()!=1)
	{
		tmpCellLog.UpdatePBCH(pbchEVM,P,bw,phichDur,phichRes,sfn,crc,evmIdx,errAvPBCH,pbchSymbCnt);
		for(cf=0;cf<pbchSymbCnt;cf++)
			ptrPBCH[cf]=d[cf];
		for(cf=0;cf<12;cf++)
			ptrPilotRaw[cf]=pilotsRaw[2*cf];
		for(cf=0;cf<12;cf++)
			ptrPilotDeSpread[cf]=pilotsDespread[2*cf];
		for(cf=0;cf<12;cf++)
			ptrPilotRaw[cf+12]=pilotsRaw[2*cf+1];
		for(cf=0;cf<12;cf++)
			ptrPilotDeSpread[cf+12]=pilotsDespread[2*cf+1];
		for(cf=0;cf<12;cf++)
			ptrPilotRaw[cf+24]=pilotsRaw[2*cf+24];
		for(cf=0;cf<12;cf++)
			ptrPilotDeSpread[cf+24]=pilotsDespread[2*cf+24];
		for(cf=0;cf<12;cf++)
			ptrPilotRaw[cf+36]=pilotsRaw[2*cf+25];
		for(cf=0;cf<12;cf++)
			ptrPilotDeSpread[cf+36]=pilotsDespread[2*cf+25];
	}
	return(crc);
}

float LTEscan::CalcSlotBndryPhShft( float ph,short foff,unsigned char NCP )
{ // extra phase shift due to CP0>CP1 when NCP=1
	float phErr=ph;
	float d=0.0f;
	phErr+=cp0mcp1dNFFT_2pi[NCP]*foff;
	d=floorf(phErr*float(rPI2));
	phErr-=d*float(PI2);
	phErr+=float(PI2)*((phErr<float(-PI))-(phErr>float(PI)));	
	phErr+=float(PI2)*((phErr<float(-PI))-(phErr>float(PI)));	
	return(phErr);
}

unsigned char LTEscan::CalcSparseIdx( unsigned char frm,unsigned char ns,unsigned char sym,unsigned char NCP )
{ // 0 SSCH, 1 PSCH, 2-5 PBCHX.0-3  X is relFrm
	unsigned char sparseIdx=12*frm; // LUT is based on 6 elements per half frame
	sparseIdx+=6*(ns>9);
	if( (ns&1)==0 )
		sparseIdx+=sym-(4+NCP);
	else
		sparseIdx+=2+sym;
	//printf("frm=%i ns=%i sym=%i NCP=%i sparseIdx=%i\n",frm,ns,sym,NCP,sparseIdx);
	return( sparseIdx );
}
