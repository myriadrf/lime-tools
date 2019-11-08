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
#include "lime.h"
#ifdef USE_C99_COMPLEX
#include <complex.h> // double _Complex, compliant with double, fftw_complex[2] Re=0, Im=1
#include <fftw3.h>  // should come after complex.h
#else // #ifdef USE_FFTW_COMPLEX
#include <fftw3.h>  // Do before complex.h to force typedef double fftw_complex[2]
#include <complex.h> 
#endif
#ifdef __cplusplus
extern "C" {
#endif
#include "turbofec/conv.h"
#include "turbofec/rate_match.h"
#include "turbofecLib.h"
#ifdef __cplusplus // If this is a C++ compiler, end C linkage
}
#endif
#include "CellItem.hpp"
#include "CellList.hpp"
#include "LTEscan.hpp"

unsigned char pbcha[40];
long nSysFrm=0;
unsigned char Pant[3]={1,2,4};
extern unsigned short NFFT; // default number of FFT bins
extern unsigned short NFFTd128;

// TS 136.211 6.6.1: b(i) data bits, xor scrabled with code c(i) --> btilda bits i=0:(Mbit-1)
// TS 136.211 6.6.2: btilda bits modulated on d(i) qpsk bits i=0:(Msymb-1)
// TS 136.211 6.6.3: layer mapping 6.3.3.1 or 6.3.3.3
// TS 136.211 6.6.3: precoding 6.3.4.1 or 6.3.4.3 y(p)(i) p=Antenna port p=0:(P-1), i=sample number i=0:(Msymb-1) P={1,2,4}
// note 6.3.3.1 => 6.3.4.1, and 6.3.3.3 => 6.3.4.3
// TS 136.211 6.3.3.1 mapping (layers v=1) x(0)(i)=d(i) - do nothing
// TS 136.211 6.3.3.3 mapping (layers v=2=P) x(0)(i)=d(i), x(1)(i)=d(i+1)
// TS 136.211 6.3.4.1 precoding (P=1)  y(p)(i)=x(0)(i)
// TS 136.211 6.3.4.3 precoding (P=2) y(0)(2i)  = Re[x(0)(i)]+jIm[x(0)(i)] => y(0)(2i)  = x(0)(i)  => y(0)(2i)  = d(i)
// TS 136.211 6.3.4.3 precoding (P=2) y(1)(2i)  =-Re[x(1)(i)]+jIm[x(1)(i)] => y(1)(2i)  =-x(1)(i)* => y(1)(2i)  =-d(i+1)*
// TS 136.211 6.3.4.3 precoding (P=2) y(0)(2i+1)= Re[x(1)(i)]+jIm[x(1)(i)] => y(0)(2i+1)= x(1)(i)  => y(0)(2i+1)= d(i+1)
// TS 136.211 6.3.4.3 precoding (P=2) y(1)(2i+1)= Re[x(0)(i)]-jIm[x(0)(i)] => y(1)(2i+1)= x(0)(i)* => y(1)(2i+1)= d(i)*
// more complex precoding exists for P=4, we do not consider that here.

// so we have three cases for MIMO: 
// (i) We only receive from p=0.  This contains PSCH, SSCH
// y(0)(0)=d(0),y(0)(1)=d(1), i.e. looks like SISO data stream.
// (ii) We only receive from p=1.  This contains PSCH, SSCH
// y(1)(0)=-d(1)*,y(1)(1)=d(0)*
// (iii) We receive from p=0 and p=1 and they are summed with some ratio {|a|,|b|} with arbitrary phase.
// y(0+1)(0)=d(0)-b*d(1)*, y(0+1)(1)=d(1)+b*d(0)*
// => d(0)*(1+b*b)

// calculate qpsk pilot sequence prior to decoding PBCH.  Is independent of BW and P, but does depend on NID,NCP.
// note the qpsk pilots (references) are unmixed,
// despreading pilots gives a single vector, so we can measure rms amplitude and phase of each antenna
// we only have 12 pilots, so only 11dB of code gain, so less SNR tolerant than SSCH and PSCH.

// original MIMO precode
// y[0][0]= d[0]
// y[0][1]= d[1]
// y[1][0]=-d[1]'
// y[1][1]= d[0]'
// both y[0][] and y[1][] are received by a single antenna p, but experience different path delay and loss
// let y[0][] be scaled by complex factor a, its value is known from pilots.
// let y[1][] be scaled by complex factor b, its value is known from pilots.
// p[0]=a*y[0][0]+b*y[1][0]=a*d[0]-b*d[1]'
// p[1]=a*y[0][1]+b*y[1][1]=a*d[1]+b*d[0]'
// note conjugate identies:
// (a+b)'=a'+b'
// (a.b)'=a'.b'
// a.a'=|a|^2
// p[0] + b/a'.p[1]' = a*d[0]-b*d[1]' + b/a'.(a*d[1]+b*d[0]')' = a*d[0]-b*d[1]'+ b/a'.(a'*d[1]'+b'*d[0]) = d[0](a+b/a'b')
// p[1] - b/a'.p[0]' = a*d[1]+b*d[0]' - b/a'.(a*d[0]-b*d[1]')' = a*d[1]+b*d[0]'- b/a'.(a'*d[0]'-b'*d[1]) = d[1](a+b/a'b')
// 1.0/(a+b/a'b')=a'/(|a|2+|b|2)
signed char LTEscan::DemodPBCH( float _Complex *eqData,unsigned char NCP,short NcellID,char *msg )
{ // pbch --> d --> pbchData,evmPBCH
	unsigned short vshift=(NcellID % 6); // cell specific frequency shift of reference pilots
	short refCnt=0;
	// d[] usually 252 for NCP=0, 264 for NCP=1 => E can be larger for better SNR? NCP=1 D=240, NCP=0 D=216
	unsigned char D=2*96; // number of QPSK symbols - was 96 (2*96=192 <216)
	unsigned char E=2*D; // number of encoded bits - defined by the Mux size 3x32 rows. (min data 3x40=120)
	float evmt=0.0;
	unsigned char csym=0;
	short cf=0;
	unsigned char ci=0;
	unsigned char rv=0; // does not affect PBCH
	unsigned char cp=0;
	// max 4 ant's p={0,1} ns=1,l=0, p={2,3} ns=1,l=1, p={0,1} repeats in l=3 NCP=0 => 6*12
	float refPhi[4]={0.0,0.0,0.0,0.0};
	float refMag[4]={0.0,0.0,0.0,0.0};
	float _Complex fac=1.0f;

	pbchSymbCnt=0;
	for(csym=0;csym<4;csym++)  // calculation independent of P antennas
	{ // always 4 OFDM symbols, regardless of NCP.  Pilots change with P,NCP
		if((csym<2) || ((csym==3)&&(NCP==0)) || ((P==4)&&(NCP==0)) ) // Assume P=={1,2}
		{	
		// v=0		
		//k1=((v+vshift) % 6);
		//k2=((v+vshift+3) % 6);
		// if csym > 2, v=3
			for(cf=0;cf<36;cf++) // negative spectrum
				if( ((cf+1+vshift)%3)>0 ) // strip pilots, was +2
					d1[pbchSymbCnt++]=eqData[csym*73+cf];
				else
					pilotsRaw[refCnt++]=eqData[csym*73+cf];
			for(cf=37;cf<73;cf++) // positive spectrum
				if( ((cf-36+vshift)%3)>0 ) // strip pilots, was -35
					d1[pbchSymbCnt++]=eqData[csym*73+cf]; // ffttmp includes DC, so 36+1
				else
					pilotsRaw[refCnt++]=eqData[csym*73+cf];
		}
		else // no pilots
		{
			for(cf=0;cf<36;cf++) // negative spectrum
				d1[pbchSymbCnt++]=eqData[csym*73+cf];
			for(cf=37;cf<73;cf++) // positive spectrum
				d1[pbchSymbCnt++]=eqData[csym*73+cf]; // ffttmp includes DC, so 36+1
		}
	}
	//printf("pbchSymbCnt=%i refCnt=%i\n",pbchSymbCnt,refCnt); // P=2 NCP=1 pbchSymbCnt=240 refCnt=48,NCP=0 216,72

	if( vshift<3 ) // calculation independent of P antennas
	{
		refPhi[0]=RefAvPh( pilotsRaw );
		refMag[0]=RefAvMg( pilotsRaw );
		refPhi[1]=RefAvPh( &pilotsRaw[1] );
		refMag[1]=RefAvMg( &pilotsRaw[1] );
		refPhi[2]=RefAvPh( &pilotsRaw[24] );
		refMag[2]=RefAvMg( &pilotsRaw[24] );
		refPhi[3]=RefAvPh( &pilotsRaw[25] );
		refMag[3]=RefAvMg( &pilotsRaw[25] );
	}
	else
	{
		refPhi[1]=RefAvPh( pilotsRaw );
		refMag[1]=RefAvMg( pilotsRaw );
		refPhi[0]=RefAvPh( &pilotsRaw[1] );
		refMag[0]=RefAvMg( &pilotsRaw[1] );
		refPhi[3]=RefAvPh( &pilotsRaw[24] );
		refMag[3]=RefAvMg( &pilotsRaw[24] );
		refPhi[2]=RefAvPh( &pilotsRaw[25] );
		refMag[2]=RefAvMg( &pilotsRaw[25] );
	}
	DespreadPilots();	
	LTESC( cprsPBCH,NcellID,E ); // descramble
	for(cp=0;cp<3;cp++)
	{
		if( cp>0 ) // create 2x2 transmit diversity 'precoded' version of d
		{
			fac=(refMag[1]/refMag[0])*cexpf(I*(refPhi[1]+refPhi[0])); // b/a'
			for(cf=0;cf<D;cf+=2)
			{ // p[cf]=d[cf]-k*d2[cf]*.  k=0 SISO, k=1 2x2 MIMO, in real life MIMO leakage may not be 1:1.
				d2[cf]=fac*conjf(d1[cf+1]);
				d2[cf+1]=-fac*conjf(d1[cf]);
			}
	
			fac=(refMag[0]*cexpf(-I*refPhi[0]))/(refMag[1]*refMag[1]+refMag[0]*refMag[0]); // a'/(|a|2+|b|2)=1.0/(a+b/a'b')
			for(cf=0;cf<D;cf++)
				d[cf]=fac*(d1[cf]+d2[cf]); // combine orthogonal MIMO channel data
		}
		else // SISO, nothing to do!
		{
			fac=1.0f;
			fac=(1.0f/refMag[0])*cexpf(-I*refPhi[0]);
			for(cf=0;cf<D;cf++)
				d[cf]=fac*(d1[cf]); // combine orthogonal MIMO channel data
		}

		for(csym=0;csym<D;csym++) // Soft QPSK decoder - vectorized
		{
			ee[csym*2]=(signed char)(-64*crealf(d[csym])); // msb
			ee[csym*2+1]=(signed char)(-64*cimagf(d[csym])); // lsb
		}
		for(csym=0;csym<D;csym++) // Fast EVM calculator - based on phase angle only
			evm[csym]=(cargf(d[csym])-qpskPh[2*(ee[csym*2]>0)+(ee[csym*2+1]>0)]);

		errAvPBCH=0.0f;
		for(csym=0;csym<D;csym++)
			errAvPBCH+=evm[csym];
		errAvPBCH/=D; // calculate average phase error of PBCH relative to PSCH/SSCH for debug

		for(csym=0;csym<D;csym++) // vectorized
			evm[csym]*=(1.0f-evm[csym]*evm[csym]*r24);
		for(csym=0;csym<D;csym++) // vectorized
			evmt+=evm[csym]*evm[csym];
		evmt*=r96;  // D=96
		pbchEVM=100.0f*sqrtf(evmt);
		//	for(csym=0;csym<D;csym++)
		//		printf("%3i %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %3i %3i\n",csym,cabsf(eqData[csym]),cargf(eqData[csym]),crealf(d[csym]),cimagf(d[csym]),cargf(d[csym]),evm[csym],ee[csym*2],ee[csym*2+1]);
		//	printf("PBCH EVM=%.2f%% %i %i\n",100*sqrt(evmt),symbCnt,D);
		SXORvec(ee,cprsPBCH,E); 
		ViterbiRateDec( ee,E,pbcha,40,rv); // 40=>2 rows of 32, E>=3*32*rows, rv has no effect
		for( ci=0; ci<16; ci++ ) // vectorizes - restore parity bits
			pbcha[24+ci]^=pbchCRCmskAnt[Pant[cp]][ci];
		if( gCRC16( pbcha,40 )==0 ) // check CRC16
		{
			P=Pant[cp];
			MIBread();
			sprintHex( pbcha,40,msg );
			return(1); // CRC OK
		}
	}
	return(-1); // CRC NOT OK
}

float LTEscan::RefAvPh( float _Complex *ptr )
{
	float refPhi=0.0;
	float rPhi[12];
	unsigned char cf;
	for(cf=0;cf<12;cf++)
		rPhi[cf]=cargf(ptr[2*cf])-cargf(pilotVal[cf]);
	for(cf=0;cf<12;cf++)
		rPhi[cf]+=PI2*((rPhi[cf]<=-PI2)-(rPhi[cf]>=PI2));
	for(cf=0;cf<12;cf++)
		rPhi[cf]+=PI2*((rPhi[cf]<-PI)-(rPhi[cf]>PI));
	for(cf=0;cf<12;cf++)
		refPhi+=rPhi[cf];
	refPhi*=r12;
	// 2nd pass on noisy data at boundaries
	for(cf=0;cf<12;cf++)
		rPhi[cf]+=PI2*(((rPhi[cf]-refPhi)<-PI)-((rPhi[cf]-refPhi)>PI));
	refPhi=0.0;
	for(cf=0;cf<12;cf++)
		refPhi+=rPhi[cf];
	refPhi*=r12;
	return(refPhi);
}

float LTEscan::RefAvMg( float _Complex *ptr )
{
	float refMag=0.0;
	unsigned char cf;
	for(cf=0;cf<12;cf++)
		refMag+=crealf(ptr[2*cf])*crealf(ptr[2*cf])+cimagf(ptr[2*cf])*cimagf(ptr[2*cf]);
	refMag*=r12;
	return(sqrtf(refMag));
}

void LTEscan::DespreadPilots(void)
{ // if we do not want to plot raw pilots, we could just use despread*conj(pilotVal)
	unsigned char cf;
	for(cf=0;cf<12;cf++)
	{ 
		pilotsDespread[cf*2]=pilotsRaw[cf*2]*conjf(pilotVal[cf]);
		pilotsDespread[cf*2+1]=pilotsRaw[cf*2+1]*conjf(pilotVal[cf]);
		pilotsDespread[cf*2+24]=pilotsRaw[cf*2+24]*conjf(pilotVal[cf]);
		pilotsDespread[cf*2+25]=pilotsRaw[cf*2+25]*conjf(pilotVal[cf]);
	}
}

void LTEscan::MIBread( void )
{ // pbcha[40]
	unsigned short ci=0;
	unsigned short cnt=0;
	char flagsTxt[12]="__________\0"; // ignore
	bw=bin2dec(pbcha,cnt,3);
	cnt+=3;
	phichDur=pbcha[cnt++];
	phichRes=bin2dec(pbcha,cnt,2);
	cnt+=2;
	sfn=(short)bin2dec(pbcha,cnt,8)*4; // 256x40ms => repeat approx 10s
	cnt+=8;
//	sprintf(sfnTxt," SFN=%i",sfn);
	for(ci=0;ci<10;++ci)
		if( pbcha[cnt+ci]==1 )
			flagsTxt[ci]='*';
	cnt+=10;
}

void LTEscan::SXORvec( signed char wtbpsk[], unsigned char bin[], unsigned int n )
{ // turbofec uses weighted signed numbers {-127:127}, descramble bpsk with binary scramble code.
	unsigned int ci;
	for(ci=0;ci<n;++ci) // vectorizes
		wtbpsk[ci]*=(1-2*bin[ci]); // SIMD friendly - convert bin to -BPSK and multiply to do XOR.
}

unsigned short LTEscan::gCRC16( unsigned char a[],int aLen )
{ // used to genereate or test crcs // TS 136.212 13.0.0 Section 5.1.1
	unsigned short p[2]={0,0}; // [0] curCRC, [1] feedforward, SIMD optimized
	int ck=0;
	for( ck=0; ck<aLen; ck++ ) // 
	{ // does not vectorize
		p[1]=a[ck]^(p[0]&1);
		p[0]>>=1;
		p[0]^=(0x0408*p[1]); // note polynomial reverse (0+0+0+0)+(0,D^5,0,0)+(0+0+0+0)+(D^12+0+0+0)
		p[0]|=(1<<15)*p[1]; // p[15]=D^0 (1<<15)=0x8000
	}
	return(p[0]);
}

unsigned int LTEscan::bin2dec( unsigned char a[], unsigned short off, unsigned char bits )
{
//	printf("bin2dec(a[],%i,%i\n",off,bits);
	unsigned int num=0;
	for( unsigned char ci=0; ci<bits; ++ci )
		num+=a[off+(bits-1)-ci]*(1<<ci); // msb first
	return( num );
}

unsigned int LTEscan::bin2decRvs( unsigned char a[], unsigned short off, unsigned char bits )
{
//	printf("bin2dec(a[],%i,%i\n",off,bits);
	unsigned int num=0;
	for( unsigned char ci=0; ci<bits; ++ci )
		num+=a[off+ci]*(1<<ci); // lsb first
	return( num );
}

void LTEscan::sprintHex( unsigned char data[], unsigned int nBits,char *msg ) // diagnostic only
{
	unsigned int nBytes=(nBits>>3)+((nBits&7)>0); // was floor, change to ceiling
	unsigned char myByte=0;
	msg[0]='\0';
	for( unsigned short cb=0; cb<nBytes; ++cb )
	{
		myByte=0;
		for( unsigned char ci=0; ci<8; ++ci)
		{
			myByte<<=1;
			if( (data[ci+8*cb]!=0) && ((ci+8*cb)<nBits) )
				myByte+=1;
		}
		sprintf(msg,"%s%02X",msg,myByte);
	}
}


