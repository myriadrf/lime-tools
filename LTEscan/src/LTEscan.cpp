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
#ifdef __cplusplus
extern "C" {
#endif
#include "turbofec/conv.h"
#include "turbofec/rate_match.h"
#include "turbofecLib.h"
#ifdef __cplusplus // If this is a C++ compiler, end C linkage
}
#endif
#include "LTEscan.hpp"

LTEscan::LTEscan( short nfft )
{ // create run time storage
	NFFT=nfft;
	NFFTm37=NFFT-37;
	NFFTd2=NFFT>>1;
	NFFTd128=NFFT>>7;
	NIQ=150*NFFT*SRCH_FRMS; // total number of samples to be collected
	short ci=0;
	short cj=0;
	float PI2dNFFTm1=float(PI2)/(NFFT-1);
	float _Complex fac;
	// intended to give a very approximate default power per MHz based on band 20.  Depends on matching network and band.
	#ifdef USE_I16
	for(ci=0;ci<256;ci++) // I16=-90, I12=-66.2, F32=0
		bnd2pwr[ci]=-167.0f; // depends on F32,I16,I12
	#endif
	#ifdef USE_I12
	for(ci=0;ci<256;ci++) // I16=-90, I12=-66.2, F32=0
		bnd2pwr[ci]=-167.0f+24.1f; // depends on F32,I16,I12
	#endif
	#ifdef USE_F32
	for(ci=0;ci<256;ci++) // I16=-90, I12=-66.2, F32=0
		bnd2pwr[ci]=-167.0f+90.3f; // depends on F32,I16,I12
	#endif
	rNFFT=1.0f/NFFT; // division slow, precompute common scale factors
	rPI2=1.0f/float(PI2);
	rPI=1.0f/float(PI);
	r62=1.0f/62; // 62 s/c used for sync
	r10=0.1f; // 72-62=10
	r24=1.0f/24; // fast sine calc
	r96=1.0f/96;
	r15=1.0f/15; // Ferr Correction
	r12=1.0f/12; // Ref averaging
	rsqrt2=1.0f/sqrtf(2.0); // ref scaling
	PI2dNFFT=float(PI2)*rNFFT;
	NFFTdPI2=rPI2*NFFT;
	cp0mcp1dNFFT_2pi[0]=0.0f;
	cp0mcp1dNFFT_2pi[1]=float(PI/64); //=2*pi/128
	cp0pNFFTdcp1pNFFT[0]=1.0f;
	cp0pNFFTdcp1pNFFT[1]=(128.0f+10.0f)/(128.0+9.0f);
	//printf("NIQ=%i SRCH_FRMS=%i NFFT=%i\n",NIQ,SRCH_FRMS,NFFT);
#ifdef USE_F32
	iqData=(float _Complex *) malloc(sizeof(float _Complex)*NIQ); // I0+jQ0,I1+jQ1,I2+jQ2,...
#else
	iqData=(int16_t *) malloc(sizeof(int16_t)*NIQ*2);  // I0,Q0,I1,Q2,I2,Q2,...
#endif
#ifdef USE_C99_COMPLEX // fftw_malloc for double, fftwf_malloc for float
	in=(double _Complex *) fftw_malloc(sizeof(double _Complex)*NFFT); // align for SIMD
	out=(double _Complex *) fftw_malloc(sizeof(double _Complex)*NFFT); // align for SIMD
	printf("FFTW with C99\n");
#else // #ifdef USE_FFTW_COMPLEX
	in=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*NFFT); // align for SIMD
	out=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*NFFT); // align for SIMD
	printf("FFTW native\n");
#endif
	pfft=fftw_plan_dft_1d(NFFT,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	pifft=fftw_plan_dft_1d(NFFT,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
//	pfft=fftw_plan_dft_1d(NFFT,in,out,FFTW_FORWARD,FFTW_MEASURE); // modest time for 2048
//	pifft=fftw_plan_dft_1d(NFFT,in,out,FFTW_BACKWARD,FFTW_MEASURE); // useful speed improvment
//	pfft=fftw_plan_dft_1d(NFFT,in,out,FFTW_FORWARD,FFTW_PATIENT); // takes too long
//	pifft=fftw_plan_dft_1d(NFFT,in,out,FFTW_BACKWARD,FFTW_PATIENT); // but results not much faster
	fftacc=(float _Complex*)malloc(sizeof(float _Complex)*NFFT); // buffer1 
	ffttmp=(float _Complex*)malloc(sizeof(float _Complex)*NFFT); // help buffer for FFTW mode
	ffttmpMg=(float *)malloc(sizeof(float)*NFFT);
	ffttmpPh=(float *)malloc(sizeof(float)*NFFT);
	for(ci=0;ci<3;ci++)
	{
		fshift[ci]=(float _Complex *) malloc(sizeof(float _Complex)*NFFT);
		fac=I*(1-ci)*0.333333f*PI2dNFFT;
		for(cj=0;cj<NFFT;cj++) // +/-5kHz offsets for aligning to 100kHz earfcn gridding
			fshift[ci][cj]=cexpf(cj*fac); // apply before fft for freq shift
	}
	for(ci=0;ci<3;ci++)
		dictPSCHzdc[ci]=(float _Complex *) malloc(sizeof(float _Complex)*64);
	for(ci=0;ci<3;ci++)
		dictPSCHzdcNFFT[ci]=(float _Complex *) malloc(sizeof(float _Complex)*NFFT);
	for(ci=0;ci<3;ci++)
		dictPSCHconjFFT[ci]=(float _Complex *) malloc(sizeof(float _Complex)*NFFT);
	for(ci=0; ci<(3*168); ci++ ) // 3*168 as SSCH depends of PSCH
	{
		dictSSCH0zdc[ci]=(float _Complex *)malloc(sizeof(float _Complex)*64);
		dictSSCH10zdc[ci]=(float _Complex *)malloc(sizeof(float _Complex)*64);
	}
	for(ci=0;ci<3;ci++)
		dictPSCHzdcPh[ci]=(float *) malloc(sizeof(float)*64);
	for(ci=0; ci<(3*168); ci++ ) // 3*168 as SSCH depends of PSCH
	{
		dictSSCH0zdcPh[ci]=(float *)malloc(sizeof(float)*64);
		dictSSCH10zdcPh[ci]=(float *)malloc(sizeof(float)*64);
	}
	LoadDictPSCH(); // 180us
	LoadDictSSCH(); // 3073us
	ViterbiInit();
	// Section of ref sequence in PBCH, note NDLRBmax-NDLRB={104,95,85,60,35,10}, ref seq independent of NDLRB
	pilotStart[0]=0+94;
	pilotStart[1]=9+85;
	pilotStart[2]=19+75;
	pilotStart[3]=44+50;
	pilotStart[4]=69+25;
	pilotStart[5]=94+0;
	for( ci=0; ci<16; ci++ ) // vectorizes
	{
		pbchCRCmskAnt[0][ci]=0; // {0,0,0,...} // not used
		pbchCRCmskAnt[1][ci]=0; // {0,0,0,...}
		pbchCRCmskAnt[2][ci]=1; // {1,1,1,...}
		pbchCRCmskAnt[3][ci]=0; // {0,0,0,...} // not used
		pbchCRCmskAnt[4][ci]=ci&1; // {0,1,0,1,...} replce %2 with &1
	}
	cp0[0]=32*NFFTd128; // extended cyclic prefix length
	cp1[0]=32*NFFTd128;
	cp0[1]=10*NFFTd128; // normal cyclic prefix length
	cp1[1]= 9*NFFTd128; // symbol0 longer than symbol1:symbol6
	FrmSize=150*NFFT;
	FrmSizeD2=FrmSize>>1;
	qpskPh[0]=float( PId4);
	qpskPh[1]=float(-PId4);
	qpskPh[2]=float(PId2+PId4);
	qpskPh[3]=float(-(PId2+PId4));

	for(ci=0;ci<63;ci++) // vectorized
		sumWt[ci]=1;
	sumWt[31]=0;
	for(ci=0;ci<63;ci++)
		mWt[ci]=ci-31;
	for(ci=0;ci<63;ci++)
		if(ci!=31)
			rmWt[ci]=1.0f/mWt[ci];
		else
			rmWt[ci]=0.0f;
	for(ci=0;ci<4;ci++)
		rotWt[ci]=(float*)malloc(sizeof(float)*63);
	//printf("NFFT=%i NFFTd2=%i NFFTd128=%i FrmSize=%i FrmSizeD2=%i\n", NFFT,NFFTd2,NFFTd128,FrmSize,FrmSizeD2);
	//printf("cp0[0]=%i cp1[0]=%i cp0[1]=%i cp1[1]=%i\n",cp0[0],cp1[0],cp0[1],cp1[1]);
}

LTEscan::~LTEscan( )
{ // delete run time storage
	short ci=0;
	for(ci=0;ci<4;ci++)
		free(rotWt[ci]);
	ViterbiFree();
	for(ci=0;ci<(3*168);ci++) // x3 as SSCH depends of PSCH
	{
		free(dictSSCH0zdcPh[ci]);
		free(dictSSCH10zdcPh[ci]);
	}
	for(ci=0;ci<3;ci++)
		free(dictPSCHzdcPh[ci]);
	for(ci=0;ci<(3*168);ci++) // x3 as SSCH depends of PSCH
	{
		free(dictSSCH0zdc[ci]);
		free(dictSSCH10zdc[ci]);
	}
	for(ci=0;ci<3;ci++)
		free(dictPSCHconjFFT[ci]);
	for(ci=0;ci<3;ci++)
		free(dictPSCHzdcNFFT[ci]);
	for(ci=0;ci<3;ci++)
		free(dictPSCHzdc[ci]);
	for(ci=0;ci<3;ci++)
		free(fshift[ci]);
	free(ffttmpPh);
	free(ffttmpMg);
	free(ffttmp);
	free(fftacc);
	fftw_destroy_plan(pifft);
	fftw_destroy_plan(pfft);
	fftw_free(out);
	fftw_free(in);
	free(iqData);
}

#ifdef USE_F32
	float _Complex *LTEscan::ShareIQptr( void )
	{
		return(iqData);
	}
#else
	int16_t *LTEscan::ShareIQptr( void )
	{
		return(iqData);
	}
#endif

// We force FFT center (NFFT/2) to be on a 300kHz grid,  Lowest frequency also on 300kHz grid is (NFFT/2)%20
// We have the additional constraint to be at least 36 subcarriers away FFT edge to prevent memory violation
// we will also shift freq -5kHz or +5kHz to align with +100kHz -100kHz with 15kHz sub carrier spacing
// {-1,0,1} * 5kHz e.g. 0=0kHz (7x15-5)=100kHz, (13*15+5)=200kHz
// => -5kHz 43 (256%20-13+20)+20 = ((256-13)%20+20)+20, ((256-13)%20)=3
// =>  0kHz 36 (256%20-0 +20)    = ((256-0 )%20+20)+20, (256%20)=16
// =>  5kHz 49 (256%20+13+20)    = ((256+13)%20+20)+20, ((256+13)%20)=9
// The upper limit is given by
// => -5kHz 475 (256-13)%20+(512-40)
// =>  0kHz 468 (256-0 )%20+(512-40)-20
// =>  5kHz 461 (256+13)%20+(512-40)-20
void LTEscan::CalcLowHighFreqs( short &fL,short &fH,signed char mode )
{
//	short ffl=((NFFT/2+mode*13)%20); // lowest freq on 300kHz grid from DC+/-5kHz (NFFT/2) [300kHz = 20*15kHz]
//	fL=((36+ffl)/20+(((36+ffl)%20)>0))*20+ffl; // lowest freq on 300kHz grid from DC (NFFT/2), meeting 36s/c keepout
//	fH=20*((NFFT-37-ffl)/20)+ffl; // highest freq on 300kHz grid from DC (NFFT/2), meeting 37s/c keepout
//	fL+=mode*13;
//	fH+=mode*13;
	//printf("fL=%i fH=%i ffl=%i NFFTd300k=%i\n",fL,fH,ffl,NFFTd300k );
	fL=((NFFTd2+mode*13)%20)+20;
	fH=((NFFTd2+mode*13)%20)+(NFFT-40);
	fL+=20*(fL<36);
	fH-=20*(fH>(NFFT-37));
}

void LTEscan::ScanIqData( const char *flog,unsigned char Bnd,double frq,unsigned char subSwpNo,short cRpt )
{
//	unsigned int HalfFrmLen=(75*NFFT+1); // FrameLength=150*NFFT, need to cover just over half a frame
	unsigned int NFFTdKBLK=NFFT/KBLK; // NFFTdKSUB=NFFT/KSUB, step size in sample points
	unsigned int cBlk=0; 
	unsigned int cOff=2*KBLK; // i.e. 2*FFT margin.  SSCH NFFT+CP before PSCH
	unsigned int NBlk=75*KBLK+cOff+1; // Number of search steps per half frame
	double tscan=0.0;
	double tcheck=0.0;
	double tident=0.0;
	double tproc=0.0;
	short fL;
	short fH;
	signed short cSh=1;
	char charBuf[128];
	for(cSh=-1;cSh<2;cSh++) // -1,0,1 for optional +/-5kHz
	{
		clock_t mytref = clock();
		CalcLowHighFreqs( fL,fH,cSh );
		//printf("%i %i %i\n",cSh,fL,fH);
		LogMsg( flog, "Read samples, processing...\n" );
		tmpCellLog.CellListInit( Bnd,frq,NFFT,bnd2pwr[Bnd],cSh,subSwpNo,cRpt );
		accCellLog.CellListInit( Bnd,frq,NFFT,bnd2pwr[Bnd],cSh,subSwpNo,cRpt );
		// define offset to ensure SSCH is never in negative time! else difficulty with test routines!
		for( cBlk=cOff; cBlk<(NBlk+cOff); ++cBlk )
			EVMpsch( (cBlk*NFFTdKBLK),cSh,fL,fH ); // 
		tscan=((double)(clock()-mytref))/CLOCKS_PER_mSEC;
		//tmpCellLog.PrintCellList();
		tmpCellLog.UpdateFreqEarfcn(); // working
		tmpCellLog.SortByFreq(); // working
		tmpCellLog.AllocTmpCid(); // working
		tmpCellLog.SortByTime();
		tmpCellLog.FprintCellListSh(flog);
		tmpCellLog.PrintCellListQwkLk();
		//printf("Time Freq optimisation\n");
		tmpCellLog.LoopStart();
		while( !tmpCellLog.LoopEnded() ) // for each element in list in freq group
		{
			CheckCellPSCHevm();
			tmpCellLog.LoopNextItemp();
		}
		tcheck=((double)(clock()-mytref))/CLOCKS_PER_mSEC;
		tmpCellLog.UpdateFreqEarfcn(); // working
		tmpCellLog.FprintCellListSh(flog);
		tmpCellLog.KeepLowestEVM(); // working. discard useless results - what criterion EVM?
		tmpCellLog.PrintCellListQwkLk();
		//LogMsg( flog, "SSCH search...\n" );
		tmpCellLog.LoopStart();
		while( !tmpCellLog.LoopEnded() ) // for each element in list in freq group
		{ // TDD/FDD and NCP
			tmpCellLog.LoopAllocSparseCor();
			tmpCellLog.LoopCalcSparseCor(); // for +/- 5kHz cases (No Ferr)
			CellSSCH();
			BulkPSCHSSCH(); // update phErr for all PSCH (+/-5kHz)
			tmpCellLog.LoopNextItemp();
		}
		tident=((double)(clock()-mytref))/CLOCKS_PER_mSEC;
		tmpCellLog.FprintCellList2(flog);
		LogMsg( flog, "Purge Bad EVM\n" );
		tmpCellLog.PurgeBadEVM(); // working
		//tmpCellLog.UpdateFreqEarfcn();
		//tmpCellLog.PurgeSimilarEARFCN();
		tmpCellLog.FprintCellList2(flog);
		tmpCellLog.LoopStart();
		while( !tmpCellLog.LoopEnded() ) // for each element in list in freq group
		{		
			tmpCellLog.LoopCalcFreqErr2(flog); 
			LogMsg( flog, "LoopCalcSparseCor+BULKPSCHSSCH\n" );
			tmpCellLog.LoopCalcSparseCor();
			BulkPSCHSSCH(); // update phErr for all PSCH (+/-5kHz)
			// do a 2nd time to remove residual freq error - not much improvment
			tmpCellLog.LoopCalcFreqErr2(flog);
			LogMsg( flog, "LoopCalcSparseCor+BULKPSCHSSCH\n" );
			tmpCellLog.LoopCalcSparseCor();
			BulkPSCHSSCH(); // update phErr for all PSCH (+/-5kHz) 
			tmpCellLog.LoopCalcdpErrdtAv(flog);
			LogMsg( flog, "FindPBCH+LoopClampPlots+LoopPlotConstAll\n" );
			FindPBCH();
			tmpCellLog.LoopClampPlots();
			#ifdef USE_GNUPLOT
				tmpCellLog.LoopPlotConstAll();
			#endif
			tmpCellLog.LoopNextItemp();
		}
		tmpCellLog.FprintCellList(flog); // keep strictly chronological
		//tmpCellLog.PrintListSummary();
		tmpCellLog.FprintListSummary(flog,1); // keep strictly chronological
	
		accCellLog.MergeLists( &tmpCellLog );
		tproc=((double)(clock()-mytref))/CLOCKS_PER_mSEC;
		sprintf(charBuf,"Frq=%6.1fMHz Scan %.1fms, Check %.1f Ident %.1f Process %.1fms cOff=%i NBlk=%i NFFTdKBLK=%i\n",frq,tscan,(tcheck-tscan),(tident-tcheck),(tproc-tident),cOff,NBlk,NFFTdKBLK);
		LogMsg( flog,(const char *)charBuf );
		printf("%s",charBuf); // <n included in string
	}
}

void LTEscan::DisplayBandSummaries( const char *flog,const char *fout,const char *fhtml,short NRpt )
{
	short CRCgood=0;
	short CRCbad=0;
	accCellLog.SortByFreq();
//	accCellLog.FprintCellList(flog);
	accCellLog.PrintListSummary();
	accCellLog.FprintListSummary(flog,1);
	accCellLog.FprintListSummary(fout,0);
	accCellLog.GenSummaryHTML(fhtml);
	CRCgood=accCellLog.CountGoodCRCs();
	CRCbad=accCellLog.CountBadCRCs();
	accCellLog.ClearList();
	printf("NRpt=%i CRCgood=%i CRCbad=%i\n",NRpt,CRCgood,CRCbad );
}

void LTEscan::LogMsg( const char *fname, const char *msg )
{
	FILE *fp;
	fp=fopen( fname, "a" );
	fprintf(fp,"%s\n",msg);
	fclose(fp);
}

void LTEscan::StartOfBandMarkerForLog(const char *flog,unsigned char Bnd,unsigned char subSwpNo,short cRpt)
{
	char msg[80];
	sprintf(msg,"B%i_SS%i_R%i_________________________________________________________\n",Bnd,subSwpNo,cRpt);
	LogMsg( flog,msg );
}

