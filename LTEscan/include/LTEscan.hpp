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

class LTEscan
{	
	short NFFT;
	short NFFTm37;
	short NFFTd128;
	int NIQ;
	float r62; // 62 s/c used for PSCH/SSCH
	float r10; // 72-62=10
	float rNFFT; // division slow, precompute common reciprocal scale factors
	float rPI2;
	float rPI;
	float r12; // PBCH Ref averaging
	float r15; // Ferr Correction
	float r24; // fast sine calc
	float r96; // PBCH evm
	float rrms;
	float rsqrt2; // refs
	float _Complex *dictPSCHzdc[3];
	float _Complex *dictPSCHzdcNFFT[3];
	float _Complex *dictPSCHconjFFT[3];
	float _Complex *dictSSCH0zdc[(3*168)];
	float _Complex *dictSSCH10zdc[(3*168)];
	float *dictPSCHzdcPh[3]; 
	float *dictSSCH0zdcPh[(3*168)];
	float  *dictSSCH10zdcPh[(3*168)];
	float bnd2pwr[256]; // convert measured digital power back to dBm/MHz (approximately)
	float PI2dNFFT;
	float NFFTdPI2;
#ifdef USE_F32
	float _Complex *iqData=NULL;
#else
	int16_t *iqData=NULL;
#endif
	fftw_plan pfft; // FFTW
	fftw_plan pifft; // FFTW
	float _Complex *fshift[3]; // -5,0,5kHz offset vectors
#ifdef USE_C99_COMPLEX
	double _Complex *in=NULL; // registered with FFTW
	double _Complex *out=NULL; // registered with FFTW
#else // #ifdef USE_FFTW_COMPLEX
	fftw_complex *in=NULL;
	fftw_complex *out=NULL;
#endif
//	float _Complex *fftbuf; // private storage
	float _Complex *fftacc; // private storage
	float _Complex *ffttmp; // Helper for FFTW mode
	float *ffttmpMg; // Helper for EVM calcs
	float *ffttmpPh; // Helper for EVM calcs
	short pilotStart[6]; // start location of PBCH for different BWs=0:5
	short pbchSymbCnt; // usually 252 for NCP=0, 264 for NCP=1
	float _Complex d[240]; // [2/3*72*2+2*72] combined MIMO signal
	float _Complex d1[240]; // raw input signal
	float _Complex d2[240]; // transmit diversity version of d1[]
	float _Complex pilotVal[12]; // calculated from PRS
	float _Complex pilotsRaw[64]; // measured, 0:23 0,1 interleaved, 24:47 2,3 interleaved,
	float _Complex pilotsDespread[48]; // measured, 0:23 0,1 interleaved, 24:47 2,3 interleaved, repeat of 0,1 for NCP=0 ignored.
	signed char ee[480]; // encoded bits
	unsigned char cprsPBCH[192];
	float evm[192];
	float errAvPBCH;
	unsigned char pbchCRCmskAnt[5][16]; // indexes {1,2,4}
	unsigned char P;
	unsigned char bw;
	unsigned char phichDur;
	unsigned char phichRes;
	short sfn;
	float pbchEVM;
	short cp0[2]; // NCP=0 => 32*NFFT/128, NCP=1 => 10*NFFT/128
	short cp1[2]; // NCP=0 => 32*NFFT/128, NCP=1 =>  9*NFFT/128
	float cp0mcp1dNFFT_2pi[2]; // NCP=0 =0, NCP=1 =1/128*2*pi=pi/64
	float cp0pNFFTdcp1pNFFT[2]; // NCP=0 =1, NCP=1 =(128+10)/(128+9)
	float qpskPh[4];
	CellList tmpCellLog;
	CellList accCellLog;
	unsigned short NFFTd2;
	unsigned short NFFTd300k;
	unsigned int NFFTdKBLK;
	unsigned int FrmSize;
	unsigned int FrmSizeD2;
	#define PTS_PRS 5 // pairs of points 4-->7?  Should get 8dB improvment moving from 2pt to 14pt
	#define r2PTS_PRS (1.0/(2*PTS_PRS))
	float *rotWt[4]; // MakeEqvVecFromXSCH3
	signed char sumWt[63]; // MakeEqvVecFromXSCH3
	float rmWt[63]; // MakeEqvVecFromXSCH3
	signed char mWt[63]; // MakeEqvVecFromXSCH3
	public:
	LTEscan( short nfft );
	~LTEscan();
	void ScanIqData( const char *flog,unsigned char Bnd,double frq,unsigned char subSwpNo,short cRpt );
	void DisplayBandSummaries( const char *flog,const char *fout,const char *fhtml,short NRpt );
	void StartOfBandMarkerForLog(const char *flog,unsigned char Bnd,unsigned char subSwpNo,short cRpt);
	void LogMsg( const char *fname, const char *msg );

#ifdef USE_F32
	float _Complex *ShareIQptr( void );
#else
	int16_t *ShareIQptr( void );
#endif	
	private:
	void CalcRMS62( short foff,float *p ); // RMS power in 62, from vec
	void CalcRMS10( short foff,float *p ); // RMS noise in 10, from vec
	void MakeEqvVecFromXSCH1(short cid,short foff,float _Complex **dict,float **dictPh,float *p,float *sd); // use two s/c to correct phase of remaining s/c
	void MakeEqvVecFromXSCH1A(short cid,short foff,float _Complex **dict,float **dictPh,float *p,float *sd); // use 4th s/c
	void MakeEqvVecFromXSCH1B(short cid,short foff,float _Complex **dict,float **dictPh,float *p,float *sd); // use 6th s/c
	void MakeEqvVecFromXSCH2(short cid,short foff,float _Complex **dict,float **dictPh,float *p,float *sd); // use all s/c to correct phase, assuming linear phase with offset
	void makeRotWt( float *ptr );
	void clampPh( float *ptr,unsigned char idx );
	void fit( float *ptr,float *prm,unsigned char idxs );
	void MakeEqvVecFromXSCH3(short cid,short foff,float _Complex **dict,float **dictPh,float *p,float *sd);
	void makeRotWtF( float *ptr );
	void clampPhF( float *ptr,unsigned char idx );
	void fitF( float *ptr,float *prm,unsigned char idxs );
	void MakeEqvVecFromXSCH3F(short cid,short foff,float _Complex **dict,float **dictPh,float *p,float *sd);
	
	float CalcFastEVM( float *phd,float dly,float phDop );
	void CheckEVM(float _Complex *ptr,float _Complex **dict,float **dictPh,short cid);

	void EVMpsch( int offset,signed char fftFIdxOff,short fL,short fH ); // take FFT of time block and scan for PSCH	
	void FftShift( float _Complex *buffer );
	void Blk2FFT( int offset,signed char fftFIdxOff,unsigned char sparseIdx ); // handles interface to and from FFTW
	void Blk2FFT2( int offset,signed char fftFIdxOff,unsigned char sparseIdx ); // handles interface to and from FFTW
	void CalcLowHighFreqs( short &fL,short &fH,signed char mode ); // calc min max search frequencies
	void CheckCellPSCHevmFoff( int offset,unsigned char evmIdx,short cf,unsigned char c2id );
	void CheckCellPSCHevm( void );

	void FastPSCHevm( int toff,unsigned char frm,unsigned char ns,unsigned char evmIdx );
	void FastSSCHevm( int toff,unsigned char frm,unsigned char ns,unsigned char evmIdx );
	signed char FastPBCHevm( int toff,unsigned char frm,unsigned char evmIdx );
	void FindPBCH( void );

	void CellSSCH( void );
	float FastEVM( float *ph );
	void EqualiseData( int toff,short foff,signed char fftFIdxOff,unsigned char NCP,unsigned char frm,unsigned char ns,unsigned char sym,float td,float phc,float rrms,float _Complex *data,unsigned char prtDiag );
	// PSCH + SSCH + Refs routines
	void LoadDictPSCH( void );
	void LoadDictSSCH(void);
	void LTEDLGold1SSSR( signed char *state );
	void LTEDLGold2SSSR( signed char *state );
	void LTEDLGold3SSSR( signed char *state );
	void Copy31( signed char *src, signed char *dest );
	void CycRot31( signed char *seq, signed short rot );
	signed char SXOR( signed char a, signed char b );
	void LoadRefs(unsigned char nSlot,unsigned char nSymbol,unsigned char NID,unsigned char NCP);
	// PBCH routines
	signed char DemodPBCH( float _Complex *eqData,unsigned char NCP,short NcellID,char *msg );
	float RefAvPh( float _Complex *ptr );
	float RefAvMg( float _Complex *ptr );
	void DespreadPilots(void);
	void MIBread( void );
	void LTESC( unsigned char *psc,unsigned int cinit,unsigned int cnt );
	void SXORvec( signed char wtbpsk[], unsigned char bin[], unsigned int n );
	unsigned short gCRC16( unsigned char a[],int aLen );
	unsigned int bin2dec( unsigned char a[], unsigned short off, unsigned char bits );
	unsigned int bin2decRvs( unsigned char a[], unsigned short off, unsigned char bits );
	void sprintHex( unsigned char data[], unsigned int nBits,char *msg );
	// test routines
	float CalcSlotBndryPhShft( float ph,short foff,unsigned char NCP );
	void CalcSparseCor( void );
	void BulkPSCHSSCH( void );
	unsigned char CalcSparseIdx( unsigned char relFrm,unsigned char ns,unsigned char sym,unsigned char NCP );
};
