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

class CellItem
{
	public:
	CellItem *nxt;
	CellItem *prv; // need bidirectional list to allow insert, deletes and moves
	int idx; // item number when list created - for debugging gives each element a unique ID
	short lteBand; // LTE band number 1-255, only 1:85 implemented
	short tmpCid; // group of frequencies within 500kHz
	double freq; // MHz - derived from FFT index and LOfrq
	short fftFIdx; // index of fft vector
	signed char fftFIdxOff; // descrete 5kHz shift used before fft. {-1,0,1}*5kHz
	int earfcn;
	short NFFT;
	short NFFTd2;
	int timeRef; // reference time location of first OFDM symbol of first detected frame
	int timePSCH[FERR_COR_HFRMS]; // time of FFT evaluation
	float evmPSCH[FERR_COR_HFRMS];
	float evm2PSCH[FERR_COR_HFRMS]; // evm^2 - useful for composite EVM check
	float pwrPSCH[FERR_COR_HFRMS]; 
	float snrPSCH[FERR_COR_HFRMS]; 
	float tdPSCH[FERR_COR_HFRMS];
	float phErrPSCH[FERR_COR_HFRMS];
	float rrmsPSCH[FERR_COR_HFRMS];

	int timeSSCH[FERR_COR_HFRMS]; // time of FFT evaluation
	float evmSSCH[FERR_COR_HFRMS];
	float pwrSSCH[FERR_COR_HFRMS]; 
	float snrSSCH[FERR_COR_HFRMS]; 
	float tdSSCH[FERR_COR_HFRMS];
	float phErrSSCH[FERR_COR_HFRMS];
	float rrmsSSCH[FERR_COR_HFRMS];
	float bnd2pwr;

	double dphErrdt[FERR_COR_HFRMS];
	double dphErrdtSD[FERR_COR_HFRMS];
	double dphErrdtAv;	// note dphdt is x2 FFTs sepparated by NFFT+cp1[NCP] points
	double dphErrdtSDav;

	double LOfrq;
	float fErr[2]; // intended for small doppler shifts <0.25 1/NFFT
	float _Complex **sparseCor; // sparseCor[FERR_COR][NFFT] - sparse matrix of frequency corrections for FFT.
	double rNFFT;
	float evmPBCH;
	float errAvPBCH; // average phase angle in PBCH symbols - debug tool to solve PBCH EVM issues
	short pbchLen;
	char pbchMsg[16]; // in hex
	signed char pbchCRC;
	unsigned char pbchIdx; // location of PBCH from start of grab in half frames
	unsigned char PhichDur;
	unsigned char PhichRes;
	

	unsigned char n2id; // psch 0:2
	unsigned char n1id; // ssch 0:127
	unsigned short nid; // 3*n1id+n2id = 0:511
	unsigned char ncp; // 0 extended, 1 normal cyclic prefix
	unsigned char ns[FERR_COR_HFRMS];
	short nSysFrm;
	unsigned char P; // number of antennas {1,2,4}
	float BW;
	short NDLRB;
	double rPI2;
	double r15;

	unsigned short cp0[2];
	unsigned short cp1[2];
	float cp0dNFFT_2pi[2]; // NCP=0 =32/128*2*pi=0.5*Pi, NCP=1 =10/128*2*pi=10/64*pi
	float cp1dNFFT_2pi[2]; // NCP=0 =32/128*2*pi=0.5*Pi, NCP=1 =9/128*2*pi=9/64*pi
	float poff;
	// some sort of time date stamp for measurements

	short pbchSymbCnt; // usually 252 for NCP=0, 264 for NCP=1
	float _Complex psch[64];
	float _Complex ssch[64];
	float _Complex pbch[292]; // 4x73, including pilots and DC
	float _Complex pilotsRaw[48]; // 0,1 interleaved 0:23, 2,3 interleaved 24:47. ignore extra 0,1 when NCP=0
	float _Complex pilotsDeSpread[48]; // 0,1 interleaved 0:23, 2,3 interleaved 24:47. ignore extra 0,1 when NCP=0
	
	char dateStr[16];
	char time24Str[16]; // yr-mnth-dy,24h time
	short cRpt;
	unsigned char subSwpNo;

	float lms7002Temp;
	char cpuTemp[256];

	CellItem();
	~CellItem();
	void PrintQwkLk( void );
	void Print( void );
	void Fprint( FILE *fp );
	int Sprint( char *str );
	void Print2( void );
	void Fprint2( FILE *fp );
	void PrintSh( void );
	void FprintSh( FILE *fp );
	void PrintSum( void );
	void FprintSum( FILE *fp );
	void GenSummaryHTML( FILE *fp );
	void ClampPSCH( void );
	void ClampSSCH( void );
	void ClampPBCH( void );
	void FprintPhVsT( FILE *fp,unsigned char idx );
	void CalcdpErrdtAv( void );
	void CalcSparseCor( void );
};
