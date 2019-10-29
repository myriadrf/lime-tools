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
class CellList
{
	public: // warning - any valid CellItem must only exist in this CellList - else memory violations
	CellItem *top;
	CellItem *cur;
	CellItem *prv;
	CellItem *nxt;
	CellItem *last;
	CellItem *tmp;
	CellItem *loop;
	CellItem *best;
	short tmpCid;
	short tmpCid2;
	double frqSampdNFFT;
	short bnd; // LTE band of FFT
	double LOfrq; // center frequency of FFT
	signed char fftFIdxOff;
	unsigned short NFFT; // default number of FFT bins
	unsigned short NFFTd2;
	unsigned short NFFTd128;
	double rNFFT;
	double r15;
	unsigned short cp0[2];
	unsigned short cp1[2];
	char dateStr[16];
	char time24Str[16]; // yr-mnth-dy,24h time
	float LTEbw[6];
	short LTEndlrb[6];
	short cRpt;
	unsigned char subSwpNo;

	CellList();
	~CellList(); // delete everything
	void ClearList( void );
	void CellListInit(short bnd,double lofrq,short nfft,float Bnd2Pwr,signed char FFTFIdxOff,unsigned char cSwp,short crpt);
	void DelCurCell( void );
	void MoveCurToTmp( void );
	void MoveTmpBeforeCur( void );
	void PrintCellListQwkLk( void );
	void PrintCellList( void );
	void FprintCellList( const char *fname );
	int SprintCellList( char *buf );
	void PrintCellList2( void );
	void FprintCellList2( const char *fname );
	void PrintCellListSh( void );
	void FprintCellListSh( const char *fname );
	void PrintListSummary( void );
	void FprintListSummary( const char *fname,unsigned char dispHdr );
	void GenSummaryHTML( const char *fname );
	void FindSmallestFreq( void );
	void FindSmallestTime( short curTmpcid );
	void SortByFreq( void );
	void AllocTmpCid( void );
	void UpdateFreqEarfcn( void );
	void RenumTmpCid( void );
	void SortByTime( void ); // within tmpCid groups
	void KeepLowestEVM( void ); // within tmpCid groups
	void LoopStart( void );
	void LoopNextItemp( void );
	short LoopEnded( void );
	int LoopGetEARFCN( void );
	short LoopGetFftFIdx( void );
	signed char LoopGetFftFIdxOff( void );
	unsigned short LoopGetNID( void );
	unsigned char LoopGetN1ID( void );
	unsigned char LoopGetN2ID( void );
	unsigned char LoopGetNCP( void );
	unsigned char LoopGetNs( unsigned char idx );
	float LoopGetTdPSCH( unsigned char idx );
	float LoopGetTdSSCH( unsigned char idx );
	float LoopGetPhErrPSCH( unsigned char idx );
	float LoopGetPhErrSSCH( unsigned char idx );
	float LoopGetRrmsPSCH( unsigned char idx );
	float LoopGetRrmsSSCH( unsigned char idx );
	float LoopGetEvmPSCH( unsigned char idx );
	float LoopGetEvmSSCH( unsigned char idx );
	int LoopGetTimePSCH( unsigned char idx );
	int LoopGetTimeSSCH( unsigned char idx );
	int LoopGetTref( void );
	float _Complex *LoopGetPSCHptr( void );
	float _Complex *LoopGetSSCHptr( void );
	float _Complex *LoopGetPBCHptr( void );
	float _Complex *LoopGetPBCHpilotRawPtr( void );
	float _Complex *LoopGetPBCHpilotDeSpreadPtr( void );
	char *LoopGetPBCHmsgPtr( void );
	float _Complex **LoopGetSparseCor( void );
	float bnd2pwr; // convert digital power (63s/c) dB into dBm + 10*log10(12*NDLRBs/62)
	signed char LoopGetPBCHcrc( void );
	//void LoopResetEVM( void );
	void LoopGetEqv( float _Complex *eqv );
	float LoopGetFerr( unsigned char idx );
	void LoopSetFerr( float fErr,unsigned char idx );
	void LoopAllocSparseCor( void );
//	void AddLastCellPSCH( int time,unsigned char c2id,short fftF,float *p,signed char fftFIdxOff );
	void AddLastCellPSCH( int time,unsigned char c2id,short fftF,float *p );
	void LoopUpdatePSCH(int time,short cf,unsigned char c2id,unsigned char evmIdx,float *p );
	void LoopUpdateSSCH(int time,unsigned char c1id,unsigned char evmIdx,float *p,unsigned char ns,unsigned char ncp );
	//void LoopAddPBCH(float _Complex *data,short dataLen,unsigned char pbchIdx);
	void LoopClampPlots( void );
	void LoopCalcFreqErr2( const char *fname );
	void LoopCalcdpErrdtAv( const char *fname );
	void LoopCalcSparseCor( void );
	double LoopGetdpErrdtAv( void );
	void UpdatePBCH(float pbchEVM,unsigned char P,unsigned char BW,unsigned char PHICHdur,unsigned char PHICHres,short nSysFrm,signed char crc,unsigned char evmIdx,float errAvPBCH,short pbchSymbCnt);
	void PurgeBadEVM( void );
	void PurgeSimilarEARFCN( void );
	short CountGoodCRCs( void );
	short CountBadCRCs( void );
	
	// test nid of 2nd sample.  Then delete ones which fail test.
	// then sort by EVM/Pwr, and delete everything else
	void MergeLists( CellList *ext ); // group later results with earlier.
	void CheckList( signed char fullList );
	void LoopPlotConstAll( void );
};
