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
#include <complex.h> // double _Complex
#include <errno.h> // ERR
#include <unistd.h> // pipes usleep
#include <fcntl.h> // file control
#include <sys/types.h> // pipes
#include <time.h> // tools for telling the time!
#include "lime.h"
#include "CellItem.hpp"
#include "CellList.hpp"
#ifdef MXD_CMPLR // needed if C compiled with gcc/clang, and cpp compiled g++/clang++
extern "C" {
#endif
#include "earfcn.h"  // .c file
#ifdef MXD_CMPLR
}
#endif

extern float temp7002;
extern char cpuTemp[256];

CellList::CellList()
{
	float bw[6]={1.4f,3.0f,5.0f,10.0f,15.0f,20.0f};
	short ndlrb[6]={6,15,25,50,75,100};
	unsigned char ci=0;
	top=NULL;
	cur=NULL;
	nxt=NULL;
	last=NULL;
	tmp=NULL;
	loop=NULL;
	best=NULL;
	tmpCid=0;
	tmpCid2=0;
	bnd=0;
	for(ci=0;ci<6;ci++)
		LTEbw[ci]=bw[ci];
	for(ci=0;ci<6;ci++)
		LTEndlrb[ci]=ndlrb[ci];
	r15=1.0/15;
}

CellList::~CellList()
{
	cur=top;
	while( cur!=NULL )
	{
		nxt=cur->nxt;
		delete cur;
		cur=nxt;
	}
}

void CellList::ClearList( void )
{
	cur=top;
	while( cur!=NULL )
	{
		nxt=cur->nxt;
		delete cur;
		cur=nxt;
	}
	top=NULL;
	cur=NULL;
	nxt=NULL;
	last=NULL;
	tmp=NULL;
	loop=NULL;
	best=NULL;
}

void CellList::CellListInit(short LTEbnd,double lofrq,short nfft,float Bnd2Pwr,signed char FFTFIdxOff,unsigned char cSwp,short crpt)
{
	bnd2pwr=Bnd2Pwr;
	time_t now=time(NULL);
	struct tm *myTime=gmtime(&now);
	sprintf(dateStr,"%04i-%02i-%02i",myTime->tm_year+1900,myTime->tm_mon+1,myTime->tm_mday);
	sprintf(time24Str,"%02i:%02i:%02i",myTime->tm_hour,myTime->tm_min,myTime->tm_sec);
	frqSampdNFFT=0.015; // constant in LTE
	bnd=LTEbnd; // LTE band being scanned
	LOfrq=lofrq; // center frequency of FFT
	fftFIdxOff=FFTFIdxOff; // -1,0,1 for optional +/-5kHz
	NFFT=nfft; // default number of FFT bins
	NFFTd2=NFFT>>1;
	rNFFT=1.0/NFFT;
	NFFTd128=NFFT>>7;
	cp0[0]=32*(NFFTd128);
	cp0[1]=10*(NFFTd128);
	cp1[0]=32*(NFFTd128);
	cp1[1]=9*(NFFTd128);
	cRpt=crpt;
	subSwpNo=cSwp;
}

void CellList::DelCurCell( void )
{
//	printf("DelCurCell\n");
	if( cur!=NULL )
	{
		prv=cur->prv;
		nxt=cur->nxt;
		if( cur==top ) // if first item of list
		{
			top=cur->nxt;
			prv=NULL;
		}
		if( cur==last ) // if last item of list
		{
			last=cur->prv;
			nxt=NULL;
		}
		if( prv != NULL ) // if there is an item before
			prv->nxt=nxt;
		if( nxt != NULL ) // if there is an item after
			nxt->prv=prv;
		delete cur;
		cur=NULL;
	}
}

void CellList::MoveCurToTmp( void )
{
//	printf("MoveCurToTmp\n");
	if( cur!=NULL )
	{
		prv=cur->prv;
		nxt=cur->nxt;
		if( cur==top ) // if first item of list
		{
			top=cur->nxt;
			prv=NULL;
		}
		if( cur==last ) // if last item of list
		{
			last=cur->prv;
			nxt=NULL;
		}
		if( prv != NULL ) // if there is an item before
			prv->nxt=nxt;
		if( nxt != NULL ) // if there is an item after
			nxt->prv=prv;
	}
	tmp=cur;
	cur=NULL;
}

void CellList::MoveTmpBeforeCur( void )
{
//	printf("MoveTmpBeforeCur\n");
	if((tmp==NULL)||(cur==NULL))
		return;
	prv=cur->prv;
	nxt=cur;
	cur=tmp;
	cur->prv=prv;
	cur->nxt=nxt;
	if( prv!=NULL )
		prv->nxt=cur;
	else
		top=cur;
	if( nxt!=NULL )
		nxt->prv=cur;
	else
		last=cur;
}

void CellList::PrintCellListQwkLk( void )
{
	cur=top;
	while( cur!=NULL )
	{
		cur->PrintQwkLk();
		cur=cur->nxt;
	}
}

void CellList::PrintCellList( void )
{
	cur=top;
	while( cur!=NULL )
	{
		cur->Print();
		cur=cur->nxt;
	}
}

void CellList::FprintCellList( const char *fname )
{
	FILE *fp;
	fp=fopen( fname, "a" );
	if( fp!=NULL )
	{
		cur=top;
		while( cur!=NULL )
		{
			cur->Fprint( fp );
			cur=cur->nxt;
		}
	}
	fclose(fp);
}

int CellList::SprintCellList( char *buf )
{
	char *ptr=buf;
	int cnt;
	cur=top;
	while( cur!=NULL )
	{
		ptr+=cur->Sprint( ptr );
		cur=cur->nxt;
	}
	cnt=ptr-buf;
	return(cnt);
}

void CellList::PrintCellListSh( void )
{
	cur=top;
	printf( "      EARFCN FRQ      EVMp   EVMp POWp   POWp SNRp SNRp    TDp    TDp   Php   Php    tp    tp Fi 5k N1ID\n");
	while( cur!=NULL )
	{
		cur->PrintSh();
		cur=cur->nxt;
	}
}

void CellList::PrintCellList2( void )
{
	cur=top;
	printf( "Bd EARFCN FRQ          EVMp   EVMp   PWRp0 PWRp1 PWRs0 PWRs1 SNRp0 SNRp1 SNRs0 SNRs1  Tp0   Tp1   Ts0   Ts1 FFTi 5k NID NCP ns0 ns1\n");
	while( cur!=NULL )
	{
		cur->Print2();
		cur=cur->nxt;
	}
}

void CellList::FprintCellList2( const char *fname )
{
	FILE *fp;
	fp=fopen( fname, "a" );
	if( fp!=NULL )
	{
		cur=top;
		fprintf( fp,"Bd EARFCN FRQ          EVMp   EVMp   PWRp0 PWRp1 PWRs0 PWRs1 SNRp0 SNRp1 SNRs0 SNRs1  Tp0   Tp1   Ts0   Ts1 FFTi 5k NID NCP ns0 ns1\n");
		while( cur!=NULL )
		{
			cur->Fprint2( fp );
			cur=cur->nxt;
		}
	}
	fclose(fp);
}

void CellList::FprintCellListSh( const char *fname )
{
	FILE *fp;
	fp=fopen( fname, "a" );
	if( fp!=NULL )
	{
		fprintf( fp,"      EARFCN FRQ    EVMp   EVMp   POWp   POWp SNRp SNRp    TDp    TDp   Php   Php    tp    tp Fi 5k N1ID\n");
		cur=top;
		while( cur!=NULL )
		{
			cur->FprintSh( fp );
			cur=cur->nxt;
		}
	}
	fclose(fp);
}

void CellList::PrintListSummary( void )
{
	cur=top;
	if(cur!=NULL)
		printf("                    BND EARFCN FreqMHZ NID NCP Ant BwMHz PHICHd PHICHr CRC PWRdB SNRdB P_EVM%% S_EVM%% B_EVM%% FerrkHz fLOMHz\n");
	while( cur!=NULL )
	{
		//printf("%s %s ",dateStr,time24Str);
		cur->PrintSum();
		cur=cur->nxt;
	}
}

void CellList::FprintListSummary( const char *fname,unsigned char dispHdr )
{
	FILE *fp;
	fp=fopen( fname,"a" );
	cur=top;
	if((dispHdr>0)&&(cur!=NULL))
		fprintf(fp,"                    BND EARFCN FreqMHZ NID NCP Ant BwMHz PHICHd PHICHr CRC PWRdB SNRdB P_EVM%% S_EVM%% B_EVM%% FerrkHz fLOMHz\n");
	while( cur!=NULL )
	{
		//fprintf(fp,"%s %s ",dateStr,time24Str);
		cur->FprintSum( fp );
		cur=cur->nxt;
	}
	fclose(fp);
}

void CellList::GenSummaryHTML( const char *fname )
{
	FILE *fp;
	fp=fopen( fname, "a" );
	cur=top;
	//fprintf(fp,"<H2>%s %s Band %i</H2>",dateStr,time24Str,bnd);
	fprintf(fp,"<H2>Band %i</H2>",bnd);
	while( cur!=NULL )
	{
		cur->GenSummaryHTML( fp );
		cur=cur->nxt;
	}
	fclose(fp);
}

void CellList::FindSmallestFreq( void )
{
//	printf("CellList::FindSmallestFreq\n");
	tmp=cur; // starting from cur...
	while( (tmp!=NULL) && (cur!=NULL) )
	{
		if(cur->freq<tmp->freq) // set tmp to lowest freq so far
			tmp=cur;
		cur=cur->nxt;
	}
	cur=tmp; // cur becomes lowest freq of list
}

void CellList::FindSmallestTime( short curTmpCid )
{
//	printf("CellList::FindSmallestTime\n");
	tmp=cur;
	while( (tmp!=NULL) && (cur!=NULL) && (cur->earfcn==tmp->earfcn ) )
	{
		if(cur->timePSCH[0]<tmp->timePSCH[0])
			tmp=cur;
		cur=cur->nxt;
	}
	cur=tmp;
}

void CellList::SortByFreq( void )
{
	//printf("CellList::SortByFreq\n");
	loop=top;
	while( loop!=NULL )
	{
		cur=loop;
		FindSmallestFreq(); // make cur the smallest freq after loop pointer, else cur=loop
		if( cur!=loop ) // move if there is smaller
		{
			MoveCurToTmp();
			cur=loop;
			MoveTmpBeforeCur();
			loop=cur;
		}
		loop=loop->nxt;	
	}
}

void CellList::AllocTmpCid( void )
{
	//printf("CellList::AllocTmpCid\n");
	prv=top;
	cur=NULL;
	tmpCid=0;
	if( prv!=NULL )
		cur=prv->nxt;
	while( cur!=NULL )
	{
		if( (prv->freq+0.65) < cur->freq ) // is upper bound more than 0.5MHz from lower bound
		{
			tmpCid++;
			prv=cur; // move lower bound to next tmpCid lower bound
		}
		cur->tmpCid=tmpCid;
		cur=cur->nxt;
	}
}

void CellList::UpdateFreqEarfcn( void )
{
	//printf("CellList::UpdateFreqEarfcn\n");
	cur=top;
	while( cur!=NULL )
	{
		cur->lteBand=bnd;
		cur->freq=(LOfrq+frqSampdNFFT*(cur->fftFIdx-NFFTd2)+0.005*cur->fftFIdxOff); // was cur->fftFIdx-NFFTd2-1
		cur->earfcn=LTEbandDLFreq2EARFCN(bnd,cur->freq);
		//printf("UpdateFreqEarfcn %f %i %i\n",cur->freq,cur->earfcn,(cur->fftFIdx-NFFTd2));
		cur=cur->nxt;
	}
}

void CellList::RenumTmpCid( void )
{ // assume numbering is discontinuous but grouped, e.g. 0,1,3,5,2,3,7
	//printf("CellList::RenumTmpCid\n");
	cur=top;
	tmpCid=cur->tmpCid;
	tmpCid2=cur->tmpCid;
	while( cur!=NULL )
	{
		if( cur->tmpCid!=tmpCid2 ) // note fail case band!=, but tmpCid==
		{
			tmpCid++;
			tmpCid2=cur->tmpCid;
		}
		cur->tmpCid=tmpCid;
		cur=cur->nxt;
	}
}

void CellList::SortByTime( void ) // most results will be in groups of time with the odd stray
{
	//printf("CellList::SortByTime\n");
	while( loop!=NULL )
	{
		cur=loop;
		tmpCid=cur->tmpCid;
		cur->Print();
//		printf("%i %i %i\n",cur->tmpCid,cur->time,cur->earfcn);
		FindSmallestTime(tmpCid);
		if( (cur!=loop) && (cur!=NULL) && (loop!=NULL) ) // only if needs moving, move it
		{
			MoveCurToTmp();
			cur=loop;
			MoveTmpBeforeCur();
			if(cur!=NULL)
				loop=cur;
		}
		loop=loop->nxt;			
	}
}

void CellList::KeepLowestEVM( void )
{ // must alloc TmpCID first
	short tmpCid=0;
	//printf("CellList::KeepLowestEVM\n");
	loop=top;
	while( loop!=NULL )
	{		
		cur=loop;
		best=cur;
		tmpCid=best->tmpCid;
		while( (cur!=NULL) && (cur->tmpCid==tmpCid) ) // find lowest evm in tmpCid 
		{
			//printf("tmpCid=%i Cur evmPSCH[0]=%.3f Best evmPSCH[0]=%.3f\n",tmpCid,cur->evmPSCH[0],best->evmPSCH[0]);
			if( cur->evmPSCH[0]<best->evmPSCH[0] )
				best=cur;
			cur=cur->nxt;
		}
//		best->Print();
		cur=loop;
		while( (cur!=NULL) && (cur->tmpCid==tmpCid) )
		{
			nxt=cur->nxt;
			if(cur!=best)
				DelCurCell(); // => nxt=cur->nxt
			cur=nxt;
		}
		loop=best->nxt;
	}
}

void CellList::PurgeBadEVM( void )
{
	float evmAv=0.0f;
	//printf("CellList::PurgeBadEVM\n");
	cur=top;
	while( cur!=NULL )
	{
		nxt=cur->nxt;
		evmAv=(cur->evmPSCH[0])*(cur->evmPSCH[0])+(cur->evmPSCH[1])*(cur->evmPSCH[1]);
		evmAv+=(cur->evmSSCH[0])*(cur->evmSSCH[0])+(cur->evmSSCH[1])*(cur->evmSSCH[1]);
		evmAv*=0.25f;
//		printf("evmPSCH[0]=%.3f evmPSCH[1]=%.3f evmAv=%.3f\n",cur->evmPSCH[0],cur->evmPSCH[1],evmAv);
//		if(cur->evmPSCH[0]>70.0) // can't use 2nd, as buffer corruption could occur.
		if(evmAv>5500.0f)	// evm _in_%^2	
			DelCurCell();  // => nxt=cur->nxt
		cur=nxt;
	}
}

void CellList::PurgeSimilarEARFCN( void )
{
	//printf("CellList::PurgeSimilarEARFCN\n");
	cur=top;
	while( cur!=NULL )
	{
		nxt=cur->nxt;
		if(nxt!=NULL)
		{
			if(abs(cur->earfcn-nxt->earfcn)<5) // if earfcn is similar, keep one with best EVM
			{
				if(cur->evmPSCH[0]>nxt->evmPSCH[0])
					DelCurCell(); // => nxt=cur->nxt
				else // cur->evmdbPSCH[0]>=nxt->evmdbPSCH[0]
				{
					cur=nxt;
					DelCurCell(); // => nxt=cur->nxt
				}
			}
		}
		cur=nxt;
	}
}

void CellList::LoopStart( void )
{
	loop=top;

}

void CellList::LoopNextItemp( void )
{
	loop=loop->nxt;

}

short CellList::LoopEnded( void )
{
	return(loop==NULL);
}

int CellList::LoopGetTref( void )
{
	return(loop->timeRef);
}

int CellList::LoopGetTimePSCH( unsigned char idx )
{
	return(loop->timePSCH[idx]);
}

int CellList::LoopGetTimeSSCH( unsigned char idx )
{
	return(loop->timeSSCH[idx]);
}

int CellList::LoopGetEARFCN( void )
{
	return(loop->earfcn);
}

short CellList::LoopGetFftFIdx( void )
{
	return(loop->fftFIdx);
}

signed char CellList::LoopGetFftFIdxOff( void )
{
	return(loop->fftFIdxOff);
}

unsigned short CellList::LoopGetNID( void )
{
	return(loop->nid);
}

unsigned char CellList::LoopGetN1ID( void )
{
	return(loop->n1id);
}

unsigned char CellList::LoopGetN2ID( void )
{
	return(loop->n2id);
}

float CellList::LoopGetEvmPSCH( unsigned char idx )
{
	return(loop->evmPSCH[idx]);// EVM in %
}

float CellList::LoopGetEvmSSCH( unsigned char idx )
{
	return(loop->evmSSCH[idx]);// EVM in %
}

float CellList::LoopGetTdPSCH( unsigned char idx )
{
	return((loop->tdPSCH[idx]));// EVM in %
}

float CellList::LoopGetTdSSCH( unsigned char idx )
{
	return((loop->tdSSCH[idx]));// EVM in %
}

float CellList::LoopGetPhErrPSCH( unsigned char idx )
{
	return(loop->phErrPSCH[idx]);// EVM in %
}

float CellList::LoopGetPhErrSSCH( unsigned char idx )
{
	return(loop->phErrSSCH[idx]);// EVM in %
}

float CellList::LoopGetRrmsPSCH( unsigned char idx )
{
	return(loop->rrmsPSCH[idx]);// EVM in %
}

float CellList::LoopGetRrmsSSCH( unsigned char idx )
{
	return(loop->rrmsSSCH[idx]);// EVM in %
}

unsigned char CellList::LoopGetNCP( void )
{
	return(loop->ncp);// EVM in %
}

unsigned char CellList::LoopGetNs( unsigned char idx )
{
	return(loop->ns[idx]);// EVM in %
}

float _Complex *CellList::LoopGetPSCHptr( void )
{
	return(loop->psch);
}

float _Complex *CellList::LoopGetSSCHptr( void )
{
	return(loop->ssch);
}

float _Complex *CellList::LoopGetPBCHptr( void )
{
	return(loop->pbch);
}

float _Complex *CellList::LoopGetPBCHpilotRawPtr( void )
{
	return(loop->pilotsRaw);
}

float _Complex *CellList::LoopGetPBCHpilotDeSpreadPtr( void )
{
	return(loop->pilotsDeSpread);
}

char *CellList::LoopGetPBCHmsgPtr( void )
{
	return(loop->pbchMsg);
}

signed char CellList::LoopGetPBCHcrc( void )
{
	return(loop->pbchCRC);
}

void CellList::LoopAllocSparseCor( void )
{ // ok for FDD, need to change this for TDD
	unsigned char cidx=0;
	// BEWARE! essentially this is a sparse matrix, we only store werr for the symbols useful to us.
	for(cidx=0;cidx<(12*SRCH_FRMS);cidx++)
		loop->sparseCor[cidx]=(float _Complex *)malloc(sizeof(float _Complex)*NFFT);
}

float _Complex **CellList::LoopGetSparseCor( void )
{
	return(loop->sparseCor);
}

float CellList::LoopGetFerr( unsigned char idx )
{
	return(loop->fErr[idx]);
}

void CellList::LoopSetFerr( float fErr,unsigned char idx )
{
	loop->fErr[idx]=fErr;
}

void CellList::LoopCalcFreqErr2( const char *fname )
{
	FILE *fp;
	fp=fopen( fname, "a" );
	loop->CalcdpErrdtAv();
	fprintf(fp,"LoopCalcFreqErr2:\n");
	loop->FprintPhVsT(fp,0);
	fclose(fp);
}

void CellList::LoopCalcSparseCor( void )
{
	loop->CalcSparseCor();
}

void CellList::LoopCalcdpErrdtAv( const char *fname )
{
	FILE *fp;
	fp=fopen( fname, "a" );
	loop->CalcdpErrdtAv();
	fprintf(fp,"LoopCalcdpErrdtAv:\n");
	loop->FprintPhVsT( fp,1 );
	fclose(fp);
}

double CellList::LoopGetdpErrdtAv( void )
{
	return(loop->dphErrdtAv);
}

void CellList::LoopClampPlots( void )
{
	loop->ClampPSCH();
	loop->ClampSSCH();
	loop->ClampPBCH();
}

/*void CellList::LoopResetEVM( void ) // not needed?
{
	loop->evmPSCH[0]=-10.0; // obsolete
}*/

/*void CellList::LoopGetEqv( float _Complex *eqv ) // not needed?
{
	unsigned char ci;
	for(ci=0;ci<64;ci++)
		eqv[ci]=loop->eqv[ci];
}*/

//add band
// convert freq to earfcn
// void CellList::AddLastCellPSCH( int time,unsigned char c2id,short fftF,float *p,signed char fftFIdxOff ) 
void CellList::AddLastCellPSCH( int time,unsigned char c2id,short fftF,float *p ) 
{ // at end of list
	//printf("CellList::AddLastCellPSCH time=%i c2id=%i fftF=%i evm=%.3f td=%.1f\n",time,c2id,fftF,100*sqrt(p[2]),p[3]);
	prv=last;
	last=new CellItem;
	last->nxt=NULL;
	last->prv=NULL;
	if(top==NULL)
		top=last;
	if( prv!=NULL )
	{
		prv->nxt=last;
		last->prv=prv;
		last->idx=prv->idx+1;
	}
	cur=last;
	cur->n2id=c2id;
	//cur->freq=?
	cur->LOfrq=LOfrq;
	cur->fftFIdx=fftF;
	cur->fftFIdxOff=fftFIdxOff; // {-1,0,1} * 5kHz e.g. 0=0kHz (7x15-5)=100kHz, (13*15+5)=200kHz
	cur->timePSCH[0]=time;
	cur->evm2PSCH[0]=p[2];
	cur->evmPSCH[0]=sqrtf(p[2])*100;
	cur->rrmsPSCH[0]=1.0f/sqrtf(p[0]);
	cur->pwrPSCH[0]=10*log10f(p[0]);
	cur->snrPSCH[0]=cur->pwrPSCH[0]-10*log10f(p[1]);
	cur->tdPSCH[0]=p[3];
	cur->phErrPSCH[0]=p[4];
	cur->cp0[0]=cp0[0];
	cur->cp0[1]=cp0[1];
	cur->cp1[0]=cp1[0];
	cur->cp1[1]=cp1[1];
	cur->r15=r15;
	cur->NFFTd2=NFFTd2;
	cur->rNFFT=rNFFT;
	memcpy(cur->dateStr,dateStr,16);
	memcpy(cur->time24Str,time24Str,16);
	cur->cRpt=cRpt;
	cur->subSwpNo=subSwpNo;
	cur->lms7002Temp=temp7002;
	sprintf(cur->cpuTemp,"%s",cpuTemp);
}

void CellList::LoopUpdatePSCH( int time,short cf,unsigned char c2id,unsigned char evmIdx,float *p )
{
//	unsigned char ci;
	//printf("LoopUpdatePSCH evm=%.1f loop->evmPSCH[%i]=%.1f %li\n",100*sqrt(p[2]),evmIdx,loop->evmPSCH[evmIdx]);
	if( (loop!=NULL) )
	{
		if( ((evmIdx==1) || (loop->fftFIdx==cf) || (loop->evm2PSCH[evmIdx]>p[2])))
		{ // when refining timing, we need to check fftFIdx+/-3 (+/-300kHz) for evmIdx==0
			//printf("CellList::LoopUpdatePSCH time=%i fftF=%i evm=%.3f td=%.1f\n",
			//	time,cf,100*sqrt(p[2]),p[3]);
			loop->bnd2pwr=bnd2pwr;
			loop->NFFT=NFFT;
			loop->fftFIdx=cf; // freq needs updating
//			loop->freq=frq+frqSampdNFFT*(f[0]-NFFTd2-1)+0.005*f[1];
//			loop->earfcn=LTEbandDLFreq2EARFCN(loop->lteBand,loop->freq);
			loop->timePSCH[evmIdx]=time;
			loop->evm2PSCH[evmIdx]=p[2];
			loop->evmPSCH[evmIdx]=sqrtf(p[2])*100;
			loop->rrmsPSCH[evmIdx]=1.0f/sqrtf(p[0]);
			loop->pwrPSCH[evmIdx]=10*log10f(p[0]);
			loop->snrPSCH[evmIdx]=loop->pwrPSCH[evmIdx]-10*log10f(p[1]);
			loop->tdPSCH[evmIdx]=p[3];
			loop->phErrPSCH[evmIdx]=p[4];
//			printf("%f %f %.2f %.2f %.2f\n",p[0],p[1],loop->pwr,loop->noise,loop->snr);
		}
	}
}

void CellList::LoopUpdateSSCH( int time,unsigned char c1id,unsigned char evmIdx,float *p,unsigned char ns,unsigned char ncp )
{
//	unsigned char ci;
	if( loop!=NULL )
	{
		//printf("evmIdx=%i ns=%i c1id=%i ncp=%i evm=%.3f\n",evmIdx,ns,c1id,ncp,p[2]);
		loop->n1id=c1id;
		loop->nid=loop->n2id+3*loop->n1id;
		if( evmIdx==0 )
		{
			loop->timeRef=time+0*(int)(p[3]+0.5f)-cp0[ncp]; // find pos of first OFDM symb of first det frame (Ns=0,Sym=0), inc. CP
			loop->timeRef-=(4+ncp)*(NFFT+cp1[ncp]);
			loop->timeRef-=ns*(6+ncp)*(NFFT+cp1[ncp]);
			loop->timeRef-=ns*(cp0[ncp]-cp1[ncp]); // calibrate time to slot0
			// correct ph for Td!=0
		}
		loop->timeSSCH[evmIdx]=time;
		loop->evmSSCH[evmIdx]=sqrtf(p[2])*100;
		loop->rrmsSSCH[evmIdx]=1.0f/sqrtf(p[0]);
		loop->pwrSSCH[evmIdx]=10*log10f(p[0]);
		loop->snrSSCH[evmIdx]=loop->pwrSSCH[evmIdx]-10*log10f(p[1]);
		loop->tdSSCH[evmIdx]=p[3];
		loop->phErrSSCH[evmIdx]=p[4];
		loop->ns[evmIdx]=ns;
		loop->ncp=ncp;
//		loop->Print();
	}
}
void CellList::UpdatePBCH(float pbchEVM,unsigned char P,unsigned char BW,unsigned char PHICHdur,unsigned char PHICHres,short nSysFrm,signed char crc,unsigned char evmIdx,float errAvPBCH,short pbchSymbCnt)
{
	if( loop!=NULL )
	{
		loop->PhichDur=PHICHdur;
		loop->PhichRes=PHICHres;
		loop->P=P;
		loop->evmPBCH=pbchEVM;
		loop->BW=LTEbw[BW];
		loop->NDLRB=LTEndlrb[BW];
		loop->nSysFrm=nSysFrm;
		loop->pbchCRC=crc;
		loop->pbchIdx=evmIdx; // location of PBCH from start of capture in half frames.
		loop->errAvPBCH=errAvPBCH;
		loop->pbchSymbCnt=pbchSymbCnt;
	}
}

short CellList::CountGoodCRCs( void )
{
	short crcs=0;
	cur=top;
	while( cur!=NULL )
	{
		crcs+=(cur->pbchCRC+1); // CRC={1,-1} => crcs+={2,0}
		cur=cur->nxt;
	}
	crcs>>=1; // divide by 2
	cur=NULL;
	return(crcs);
}

short CellList::CountBadCRCs( void )
{
	short crcs=0;
	cur=top;
	while( cur!=NULL )
	{
		crcs+=(1-cur->pbchCRC); // CRC={1,-1} => crcs+={2,0}
		cur=cur->nxt;
	}
	crcs>>=1; // divide by 2
	cur=NULL;
	return(crcs);
}

void CellList::MergeLists( CellList *ext )
{
	int curTmpCid=0;
	//printf("CellList::MergeLists\n");
	prv=last;
	nxt=ext->top;
	if( prv!=NULL )
	{
		if( nxt!=NULL )
		{
			last=ext->last;
			prv->nxt=nxt;
			nxt->prv=prv;
			curTmpCid=prv->tmpCid;
			cur=nxt;
			while(cur!=NULL) // update tmpCid's of new list so do not duplicate.
			{
				cur->tmpCid+=curTmpCid;
				cur=cur->nxt;
			}
		} // else adding empty list, nothing to do.
	}
	else
	{ // existing list empty, just add new list
		top=nxt;
		last=ext->last;
	}	
	ext->top=NULL; // clear ext list, to allow reuse.
	ext->last=NULL;
}

void CellList::CheckList( signed char fullList )
{
	if( top!=NULL )
		printf("TOP=%lX  TOP.PREV=%lX  TOP.NEXT=%lX\n",(unsigned long)top,(unsigned long)top->prv,(unsigned long)top->nxt);
	else
		printf("TOP=NULL\n");

	if( cur!=NULL )
		printf("CUR=%lX  CUR.PREV=%lX  CUR.NEXT=%lX\n",(unsigned long)cur,(unsigned long)cur->prv,(unsigned long)cur->nxt);
	else
		printf("CUR=NULL\n");

	if( prv!=NULL )
		printf("PRV=%lX  PRV.PREV=%lX  PRV.NEXT=%lX\n",(unsigned long)prv,(unsigned long)prv->prv,(unsigned long)prv->nxt);
	else
		printf("PRV=NULL\n");

	if( nxt!=NULL )
		printf("NXT=%lX  NXT.PREV=%lX  NXT.NEXT=%lX\n",(unsigned long)nxt,(unsigned long)nxt->prv,(unsigned long)nxt->nxt);
	else
		printf("NXT=NULL\n");

	if( best!=NULL )
		printf("BEST=%lX BEST.PREV=%lX BEST.NEXT=%lX\n",(unsigned long)best,(unsigned long)best->prv,(unsigned long)best->nxt);
	else
		printf("BEST=NULL\n");

	if( loop!=NULL )
		printf("LOOP=%lX LOOP.PREV=%lX LOOP.NEXT=%lX\n",(unsigned long)loop,(unsigned long)loop->prv,(unsigned long)loop->nxt);
	else
		printf("LOOP=NULL\n");

	if( last!=NULL )
		printf("LAST=%lX LAST.PREV=%lX LAST.NEXT=%lX\n",(unsigned long)last,(unsigned long)last->prv,(unsigned long)last->nxt);
	else
		printf("LAST=NULL\n");

	if( tmp!=NULL )
		printf("TMP=%li  TMP.PREV=%li  TMP.NEXT=%li\n",(unsigned long)tmp,(unsigned long)tmp->prv,(unsigned long)tmp->nxt);
	else
		printf("TEMP=NULL\n");
	if(fullList>0)
	{
		printf("Full List...\n");
		cur=top;
		while(cur!=NULL)
		{
			printf("IDX=%i(%lX,%lX,%lX)\n",cur->idx,(unsigned long)cur,(unsigned long)cur->prv,(unsigned long)cur->nxt);
			cur=cur->nxt;
		}
	}
}

void CellList::LoopPlotConstAll( void )
{
	unsigned char cv=0;
	unsigned char cj=0;
	char fName[256];
	short pbchCnt=0;
	sprintf(fName,"output/%i_%i_constALL_%i_%i",loop->lteBand,loop->earfcn,subSwpNo,cRpt);
	FILE *GpPipe;
	GpPipe=popen("gnuplot","w"); // init GnuPlot Pipe
	fprintf(GpPipe,"set term gif\n");
	fprintf(GpPipe,"set output \"%s.gif\"\n",fName);
	fprintf(GpPipe,"set title \"PBCH Bnd=%i Earfcn=%i\"\n",loop->lteBand,loop->earfcn);
//	fprintf(GpPipe,"set ylabel \"%s\"\n",yLab);
//	fprintf(GpPipe,"set xlabel \"%s\"\n",xLab);
	//fprintf(GpPipe,"set xrange [-1.5:1.5]\n");
	//fprintf(GpPipe,"set yrange [-1.5:1.5]\n");
	// pt 4 open sq, pt 5 sq, pt 6 open circ, pt 7 circ pt, 12 open dia
	fprintf(GpPipe,"set style line 1 lt 2 lc rgb \"orange\" lw 1 pt 12\n");
	fprintf(GpPipe,"set style line 2 lt 2 lc rgb \"red\" lw 1 pt 12\n");
	fprintf(GpPipe,"set style line 3 lt 2 lc rgb \"magenta\" lw 1 pt 12\n");
	fprintf(GpPipe,"set style line 4 lt 2 lc rgb \"cyan\" lw 1 pt 12\n");
	fprintf(GpPipe,"set style line 5 lt 2 lc rgb \"blue\" lw 1 pt 12\n");
	fprintf(GpPipe,"set style line 6 lt 2 lc rgb \"green\" lw 1 pt 12\n");
	fprintf(GpPipe,"set style line 7 lt 2 lc rgb \"brown\" lw 1 pt 12\n");
	fprintf(GpPipe,"set style line 8 lt 2 lc rgb \"purple\" lw 1 pt 12\n");
	fprintf(GpPipe,"set grid\n");
	fprintf(GpPipe,"plot");
	for( cv=0;cv<4;cv++ )
		fprintf(GpPipe," '-' using 1:2 with points ls %i title \"pbch%i\",\\\n",(cv+1),cv);
	fprintf(GpPipe," '-' using 1:2 with points ls %i title \"psch\",\\\n",5);
	fprintf(GpPipe," '-' using 1:2 with points ls %i title \"ssch\",\\\n",6);
	//fprintf(GpPipe," '-' using 1:2 with points ls %i title \"refRaw\",\\\n",7);
	fprintf(GpPipe," '-' using 1:2 with points ls %i title \"refDsprd\"\n",8);
	for( cv=0;cv<4;cv++ )
	{
		if( (cv<2) || ((cv==3) && (loop->ncp==0)) )
		{
			for( cj=0; cj<60; cj++ ) // includes 36+DC+36
				fprintf( GpPipe, "%f %f\n",crealf(loop->pbch[pbchCnt+cj]),cimagf(loop->pbch[pbchCnt+cj]) );
			fprintf(GpPipe,"e\n");
			pbchCnt+=60;
		}
		else
		{
			for( cj=0; cj<72; cj++ ) // includes 36+DC+36
				fprintf( GpPipe, "%f %f\n",crealf(loop->pbch[pbchCnt+cj]),cimagf(loop->pbch[pbchCnt+cj]) );
			fprintf(GpPipe,"e\n");
			pbchCnt+=72;
		}
	}
	for( cj=0; cj<63; cj++ ) // includes DC
		fprintf( GpPipe, "%f %f\n",crealf(loop->psch[cj]),cimagf(loop->psch[cj]) );
	fprintf(GpPipe,"e\n");
	for( cj=0; cj<63; cj++ ) // includes DC
		fprintf( GpPipe, "%f %f\n",crealf(loop->ssch[cj]),cimagf(loop->ssch[cj]) );
	fprintf(GpPipe,"e\n");
	short cp=12*(loop->P);
	if( cp==0 )
		cp=24;
	//for( cj=0; cj<cp; cj++ )
	//	fprintf( GpPipe, "%f %f\n",crealf(loop->pilotsRaw[cj]),cimagf(loop->pilotsRaw[cj]) );
	//fprintf(GpPipe,"e\n");
	for( cj=0; cj<cp; cj++ )
		fprintf( GpPipe, "%f %f\n",crealf(loop->pilotsDeSpread[cj]),cimagf(loop->pilotsDeSpread[cj]) );
	fprintf(GpPipe,"e\n");
	fflush(GpPipe);
	pclose(GpPipe); // kill gnuplot process!
}

