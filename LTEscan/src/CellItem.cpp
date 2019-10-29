/*
 Copyright 2018 Lime Microsystems Ltd.

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
#include <complex.h> // double _Complex, compliant with double, fftw_complex[2] Re=0, Im=1
#include "lime.h"
#include "CellItem.hpp"

extern short NFFT;

CellItem::CellItem()
{
	unsigned short ci=0;
	nxt=NULL;
	prv=NULL;
	idx=1;
	tmpCid=0;
	lteBand=0;
	earfcn=0;
	freq=0.0;
	fErr[0]=0.0f;
	fErr[1]=0.0f;
	sparseCor=(float _Complex **)malloc(sizeof(float _Complex *)*(12*SRCH_FRMS));
	for(ci=0;ci<(12*SRCH_FRMS);ci++)
		sparseCor[ci]=NULL;
	for(ci=0;ci<FERR_COR_HFRMS;ci++)
	{
		timePSCH[ci]=0;
		evmPSCH[ci]=200.0f;
		evm2PSCH[ci]=2.0f;
		pwrPSCH[ci]=-101.0f;
		snrPSCH[ci]=-101.0f;
		tdPSCH[ci]=0.0f;
		phErrPSCH[ci]=0.0f;
		rrmsPSCH[ci]=0.0f;

		timeSSCH[ci]=0;
		evmSSCH[ci]=200.0f;
		pwrSSCH[ci]=-101.0f;
		snrSSCH[ci]=-101.0f;
		tdSSCH[ci]=0.0f;
		phErrSSCH[ci]=0.0f;
		rrmsSSCH[ci]=0.0f;
		ns[ci]=0;
	}
	nSysFrm=0;
	pbchCRC=0;
	pbchIdx=0;
	errAvPBCH=0.0;
	sprintf(pbchMsg,"-----------");

	fftFIdxOff=0;
	fftFIdx=0;
	n1id=0;
	n2id=0;
	nid=0;
	ncp=0;
	P=0;
	PhichDur=0;
	PhichRes=0;
	for(ci=0;ci<64;ci++)
		psch[ci]=0.0f;
	for(ci=0;ci<64;ci++)
		ssch[ci]=0.0f;
	for(ci=0;ci<292;ci++)
		pbch[ci]=0.0f;
	rPI2=1.0f/float(PI2);
	cp0dNFFT_2pi[0]=float(0.5*PI); // phase shift of CP per subcarrier
	cp0dNFFT_2pi[1]=float((10*PI)/64); // =10/128*2*pi
	cp1dNFFT_2pi[0]=float(0.5*PI); // =32/128*2*pi
	cp1dNFFT_2pi[1]=float((9*PI)/64); // =9/128*2*pi
	poff=10.0f*log10(0.93f);
	cRpt=0;
	subSwpNo=0;
}

CellItem::~CellItem()
{
	unsigned char ci=0;
	for(ci=0;ci<(12*SRCH_FRMS);ci++)
		if(sparseCor[ci]!=NULL)
			free(sparseCor[ci]);
	free(sparseCor);
}

void CellItem::PrintQwkLk( void )
{ // quick look at early results
	printf( "%3i %2i %2i %5i %.3f EVMp %5.1f TDp %6.1f PhErrp %5.2f Timep %6i [%i %i]\n",
		idx,tmpCid,lteBand,earfcn,freq,evmPSCH[0],tdPSCH[0],phErrPSCH[0],timePSCH[0],fftFIdx,fftFIdxOff);
}

void CellItem::Print( void )
{
	printf( "%3i %2i %2i %5i %.3f EVM{p%.1f p%.1f s%.1f s%.1f b%.1f} PWR{p%.1f p%.1f s%.1f s%.1f} SNR{p%.1f p%.1f s%.1f s%.1f} TD{p%.1f p%.1f s%.1f s%.1f} PhErr{p%.2f p%.2f s%.2f s%.2f} TIME/p%i p%i s%i s%i/ rrms{p%.2g p%.2g s%.2g s%.2g} [%i %i] (NID %i %i %i NCP %i [%i %i])\n",
		idx,tmpCid,lteBand,earfcn,freq,
		evmPSCH[0],evmPSCH[1],evmSSCH[0],evmSSCH[1],evmPBCH, // { }
		pwrPSCH[0],pwrPSCH[1],pwrSSCH[0],pwrSSCH[1], // { }
		snrPSCH[0],snrPSCH[1],snrSSCH[0],snrSSCH[1], // { }
		tdPSCH[0],tdPSCH[1],tdSSCH[0],tdSSCH[1], // { }
		phErrPSCH[0],phErrPSCH[1],phErrSSCH[0],phErrSSCH[1], // { }
		timePSCH[0],timePSCH[1],timeSSCH[0],timeSSCH[1], // / /
		rrmsPSCH[0],rrmsPSCH[1],rrmsSSCH[0],rrmsSSCH[1], // { }
		fftFIdx,fftFIdxOff, // []
		nid,n1id,n2id,ncp,ns[0],ns[1] ); // ( )
}

void CellItem::Fprint( FILE *fp )
{
	fprintf( fp,"%3i %2i %2i %5i %.3f EVM{p%.1f p%.1f s%.1f s%.1f b%.1f} PWR{p%.1f p%.1f s%.1f s%.1f} SNR{p%.1f p%.1f s%.1f s%.1f} TD{p%.1f p%.1f s%.1f s%.1f} PhErr{p%.2f p%.2f s%.2f s%.2f} TIME/p%i p%i s%i s%i/ rrms{p%.2g p%.2g s%.2g s%.2g} [%i %i] (NID %i %i %i NCP %i [%i %i])\n",
		idx,tmpCid,lteBand,earfcn,freq,
		evmPSCH[0],evmPSCH[1],evmSSCH[0],evmSSCH[1],evmPBCH, // {}
		pwrPSCH[0],pwrPSCH[1],pwrSSCH[0],pwrSSCH[1], // {}
		snrPSCH[0],snrPSCH[1],snrSSCH[0],snrSSCH[1], // {}
		tdPSCH[0],tdPSCH[1],tdSSCH[0],tdSSCH[1], // {}
		phErrPSCH[0],phErrPSCH[1],phErrSSCH[0],phErrSSCH[1], // {}
		timePSCH[0],timePSCH[1],timeSSCH[0],timeSSCH[1], // / /
		rrmsPSCH[0],rrmsPSCH[1],rrmsSSCH[0],rrmsSSCH[1], // { }
		fftFIdx,fftFIdxOff, // []
		nid,n1id,n2id,ncp,ns[0],ns[1] ); // ( )
}

void CellItem::Print2( void )
{
	printf( "%2i %5i %.3f %2i %2i %6.1f %6.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5i %5i %5i %5i %3i %2i %3i %3i %3i %3i\n",
		lteBand,earfcn,freq,idx,tmpCid,
		evmPSCH[0],evmPSCH[1],
		pwrPSCH[0],pwrPSCH[1],pwrSSCH[0],pwrSSCH[1],
		snrPSCH[0],snrPSCH[1],snrSSCH[0],snrSSCH[1],
		timePSCH[0],timePSCH[1],timeSSCH[0],timeSSCH[1],
		fftFIdx,fftFIdxOff,
		nid,ncp,ns[0],ns[1] );
}

void CellItem::Fprint2( FILE *fp )
{
	fprintf( fp,"%2i %5i %.3f %2i %2i %6.1f %6.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5i %5i %5i %5i %3i %2i %3i %3i %3i %3i\n",
		lteBand,earfcn,freq,idx,tmpCid,
		evmPSCH[0],evmPSCH[1],
		pwrPSCH[0],pwrPSCH[1],pwrSSCH[0],pwrSSCH[1],
		snrPSCH[0],snrPSCH[1],snrSSCH[0],snrSSCH[1],
		timePSCH[0],timePSCH[1],timeSSCH[0],timeSSCH[1],
		fftFIdx,fftFIdxOff,
		nid,ncp,ns[0],ns[1] );
}

int CellItem::Sprint( char *str )
{
	int len=0;
	len=sprintf( str,"%3i %2i %2i %5i %.3f EVM{p%.1f p%.1f s%.1f s%.1f b%.1f} PWR{p%.1f p%.1f s%.1f s%.1f} SNR{p%.1f p%.1f s%.1f s%.1f} TD{p%.1f p%.1f s%.1f s%.1f} PhErr{p%.2f p%.2f s%.2f s%.2f} TIME/p%i p%i s%i s%i/ rrms{p%.2g p%.2g s%.2g s%.2g} [%i %i] (NID %i %i %i NCP %i [%i %i])\n",
		idx,tmpCid,lteBand,earfcn,freq,
		evmPSCH[0],evmPSCH[1],evmSSCH[0],evmSSCH[1],evmPBCH, // {}
		pwrPSCH[0],pwrPSCH[1],pwrSSCH[0],pwrSSCH[1], // {}
		snrPSCH[0],snrPSCH[1],snrSSCH[0],snrSSCH[1], // {}
		tdPSCH[0],tdPSCH[1],tdSSCH[0],tdSSCH[1], // {}
		phErrPSCH[0],phErrPSCH[1],phErrSSCH[0],phErrSSCH[1], // {}
		timePSCH[0],timePSCH[1],timeSSCH[0],timeSSCH[1], // / /
		rrmsPSCH[0],rrmsPSCH[1],rrmsSSCH[0],rrmsSSCH[1], // { }
		fftFIdx,fftFIdxOff, // []
		nid,n1id,n2id,ncp,ns[0],ns[1] ); // ( )
	return(len);
}

void CellItem::PrintSh( void )
{
	printf( "%2i %2i %5i %.1f %6.1f %6.1f %6.1f %6.1f %4.0f %4.0f %6.1f %6.1f %5.2f %5.2f %5i %5i %i %i %i\n",
		idx,tmpCid,earfcn,freq,
		evmPSCH[0],evmPSCH[1],
		pwrPSCH[0],pwrPSCH[1],
		snrPSCH[0],snrPSCH[1],
		tdPSCH[0],tdPSCH[1],
		phErrPSCH[0],phErrPSCH[1],
		timePSCH[0],timePSCH[1],
		fftFIdx,fftFIdxOff,
		n1id );
}

void CellItem::FprintSh( FILE *fp )
{
	fprintf( fp,"%2i %2i %5i %.1f %6.1f %6.1f %6.1f %6.1f %4.0f %4.0f %6.1f %6.1f %5.2f %5.2f %5i %5i %i %i %i\n",
		idx,tmpCid,earfcn,freq,
		evmPSCH[0],evmPSCH[1],
		pwrPSCH[0],pwrPSCH[1],
		snrPSCH[0],snrPSCH[1],
		tdPSCH[0],tdPSCH[1],
		phErrPSCH[0],phErrPSCH[1],
		timePSCH[0],timePSCH[1],
		fftFIdx,fftFIdxOff,
		n1id );
}
void CellItem::PrintSum( void )
{
	float pwr=(pwrPSCH[0]+pwrPSCH[1]+pwrSSCH[0]+pwrSSCH[1])/4; // approx 0.1dB error using this method
	float snr=(snrPSCH[0]+snrPSCH[1]+snrSSCH[0]+snrSSCH[1])/4;
	float evmpsch=sqrt((evmPSCH[0]*evmPSCH[0]+evmPSCH[1]*evmPSCH[1])/2);
	float evmssch=sqrt((evmSSCH[0]*evmSSCH[0]+evmSSCH[1]*evmSSCH[1])/2);
	float pwr2=pwr+bnd2pwr+poff; // 62*15kHz = 0.93Mhz
	printf( "%s %s ",dateStr,time24Str);
	printf( " %2i %6i %6.1f %3i %3i %3i %5.1f ",lteBand,earfcn,freq,nid,ncp,P,BW );
	if(PhichDur==0)
		printf( "  NORM  ");
	else
		printf( "  EXTD  ");
	if(PhichRes==0)
		printf( "  1/6 ");
	if(PhichRes==1)
		printf( "  1/2 ");
	if(PhichRes==2)
		printf( "   1  ");
	if(PhichRes==3)
		printf( "   2  ");
	printf( "%3i  %5.1f %5.1f %6.1f %6.1f %6.1f %6.3f  %6.1f %2i %2i\n",pbchCRC,pwr2,snr,evmpsch,evmssch,evmPBCH,fErr[0],LOfrq,subSwpNo,cRpt );
}

void CellItem::FprintSum( FILE *fp )
{
	unsigned char ci=0;
	float pwr=0.0f;
	float snr=0.0f;
	for(ci=0;ci<FERR_COR_HFRMS;ci++)
		pwr+=pwrPSCH[ci]+pwrSSCH[ci]; // approx 0.1dB error using this method
	pwr*=0.5f*rFERR_COR_HFRMS;
	for(ci=0;ci<FERR_COR_HFRMS;ci++)
		snr+=snrPSCH[0]+snrSSCH[0];
	snr*=0.5f*rFERR_COR_HFRMS;
	float evmpsch=sqrt((evmPSCH[0]*evmPSCH[0]+evmPSCH[1]*evmPSCH[1])/2);
	float evmssch=sqrt((evmSSCH[0]*evmSSCH[0]+evmSSCH[1]*evmSSCH[1])/2);
	float pwr2=pwr+bnd2pwr+poff; // 62*15kHz = 0.93Mhz
	fprintf( fp,"%s %s ",dateStr,time24Str);
	fprintf( fp," %2i %6i %6.1f %3i %3i %3i %5.1f ",lteBand,earfcn,freq,nid,ncp,P,BW );
	if(PhichDur==0)
		fprintf( fp,"  NORM  ");
	else
		fprintf( fp,"  EXTD  ");
	if(PhichRes==0)
		fprintf( fp,"  1/6 ");
	if(PhichRes==1)
		fprintf( fp,"  1/2 ");
	if(PhichRes==2)
		fprintf( fp,"   1  ");
	if(PhichRes==3)
		fprintf( fp,"   2  ");
	fprintf( fp,"%3i  %5.1f %5.1f %6.1f %6.1f %6.1f %6.3f  %6.1f %2i %2i\n",pbchCRC,pwr2,snr,evmpsch,evmssch,evmPBCH,fErr[0],LOfrq,subSwpNo,cRpt );
}

void CellItem::GenSummaryHTML( FILE *fp )
{
	unsigned char ci=0;
	float pwr=0.0f;
	float snr=0.0f;
	for(ci=0;ci<FERR_COR_HFRMS;ci++)
		pwr+=pwrPSCH[ci]+pwrSSCH[ci]; // approx 0.1dB error using this method
	pwr*=0.5f*rFERR_COR_HFRMS;
	for(ci=0;ci<FERR_COR_HFRMS;ci++)
		snr+=snrPSCH[0]+snrSSCH[0];
	snr*=0.5f*rFERR_COR_HFRMS;
	float evmpsch=sqrt((evmPSCH[0]*evmPSCH[0]+evmPSCH[1]*evmPSCH[1])/2);
	float evmssch=sqrt((evmSSCH[0]*evmSSCH[0]+evmSSCH[1]*evmSSCH[1])/2);
	float pwr2=pwr+bnd2pwr+poff; // 62*15kHz = 0.93Mhz
	fprintf( fp,"<H3>EARFCN=%i (%.3fMHz)</H3>\n",earfcn,freq);
	fprintf( fp,"%s %s ss%i r%i<P>\n",dateStr,time24Str,subSwpNo,cRpt);
	fprintf( fp,"NID=%i NCP=%i P=%i BW=%.1fMHz ",nid,ncp,P,BW);
	if(PhichDur==0)
		fprintf( fp,"PHICHd=NORM ");
	else
		fprintf( fp,"PHICHd=EXTD ");
	if(PhichRes==0)
		fprintf( fp," PHICHr=1/6 ");
	if(PhichRes==1)
		fprintf( fp," PHICHr=1/2 ");
	if(PhichRes==2)
		fprintf( fp," PHICHr=1 ");
	if(PhichRes==3)
		fprintf( fp," PHICHr=2 ");
	fprintf( fp," SFN=%i CRC[%i]=%i PBCHmsg={0x%s}<P>PWR=%.1fdBm/MHz SNR=%.1fdB EVMp=%.1f%% EVMs=%.1f%% EVMb=%.1f%% FreqErr=%.3fkHz (%.3fkHz) PBCHav=%.3f<P>LO=%.3fMHz NFFT=%i FFTidx-DC=%i Fshift=%4.1fkHz<P>", nSysFrm,pbchIdx,pbchCRC,pbchMsg,pwr2,snr,evmpsch,evmssch,evmPBCH,fErr[0],fErr[1],errAvPBCH,LOfrq,NFFT,(fftFIdx-(NFFT>>1)),(5.0*fftFIdxOff));
	fprintf( fp,"TDpsch[0]=%.1f TDssch[0]=%.1f Tref=%i timePSCH[0]=%i timeSSCH[0]=%i ns[0]=%i evmPSCH[0]=%.2f evmSSCH[0]=%.2f dPhdt=%.3f<P>\n",tdPSCH[0],tdSSCH[0],timeRef,timePSCH[0],timeSSCH[0],ns[0],evmPSCH[0],evmSSCH[0],dphErrdtAv);
	fprintf( fp,"lms7002=%.1foC<P>\n%s",lms7002Temp,cpuTemp );
	fprintf( fp,"<img src=\"./output/%i_%i_constALL_%i_%i.gif\"><P>",lteBand,earfcn,subSwpNo,cRpt);
}

void CellItem::CalcdpErrdtAv( void )
{
	unsigned char cidx=0;
	unsigned char cnt=0;
	double oldAv=0.0; 
	//printf("CellItem::CalcdpErrdtAv\n");
	dphErrdtAv=0.0;
	for(cidx=0;cidx<FERR_COR_HFRMS;cidx++)
		dphErrdt[cidx]=phErrPSCH[cidx];
	for(cidx=0;cidx<FERR_COR_HFRMS;cidx++) // vectorizes
		dphErrdt[cidx]-=phErrSSCH[cidx];
	for(cidx=0;cidx<FERR_COR_HFRMS;cidx++) // vectorizes
		dphErrdt[cidx]+=PI2*((dphErrdt[cidx]<-PI)-(PI<dphErrdt[cidx]));
	for(cidx=0;cidx<FERR_COR_HFRMS;cidx++) // vectorizes
		dphErrdt[cidx]+=PI2*((dphErrdt[cidx]<-PI)-(PI<dphErrdt[cidx]));
	for(cidx=0;cidx<FERR_COR_HFRMS;cidx++)
		dphErrdtAv+=dphErrdt[cidx];
	dphErrdtAv*=rFERR_COR_HFRMS; // this average might be complete rubbish if data noisy.
	oldAv=dphErrdtAv;
	for(cidx=0;cidx<FERR_COR_HFRMS;cidx++)
		dphErrdtSD[cidx]=(dphErrdt[cidx]-dphErrdtAv)*(dphErrdt[cidx]-dphErrdtAv); // get standard deviation
	dphErrdtSDav=0.0; // value <1.0 is good, >6.0 bad
	for(cidx=0;cidx<FERR_COR_HFRMS;cidx++)
		dphErrdtSDav+=dphErrdtSD[cidx];
	dphErrdtSDav*=rFERR_COR_HFRMS; // for SD to be useful must be less than Pi
	dphErrdtAv=0.0;
	for(cidx=0;cidx<FERR_COR_HFRMS;cidx++)
		if( dphErrdtSD[cidx]<4.7 )
		{
			dphErrdtAv+=dphErrdt[cidx];
			cnt++;
		}
	dphErrdtAv/=cnt;
	if(cnt==0)
		dphErrdtAv=oldAv; // ensure we always have an answer!

	dphErrdtAv-=PI2*(dphErrdtAv>PI);
	dphErrdtAv+=PI2*(dphErrdtAv<-PI);
	// reject highest standard deviation(s) from average
}

void CellItem::FprintPhVsT( FILE *fp,unsigned char idx )
{
	unsigned char cidx=0;
	short foff=fftFIdx-NFFTd2; // signed relative to DC
	double freqErr=0.0;
	float d=0.0f;
	fprintf(fp,"   dphdt dphdtSD phP    phS     evmP    evmS     rmsP       rmsS     tdP  tdS    ns pwrP  pwrS   snrP  snrS     \n");
	for(cidx=0;cidx<FERR_COR_HFRMS;cidx++)
		fprintf(fp,"%i %6.3f %6.5f %6.3f %6.3f %7.3f %7.3f %10.4g %10.4g %5.1f %5.1f %4i %5.1f %5.1f %5.1f %5.1f %6i %6i\n",cidx,dphErrdt[cidx],dphErrdtSD[cidx],phErrPSCH[cidx],phErrSSCH[cidx],evmPSCH[cidx],evmSSCH[cidx],rrmsPSCH[cidx],rrmsSSCH[cidx],tdPSCH[cidx],tdSSCH[cidx],ns[cidx],pwrPSCH[cidx],pwrSSCH[cidx],snrPSCH[cidx],snrSSCH[cidx],timePSCH[cidx],timeSSCH[cidx]);
	double kzk=foff*cp1dNFFT_2pi[ncp];
	kzk*=(rPI2);
	kzk-=floor(kzk);
	kzk*=PI2;
	freqErr=dphErrdtAv-foff*cp1dNFFT_2pi[ncp]; // radians
	d=floor(freqErr*rPI2);
	freqErr-=d*PI2;
	freqErr-=PI2*(freqErr>PI);
	freqErr+=PI2*(freqErr<-PI);	
	freqErr-=PI2*(freqErr>PI);
	freqErr+=PI2*(freqErr<-PI);
	freqErr*=(6+ncp)/PI; // convert to kHz, dph/(2*Pi)*14/1ms=7*dph/Pi kHz
	fErr[idx]+=freqErr; // accumulate with previous result (inital value 0.0)
	fprintf(fp,"Ferr[%i]=%.6fkHz Ferr=%.6fkHz dPh=%.4f foff=%i kzk=%f 2*PI*cp1/NFFT=%f SD=%f\n",idx,fErr[idx],freqErr,dphErrdtAv,foff,kzk,PI2*(9.0/128),dphErrdtSDav );
}

void CellItem::CalcSparseCor( void )
{
// applying doppler frequency and +/-5kHz shift
// werr*t+werr*(cp0+csym*(NFFT+cp1[NCP])+ns*(cp0[NCP]-cp1[NCP])+(6+NCP)*(NFFT+cp1[NCP]))
// since we only need werr*...%NFFT, this simplifies to 
// werr*t+werr*(cp0+csym*cp1[NCP]+ns*(cp0[NCP]-cp1[NCP])+(6+NCP)*(cp1[NCP])
// initial cp0 term is optional, and can be ignored
	unsigned char cidx=0;
	unsigned char chFrms=0;
	unsigned char cs=0;
	unsigned char csym=0;
	unsigned short ct=0;
	float Ferr=fErr[0]+5.0f*fftFIdxOff; // khz
	// Ferr*=PI/(6+NCP);
	float _Complex werr=-I*PI2*rNFFT*(Ferr*r15); // 1.92MS/s FFT=128 => tFFT = 1.0/15kHz
//	float _Complex wPh=-I*PI2*rNFFT;
	int ph=0;
//	printf("CalcSparseCor %6.3f %6.3f %6.3f %6.5f %6.5f %6.5f %6.5f\n",Ferr,fErr[0],5.0*fftFIdxOff,cimagf(werr),(Ferr*r15),rNFFT,cargf(cexpf(werr*NFFT)));
	cidx=0;
	for(chFrms=0;chFrms<(2*SRCH_FRMS);chFrms++) // for each half frame
	{
		cs=10*chFrms;
		for(csym=(4+ncp);csym<(6+ncp);csym++) // SSCH,PSCH
		{
			ph=cp0[ncp]+csym*(NFFT+cp1[ncp])+cs*((cp0[ncp]-cp1[ncp])+(6+ncp)*(NFFT+cp1[ncp])); // full
			//ph=cp0[ncp]+csym*cp1[ncp]+cs*((cp0[ncp]-cp1[ncp])+(6+ncp)*cp1[ncp]); // %NFFT ver, error for +/-5kHz
			//ph&=(NFFT-1); // ph%NFFT, note do not need to consider frames and slot+10 as repeats.
//			if(ncp==1)
//				printf("cidx=%i chFrms=%i cs=%i csym=%i ph=%i cp0=%i cp1=%i NCP=%i\n",cidx,chFrms,cs,csym,ph,cp0[ncp],cp1[ncp],ncp);
			for(ct=0;ct<NFFT;ct++)
				sparseCor[cidx][ct]=cexpf(werr*(ct+ph)); // works using full definition of ph (no NFFT truncation)
//				sparseCor[cidx][ct]=cexpf(werr*(ct+15*ph)); // |NFFT| multiply phase by x15, why?  error for +/-5kH
//				sparseCor[cidx][ct]=cexpf(werr*(ct)+wPh*ph);
			cidx++;
		}
		cs=10*chFrms+1;
		for(csym=0;csym<4;csym++) // x4 PBCH
		{
			ph=cp0[ncp]+csym*(NFFT+cp1[ncp])+cs*((cp0[ncp]-cp1[ncp])+(6+ncp)*(NFFT+cp1[ncp]));
			//ph=cp0[ncp]+csym*cp1[ncp]+cs*((cp0[ncp]-cp1[ncp])+(6+ncp)*cp1[ncp]);
			//ph&=(NFFT-1); // ph%NFFT
			//if(ncp==1)
			//	printf("cidx=%i chFrms=%i cs=%i csym=%i ph=%i cp0=%i cp1=%i NCP=%i\n",cidx,chFrms,cs,csym,ph,cp0[ncp],cp1[ncp],ncp);
			for(ct=0;ct<NFFT;ct++)
				sparseCor[cidx][ct]=cexpf(werr*(ct+ph)); // works using full definition of ph (no NFFT truncation)
//				sparseCor[cidx][ct]=cexpf(werr*(ct+15*ph)); // %NFFT multiply phase by x15, why?
//				sparseCor[cidx][ct]=cexpf(werr*(ct)+wPh*ph);
			cidx++;
		}
	}
}

void CellItem::ClampPSCH( void )
{
	short cp=0;
	for(cp=0;cp<64;cp++)
		if(cabsf(psch[cp])>2.5f)
			psch[cp]=2.5f*cexpf(I*cargf(psch[cp]));
}

void CellItem::ClampSSCH( void )
{
	short cp=0;
	for(cp=0;cp<64;cp++)
		if(cabsf(ssch[cp])>2.5f)
			ssch[cp]=2.5f*cexpf(I*cargf(ssch[cp]));
}

void CellItem::ClampPBCH( void )
{
	short cp=0;
	for(cp=0;cp<292;cp++)
		if(cabsf(pbch[cp])>2.5f)
			pbch[cp]=2.5f*cexpf(I*cargf(pbch[cp]));
}

