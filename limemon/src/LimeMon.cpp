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
// test case - pipe output, boundary post c
// ./LimeScan -f 800M:832M -C 0 -A "LNAW" -w 35M -r 16M -OSR 16 -b 256 -g 48.0 -n 16 | tail -10 > log.txt
// ./LimeScan -f 800M:831M -C 0 -A "LNAW" -w 35M -r 16M -OSR 16 -b 256 -g 48.0 -n 16 | tail -10 > log.txt
// ./LimeScan -f 800M:833M -C 0 -A "LNAW" -w 35M -r 16M -OSR 16 -b 256 -g 48.0 -n 16 | tail -15 > log.txt
// ./LimeScan -f 834M:883M -C 0 -A "LNAW" -w 35M -r 16M -OSR 8 -b 256 -g 48.0 -n 16 -Tst -50 | tail -15 > log.txt

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h> // double _Complex compliant with double, fftw_complex[2]; // [0] real, [1] imaginary
#include <fftw3.h>
#include <unistd.h> // pipes usleep
#include <fcntl.h> // file control
#include <sys/types.h> // pipes
#include <sys/stat.h> // mkdir
#include <time.h> // tools for telling the time!
#include <errno.h> // EPERM etc
#include "lime/LimeSuite.h"

using namespace std;

void DecodeArg( int argc, char *argv[] );
void Init( void );
void OpenSDR( void );
void ScanMode1( void );
void ScanMode2( void );
void GnuPlotDispAndSave( char *fName,double **zz );
void DisplayBinFile( char *fname );
void CloseSDR( void );
int error( void );
void ShutDown( void );
int ReadPwl( char fname[],double xlim[],double mPwl[],double cPwl[] );
void ApplyPwl( double *vecy[],double xlim[],double mPwl[],double cPwl[],unsigned char nPwl );
void ReadNMEAtraffic( char *gpsPath );
char* ReadTil( char *buffer, char word[] );
void ReadSubWord(char word[],char subWord[],int pos,int len);
void ReadGPRMC( char *buffer );

// global control variables
float frq_cen=70.0e6; // 800 MHz
float frq_step=16.0e6; // 16 MHz
float frq_LPF=10.0E6; // ifbandwidth of LimeSDR
float frq_Samp=frq_step; // USB Sample Rate for LimeSDR
float t_Tic=0.12; // ms
unsigned char OSR=4;
unsigned int gaindB=48; //gain in dB - has to be int for library compatibility
unsigned char Ch=0;
unsigned int NFFT=256;
unsigned int NRpt=32;
char LNAnum=3; // off - use default 1=LNAH 2=LNAL 3=LNAW
// private global variables
unsigned int NUSB=1420*10; // number of samples read per USB packet
unsigned int Nlat=NUSB/NFFT+((NUSB%NFFT)>0); // ceil(USBsize/FFTsize)
unsigned int tCnt;
lms_device_t* device=NULL; // SDR device
lms_stream_t streamId; // SDR stream
char fNameStem[100]="output";
fftw_plan pfft;
double _Complex *in=NULL; // registered with FFTW
double _Complex *out=NULL; // registered with FFTW
double *win=NULL; // hamming window coefficients
time_t *swpTime;
double time_len=500e-3;
double **fmat; // freq values for interpolation with pwl
double *tvec; // time values for plotting
double **z; // average data
double RNFFTdB;
char gpsPath[20]="/dev/ttyACM0";
unsigned char doGPS=0; // 0=don't 1=do
unsigned int tstLvl=-50;
double tstFrq=860e6; // 860MHz is an amature band, dont want to jam sensitive bands
unsigned char doTst=0; // 0=don't 1=do

char doPwlLNAW=0;
char doPwlLNAH=0;
char doPwlLNAL=0;
char doPwlAnt=0;
char fPwlLNAW[70]="pwl/lnaw.txt";
char fPwlLNAH[70]="pwl/lnah.txt";
char fPwlLNAL[70]="pwl/lnal.txt";
char fPwlAnt[70]="pwl/ant.txt";

int main( int argc, char *argv[] )
{
	unsigned char nPwl;
	double xPwl[255];
	double mPwl[255];
	double cPwl[255];
	char fname[100];
	DecodeArg(argc,argv);
	if( doGPS>0 )
		ReadNMEAtraffic( gpsPath );
	Init();
	OpenSDR();
	ScanMode1(); // long non contiguous, average option - seg fault
//	ScanMode2(); // short contiguous, no average
	CloseSDR();
/*	if( doPwlLNAW>0 )
	{
		nPwl=ReadPwl( fPwlLNAW,xPwl,mPwl,cPwl );
		ApplyPwl( zrms,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmin,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmax,xPwl,mPwl,cPwl,nPwl );
	}
	if( doPwlLNAH>0 )
	{
		nPwl=ReadPwl( fPwlLNAH,xPwl,mPwl,cPwl );
		ApplyPwl( zrms,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmin,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmax,xPwl,mPwl,cPwl,nPwl );
	}
	if( doPwlLNAL>0 )
	{
		nPwl=ReadPwl( fPwlLNAL,xPwl,mPwl,cPwl );
		ApplyPwl( zrms,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmin,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmax,xPwl,mPwl,cPwl,nPwl );
	}
	if( doPwlAnt>0 )
	{
		nPwl=ReadPwl( fPwlAnt,xPwl,mPwl,cPwl );
		ApplyPwl( zrms,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmin,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmax,xPwl,mPwl,cPwl,nPwl );
	}*/
	sprintf(fname,"%sRMS",fNameStem);
	GnuPlotDispAndSave( fname,z ); // average
	sprintf(fname,"%s/%sRMS.gif",fNameStem,fNameStem);
	DisplayBinFile( fname );// do we get multiple windows, or updated window?
//	sprintf(fname,"%s/%sPk.gif",fNameStem,fNameStem);
//	DisplayBinFile( fname );// do we get multiple windows, or updated window?
//	sprintf(fname,"%s/%sMin.gif",fNameStem,fNameStem);
//	DisplayBinFile( fname );// do we get multiple windows, or updated window?
	ShutDown();
	return(0);
}

void Init( void )
{
	unsigned int ct,ci;
	swpTime=(time_t*)malloc(sizeof(time_t)*NFFT);
	in=(double _Complex *) fftw_malloc(sizeof(double _Complex)*NFFT); // align for SIMD
	out=(double _Complex *) fftw_malloc(sizeof(double _Complex)*NFFT); // align for SIMD
	win=(double *)malloc(sizeof(double)*NFFT);
	frq_step=frq_Samp;
	tCnt=ceil((time_len*frq_Samp)/(NFFT*NRpt)); // tfft=NFFT/Frq_Samp
	printf("time_len=%.3fms tstep=%.3fms tCnt=%i frq_Samp=%.3fMHz\n",time_len*1e3,((NFFT*NRpt)/frq_Samp)*1e3,tCnt,frq_Samp*1e-6);
	tvec=(double*)malloc(sizeof(double)*tCnt);
	fmat=(double**)malloc(sizeof(double*)*tCnt);
	z=(double**)malloc(sizeof(double*)*tCnt);
	double RFstep=frq_Samp/NFFT;
	for(ct=0;ct<tCnt;ct++)
		tvec[ct]=ct*NFFT/frq_Samp*NRpt;
	for(ct=0;ct<tCnt;ct++)
	{
		fmat[ct]=(double*)malloc(sizeof(double)*NFFT);
		for(ci=0;ci<NFFT;ci++) // freq matrix for pwl scaling
			fmat[ct][ci]=frq_cen+ci*RFstep-(frq_Samp/2);
		z[ct]=(double*)malloc(sizeof(double)*NFFT);
		for(ci=0;ci<NFFT;ci++) // freq matrix for pwl scaling
			z[ct][ci]=-1.0*(ct+ci); // debug test pattern
	}
	for(ci=0;ci<NFFT;ci++)
		win[ci]=1-cos((2*3.141592*ci)/(NFFT-1)); // Hann
//		win[ci]=25.0/46.0-(1-25.0/46.0)*cos((2*3.141592*ci)/(NFFT-1)); // Hamming
//		win[ci]=0.355768-0.487396*cos((2*3.141592*ci)/(NFFT-1))+0.144232*cos((4*3.141592*ci)/(NFFT-1))-0.012604*cos((6*3.141592*ci)/(NFFT-1)); // Nutall
//		win[ci]=0.3635819-0.4891775*cos((2*3.141592*ci)/(NFFT-1))+0.1365995*cos((4*3.141592*ci)/(NFFT-1))-0.0106411*cos((6*3.141592*ci)/(NFFT-1)); // Blackman Nutall
//		win[ci]=0.35875-0.48829*cos((2*3.141592*ci)/(NFFT-1))+0.14128*cos((4*3.141592*ci)/(NFFT-1))-0.01168*cos((6*3.141592*ci)/(NFFT-1)); // Blackman Harris
	pfft=fftw_plan_dft_1d(NFFT,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
}

void ShutDown( void )
{ // clean out the trash
	unsigned int ct;
	for(ct=0;ct<tCnt;ct++)
	{
		free(fmat[ct]);
		free(z[ct]);
	}
	free(win);
	free(swpTime);
	free(tvec);
	fftw_destroy_plan(pfft);
	fftw_free(in);
	fftw_free(out);
}

void ScanMode1( void )
{ // successive Read-FFT - non contiguous monitoring
	int samplesRead;
	unsigned int ci,ct,crpt;
	float _Complex *buf=(float _Complex*)malloc(sizeof(float _Complex)*NFFT); // for LimeSDR interface
	float *mb[4]; // 0 mag, 1 min mag, 2 max mag, 3 rms mag, vector form for speed
//	char fname[100];
	unsigned int NFFTd2=NFFT>>1;
	float RNRptdB=-10*log10(NRpt); // divide Sum( mag^2 )/NRpt in dBs
	float RNFFTdB=-20*log10(NFFT); // divide by NFFT in dBs
	double Tstep=(NFFT*NRpt)/frq_Samp;
	for( ci=0;ci<4;ci++)
		mb[ci]=(float *)malloc(sizeof(float)*NFFT);
	if(LMS_SetLOFrequency(device,LMS_CH_RX,0,frq_cen)!= 0)
		error();
	usleep(100000); // flush old data was 100us
	printf("tCnt=%i\n",tCnt);
	for( ct=0; ct<tCnt; ct++ )
	{
		for( ci=0;ci<NFFT;ci++) // vectorized
		{ // SIMD vector reset
			mb[0][ci]=0.0; // now
			mb[1][ci]=0.0; // rms
			mb[2][ci]=1.0; // min
			mb[3][ci]=0.0; // max
		}		
		swpTime[ct]=time(NULL);
		for( crpt=0; crpt<NRpt; crpt++ )
		{
//			for(ci=0;ci<Nlat;ci++) // discard latency data from prev freq - NOT WORKING
//				samplesRead=LMS_RecvStream(&streamId,buf,NFFT,NULL,NFFT);
			samplesRead=LMS_RecvStream(&streamId,buf,NFFT,NULL,NFFT);
			for(ci=0;ci<NFFT;ci++) // copy SDR buffer to FFTW and window (Vectorized)
				in[ci]=buf[ci]*win[ci];
			fftw_execute(pfft);
			for( ci=0;ci<NFFT;ci++) // vectorized
			{
				mb[0][ci]=creal(out[ci]*conj(out[ci]));
				mb[1][ci]+=mb[0][ci];
			}
			for( ci=0;ci<NFFT;ci++) // vectorized
			{
				mb[2][ci]=mb[0][ci]*(mb[0][ci]<mb[2][ci])+mb[2][ci]*(mb[0][ci]>=mb[2][ci]); // SIMD min
				mb[3][ci]=mb[0][ci]*(mb[0][ci]>mb[3][ci])+mb[3][ci]*(mb[0][ci]<=mb[3][ci]); // SIMD max
			}
		}
		tvec[ct]=ct*Tstep; // should derive from tstep
		for(ci=0;ci<NFFTd2;ci++) // Include FFT pos/neg swap
		{
			z[ct][ci+NFFTd2]=10*log10(mb[1][ci]+1e-20)+RNFFTdB+RNRptdB; // divide by NFFT
			z[ct][ci]=10*log10(mb[1][ci+NFFTd2]+1e-20)+RNFFTdB+RNRptdB; // divide by NFFT
		}
	}
	for( ci=0;ci<4;ci++)
		free(mb[ci]);
	free(buf);
}

void ScanMode2( void )
{ // Bulk Read, Bulk FFT - shorter contiguous monitoring, no average.

}

void GnuPlotDispAndSave( char *fName,double **zz )
{ // zz[NY][NX], fName != fNameStem, as we have several sets of files to plot/print
	unsigned int ci,cj,ct;
	char fname[100];
	FILE *ftdv, *fcsv;
	FILE *GpPipe;
	
	float fftsf=(frq_Samp/1.0e6)/NFFT;
	GpPipe=popen("gnuplot","w"); // init GnuPlot Pipe
	fprintf(GpPipe,"set term gif\n");
	fprintf(GpPipe,"set output \"%s/%s.gif\"\n",fNameStem,fName);
	fprintf(GpPipe,"set grid\n");
	fprintf(GpPipe,"set surface\n"); // not required?
//	fprintf(GpPipe,"set view 10,340\n");
//	fprintf(GpPipe,"set xrange [0:10]\n");
//	fprintf(GpPipe,"set yrange [0:10]\n");
	fprintf(GpPipe,"set zrange [-90:0]\n");
	fprintf(GpPipe,"set pm3d\n"); // add map for heat map
	fprintf(GpPipe,"set xlabel \"FFT MHz (Cen %.1f MHz)\"\n",frq_cen*1.0e-6);
	fprintf(GpPipe,"set ylabel \"Time ms\"\n");
	fprintf(GpPipe,"set xtics 2\n");
	fprintf(GpPipe,"set ytics %f\n",t_Tic);
//	fprintf(GpPipe,"set palette defined (-1 \"blue\", 0 \"green\", 1 \"yellow\", 2 \"red\")\n");
//	fprintf(GpPipe,"set palette defined (-3 \"dark-green\", -2 \"green\", -1 \"cyan\", 0 \"blue\", 1 \"purple\", 2 \"red\", 3 \"orange\", 4 \"white\")\n");
	fprintf(GpPipe,"set palette defined (-5 \"black\",-4 \"dark-green\",-3 \"forest-green\", -2 \"green\", -1 \"sea-green\", 0 \"cyan\", 1 \"blue\", 2 \"purple\", 3 \"pink\",4 \"red\", 5 \"orange\", 6 \"yellow\", 7 \"white\")\n");
//	fprintf(GpPipe,"set palette defined ( 0 0 0 0, 1 1 1 1 )\n"); // grey scale
//	fprintf(GpPipe,"set palette defined ( 0 0 0 0, 1 0 1 1 )\n"); // cyan scale
	fprintf(GpPipe,"set hidden3d\n");
	fprintf(GpPipe,"set contour base\n"); // base surface both
	fprintf(GpPipe,"splot '-' using 1:2:3 with pm3d title 'ch 1'\n");
	for(ct=0;ct<tCnt;++ct)
	{ // 1-D x y z format
		for(cj=0;cj<NFFT;++cj)
			fprintf(GpPipe,"%f %f %f\n",(1.0*cj-(NFFT/2))*fftsf,tvec[ct]*1.0e3,zz[ct][cj]);
		fprintf(GpPipe,"\n");
	}
	fprintf(GpPipe,"e\n");
	fflush(GpPipe);
	pclose(GpPipe); // kill gnuplot process!
	
	sprintf( fname, "%s/%s.xls",fNameStem,fName); 
	ftdv=fopen(fname,"w"); // standard xls tab delimited variable "\t"
	fprintf(ftdv,"\t\t");
	for(cj=0;cj<NFFT;++cj)
		fprintf(ftdv,"\t%i",cj);
	fprintf(ftdv,"\n");
	for(ct=0;ct<tCnt;++ct)
	{
		fprintf(ftdv,"%.3f\t%.3f\t%i\t",frq_cen*1e-6,tvec[ct]*1.0e3,ct);
		for(cj=0;cj<NFFT;++cj)
			fprintf(ftdv,"%.2f\t",zz[ct][cj]);
		fprintf(ftdv,"\n");
	}
	fclose(ftdv);

	sprintf( fname, "%s/%s.csv",fNameStem,fName);
	if( strcmp( fNameStem,"output")!=0 )
		fcsv=fopen(fname,"w"); // "soapy power" style csv ", "
	else
		fcsv=stdout;
	fprintf(fcsv,"%s\n",fname);
	for(ci=0;ci<tCnt;++ci)
	{
		struct tm *myTime=gmtime(&swpTime[ci]); // utc time
//		struct tm *myTime=localtime(&swpTime[ci]);
		char dateStr[100];
		char time24Str[100]; // yr-mnth-dy,24h time
		sprintf(dateStr,"%i-%i-%i",myTime->tm_year+1900,myTime->tm_mon,myTime->tm_mday);
		sprintf(time24Str,"%i:%i:%i",myTime->tm_hour,myTime->tm_min,myTime->tm_sec);
		fprintf(fcsv,"%s, %s, ",dateStr,time24Str);
		fprintf(fcsv,"%.3f, %.3f, %.1f, %i",(frq_cen+frq_step/2)*1e-6,(frq_cen+frq_step)*1e-6,NRpt*NFFT/frq_step*1e3,NFFT);
		for(cj=0;cj<NFFT;++cj)
			fprintf(fcsv,", %.2f",zz[ci][cj]);
		fprintf(fcsv,"\n");
	}
	if( fcsv!=stdout )
		fclose(fcsv);
}

void OpenSDR( void )
{
	int n; //Find devices
	int ci=0;
	float_type rate=frq_Samp;
	float_type rf_rate=frq_Samp*OSR;
	float_type gain; 
    lms_name_t antenna_list[10]; //large enough list for antenna names.
	if((n=LMS_GetDeviceList(NULL))<0) // Pass NULL to only obtain number of devices
		error();
	printf("Devices found: %i \n",n);
	if(n<1)
		error();
	lms_info_str_t* list = new lms_info_str_t[n]; // allocate device list storage
	if(LMS_GetDeviceList(list)<0)
		error();
	for(ci=0;ci<n;ci++) // print device list
		printf("%i:%s\n",ci,list[ci]);
	if(LMS_Open(&device,list[0],NULL)) //Open the first device
		error();
	delete [] list;
    //Initialize device with default configuration
    //Use LMS_LoadConfig(device, "/path/to/file.ini") to load config from INI
    if(LMS_Init(device) != 0)
        error();
    if(LMS_EnableChannel(device, LMS_CH_RX,Ch,true)!=0) // Rx, Channels 0:(N-1)
        error();
	printf("LOfreq=%.3fMHz\n",frq_cen/1.0e6);
	if(LMS_SetLOFrequency(device, LMS_CH_RX,Ch,frq_cen)!=0)
        error();
	//Alternatively, NULL can be passed to LMS_GetAntennaList() to obtain number of antennae
	if((n = LMS_GetAntennaList(device, LMS_CH_RX, 0, antenna_list)) < 0)
		error();
	printf("Ae:\n"); //print available antennae names
	for(ci = 0; ci < n; ci++)
		printf(" %i:%s",ci,antenna_list[ci]);
	if((n = LMS_GetAntenna(device, LMS_CH_RX, Ch)) < 0) //get currently selected antenna index
		error();
	printf("\nDefault: %i:%s, ",n,antenna_list[n]);
	if( LNAnum==1 ) // manually select antenna
		if(LMS_SetAntenna(device, LMS_CH_RX, Ch, LMS_PATH_LNAH) != 0) 
			error();
	if( LNAnum==2 ) // manually select antenna
		if(LMS_SetAntenna(device, LMS_CH_RX, Ch, LMS_PATH_LNAL) != 0) 
			error();
	if( LNAnum==3 ) // manually select antenna
		if(LMS_SetAntenna(device, LMS_CH_RX, Ch, LMS_PATH_LNAW) != 0) 
			error();
	if((n = LMS_GetAntenna(device, LMS_CH_RX, Ch)) < 0) //get currently selected antenna index
		error();
	printf("Selected: %i:%s\n",n,antenna_list[n]);
    if(LMS_SetSampleRate(device,frq_Samp,OSR) != 0) // oversampling in RF OSR x sample rate
        error();
	if(LMS_GetSampleRate(device, LMS_CH_RX, 0, &rate, &rf_rate) != 0)  // NULL can be passed
		error();
	printf("\nUSB rate: %f MS/s, ADC rate: %fMS/s\n", rate/1.0e6, rf_rate/1.0e6);
    //Example of getting allowed parameter value range
    //There are also functions to get other parameter ranges (check LimeSuite.h)
	lms_range_t range; //Get allowed LPF bandwidth range
	if(LMS_GetLPFBWRange(device,LMS_CH_RX,&range)!=0)
		error(); // RX LPF range: 1.400100 - 130.000000 MHz
//	printf("RX LPF bandwitdh range: %f - %f MHz\n", range.min/1.0e6,range.max/1.0e6);   
	if(LMS_SetLPFBW(device, LMS_CH_RX, Ch, frq_step)!=0)  //Configure LPF, bandwidth 8 MHz
		error(); 
	if(LMS_SetNormalizedGain(device,LMS_CH_RX,Ch,gaindB/73.0)!=0) //Set RX gain
		error();
	int err=0;
	if(LMS_SetGaindB(device,LMS_CH_RX,Ch,gaindB)!=0) // 0:73 breaks
	{
		printf("Ch=%i gaindB=%i err=%i\n",Ch,gaindB,err);
		error();
	}
	if(LMS_GetNormalizedGain(device,LMS_CH_RX,Ch,&gain)!=0) //normalized gain
		error();
	if(LMS_GetGaindB(device,LMS_CH_RX,Ch,&gaindB)!=0)
		error();
	printf("Normalized RX Gain: %f, RX Gain: %i dB\n",gain,gaindB);
	if(LMS_Calibrate(device,LMS_CH_RX,Ch,frq_step,0)!=0)
		error();
	// optional TX test signal
	if( doTst>0 )
	{
		if(LMS_EnableChannel(device,LMS_CH_TX,Ch,true)!=0)
			error();
		if(LMS_SetAntenna(device,LMS_CH_TX,Ch,LMS_PATH_TX1)!=0) 
			error();
		if(LMS_SetLOFrequency(device,LMS_CH_TX,Ch,tstFrq)!=0)
			error();
		if(LMS_SetTestSignal(device,LMS_CH_TX,Ch,LMS_TESTSIG_NCODIV8,0,0)!=0)
			error(); // freq offset = rate/NCODIV
		if(LMS_SetLPFBW(device,LMS_CH_TX,Ch,16.0E6)!=0)
			error(); // TX LPF range: 1.400100 - 130.000000 MHz
		if(LMS_SetGaindB(device,LMS_CH_TX,Ch,(79+tstLvl))!= 0) // 0:73
			error(); //tstLvl
		if(LMS_GetNormalizedGain(device,LMS_CH_TX,Ch,&gain)!=0) //normalized gain
			error();
		if(LMS_GetGaindB(device,LMS_CH_TX,Ch,&gaindB)!=0)
			error();
		printf("Normalized TX Gain: %f, TX Gain: %i dB\n",gain,gaindB);
		if(LMS_Calibrate(device,LMS_CH_TX,Ch,frq_step,0)!=0)
			error();
	}
	streamId.channel=Ch; //channel number
	streamId.fifoSize=1024*1024; //fifo size in samples
	streamId.throughputVsLatency=1.0; //optimize for max throughput
	streamId.isTx=false; //RX channel
	streamId.dataFmt=lms_stream_t::LMS_FMT_F32; //32-bit floats
	if(LMS_SetupStream(device,&streamId)!=0)
		error();
	LMS_StartStream(&streamId);
}

void CloseSDR( void )
{
	if(LMS_EnableChannel(device,LMS_CH_TX,Ch,false)!=0)
		error();
	if(LMS_SetGaindB(device,LMS_CH_TX,Ch,0)!= 0) // 0:73
		error(); // switch off Tx so not to interfere.
	LMS_StopStream(&streamId); //stream is stopped but can be started again with LMS_StartStream()
	LMS_DestroyStream(device, &streamId); //stream is deallocated and can no longer be used
	LMS_Close(device);
}

int error( void )
{
	printf("LimeSDR: ERROR\n");
	if(device!=NULL)
		LMS_Close(device);
	exit(-1);
}

void DecodeArg( int argc, char *argv[] )
{ // support a subset of Soapy Power for compatibility
	int ci=0;
	char units1;
	int scanLen;
	char fname[100];
	FILE *fp;
	while( ci<argc )
	{
		if( strcmp(argv[ci],"-A" )==0 ) // antenna/LNA name
		{
			ci++;
			if( strcmp(argv[ci],"LNAH" )==0 )
				LNAnum=1;
			if( strcmp(argv[ci],"LNAL" )==0 )
				LNAnum=2;
			if( strcmp(argv[ci],"LNAW" )==0 )
				LNAnum=3;
		}
		if( strcmp(argv[ci],"-b" )==0 ) // number of FFT bins
		{
			ci++;
			NFFT=atoi(argv[ci]);
		}
/*		if( strcmp(argv[ci],"-C" )==0 ) // channel number - always 0 for LimeSDRmini
		{
			ci++;
			Ch=atoi(argv[ci]);
		}*/
		if( strcmp(argv[ci],"-n" )==0 ) // repeats
		{
			ci++;
			NRpt=atoi(argv[ci]);
		}
		if( strcmp(argv[ci],"-w" )==0 ) // RF filter width Hz
		{
			ci++;
			scanLen=sscanf(argv[ci],"%f%c",&frq_LPF,&units1);
			if( scanLen>1 ) // detect ifunits used
			{
				if( units1=='k' )
					frq_LPF*=1e3;
				if( units1=='M' )
					frq_LPF*=1e6;
				if( units1=='G' )
					frq_LPF*=1e9;
			}
		}
		if( strcmp(argv[ci],"-OSR" )==0 ) // channel number
		{
			ci++;
			OSR=atoi(argv[ci]);
		}
		if( strcmp(argv[ci],"-r" )==0 ) // Sample Rate at USB S/s
		{
			ci++;
			scanLen=sscanf(argv[ci],"%f%c",&frq_Samp,&units1);
			if( scanLen>1 ) // detect ifunits used
			{
				if( units1=='k' )
					frq_Samp*=1e3;
				if( units1=='M' )
					frq_Samp*=1e6;
				if( units1=='G' )
					frq_Samp*=1e9;
			}
		}
		if( strcmp(argv[ci],"-g" )==0 ) // Gain dB
		{
			ci++;
			gaindB=atoi(argv[ci]);
		}
		if( (strcmp(argv[ci],"-O" )==0) || (strcmp(argv[ci],"-o" )==0) ) // Output Filename
		{
			ci++;
			strcpy(fNameStem,argv[ci]);
		}
		if( strcmp(argv[ci],"-f" )==0 ) // center frequency 
		{
			ci++;
			units1=' ';
			scanLen=sscanf(argv[ci],"%f%c",&frq_cen,&units1);
			if( scanLen>1 ) // detect ifunits used
			{
				if( units1=='k' )
					frq_cen*=1.0e3;
				if( units1=='M' )
					frq_cen*=1.0e6;
				if( units1=='G' )
					frq_cen*=1.0e9;
			}
		}
		if( (strcmp(argv[ci],"-v" )==0) || (strcmp(argv[ci],"--version" )==0) ) // Version number
			printf("Version=%f\n",0.0);
		if( (strcmp(argv[ci],"-h" )==0) || (strcmp(argv[ci],"--help" )==0) )
			printf("LimeScan [-h] [--help] [-f Hz:Hz] [-O FILEstub] [-o FILEstub] [--info] [-v] [--version] [-b BINS] [-w Hz] [-r S/s] [-n REPEATS] [-t SECONDS] [-A ANTENNA] [-C CHANNEL]\n");
		if( strcmp(argv[ci],"--info" )==0 ) // Soapy device info - not relevant
			printf("No information available, please consult your system administrator\n");
		if( strcmp(argv[ci],"-gps" )==0 )
		{
			ci++;
			strcpy(gpsPath,argv[ci]);
			doGPS=1;
		}
		if( strcmp(argv[ci],"-Tst" )==0 )
		{
			ci++;
			tstLvl=atoi(argv[ci]);
			doTst=1;
			printf("-Tst %i\n",tstLvl);
		}
		if( strcmp(argv[ci],"-T" )==0 ) // should scan for unit e.g. 100m
		{
			ci++;
			time_len=atof(argv[ci]);
			doTst=1;
			printf("-T %f\n",time_len);
		}
		ci++;
	}
	mkdir( fNameStem, S_IRWXU ); // S_IROTH Note: default ./output - for graphics and .xls; .csv to stdout
	sprintf(fname,"%s/%s_cmd_params.txt",fNameStem,fNameStem);
	fp=fopen(fname,"w");
	for(ci=0;ci<argc;ci++) // make copy of command line with data.
		fprintf(fp,"%s ",argv[ci]);
	fprintf(fp,"\n\n");
	printf("File=%s\n",fNameStem);
	if( LNAnum==1 )
		fprintf(fp,"LNAH Ch=%i\n",Ch);
	if( LNAnum==2 )
		fprintf(fp,"LNAL Ch=%i\n",Ch);
	if( LNAnum==3 )
		fprintf(fp,"LNAW Ch=%i\n",Ch);	
	fprintf(fp,"Freq %.3f MHz\n", frq_cen*1.0e-6 );
	if(time_len>=0.003) // s
		t_Tic=0.2; // ms
	if(time_len>=0.01) //s
		t_Tic=2; // ms
	if(time_len>=0.03) // s
		t_Tic=6; // ms
	if(time_len>=0.1) // s
		t_Tic=20; // ms
	if(time_len>=0.3) // s
		t_Tic=60; // ms
	if(time_len>=1) // s
		t_Tic=200; // ms
	fprintf(fp,"time=%.3fms tsamp=%.3fms\n",time_len*1e3,(NRpt*NFFT)/frq_Samp*1e3);
	fprintf(fp,"Rate=%.3f MS/s OSR=%i ADC=%.3fMS/s CLKGEN=%.3fMHz\n",frq_Samp/1.0e6,OSR,OSR*(frq_Samp*1.0e-6),4*OSR*(frq_Samp*1.0e-6));
	if( (4*OSR*(frq_Samp*1.0e-6))>640 )
		fprintf(fp,"ERROR: CLK_GEN 4*Rate*OSR must be < 640MS/s\n");
	fprintf(fp,"Gain=%idB\n",gaindB);
	fprintf(fp,"LPF=%.3f MHz (RF)\n",frq_LPF*1.0e-6);
	fprintf(fp,"FFTW %i bins, with Hann Window, NRpt=%i\n",NFFT,NRpt);
	fclose(fp);
}

void DisplayBinFile( char *fname ) // xdg-open can be used to open txt,pdf,png,jpg,gif,mpg,avi
{
	char cmd[255];
	int err;
	sprintf(cmd,"xdg-open %s",fname); // xdg-open preinstalled in Ubuntu 16.04
	printf("SYSTEM(\"%s\")\n",cmd);
	err=system(cmd); // alternatives numerous include feh,eog,shotwell,tycat,tiv,gnome-open,pxl
	if( err!=0 )
		printf("system() %i %s\n",err,strerror(errno));
}

int ReadPwl( char fname[],double xlim[],double mPwl[],double cPwl[] )
{
	double ylim[256];
	short ci;
	unsigned char limCnt=0;
	FILE *fp=NULL;
	if( (fp=fopen(fname,"r"))==NULL )
	{
		printf("PWL file %s not found\n",fname);
		exit(1);
	}
	while( 1 )
	{
		if( fscanf(fp,"%lf %lf",&xlim[limCnt],&ylim[limCnt])<2 )
			break;
		limCnt++;
	}
	if( limCnt<2 )
	{
		printf("ERROR - too few limits\n");
		exit(1);
	}
	for(ci=0;ci<(limCnt-1);ci++)
	{
		mPwl[ci]=(ylim[ci+1]-ylim[ci])/(xlim[ci+1]-xlim[ci]);
		cPwl[ci]=ylim[ci]-xlim[ci]*mPwl[ci];
	}
	return(limCnt-1);
}

void ApplyPwl( double *vecy[],double xlim[],double mPwl[],double cPwl[],unsigned char nPwl )
{
	int ci;
	int cj;
	unsigned char ck=0;
	for(cj=0;cj<tCnt;cj++)
		for(ci=0;ci<NFFT;ci++)
		{
			while( (xlim[ck+1]<fmat[cj][ci]) && (ck<nPwl) ) // check in interval range
				ck++;
			vecy[cj][ci]=fmat[cj][ci]*mPwl[ck]+cPwl[ck];
		}
}

void ReadGPRMC( char *buffer )
{ // time UTC, status, lat, long, speed, track, date, 
	float gpsTime[3]={0,0,0}; // 0 hr, 1 min, 2 sec
	int gpsDate[3]={0,0,0}; // 0 day, 1 mnth, 2 yr
//	char gpsStatus='N';
	short gpsLatDeg=0; // deg
	float gpsLatMin=0; // min
	char gpsLatRef=' ';
	short gpsLngDeg=0; // deg
	float gpsLngMin=0; // min
	char gpsLngRef=' ';
//	float gpsAlt=0; // m
	float gpsTrackT=0; // heading relative to true north
	float gpsSpeedK=0; // km/h

	char *buf=buffer;
	char word[20];
	char subWord[20];
	if( buf!=NULL )
	{ // hh mm ss.ss
		buf=ReadTil( buf, word );
		ReadSubWord(word,subWord,0,2);
		gpsTime[0]=atof(subWord);
		ReadSubWord(word,subWord,2,2);
		gpsTime[1]=atof(subWord);
		ReadSubWord(word,subWord,4,5);
		gpsTime[2]=atof(subWord);
	}
	if( buf!=NULL )
	{ // hh mm ss.ss
		buf=ReadTil( buf, word );
		if( word[0]=='V' )
			printf("WARNING\n");
	}
	if( buf!=NULL )
	{ // latitude
		buf=ReadTil( buf, word );
		ReadSubWord(word,subWord,0,2);
		gpsLatDeg=atoi(subWord);
		ReadSubWord(word,subWord,2,10);
		gpsLatMin=atof(subWord);
	}
	if( buf!=NULL )
	{ // N or S
		buf=ReadTil( buf, word );
		gpsLatRef=word[0];
	}
	if( buf!=NULL )
	{ // long
		buf=ReadTil( buf, word );
		ReadSubWord(word,subWord,0,2);
		gpsLngDeg=atoi(subWord);
		ReadSubWord(word,subWord,2,10);
		gpsLngMin=atof(subWord);
	}
	if( buf!=NULL )
	{ // E or W
		buf=ReadTil( buf, word );
		gpsLngRef=word[0];
	}
	if( buf!=NULL )
	{ // speed
		buf=ReadTil( buf, word );
		gpsSpeedK=atof(word); // knots
	}
	if( buf!=NULL )
	{ // track
		buf=ReadTil( buf, word );
		gpsTrackT=atof(word);
	}
	if( buf!=NULL )
	{ // data 
		buf=ReadTil( buf, word );
		ReadSubWord(word,subWord,0,2);
		gpsDate[0]=atoi(subWord);
		ReadSubWord(word,subWord,2,2);
		gpsDate[1]=atoi(subWord);
		ReadSubWord(word,subWord,4,2);
		gpsDate[2]=atoi(subWord);
	}
	printf("GPRMC %.0f:%.0f:%.2f UTC ",gpsTime[0],gpsTime[1],gpsTime[2]);
	printf("%i/%i/%i UK ",gpsDate[0],gpsDate[1],gpsDate[2]);
	printf("%io %.3f'%c %io %.3f'%c ",gpsLatDeg,gpsLatMin,gpsLatRef, gpsLngDeg,gpsLngMin,gpsLngRef);
	printf("vel=%.1fkt Hdg=%.3foT\n",gpsSpeedK,gpsTrackT);
}

char* ReadTil( char buffer[], char word[] )
{ // read til comma, newline or /0, occasionally get ,, * is start of checksum
	unsigned char ci=0;
	char *buf=buffer;
	word[0]='\0';
	while( ((*buf)!=',') && ((*buf)!='\n') && ((*buf)!='\0') && ((*buf)!='*') )
		word[ci++]=(*buf++);
	word[ci]='\0';
	if( ((*buf)=='\n') || ((*buf)=='\0') || ((*buf)=='*') )
		buf=NULL;
	else
		buf++;
	return(buf);
}

void ReadSubWord(char word[],char subWord[],int pos,int len)
{ // read len chars, or until end of string
	int ci=0;
	while( (ci<len) && (word[pos+ci]!='\0') )
	{
		subWord[ci]=word[pos+ci];
		ci++;
	}
	subWord[ci]='\0';
}

void ReadNMEAtraffic( char *gpsPath )
{
	int ci=0;
	int cj=0;
	int ck=0;
	char buf[256];
	size_t size;
	char byte='\0';
	char header[10];
	int fd=0;
	if( (fd=open(gpsPath, O_RDWR | O_NOCTTY | O_NDELAY ))<0 )
	{
		printf("Error - cannot open Ublox GPS %s\n",gpsPath );
		exit(1);
	}
	fcntl(fd,F_SETFL,0);
	for(ci=0; ci<150;ci++)
	{
		ck=0;
		byte='\0';
		while(byte!='\n')
		{
			while( (size=read(fd,&byte,1))<1 ); // wait til ready
			buf[ck++]=byte;
		}
		ck--;
		byte='\0';
		buf[ck]='\0';
		if( ck>0 ) // ignore empty lines, seems to send \n\n
		{
			for(cj=0;cj<6;cj++)
				header[cj]=buf[cj];
			header[6]='\0';
			if(strcmp(header,"$GPRMC")==0)
				ReadGPRMC( buf+7 );
		}
	}
	close(fd);
}

