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
#include <fftw3.h>  // should come after complex.h
#include <unistd.h> // pipes usleep
#include <fcntl.h> // file control
#include <sys/types.h> // pipes
#include <sys/stat.h> // mkdir
#include <time.h> // tools for telling the time!
#include <errno.h> // EPERM etc
#include "LimeSuite.h"

#define PI 3.14159265
#define PI2 (2*PI)


using namespace std;

void DecCmdLine( int argc, char *argv[] );
void Init( void );
void OpenSDR( void );
void Scan( void );
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

// default global control variables
float frq_min=800.0e6; // 800 MHz
float frq_max=1000.0e6; // 1000 MHz
float frq_step=16.0e6; // 16 MHz
float frq=0.0;
float frq_LPF=3.0e6; // IF bandwidth of LimeSDR
float frq_Samp=16.0e6; // USB Sample Rate for LimeSDR
float frq_Tic=10.0;
unsigned char OSR=4; // ADC runs faster than USB rate.
unsigned int gaindB=48; // gain in dB 0:73 - int for library compatibility
unsigned char Ch=0; // SDR Channel number. 0:1
unsigned int NFFT=256; // default number of FFT bins
unsigned int NRpt=16; // number of repeat reads for RMS integration
char LNAnum=3; // off - use default
// private global variables
char fNameStem[50]="output"; // default output location
unsigned int NUSB=1420*10; // number of samples read per USB packet
unsigned int Nlat=NUSB/NFFT+((NUSB%NFFT)>0); // ceil(NUSB/NFFT)
unsigned int fCnt;
lms_device_t* device=NULL; // SDR device
lms_stream_t streamId; // SDR stream
double *win=NULL; // hamming window coefficients
double _Complex *in=NULL; // registered with FFTW
double _Complex *out=NULL; // registered with FFTW
fftw_plan pfft; // FFTW
time_t *swpTime;
double **fmat; // freq values for interpolation with pwl (MHz)
double **zrms; // average data
double **zmin; // min data
double **zmax; // max data
double RNFFTdB;
char gpsPath[20]="/dev/ttyACM0";
unsigned char doGPS=0; // 0=don't 1=do
unsigned int tstLvl=-50; // for optional test signal
double tstFrq=860e6; // 860MHz is an amature band, dont want to jam sensitive bands
unsigned char doTst=0; // 0=don't 1=do, test signal defaul off.

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
	DecCmdLine(argc,argv);
	if( doGPS>0 )
		ReadNMEAtraffic( gpsPath );
	Init();
	OpenSDR();
	Scan();
	CloseSDR();
	if( (doPwlLNAH>0) && (LNAnum==1) ) // Apply PWL corrections
	{ // block pwl if for wrong LNA
		nPwl=ReadPwl( fPwlLNAH,xPwl,mPwl,cPwl );
		ApplyPwl( zrms,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmin,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmax,xPwl,mPwl,cPwl,nPwl );
	}
	if( (doPwlLNAL>0) && (LNAnum==2) )
	{
		nPwl=ReadPwl( fPwlLNAL,xPwl,mPwl,cPwl );
		ApplyPwl( zrms,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmin,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmax,xPwl,mPwl,cPwl,nPwl );
	}
	if( (doPwlLNAW>0) && (LNAnum==3) )
	{
		nPwl=ReadPwl( fPwlLNAW,xPwl,mPwl,cPwl );
		ApplyPwl( zrms,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmin,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmax,xPwl,mPwl,cPwl,nPwl );
	}
	if( (doPwlAnt>0) )
	{
		nPwl=ReadPwl( fPwlAnt,xPwl,mPwl,cPwl );
		ApplyPwl( zrms,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmin,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmax,xPwl,mPwl,cPwl,nPwl );
	} // display results	sprintf(fname,"%sRMS",fNameStem);
	sprintf(fname,"%sRMS",fNameStem);
	GnuPlotDispAndSave( fname,zrms ); // average
	sprintf(fname,"%sMin",fNameStem);
	GnuPlotDispAndSave( fname,zmin ); // max
	sprintf(fname,"%sPk",fNameStem);
	GnuPlotDispAndSave( fname,zmax ); // min
#ifdef USE_GNUPLOT
	sprintf(fname,"%s/%sRMS.gif",fNameStem,fNameStem);
	DisplayBinFile( fname ); // provides updated window
//	sprintf(fname,"%s/%sPk.gif",fNameStem,fNameStem);
//	DisplayBinFile( fname );
//	sprintf(fname,"%s/%sMin.gif",fNameStem,fNameStem);
//	DisplayBinFile( fname );
#endif
	ShutDown();
	return(0);
}

void Init( void )
{ // memory allocation and frequency settings
	unsigned int ci;
	in=(double _Complex *) fftw_malloc(sizeof(double _Complex)*NFFT); // align for SIMD
	out=(double _Complex *) fftw_malloc(sizeof(double _Complex)*NFFT); // align for SIMD
	win=(double *)malloc(sizeof(double)*NFFT);
	frq_step=frq_Samp;
	frq_min+=frq_step/2;
	frq_max-=frq_step/2;
	fCnt=ceil((frq_max-frq_min)/frq_step)+1;
	printf("fCnt=%i\n",fCnt);
	swpTime=(time_t*)malloc(sizeof(time_t)*fCnt); // was NFFT
	fmat=(double**)malloc(sizeof(double*)*fCnt);
	zrms=(double**)malloc(sizeof(double*)*fCnt);
	zmin=(double**)malloc(sizeof(double*)*fCnt);
	zmax=(double**)malloc(sizeof(double*)*fCnt);
	for(ci=0;ci<fCnt;ci++)
	{
		fmat[ci]=(double*)malloc(sizeof(double)*NFFT);
		zrms[ci]=(double*)malloc(sizeof(double)*NFFT);
		zmin[ci]=(double*)malloc(sizeof(double)*NFFT);
		zmax[ci]=(double*)malloc(sizeof(double)*NFFT);
	}
	for(ci=0;ci<NFFT;ci++) // calcualte Windowing function coefficients
		win[ci]=1-cos((PI2*ci)/(NFFT-1)); // Hann
//		win[ci]=25.0/46.0-(1-25.0/46.0)*cos((PI2*ci)/(NFFT-1)); // Hamming
//		win[ci]=0.355768-0.487396*cos((PI2*ci)/(NFFT-1))+0.144232*cos((2*PI2*ci)/(NFFT-1))-0.012604*cos((3*PI2*ci)/(NFFT-1)); // Nutall
//		win[ci]=0.3635819-0.4891775*cos((PI2*ci)/(NFFT-1))+0.1365995*cos((2*PI2*ci)/(NFFT-1))-0.0106411*cos((3*PI2*ci)/(NFFT-1)); // Blackman Nutall
//		win[ci]=0.35875-0.48829*cos((PI2*ci)/(NFFT-1))+0.14128*cos((2*PI2*ci)/(NFFT-1))-0.01168*cos((3*PI2*ci)/(NFFT-1)); // Blackman Harris
	pfft=fftw_plan_dft_1d(NFFT,
                reinterpret_cast<fftw_complex*>(in),
                reinterpret_cast<fftw_complex*>(out),
                FFTW_FORWARD,
                FFTW_ESTIMATE);
}

void ShutDown( void )
{ // deallocate memory, and other good house keeping
	unsigned int ci;
	for(ci=0;ci<fCnt;ci++)
	{
		free(fmat[ci]);
		free(zrms[ci]);
		free(zmin[ci]);
		free(zmax[ci]);
	}
	free(win);
	free(swpTime);
	fftw_destroy_plan(pfft);
	fftw_free(in);
	fftw_free(out);
}

void OpenSDR( void )
{ // software based SDR set up, allows use with LimeSDR and LimeSDRmini.
	int n; // Find devices
	int ci=0;
	unsigned int gaindBTx;
	float_type rate=frq_Samp;
	float_type rf_rate=frq_Samp*OSR;
	float_type gain;
    lms_name_t antenna_list[10]; // large enough list for antenna names.
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
    // Initialize device with default configuration
    // Use LMS_LoadConfig(device, "/path/to/file.ini") to load config from INI
    if(LMS_Init(device) != 0)
        error();
    if(LMS_EnableChannel(device, LMS_CH_RX,Ch,true)!=0) // Rx, Channels 0:(N-1)
        error();
     if(LMS_SetLOFrequency(device, LMS_CH_RX,Ch,frq_min)!=0)
        error();
	// Alternatively, NULL can be passed to LMS_GetAntennaList() to obtain antennae num
	if((n = LMS_GetAntennaList(device, LMS_CH_RX, 0, antenna_list)) < 0)
		error();
	printf("Ae:\n"); // print available antennae names
	for(ci = 0; ci < n; ci++)
		printf(" %i:%s",ci,antenna_list[ci]);
	if((n = LMS_GetAntenna(device, LMS_CH_RX, Ch)) < 0) // get selected antenna index
		error();
	printf("\nDefault: %i:%s, ",n,antenna_list[n]);
	if( LNAnum==1 )
		if(LMS_SetAntenna(device, LMS_CH_RX, Ch, LMS_PATH_LNAH) != 0)
			error();
	if( LNAnum==2 )
		if(LMS_SetAntenna(device, LMS_CH_RX, Ch, LMS_PATH_LNAL) != 0)
			error();
	if( LNAnum==3 )
		if(LMS_SetAntenna(device, LMS_CH_RX, Ch, LMS_PATH_LNAW) != 0)
			error();
	if((n = LMS_GetAntenna(device, LMS_CH_RX, Ch)) < 0) // get selected antenna index
		error();
	printf("Selected: %i:%s\n",n,antenna_list[n]);
    if(LMS_SetSampleRate(device,frq_Samp,OSR) != 0) // oversampling in RF OSR x sample rate
        error();
	if(LMS_GetSampleRate(device, LMS_CH_RX, 0, &rate, &rf_rate) != 0)  // can pass NULL
		error();
	printf("\nUSB rate: %f MS/s, ADC rate: %fMS/s\n", rate*1.0e-6, rf_rate*1.0e-6);
    //Example of getting allowed parameter value range
    //There are also functions to get other parameter ranges (check LimeSuite.h)
	lms_range_t range; //Get allowed LPF bandwidth range
	if(LMS_GetLPFBWRange(device,LMS_CH_RX,&range)!=0)
		error(); // RX LPF range: 1.400100 - 130.000000 MHz
//	printf("RX LPF bandwitdh range: %f - %f MHz\n", range.min*1.0e-6,range.max*1.0e-6);
	if(LMS_SetLPFBW(device, LMS_CH_RX, Ch, frq_step)!=0)  //Configure LPF, bandwidth 8 MHz
		error();
	if(LMS_SetNormalizedGain(device,LMS_CH_RX,Ch,gaindB/73.0)!=0) //Set RX gain
		error();
	int err=0;
	if(LMS_SetGaindB(device,LMS_CH_RX,Ch,gaindB)!=0) // 0:73
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
	if( doTst>0 )
	{ // optional TX test signal
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
		if(LMS_GetGaindB(device,LMS_CH_TX,Ch,&gaindBTx)!=0)
			error();
		printf("Normalized TX Gain: %f, TX Gain: %i dB\n",gain,gaindBTx);
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
	LMS_StopStream(&streamId); // stream is stopped, start again with LMS_StartStream()
	LMS_DestroyStream(device, &streamId); //stream can no longer be used
	LMS_Close(device);
}

int error( void )
{
	printf("LimeSDR: ERROR\n");
	if(device!=NULL)
		LMS_Close(device);
	exit(-1);
}
void Scan( void )
{
	int samplesRead;
	unsigned int ci,cf,ct;
	float _Complex *buf=(float _Complex*)malloc(sizeof(float _Complex)*NFFT); // LimeSDR
	float *mb[4]; // 0 mag, 1 min mag, 2 max mag, 3 rms mag, do in vector form for speed
	unsigned int NFFTd2=NFFT>>1;
	float RNRptdB=-10*log10(NRpt); // divide Sum( mag^2 )/NRpt in dBs
	float RNFFTdB=-20*log10(NFFT); // divide by NFFT in dBs
	double RFstep=frq_Samp/NFFT;
	for( ci=0;ci<4;ci++)
		mb[ci]=(float *)malloc(sizeof(float)*NFFT);
	for( cf=0; cf<fCnt; cf++ )
	{
		frq=frq_min+cf*frq_step;
		if(LMS_SetLOFrequency(device,LMS_CH_RX,0,frq)!= 0)
			error();
		usleep(100000); // flush old data, was 100us
		for( ci=0;ci<NFFT;ci++) // vectorized
		{ // SIMD vector reset
			mb[0][ci]=0.0; // now
			mb[1][ci]=0.0; // rms
			mb[2][ci]=1.0; // min
			mb[3][ci]=0.0; // max
		}
		for( ct=0; ct<NRpt; ct++ )
		{
			swpTime[cf]=time(NULL);
			samplesRead=LMS_RecvStream(&streamId,buf,NFFT,NULL,NFFT);
			for(ci=0;ci<NFFT;ci++) // copy SDR buffer to FFTW and window (Vectorized)
				in[ci]=buf[ci]*win[ci];
			fftw_execute(pfft);
			for( ci=0;ci<NFFT;ci++) // vectorized
			{
				mb[0][ci]=creal(out[ci]*conj(out[ci]));
				mb[1][ci]+=creal(out[ci]*conj(out[ci]));
			}
			for( ci=0;ci<NFFT;ci++) // vectorized
			{
				mb[2][ci]=mb[0][ci]*(mb[0][ci]<mb[2][ci])+mb[2][ci]*(mb[0][ci]>=mb[2][ci]); // SIMD min
				mb[3][ci]=mb[0][ci]*(mb[0][ci]>mb[3][ci])+mb[3][ci]*(mb[0][ci]<=mb[3][ci]); // SIMD max
			}
		}
		for(ci=0;ci<NFFT;ci++) // freq matrix for pwl scaling
			fmat[cf][ci]=(frq+ci*RFstep-(frq_Samp/2))*1.0e-6; // MHz
		for(ci=0;ci<NFFTd2;ci++) // Include FFT pos/neg swap
		{
			zrms[cf][ci+NFFTd2]=10*log10(mb[1][ci]+1e-20)+RNFFTdB+RNRptdB-gaindB; // /NFFT
			zrms[cf][ci]=10*log10(mb[1][ci+NFFTd2]+1e-20)+RNFFTdB+RNRptdB-gaindB; // /NFFT
			zmin[cf][ci+NFFTd2]=10*log10(mb[2][ci]+1e-20)+RNFFTdB-gaindB; // /NFFT
			zmin[cf][ci]=10*log10(mb[2][ci+NFFTd2]+1e-20)+RNFFTdB-gaindB; // /NFFT
			zmax[cf][ci+NFFTd2]=10*log10(mb[3][ci]+1e-20)+RNFFTdB-gaindB; // /NFFT
			zmax[cf][ci]=10*log10(mb[3][ci+NFFTd2]+1e-20)+RNFFTdB-gaindB; // /NFFT
		}
	}
	for( ci=0;ci<4;ci++)
		free(mb[ci]);
	free(buf);
}

void GnuPlotDispAndSave( char *fName,double **zz )
{ // zz[NY][NX], fName != fNameStem, as we have several sets of files to plot/print
	unsigned int ci,cj;
	char fname[100];
	FILE *ftdv, *fcsv;
	float fftsf=(frq_Samp*1.0e-6)/NFFT;
	float frq=frq_min;
	FILE *GpPipe;
#ifdef USE_GNUPLOT
	GpPipe=popen("gnuplot","w"); // init GnuPlot Pipe
	fprintf(GpPipe,"set term gif\n");
	fprintf(GpPipe,"set output \"%s/%s.gif\"\n",fNameStem,fName);
	fprintf(GpPipe,"set grid\n");
	fprintf(GpPipe,"set surface\n"); // not required?
//	fprintf(GpPipe,"set view 10,340\n");
//	fprintf(GpPipe,"set xrange [0:10]\n");
//	fprintf(GpPipe,"set yrange [0:10]\n");
//	fprintf(GpPipe,"set zrange [-90:0]\n");
	fprintf(GpPipe,"set pm3d\n"); // add map for heat map
	fprintf(GpPipe,"set xlabel \"FFT Band MHz\"\n");
	fprintf(GpPipe,"set ylabel \"FFT Center MHz\"\n");
	fprintf(GpPipe,"set xtics 2\n");
	fprintf(GpPipe,"set ytics %f\n",frq_Tic);
//	fprintf(GpPipe,"set palette defined (-1 \"blue\", 0 \"green\", 1 \"yellow\", 2 \"red\")\n");
//	fprintf(GpPipe,"set palette defined (-3 \"dark-green\", -2 \"green\", -1 \"cyan\", 0 \"blue\", 1 \"purple\", 2 \"red\", 3 \"orange\", 4 \"white\")\n");
	fprintf(GpPipe,"set palette defined (-5 \"black\",-4 \"dark-green\",-3 \"forest-green\", -2 \"green\", -1 \"sea-green\", 0 \"cyan\", 1 \"blue\", 2 \"purple\", 3 \"pink\",4 \"red\", 5 \"orange\", 6 \"yellow\", 7 \"white\")\n");
//	fprintf(GpPipe,"set palette defined ( 0 0 0 0, 1 1 1 1 )\n"); // grey scale
//	fprintf(GpPipe,"set palette defined ( 0 0 0 0, 1 0 1 1 )\n"); // cyan scale
	fprintf(GpPipe,"set hidden3d\n");
	fprintf(GpPipe,"set contour base\n"); // base surface both
	fprintf(GpPipe,"splot '-' using 1:2:3 with pm3d title 'ch 1'\n");
	for(ci=0;ci<fCnt;++ci)
	{ // 1-D x y z format
		for(cj=0;cj<NFFT;++cj)
			fprintf(GpPipe,"%f %f %f\n",(1.0*cj-(NFFT/2))*fftsf,frq*1.0e-6,zz[ci][cj]);
		fprintf(GpPipe,"\n");
		frq+=frq_step;
	}
	fprintf(GpPipe,"e\n");
	fflush(GpPipe);
	pclose(GpPipe); // kill gnuplot process!
#endif
	sprintf( fname, "%s/%s.xls",fNameStem,fName);
	ftdv=fopen(fname,"w"); // standard xls tab delimited variable "\t" format
	fprintf(ftdv,"\t\t");
	for(cj=0;cj<NFFT;++cj)
		fprintf(ftdv,"\t%i",cj);
	fprintf(ftdv,"\n");
	for(ci=0;ci<fCnt;++ci)
	{
		fprintf(ftdv,"%.3f\t%.3f\t%i\t",(frq_min+frq_step*(ci-0.5))/1e6,(frq_min+frq_step*(ci+0.5))/1e6,ci);
		for(cj=0;cj<NFFT;++cj)
			fprintf(ftdv,"%.2f\t",zz[ci][cj]);
		fprintf(ftdv,"\n");
	}
	fclose(ftdv);

	sprintf( fname, "%s/%s.csv",fNameStem,fName);
	if( strcmp( fNameStem,"output")!=0 )
		fcsv=fopen(fname,"w"); // "soapy power" style csv ", " format
	else
		fcsv=stdout;
	fprintf(fcsv,"%s\n",fname);
	for(ci=0;ci<fCnt;++ci)
	{
		struct tm *myTime=gmtime(&swpTime[ci]); // utc time
//		struct tm *myTime=localtime(&swpTime[ci]);
		char dateStr[100];
		char time24Str[100]; // yr-mnth-dy,24h time
		sprintf(dateStr,"%i-%i-%i",myTime->tm_year+1900,myTime->tm_mon,myTime->tm_mday);
		sprintf(time24Str,"%i:%i:%i",myTime->tm_hour,myTime->tm_min,myTime->tm_sec);
		fprintf(fcsv,"%s, %s, ",dateStr,time24Str);
		fprintf(fcsv,"%.3f, %.3f, %.1f, %i",(frq_min+frq_step*(ci-0.5))/1e6,(frq_min+frq_step*(ci+0.5))/1e6,frq_step/NFFT,NFFT); // frq_min, frq_max, fftbin_frq_width, buffer_length_samples
		for(cj=0;cj<NFFT;++cj)
			fprintf(fcsv,", %.2f",zz[ci][cj]);
		fprintf(fcsv,"\n");
	}
	if( fcsv!=stdout )
		fclose(fcsv);
}

void DecCmdLine( int argc, char *argv[] )
{ // support a subset of Soapy Power for compatibility
	int ci=0;
	char units1,units2,sep;
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
		if( strcmp(argv[ci],"-C" )==0 ) // channel number - always 0 for LimeSDRmini
		{
			ci++;
			Ch=atoi(argv[ci]);
		}
		if( strcmp(argv[ci],"-n" )==0 ) // repeats
		{
			ci++;
			NRpt=atoi(argv[ci]);
		}
		if( strcmp(argv[ci],"-w" )==0 ) // RF filter width Hz
		{
			ci++;
			scanLen=sscanf(argv[ci],"%f%c",&frq_LPF,&units1);
			if( scanLen>1 ) // detect if units used
			{
				if( units1=='k' )
					frq_LPF*=1e3;
				if( units1=='M' )
					frq_LPF*=1e6;
				if( units1=='G' )
					frq_LPF*=1e9;
			}
		}
		if( strcmp(argv[ci],"-OSR" )==0 ) // Oversample rate
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
		if( strcmp(argv[ci],"-f" )==0 ) // frequency range, FORTRAN style fmin:fmax
		{ // sopay power also supports a center frequency mode, we do not!
			ci++;
			units1=' ';
			units2=' ';
			scanLen=sscanf(argv[ci],"%f%c%c%f%c",&frq_min,&units1,&sep,&frq_max,&units2);
			if( scanLen>3 ) // detect if units used
			{
				if( units1=='k' )
					frq_min*=1.0e3;
				if( units1=='M' )
					frq_min*=1.0e6;
				if( units1=='G' )
					frq_min*=1.0e9;
				if( units2=='k' )
					frq_max*=1.0e3;
				if( units2=='M' )
					frq_max*=1.0e6;
				if( units2=='G' )
					frq_max*=1.0e9;
			}
		}
		if( (strcmp(argv[ci],"-v" )==0) || (strcmp(argv[ci],"--version" )==0) )
			printf("Version=%f\n",0.0); // Version number
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
		if( strcmp(argv[ci],"-pwla" )==0 )
			doPwlAnt=1;
		if( strcmp(argv[ci],"-pwlw" )==0 )
			doPwlLNAW=1;
		if( strcmp(argv[ci],"-pwlh" )==0 )
			doPwlLNAH=1;
		if( strcmp(argv[ci],"-pwll" )==0 )
			doPwlLNAL=1;
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
	fprintf(fp,"Freq Range %.3f:%.3f MHz\n", frq_min*1.0e-6,frq_max*1.0e-6 );
	if(((frq_max-frq_min)/1.0E6)>=10)
		frq_Tic=2.0;
	if(((frq_max-frq_min)/1.0E6)>=30)
		frq_Tic=10.0;
	if(((frq_max-frq_min)/1.0E6)>=100)
		frq_Tic=20.0;
	if(((frq_max-frq_min)/1.0E6)>=300)
		frq_Tic=100.0;
	if(((frq_max-frq_min)/1.0E6)>=1000)
		frq_Tic=200.0;
	if(((frq_max-frq_min)/1.0E6)>=2000)
		frq_Tic=500.0;
	fprintf(fp,"Rate=%.3f Ms/s OSR=%i\n",frq_Samp*1.0e-6,OSR);
	fprintf(fp,"%f\n",4*OSR*(frq_Samp*1.0e-6));
	if( (4*OSR*(frq_Samp*1.0e-6))>640 )
		fprintf(fp,"ERROR: Rate*OSR must be < 640MS/s\n");
	fprintf(fp,"Gain=%idB\n",gaindB);
	fprintf(fp,"LPF=%.3f MHz (RF)\n",frq_LPF*1.0e-6);
	fprintf(fp,"FFTW %i bins, with Hann Window, NRpt=%i\n",NFFT,NRpt);
	fclose(fp);
}

void DisplayBinFile( char *fname )
{ // xdg-open can open txt,pdf,png,jpg,gif,mpg,avi. Preinstalled in Ubuntu 16.04
	char cmd[255];
	int err;
	sprintf(cmd,"xdg-open %s",fname);
	printf("SYSTEM(\"%s\")\n",cmd);
	err=system(cmd); // alternatives include: feh,eog,shotwell,tycat,tiv,gnome-open,pxl
	if( err!=0 )
		printf("system() %i %s\n",err,strerror(errno));
}

int ReadPwl( char fname[],double xlim[],double mPwl[],double cPwl[] )
{ // file format 'freq MHz', 'level dB'
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
{ // linear interpolation between pwl table
	int ci;
	int cj;
	unsigned char ck=0;
	for(cj=0;cj<fCnt;cj++)
		for(ci=0;ci<NFFT;ci++)
		{
			while( (xlim[ck+1]<fmat[cj][ci]) && (ck<nPwl) ) // check in interval range
				ck++;
			vecy[cj][ci]-=fmat[cj][ci]*mPwl[ck]+cPwl[ck];
		}
}

void ReadGPRMC( char *buffer )
{ // time UTC, status, lat, long, speed, track, date,
	float gpsTime[3]={0,0,0}; // 0 hr, 1 min, 2 sec
	int gpsDate[3]={0,0,0}; // 0 day, 1 mnth, 2 yr
	short gpsLatDeg=0; // deg
	float gpsLatMin=0; // min
	char gpsLatRef=' ';
	short gpsLngDeg=0; // deg
	float gpsLngMin=0; // min
	char gpsLngRef=' ';
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

