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

/* Algorithm

(A) Coarse Scan Search
(1) Do coarse scan 30.72MHz wide 50ms (5 frames, to catch SIB)
(2) Search for possible PSCHs
(3) Verify 5ms PSCH repeat
(4) Search for SSCH using EVM {N2ID, NCP} --> Log
(5) grab PBCH without a fresh data grab. --> most of the settings we need before calling LTE decoder?

Note min bandwidth 1.08MHz.  So anything extra within this is false alarms or contain freq errors

(B) Detailed Search
(1) For each EARFCN NID 
(2) Equalise
(3) Obtain PBCH
(4) Sniff for SIB1
(5) Sniff for RNTIs
(6) Decode SIB1 {EARFCN, MCC, MNC, c2id, TAC, NID, NCP, nRNTIs}--> Log 

Use multiple data reads to avoid buffer overwrite.
Use 16bit reads for fastest data transfer.
Beware of mechanical connector issues with big antenna.

How to take advantage of multi threading?

UK LTE Bands
Band 20 DL 791-821 FDD (30MHz) EE 02 VF 3 try 800
Band 8 DL 925-960 FDD (35MHz) 02 VF  try 940
Band 3 DL 1805-1880 FDD (75MHz) EE 02 VF 3 try 1830 -cal fail
Band 1 DL 2110-2170 FDD (60MHz) EE 02 VF 3
Band 7 DL 2620-2690 FDD (70MHz) EE VF try 2640
// Looks like we have 3 channels locally, 793.5MHz 5MHz 798? 5MHz, 806MHz 10MHz, 816MHz 10MHz
// we have log periodic antenna with more gain
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lime.h"
#ifdef USE_C99_COMPLEX
#include <complex.h> // double _Complex, compliant with double, fftw_complex[2] Re=0, Im=1
#include <fftw3.h>  // should come after complex.h
#else //#ifdef USE_FFTW_COMPLEX
#include <fftw3.h>  // Do before complex.h to force typedef double fftw_complex[2]
#include <complex.h> 
#endif

#include <sys/stat.h> // mkdir
#include <time.h> // tools for telling the time!
#include <errno.h> // EPERM etc
#include <sensors/sensors.h> // lm-sensors library
//#include "LimeSuite.h"
#include "CellItem.hpp"
#include "CellList.hpp"
#include "LTEscan.hpp"
#include "SDRlib.hpp"
#ifdef MXD_CMPLR // needed if C compiled with gcc/clang, and cpp compiled g++/clang++
extern "C" {
#endif
#include "earfcn.h"  // .c file
#ifdef MXD_CMPLR
}
#endif
#define USE_SDR

#ifdef USE_F32
	float _Complex *iqDataPtr=NULL;
#else
	int16_t *iqDataPtr=NULL;
#endif

unsigned short NFFT=512; // default number of FFT bins
short NRpt=2; // default number of repeat scans for each band

//unsigned short NFFTd300k=((NFFT-72)*15)/300+((((NFFT-72)*15)%300)>0); // 300k blocks in NFFT (ceiling), note 6RBs keep out

// NOTE: LimeSuite.h typedef double float_type;
 // LTE Band 20 DL 791-821 (30MHz)
// private global variables
char fNameStem[250]="./output/LTEscan"; // default output location - will add extention .txt etc, so less than 256.
//unsigned int NUSB=1420*10; // number of samples read per USB packet
//unsigned int Nlat=NUSB/NFFT+((NUSB%NFFT)>0); // ceil(NUSB/NFFT)
unsigned char bandList[256]; // currently less than 256 bands

unsigned int fCnt;
extern double frq_LPF; // MHz, IF bandwidth of LimeSDR
extern unsigned char Ch; // SDR Channel number. 0:1
extern char LNAnum; // off - use default
extern unsigned char OSR; // ADC runs faster than USB rate.  was 4
extern unsigned int gaindB; // 48, gain in dB 0:73 - int for library compatibility
extern double frq_Samp; // MS/s
extern float fUserCor; // User supplied frequency correction.  Usually -800:800 Hz

char cpuName[256]; // includes \n
char ramSize[128]; // includes \n
char linuxName[128]; // includes \n
char cpuTemp[256]; // includes \n
unsigned char cpuType=-1; // 0=Intel, 1=AMD, 2=ARM
extern char sdrDetails[128];
extern float temp7002;


LTEscan *MyLTEscan; // ( NFFT )

void DecCmdLine( int argc, char *argv[] );
void Scan( const char *flog,const char *fout,const char *fhtml );

void myCPUsys( void );
void myCPUtempARMpi( void );
void myCPUtempIntelAMD( void );

int main( int argc, char *argv[] )
{
	clock_t mytref = clock();
	char flog[256]; // our diagnostic
	char fout[256]; // customer file
	char fhtml[256]; // our diagnostic
	unsigned char ci=0;
	bandList[0]=8;
	bandList[1]=20;
	bandList[2]=0; // stop
	sprintf(flog,"./output/LTElog.txt"); // true linux, user should be able to specify full path
	sprintf(fout,"./output/LTEscan.txt");
	sprintf(fhtml,"./LTEscan.html"); // issue with imag paths
	//char fname[100];
	DecCmdLine(argc,argv); // Read Settings
	myCPUsys();
	if(cpuType<2)
		myCPUtempIntelAMD();
	else
		myCPUtempARMpi();
	sprintf(fout,"%s.txt",fNameStem); // true linux, user should be able to specify full path
	frq_Samp=0.015*NFFT; // frqSamp and NFFT locked together in LTE.
	MyLTEscan=new LTEscan( NFFT ); // LTEscan(NFFT)
	iqDataPtr=MyLTEscan->ShareIQptr();
#ifdef USE_SDR
	OpenSDR(flog); // updates LPF value for html, get SDR details
#endif
	FILE *fp=fopen(flog,"w");
	fprintf(fp,"CMD LN: ");
	for(ci=0;ci<argc;ci++) // make copy of command line with data.
		fprintf(fp,"%s ",argv[ci]);
	fprintf(fp,"\n");
	fprintf(fp,"SDR: %s\n",sdrDetails);
	fprintf(fp,"LINUX: %sCPU: %s%s",linuxName,cpuName,ramSize); // record one, does not change!
	fprintf(fp,"Initial Temperatures:\nLMS7002M=%.1foC\n",(double)temp7002);
	fprintf(fp,"%s\n",cpuTemp);
	fclose(fp);
	fp=fopen(fout,"w");
		fprintf(fp,"                    BND EARFCN FreqMHZ NID NCP Ant BwMHz PHICHd PHICHr CRC PWRdB SNRdB P_EVM%% S_EVM%% B_EVM%% FerrkHz fLO_MHz\n");
	fclose(fp); // Reset log file
	fp=fopen(fhtml,"w");
	fprintf(fp,"<H1>LTE Scan</H1>\n");
	fprintf(fp,"CMD LN: ");
	for(ci=0;ci<argc;ci++) // make copy of command line with data.
		fprintf(fp,"%s ",argv[ci]);
	fprintf(fp,"<P>\n");
	fprintf(fp,"LINUX: %s<P>\nCPU: %s<P>\n%s<P>\n",linuxName,cpuName,ramSize); // record one, does not change!
	fprintf(fp,"SDR: %s<P>\n",sdrDetails);
	fprintf(fp,"SDR: ADC=%.2fMS/s USB=%.2fMS/s OSR=%i CLKGEN=%.2f LPF=%.2fMHz Gain=%idB LNA=%i<P>LTE: NFFT=%i SrchFrms=%i<P>\n",OSR*frq_Samp,frq_Samp,OSR,(4*OSR*frq_Samp),frq_LPF,gaindB,LNAnum,NFFT,SRCH_FRMS);
	fprintf(fp,"<H3>Initial Temperatures:</H3>\nLMS7002M=%.1foC<P>\n",(double)temp7002);
	fprintf(fp,"%s\n",cpuTemp);
	fclose(fp); // Reset html file
#ifdef USE_SDR
	Scan((const char*)flog,(const char*)fout,(const char*)fhtml); // (const char*)
	CloseSDR();
#else
	MyLTEscan->MakeTestSignal( 7,1 ); // NID,NCP
	MyLTEscan->CheckTestSignal( 1 ); // NCP
	MyLTEscan->UseWFM();
	MyLTEscan->CheckTestSignal( 1 ); // NCP
#endif
	if(cpuType<2)
		myCPUtempIntelAMD();
	else
		myCPUtempARMpi();
	fp=fopen(fhtml,"a"); // final temperature readout in case no successful detection made
	fprintf(fp,"<H3>Final Temperatures:</H3>\nLMS7002M=%.1foC<P>\n",(double)temp7002);
	fprintf(fp,"%s\n",cpuTemp);
	fprintf(fp,"Total Scan Time %.1fs<P>\n",((double)(clock()-mytref))/CLOCKS_PER_SEC);
	fclose(fp); // Reset html file
	fp=fopen(flog,"a"); // final temperature readout in case no successful detection made
	fprintf(fp,"Final Temperatures:\nLMS7002M=%.1foC\n",(double)temp7002);
	fprintf(fp,"%s\n",cpuTemp);
	fprintf(fp,"Total Scan Time %.1fs\n",((double)(clock()-mytref))/CLOCKS_PER_SEC);
	printf("Done!  Total Scan Time %.1fs\n",((double)(clock()-mytref))/CLOCKS_PER_SEC);
	fclose(fp); // Reset html file
	delete MyLTEscan;
	printf("Done!\n");
}

void Scan( const char *flog,const char *fout,const char *fhtml )
{
	double frq=0.0;
	double frq_min=791.0; // MHz
	double frq_max=821.0;  // MHz
	double frq_Step=30.0;
	double frq_HalfStep=0.0;
	unsigned char cf=0;
	unsigned char cBndIdx=0;
	short cRpt=0;
	int FrameSize=150*NFFT;
	short good=0;
	FILE *fp=NULL;
//	struct tm *myTime=NULL;
	printf("Scan starting...\n");
	while( bandList[cBndIdx]!=0 )
	{ // need to convert band number to frq_min and center frequency
		frq_min=LTEbandDLLow( bandList[cBndIdx] )+0.7; // EARFCNs close to band edge not used 136.104 5.7.3 Note 1
		frq_max=LTEbandDLHigh( bandList[cBndIdx] )-0.6; // EARFCNs close to band edge not used 136.104 5.7.3 Note 1
		// note, due to decimation, FIR filters cause amplitude damage to signal, reduce step size for good EVM
		fCnt=(int)ceil(LTEbandBW( bandList[cBndIdx] )/(0.75*frq_Samp)); // probably want sub-bands to overlap
		frq_Step=LTEbandBW( bandList[cBndIdx] )/fCnt;
		frq_HalfStep=floor(10*frq_Step/2)/10; // align center frequency on 100kHz grid regardless of NFFT
		frq_Step=2*frq_HalfStep;
// 30.72
// (frq_HalfStep-frq_min)/0.1
		frq=0.1*floor((frq_min+frq_HalfStep)*10);// must be on 100kHz earfcn grid
		printf("Band %i %.1fMHz:%.1fMHz:%.1fMHz %i %.3f\n",bandList[cBndIdx],frq_min,frq_Step,frq_max,fCnt,LTEbandBW( bandList[cBndIdx] ));
		for( cRpt=0; cRpt<NRpt; cRpt++ )
		{
			for( cf=0; cf<fCnt; cf++ )
			{
				if(cpuType<2)
					myCPUtempIntelAMD();
				else
					myCPUtempARMpi();
				MyLTEscan->StartOfBandMarkerForLog(flog,bandList[cBndIdx],cf,cRpt);
				frq=frq_min+frq_HalfStep+cf*frq_Step;
				good=SetFreqGrabFrames( flog,frq,FrameSize,SRCH_FRMS ); // 10ms Frame = NFFT*150 Samples
				if(good==1)
					MyLTEscan->ScanIqData( flog,bandList[cBndIdx],frq,cf,cRpt );
				else
				{
					MyLTEscan->LogMsg( flog,"\nSDR: Bad SDR read\n\n" );// else reject samples
					fp=fopen(fhtml,"a");
					fprintf(fp,"<B>ERROR: Scan at %6.1fMHz ss%i r%i failed after %i tries </B><P>\n",frq,cf,cRpt,MAX_SDR_READS);
					fclose(fp); // Reset html file
				}
			}
		}
		MyLTEscan->DisplayBandSummaries( flog,fout,fhtml,NRpt );
		cBndIdx++;
	}
	delete MyLTEscan;
	MyLTEscan=NULL;
}
// note offgridness is only a few possible values of phase shift
//B.000 * spot on
//B.090 * -10kHz error
//B.105 * +5kHz error
//B.195 * -5kH
//B.210 * + 10KHz
// 2048 (NFFT) pts, 30.72MS/s 2*pi*0.005/30.72
// how does evm phase align - through using tone +1/-1 as pilot


void DecCmdLine( int argc, char *argv[] )
{ // support a subset of Soapy Power for compatibility
	int ci=0;
	//char units1,units2,sep;
	//int scanLen;
	char csv[256];
	char csvStr[5];
	unsigned char csvStrCnt=0;
	unsigned char ccsv=0;
	unsigned char csvl=0;
	signed char quitOnHelp=0;
	unsigned char nfftok=0;
	//unsigned short osr[5]={32,32,16,8,4}; // gave calibration error failure
	unsigned short osr[5]={16,16,8,4,2};
	unsigned short nfft[5]={128,256,512,1024,2048};
	printf( "DecCmdLine\n" );
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
		if( strcmp(argv[ci],"-B" )==0 ) // list of LTE bands  1,2,3
		{ // use csv object to fool argc/argv to give bands as a single object
			ci++;
			sprintf(csv,"%s",argv[ci]);
			csvl=0;
			ccsv=0;
			while(csv[ccsv]!='\0') // parse csv object into individual bands
			{
				csvStrCnt=0;
				while((csv[ccsv]!='\0')&&(csv[ccsv]!=','))
					csvStr[csvStrCnt++]=csv[ccsv++];
				csvStr[csvStrCnt]='\0';
				bandList[csvl++]=(unsigned char)atoi(csvStr);
				if(csv[ccsv]==',')
					ccsv++;
			}
			bandList[csvl++]=0;
			csvl=0;
			//while(bandList[csvl]!=0)
			//	printf("%i ",bandList[csvl++]);
			//printf("\n");
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
		// IF filter defined by NFFT
		if( strcmp(argv[ci],"-OSR" )==0 ) // Oversample rate
		{
			ci++;
			OSR=atoi(argv[ci]);
		}
		if( strcmp(argv[ci],"-f" )==0 ) // User supplied frequency correction, usually -800:800 Hz
		{
			ci++;
			fUserCor=atof(argv[ci]);
		}
		if( strcmp(argv[ci],"-g" )==0 ) // Gain dB, Rx gain, default is 48dB for best dynamic range
		{
			ci++;
			gaindB=atoi(argv[ci]);
		}
		if( (strcmp(argv[ci],"-O" )==0) || (strcmp(argv[ci],"-o" )==0) ) // Output Filename
		{
			ci++;
			strcpy(fNameStem,argv[ci]);
		}
		// frequency ranges defined by bands
		if( (strcmp(argv[ci],"-v" )==0) || (strcmp(argv[ci],"--version" )==0) ) 
			printf("Version=%f\n",0.0); // Version number
		if( (strcmp(argv[ci],"-h" )==0) || (strcmp(argv[ci],"--help" )==0) )
		{
			printf("LimeScan [-h] [--help] [--info] [-v] [--version] [-b BINS] [-B BANDS] [-O FILEstub] [-o FILEstub] [-A ANTENNA] [-C CHANNEL] [-n REPEATS]\n"); // [-t SECONDS]
			quitOnHelp=1;
		}
		if( strcmp(argv[ci],"--info" )==0 ) // Soapy device info - not relevant
			printf("No information available, please consult your system administrator\n");
		ci++;
	}
	for( ci=0;ci<5;ci++)
		if(nfft[ci]==NFFT)
		{
			nfftok=1;
			OSR=osr[ci];
		}
	if( nfftok!=1 )
	{
		printf("ERROR: -b must be 128<= 2^N <=2048.  Forcing NFFT=256\n");
		NFFT=256;
	}
	frq_Samp=1.92*(NFFT>>7);
	//mkdir( fNameStem, S_IRWXU ); // S_IROTH Note: default ./output - for graphics and .xls; .csv to stdout

	if(quitOnHelp>0)
		exit(0);
}

void myCPUsys( void )
{ // discover about computer type and capabilities, affects how we measure CPU temp.
	FILE *fp; // create pipe from system command
	char sysline[256];
	// get processor description
	fp=popen("grep -m 1 'model name' /proc/cpuinfo | sed 's/model name\\s*:\\s*//' | sed 's/$/ /' | tr -d '\\n'","r");
	fgets(cpuName,sizeof(cpuName)-1,fp);
	pclose(fp);
//	fp=popen("grep -c 'cpu cores' /proc/cpuinfo | sed 's/$/ cores /' | tr -d '\\n'","r"); // get number of cores
//	fgets(sysline,sizeof(sysline)-1,fp);
//	sprintf(cpuName,"%s%s",cpuName,sysline); // BEWARE this can fail if cpuName is altered.  Fails on ARM
//	pclose(fp);
	fp=popen("lscpu | grep -m 1 'Stepping:' | sed 's/Stepping:\\s*//' | sed 's/$/ /'","r"); // get number of threads/cores
	fgets(sysline,sizeof(sysline)-1,fp);
	sprintf(cpuName,"%s stepping=%s",cpuName,sysline); // includes \n
	pclose(fp);
	fp=popen("lscpu | grep -m 1 'Socket(s):' | sed 's/Socket(s):\\s*//' | sed 's/$/ sockets,/'","r"); // get number of threads/cores
	fgets(sysline,sizeof(sysline)-1,fp);
	sprintf(cpuName,"%s. %s",cpuName,sysline); // includes \n
	pclose(fp);
	fp=popen("lscpu | grep -m 1 'Core(s) per socket:' | sed 's/Core(s) per socket:\\s*//' | sed 's/$/ cores\\/socket,/'","r"); // get number of threads/cores
	fgets(sysline,sizeof(sysline)-1,fp);
	sprintf(cpuName,"%s%s",cpuName,sysline); // includes \n
	pclose(fp);
	fp=popen("lscpu | grep -m 1 'Thread(s) per core:' | sed 's/Thread(s) per core:\\s*//' | sed 's/$/ threads\\/core./'","r"); // get number of threads/cores
	fgets(sysline,sizeof(sysline)-1,fp);
	sprintf(cpuName,"%s%s",cpuName,sysline); // includes \n
	pclose(fp);	
	char INTELstr[16]="Intel";
	char AMDstr[16]="AMD";
	char ARMstr[16]="ARM";
	if( strstr(cpuName,INTELstr)!=NULL )
		cpuType=0;
	if( strstr(cpuName,AMDstr)!=NULL )
		cpuType=1;
	if( strstr(cpuName,ARMstr)!=NULL )
		cpuType=2;
	printf("CPU=%i\n",cpuType);
	if( cpuType<2 ) // Intel and AMD have Cache
	{ // from lspcu we can also extract cache sizes, but only applies to Intel Machines
		fp=popen("lscpu | grep -m 1 'L1d cache:' | sed 's/L1d cache:\\s*//' | sed 's/$/ ,/'","r"); // get number of threads/cores
		sprintf(sysline,"0K\0"); // expect empty line
		fgets(sysline,sizeof(sysline)-1,fp);
		sprintf(cpuName,"%s L1d=%s",cpuName,sysline);
		pclose(fp);
		fp=popen("lscpu | grep -m 1 'L1i cache:' | sed 's/L1i cache:\\s*//' | sed 's/$/ ,/'","r"); // get number of threads/cores
		sprintf(sysline,"0K\0"); // expect empty line
		fgets(sysline,sizeof(sysline)-1,fp);
		sprintf(cpuName,"%s L1i=%s",cpuName,sysline);
		pclose(fp);
		fp=popen("lscpu | grep -m 1 'L2 cache:' | sed 's/L2 cache:\\s*//' | sed 's/$/ ,/'","r"); // get number of threads/cores
		sprintf(sysline,"0K\0"); // expect empty line
		fgets(sysline,sizeof(sysline)-1,fp); 
		sprintf(cpuName,"%s L2=%s",cpuName,sysline);
		pclose(fp);  // Intel atom and AMD Athalon has no L3 Cache
		fp=popen("lscpu | grep -m 1 'L3 cache:' | sed 's/L3 cache:\\s*//' | sed 's/$/ ./'","r"); // get number of threads/cores
		sprintf(sysline,"0K\0"); // expect empty line
		fgets(sysline,sizeof(sysline)-1,fp);
		sprintf(cpuName,"%s L3=%s",cpuName,sysline);
		pclose(fp);
	}
	fp=popen("grep -m 1 'MemTotal:' /proc/meminfo | sed 's/MemTotal:\\s*/RAM: /'","r"); // get total memory
	fgets(ramSize,sizeof(ramSize)-1,fp); // includes \n
	pclose(fp);
	fp=popen("grep -m 1 'MemAvailable:' /proc/meminfo | sed 's/MemAvailable:\\s*/ Available: /'","r"); // get available memory
	fgets(sysline,sizeof(sysline)-1,fp);
	sprintf(ramSize,"%s %s",ramSize,sysline);
	pclose(fp);
	fp=popen("lsb_release -d | grep -o ':.*' | sed 's/\\s*:\\s*//'","r"); // get linux version - included in the name!
	fgets(linuxName,sizeof(linuxName)-1,fp);
	pclose(fp);
	// ideally we want to monitor hard disk size too!  Detect disk flooding problems.
}

void myCPUtempARMpi( void )
{ // Rpi does not use lm-sensors
	FILE *fp; // create pipe from system command
	char pitemp[32];
	sprintf(cpuTemp,"");
	fp=popen("/opt/vc/bin/vcgencmd measure_temp","r");
	fgets(pitemp,sizeof(pitemp)-1,fp);
	sprintf(cpuTemp,"PI/CM%s<P>\n",pitemp); // BEWARE! Cannot append front of string with same name!!!
	pclose(fp);
}

void myCPUtempIntelAMD( void )
{ // uses lm-sensors library
	int cc=0;
	int err=0;
	double val;
	sensors_chip_name const *scname;
	sensors_cleanup();
	err=sensors_init(NULL);
	sprintf(cpuTemp,"");
	while( (scname=sensors_get_detected_chips(0,&cc)) != 0 )
	{
		sensors_feature const *sfeat;
		int cf=0;
		while( (sfeat=sensors_get_features(scname,&cf)) != 0 )
		{
			int csf=0;
			sensors_subfeature const *ssubfeat;
			if( (ssubfeat=sensors_get_all_subfeatures(scname,sfeat,&csf)) != 0 ) 
			{ // only need first subfeat for temp, rest critical thresholds
				if( ssubfeat->flags & SENSORS_MODE_R )
				{
					if( (err=sensors_get_value(scname,ssubfeat->number,&val))>=0)
						sprintf(cpuTemp,"%s%s(%i)=%.1foC ",cpuTemp,scname->prefix,cf,val);
				}
			}
		}
	}
	sprintf(cpuTemp,"%s<P>\n",cpuTemp);
}

