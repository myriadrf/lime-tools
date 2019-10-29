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
#include <unistd.h>
#include "lime.h" // local Lime specfic defs
#include "lime/LimeSuite.h" // LimeSuite library defs, API etc
#include "SDRlib.hpp"

using namespace std;

#ifdef USE_F32 // store data from SDR
extern float _Complex *iqDataPtr;
#else
extern int16_t *iqDataPtr;
#endif

lms_device_t* device=NULL; // SDR device
lms_stream_t streamId; // SDR stream
double frq_LPF=30.0; // MHz, IF bandwidth of LimeSDR - auto calculated from sample rate.
unsigned char Ch=0; // Default SDR Channel number. 0:1
char LNAnum=2; // default LNAL {1,2,3} => {LNAH,LNAL,LNAW} Note LimeNet micro has no LNAW, LimeSDR mini has no LNAL
unsigned char OSR=4; // ADC runs faster than USB rate.  was 4
unsigned int gaindB=48; // gain in dB 0:73 - int for library compatibility, was 48 for LNA=0,TIA=-3,PGA=-3
float fUserCor=0.0; // User supplied frequency correction.
double frq_Samp=30.72; // MS/s
float temp7002=0.0;
char sdrDetails[128];
extern char cpuTemp[256];

void SDRupdateSPIreg(uint32_t addr,unsigned char bitsH,unsigned char bitsL,unsigned char val); // private
void SDRreadStatus( FILE *fp,int FrmSize,double frq,unsigned char frms, clock_t mytref2,lms_stream_status_t status,uint64_t *timeStamp,float *fifosizelog,uint32_t *droplog );

void SDRupdateSPIreg(uint32_t addr,unsigned char bitsH,unsigned char bitsL,unsigned char val)
{ // conform to LimeSuite help format, e.g. 0x0113 [1:0] = 2
	uint16_t reg;
	uint16_t mask=((1<<(bitsH-bitsL+1))-1)<<bitsL;
	uint16_t val2=((uint16_t)val)<<bitsL;
	int err; // ignore
	val2&=mask;
	err=LMS_ReadLMSReg(device,addr,&reg);
	reg&=(65535-mask);
	reg|=val2;
	err=LMS_WriteLMSReg(device,addr,reg);
}

void OpenSDR( const char *fname )
{ // software based SDR set up, allows use with LimeSDR and LimeSDRmini.
	FILE *fp;
	fp=fopen( fname, "a" );
	int n; // Find devices
	int ci=0;
	double frq=850; // MHz, our start default frequency
	float_type temp=0.0; // temperature 'C
	float_type rate=frq_Samp*1.0e6; // USB rate
	float_type rf_rate=frq_Samp*OSR*1.0e6; // DAC rate
	float_type gain; 
    lms_name_t antenna_list[10]; // large enough list for antenna names.
	if((n=LMS_GetDeviceList(NULL))<0) // Pass NULL to only obtain number of devices
		error();
	fprintf(fp,"SDR Devices found: %i \n",n);
	if(n<1)
		error();
	lms_info_str_t* list = new lms_info_str_t[n]; // allocate device list storage
	// LimeSDR Mini, media=USB 3.0, module=FT601, addr=...
	// LimeSDR-USB, media=USB 3.0, module=FX3, addr=...
	// LimeNet-Micro, media=USB 2.0, module=FT601, addr=...
	if(LMS_GetDeviceList(list)<0)
		error();
	for(ci=0;ci<n;ci++) // print device list
	{
		printf("%i:%s\n",ci,list[ci]); // char[256]
		fprintf(fp,"%i:%s\n",ci,list[ci]); // char[256]
	}
	if(LMS_Open(&device,list[0],NULL)) //Open the first device
		error();
	delete [] list;
 
	const lms_dev_info_t *devInfo=LMS_GetDeviceInfo(device);
	sprintf(sdrDetails,"DEV: %s HW: %s FW: %s GW: %s S/N: %llu LIB: %s",devInfo->deviceName,devInfo->hardwareVersion,devInfo->firmwareVersion,devInfo->gatewareVersion,(unsigned long long)devInfo->boardSerialNumber,LMS_GetLibraryVersion());
	printf("%s\n",sdrDetails);
	char SDRminiStr[32]="LimeSDR-Mini"; // LNAH + LNAW
	char SDRmicroStr[32]="LimeNET-Micro"; // LNAH + LNAL
	if( (LNAnum==2) && (strstr(devInfo->deviceName,SDRminiStr)!=NULL) )
		LNAnum=3; // LNAW, LNAL not supported on this model
	if( (LNAnum==3) && (strstr(devInfo->deviceName,SDRmicroStr)!=NULL) )
		LNAnum=2; // LNAW, LNAL not supported on this model 
	
	if(LMS_Reset(device))
		error();
	if(LMS_Init(device) != 0) // Initialize device with default configuration
		error();
	// Use LMS_LoadConfig(device, "/path/to/file.ini") to load config from INI
	usleep(10000);

	if(LMS_EnableChannel(device, LMS_CH_RX,Ch,true)!=0) // Rx, Channels 0:(N-1)
		error();
	if(LMS_SetLOFrequency(device, LMS_CH_RX,Ch,(frq*1.0e6))!=0)
		error();
	if(((n=LMS_GetAntennaList(device, LMS_CH_RX, 0, antenna_list))) < 0)
		error(); // Alternatively, NULL can be passed to LMS_GetAntennaList() to obtain antennae num
	fprintf(fp,"SDR Ae: "); // print available antennae names
	for(ci = 0; ci < n; ci++)
		fprintf(fp," %i:%s",ci,antenna_list[ci]);
	if((n = LMS_GetAntenna(device, LMS_CH_RX, Ch)) < 0) // get selected antenna index
		error();
	fprintf(fp,", Default: %i:%s, ",n,antenna_list[n]);
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
	fprintf(fp,"Requested:%i Selected: %i:%s\n",LNAnum,n,antenna_list[n]);
	if(LMS_SetSampleRate(device,frq_Samp*1.0e6,OSR) != 0) // oversampling in RF OSR x sample rate, typedef double float_type;
		error();
	if(LMS_GetSampleRate(device, LMS_CH_RX, 0, &rate, &rf_rate) != 0)  // can pass NULL, typedef double float_type;
		error();
	frq_LPF=1.6*frq_Samp; // frq_LPF=frq_Samp does not give enough flat bandwidth for a scan. 1.6 gives LPF flat = SampRate
	if(LMS_SetLPFBW(device, LMS_CH_RX,Ch,frq_LPF*1.0e6)!=0)  //Configure LPF, bandwidth 8 MHz
		error(); 
	int err=0;
	if(LMS_SetGaindB(device,LMS_CH_RX,Ch,gaindB)!=0) // 0:73
	{
		fprintf(fp,"Ch=%i gaindB=%i err=%i\n",Ch,gaindB,err);
		error();
	}
	if(LMS_GetNormalizedGain(device,LMS_CH_RX,Ch,&gain)!=0) //normalized gain
		error();
	if(LMS_GetGaindB(device,LMS_CH_RX,Ch,&gaindB)!=0)
		error();
	// minimise noise figure
	SDRupdateSPIreg( (uint32_t)0x0113,1,0,3 ); // RFE:TIA Gain=2
	SDRupdateSPIreg( (uint32_t)0x0110,9,5,18 ); // RFE:LNAbias=18 (750uA)
	// minimise phase noise
	SDRupdateSPIreg( (uint32_t)0x0123,11,8,3 ); // SXR:CP2=3 (7.0pF)
	SDRupdateSPIreg( (uint32_t)0x0123,7,4,2 ); // SXR:CP3=2 (11.7pF)
	SDRupdateSPIreg( (uint32_t)0x0122,11,6,10 ); // SXR:IOFF=10 (2.4uA)
	SDRupdateSPIreg( (uint32_t)0x0122,5,0,32 ); // SXR:ICP=40 (76uA)
	SDRupdateSPIreg( (uint32_t)0x0120,7,0,150 ); // SXR:VCObias=150

	if( LMS_GetChipTemperature(device,0,&temp)!=0 )
		error();
	temp7002=(float)temp;
	fprintf(fp,"SDR USB rate: %.3f MS/s, ADC rate: %.3fMS/s OSR: %i CLKGEN=%.3f LPF=%.3fMHz\n", (rate*1.0e-6), (rf_rate*1.0e-6),OSR,(4*rf_rate*1.0e-6),frq_LPF);
	fprintf(fp,"SDR Normalized RX Gain: %f, RX Gain: %i dB\n\n",gain,gaindB);
	streamId.channel=Ch; //channel number
	streamId.fifoSize=FIFO_SIZE; //fifo size in samples uint32_t
	streamId.throughputVsLatency=1.0; //optimize for max throughput
	streamId.isTx=false; //RX channel
#ifdef USE_F32
	streamId.dataFmt=lms_stream_t::LMS_FMT_F32; //32-bit floats
#endif
#ifdef USE_I16
	streamId.dataFmt=lms_stream_t::LMS_FMT_I16; //16-bit int
#endif
#ifdef USE_I12
	streamId.dataFmt=lms_stream_t::LMS_FMT_I12; //12-bit int
#endif
	if(LMS_SetupStream(device,&streamId)!=0)
		error();
	LMS_StartStream(&streamId);
	fclose(fp);
}

void CloseSDR( void )
{
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

unsigned int GrabSamples( FILE *fp,int FrmSize,unsigned char frms,double frq ) // <= SRCH_FRMS
{
#ifdef USE_F32
	float _Complex *ptr=iqDataPtr;
#else
	int16_t *ptr=iqDataPtr;
#endif
	float fifosizelog[SRCH_FRMS];
	uint32_t droplog[SRCH_FRMS];
	uint64_t timeStamp[SRCH_FRMS]; // make ARM resistant! Intel unsigned long, ARM unsigned long long
	lms_stream_status_t status;
	clock_t mytref2=clock();	
	unsigned char cp=0;
	int err=0;
	int samplesRead; // size_t
	uint32_t bad=0; // uint32_t same as unsigned int for both ARM and Intel64
	for( cp=0; cp<frms; cp++ ) // avoid buffer overwrite, read one frame at a time (SRCH_FRMS=16)
	{ // each frame is 10ms, need at least 5 frames for complete SIB1
		samplesRead=LMS_RecvStream(&streamId,ptr,FrmSize,NULL,1000);
		err=LMS_GetStreamStatus(&streamId,&status); // most of the fields are uint32, except timestamp
		timeStamp[cp]=status.timestamp;
		droplog[cp]=status.droppedPackets;
		fifosizelog[cp]=(100.0*status.fifoFilledCount)/status.fifoSize;
		bad=status.droppedPackets;
		if(bad!=((uint32_t)0))
		{
			SDRreadStatus(fp,FrmSize,frq,(cp+1),mytref2,status,timeStamp,fifosizelog,droplog);
			if(bad>10000) // ARM/Pi fault, restart stream
			{ // timestamp reset to 0 when stream restarted.
				printf("ERtempROR: Stream/USB fault - restart stream\n");
				fprintf(fp,"ERROR: Stream/USB fault - restart stream\n");
				LMS_StopStream(&streamId);
				usleep(1000); // 1ms
				LMS_StartStream(&streamId);
			}
			return((unsigned int)bad);
		}	
	#ifdef USE_F32
		ptr+=FrmSize;
	#else
		ptr+=2*FrmSize;
	#endif
	}
	SDRreadStatus(fp,FrmSize,frq,frms,mytref2,status,timeStamp,fifosizelog,droplog);
	return((unsigned int)bad);
}

void SDRreadStatus( FILE *fp,int FrmSize,double frq,unsigned char frms, clock_t mytref2,lms_stream_status_t status,uint64_t *timeStamp,float *fifosizelog,uint32_t *droplog )
{
	unsigned char cp=0;
#ifdef USE_F32
	fprintf(fp"%8.3fMHz F32 Time to read %i*%i samples is %.3fms, expected is %.3fms FIFO=%i\n",frq,SRCH_FRMS,FrmSize,((double)(clock()-mytref2))/CLOCKS_PER_mSEC,10.0*frms,status.fifoSize);
#endif
#ifdef USE_I16
	fprintf(fp,"%8.3fMHz I16 Time to read %i*%i samples is %.3fms, expected is %.3fms FIFO=%i\n",frq,SRCH_FRMS,FrmSize,((double)(clock()-mytref2))/CLOCKS_PER_mSEC,10.0*frms,status.fifoSize);
#endif
#ifdef USE_I12
	fprintf(fp,"%8.3fMHz I12 Time to read %i*%i samples is %.3fms, expected is %.3fms FIFO=%i\n",frq,SRCH_FRMS,FrmSize,((double)(clock()-mytref2))/CLOCKS_PER_mSEC,10.0*frms,status.fifoSize);
#endif
	fprintf(fp,"%8.3fMHz FIFO: ",frq);
	for( cp=0; cp<frms; cp++ )
		fprintf(fp,"%.3f%% ",fifosizelog[cp]);
	fprintf(fp," DROP: ");
	for( cp=0; cp<frms; cp++ )
		fprintf(fp,"%u ",droplog[cp]);
	fprintf(fp," TIME: ");
	for( cp=0; cp<frms; cp++ )
		fprintf(fp,"%llu ",(unsigned long long)timeStamp[cp]); // make ARM compatible, coerce uint_64t to unsigned long long
	fprintf(fp,"\n");
}

short SetFreqGrabFrames( const char *fname,double frq,int FrmSize,unsigned char frms )
{
	unsigned int bad=0;
	unsigned short itt=0;
	FILE *fp;
	fp=fopen( fname,"a" );
	clock_t mytref = clock();
	float_type temp=0.0;
	LMS_StopStream(&streamId); // stream is stopped, start again with LMS_StartStream()
	if(LMS_SetLOFrequency(device,LMS_CH_RX,0,(frq*1.0e6+fUserCor))!= 0) // typedef double float_type;
		error(); // takes about 100us for PLL to settle.  Frames are 10ms long.
	usleep(1000); // 1ms
	if(LMS_Calibrate(device,LMS_CH_RX,Ch,floor(frq_Samp)*1.0e6,0)!=0) // => Rx calibration finished
		error();
	usleep(10000); // Tried 100ms
	if( LMS_GetChipTemperature(device,0,&temp)!=0 )
		error();
	temp7002=(float)temp;
	fprintf(fp,"LMS7002M=%.1foC\n",(double)temp7002);
	fprintf(fp,"%s\n",cpuTemp);
	LMS_StartStream(&streamId);
	bad=GrabSamples(fp,FrmSize,SRCH_FRMS,frq);
	bad=1;
	while((bad>0) && ((itt++)<MAX_SDR_READS) )
		bad=GrabSamples(fp,FrmSize,SRCH_FRMS,frq);
	if(itt>=MAX_SDR_READS)
	{
		fprintf(fp,"SDR: MAX_SDR_READS exceeded %i/%i.. Bad=%u\n",(itt-1),MAX_SDR_READS,bad);
		printf("SDR: MAX_SDR_READS exceeded %i/%i.. Bad=%u\n",(itt-1),MAX_SDR_READS,bad);
		fclose(fp);
		return(-1); // was -1
	}		
	fprintf(fp,"SDR cal/read took %.1fms\n",((double)(clock()-mytref))/CLOCKS_PER_mSEC);
	fclose(fp);
	return(1);
}

