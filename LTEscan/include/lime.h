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
#define MXD_CMPLR 1 // e.g. gcc,g++ or clang,clang++

// comment out to force test
#if __GNUC__ <6
#define USE_C99_COMPLEX 1 // assumes fftw_complex maps to C99 double _Complex (Ubuntu 16.04.5 gcc 5.4.0)
// #warning "GNU <6"
#else
#define USE_FFTW_COMPLEX 2 // assumes typedef double fftw_complex[2]; (Raspbian 6.0.3 gcc 6.3.0)
// #warning "GNU >=6"
#endif
// Also see same problem on Ubuntu 18.04.1 with gcc 7.3.0.  Version 8.1.0 available
// cannot easily see how we can automatically detect correct type from <fftw3.h>
// - see 4.1.1 of fftw3.pdf on C99 compatibility with FFTW3, only honoured in gcc 5?

// SDR data format.  For LimeNet Micro (USB2/ARM) use I16 or I12
//#define USE_F32 // bigest USB packets
#define USE_I16 // minimises post processing USB packets
//#define USE_I12 // smallest USB packets

#define USE_GNUPLOT 1 // comment out to disable GNUPLOT graphs

#define CLOCKS_PER_uSEC (CLOCKS_PER_SEC/1000000)
#define CLOCKS_PER_mSEC (CLOCKS_PER_SEC/1000)

#define KBLK 16 // was 8, how coarse to search Frame  OFDM symbol is CP+NFFT, approx 1.1*NFFT

#define PI 3.14159265f
#define PI2 (2*PI)
#define PId2 (PI/2)
#define PId4 (PI/4)
#define PI3d2 (3*PId2)
#define FIFO_SIZE (1048*1024) // was 1024 x1024, uint32_t
#define SRCH_FRMS 5 // number of frames to search for PBCH, was 16 for SIB1
// need to access up to 4.5 frames of psch,ssch to find a valid PBCH0 e.g.
//..SSC,PSC..SSC,PSCH,PBCH1..SSC,PSC..SSC,PSCH,PBCH2..SSC,PSC..SSC,PSCH,PBCH3..SSC,PSC..SSC,PSCH,PBCH0..
// as pbch and ssch are used to equalise pbch in the presence of a frequency error
#define SRCH_HFRMS (SRCH_FRMS*2-1)
#define FERR_COR_HFRMS (SRCH_FRMS*2-1)
#define rFERR_COR_HFRMS (1.0/FERR_COR_HFRMS)
#define rFERR_COR_HFRMSm1 (1.0/(FERR_COR_HFRMS-1))
#define FERR_COR_SYMS (6*FERR_COR_HFRMS) // 9*(PSCH+SSCH+4*PBCH)
#define MAX_SDR_READS 5 // maxinum number of tries to try to read SDR

