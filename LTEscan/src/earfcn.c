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
#include <math.h>

// TS 136.101 13.0.0 Table 5.5-1, does not contain EARFCN numbers
// EARFCN = E-UTRA Absolute Radio Frequency Channel Number!
// TS 136.104 13.5.0 Table 5.7.3-1 contains EARFCN numbers.  100kHz steps
// Note 1, band edge, following Earfcns cannot be used
//  1.4MHz First   7, Last  6
//  3.0Mhz First  15, Last 14
//  5.0MHz First  25, Last 24
// 10.0MHz First  50, Last 49
// 15.0MHz First  75, Last 74
// 20.0MHz First 100, Last 99
// note band definitions may change slightly between various editions of the standard
// Note NBIOT channel numbers different
// Note band 46, and Carrier aggregation have further restrictions.
// Similar to ARFCN GSM 45.005, which uses 200kHz steps.
// DLmin, DLmax,DL-UL,DLBW,DLEarfcnMin,DLEarfcnMax,ULEarfcnMin,ULEarfcnMax,{FDD=1 TDD=2, DLonly=3}
double LTEbands[86][9]={
{0,0,0,0,0,0,0,0,0}, // no channel 0, pad to avoid silly subtractions
{2110,2170,190,60,0,599,18000,18599,1}, // 1
{1930,1990,80,60,600,1199,18600,19199,1}, // 2
{1805,1880,95,75,1200,1949,19200,19949,1}, // 3
{2110,2155,400,45,1950,2399,19950,20399,1}, // 4
{869,894,45,25,2400,2649,20400,20649,1}, // 5
{875,885,45,10,2650,2749,20650,20749,1}, // 6
{2620,2690,120,70,2750,3449,20750,21449,1}, // 7
{925,960,45,35,3450,3799,21450,21799,1}, // 8
{1844.9,1879.9,95,35,3800,4149,21800,22149,1}, // 9
{2110,2170,-400,60,4150,4749,22150,22749,1}, // 10
{1475.9,1495.9,48,20,4750,4949,22750,22949,1}, // 11
{729,746,30,17,5010,5179,23010,23179,1}, // 12
{746,768,-31,10,5180,5279,23180,23279,1}, // 13
{758,768,-30,10,5280,5379,23280,23379,1}, // 14

{0,0,0,0,0,0,0,0,0}, // 15 Reserved
{0,0,0,0,0,0,0,0,0}, // 16 Reserved

{734,746,30,12,5730,5849,23730,23849,1}, // 17
{860,875,45,15,5850,5999,23850,23999,1}, // 18
{875,890,45,15,6000,6149,24000,24149,1}, // 19
{791,821,-41,30,6150,6449,24150,24449,1}, // 20
{1495.9,1510.9,48,15,6450,6599,24450,24599,1}, // 21
{3510,3590,100,80,6600,7399,24600,25399,1}, // 22
{2180,2200,180,20,7500,7699,25500,25699,1}, // 23
{1525,1559,-101.5,34,7700,8039,25700,26039,1}, // 24
{1930,1995,80,65,8040,8689,26040,26689,1}, // 25
{859,894,45,35,8690,9039,26690,27039,1}, // 26
{852,869,45,17,9040,9209,27040,27209,1}, // 27
{758,803,55,45,9210,9659,27210,27659,1}, // 28
{717,728,0,11,9660,9769,0,0,3}, // 29
{2350,2360,45,10,9770,9869,27660,27759,1}, // 30
{462.5,467.5,10,5,9870,9919,27760,27809,1}, // 31
{1452,1496,0,44,9920,10359,0,0,3}, // 32

{1900,1920,0,20,36000,36199,2}, // 33
{2010,2025,0,15,36200,36349,2}, // 34
{1850,1910,0,60,36350,36949,2}, // 35
{1930,1990,0,60,36950,37549,2}, // 36
{1910,1930,0,20,37550,37749,2}, // 37
{2570,2620,0,50,37750,38249,2}, // 38
{1880,1920,0,40,38250,38649,2}, // 39
{2300,2400,0,100,38650,39649,2}, // 40
{2496,2690,0,194,39650,41589,2}, // 41
{3400,3600,0,200,41590,43589,2}, // 42
{3600,3800,0,200,43590,45589,2}, // 43
{703,803,0,100,45590,46589,2}, // 44
{1447,1467,0,20,46590,46789,2}, // 45
{5150,5925,0,775,46790,54539,2}, // 46
{5855,5925,0,70,54540,55239,2}, // 47
{3550,3700,0,150,55240,56739,2}, // 48
{3550,3700,0,150,56740,58239,2}, // 49
{1432,1517.5,0,85,58240,59089,2}, // 50
{1427,1432,0,5,59090,59139,2}, // 51
{3300,3400,0,100,59140,60139,2}, // 52

{0,0,0,0,0,0,0,0,0}, // 53
{0,0,0,0,0,0,0,0,0}, // 54
{0,0,0,0,0,0,0,0,0}, // 55
{0,0,0,0,0,0,0,0,0}, // 56
{0,0,0,0,0,0,0,0,0}, // 57
{0,0,0,0,0,0,0,0,0}, // 58
{0,0,0,0,0,0,0,0,0}, // 59
{0,0,0,0,0,0,0,0,0}, // 60
{0,0,0,0,0,0,0,0,0}, // 61
{0,0,0,0,0,0,0,0,0}, // 62
{0,0,0,0,0,0,0,0,0}, // 63
{0,0,0,0,0,0,0,0,0}, // 64 Reserved

{2110,2200,190,90,65536,66435,131072,131971,1}, // 65
{2110,2200,400,90,66436,67335,131972,132671,1}, // 66
{738,758,0,20,67336,67535,0,0,3}, // 67
{753,783,55,30,67536,67835,132672,132971,1}, // 68
{2570,2620,0,50,67836,68335,0,0,3}, // 69
{1995,2020,300,25,68336,68585,132972,133121,1}, // 70
{617,652,-46,35,68586,68935,133122,133471,1}, // 71
{461,466,10,5,68936,68985,133472,133521,1}, // 72
{460,465,10,5,68986,69035,133522,133571,1}, // 73
{1475,1518,48,43,69036,69465,133572,134001,1}, // 74
{1432,1517,0,85,69466,70315,0,0,3}, // 75
{1427,1432,0,5,70316,70365,0,0,3}, // 76

{0,0,0,0,0,0,0,0,0}, // 77
{0,0,0,0,0,0,0,0,0}, // 78
{0,0,0,0,0,0,0,0,0}, // 79
{0,0,0,0,0,0,0,0,0}, // 80
{0,0,0,0,0,0,0,0,0}, // 81
{0,0,0,0,0,0,0,0,0}, // 82
{0,0,0,0,0,0,0,0,0}, // 83
{0,0,0,0,0,0,0,0,0}, // 84

{728,746,30,18,70366,70545,134002,134181,1}}; // 85

double LTEbandDLLow(unsigned char band )
{
	if(band>85)
		return(0.0);
	return(LTEbands[band][0]);
}

double LTEbandDLHigh(unsigned char band )
{
	if(band>85)
		return(0.0);
	return(LTEbands[band][1]);
}

double LTEbandULLow(unsigned char band )
{
	if(band>85)
		return(0.0);
	return(LTEbands[band][0]-LTEbands[band][2]);
}

double LTEbandULHigh(unsigned char band )
{
	if(band>85)
		return(0.0);
	return(LTEbands[band][1]-LTEbands[band][2]);
}

double LTEbandBW(unsigned char band )
{
	if(band>85)
		return(0.0);
	return(LTEbands[band][3]);
}

int LTEbandDLFreq2EARFCN(unsigned char band, double freq)
{
	if(band>85)
		return(0);
	if( freq<LTEbands[band][0] )
		return(0);
	if( LTEbands[band][1]<freq )
		return(0);
	return( LTEbands[band][4]+(int)floor((freq-LTEbands[band][0]+0.05)/0.1) );
}

double LTEbandDLEARFCN2Freq(unsigned char band, int earfcn)
{
	if(band>85)
		return(0.0);
	if( earfcn<LTEbands[band][4] )
		return(0.0);
	if( LTEbands[band][5]<earfcn )
		return(0.0);
	return( LTEbands[band][0]+0.1*(earfcn-LTEbands[band][4]) );
}

