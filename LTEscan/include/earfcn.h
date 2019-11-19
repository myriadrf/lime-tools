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
double LTEbandDLLow(unsigned char band );
double LTEbandDLHigh(unsigned char band );
double LTEbandULLow(unsigned char band );
double LTEbandULHigh(unsigned char band );
double LTEbandBW(unsigned char band );
int LTEbandDLFreq2EARFCN(unsigned char band, double freq);
double LTEbandDLEARFCN2Freq(unsigned char band, int earfcn);


