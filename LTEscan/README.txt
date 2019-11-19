Introdcution:
=============

Welcome to the LTEScan software.  LTEScan produces automated scan of LTE bands.
Software has command line interface compatible with LimeScan.

Dependecies:
============

autoconf (turbofec)
libtools (turbofec)
lm-sensors (intel temperature)
libsensors4-dev (intel temperature)
libfftw3-3 (fft)
xgnuplot (graphs)
xdg-open
turbofec
limesuite

Turbofec and LimeSuite gits should be cloned into the home directory.

To obtain limesuite, follow instructions on https://wiki.myriadrf.org/Lime_Suite
If you choose the repository version of LimeSuite, be sure to also include liblimesuite-dev 

To obtain turbofec
git clone https://github.com/ttsou/turbofec.git
cd turbofec
autoreconf -i
./configure
make
make check
sudo make install
sudo ldconfig

To compile:
===========
use lime-tools cmake build

To run:
=======

For best results use an external antenna.  

Fit antenna to RX input, and run the following script.  
./lteScan.sh

Else run as a command line with parameters, e.g.
./LTEscan -B 20 -n 2

This script contains various example commands.  The parameters for the LTEScan command are as follows:

-B List of LTE bands, separated by commas and no spaces e.g. 1,5,7,20  Not all bands may be available in your area.
-n Number of repeats {1,2,3,...} (default 2) (useful for low signal strength)
-b NFFT size {128,256,512,1024,512} (default 512)
-A Antenna name {LNAH,LNAL,LNAW} (default LNAL)
-C channel number {0,1} (default 0)
-f frequency correction in Hz.  Typically between -800Hz and 800Hz, correcting for TXCO error.
-h help
-v version
// -OSR Oversample rate for ADC {1,2,4,8,16,32} (do not alter). Max value depends on NFFT, 2048 OSR<=4, 1024 OSR<=8, 512 OSR<=16 256 OSR<=32 128 OSR<=32
// -g SDR gain (default 48 for best dynamic range, LNAgain=max, PGAgain=-3dB, TIA=max)
// --info Not implemented.
// -gps "/dev/ttyACM0" (supports UBLOX-7 NMEA GPRMC messages) - not currently included

Compatible with LimeSDR-USB, LimeSDR-mini and LimeNet Micro..
Note: LimeNet micro contains an RF switch that selects LNAH or LNAL
Note: LimeSDRmini contains an RF switch that selects LNAH or LNAW
Note: LimeSDRmini and LimeNet Micro only has channel 0 (-C 0)
Note: Aluminium LimeSDR also requires internal cable changes to select LNAH,LNAL,LNAW

Note: GPS Ublox-7 USB sticks, sudo adduser $USER dialout # then reboot
Note: GPS Ublox-7 USB sticks work best outside.  Accuracy will be degraded using indoors, and may become unusable.

Results in LTEscan.html and ./output/LTEscan.txt
Detailed diagnostics in ./output/LTElog.txt

