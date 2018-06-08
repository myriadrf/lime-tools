Introdcution:

Welcome to the LimeScan software.  Produces automated min,max,rms power sweep measurements vs frequency and time in .csv, .xls and .gif formats.  Software has compatible command line to Soapy Power.

Dependecies:

libfftw3-3 xgnuplot xdg-open limesuite
To obtain limesuite, follow instructions on https://wiki.myriadrf.org/Lime_Suite
If you choose the repository version of LimeSuite, be sure to also include liblimesuite-dev 

To compile:

make

To run:

./testScript.sh

This script contains various example commands, and the results are stored in separate directories.  The parameters for the LimeScan command are as follows.

-f Frequency range fmin:fmax (Hz) e.g. 800M:1000M
-C Channel number 0:1
-A Antenna name LNAH,LNAL,LNAW
-w IF filter bandwidth, e.g. 35M (Hz), should be twice as large as -r
-r Sample rate, e.g. 16M (S/s)
-OSR Oversample rate for ADC e.g. 8
-b FFT bins
-g SDR gain
-n Number of repeats for RMS FFT
-o -O Output folder name. e.g. -o TESTmini.
-h help
-v version
-T Measurement time (s)
-gps "/dev/ttyACM0" (supports UBLOX-7 NMEA GPRMC messages)
-tst LVL Generate test signal at 860MHz LVL dBm -50:-6 e.g. -test -50

Compatible with LimeSDR and LimeSDRmini.
Note: LimeSDRmini contains an RF switch that selects LNAH or LNAW
Note: LimeSDRmini only has channel 0 (-C 0)
Note: Aluminium LimeSDR requires cables to change to select LNAH,LNAL,LNAW

Note: GPS Ublox-7 USB sticks, sudo adduser $USER dialout # then reboot
Note: GPS Ublox-7 USB sticks work best outside.  Accuracy will be degraded using indoors.

Note: Testfrequency is 860MHz+SampleRate/8
