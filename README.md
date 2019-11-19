# Lime Tools

A collection of simple, useful tools for use with LimeSDR hardware.

## Contents

### LimeMon

Min, max, rms power sweep measurements vs. frequency and time in .csv, .xls and .gif format.

### LimeScan

Power spectrum tool that aims to be compatible with rtl_power and soapy_power.

### LTEScan

LTEScan produces automated scan of LTE bands.

## Dependencies

* [LimeSuite](http://wiki.myriadrf.org/Lime_Suite)
* GNUplot (`sudo apt-get install gnuplot`)
* FFTW3 (`sudo apt-get install libfftw3-dev`)
* libsensors (`sudo apt-get install libsensors4-dev`)
* [turbofec](https://github.com/ttsou/turbofec)

## CMake build instructions

    mkdir build
    cd build
    cmake ..
    make
    sudo make install

## Licensing

This software is published under the Apache License 2.0.
