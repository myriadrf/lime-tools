# Lime Tools

A collection of simple, useful tools for use with LimeSDR hardware.

## Contents

### LimeMon

Min, max, rms power sweep measurements vs. frequency and time in .csv, .xls and .gif format.

### LimeScan

Power spectrum tool that aims to be compatible with rtl_power and soapy_power.

## Dependencies

* [LimeSuite](http://wiki.myriadrf.org/Lime_Suite)
* GNUplot (`sudo apt-get install gnuplot`)
* FFTW3 (`sudo apt-get install libfftw3-dev`)

## CMake build instructions

    mkdir build
    cd build
    cmake ..
    make
    sudo make install

## Licensing

This software is published under the Apache License 2.0.
