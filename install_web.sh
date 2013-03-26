#!/bin/bash

cd model_builder/C/core/

sed -ie "s/#define FLAG_JSON \([0-9]*\)/#define FLAG_JSON 1/" plom.h

cd ..

##compile with gcc-4.7 (on OSX, gcc shipped with Xcode is broken so brew install gcc-4.7)
mygcc=gcc-4.7
#if gcc-4.7 doesn't exists, then use gcc
type gcc-4.7 >/dev/null 2>&1 || mygcc=gcc;

make -f Makefile_web CC=$mygcc
make -f Makefile_web install
