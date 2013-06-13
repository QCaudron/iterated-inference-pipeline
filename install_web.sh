#!/bin/bash

cd model_builder/C/core/

sed -ie "s/#define FLAG_JSON \([0-9]*\)/#define FLAG_JSON 1/" plom.h

cd ..

make -f Makefile_web
make -f Makefile_web install
