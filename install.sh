#!/bin/bash

##NOTES: for the webApp: set FLAG_JSON to 1
## on OSX, gcc shipped with Xcode is broken so brew install gcc-4.7

cd model_builder/C/core

sed -ie "s/#define FLAG_JSON \([0-9]*\)/#define FLAG_JSON 0/" plom.h

make CC=gcc-4.7 ##remove CC=gcc-4.7 on linux
make install
cd ../

paths=( "simplex" "smc" "worker" "mif" "pmcmc" "simulation" "kalman" )

for i in "${paths[@]}"
do
    echo building $i
    cd $i
    make CC=gcc-4.7 obj ##remove CC=gcc-4.7 on linux
    cd ..
done

cd ../../


if [ -d "dist" ]; then
    rm -rf dist
fi

python setup.py sdist

cd dist
tar -zxvf plom-0.1.0.tar.gz
cd plom-0.1.0
python setup.py install
