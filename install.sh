#!/bin/bash


./clean.sh

cd model_builder/C/core

##NOTES: for the webApp: set FLAG_JSON to 1
sed -ie "s/#define FLAG_JSON \([0-9]*\)/#define FLAG_JSON 0/" plom.h

##compile with gcc-4.7  (on OSX, gcc shipped with Xcode is broken so brew install gcc-4.7)
mygcc=gcc-4.7
#if gcc-4.7 doesn't exists, then use gcc
type gcc-4.7 >/dev/null 2>&1 || mygcc=gcc;

make CC=$mygcc
make install
cd ../

paths=( "simplex" "smc" "worker" "mif" "pmcmc" "simulation" "kalman" )

for i in "${paths[@]}"
do
    echo building $i
    cd $i
    make CC=$mygcc obj
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
