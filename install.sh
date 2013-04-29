#!/bin/bash

cd model_builder/C/core/

sed -ie "s/#define FLAG_JSON \([0-9]*\)/#define FLAG_JSON 0/" plom.h

cd ..

##compile with gcc-4.7 (on OSX, gcc shipped with Xcode is broken so brew install gcc-4.7)
mygcc=gcc-4.7
#if gcc-4.7 doesn't exists, then use gcc
type gcc-4.7 >/dev/null 2>&1 || mygcc=gcc;

make CC=$mygcc
make install

cd ../../

if [ -d "dist" ]; then
    rm -rf dist
fi

python setup.py sdist

cd dist
tar -zxvf plom-0.6.0.tar.gz
cd plom-0.6.0
python setup.py install
