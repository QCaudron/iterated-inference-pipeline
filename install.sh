#!/bin/bash

cd model_builder/C/core/

sed -ie "s/#define FLAG_JSON \([0-9]*\)/#define FLAG_JSON 0/" plom.h
#NOTE on OSX openMP require to compile with gcc-4.7 or gcc-4.8
sed -ie "s/#define FLAG_OMP \([0-9]*\)/#define FLAG_OMP 0/" plom.h 

cd ..

##compile with gcc-4.8 (on OSX, gcc shipped with Xcode is broken so brew tap homebrew/versions && brew install gcc48)
mygcc=gcc-4.8
#if gcc-4.8 doesn't exists, then use gcc
type gcc-4.8 >/dev/null 2>&1 || mygcc=gcc;

make CC=$mygcc
make install

cd ../../

##sudo is necessary if using python 2.7 shipping with OSX 

if [ -d "dist" ]; then
    sudo rm -rf dist
fi

python setup.py sdist

cd dist
tar -zxvf plom-0.9.0.tar.gz
cd plom-0.9.0
sudo python setup.py install
