#!/bin/bash

cd model_builder/C

paths=("core" "simplex" "smc" "worker" "mif" "pmcmc" "simulation" "kalman" )

for i in "${paths[@]}"
do
    echo cleaning $i
    cd $i
    make clean
    cd ..
done

rm -r lib/

cd ../../
