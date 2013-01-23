#!/bin/bash


echo "building noise example..."
cd noise
pmbuilder context.json process.json link.json -o model && fit adapt theta.json context.json process.json link.json -o model/theta.json
cd model
fit theta | ./smc deter
fit theta | ./smc sto
fit theta | ./simplex -M 1000 --no_trace
fit theta | ./kalman deter

cd ../../

echo "building drift example..."
cd drift
pmbuilder context.json process.json link.json -o model && fit adapt theta.json context.json process.json link.json -o model/theta.json
cd model
fit theta | ./smc deter
fit theta | ./smc sto
fit theta | ./ksimplex deter -M 100 --no_trace
fit theta | ./kalman deter

cd ../../


echo "clean"

rm -r noise/model drift/model
