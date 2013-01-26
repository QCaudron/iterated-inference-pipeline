#!/bin/bash


echo "building noise example..."
cd noise
pmbuilder context.json process.json link.json -o model && fit adapt theta.json context.json process.json link.json -o model/theta.json
cd model
fit theta | ./smc deter
fit theta | ./smc sto

fit theta | ./simplex -M 1000 --no_trace

fit theta -B | ./pmcmc deter -J 1 -M 1000
fit theta -B -C | ./pmcmc deter -J 1 -M 1000 -c --full


fit theta -B -o mle.json
fit theta mle.json | ./worker deter -J 50 -P 2 &
fit theta mle.json | ./worker deter -J 50 -P 2 &
fit theta mle.json | ./pmcmc deter -J 100 -M 100 -P 1 -C 50 -Z


fit theta -B | ./mif deter -J 1000 -M 50 -b 4

fit theta | ./kalman deter
fit theta | ./kalman sto
fit theta | ./ksimplex deter -M 100 --no_trace

cd ../../

echo "building drift example..."
cd drift
pmbuilder context.json process.json link.json -o model && fit adapt theta.json context.json process.json link.json -o model/theta.json
cd model

fit theta | ./smc deter
fit theta | ./smc sto
fit theta | ./kalman deter
fit theta | ./kalman sto

fit theta | ./ksimplex deter -M 100 --no_trace

fit theta -B | ./pmcmc deter -J 1 -M 1000
fit theta -B -C | ./pmcmc deter -J 1 -M 1000 -c --full

fit theta -B | ./mif deter -J 1000 -M 50 -b 4

cd ../../


echo "clean"

#rm -r noise/model drift/model
