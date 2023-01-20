#!/bin/bash
# This is the code used for parallel multiple imputation
PID=()
for i in `seq 1 10` # create 10 threads of same R program
do
RNDseed=$RANDOM # get random number for the seed of multiple imputation algorithm
	echo $RNDseed
	Rscript Multiple_imputation_code.R $RNDseed $1 & # input 2 arguments: random seed and k (the offset) for sensitivity analysis
	PID[i]=$!
done

# Ensure that every threads created return their results
# before closing of the current program
for i in `seq 0 9`
do
	echo wait
	echo ${PID[i]}
        wait ${PID[i]}
done

