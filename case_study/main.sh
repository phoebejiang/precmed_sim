#!/bin/bash

# Allow each combination of arguments to run separately

## Argument 1: PM Method Group
## default=('allDMF' 'allGA' 'linear' 'negBin' 'dWOLS' 'poisson' 'boosting' 'twoReg' 'contrastReg' 'listDTR2' 'listDTR3')
methods=('allDMF' 'allGA' 'linear' 'negBin' 'dWOLS' 'poisson' 'boosting' 'twoReg' 'contrastReg' 'listDTR2' 'listDTR3')
len=${#methods[@]}
stop=$(($len-1))

## Argument 2: Total number of CV repetition batches
batch=5

# Loop through each PM method
for ((j=0;j<=$stop;j+=1)); do
	m=${methods[$j]} 
	echo $m

	# Loop through each CV repetition batch
	for ((b=1;b<=$batch;b+=1)); do
		echo $b
		sed "s/method/$m/g" main.sl.temp > main.sl
		sed -i "s/batch/$b/g" main.sl
  		sbatch main.sl $m $b
	done
	
done