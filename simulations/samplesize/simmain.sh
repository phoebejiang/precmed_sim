#!/bin/bash

# Allow each combination of arguments to run separately

## Argument 1: PM Method
#default=('allA1' 'allA0' 'linear' 'weightedLinear' 'weightedPoisson' 'negBin' 'weightedNegBin' 'dWOLS' 'listDTR2' 'listDTR3' 'poisson' 'boosting' 'twoReg' 'contrastReg')
methods=('allA1' 'allA0' 'linear' 'negBin' 'dWOLS' 'listDTR2' 'poisson' 'boosting' 'twoReg' 'contrastReg')
len=${#methods[@]}
stop=$(($len-1))
  
## Argument 2: transformed y variable for PM method that cannot accept count outcomes
## default=("logarr0001" "logarr1" "logarr01r")
contiyvars=("logarr0001")

## Argument 3: Total number of CV repetition batches
## Should match with the order of `ns` below
## default = 5 for n = 100, 250, 500, 1000, 2000 with batch_size = 5, so a total of 25 CV repetitions
## new = 25, for n = 5000, 10000, with batch_size = 1, still a total of 25 CV repetitions
## default=(5 5 25 25 25)
batches=(5 5 25 25 25)

## Argument 4: Size of data (if n < 8599, then randomly select a subset of the full MarketScan data)
## default=(500 1000 2500 5000 10000)
ns=(500 1000 2500 5000 10000)
len3=${#ns[@]}
stop3=$(($len3-1))

## Argument 5: level of treatment effect heterogeneity
## default=("c(-0.2,-0.2,-0.2,-0.2,-0.2)" "c(-0.92,-0.69,0,0.1,0.18)" "c(-0.36,-0.29,0,0.05,0.1)" "c(-1.2,-0.69,0,0.1,0.41)")  
betas=("c(-0.2,-0.2,-0.2,-0.2,-0.2)" "c(-0.92,-0.69,0,0.1,0.18)" "c(-0.36,-0.29,0,0.05,0.1)" "c(-1.2,-0.69,0,0.1,0.41)")
# ("c(-0.2,-0.2,-0.2,-0.2,-0.2)" No 
#  "c(log(0.7), log(0.75), log(1), log(1.05), log(1.1)" "c(-0.36,-0.29,0,0.05,0.1)" Low 
#  "c(log(0.4), log(0.5), log(1), log(1.1), log(1.2))" "c(-0.92,-0.69,0,0.1,0.18)" Medium
#  "c(log(0.3), log(0.5), log(1), log(1.1), log(1.5))" "c(-1.2,-0.69,0,0.1,0.41)" High
#  ")
echo $betas 
len2=${#betas[@]}
stop2=$(($len2-1))
  
## Argument 6: responder subgroup percentile
perc=("seq(0,1,by=0.2)") # c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
echo $perc

# Loop through each PM method
for ((i=0;i<=$stop;i+=1)); do
	m=${methods[$i]} 
	
	echo $m
  
	# Determine outcome candidates for the current PM method
	if [[ $m = "weightedNegBin" ]] || [[ $m = "negBin" ]] || [[ $m = "weightedPoisson" ]] || \
	   [[ $m = "poisson" ]] || [[ $m = "boosting" ]] || [[ $m = "twoReg" ]] || \
	   [[ $m = "contrastReg" ]] || [[ $m = "allA1" ]] || [[ $m = "allA0" ]]; then 
		
		yvar=("postrelapse_num") # these methods can accept count outcomes
	
	else
	
		yvar=$contiyvars
	
	fi
 
 	echo $yvar

	# Loop through each sample size 
	for ((k=0;k<=$stop3;k+=1)); do
		n=${ns[$k]}
		batch=${batches[$k]}

		echo $n
		echo $batch
			
		# Loop through each CV repetition batch
		for ((b=1;b<=$batch;b+=1)); do
			echo $b
					
			# Loop through each beta variable
			for ((j=0;j<=$stop2;j+=1)); do
				beta=${betas[$j]}
				echo $beta
		
				sed "s/method/$m/g" simmain.sl.temp > simmain.sl
				sed -i "s/yvar/$yvar/g" simmain.sl
				sed -i "s/batch/$b/g" simmain.sl
				sed -i "s/sized/$n/g" simmain.sl
				sed -i "s/beta/$beta/g" simmain.sl
				sed -i "s/perc/$perc/g" simmain.sl
				
				sbatch simmain.sl $m $yvar $b $n $beta $perc
			
			done
		done
	done
done