#!/bin/bash

# Allow each combination of arguments to run separately

## Argument 1: PM Method
#default=('allA1' 'allA0' 'linear' 'weightedLinear' 'weightedPoisson' 'negBin' 'weightedNegBin' 'dWOLS' 'listDTR2' 'listDTR3' 'poisson' 'boosting' 'twoReg' 'contrastReg')
methods=('allA1' 'allA0' 'linear' 'negBin' 'dWOLS' 'listDTR2' 'listDTR3' 'poisson' 'boosting' 'twoReg' 'contrastReg')
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
# betas=("c(-0.92,-0.69,0,0.1,0.18)" "c(-0.92,-0.69,0,0.1,0.18)" "c(-0.92,-0.69,-0.69,-0.69,0)" "c(-0.92,-0.69,-0.69,0,0.1)")  
betas=("c(-0.2,-0.2,-0.2,-0.2,-0.2)" "c(-0.92,-0.69,0,0.1,0.18)" "c(-0.92,-0.69,0,0.1,0.18)" "c(-0.92,-0.69,-0.69,-0.69,0)" "c(-0.92,-0.69,-0.69,0,0.1)")  
#  "c(log(0.7), log(0.75), log(1), log(1.05), log(1.1)" "c(-0.36,-0.29,0,0.05,0.1)" Low
#  "c(log(0.4), log(0.5), log(1), log(1.1), log(1.2))" "c(-0.92,-0.69,0,0.1,0.18)" Medium + symmetric
#  "c(log(0.4), log(0.5), log(0.5), log(0.5), log(1))" "c(-0.92,-0.69,-0.69,-0.69,0)" Medium + asymmetric1 (55%, 15%, 15%, 15%)
#  "c(log(0.4), log(0.5), log(0.5), log(1), log(1.1))" "c(-0.92,-0.69,-0.69,0,0.1)" Medium + asymmetric2  (55%, 30%, 15%)
#  "c(log(0.3), log(0.5), log(1), log(1.1), log(1.5))" "c(-1.2,-0.69,0,0.1,0.41)" High
# )
echo $betas 
  
## Argument 6: responder subgroup percentile
# percs=("c(0,0.1,0.25,0.75,0.9,1)" "c(0,0.1,0.4,0.6,0.9,1)" "c(0,0.55,0.65,0.75,0.85,1)" "c(0,0.55,0.65,0.7,0.85,1)") 
percs=("c(0,0.1,0.25,0.75,0.9,1)" "c(0,0.1,0.4,0.6,0.9,1)" "c(0,0.55,0.65,0.75,0.85,1)" "c(0,0.55,0.65,0.7,0.85,1)")
# c(0,0.2,0.4,0.6,0.8,1) 
# c(0,0.1,0.25,0.75,0.9,1) => symmetric (10%, 15%, 50%, 15%, 10%)
# c(0,0.1,0.4,0.6,0.9,1) => symmetric (10%, 30%, 20%, 30%, 10%)
# c(0,0.55,0.85,1,1,1) => c(0,0.55,0.65,0.75,0.85,1) => asymmetric1 (55%, 15%, 15%, 15%)
# c(0,0.55,0.7,0.85,1,1) => c(0,0.55,0.65,0.7,0.85,1) => asymmetric2  (55%, 30%, 15%)
echo $percs
len4=${#percs[@]}
stop4=$(($len4-1))


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
  
  echo yvar

	# Loop through each sample size 
	for ((k=0;k<=$stop3;k+=1)); do
		n=${ns[$k]}
		batch=${batches[$k]}

		echo $n
		echo $batch
			
		# Loop through each CV repetition batch
		for ((b=1;b<=$batch;b+=1)); do

			echo $b
			
			# Loop through each proportion
			for ((p=0;p<=$stop4;p+=1)); do
				perc=${percs[$p]}
				beta=${betas[$p]}
				echo $perc $beta

				sed "s/method/$m/g" simmain.sl.temp > simmain.sl
				sed -i "s/yvar/$yvar/g" simmain.sl
				sed -i "s/batch/$b/g" simmain.sl
				sed -i "s/sized/$n/g" simmain.sl
				sed -i "s/beta/$beta/g" simmain.sl
				sed -i "s/perc/$perc/g" simmain.sl

				sbatch simmain.sl $m $yvar $b $n $beta $perc

				if [ $? -ne 0 ]; then
					echo "cannot do more sbatch, first failed submission: simmain.sl $m $y $b $n $beta $perc" 
					echo "Terminating at this point, please re-run later from this step"
					exit 1
			 	fi
			 	
		    done
		 done
	done
done