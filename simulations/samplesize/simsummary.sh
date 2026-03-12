#!/bin/bash

# Allow each combination of arguments to run separately

## Argument 1: Size of data 
## default=(500 1000 2500 5000 10000)
ns=(500 1000 2500 5000 10000)
len1=${#ns[@]}
stop1=$(($len1-1))

## Argument 2: level of treatment effect heterogeneity
betas=("c(-0.2,-0.2,-0.2,-0.2,-0.2)" "c(-0.36,-0.29,0,0.05,0.1)" "c(-0.92,-0.69,0,0.1,0.18)" "c(-1.2,-0.69,0,0.1,0.41)")  # ("c(-1,-0.5,0,0.5,1)")
  # ("c(-1,-0.5,0,0.5,1)")
# ("c(-0.2,-0.2,-0.2,-0.2,-0.2)" 
#  "c(log(0.7), log(0.75), log(1), log(1.05), log(1.1)" "c(-0.36,-0.29,0,0.05,0.1)"
#  "c(log(0.4), log(0.5), log(1), log(1.1), log(1.2))" "c(-0.92,-0.69,0,0.1,0.18)" 
#  "c(log(0.3), log(0.5), log(1), log(1.1), log(1.5))" "c(-1.2,-0.69,0,0.1,0.41)"
#  ")
len2=${#betas[@]}
stop2=$(($len2-1))

## Argument 3: responder subgroup percentile
percs=("seq(0,1,by=0.2)") # c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
len3=${#percs[@]}
stop3=$(($len3-1))

# Loop through each sample size 
for ((k=0;k<=$stop1;k+=1)); do
	n=${ns[$k]}
	echo $n

    # Loop through each level of heterogeneity
	for ((j=0;j<=$stop2;j+=1)); do
		beta=${betas[$j]}
    echo $beta

    # Loop through each subgroup percentiles
    for ((i=0;i<=$stop3;i+=1)); do
    	perc=${percs[$i]}
    	echo $perc

			sed "s/sized/$n/" simsummary.sl.temp > simsummary.sl
    	sed -i "s/beta/$beta/" simsummary.sl
    	sed -i "s/perc/$perc/" simsummary.sl
    	
	    sbatch simsummary.sl $n $beta $perc
   
    done
  done
done
