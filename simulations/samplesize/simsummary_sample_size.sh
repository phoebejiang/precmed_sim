#!/bin/bash

# Allow each combination of arguments to run separately

## Argument 1: level of treatment effect heterogeneity
betas=("c(-0.2,-0.2,-0.2,-0.2,-0.2)" "c(-0.36,-0.29,0,0.05,0.1)" "c(-0.92,-0.69,0,0.1,0.18)" "c(-1.2,-0.69,0,0.1,0.41)")  
len1=${#betas[@]}
stop1=$(($len1-1))

## Argument 2: responder subgroup percentile
percs=("seq(0,1,by=0.2)") 
len2=${#percs[@]}
stop2=$(($len2-1))

# Loop through each level of heterogeneity
for ((j=0;j<=$stop1;j+=1)); do
  beta=${betas[$j]}
  echo $beta

  # Loop through each subgroup percentiles
  for ((i=0;i<=$stop2;i+=1)); do
    perc=${percs[$i]}
 	  echo $perc

    sed "s/beta/$beta/" simsummary_sample_size.sl.temp > simsummary_sample_size.sl
    sed -i "s/perc/$perc/" simsummary_sample_size.sl 
    
    sbatch simsummary_sample_size.sl $beta $perc
   
  done
done

