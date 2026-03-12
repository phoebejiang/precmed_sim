#!/bin/bash

# Run this after simsummary.sh

# Do this once to copy the all 20% proportion group results ran from the samplesize simulations to the proportion simulations
# Samplesize simulations already did the medium HTE with 20% all equal proportion so no need to regenerate the same results
# START
# cp /home/pjiang/pmms/simulations/samplesize/outputs/simulations_n500_stratified10foldCV/*-0.92*.csv /home/pjiang/pmms/simulations/proportion/outputs/
# cp /home/pjiang/pmms/simulations/samplesize/outputs/simulations_n1000_stratified10foldCV/*-0.92*.csv /home/pjiang/pmms/simulations/proportion/outputs/
# cp /home/pjiang/pmms/simulations/samplesize/outputs/simulations_n2500_stratified10foldCV/*-0.92*.csv /home/pjiang/pmms/simulations/proportion/outputs/
# cp /home/pjiang/pmms/simulations/samplesize/outputs/simulations_n5000_stratified10foldCV/*-0.92*.csv /home/pjiang/pmms/simulations/proportion/outputs/
# cp /home/pjiang/pmms/simulations/samplesize/outputs/simulations_n10000_stratified10foldCV/*-0.92*.csv /home/pjiang/pmms/simulations/proportion/outputs/
# END

ns=(500 1000 2500 5000 10000)
len=${#ns[@]}
stop=$(($len-1))

# Loop through each sample size
for ((i=0;i<=$stop;i+=1)); do
  
  n=${ns[$i]} 
  echo $n
  sed "s/size/$n/g" summary_proportion.sl.temp > simsummary_proportion.sl
  sbatch simsummary_proportion.sl $n
         
done 
            
