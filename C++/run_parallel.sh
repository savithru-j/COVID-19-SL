#!/bin/bash

#Copy this file to the build/release directory and run from there.

if [ -z $1 ]
then
  echo "Required arguments: country_name [num_threads = 1]"
  exit
else
  echo "Country: $1"
fi

country=$1
wt_conf=1;
wt_recov=1;
wt_fatal=1;

node_interval=5
num_iter=1000;
num_passes=10;
num_threads=1;
if [ ! -z $2 ]
then
  num_threads=$(($2))
fi

echo "No. of threads: $num_threads"

for seed in $(seq 1 $num_threads)
do
  echo "Launching thread $seed"
  ./CovidOptPiecewise $country $wt_conf $wt_recov $wt_fatal $node_interval $num_iter $num_passes $seed > results/output$seed.txt &
done


