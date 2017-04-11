#!/usr/bin/env bash

. utils.sh

RESULTS_CSV=test_case_4_results.csv

echo 'Command,Number of executors,Executor cores,Executor memory,Exit code,Time (mins)' > $RESULTS_CSV

for exec_cores in 1 2 4; do
  for total_cores in 4 8 16 32 64; do
    num_exec=$((total_cores * exec_cores))
    time_gatk "CountReadsSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878.bam -L Broad.human.exome.b37.interval_list" $num_exec $exec_cores 1g 1g >> $RESULTS_CSV
  done
done
