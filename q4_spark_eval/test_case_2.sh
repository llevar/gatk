#!/usr/bin/env bash

. utils.sh

RESULTS_CSV=test_case_2_results.csv

echo 'Command,Number of executors,Executor cores,Executor memory,Exit code,Time (mins)' > $RESULTS_CSV

run_id=0
for exec_cores in 1 2 4; do
  for total_cores in 4 8 16 32 64; do
    num_exec=$((total_cores * exec_cores))
    time_gatk "CountReadsSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878.bam" $num_exec $exec_cores 1g 1g >> $RESULTS_CSV
    run_id=$((run_id+1))
  done
done
