#!/usr/bin/env bash

. utils.sh

RESULTS_CSV=test_case_2_results.csv

echo 'Command,Number of executors,Executor cores,Executor memory,Exit code,Time (mins)' > $RESULTS_CSV

for num_exec in 4 8 16 32 64
do
  time_gatk "CountReadsSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878.bam" $num_exec 1 1g >> $RESULTS_CSV
done

for num_exec in 2 4 8 16 32
do
  time_gatk "CountReadsSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878.bam" $num_exec 2 1g >> $RESULTS_CSV
done

for num_exec in 1 2 4 8 16
do
  time_gatk "CountReadsSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878.bam" $num_exec 4 1g >> $RESULTS_CSV
done